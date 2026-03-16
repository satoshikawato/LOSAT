//! BLASTP search coordination and execution
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c

use anyhow::{bail, Context, Result};
use bio::io::fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::blast_gapalign::{get_ungapped_hsp_list, InitHSP};
use crate::algorithm::tblastx::chaining::UngappedHit;
use crate::algorithm::tblastx::constants::{
    GAP_TRIGGER_BIT_SCORE, X_DROP_GAPPED_FINAL, X_DROP_GAPPED_PRELIM, X_DROP_UNGAPPED_BITS,
};
use crate::algorithm::tblastx::lookup::{build_ncbi_lookup, encode_kmer, QueryContext};
use crate::algorithm::tblastx::ncbi_cutoffs::{
    cutoff_score_from_evalue, gap_trigger_raw_score, x_drop_raw_score,
};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::common::Hit;
use crate::config::ScoringMatrix;
use crate::core::blast_diagnostics::diagnostics_enabled;
use crate::core::blast_stat::compute_blosum62_ideal_karlin_params;
use crate::core::composition_adjustment::adjust_scores::build_matrix_info;
use crate::core::composition_adjustment::redo_alignment::{
    BlastCompoAdjustMode, BlastCompoGappingParams,
};
use crate::report::{
    format_bitscore_ncbi, format_evalue_ncbi_tabular, write_pairwise, OutputFormat, PairwiseConfig,
    PairwiseHit, ReportContext,
};
use crate::stats::search_space::SearchSpace;
use crate::stats::{
    lookup_protein_gumbel_params, lookup_protein_params, lookup_protein_params_ungapped,
    protein_scoring_supported,
};
use crate::utils::matrix::protein_score;
use crate::utils::matrix::BLASTAA_SIZE;

use super::alignment::build_alignment_view_with_matrix;
use super::args::{BlastpArgs, BlastpCompositionMode, ResolvedBlastpArgs};
use super::encoding::{
    encode_protein_query_frame_with_seg, encode_protein_sequence, ncbistdaa_to_ascii,
    EncodedProtein,
};
use super::extension::{extend_one_hit, extend_two_hit};
use super::gapalign::{
    blastp_get_start_for_gapped_alignment, blastp_score_only_gapped_alignment,
    BlastpGappedAlignmentMode, BlastpPreliminaryHsp,
};
use super::hsp::{
    collect_hits_from_hit_lists, get_prelim_hitlist_size, trim_by_max_hsps, BlastpHitList,
    BlastpHsp, BlastpHspList,
};
use super::kappa::{postprocess_preliminary_hits, BlastRedoAlignParams};

#[inline]
fn fasta_id(record: &fasta::Record) -> Arc<str> {
    Arc::<str>::from(record.id().split_whitespace().next().unwrap_or("unknown"))
}

#[inline]
fn path_label(path: &Path) -> String {
    path.file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string())
}

#[inline]
fn fasta_defline(record: &fasta::Record) -> String {
    match record.desc() {
        Some(desc) => format!("{} {}", record.id(), desc),
        None => record.id().to_string(),
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:81-83
// ```c
// NCBI_XBLAST_EXPORT const int   kBlastMajorVersion = 2;
// NCBI_XBLAST_EXPORT const int   kBlastMinorVersion = 17;
// NCBI_XBLAST_EXPORT const int   kBlastPatchVersion = 0;
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/api/version.hpp:57-60
// ```c
// virtual string Print(void) const {
//     return CVersionInfo::Print() + "+";
// }
// ```
const NCBI_BLASTP_VERSION: &str = "2.17.0+";
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3698-3699
// ```c
// const double kRestrictedMult = 0.68;/* fraction of the ordinary cutoff score
// ```
const BLASTP_RESTRICTED_MULT: f64 = 0.68;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:672-676
// ```c
// /** SCALING_FACTOR is a multiplicative factor used to get more bits of
//  * precision in the integer matrix scores.
//  */
// #define SCALING_FACTOR 32
// ```
const BLASTP_KAPPA_SCALING_FACTOR: f64 = 32.0;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2402-2404
// ```c
// /* Bit score per alignment position threshold for preliminaru near identical
//    test */
// #define NEAR_IDENTICAL_BITS_PER_POSITION (1.74)
// ```
const BLASTP_NEAR_IDENTICAL_BITS_PER_POSITION: f64 = 1.74;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2387-2415
// ```c
// gapping_params->x_dropoff = (Int4)
//     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
//         extendParams->gap_x_dropoff_final);
// ```
fn blastp_redo_x_dropoff(
    scaled_gapped_lambda: f64,
    gapped_params: &crate::stats::KarlinParams,
    x_drop_gapped_final: i32,
) -> i32 {
    x_drop_raw_score(
        X_DROP_GAPPED_FINAL as f64,
        &crate::stats::KarlinParams {
            lambda: scaled_gapped_lambda,
            k: gapped_params.k,
            h: gapped_params.h,
            alpha: gapped_params.alpha,
            beta: gapped_params.beta,
        },
        1.0,
    )
    .max(x_drop_gapped_final)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2438-2448
// ```c
// near_identical_cutoff =
//     (near_identical_cutoff_bits * NCBIMATH_LN2)
//         / context->sbp->kbp_gap[index]->Lambda;
// ```
fn blastp_redo_near_identical_cutoff(scaled_gapped_lambda: f64) -> f64 {
    (BLASTP_NEAR_IDENTICAL_BITS_PER_POSITION * std::f64::consts::LN_2) / scaled_gapped_lambda
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2453-2462
// ```c
// if (do_link_hsps) {
//     cutoff_s =
//         (int) (hitParams->cutoff_score_min * context->localScalingFactor);
// } else {
//     cutoff_s = 1;
// }
// ```
fn blastp_redo_cutoff_score(
    do_link_hsps: bool,
    cutoff_score: i32,
    local_scaling_factor: f64,
) -> i32 {
    if do_link_hsps {
        ((cutoff_score as f64) * local_scaling_factor) as i32
    } else {
        1
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:825-886
// ```c
// switch (comp_stat_string[0]) {
// case '0': ... eNoCompositionBasedStats;
// case '1': ... eCompositionBasedStats;
// case '2': ... eCompositionMatrixAdjust;
// case '3': ... eCompoForceFullMatrixAdjust;
// }
// ```
fn blastp_compo_adjust_mode(mode: BlastpCompositionMode) -> BlastCompoAdjustMode {
    match mode {
        BlastpCompositionMode::NoCompositionBasedStats => {
            BlastCompoAdjustMode::NoCompositionBasedStats
        }
        BlastpCompositionMode::CompositionBasedStats => BlastCompoAdjustMode::CompositionBasedStats,
        BlastpCompositionMode::CompositionMatrixAdjust => {
            BlastCompoAdjustMode::CompositionMatrixAdjust
        }
        BlastpCompositionMode::ForceFullMatrixAdjust => BlastCompoAdjustMode::ForceFullMatrixAdjust,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3099-3107
// ```c
// if (compo_adjust_mode != eNoCompositionBasedStats) {
//     if((0 == strcmp(scoringParams->options->matrix, "BLOSUM62_20"))) {
//         localScalingFactor = SCALING_FACTOR / 10;
//     } else {
//         localScalingFactor = SCALING_FACTOR;
//     }
// } else {
//     localScalingFactor = 1.0;
// }
// ```
fn blastp_local_scaling_factor(mode: BlastpCompositionMode, matrix: ScoringMatrix) -> f64 {
    if mode == BlastpCompositionMode::NoCompositionBasedStats {
        return 1.0;
    }
    let _ = matrix;
    BLASTP_KAPPA_SCALING_FACTOR
}
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_options.c:1163-1169
// ```c
// if ((program_number == eBlastTypeTblastn ||
//      program_number == eBlastTypeBlastp ||
//      program_number == eBlastTypeBlastx) &&
//     word_size > 5)
//     options->lut_type = eCompressedAaLookupTable;
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3036-3062
// ```c
// ECompoAdjustModes compo_adjust_mode =
//     (ECompoAdjustModes) extendParams->options->compositionBasedStats;
// if (program_number != eBlastTypeBlastp ||
//     compo_adjust_mode == eNoCompositionBasedStats) {
//     ...
// }
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319
// ```c
// s_BlastProtGappedAlignment(..., const BlastScoringParameters* score_params,
//                            BlastGapAlignStruct* gap_align, ...)
// {
//     ...
//     score_left = Blast_SemiGappedAlign(..., score_params, ...);
//     ...
//     score_right = Blast_SemiGappedAlign(..., score_params, ...);
// }
// ```
fn validate_requested_blastp_support(args: &ResolvedBlastpArgs) -> Result<()> {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:917-918
    // ```c
    // arg_desc.AddFlag(kArgUngapped, "Perform ungapped alignment only?", true);
    // ```
    if args.ungapped {
        bail!(
            "unsupported pure-Rust blastp --ungapped: ungapped-only blastp search is not yet ported"
        );
    }
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:113
    // ```c
    // const string kArgUseSWTraceback(\"use_sw_tback\");
    // ```
    if args.use_sw_tback {
        bail!(
            "unsupported pure-Rust blastp -use_sw_tback: Smith-Waterman traceback is not yet ported"
        );
    }
    if args.word_size != 3 {
        bail!(
            "unsupported pure-Rust blastp word_size={}: compressed amino-acid lookup for word sizes 5-7 is not yet ported; only --word-size 3 is currently supported",
            args.word_size
        );
    }
    if args.scoring.matrix != ScoringMatrix::Blosum62
        || args.scoring.gap_open != 11
        || args.scoring.gap_extend != 1
    {
        bail!(
            "unsupported pure-Rust blastp scoring {:?} with gap penalties {}/{}: matrix-aware protein gap alignment is not yet ported; only BLOSUM62 with gap open 11 and gap extend 1 is currently supported",
            args.scoring.matrix,
            args.scoring.gap_open,
            args.scoring.gap_extend
        );
    }
    if args.comp_based_stats.mode != BlastpCompositionMode::CompositionMatrixAdjust
        || args.comp_based_stats.unified_p
    {
        bail!(
            "unsupported pure-Rust blastp composition mode {:?}{}: only --comp_based_stats 2 is currently supported",
            args.comp_based_stats.mode,
            if args.comp_based_stats.unified_p {
                " with unified P-values"
            } else {
                ""
            }
        );
    }
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_options.c:908-936
    // ```c
    // if ((status=Blast_KarlinBlkGappedLoadFromTables(NULL, options->gap_open,
    //       options->gap_extend, options->matrix, std_matrix_only)) != 0)
    // {
    //     ...
    //     return BLASTERR_OPTION_VALUE_INVALID;
    // }
    // ```
    if !protein_scoring_supported(&args.scoring) {
        bail!(
            "blastp matrix {:?} gap penalties {}/{} are not present in the NCBI protein statistics tables used by LOSAT",
            args.scoring.matrix,
            args.scoring.gap_open,
            args.scoring.gap_extend
        );
    }
    Ok(())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:39-45
// ```c
// const char* kDfltArgTabularOutputFmt =
//     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
//     "evalue bitscore";
// const char* kDfltArgTabularOutputFmtTag("std");
// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BlastpTabularField {
    QuerySeqId,
    QueryAccessionVersion,
    QueryLength,
    SubjectSeqId,
    SubjectAccessionVersion,
    SubjectLength,
    QueryStart,
    QueryEnd,
    SubjectStart,
    SubjectEnd,
    QuerySeq,
    SubjectSeq,
    Evalue,
    BitScore,
    Score,
    AlignmentLength,
    PercentIdentical,
    NumIdentical,
    Mismatches,
    Positives,
    GapOpenings,
    Gaps,
    PercentPositives,
    QueryFrame,
    SubjectFrame,
    Frames,
    Btop,
    SubjectTitle,
}

impl BlastpTabularField {
    #[inline]
    fn name(self) -> &'static str {
        match self {
            Self::QuerySeqId => "qseqid",
            Self::QueryAccessionVersion => "qaccver",
            Self::QueryLength => "qlen",
            Self::SubjectSeqId => "sseqid",
            Self::SubjectAccessionVersion => "saccver",
            Self::SubjectLength => "slen",
            Self::QueryStart => "qstart",
            Self::QueryEnd => "qend",
            Self::SubjectStart => "sstart",
            Self::SubjectEnd => "send",
            Self::QuerySeq => "qseq",
            Self::SubjectSeq => "sseq",
            Self::Evalue => "evalue",
            Self::BitScore => "bitscore",
            Self::Score => "score",
            Self::AlignmentLength => "length",
            Self::PercentIdentical => "pident",
            Self::NumIdentical => "nident",
            Self::Mismatches => "mismatch",
            Self::Positives => "positive",
            Self::GapOpenings => "gapopen",
            Self::Gaps => "gaps",
            Self::PercentPositives => "ppos",
            Self::QueryFrame => "qframe",
            Self::SubjectFrame => "sframe",
            Self::Frames => "frames",
            Self::Btop => "btop",
            Self::SubjectTitle => "stitle",
        }
    }
}

const DEFAULT_BLASTP_TABULAR_FIELDS: [BlastpTabularField; 12] = [
    BlastpTabularField::QueryAccessionVersion,
    BlastpTabularField::SubjectAccessionVersion,
    BlastpTabularField::PercentIdentical,
    BlastpTabularField::AlignmentLength,
    BlastpTabularField::Mismatches,
    BlastpTabularField::GapOpenings,
    BlastpTabularField::QueryStart,
    BlastpTabularField::QueryEnd,
    BlastpTabularField::SubjectStart,
    BlastpTabularField::SubjectEnd,
    BlastpTabularField::Evalue,
    BlastpTabularField::BitScore,
];

#[inline]
fn default_blastp_tabular_fields() -> &'static [BlastpTabularField] {
    &DEFAULT_BLASTP_TABULAR_FIELDS
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:39-145
// ```c
// const char* kDfltArgTabularOutputFmtTag("std");
// const SFormatSpec sc_FormatSpecifiers[] = {
//     SFormatSpec("qseqid", ...),
//     SFormatSpec("qaccver", ...),
//     ...
//     SFormatSpec("btop", ...),
//     SFormatSpec("stitle", ...)
// };
// ```
fn parse_blastp_tabular_fields(spec: &str) -> Result<Vec<BlastpTabularField>> {
    let mut fields = Vec::new();

    for token in spec.split_whitespace() {
        let expanded = if token.eq_ignore_ascii_case("std") {
            DEFAULT_BLASTP_TABULAR_FIELDS.to_vec()
        } else {
            vec![match token {
                "qseqid" => BlastpTabularField::QuerySeqId,
                "qacc" | "qaccver" => BlastpTabularField::QueryAccessionVersion,
                "qlen" => BlastpTabularField::QueryLength,
                "sseqid" => BlastpTabularField::SubjectSeqId,
                "sacc" | "saccver" => BlastpTabularField::SubjectAccessionVersion,
                "slen" => BlastpTabularField::SubjectLength,
                "qstart" => BlastpTabularField::QueryStart,
                "qend" => BlastpTabularField::QueryEnd,
                "sstart" => BlastpTabularField::SubjectStart,
                "send" => BlastpTabularField::SubjectEnd,
                "qseq" => BlastpTabularField::QuerySeq,
                "sseq" => BlastpTabularField::SubjectSeq,
                "evalue" => BlastpTabularField::Evalue,
                "bitscore" => BlastpTabularField::BitScore,
                "score" => BlastpTabularField::Score,
                "length" => BlastpTabularField::AlignmentLength,
                "pident" => BlastpTabularField::PercentIdentical,
                "nident" => BlastpTabularField::NumIdentical,
                "mismatch" => BlastpTabularField::Mismatches,
                "positive" => BlastpTabularField::Positives,
                "gapopen" => BlastpTabularField::GapOpenings,
                "gaps" => BlastpTabularField::Gaps,
                "ppos" => BlastpTabularField::PercentPositives,
                "qframe" => BlastpTabularField::QueryFrame,
                "sframe" => BlastpTabularField::SubjectFrame,
                "frames" => BlastpTabularField::Frames,
                "btop" => BlastpTabularField::Btop,
                "stitle" => BlastpTabularField::SubjectTitle,
                other => bail!("unsupported blastp tabular field '{other}'"),
            }]
        };

        for field in expanded {
            if !fields.contains(&field) {
                fields.push(field);
            }
        }
    }

    if fields.is_empty() {
        bail!("blastp custom tabular field list cannot be empty");
    }

    Ok(fields)
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:985-1027
// ```c
// if (m_QuerySeq[i] == m_SubjectSeq[i]) {
//     ++num_ident;
//     ++num_positives;
//     ++num_matches;
// } else {
//     if (num_matches > 0) {
//         btop_string += NStr::Int8ToString(num_matches);
//         num_matches = 0;
//     }
//     btop_string += m_QuerySeq[i];
//     btop_string += m_SubjectSeq[i];
// }
// ```
fn build_btop(query_seq: &str, subject_seq: &str) -> String {
    let mut btop = String::new();
    let mut matches = 0usize;

    for (q, s) in query_seq.chars().zip(subject_seq.chars()) {
        if q == s {
            matches += 1;
        } else {
            if matches > 0 {
                btop.push_str(&matches.to_string());
                matches = 0;
            }
            btop.push(q);
            btop.push(s);
        }
    }

    if matches > 0 {
        btop.push_str(&matches.to_string());
    }

    btop
}

fn blastp_frame_value(frame: Option<i8>) -> i8 {
    frame.unwrap_or(1)
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1264-1283
// ```c
// if (num_hits != 0) {
//     PrintFieldNames(is_csv);
// }
// m_Ostream << "# " << num_hits << " hits found" << "\n";
// ```
fn write_blastp_outfmt7_header<W: Write>(
    writer: &mut W,
    context: &ReportContext,
    fields: &[BlastpTabularField],
    num_hits: usize,
) -> std::io::Result<()> {
    let version_str = context.version.as_deref().unwrap_or("0.1.0");
    writeln!(
        writer,
        "# {} {}",
        context.program.to_uppercase(),
        version_str
    )?;
    if let Some(ref query) = context.query_name {
        writeln!(writer, "# Query: {}", query)?;
    }
    if let Some(ref db) = context.subject_name {
        writeln!(writer, "# Database: {}", db)?;
    }
    if num_hits > 0 {
        let field_list = fields
            .iter()
            .map(|field| field.name())
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(writer, "# Fields: {}", field_list)?;
    }
    writeln!(writer, "# {} hits found", num_hits)?;
    Ok(())
}

fn blastp_tabular_field_text(
    field: BlastpTabularField,
    hit: &PairwiseHit,
    query_id: &str,
    subject_id: &str,
) -> String {
    let h = &hit.hit;
    let positives = hit.positives.unwrap_or(h.num_positives);
    let gaps = hit.gaps.unwrap_or_else(|| h.gap_letters());
    let query_frame = blastp_frame_value(hit.query_frame);
    let subject_frame = blastp_frame_value(hit.subject_frame);

    match field {
        BlastpTabularField::QuerySeqId | BlastpTabularField::QueryAccessionVersion => {
            query_id.to_string()
        }
        BlastpTabularField::QueryLength => h.query_length.to_string(),
        BlastpTabularField::SubjectSeqId | BlastpTabularField::SubjectAccessionVersion => {
            subject_id.to_string()
        }
        BlastpTabularField::SubjectLength => hit.subject_length.unwrap_or_default().to_string(),
        BlastpTabularField::QueryStart => h.q_start.to_string(),
        BlastpTabularField::QueryEnd => h.q_end.to_string(),
        BlastpTabularField::SubjectStart => h.s_start.to_string(),
        BlastpTabularField::SubjectEnd => h.s_end.to_string(),
        BlastpTabularField::QuerySeq => hit.query_seq.clone().unwrap_or_default(),
        BlastpTabularField::SubjectSeq => hit.subject_seq.clone().unwrap_or_default(),
        BlastpTabularField::Evalue => format_evalue_ncbi_tabular(h.e_value),
        BlastpTabularField::BitScore => format_bitscore_ncbi(h.bit_score),
        BlastpTabularField::Score => h.raw_score.to_string(),
        BlastpTabularField::AlignmentLength => h.length.to_string(),
        BlastpTabularField::PercentIdentical => format!("{:.3}", h.identity),
        BlastpTabularField::NumIdentical => h.num_ident.to_string(),
        BlastpTabularField::Mismatches => h.mismatch.to_string(),
        BlastpTabularField::Positives => positives.to_string(),
        BlastpTabularField::GapOpenings => h.gapopen.to_string(),
        BlastpTabularField::Gaps => gaps.to_string(),
        BlastpTabularField::PercentPositives => {
            let ppos = if h.length > 0 {
                (positives as f64 / h.length as f64) * 100.0
            } else {
                0.0
            };
            format!("{:.2}", ppos)
        }
        BlastpTabularField::QueryFrame => query_frame.to_string(),
        BlastpTabularField::SubjectFrame => subject_frame.to_string(),
        BlastpTabularField::Frames => format!("{query_frame}/{subject_frame}"),
        BlastpTabularField::Btop => {
            let query_seq = hit.query_seq.as_deref().unwrap_or("");
            let subject_seq = hit.subject_seq.as_deref().unwrap_or("");
            build_btop(query_seq, subject_seq)
        }
        BlastpTabularField::SubjectTitle => hit.subject_title.clone().unwrap_or_default(),
    }
}

fn open_output_writer(out_path: Option<&PathBuf>) -> Result<Box<dyn Write>> {
    let stdout = std::io::stdout();
    Ok(if let Some(path) = out_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout))
    })
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// void CBlastTabularInfo::Print()
// {
//     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//         if (iter != m_FieldsToShow.begin())
//             m_Ostream << m_FieldDelimiter;
//         x_PrintField(*iter);
//     }
//     m_Ostream << "\n";
// }
// ```
fn write_blastp_tabular_output(
    hits: &[PairwiseHit],
    fields: &[BlastpTabularField],
    outfmt: OutputFormat,
    out_path: Option<&PathBuf>,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    context: &ReportContext,
) -> Result<()> {
    let mut writer = open_output_writer(out_path)?;

    if outfmt == OutputFormat::TabularWithComments {
        write_blastp_outfmt7_header(&mut writer, context, fields, hits.len())?;
    }

    for hit in hits {
        let (query_id, subject_id) = hit.hit.resolve_ids(query_ids, subject_ids);
        for (index, field) in fields.iter().enumerate() {
            if index > 0 {
                writer.write_all(b"\t")?;
            }
            let value = blastp_tabular_field_text(*field, hit, query_id, subject_id);
            writer.write_all(value.as_bytes())?;
        }
        writer.write_all(b"\n")?;
    }

    Ok(())
}

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:2122-2149
// ```c
// if(sequence_standard[i]==sequence[i]){
//     match ++;
// } else {
//     if ((m_AlignType&eProt)
//         && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){
//         positive ++;
//     }
// }
// ```
fn count_identities_and_positives(
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    len: usize,
    matrix: ScoringMatrix,
) -> Option<(usize, usize)> {
    if q_off + len > query.len() || s_off + len > subject.len() {
        return None;
    }

    let mut identities = 0usize;
    let mut positives = 0usize;

    for i in 0..len {
        let q = query[q_off + i];
        let s = subject[s_off + i];
        if q == s {
            identities += 1;
            positives += 1;
        } else if protein_score(matrix, q, s) > 0 {
            positives += 1;
        }
    }

    Some((identities, positives))
}

fn decode_alignment(seq: &[u8], start: usize, end: usize) -> Option<String> {
    if start > end || end > seq.len() {
        return None;
    }
    Some(
        seq[start..end]
            .iter()
            .map(|&aa| ncbistdaa_to_ascii(aa))
            .collect(),
    )
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:518-523
// ```c
// query = query_blk->sequence + query_info->contexts[hsp->context].query_offset;
// query_nomask = query_blk->sequence_nomask +
//     query_info->contexts[hsp->context].query_offset;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:649-653
// ```c
// query = query_blk->sequence + query_info->contexts[hsp->context].query_offset;
// query_nomask = query_blk->sequence_nomask +
//     query_info->contexts[hsp->context].query_offset;
// ```
#[inline]
fn query_nomask_sequence(ctx: &QueryContext) -> &[u8] {
    let query_nomask_full = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
    &query_nomask_full[1..query_nomask_full.len() - 1]
}

fn build_pairwise_hits(
    hits: &[Hit],
    query_contexts: &[QueryContext],
    subjects: &[EncodedProtein],
    subject_titles: &[Option<String>],
    matrix: ScoringMatrix,
) -> Vec<PairwiseHit> {
    hits.iter()
        .cloned()
        .map(|hit| {
            let q_idx = hit.q_idx as usize;
            let s_idx = hit.s_idx as usize;
            let query = query_nomask_sequence(&query_contexts[q_idx]);
            let subject = &subjects[s_idx].aa_seq[1..subjects[s_idx].aa_seq.len() - 1];

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
            // ```c
            // s_Blast_HSPGetNumIdentitiesAndPositives(...)
            // {
            //    if (!hsp->gap_info) { ... } else { ... }
            // }
            // ```
            let rendered = build_alignment_view_with_matrix(query, subject, &hit, matrix);
            let (query_seq, subject_seq, positives, gaps) = if let Some(view) = rendered {
                (
                    Some(view.query_aln),
                    Some(view.subject_aln),
                    Some(view.num_positives),
                    Some(view.gap_letters),
                )
            } else {
                let q_start = hit.q_start.saturating_sub(1);
                let q_end = hit.q_end.min(query.len());
                let s_start = hit.s_start.saturating_sub(1);
                let s_end = hit.s_end.min(subject.len());
                let ungapped_span = q_end
                    .saturating_sub(q_start)
                    .min(s_end.saturating_sub(s_start));
                if ungapped_span == hit.length {
                    (
                        decode_alignment(query, q_start, q_start + ungapped_span),
                        decode_alignment(subject, s_start, s_start + ungapped_span),
                        Some(hit.num_positives),
                        Some(hit.gap_letters()),
                    )
                } else {
                    (None, None, Some(hit.num_positives), Some(hit.gap_letters()))
                }
            };

            PairwiseHit {
                hit,
                query_seq,
                subject_seq,
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1049-1093
                // ```c
                // int query_frame = 1, subject_frame = 1;
                // ...
                // SetCounts(num_ident, align_length, num_gaps, num_gap_opens,
                //           num_positives, query_frame, subject_frame);
                // ```
                query_frame: Some(1),
                subject_frame: Some(1),
                positives,
                gaps,
                subject_length: Some(subjects[s_idx].aa_len),
                subject_title: subject_titles[s_idx].clone(),
            }
        })
        .collect()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
// ```c
// diag_array_length = 1;
// while (diag_array_length < (qlen+window_size))
//     diag_array_length = diag_array_length << 1;
// diag_table->diag_mask = diag_array_length-1;
// ```
fn compute_diag_table_shape(contexts: &[QueryContext], window: i32) -> (i32, i32) {
    let query_length = contexts
        .last()
        .map(|context| context.frame_base + context.aa_len as i32)
        .unwrap_or(0);

    let mut diag_array_size = 1i32;
    while diag_array_size < query_length + window {
        diag_array_size <<= 1;
    }

    (diag_array_size, diag_array_size - 1)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3748-3786
// ```c
// /* find the first ungapped alignment for each query and determine
//    whetherto use approximate gapped alignment for the query */
// for (i = 0;i < init_hitlist->total;i++) {
//     ...
//     if (found[query_index]) {
//         continue;
//     }
//     found[query_index] = TRUE;
//     if (init_hsp->ungapped_data && init_hsp->ungapped_data->score <
//         (Int4)(kRestrictedMult * hit_params->cutoffs[contxt].cutoff_score)) {
//         restricted_align_array[query_index] = TRUE;
//     } else {
//         restricted_align_array[query_index] = FALSE;
//     }
// }
// ```
fn build_restricted_align_array(
    ungapped_hits: &[UngappedHit],
    cutoff_scores: &[i32],
    num_queries: usize,
) -> Vec<bool> {
    let mut found = vec![false; num_queries];
    let mut restricted = vec![false; num_queries];
    for hit in ungapped_hits {
        let query_index = hit.q_idx as usize;
        if query_index >= num_queries || found[query_index] {
            continue;
        }
        found[query_index] = true;
        restricted[query_index] =
            hit.raw_score < (BLASTP_RESTRICTED_MULT * cutoff_scores[hit.ctx_idx] as f64) as i32;
    }
    restricted
}

fn preliminary_tree_hsp(
    preliminary_hsp: &BlastpPreliminaryHsp,
    query_context_offset: i32,
) -> TreeHsp {
    TreeHsp {
        query_offset: preliminary_hsp.query_start as i32,
        query_end: preliminary_hsp.query_end as i32,
        subject_offset: preliminary_hsp.subject_start as i32,
        subject_end: preliminary_hsp.subject_end as i32,
        score: preliminary_hsp.raw_score,
        query_frame: 1,
        query_length: preliminary_hsp.query_length as i32,
        query_context_offset,
        subject_frame_sign: 1,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
// ```c
// int ScoreCompareHSPs(const void* h1, const void* h2) {
//    ...
// }
// ```
fn preliminary_score_compare(
    a: &BlastpPreliminaryHsp,
    b: &BlastpPreliminaryHsp,
) -> std::cmp::Ordering {
    b.raw_score
        .cmp(&a.raw_score)
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| b.query_end.cmp(&a.query_end))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
// ```c
// Int4
// Blast_HSPListPurgeHSPsWithCommonEndpoints(EBlastProgramType program,
//                                           BlastHSPList* hsp_list,
//                                           Boolean purge)
// {
//    ...
// }
// ```
fn purge_preliminary_hits_with_common_endpoints(
    mut hits: Vec<BlastpPreliminaryHsp>,
) -> Vec<BlastpPreliminaryHsp> {
    if hits.len() <= 1 {
        return hits;
    }

    hits.sort_by(|a, b| {
        (a.q_idx, a.query_frame)
            .cmp(&(b.q_idx, b.query_frame))
            .then_with(|| a.query_start.cmp(&b.query_start))
            .then_with(|| a.subject_start.cmp(&b.subject_start))
            .then_with(|| b.raw_score.cmp(&a.raw_score))
            .then_with(|| b.query_end.cmp(&a.query_end))
            .then_with(|| b.subject_end.cmp(&a.subject_end))
    });

    let mut i = 0usize;
    while i + 1 < hits.len() {
        let j = i + 1;
        let same_context =
            hits[i].q_idx == hits[j].q_idx && hits[i].query_frame == hits[j].query_frame;
        let same_q_start = hits[i].query_start == hits[j].query_start;
        let same_s_start = hits[i].subject_start == hits[j].subject_start;
        if same_context && same_q_start && same_s_start {
            hits.remove(j);
        } else {
            i += 1;
        }
    }

    hits.sort_by(|a, b| {
        (a.q_idx, a.query_frame)
            .cmp(&(b.q_idx, b.query_frame))
            .then_with(|| a.query_end.cmp(&b.query_end))
            .then_with(|| a.subject_end.cmp(&b.subject_end))
            .then_with(|| b.raw_score.cmp(&a.raw_score))
            .then_with(|| b.query_start.cmp(&a.query_start))
            .then_with(|| b.subject_start.cmp(&a.subject_start))
    });

    let mut i = 0usize;
    while i + 1 < hits.len() {
        let j = i + 1;
        let same_context =
            hits[i].q_idx == hits[j].q_idx && hits[i].query_frame == hits[j].query_frame;
        let same_q_end = hits[i].query_end == hits[j].query_end;
        let same_s_end = hits[i].subject_end == hits[j].subject_end;
        if same_context && same_q_end && same_s_end {
            hits.remove(j);
        } else {
            i += 1;
        }
    }

    hits
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:540-555
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
// ...
// Blast_HSPListSortByScore(hsp_list);
// ```
fn prepare_preliminary_hits_for_kappa(
    mut hits: Vec<BlastpPreliminaryHsp>,
) -> Vec<BlastpPreliminaryHsp> {
    hits = purge_preliminary_hits_with_common_endpoints(hits);
    hits.sort_by(preliminary_score_compare);
    hits
}

fn rebuild_preliminary_interval_tree(
    interval_tree: &mut BlastIntervalTree,
    preliminary_hits_by_query: &HashMap<u32, Vec<BlastpPreliminaryHsp>>,
    contexts: &[QueryContext],
) {
    interval_tree.reset();
    for preliminary_hits in preliminary_hits_by_query.values() {
        for preliminary_hsp in preliminary_hits {
            interval_tree.add_hsp(
                preliminary_tree_hsp(
                    preliminary_hsp,
                    contexts[preliminary_hsp.q_idx as usize].frame_base,
                ),
                contexts[preliminary_hsp.q_idx as usize].frame_base,
                IndexMethod::QueryAndSubject,
            );
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3777-3786
// ```c
// /* use approximate gapped alignment if the highest scoring ungapped
//    alignment scores below the reduced cutoff */
// if (init_hsp->ungapped_data && init_hsp->ungapped_data->score <
//     (Int4)(kRestrictedMult * hit_params->cutoffs[contxt].cutoff_score)) {
//     restricted_align_array[query_index] = TRUE;
// } else {
//     restricted_align_array[query_index] = FALSE;
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:635-688
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
// ...
// Blast_HSPListSortByScore(hsp_list);
// ...
// if (BlastIntervalTreeContainsHSP(tree, hsp, query_info,
//                                  hit_options->min_diag_separation)) {
//     hsp_array[index] = Blast_HSPFree(hsp);
// }
// ```
fn extend_preliminary_blastp_hit(
    hit: &UngappedHit,
    contexts: &[QueryContext],
    subject_raw: &[u8],
    cutoff_score: i32,
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    mode: BlastpGappedAlignmentMode,
) -> Option<BlastpPreliminaryHsp> {
    let ctx = &contexts[hit.ctx_idx];
    let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
    let ungapped_len = hit.q_aa_end.saturating_sub(hit.q_aa_start);
    if ungapped_len == 0 {
        return None;
    }
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3937-3945
    // ```c
    // max_offset =
    //    BlastGetStartForGappedAlignment(query_tmp.sequence,
    //       subject->sequence, gap_align->sbp,
    //       init_hsp->ungapped_data->q_start,
    //       init_hsp->ungapped_data->length,
    //       init_hsp->ungapped_data->s_start,
    //       init_hsp->ungapped_data->length);
    // init_hsp->offsets.qs_offsets.s_off +=
    //     max_offset - init_hsp->offsets.qs_offsets.q_off;
    // init_hsp->offsets.qs_offsets.q_off = max_offset;
    // ```
    let gapped_q_start = blastp_get_start_for_gapped_alignment(
        query_raw,
        subject_raw,
        hit.q_aa_start,
        ungapped_len,
        hit.s_aa_start,
        hit.s_aa_end.saturating_sub(hit.s_aa_start),
        matrix,
    );
    let gapped_s_start = hit
        .s_aa_start
        .saturating_add(gapped_q_start.saturating_sub(hit.q_aa_start));
    if gapped_q_start >= query_raw.len() || gapped_s_start >= subject_raw.len() {
        return None;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3949-3978
    // ```c
    // status = s_BlastProtGappedAlignment(program_number, &query_tmp,
    //                                     subject, gap_align,
    //                                     score_params, init_hsp,
    //                                     restricted_alignment,
    //                                     fence_hit);
    // if (restricted_alignment &&
    //     gap_align->score < cutoff &&
    //     gap_align->score >= restricted_cutoff) {
    //     ...
    // }
    // ```
    let aligned = blastp_score_only_gapped_alignment(
        query_raw,
        subject_raw,
        gapped_q_start,
        gapped_s_start,
        matrix,
        gap_open,
        gap_extend,
        x_dropoff,
        mode,
    )?;

    if aligned.score < cutoff_score {
        return None;
    }

    if aligned.query_stop < aligned.query_start || aligned.subject_stop < aligned.subject_start {
        return None;
    }

    Some(BlastpPreliminaryHsp {
        query_start: aligned.query_start,
        query_end: aligned.query_stop,
        subject_start: aligned.subject_start,
        subject_end: aligned.subject_stop,
        gapped_query_start: gapped_q_start,
        gapped_subject_start: gapped_s_start,
        raw_score: aligned.score,
        query_frame: 0,
        query_length: ctx.aa_len,
        q_idx: hit.q_idx,
        s_idx: hit.s_idx,
    })
}

pub fn run(args: BlastpArgs) -> Result<()> {
    let args = args.resolve()?;
    validate_requested_blastp_support(&args)?;

    let (outfmt, custom_fields) = OutputFormat::parse(&args.outfmt).map_err(anyhow::Error::msg)?;
    if outfmt == OutputFormat::Pairwise && custom_fields.is_some() {
        bail!("blastp outfmt 0 does not accept custom field lists");
    }
    let custom_fields = custom_fields
        .as_deref()
        .map(parse_blastp_tabular_fields)
        .transpose()?;

    let query_records: Vec<fasta::Record> = fasta::Reader::from_file(&args.query)
        .with_context(|| format!("failed to open query FASTA {}", args.query.display()))?
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("failed to read query FASTA {}", args.query.display()))?;

    let subject_records: Vec<fasta::Record> = fasta::Reader::from_file(&args.subject)
        .with_context(|| format!("failed to open subject FASTA {}", args.subject.display()))?
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("failed to read subject FASTA {}", args.subject.display()))?;

    let query_ids: Vec<Arc<str>> = query_records.iter().map(fasta_id).collect();
    let subject_ids: Vec<Arc<str>> = subject_records.iter().map(fasta_id).collect();
    let subject_titles: Vec<Option<String>> = subject_records
        .iter()
        .map(|record| record.desc().map(str::to_string))
        .collect();

    let seg_params = args.seg.params();
    let query_frames: Vec<Vec<QueryFrame>> = query_records
        .iter()
        .map(|record| {
            vec![encode_protein_query_frame_with_seg(
                record.seq(),
                seg_params.as_ref(),
            )]
        })
        .collect();
    let subjects: Vec<EncodedProtein> = subject_records
        .iter()
        .map(|record| encode_protein_sequence(record.seq()))
        .collect();

    let ungapped_params = lookup_protein_params_ungapped(args.scoring.matrix);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2833-2852
    // ```c
    // Blast_ScoreBlkKbpIdealCalc(BlastScoreBlk* sbp)
    // {
    //    Blast_ResFreqStdComp(sbp, stdrfp);
    //    BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
    //    Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);
    // }
    // ```
    let ideal_ungapped_params = compute_blosum62_ideal_karlin_params()
        .map_err(anyhow::Error::msg)
        .context("failed to calculate NCBI blastp kbp_ideal parameters")?;
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1833-1839
    // ```c
    // kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    // hsp->evalue = BLAST_KarlinStoE_simple(hsp->score, kbp[hsp->context],
    //                                       hsp_list->hspcnt ? searchsp : 1);
    // ```
    let gapped_params = lookup_protein_params(&args.scoring);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:457-463
    // ```c
    // params->gap_x_dropoff = (Int4)
    //     (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
    // params->gap_x_dropoff_final = (Int4)
    //     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
    //         params->gap_x_dropoff);
    // ```
    let x_drop_gapped = x_drop_raw_score(X_DROP_GAPPED_PRELIM as f64, &gapped_params, 1.0);
    let x_drop_gapped_final =
        x_drop_raw_score(X_DROP_GAPPED_FINAL as f64, &gapped_params, 1.0).max(x_drop_gapped);
    let (lookup, contexts) =
        build_ncbi_lookup(&query_frames, args.threshold as i32, &ungapped_params);

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_setup.c:766-769
    // ```c
    // kbp_ptr = (scoring_options->gapped_calculation ? sbp->kbp_gap_std : sbp->kbp);
    // ```
    let total_db_len = subjects.iter().map(|subject| subject.aa_len).sum::<usize>();
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_setup.c:914-921
    // ```c
    // /* Set the database length for new FSC */
    // if (sbp->gbp) {
    //     sbp->gbp->db_length =
    //         (Blast_SubjectIsTranslated(program_number)) ? dbl/3 : dbl;
    // }
    // ```
    let gapped_gumbel = lookup_protein_gumbel_params(&args.scoring, total_db_len as i64)
        .context("missing NCBI gumbel parameters for supported blastp scoring")?;
    let num_subjects = subjects.len().max(1);
    let search_spaces: Vec<SearchSpace> = contexts
        .iter()
        .map(|ctx| {
            SearchSpace::for_database_search(
                ctx.aa_len,
                total_db_len,
                num_subjects,
                &gapped_params,
                true,
            )
        })
        .collect();
    let x_dropoff_per_context: Vec<i32> = contexts
        .iter()
        .map(|ctx| x_drop_raw_score(X_DROP_UNGAPPED_BITS, &ctx.karlin_params, 1.0))
        .collect();
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:336-373
    // ```c
    // if (sbp->kbp_std) {
    //    gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 +
    //                             kbp->logK) / kbp->Lambda);
    // }
    // ...
    // if (!gapped_calculation || sbp->matrix_only_scoring) {
    //    BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp, searchsp, TRUE, gap_decay_rate);
    //    new_cutoff = MIN(new_cutoff, gap_trigger);
    // } else {
    //    new_cutoff = gap_trigger;
    // }
    // new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
    // ```
    let cutoff_scores: Vec<i32> = contexts
        .iter()
        .zip(search_spaces.iter())
        .map(|(ctx, search_space)| {
            let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, &ctx.karlin_params);
            let cutoff_score_max = cutoff_score_from_evalue(
                args.evalue,
                search_space.effective_space.max(1.0) as i64,
                &gapped_params,
            );
            gap_trigger.min(cutoff_score_max)
        })
        .collect();

    let prelim_hitlist_size = get_prelim_hitlist_size(
        args.max_target_seqs,
        args.comp_based_stats.is_enabled(),
        true,
    );
    let mut hit_lists: Vec<Option<BlastpHitList>> = std::iter::repeat_with(|| None)
        .take(query_ids.len())
        .collect();

    let window = args.window_size as i32;
    let word_size = args.word_size as i32;
    let (diag_array_size, diag_mask) = compute_diag_table_shape(&contexts, window);
    let total_query_span = contexts
        .last()
        .map(|context| context.frame_base + context.aa_len as i32 + 1)
        .unwrap_or(1);

    for (s_idx, subject) in subjects.iter().enumerate() {
        let subject_frame = QueryFrame {
            frame: 0,
            aa_seq: subject.aa_seq.clone(),
            aa_seq_nomask: None,
            aa_len: subject.aa_len,
            orig_len: subject.aa_len,
            seg_masks: Vec::new(),
        };
        let subject_frames = [subject_frame];
        let subject_raw = &subject.aa_seq[1..subject.aa_seq.len() - 1];
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3826-3837
        // ```c
        // tree = Blast_IntervalTreeInit(0, query->length+1,
        //                               0, subject->length+1);
        // ```
        let mut interval_tree =
            BlastIntervalTree::new(0, total_query_span, 0, subject.aa_len as i32 + 1);

        let mut init_hsps: Vec<InitHSP> = Vec::new();
        let mut diag_array = vec![DiagStruct::clear(window); diag_array_size as usize];
        let diag_offset = window;

        if subject_raw.len() >= args.word_size {
            for s_off in 0..=subject_raw.len() - args.word_size {
                let index = encode_kmer(
                    subject_raw[s_off] as usize,
                    subject_raw[s_off + 1] as usize,
                    subject_raw[s_off + 2] as usize,
                );

                for &query_offset_i32 in lookup.get_hits(index) {
                    let query_offset = query_offset_i32 as u32;
                    let subject_offset = s_off as u32;
                    let diag_coord =
                        (query_offset.wrapping_sub(subject_offset) & (diag_mask as u32)) as usize;
                    let diag_entry = &mut diag_array[diag_coord];
                    let ctx_idx = lookup.get_context_idx(query_offset as i32);
                    let ctx = &contexts[ctx_idx];
                    let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
                    let q_right_off = query_offset.wrapping_sub(ctx.frame_base as u32) as usize;
                    let x_drop = x_dropoff_per_context[ctx_idx];
                    let cutoff = cutoff_scores[ctx_idx];

                    if args.window_size == 0 {
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:713-775
                        // ```c
                        // diff = subject_offset -
                        //     (diag_array[diag_coord].last_hit - diag_offset);
                        // if (diff >= 0) {
                        //     score = s_BlastAaExtendOneHit(...);
                        //     if (score >= cutoffs->cutoff_score) {
                        //         BlastSaveInitHsp(...);
                        //     }
                        //     diag_array[diag_coord].last_hit =
                        //         s_last_off - (wordsize - 1) + diag_offset;
                        // }
                        // ```
                        let diff = subject_offset as i32 - (diag_entry.last_hit() - diag_offset);
                        if diff < 0 {
                            continue;
                        }

                        let Some((q_start, q_end, s_start, s_end, score, s_last_off)) =
                            extend_one_hit(
                                args.scoring.matrix,
                                query_raw,
                                subject_raw,
                                q_right_off,
                                subject_offset as usize,
                                x_drop,
                                args.word_size,
                            )
                        else {
                            continue;
                        };

                        diag_entry.set_last_hit(s_last_off as i32 - (word_size - 1) + diag_offset);

                        if score >= cutoff {
                            init_hsps.push(InitHSP {
                                q_start_absolute: ctx.frame_base + q_start as i32,
                                q_end_absolute: ctx.frame_base + q_end as i32,
                                s_start: s_start as i32,
                                s_end: s_end as i32,
                                score,
                                ctx_idx,
                                s_f_idx: 0,
                                q_idx: ctx.q_idx,
                                s_idx: s_idx as u32,
                                q_frame: 0,
                                s_frame: 0,
                                q_orig_len: ctx.aa_len,
                                s_orig_len: subject.aa_len,
                            });
                        }
                        continue;
                    }

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:488-574
                    // ```c
                    // if (diag_array[diag_coord].flag) {
                    //     ...
                    // } else {
                    //     last_hit = diag_array[diag_coord].last_hit - diag_offset;
                    //     diff = subject_offset - last_hit;
                    //     if (diff >= window) { ... }
                    //     if (diff < wordsize) continue;
                    //     ...
                    //     score = s_BlastAaExtendTwoHit(...);
                    //     if (right_extend) {
                    //         diag_array[diag_coord].flag = 1;
                    //         diag_array[diag_coord].last_hit =
                    //             s_last_off - (wordsize - 1) + diag_offset;
                    //     } else {
                    //         diag_array[diag_coord].last_hit =
                    //             subject_offset + diag_offset;
                    //     }
                    // }
                    // ```
                    if diag_entry.flag() != 0 {
                        let subject_plus_offset = subject_offset.wrapping_add(diag_offset as u32);
                        if subject_plus_offset < diag_entry.last_hit() as u32 {
                            continue;
                        }
                        diag_entry.set_last_hit(subject_plus_offset as i32);
                        diag_entry.set_flag(0);
                        continue;
                    }

                    let last_hit = diag_entry.last_hit() - diag_offset;
                    let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;
                    if diff >= window {
                        diag_entry
                            .set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32);
                        continue;
                    }
                    if diff < word_size {
                        continue;
                    }

                    let q_minus_diff = query_offset.wrapping_sub(diff as u32);
                    if q_minus_diff < ctx.frame_base as u32 {
                        diag_entry
                            .set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32);
                        continue;
                    }

                    let Some((q_start, q_end, s_start, s_end, score, right_extended, s_last_off)) =
                        extend_two_hit(
                            args.scoring.matrix,
                            query_raw,
                            subject_raw,
                            (last_hit + word_size) as usize,
                            subject_offset as usize,
                            q_right_off,
                            x_drop,
                            args.word_size,
                        )
                    else {
                        continue;
                    };

                    if right_extended {
                        diag_entry.set_flag(1);
                        diag_entry.set_last_hit(s_last_off as i32 - (word_size - 1) + diag_offset);
                    } else {
                        diag_entry
                            .set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32);
                    }

                    if score >= cutoff {
                        init_hsps.push(InitHSP {
                            q_start_absolute: ctx.frame_base + q_start as i32,
                            q_end_absolute: ctx.frame_base + q_end as i32,
                            s_start: s_start as i32,
                            s_end: s_end as i32,
                            score,
                            ctx_idx,
                            s_f_idx: 0,
                            q_idx: ctx.q_idx,
                            s_idx: s_idx as u32,
                            q_frame: 0,
                            s_frame: 0,
                            q_orig_len: ctx.aa_len,
                            s_orig_len: subject.aa_len,
                        });
                    }
                }
            }
        }
        let ungapped_hits = get_ungapped_hsp_list(init_hsps, &contexts, &subject_frames);
        let mut preliminary_hits_by_query: HashMap<u32, Vec<BlastpPreliminaryHsp>> = HashMap::new();
        let mut restricted_align_array =
            build_restricted_align_array(&ungapped_hits, &cutoff_scores, query_ids.len());
        let mut redo_index: Option<usize> = None;
        let mut redo_query: Option<u32> = None;
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1481-1497
        // ```c
        // if (!gapped_calculation) {
        //     /* The following must be performed for any ungapped
        //        search with a nucleotide database. */
        //     status = Blast_HSPListReevaluateUngapped(...);
        // }
        // ```
        let mut hit_index = 0usize;
        while hit_index < ungapped_hits.len() {
            let hit = &ungapped_hits[hit_index];
            let ctx = &contexts[hit.ctx_idx];
            let query_index = hit.q_idx as usize;
            if let (Some(redo_start), Some(redo_query_index)) = (redo_index, redo_query) {
                if hit_index < redo_start && hit.q_idx != redo_query_index {
                    hit_index += 1;
                    continue;
                }
            }
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3918
            // ```c
            // tmp_hsp.query.offset = q_start;
            // tmp_hsp.query.end = q_end;
            // tmp_hsp.subject.offset = s_start;
            // tmp_hsp.subject.end = s_end;
            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                   hit_options->min_diag_separation))
            // ```
            let ungapped_tree_hsp = TreeHsp {
                query_offset: hit.q_aa_start as i32,
                query_end: hit.q_aa_end as i32,
                subject_offset: hit.s_aa_start as i32,
                subject_end: hit.s_aa_end as i32,
                score: hit.raw_score,
                query_frame: 1,
                query_length: ctx.aa_len as i32,
                query_context_offset: ctx.frame_base,
                subject_frame_sign: 1,
            };
            if interval_tree.contains_hsp(&ungapped_tree_hsp, ctx.frame_base, 0) {
                hit_index += 1;
                continue;
            }
            let restricted_alignment = restricted_align_array
                .get(query_index)
                .copied()
                .unwrap_or(false);
            let cutoff = cutoff_scores[hit.ctx_idx];
            let restricted_cutoff = (BLASTP_RESTRICTED_MULT * cutoff as f64) as i32;
            let mode = if restricted_alignment {
                BlastpGappedAlignmentMode::Restricted
            } else {
                BlastpGappedAlignmentMode::Exact
            };
            let Some(preliminary_hsp) = extend_preliminary_blastp_hit(
                hit,
                &contexts,
                subject_raw,
                cutoff,
                args.scoring.matrix,
                args.scoring.gap_open,
                args.scoring.gap_extend,
                x_drop_gapped,
                mode,
            ) else {
                hit_index += 1;
                continue;
            };
            if restricted_alignment
                && preliminary_hsp.raw_score < cutoff
                && preliminary_hsp.raw_score >= restricted_cutoff
            {
                if let Some(restricted_slot) = restricted_align_array.get_mut(query_index) {
                    *restricted_slot = false;
                }
                preliminary_hits_by_query.remove(&hit.q_idx);
                rebuild_preliminary_interval_tree(
                    &mut interval_tree,
                    &preliminary_hits_by_query,
                    &contexts,
                );
                redo_index = Some(hit_index);
                redo_query = Some(hit.q_idx);
                hit_index = 0;
                continue;
            }
            interval_tree.add_hsp(
                preliminary_tree_hsp(&preliminary_hsp, ctx.frame_base),
                ctx.frame_base,
                IndexMethod::QueryAndSubject,
            );
            preliminary_hits_by_query
                .entry(hit.q_idx)
                .or_default()
                .push(preliminary_hsp);
            hit_index += 1;
        }

        for (q_idx, preliminary_hits) in preliminary_hits_by_query {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3658-3691
            // ```c
            // s_HSPListFromDistinctAlignments(hsp_list, ...);
            // if (hsp_list->hspcnt > 1) {
            //     s_HitlistReapContained(hsp_list->hsp_array, &hsp_list->hspcnt);
            // }
            // s_HitlistEvaluateAndPurge(...);
            // if (best_evalue <= hitParams->options->expect_value) {
            //     s_HSPListNormalizeScores(...);
            //     s_ComputeNumIdentities(...);
            // }
            // ```
            let local_scaling_factor =
                blastp_local_scaling_factor(args.comp_based_stats.mode, args.scoring.matrix);
            let scaled_gapped_lambda = gapped_params.lambda / local_scaling_factor;
            let do_link_hsps = false;
            let matrix_info = build_matrix_info(
                args.scoring.matrix,
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2216-2231
                // ```c
                // self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
                // status = s_GetStartFreqRatios(self->startFreqRatios, matrixName);
                // Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
                //                          self->startFreqRatios, self->ungappedLambda);
                // ```
                ideal_ungapped_params.lambda / local_scaling_factor,
            )?;
            let redo_align_params = BlastRedoAlignParams {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2216-2231
                // ```c
                // self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
                // status = s_GetStartFreqRatios(self->startFreqRatios, matrixName);
                // Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
                //                          self->startFreqRatios, self->ungappedLambda);
                // ```
                matrix_info,
                // NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:97-106
                // ```c
                // typedef struct BlastCompo_GappingParams {
                //     int gap_open;
                //     int gap_extend;
                //     int decline_align;
                //     int x_dropoff;
                //     void * context;
                // } BlastCompo_GappingParams;
                // ```
                //
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2117-2133
                // ```c
                // kbp->Lambda /= scale_factor;
                // sp->gap_open = (Int4)BLAST_Nint(sp->gap_open  * scale_factor);
                // sp->gap_extend = (Int4)BLAST_Nint(sp->gap_extend * scale_factor);
                // ```
                //
                gapping_params: BlastCompoGappingParams {
                    gap_open: ((args.scoring.gap_open as f64) * local_scaling_factor).round()
                        as i32,
                    gap_extend: ((args.scoring.gap_extend as f64) * local_scaling_factor).round()
                        as i32,
                    decline_align: 0,
                    x_dropoff: blastp_redo_x_dropoff(
                        scaled_gapped_lambda,
                        &gapped_params,
                        x_drop_gapped_final,
                    ),
                    context: None,
                },
                compo_adjust_mode: blastp_compo_adjust_mode(args.comp_based_stats.mode),
                // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1107-1118
                // ```c
                // Blast_RedoOneMatch(..., int ** matrix, int alphsize,
                //                    ..., int compositionTestIndex, ...)
                // ```
                alphsize: BLASTAA_SIZE as i32,
                composition_test_index: i32::from(args.comp_based_stats.unified_p),
                unified_p: args.comp_based_stats.unified_p,
                log_k: gapped_params.k.ln(),
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3584-3690
                // ```c
                // localScalingFactor      /* i */
                // ...
                // s_HSPListNormalizeScores(hsp_list, kbp->Lambda, kbp->logK,
                //                          localScalingFactor);
                // ```
                score_divisor: local_scaling_factor,
                restricted_alignment: false,
                smith_waterman: false,
                near_identical_cutoff: blastp_redo_near_identical_cutoff(scaled_gapped_lambda),
                position_based: false,
                re_matrix_adjustment_pseudocounts: 20,
                ccat_query_length: contexts[q_idx as usize].aa_len as i32,
                query_is_translated: false,
                subject_is_translated: false,
                cutoff_score: blastp_redo_cutoff_score(
                    do_link_hsps,
                    cutoff_scores[q_idx as usize],
                    local_scaling_factor,
                ),
                cutoff_evalue: args.evalue,
                do_link_hsps,
                // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:382-406
                // ```c
                // Nearly identical alignments are computed with exact subjects,
                // and others with segged subjects; this makes comparing end
                // points more difficult.
                // ```
                is_same_adjustment: false,
            };
            let ctx = &contexts[q_idx as usize];
            let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
            let query_nomask = query_nomask_sequence(ctx);
            let preliminary_hits = prepare_preliminary_hits_for_kappa(preliminary_hits);
            let diagnostics_preliminary_hits =
                diagnostics_enabled().then(|| preliminary_hits.clone());
            let first_preliminary_hit = diagnostics_preliminary_hits
                .as_ref()
                .and_then(|hits| hits.first().copied());
            let postprocessed = postprocess_preliminary_hits(
                preliminary_hits,
                query_raw,
                query_nomask,
                subject_raw,
                args.scoring.matrix,
                &gapped_params,
                &gapped_gumbel,
                ctx.aa_len as i32,
                search_spaces[q_idx as usize].length_adjustment as i32,
                search_spaces[q_idx as usize].effective_space,
                &redo_align_params,
                args.evalue,
            )?;
            let hits = postprocessed.hits;
            if diagnostics_enabled() {
                if let (Some(prelim), Some(hit)) = (first_preliminary_hit, hits.first()) {
                    let query_id = &query_ids[q_idx as usize];
                    let subject_id = &subject_ids[prelim.s_idx as usize];
                    eprintln!(
                        "[BLASTP KAPPA] query_id={} subject_id={} q_idx={} s_idx={} prelim q={}..{} s={}..{} gapped={}..{} -> redo q={}..{} s={}..{} score={} prelim_count={} hit_count={}",
                        query_id,
                        subject_id,
                        q_idx,
                        prelim.s_idx,
                        prelim.query_start,
                        prelim.query_end,
                        prelim.subject_start,
                        prelim.subject_end,
                        prelim.gapped_query_start,
                        prelim.gapped_subject_start,
                        hit.q_start,
                        hit.q_end,
                        hit.s_start,
                        hit.s_end,
                        hit.raw_score,
                        diagnostics_preliminary_hits.as_ref().map(|hits| hits.len()).unwrap_or(0),
                        hits.len(),
                    );
                }
            }
            if hits.is_empty() {
                continue;
            }

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3658-3724
            // ```c
            // s_HSPListFromDistinctAlignments(hsp_list, ...);
            // if (hsp_list->hspcnt > 1) {
            //     s_HitlistReapContained(hsp_list->hsp_array, &hsp_list->hspcnt);
            // }
            // s_HitlistEvaluateAndPurge(...);
            // if (best_evalue <= hitParams->options->expect_value) {
            //     s_HSPListNormalizeScores(hsp_list, kbp->Lambda, kbp->logK,
            //                              localScalingFactor);
            //     s_ComputeNumIdentities(...);
            //     if (BlastCompo_HeapWouldInsert(...)) {
            //         BlastCompo_HeapInsert(...);
            //     }
            // }
            // ```
            let hsp_list = BlastpHspList {
                oid: s_idx as u32,
                query_index: q_idx,
                hsps: hits.into_iter().map(BlastpHsp::from_hit).collect(),
                best_evalue: i32::MAX as f64,
            };

            let hit_list = hit_lists[q_idx as usize]
                .get_or_insert_with(|| BlastpHitList::new(prelim_hitlist_size));
            hit_list.update(hsp_list);
        }
    }

    for hit_list_opt in &mut hit_lists {
        if let Some(hit_list) = hit_list_opt.as_mut() {
            hit_list.sort_by_evalue();
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:845-852
            // ```c
            // for (subject_index = 0; subject_index < hit_list->hsplist_count; ++subject_index) {
            //    BlastHSPList * hsp_list = hit_list->hsplist_array[subject_index];
            //    if(hit_options->max_hsps_per_subject) {
            //       Blast_TrimHSPListByMaxHsps(hsp_list, hit_options);
            //    }
            // }
            // ```
            if args.max_hsps_per_subject > 0 {
                for hsp_list in hit_list
                    .hsplist_array
                    .iter_mut()
                    .take(hit_list.hsplist_count)
                {
                    trim_by_max_hsps(hsp_list, args.max_hsps_per_subject);
                }
            }
            hit_list.prune_by_size(args.max_target_seqs);
        }
    }

    let final_hits = collect_hits_from_hit_lists(&hit_lists);
    let query_header = query_records.first().map(fasta_defline);
    let context = ReportContext {
        query_name: query_header,
        query_length: query_records.first().map(|record| record.seq().len()),
        subject_name: Some(path_label(&args.subject)),
        database_name: Some(format!(
            "User specified sequence set (Input: {})",
            args.subject.display()
        )),
        database_num_sequences: Some(subject_records.len()),
        database_total_letters: Some(total_db_len),
        program: "blastp".to_string(),
        version: Some(NCBI_BLASTP_VERSION.to_string()),
    };

    let pairwise_hits = if matches!(
        outfmt,
        OutputFormat::Pairwise | OutputFormat::Tabular | OutputFormat::TabularWithComments
    ) {
        Some(build_pairwise_hits(
            &final_hits,
            &contexts,
            &subjects,
            &subject_titles,
            args.scoring.matrix,
        ))
    } else {
        None
    };

    match outfmt {
        OutputFormat::Pairwise => {
            let mut writer = open_output_writer(args.out.as_ref())?;
            let pairwise_config = PairwiseConfig {
                line_length: 60,
                show_gi: false,
                show_frame: false,
                program: "blastp".to_string(),
                protein_matrix: args.scoring.matrix,
            };
            write_pairwise(
                pairwise_hits
                    .as_ref()
                    .expect("pairwise hits prepared for blastp outfmt 0"),
                &mut writer,
                &pairwise_config,
                &query_ids,
                &subject_ids,
                &context,
            )?;
        }
        OutputFormat::Tabular | OutputFormat::TabularWithComments => {
            let fields: &[BlastpTabularField] = if let Some(fields) = custom_fields.as_ref() {
                fields.as_slice()
            } else {
                default_blastp_tabular_fields()
            };
            write_blastp_tabular_output(
                pairwise_hits
                    .as_ref()
                    .expect("pairwise hits prepared for blastp tabular outfmt"),
                fields,
                outfmt,
                args.out.as_ref(),
                &query_ids,
                &subject_ids,
                &context,
            )?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::args::{BlastpCompBasedStats, BlastpSegSpec};

    fn make_blastp_args() -> BlastpArgs {
        BlastpArgs {
            query: PathBuf::from("query.faa"),
            subject: PathBuf::from("subject.faa"),
            task: "blastp".to_string(),
            evalue: None,
            threshold: None,
            word_size: None,
            num_threads: 1,
            out: None,
            max_target_seqs: 500,
            max_hsps_per_subject: 0,
            ungapped: false,
            window_size: None,
            matrix: None,
            gap_open: None,
            gap_extend: None,
            comp_based_stats: None,
            seg: None,
            use_sw_tback: false,
            outfmt: "6".to_string(),
        }
    }

    fn make_pairwise_hit() -> PairwiseHit {
        PairwiseHit {
            hit: Hit {
                identity: 75.0,
                length: 4,
                mismatch: 1,
                gapopen: 1,
                q_start: 2,
                q_end: 5,
                s_start: 7,
                s_end: 10,
                e_value: 1.0e-25,
                bit_score: 42.5,
                num_ident: 3,
                query_frame: 0,
                query_length: 12,
                q_idx: 0,
                s_idx: 0,
                raw_score: 80,
                gap_info: None,
                num_positives: 3,
            },
            query_seq: Some("AB-C".to_string()),
            subject_seq: Some("ABDC".to_string()),
            query_frame: Some(1),
            subject_frame: Some(1),
            positives: Some(3),
            gaps: Some(1),
            subject_length: Some(30),
            subject_title: Some("subject description".to_string()),
        }
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:39-45
    // ```c
    // const char* kDfltArgTabularOutputFmt =
    //     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
    //     "evalue bitscore";
    // const char* kDfltArgTabularOutputFmtTag("std");
    // ```
    #[test]
    fn test_parse_blastp_tabular_fields_expands_std_and_dedups() {
        let fields = parse_blastp_tabular_fields("std qlen qaccver").expect("parsed fields");
        assert_eq!(fields[0], BlastpTabularField::QueryAccessionVersion);
        assert_eq!(fields[1], BlastpTabularField::SubjectAccessionVersion);
        assert!(fields.contains(&BlastpTabularField::QueryLength));
        assert_eq!(
            fields
                .iter()
                .filter(|field| **field == BlastpTabularField::QueryAccessionVersion)
                .count(),
            1
        );
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:985-1027
    // ```c
    // if (num_matches > 0) {
    //     btop_string += NStr::Int8ToString(num_matches);
    //     num_matches = 0;
    // }
    // btop_string += m_QuerySeq[i];
    // btop_string += m_SubjectSeq[i];
    // ```
    #[test]
    fn test_build_btop_matches_ncbi_style_match_run_encoding() {
        assert_eq!(build_btop("AB-C", "ABDC"), "2-D1");
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
    // ```c
    // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
    //     if (iter != m_FieldsToShow.begin())
    //         m_Ostream << m_FieldDelimiter;
    //     x_PrintField(*iter);
    // }
    // ```
    #[test]
    fn test_blastp_tabular_field_text_renders_requested_metadata() {
        let pairwise_hit = make_pairwise_hit();
        assert_eq!(
            blastp_tabular_field_text(
                BlastpTabularField::QueryAccessionVersion,
                &pairwise_hit,
                "query1",
                "subject1"
            ),
            "query1"
        );
        assert_eq!(
            blastp_tabular_field_text(
                BlastpTabularField::SubjectTitle,
                &pairwise_hit,
                "query1",
                "subject1"
            ),
            "subject description"
        );
        assert_eq!(
            blastp_tabular_field_text(
                BlastpTabularField::Btop,
                &pairwise_hit,
                "query1",
                "subject1"
            ),
            "2-D1"
        );
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_options.c:1163-1169
    // ```c
    // if ((program_number == eBlastTypeTblastn ||
    //      program_number == eBlastTypeBlastp ||
    //      program_number == eBlastTypeBlastx) &&
    //     word_size > 5)
    //     options->lut_type = eCompressedAaLookupTable;
    // ```
    #[test]
    fn test_validate_requested_blastp_support_rejects_compressed_lookup_path() {
        let mut args = make_blastp_args();
        args.word_size = Some(5);
        let resolved = args.resolve().expect("resolved blastp args");
        let err = validate_requested_blastp_support(&resolved).expect_err("unsupported word size");
        assert!(err.to_string().contains("compressed amino-acid lookup"));
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3036-3062
    // ```c
    // ECompoAdjustModes compo_adjust_mode =
    //     (ECompoAdjustModes) extendParams->options->compositionBasedStats;
    // ```
    #[test]
    fn test_validate_requested_blastp_support_rejects_unified_p_mode() {
        let mut args = make_blastp_args();
        args.comp_based_stats = Some(BlastpCompBasedStats {
            mode: BlastpCompositionMode::CompositionMatrixAdjust,
            unified_p: true,
        });
        let resolved = args.resolve().expect("resolved blastp args");
        let err = validate_requested_blastp_support(&resolved).expect_err("unsupported unified p");
        assert!(err.to_string().contains("only --comp_based_stats 2"));
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319
    // ```c
    // s_BlastProtGappedAlignment(..., const BlastScoringParameters* score_params,
    //                            BlastGapAlignStruct* gap_align, ...)
    // ```
    #[test]
    fn test_validate_requested_blastp_support_rejects_non_default_scoring() {
        let mut args = make_blastp_args();
        args.matrix = Some("PAM70".to_string());
        args.gap_open = Some(10);
        args.gap_extend = Some(1);
        args.seg = Some(BlastpSegSpec::No);
        let resolved = args.resolve().expect("resolved blastp args");
        let err = validate_requested_blastp_support(&resolved).expect_err("unsupported scoring");
        assert!(err
            .to_string()
            .contains("matrix-aware protein gap alignment"));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2387-2415
    // ```c
    // gapping_params->x_dropoff = (Int4)
    //     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
    //         extendParams->gap_x_dropoff_final);
    // ```
    #[test]
    fn test_blastp_redo_x_dropoff_uses_scaled_final_vs_unscaled_extend_max() {
        let gapped = crate::stats::KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.0,
            beta: -1.0,
        };
        let scaled_lambda = gapped.lambda / 32.0;
        let unscaled_final = x_drop_raw_score(X_DROP_GAPPED_FINAL as f64, &gapped, 1.0)
            .max(x_drop_raw_score(X_DROP_GAPPED_PRELIM as f64, &gapped, 1.0));

        assert_eq!(
            blastp_redo_x_dropoff(scaled_lambda, &gapped, unscaled_final),
            x_drop_raw_score(
                X_DROP_GAPPED_FINAL as f64,
                &crate::stats::KarlinParams {
                    lambda: scaled_lambda,
                    k: gapped.k,
                    h: gapped.h,
                    alpha: gapped.alpha,
                    beta: gapped.beta,
                },
                1.0,
            )
            .max(unscaled_final)
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2438-2462
    // ```c
    // near_identical_cutoff =
    //     (near_identical_cutoff_bits * NCBIMATH_LN2)
    //         / context->sbp->kbp_gap[index]->Lambda;
    // ...
    // if (do_link_hsps) {
    //     cutoff_s =
    //         (int) (hitParams->cutoff_score_min * context->localScalingFactor);
    // } else {
    //     cutoff_s = 1;
    // }
    // ```
    #[test]
    fn test_blastp_redo_params_follow_ncbi_cutoff_formulas() {
        let scaled_lambda = 0.267 / 32.0;
        let near_identical_cutoff = blastp_redo_near_identical_cutoff(scaled_lambda);
        assert!(near_identical_cutoff > 0.0);
        assert_eq!(blastp_redo_cutoff_score(false, 42, 32.0), 1);
        assert_eq!(blastp_redo_cutoff_score(true, 42, 32.0), 1344);
        assert_eq!(blastp_redo_cutoff_score(true, 7, 3.2), 22);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:540-555
    // ```c
    // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    // Blast_HSPListSortByScore(hsp_list);
    // ```
    #[test]
    fn test_prepare_preliminary_hits_for_kappa_purges_shared_endpoints_then_sorts() {
        let hits = vec![
            BlastpPreliminaryHsp {
                query_start: 10,
                query_end: 60,
                subject_start: 20,
                subject_end: 70,
                gapped_query_start: 30,
                gapped_subject_start: 40,
                raw_score: 80,
                query_frame: 0,
                query_length: 100,
                q_idx: 0,
                s_idx: 0,
            },
            BlastpPreliminaryHsp {
                query_start: 10,
                query_end: 80,
                subject_start: 20,
                subject_end: 90,
                gapped_query_start: 35,
                gapped_subject_start: 45,
                raw_score: 90,
                query_frame: 0,
                query_length: 100,
                q_idx: 0,
                s_idx: 0,
            },
            BlastpPreliminaryHsp {
                query_start: 5,
                query_end: 40,
                subject_start: 7,
                subject_end: 42,
                gapped_query_start: 20,
                gapped_subject_start: 22,
                raw_score: 120,
                query_frame: 0,
                query_length: 100,
                q_idx: 0,
                s_idx: 0,
            },
        ];

        let prepared = prepare_preliminary_hits_for_kappa(hits);
        assert_eq!(prepared.len(), 2);
        assert_eq!(prepared[0].raw_score, 120);
        assert_eq!(prepared[1].raw_score, 90);
        assert_eq!(prepared[1].query_start, 10);
        assert_eq!(prepared[1].subject_start, 20);
        assert_eq!(prepared[1].query_end, 80);
    }
}
