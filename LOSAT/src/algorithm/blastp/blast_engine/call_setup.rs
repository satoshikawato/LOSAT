//! The setup of a BLASTP call before its search (`run_call`, from its requested formats to its
//! lookup table): the reports that the call's formats ask for (`requested_reports`), the
//! queries' IDs, titles and lengths (`query_meta`), the query frames (`query_frames_of`), the
//! encoded subjects and their titles (`call_subjects`), the scoring values and X-drops
//! (`call_scoring`), and the query batches, the lookup table and the query contexts
//! (`call_lookup`). Each function is the statements of `run_call` at its place, in their order;
//! `run_call` calls them once each, in that order, with its timing clocks between them, and its
//! `?` passes their errors on.
//!
//! NCBI sets a search up in one function, the query data first, then the score block and the
//! lookup table:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:132-136
//! ```c++
//! BlastSetupPreliminarySearchEx(CRef<IQueryFactory> qf,
//!                               CRef<CBlastOptions> options,
//!                               CConstRef<CPssmWithParameters> pssm,
//!                               BlastSeqSrc* seqsrc,
//!                               size_t num_threads)
//! ```
//!
//! This module is a child of `blast_engine` for the private items of that file it reads
//! (`BlastpReport`, `BlastpRanges`, `BlastpReports`, `PrelimOnly`, `BlastpRuntimeLookup`,
//! `BlastpContextLookupBounds`, the format kinds' functions, the query frames' and subjects'
//! encoders, the lookup builders, `x_prelim_block_plan`) and its sibling modules `psi_round`,
//! `pssm_hooks`, `rps_hooks`, `x_lookup_ahead` and `x_subject_meta`.

use std::ops::Range;
use std::sync::Arc;

use anyhow::{bail, Context, Result};

use super::x_lookup_ahead::{XAheadBlock, XAheadLookup};
use super::x_subject_meta::XSubjectMeta;
// `compute_blosum62_ideal_karlin_params` is `blast_engine`'s import (from `stats`), which only this
// module reads.
use super::{
    blastp_aa_lookup, blastp_aa_lookup_on_block_spare, blastp_compressed_lookup,
    blastp_context_lookup_bounds, blastp_is_archive, blastp_is_xml, blastp_query_anchored_view,
    blastp_query_frames_marked, blastp_seqalign_format, blastp_xml2_kind,
    build_blastp_concatenated_query, compute_blosum62_ideal_karlin_params,
    prepare_blastp_subjects_preserving_order, psi_round, pssm_hooks, rps_hooks, x_blocks_contexts,
    x_lookup_ahead, x_prelim_block_plan, BlastpContextLookupBounds, BlastpRanges, BlastpReport,
    BlastpReports, BlastpRuntimeLookup, PrelimOnly,
};
use crate::algorithm::blastp::args::{
    BlastpCompositionMode, BlastpLookupTableType, ResolvedBlastpArgs,
};
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::x_seed_fused::XContextTable;
use crate::algorithm::blastp::{tabular, x_query_blocks, x_seed_fused};
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::x_bpclut_on;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::api::local_blast::ReportOutputs;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::RangedBatch;
use crate::config::ScoringMatrix;
use crate::report::record::ReportNeeds;
use crate::report::OutputFormat;
use crate::stats::{lookup_protein_params, KarlinParams};
use crate::utils::threading::SearchPool;

/// What `query_frames_of`, `call_subjects` and `call_lookup` read of `run_call`: each field is
/// `run_call`'s value of the same name, by reference where it is not `Copy` (`subject_records`
/// are those of a FASTA subject set, none for a database; `x_block_inputs` are a query block's
/// subjects, query frames and subject scan ranges, LOSAT_X_BPQBLOCK). `'a` is the lifetime of
/// what the query frames and the subjects borrow; `'s` that of `run_call`'s own values (`args`,
/// `prelim`), whose borrow ends with the last call.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:132-134
/// ```c++
///     CRef<SBlastSetupData> setup_data =
///         BlastSetupPreliminarySearchEx(query_factory, options, pssm, seqsrc,
///                                       num_threads);
/// ```
pub(super) struct CallSetup<'s, 'a> {
    pub(super) args: &'s ResolvedBlastpArgs,
    pub(super) query_records: &'a [FastaRecord],
    pub(super) subject_records: &'a [FastaRecord],
    pub(super) prelim: &'s Option<PrelimOnly<'a>>,
    pub(super) query_chunk: bool,
    #[allow(clippy::type_complexity)]
    pub(super) x_block_inputs: Option<(
        &'a [EncodedProtein],
        &'a [Vec<QueryFrame>],
        Option<&'a [Vec<(i32, i32)>]>,
    )>,
    pub(super) blastp_parallel_pool: &'s SearchPool<'a>,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) blastp_num_threads: usize,
    pub(super) ranges: &'s BlastpRanges<'a>,
    pub(super) reports: &'s BlastpReports<'a>,
}

/// The reports that the call's output formats ask for, in the formats' order (`BlastpReport`:
/// each `-outfmt` parsed as the single CLI `-outfmt` is; the query-anchored outfmt 1 to 4 are
/// BLASTP's only), and what they need beyond outfmt 0, 6 and 7 (`ReportNeeds`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2795-2798
/// ```c++
/// CFormattingArgs::ParseFormattingString(const CArgs& args,
///                                        EOutputFormat& fmt_type,
///                                        string& custom_fmt_spec,
///                                        string& custom_delim) const
/// ```
pub(super) fn requested_reports(
    args: &ResolvedBlastpArgs,
    outputs: &ReportOutputs<'_>,
) -> Result<(Vec<BlastpReport>, ReportNeeds)> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
    // ```c
    // if (args[kArgOutputFormat]) {
    //     string fmt_choice =
    //         NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
    // ```
    // Each requested output format is parsed exactly as the single CLI `-outfmt` was.
    let requested_formats = outputs
        .formats
        .iter()
        .map(|format| -> Result<BlastpReport> {
            // The tabular reports with a field list or a delimiter, and outfmt 10 and 20.
            if let Some(spec) = tabular::custom_tabular_spec(format.outfmt) {
                return Ok(BlastpReport::CustomTabular(spec));
            }
            if let Some(serial) = blastp_seqalign_format(format.outfmt) {
                return Ok(BlastpReport::Seqalign(serial));
            }
            if blastp_is_xml(format.outfmt) {
                return Ok(BlastpReport::Xml);
            }
            if blastp_is_archive(format.outfmt) {
                return Ok(BlastpReport::Archive);
            }
            if let Some(kind) = blastp_xml2_kind(format.outfmt) {
                return Ok(BlastpReport::Xml2(kind));
            }
            // The query-anchored outfmt 1 to 4: BLASTP's (`Program::writes`; psiblast and rpsblast
            // refuse them, line RG PS15).
            if let Some(view) = blastp_query_anchored_view(args.app, format.outfmt)? {
                return Ok(BlastpReport::Pairwise(view));
            }
            // Outfmt 0, and outfmt 6 and 7 of the `std` fields and a tab (`TabularSpec::of`).
            let (outfmt, custom_fields) =
                OutputFormat::parse(format.outfmt).map_err(anyhow::Error::msg)?;
            if outfmt == OutputFormat::Pairwise && custom_fields.is_some() {
                bail!("blastp outfmt 0 does not accept custom field lists");
            }
            Ok(match outfmt {
                OutputFormat::Pairwise => {
                    BlastpReport::Pairwise(crate::report::pairwise::PairwiseView::Pairwise)
                }
                _ => BlastpReport::Text(outfmt),
            })
        })
        .collect::<Result<Vec<_>>>()?;
    // What the requested formats need beyond outfmt 0, 6 and 7, decided once for the search
    // (`ReportNeeds`): the engine's positives of the composition-adjusted matches
    // (`engine_positives`). The engine's counts of the HSPs that the pyrrolysine recount
    // changes are kept for every format (`EngineCounts`: the Seq-aligns', which outfmt 0 reads
    // back with `-sorthits` and `-sorthsps`).
    let needs = crate::report::record::ReportNeeds {
        // The archive's Seq-aligns are those of outfmt 8.
        seqalign: requested_formats
            .iter()
            .any(|report| matches!(report, BlastpReport::Seqalign(_) | BlastpReport::Archive)),
    };
    Ok((requested_formats, needs))
}

/// The IDs, titles and lengths of the call's queries by query index, as the reports write them
/// (`query_meta`).
///
/// No NCBI counterpart: NCBI's formatter reads each query's ID, title and length from the
/// query's Bioseq where it writes them; LOSAT makes them once for the call.
pub(super) struct BlastpQueryMeta {
    pub(super) query_ids: Vec<Arc<[u8]>>,
    pub(super) query_titles: Vec<Arc<[u8]>>,
    pub(super) query_lengths: Vec<usize>,
}

/// The IDs, titles and lengths of `query_records` (`BlastpQueryMeta`), which NCBI writes in each
/// query's preamble.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1492-1495
/// ```c++
///     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
///                                             m_Outfile, m_BelieveQuery,
///                                             m_IsHTML, kIsTabularOutput,
///                                             results.GetRID());
/// ```
pub(super) fn query_meta(query_records: &[FastaRecord]) -> BlastpQueryMeta {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:474-504
    // (`s_ReplaceLocalId`)
    // The query ID of the tabular formats is the title's first word, or the local ID
    // (`FastaRecord::shown_id`); the reports write the titles' bytes.
    let query_ids: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.shown_id()))
        .collect();
    let query_titles: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.query_line().as_ref()))
        .collect();
    let query_lengths: Vec<usize> = query_records
        .iter()
        .map(|record| record.seq().len())
        .collect();
    BlastpQueryMeta {
        query_ids,
        query_titles,
        query_lengths,
    }
}

/// The call's query frames (encoded, with SEG where it applies), in `x_owned_frames` where the
/// call makes them: a query block's from the search of all queries (LOSAT_X_BPQBLOCK), those
/// built ahead of the block (LOSAT_X_BPLKPAHEAD, taken from `x_prebuilt`), psiblast's later
/// rounds' (`psi_round`), or the queries' own.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:149-152
/// ```c++
///     // 1. Initialize the query data (borrow it from the factory)
///     CRef<ILocalQueryData> query_data(qf->MakeLocalQueryData(&*options));
///     retval->m_InternalData->m_Queries = query_data->GetSequenceBlk();
///     retval->m_InternalData->m_QueryInfo = query_data->GetQueryInfo();
/// ```
pub(super) fn query_frames_of<'a>(
    setup: &CallSetup<'_, 'a>,
    x_prebuilt: &mut Option<XAheadLookup>,
    x_owned_frames: &'a mut Vec<Vec<QueryFrame>>,
) -> &'a [Vec<QueryFrame>] {
    let CallSetup {
        args,
        query_records,
        x_block_inputs,
        ..
    } = *setup;
    match (x_block_inputs, x_prebuilt.as_mut()) {
        // EXPERIMENT (LOSAT_X_BPQBLOCK): a block of the preliminary stage takes its frames from the
        // search of all queries (its records, encoded with the same parameters).
        (Some((_, frames, _)), _) => frames,
        (None, Some(prebuilt)) => {
            *x_owned_frames = prebuilt.take_query_frames();
            x_owned_frames
        }
        // psiblast's later rounds: the query as the rounds before left it
        // (`psi_round::RoundQuery`).
        (None, None) if psi_round::searched_other_letters(args) => {
            *x_owned_frames = psi_round::round_query_frames(args);
            x_owned_frames
        }
        // Diagnostic (LOSAT_TIMING): this call's encoding and SEG go to the common line
        // (utils/xtiming.rs); frames built ahead (LOSAT_X_BPLKPAHEAD) are not marked. No NCBI
        // counterpart.
        (None, None) => {
            *x_owned_frames = blastp_query_frames_marked(args, query_records, true);
            x_owned_frames
        }
    }
}

/// The call's encoded subjects by subject index, the block search's subject IDs and titles
/// (`x_meta`, LOSAT_X_BPQBLOCK) and the subjects' titles as the reports write them
/// (`call_subjects`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/setup_factory.hpp:337-338
/// ```c++
///     /// The source of subject sequence data
///     CRef<TBlastSeqSrc> m_SeqSrc;
/// ```
pub(super) struct BlastpCallSubjects<'a> {
    pub(super) subjects: &'a [EncodedProtein],
    pub(super) x_meta: Option<&'a XSubjectMeta>,
    pub(super) subject_titles: std::borrow::Cow<'a, [Arc<[u8]>]>,
}

/// The subjects of a FASTA subject set encoded once, in `x_owned_subjects` where the call encodes
/// them (none for a database, whose sequences are read in blocks; the block search's encoded
/// once for all blocks, LOSAT_X_BPQBLOCK), the block search's subject IDs and titles, and the
/// subjects' titles (none for a call of the preliminary stage only).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_query_factory.cpp:101-102
/// ```c++
///     SetupSubjects_OMF(*m_QuerySource, program, &m_SeqBlkVector, &m_MaxLength);
///     m_NumSeqs = static_cast<Uint4>(m_QuerySource->Size());
/// ```
pub(super) fn call_subjects<'a>(
    setup: &CallSetup<'_, 'a>,
    x_owned_subjects: &'a mut Vec<EncodedProtein>,
) -> Result<BlastpCallSubjects<'a>> {
    let CallSetup {
        args,
        subject_records,
        prelim,
        x_block_inputs,
        blastp_parallel_pool,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        blastp_num_threads,
        ..
    } = *setup;
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1426
    // ```c
    // /* iterate over all subject sequences */
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //        continue;
    //    }
    // ```
    // EXPERIMENT (LOSAT_X_BPQBLOCK): the block search encodes the subjects once
    // (`x_query_blocks::XBlockSubjects`), and so does the search of all queries whose preliminary
    // stage runs in blocks (`PrelimOnly::Block`); the same values.
    let subjects: &[EncodedProtein] = if let Some(encoded) = x_block_inputs
        .map(|(subjects, _, _)| subjects)
        .or_else(|| x_query_blocks::x_block_subjects(blastp_parallel_pool, subject_records))
    {
        encoded
    } else {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        let prepared_subjects = prepare_blastp_subjects_preserving_order(
            subject_records,
            args.num_threads,
            blastp_num_threads,
            blastp_parallel_pool,
        )?;
        #[cfg(any(
            not(feature = "parallel"),
            all(target_arch = "wasm32", not(feature = "wasm-threads"))
        ))]
        let prepared_subjects =
            prepare_blastp_subjects_preserving_order(subject_records, args.num_threads)?;
        *x_owned_subjects = prepared_subjects
            .into_iter()
            .map(|prepared_subject| prepared_subject.encoded)
            .collect();
        x_owned_subjects
    };
    // The subjects' titles by subject index, as the reports write them. EXPERIMENT
    // (LOSAT_X_BPQBLOCK): the block search's titles and tabular IDs, made once for all blocks
    // (`blast_engine/x_subject_meta.rs`); the same values, borrowed.
    let x_meta = x_query_blocks::x_block_subject_meta(blastp_parallel_pool, subject_records);
    let subject_titles: std::borrow::Cow<'_, [Arc<[u8]>]> = match x_meta {
        // A call of the preliminary stage only (`PrelimOnly`) writes no report.
        _ if prelim.is_some() => std::borrow::Cow::Owned(Vec::new()),
        Some(meta) => std::borrow::Cow::Borrowed(meta.titles()),
        None => std::borrow::Cow::Owned(
            subject_records
                .iter()
                .map(|record| Arc::from(record.title()))
                .collect(),
        ),
    };
    Ok(BlastpCallSubjects {
        subjects,
        x_meta,
        subject_titles,
    })
}

/// The call's ideal ungapped Lambda (the redo's), its matrix's gapped Karlin-Altschul parameters
/// and its gapped X-drops in raw score, preliminary and final (`call_scoring`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/setup_factory.hpp:324-325
/// ```c++
///     /// BLAST score block structure
///     CRef<TBlastScoreBlk> m_ScoreBlk;
/// ```
pub(super) struct BlastpScoring {
    pub(super) ideal_ungapped_lambda: f64,
    pub(super) gapped_params: KarlinParams,
    pub(super) x_drop_gapped: i32,
    pub(super) x_drop_gapped_final: i32,
}

/// The call's `BlastpScoring` of `args`, with the refusal of an X-drop that NCBI's gapped kernels
/// cannot size before the search.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:422-424
/// ```c
/// Int2 BlastExtensionParametersNew(EBlastProgramType program_number,
///         const BlastExtensionOptions* options, BlastScoreBlk* sbp,
///         BlastQueryInfo* query_info, BlastExtensionParameters* *parameters)
/// ```
pub(super) fn call_scoring(args: &ResolvedBlastpArgs) -> Result<BlastpScoring> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2833-2852
    // ```c
    // Blast_ScoreBlkKbpIdealCalc(BlastScoreBlk* sbp)
    // {
    //    Blast_ResFreqStdComp(sbp, stdrfp);
    //    BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
    //    Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);
    // }
    // ```
    // BLOSUM62's block is computed once per process; another matrix's is the one its lookup
    // table's query preparation computes (`ideal_karlin_params_of`). The redo reads its λ only.
    let ideal_ungapped_lambda = if args.scoring.matrix == ScoringMatrix::Blosum62 {
        compute_blosum62_ideal_karlin_params()
            .map_err(anyhow::Error::msg)
            .context("failed to calculate NCBI blastp kbp_ideal parameters")?
            .lambda
    } else {
        crate::algorithm::tblastx::lookup::ideal_karlin_params_of(args.scoring.matrix).lambda
    };
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
    // The X-drops in bits are `-xdrop_gap` and `-xdrop_gap_final` or BLASTP's defaults
    // (`ResolvedBlastpArgs`), converted as the protein programs convert them (`core::x_drop`); the
    // contexts of a batch share the matrix's gapped Lambda, the smallest one. An X-drop that NCBI's
    // gapped kernels cannot size is refused before the search.
    let x_drop_bits = (args.x_drop_gapped_bits, args.x_drop_final_bits);
    let (x_drop_gapped, x_drop_gapped_final) =
        crate::core::x_drop::gapped_x_dropoffs(x_drop_bits.0, x_drop_bits.1, gapped_params.lambda);
    // The final X-drop reaches the unscaled `ALIGN_EX` only on the ordinary traceback (mode 0,
    // after IDENTITY's reset too; the same choice as `ordinary_traceback` below). A PSSM search
    // converts its X-drops with the PSSM's Lambda and checks them there (`pssm_hooks::pssm_prelim`).
    if !args.ungapped && args.psi_search.is_none() {
        crate::core::x_drop::check_gapped_x_dropoffs(
            args.app.name(),
            x_drop_bits,
            (x_drop_gapped, x_drop_gapped_final),
            (args.scoring.gap_open, args.scoring.gap_extend),
            args.comp_based_stats.mode == BlastpCompositionMode::NoCompositionBasedStats
                || args.rps_search.is_some(),
        )?;
    }
    Ok(BlastpScoring {
        ideal_ungapped_lambda,
        gapped_params,
        x_drop_gapped,
        x_drop_gapped_final,
    })
}

/// The call's query batches and its plan of query blocks (LOSAT_X_BPQBLOCK), its lookup table and
/// query contexts, the contexts' length bounds, the concatenated query, the context of each query
/// offset, which queries are valid and which batches are not searched (`call_lookup`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/setup_factory.hpp:327-329
/// ```c++
///     /// Lookup table, usually only needed in the preliminary stage of the
///     /// search, but for PHI-BLAST it's also needed in the traceback stage.
///     CRef<TLookupTableWrap> m_LookupTable;
/// ```
pub(super) struct BlastpCallLookup {
    pub(super) query_batches: Vec<RangedBatch>,
    pub(super) x_blocks: Option<Vec<Range<usize>>>,
    pub(super) lookup: BlastpRuntimeLookup,
    pub(super) contexts: Vec<QueryContext>,
    pub(super) context_lookup_bounds: BlastpContextLookupBounds,
    pub(super) concatenated_query: Vec<u8>,
    pub(super) x_context_table: XContextTable,
    pub(super) query_valid: Vec<bool>,
    pub(super) query_batch_skipped: Vec<bool>,
}

/// The call's `BlastpCallLookup`: its query batches and block plan, then its lookup table and
/// contexts (rpsblast's contexts for the database's table, psiblast's PSSM table, a block plan's
/// contexts, the table of `-lookup_table_type`, or one built ahead of the block,
/// LOSAT_X_BPLKPAHEAD, taken from `x_prebuilt`), after which the build of the next block's table
/// starts (`x_ahead`), then the bounds, the concatenated query and the context table, and the
/// queries' validity.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:47-53
/// ```c
/// Int2 LookupTableWrapInit(BLAST_SequenceBlk* query,
///         const LookupTableOptions* lookup_options,
///         const QuerySetUpOptions* query_options,
///         BlastSeqLoc* lookup_segments, BlastScoreBlk* sbp,
///         LookupTableWrap** lookup_wrap_ptr, const BlastRPSInfo *rps_info,
///         Blast_Message* *error_msg,
///         BlastSeqSrc* seqsrc)
/// ```
pub(super) fn call_lookup(
    setup: &CallSetup<'_, '_>,
    query_frames: &[Vec<QueryFrame>],
    query_lengths: &[usize],
    x_prebuilt: &mut Option<XAheadLookup>,
    x_ahead: &XAheadBlock<'_>,
) -> Result<BlastpCallLookup> {
    let CallSetup {
        args,
        query_records,
        prelim,
        query_chunk,
        blastp_parallel_pool,
        ranges,
        reports,
        ..
    } = *setup;
    let batch_size = args.query_batch_size()?;
    // The batches count whole records (`seq_range::QueryInput::batches`); the batches
    // here are those that NCBI searches (`read_query_batches`).
    let query_batches = ranges.query_input.batches(batch_size);
    // EXPERIMENT (LOSAT_X_BPQBLOCK): the preliminary stage of the search of all queries in blocks of
    // whole batches, each block a call of its own with its own lookup table (`PrelimOnly::Block`,
    // `x_prelim_block_plan`); `None`: one table for all queries. Planned only by the search of all
    // queries (`BlastpReports::x_plan_blocks`).
    let x_blocks = if prelim.is_none() && reports.x_plan_blocks {
        x_prelim_block_plan(
            args,
            query_records,
            ranges.query_input,
            batch_size,
            blastp_parallel_pool,
        )
    } else {
        None
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2746-2748
    // ```c
    // Boolean check_ideal =
    //    (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
    //     program == eBlastTypeRpsTblastn);
    // ```
    // Exact blastp (`AaLookupTable`) builds no compressed table: LOSAT_X_BPCLUT does not act.
    if matches!(args.lookup_table_type, BlastpLookupTableType::AaLookupTable) && x_bpclut_on() {
        crate::utils::xreach::x_reach!("LOSAT_X_BPCLUT", false);
    }
    let (lookup, mut contexts) = match args.lookup_table_type {
        // rpsblast: the queries' contexts; the search scans with the database's prebuilt table
        // (`rps_hooks::rps_contexts`).
        _ if args.rps_search.is_some() => rps_hooks::rps_contexts(args, query_frames),
        // psiblast `-in_pssm`: the lookup table of the PSSM (`pssm_hooks::pssm_lookup`).
        _ if args.psi_search.is_some() => pssm_hooks::pssm_lookup(args, query_frames),
        // EXPERIMENT (LOSAT_X_BPQBLOCK): the preliminary stage in blocks (`x_blocks`): each block
        // builds its own table; this call has the contexts of all queries (`x_blocks_contexts`).
        _ if x_blocks.is_some() => x_blocks_contexts(args, query_frames, blastp_parallel_pool),
        // `blastp_aa_lookup`. EXPERIMENT (LOSAT_X_BPLKPAHEAD): or the table and contexts built
        // ahead of this block by the same call on the block search's spare vectors
        // (`blast_engine/x_lookup_ahead.rs`). EXPERIMENT (LOSAT_X_BPLKREUSE): or, where nothing is
        // built ahead, built here on the vectors of the block before
        // (`x_query_blocks::x_bp_aa_block_spare`).
        BlastpLookupTableType::AaLookupTable => match x_prebuilt.take() {
            Some(prebuilt) => x_lookup_ahead::x_use_prebuilt(prebuilt),
            None => match x_query_blocks::x_bp_aa_block_spare(blastp_parallel_pool) {
                Some(spare) => blastp_aa_lookup_on_block_spare(
                    args,
                    query_frames,
                    spare,
                    query_chunk,
                    blastp_parallel_pool,
                ),
                None => blastp_aa_lookup(args, query_frames, query_chunk, blastp_parallel_pool),
            },
        },
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1349-1352
        // ```c
        // s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
        //                                 query, locations);
        // s_CompressedLookupFinalize(lookup);
        // ```
        // `blastp_compressed_lookup`. EXPERIMENT (LOSAT_X_BPLKPAHEAD): or the table and contexts built
        // ahead of this block by the same calls (`blast_engine/x_lookup_ahead.rs`).
        BlastpLookupTableType::CompressedAaLookupTable => match x_prebuilt.take() {
            Some(prebuilt) => x_lookup_ahead::x_use_prebuilt(prebuilt),
            None => blastp_compressed_lookup(
                args,
                query_frames,
                &query_batches,
                query_lengths,
                blastp_parallel_pool,
                if query_chunk {
                    x_query_blocks::XLookupFor::QueryChunk
                } else {
                    x_query_blocks::XLookupFor::Search
                },
                if query_chunk {
                    0
                } else {
                    x_query_blocks::x_block_appends(blastp_parallel_pool, query_records)
                },
            )?,
        },
    };
    // EXPERIMENT (LOSAT_X_BPLKPAHEAD): this block's table is ready; the build of the next block's
    // starts now, on the pool, beside this block's search (`blast_engine/x_lookup_ahead.rs`).
    x_ahead.lookup_ready(&lookup);
    let context_lookup_bounds = blastp_context_lookup_bounds(&contexts);
    // EXPERIMENT (LOSAT_X_BPQBLOCK): with the preliminary stage in blocks (`x_blocks`) this call
    // scans no subject (each block's call makes its own query and table), so it makes neither
    // the concatenated query nor the context table (both empty).
    let x_scans = x_blocks.is_none();
    let concatenated_query = if x_scans {
        build_blastp_concatenated_query(&contexts)
    } else {
        Vec::new()
    };
    // The context of every query offset, for the candidates of the fused seed loop (the binary
    // search gives the same index).
    // No NCBI counterpart: a table of BSearchContextInfo's results; it does not change any value
    // NCBI computes.
    let x_context_table = {
        let frame_bases: Vec<i32> = if x_scans {
            contexts.iter().map(|ctx| ctx.frame_base).collect()
        } else {
            Vec::new()
        };
        x_seed_fused::XContextTable::new(&frame_bases, concatenated_query.len())
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:177-180
    // ```c
    //     int status = m_PrelimSearch->CheckInternalData();
    //     if (status != 0)
    //     {
    //          // Search was not run, but we send back an empty CSearchResultSet.
    // ```
    // A query whose Karlin-Altschul parameters cannot be computed is invalid
    // (`is_valid`); a query batch (`GetQueryBatchSize`) of invalid queries only is not
    // searched, and its reports show -1 parameters.
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-640
    // ```c
    //         } catch (const CException& e) {
    //             ...
    //             CRef<CSearchMessage> m
    //                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
    //             messages[index].push_back(m);
    //             s_InvalidateQueryContexts(qinfo, index);
    //         }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2771-2772
    // ```c
    //       if ( !contexts[context].is_valid )
    //           continue;
    // ```
    // A query without letters (`Sequence contains no data`, blast_setup.hpp:190-198) is
    // invalidated when the queries are set up, before its Karlin-Altschul parameters would
    // be computed: it is never searched and has no message of its own about them.
    for context in &mut contexts {
        if context.aa_len == 0 {
            context.is_valid = false;
        }
    }
    let mut query_valid = vec![true; query_records.len()];
    for context in &contexts {
        if !context.is_valid {
            query_valid[context.q_idx as usize] = false;
        }
    }
    let mut query_batch_skipped = vec![false; query_records.len()];
    for batch in &query_batches {
        let range = batch.searched.clone();
        let skipped = query_valid[range.clone()].iter().all(|valid| !valid);
        for flag in &mut query_batch_skipped[range] {
            *flag = skipped;
        }
    }
    Ok(BlastpCallLookup {
        query_batches,
        x_blocks,
        lookup,
        contexts,
        context_lookup_bounds,
        concatenated_query,
        x_context_table,
        query_valid,
        query_batch_skipped,
    })
}
