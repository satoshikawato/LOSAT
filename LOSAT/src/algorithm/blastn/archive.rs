//! BLASTN's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:302-306
//! ```c++
//!             if (isArchiveFormat) {
//!                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
//!                 m_Bah.ResetMessages();
//!             } else {
//! ```
//! One archive per query batch, written after the batch's search; `report::archive` writes
//! the document, this module gathers what it holds from the hit lists, the query masks and
//! the Karlin-Altschul blocks that BLASTN keeps for its other formats.

use std::io::Write;
use std::ops::Range;

use anyhow::Result;

use super::args::BlastnArgs;
use super::hsp::BlastnHitList;
use super::pairwise::shown_query_masks;
use super::scoring::ContextKarlin;
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::blastinput::value_parsers::{parse_dust_filtering, DustSpec};
use crate::report::archive::bioseq::BioseqInput;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::{Filter, OptionInputs};
use crate::report::archive::{
    write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive, KaBlock,
};
#[cfg(test)]
use crate::report::diag_log::ERROR;
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryNames, RunRecord, SubjectNames, SubjectRecord};
use crate::report::serial::seqalign::{write_seq_aligns, Segments};
use crate::utils::dust::MaskedInterval;

/// The recorded options of a BLASTN run (`options::algorithm_options`'s inputs). `dust` is
/// the `-dust` value as given (`BlastnArgs::resolve_dust` consumes it).
pub(crate) fn option_inputs(args: &BlastnArgs, dust: Option<&str>) -> OptionInputs {
    let mut inputs = OptionInputs::new(Program::Blastn, &args.task);
    inputs.dbsize = args.dbsize;
    inputs.evalue = args.evalue;
    inputs.gapopen = args.gap_open;
    inputs.gapextend = args.gap_extend;
    inputs.word_size = args.word_size.map(|size| size as i32);
    inputs.perc_identity = args.percent_identity;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    inputs.max_hsps = args.max_hsps_per_subject.map(|count| count as i32);
    inputs.penalty = args.penalty;
    inputs.reward = args.reward;
    inputs.template_type = args.template_type.map(|kind| {
        use crate::algorithm::blastn::disc_lookup::DiscWordType;
        match kind {
            DiscWordType::Coding => 0,
            DiscWordType::Optimal => 1,
            DiscWordType::TwoTemplates => 2,
        }
    });
    inputs.template_length = args.template_length.map(i32::from);
    inputs.dust = dust
        .and_then(|value| parse_dust_filtering(value).ok())
        .map(|spec| match spec {
            DustSpec::No => Filter::No,
            DustSpec::Yes => Filter::Yes,
            DustSpec::Parameters {
                level,
                window,
                linker,
            } => Filter::Dust {
                level: level as i32,
                window: window as u32 as i32,
                linker: linker as u32 as i32,
            },
        });
    inputs.culling_limit = args.culling_limit;
    inputs.subject_besthit = args.subject_besthit;
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(500, |size| size as i32);
    inputs
}

/// The handle of a run that writes the archive, or `None` when no requested format is the
/// archive (`ArchiveHandle::for_program`, with BLASTN's recorded options).
pub(crate) fn handle_for(
    args: &BlastnArgs,
    dust: Option<&str>,
    archive_requested: bool,
) -> Result<Option<ArchiveHandle>> {
    if !archive_requested {
        return Ok(None);
    }
    ArchiveHandle::for_program(&option_inputs(args, dust), args.db.as_deref(), "BLASTN").map(Some)
}

/// `report::archive::check_titles` for BLASTN's records.
pub(crate) fn check_titles(records: &[FastaRecord], role: &str) -> Result<()> {
    crate::report::archive::check_titles(records, role, "BLASTN")
}

pub(crate) use crate::report::archive::error_lines;

/// What the archive of a call needs: the call's queries (one batch, or several when the
/// reports of every batch are written after the search) and their data.
pub(crate) struct ArchiveCall<'a> {
    pub handle: &'a ArchiveHandle,
    pub queries: &'a [FastaRecord],
    pub query_lengths: &'a [usize],
    pub query_placements: &'a Placements,
    /// `-subject`: every subject record, whole (the Bioseqs of the request).
    pub subject_records: Option<&'a [FastaRecord]>,
    /// The subject names of the Seq-aligns.
    pub subjects: SubjectNames<'a>,
    pub subject_placements: &'a Placements,
    pub hit_lists: &'a [Option<BlastnHitList>],
    pub query_masks: &'a [Vec<MaskedInterval>],
    pub unsearched: &'a [bool],
    pub lcase_masking: bool,
    pub both_strands: bool,
    /// The queries of each batch of the call, in order; empty for a call that is one batch.
    pub spans: &'a [Range<usize>],
}

/// The values of a batch's first query that `search-stats`, `ka-blocks` and the
/// `EffectiveSearchSpace` option take.
pub(crate) struct FirstContext {
    /// The ungapped and gapped blocks of the first valid context; `None` when the first
    /// query has no valid context.
    pub karlin: Option<ContextKarlin>,
    /// The effective search space and length adjustment of that context (0 without one).
    pub search_space: i64,
    pub length_adjustment: i64,
    /// The search space of the batch's context 0 (the option).
    pub context0_space: i64,
}

fn ka_block(params: &crate::stats::tables::KarlinParams, gapped: bool) -> KaBlock {
    KaBlock {
        lambda: params.lambda,
        k: params.k,
        h: params.h,
        gapped,
    }
}

fn defline_of(record: &FastaRecord) -> Vec<u8> {
    if record.source.defline.is_empty() {
        let mut line = vec![b'>'];
        line.extend_from_slice(&record.title);
        line
    } else {
        record.source.defline.clone()
    }
}

/// The Seq-aligns of one batch of a call.
struct BlastnAligns<'a> {
    run: RunRecord<'a>,
    queries: QueryNames<'a>,
    subjects: SubjectNames<'a>,
    subject_placements: &'a Placements,
    hit_lists: &'a [Option<BlastnHitList>],
    /// The index of the batch's first query in the call.
    first: usize,
}

impl<'a> AlignSource for BlastnAligns<'a> {
    fn write<E: crate::report::serial::encoder::Encoder>(
        &mut self,
        out: &mut ArchiveOut<'_, E>,
    ) -> Result<()> {
        let mut segments = Segments::default();
        let mut subject_records: Vec<SubjectRecord<'a>> = Vec::new();
        let mut hsp_records: Vec<HspRecord<'a>> = Vec::new();
        let hit_lists: &'a [Option<BlastnHitList>] = self.hit_lists;
        for (offset, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            super::seqalign::query_records(
                hit_list,
                &self.subjects,
                self.subject_placements,
                &mut subject_records,
                &mut hsp_records,
            );
            let query = self.queries.record(self.first + offset);
            out.write_with(|enc, flush| {
                write_seq_aligns(
                    enc,
                    &self.run,
                    &query,
                    &subject_records,
                    &hsp_records,
                    &mut segments,
                    |enc| flush(enc),
                );
            });
        }
        Ok(())
    }
}

/// Writes the archive of each batch of `call` to `writer`; `warnings` writes the warnings
/// of a batch's queries (to the diagnostics, which saves them as the batch's messages)
/// before its archive. `first_context(first_query)` gives the values of a batch's first
/// query.
pub(crate) fn write_blastn_archives(
    call: &ArchiveCall<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    first_context: &dyn Fn(usize) -> Result<FirstContext>,
    writer: &mut dyn Write,
) -> Result<()> {
    let whole = [0..call.queries.len()];
    let spans: &[Range<usize>] = if call.spans.is_empty() {
        &whole
    } else {
        call.spans
    };
    let run = RunRecord {
        program: Program::Blastn,
        db_title: call.subjects.db_title(),
        gapped: true,
    };
    // The Bioseq inputs borrow their deflines: keep them alive for the whole call.
    let subject_deflines: Vec<Vec<u8>> = call
        .subject_records
        .map(|records| records.iter().map(defline_of).collect())
        .unwrap_or_default();
    let subjects_input: Vec<BioseqInput<'_>> = call
        .subject_records
        .map(|records| {
            records
                .iter()
                .zip(&subject_deflines)
                .map(|(record, defline)| BioseqInput {
                    id: &record.local_id,
                    defline,
                    title: &record.title,
                    letters: record.archive_letters(),
                    protein: false,
                })
                .collect()
        })
        .unwrap_or_default();
    for span in spans.iter().cloned() {
        if span.is_empty() {
            continue;
        }
        if let Some(warnings) = warnings.as_deref_mut() {
            for index in span.clone() {
                warnings.before_query(index, &mut *writer)?;
            }
        }
        let records = &call.queries[span.clone()];
        let deflines: Vec<Vec<u8>> = records.iter().map(defline_of).collect();
        let queries: Vec<BioseqInput<'_>> = records
            .iter()
            .zip(&deflines)
            .map(|(record, defline)| BioseqInput {
                id: &record.local_id,
                defline,
                title: &record.title,
                letters: record.archive_letters(),
                protein: false,
            })
            .collect();
        let first = span.start;
        let searched = !call.unsearched[first];
        let context = first_context(first)?;

        let required = if call.query_placements.ranged() {
            let placement = &call.query_placements.0[first];
            let cut = call.queries[first].seq().len();
            (cut != placement.length)
                .then_some((placement.offset as i32, (placement.offset + cut) as i32 - 1))
        } else {
            None
        };

        let mut masks: Vec<MaskEntry> = Vec::new();
        let mut lcase_masks: Vec<MaskEntry> = Vec::new();
        if searched {
            for index in span.clone() {
                let record = &call.queries[index];
                let offset = call.query_placements.offset(index);
                // The results' masks leave out a mask over the whole query (`shown_query_masks`,
                // as the report's masks); the `LCaseMask` parameters keep every interval of
                // the query's masks (probe `multi_err`: a query that DUST masks whole).
                let raw: Vec<(i32, i32)> = call.query_masks[index]
                    .iter()
                    .map(|mask| ((mask.start + offset) as i32, (mask.end + offset) as i32 - 1))
                    .collect();
                if raw.is_empty() {
                    continue;
                }
                let shown: Vec<(i32, i32)> =
                    shown_query_masks(&call.query_masks[index], call.query_lengths[index], offset)
                        .iter()
                        .map(|mask| (mask.start as i32, mask.end as i32 - 1))
                        .collect();
                // Both strands' contexts record the intervals, and the lower-case masks of
                // a `-lcase_masking` run are kept for both strands (probes `in_n_dust*`,
                // `t_n_strandplus_lc`, `m_n_*`).
                let has_lower_case = call.lcase_masking
                    && record
                        .archive_letters()
                        .iter()
                        .any(|letter| letter.is_ascii_lowercase());
                let times = if call.both_strands || has_lower_case {
                    2
                } else {
                    1
                };
                if !shown.is_empty() {
                    masks.push(MaskEntry {
                        id: record.local_id.clone(),
                        frame: Frame::Plus1,
                        intervals: shown,
                    });
                }
                lcase_masks.push(MaskEntry {
                    id: record.local_id.clone(),
                    frame: Frame::Plus1,
                    intervals: raw
                        .iter()
                        .flat_map(|&interval| std::iter::repeat_n(interval, times))
                        .collect(),
                });
            }
        }

        let ka_blocks = match &context.karlin {
            Some(blocks) => [
                ka_block(&blocks.ungapped, false),
                ka_block(&blocks.gapped, true),
            ],
            None => [KaBlock::missing(false), KaBlock::missing(true)],
        };
        let batch = BatchArchive {
            queries: &queries,
            subjects: &subjects_input,
            required,
            lcase_masks: &lcase_masks,
            masks: &masks,
            ka_blocks,
            search_space: context.search_space,
            length_adjustment: context.length_adjustment,
            searched_space: searched.then_some(context.context0_space),
        };
        let mut aligns = BlastnAligns {
            run,
            queries: QueryNames {
                records: call.queries,
                searched_lengths: call.query_lengths,
                placements: call.query_placements,
            },
            subjects: call.subjects,
            subject_placements: call.subject_placements,
            hit_lists: &call.hit_lists[span.clone()],
            first,
        };
        write_archive(call.handle, &batch, &mut aligns, writer)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{try_parse_from, Cli, Commands};
    use crate::report::archive::options::{self, algorithm_options, ParamValue};
    use crate::report::serial::real::double_to_string_posix;

    // The `algorithm-options` that NCBI BLAST+ 2.17.0 wrote in the OF8 design's probes
    // (S/probe8/out, `check_archive_layout.py opts`), without the `EffectiveSearchSpace` that
    // each archive adds: the probe's name, the options given besides the inputs, and the
    // list as `Name=variant:value`.
    const CASES: &[(&str, &[&str], &str)] = &[
    ("b_n_allN_db_fmt11", &[], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("e_n_badpair_fmt11", &["-reward", "2", "-penalty", "-7"], "EvalueThreshold=cutoff:e-value:10 MismatchPenalty=integer:-7 MatchReward=integer:2 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("in_n_dust_minus_fmt11", &["-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("in_n_lcase_mask_fmt11", &["-lcase_masking"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("in_n_mts1_fmt11", &["-task", "blastn", "-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:1"),
    ("in_n_qloc_fmt11", &["-query_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("in_n_qloc_full_fmt11", &["-query_loc", "1-480"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("m_n_dust_plus_fmt11", &["-strand", "plus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("m_n_lc_minus_fmt11", &["-lcase_masking", "-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g1", &["-word_size", "20", "-evalue", "1e-5", "-perc_identity", "90", "-max_hsps", "1", "-qcov_hsp_perc", "50"], "EvalueThreshold=cutoff:e-value:1e-05 WordSize=integer:20 PercentIdentity=real:90 QueryCovHspPerc=real:50 MaxHspsPerSubject=integer:1 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g10", &["-gapopen", "2", "-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:2 GapExtensionCost=integer:2 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g1rev", &["-qcov_hsp_perc", "50", "-max_hsps", "1", "-perc_identity", "90", "-evalue", "1e-5", "-word_size", "20"], "EvalueThreshold=cutoff:e-value:1e-05 WordSize=integer:20 PercentIdentity=real:90 QueryCovHspPerc=real:50 MaxHspsPerSubject=integer:1 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g2", &["-penalty", "-3", "-culling_limit", "1", "-strand", "plus", "-dust", "no"], "EvalueThreshold=cutoff:e-value:10 MismatchPenalty=integer:-3 MaskAtHash=boolean:True DustFiltering=boolean:False UngappedMode=boolean:False Culling=integer:1 HitlistSize=integer:500"),
    ("opt_blastn_g3", &["-reward", "2", "-penalty", "-3", "-gapopen", "5", "-gapextend", "2", "-subject_besthit", "-lcase_masking"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:5 GapExtensionCost=integer:2 MismatchPenalty=integer:-3 MatchReward=integer:2 MaskAtHash=boolean:True UngappedMode=boolean:False UseSubjectBestHit=boolean:True HitlistSize=integer:500"),
    ("opt_blastn_g4", &["-max_target_seqs", "1", "-dust", "20 64 1", "-query_loc", "10-300", "-subject_loc", "10-300", "-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 MaskAtHash=boolean:True DustFilteringLevel=integer:20 DustFilteringWindow=integer:64 DustFilteringLinker=integer:1 UngappedMode=boolean:False HitlistSize=integer:1"),
    ("opt_blastn_g5", &["-num_descriptions", "5", "-num_alignments", "5", "-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:5"),
    ("opt_blastn_g6", &["-task", "dc-megablast", "-template_type", "coding", "-template_length", "16", "-word_size", "12"], "EvalueThreshold=cutoff:e-value:10 MBTemplateType=integer:0 MBTemplateLength=integer:16 WordSize=integer:12 WindowSize=integer:40 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g7", &["-task", "blastn-short"], "EvalueThreshold=cutoff:e-value:1000 MatchReward=integer:1 MismatchPenalty=integer:-3 WordSize=integer:7 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g8", &["-task", "blastn", "-word_size", "7", "-dust", "yes"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:7 MaskAtHash=boolean:True DustFiltering=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("opt_blastn_g9", &["-task", "megablast", "-word_size", "28", "-evalue", "10", "-strand", "both", "-dust", "yes"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:28 MaskAtHash=boolean:True DustFiltering=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("optdb_blastn_g4", &["-max_target_seqs", "1", "-dust", "20 64 1", "-query_loc", "10-300", "-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 MaskAtHash=boolean:True DustFilteringLevel=integer:20 DustFilteringWindow=integer:64 DustFilteringLinker=integer:1 UngappedMode=boolean:False HitlistSize=integer:1"),
    ("t_n_cull0_fmt11", &["-culling_limit", "0"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False Culling=integer:0 HitlistSize=integer:500"),
    ("t_n_dmb_opt_fmt11", &["-task", "dc-megablast", "-template_type", "optimal", "-template_length", "21"], "EvalueThreshold=cutoff:e-value:10 MBTemplateType=integer:1 MBTemplateLength=integer:21 WordSize=integer:11 WindowSize=integer:40 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("t_n_naln_fmt11", &["-num_alignments", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:3"),
    ("t_n_ndesc_fmt11", &["-num_descriptions", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("t_n_pid0_fmt11", &["-perc_identity", "0"], "EvalueThreshold=cutoff:e-value:10 PercentIdentity=real:0 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
    ("t_n_strandplus_lc_fmt11", &["-lcase_masking", "-strand", "plus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
];

    fn posix(value: f64) -> String {
        let mut digits = Vec::new();
        double_to_string_posix(value, 15, &mut digits);
        String::from_utf8(digits).unwrap()
    }

    fn render(params: &[crate::report::archive::options::Param]) -> String {
        params
            .iter()
            .map(|param| {
                let value = match &param.value {
                    ParamValue::Evalue(value) => format!("cutoff:e-value:{}", posix(*value)),
                    ParamValue::Integer(value) => format!("integer:{value}"),
                    ParamValue::BigInteger(value) => format!("big-integer:{value}"),
                    ParamValue::Boolean(true) => "boolean:True".to_string(),
                    ParamValue::Boolean(false) => "boolean:False".to_string(),
                    ParamValue::Real(value) => format!("real:{}", posix(*value)),
                    ParamValue::Str(value) => format!("string:{value}"),
                };
                format!("{}={value}", param.name)
            })
            .collect::<Vec<_>>()
            .join(" ")
    }

    fn parse(extra: &[&str]) -> Option<BlastnArgs> {
        let mut argv = vec!["losat", "blastn", "-query", "q.fa", "-subject", "s.fa"];
        argv.extend_from_slice(extra);
        argv.extend_from_slice(&["-outfmt", "11"]);
        let cli: Cli = try_parse_from(argv).ok()?;
        match cli.command {
            Commands::Blastn(args) => Some(args),
            _ => None,
        }
    }

    #[test]
    fn the_option_lists_equal_those_of_ncbis_probes() {
        let mut compared = 0;
        for &(name, extra, expected) in CASES {
            let args = parse(extra).unwrap_or_else(|| panic!("{name} {extra:?} does not parse"));
            let dust = args.dust_filtering.clone();
            let inputs = option_inputs(&args, dust.as_deref());
            assert_eq!(
                render(&algorithm_options(&inputs)),
                expected,
                "{name} {extra:?}"
            );
            compared += 1;
        }
        assert!(compared >= 25, "{compared}");
    }

    #[test]
    fn the_service_follows_the_task() {
        for (task, service) in [
            ("megablast", "megablast"),
            ("dc-megablast", "megablast"),
            ("blastn", "plain"),
            ("blastn-short", "plain"),
        ] {
            assert_eq!(options::service(Program::Blastn, task), service);
        }
    }

    #[test]
    fn error_lines_are_ncbis_error_texts_only() {
        let native = |message: &str| {
            anyhow::Error::from(crate::cli::NativeError {
                exit: 1,
                message: message.to_string(),
            })
        };
        let lines = error_lines(
            &native(
                "BLAST query/options error: bad  option\nPlease refer to the BLAST+ user manual.\n",
            ),
            false,
        )
        .unwrap();
        assert_eq!(
            lines,
            vec![
                (ERROR, b"BLAST query/options error: bad option".to_vec()),
                (ERROR, b"Please refer to the BLAST+ user manual.".to_vec()),
            ]
        );
        assert!(error_lines(&native("Command line argument error: x\n"), false).is_none());
        let subject = error_lines(
            &native("Command line argument error: Argument \"subject\". File is not accessible:  `a.fa'\n"),
            true,
        )
        .unwrap();
        assert_eq!(
            subject,
            vec![(
                ERROR,
                b"Command line argument error: Argument \"subject\". File is not accessible: `a.fa'".to_vec()
            )]
        );
        assert!(error_lines(&anyhow::anyhow!("unsupported by LOSAT"), false).is_none());
    }
}
