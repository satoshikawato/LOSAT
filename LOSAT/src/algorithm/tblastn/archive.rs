//! TBLASTN's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:295
//! ```c++
//!                 if (fmt_args->ArchiveFormatRequested(args)) {
//!                     formatter.WriteArchive(*query_factory, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
//!                     m_Bah.ResetMessages();
//!                 } else {
//! ```
//! One archive per query batch (the reader's batch, `BATCH_SIZE` letters, 20000 by default),
//! written after the batch's search. LOSAT's TBLASTN searches the batches first and writes every
//! format at the end, from the hit lists it holds (`stage_e_report::render`), so the archives are
//! written there, one per retained batch with that batch's searched flag; `report::archive` writes
//! the document and this module gathers what TBLASTN holds for its other formats: the Seq-aligns of
//! outfmt 8, the SEG and lower-case masks of the queries, and the Karlin-Altschul blocks and search
//! space of the batch's first query. A query is a protein with one context, so the first valid
//! context of a batch is its first query (a first query without one gives both blocks the `-1`
//! triple), and the batch's `EffectiveSearchSpace` option is that query's search space.

use std::borrow::Cow;
use std::io::Write;
use std::ops::Range;

use anyhow::Result;

use super::seqalign::query_records_of;
use super::stage_d_results::KappaResultHitList;
use super::stage_d_stats::LocalSubjectParameters;
use super::stage_e_report::{ReportRanges, ShownSubjects};
use super::x_subject_store::SubjectLetters;
use super::TblastnArgs;
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::value_parsers::SegSpec;
use crate::report::archive::bioseq::BioseqInput;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::{Filter, OptionInputs};
use crate::report::archive::{
    defline_of, ka_block, write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive,
    KaBlock,
};
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::seqid::SeqId;
use crate::report::serial::seqalign::{write_seq_aligns, Segments};
use crate::stats::tables::KarlinParams;
use crate::utils::seg::SegParams;

/// The recorded options of a TBLASTN run (`options::algorithm_options`'s inputs), from the
/// arguments as given: a value is recorded only when it was given or the class of its argument
/// always sets it (`S/DESIGN_OF8.md` 2.2).
pub(crate) fn option_inputs(args: &TblastnArgs) -> OptionInputs {
    let mut inputs = OptionInputs::new(Program::Tblastn, &args.task);
    inputs.dbsize = args.dbsize;
    inputs.searchsp = args.searchsp;
    inputs.evalue = args.evalue;
    inputs.gapopen = args.gap_open;
    inputs.gapextend = args.gap_extend;
    inputs.matrix = args.matrix.clone();
    inputs.xdrop_gap = args.xdrop_gap;
    inputs.xdrop_gap_final = args.xdrop_gap_final;
    inputs.word_size = args.word_size;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    inputs.max_hsps = args.max_hsps.map(|count| count as i32);
    inputs.culling_limit = args.culling_limit;
    inputs.sum_stats = args.sum_stats;
    // `-threshold` is read as a double and recorded as an `int` (`static_cast<int>`).
    inputs.threshold = args.threshold.map(|threshold| threshold as i32);
    inputs.window_size = args.window_size;
    inputs.ungapped = args.ungapped;
    inputs.soft_masking = args.soft_masking;
    inputs.seg = args
        .seg
        .as_deref()
        .and_then(|value| crate::blastinput::app::parse_seg_option(value, "TBLASTN").ok())
        .map(|spec| match spec {
            SegSpec::No => Filter::No,
            SegSpec::Yes => Filter::Yes,
            SegSpec::WindowLocutHicut {
                window,
                locut,
                hicut,
            } => Filter::Seg {
                window,
                locut,
                hicut,
            },
        });
    inputs.max_intron_length = args.max_intron_length.unwrap_or(0);
    inputs.db_gencode = i32::from(args.db_gencode);
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(500, |size| size as i32);
    if let Ok((mode, _)) =
        crate::blastinput::app::parse_comp_based_stats(&args.comp_based_stats, false, false)
    {
        inputs.comp_based_stats = mode as i32;
    }
    inputs
}

/// The handle of a run that writes the archive, or `None` when no requested format is the
/// archive (`ArchiveHandle::for_program`, with TBLASTN's recorded options).
pub(crate) fn handle_for(
    args: &TblastnArgs,
    archive_requested: bool,
) -> Result<Option<ArchiveHandle>> {
    if !archive_requested {
        return Ok(None);
    }
    // NCBI keeps the name without the white space at its ends (`NStr::TruncateSpaces`,
    // blast_args.cpp:2472-2476), and the request's `database` is that name.
    ArchiveHandle::for_program(&option_inputs(args), args.database_name(), "TBLASTN").map(Some)
}

/// `report::archive::check_titles` for TBLASTN's records.
pub(crate) fn check_titles(records: &[FastaRecord], role: &str) -> Result<()> {
    crate::report::archive::check_titles(records, role, "TBLASTN")
}

pub(crate) use crate::report::archive::error_lines;

/// What the archives of a run need besides the report's own data.
pub(super) struct ArchiveRender<'a> {
    pub handle: &'a ArchiveHandle,
    /// The searched queries of each batch NCBI read and searched, in order (the batches before
    /// a stop; `read_query_batches`).
    pub spans: &'a [Range<usize>],
}

/// The values of a batch's first query that `search-stats`, `ka-blocks` and the
/// `EffectiveSearchSpace` option take.
struct FirstContext {
    /// The ungapped and the gapped block; `None` when the first query has no valid context
    /// (both blocks are then missing).
    karlin: Option<(KarlinParams, KarlinParams)>,
    /// The effective search space and length adjustment of that context (0 without one).
    search_space: i64,
    length_adjustment: i64,
    /// The `EffectiveSearchSpace` option: `-searchsp` in an invalid context too.
    context0_space: i64,
}

/// The runs of lower-case letters of a record (inclusive 0-based positions).
pub(crate) fn lowercase_runs(letters: &[u8]) -> Vec<(i32, i32)> {
    let mut runs = Vec::new();
    let mut start = None;
    for (index, letter) in letters.iter().enumerate() {
        if letter.is_ascii_lowercase() {
            start.get_or_insert(index);
        } else if let Some(from) = start.take() {
            runs.push((from as i32, index as i32 - 1));
        }
    }
    if let Some(from) = start {
        runs.push((from as i32, letters.len() as i32 - 1));
    }
    runs
}

/// The masked regions of each query as NCBI's `TMaskedQueryRegions` holds them (`Map`,
/// api/blast_aux.cpp:825-842): the SEG intervals and, with `-lcase_masking`, the lower-case runs of
/// the searched letters, merged where they overlap, a mask that covers the whole query dropped,
/// and the rest moved into the query's record by the start of its interval; inclusive 0-based
/// positions.
///
/// The masks are computed here, per batch, at report time: SEG runs again on every query of the
/// batch (and the lower-case runs are read again), instead of the search keeping each query's
/// mask until the last archive is written. That is the price of not holding the masks of all
/// batches (memory); it is CPU only, for outfmt 11 only, and below
/// NCBI's cost (200,000 peptides: 24 s user, NCBI 87 s).
fn query_record_masks(
    queries: &[FastaRecord],
    first: usize,
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    ranges: &ReportRanges<'_>,
) -> Vec<Vec<(i32, i32)>> {
    queries
        .iter()
        .enumerate()
        .map(|(offset, query)| {
            let q_idx = first + offset;
            let letters = query.seq();
            let intervals =
                super::search_seed::tblastn_report_mask_intervals(letters, seg, mask_lowercase)
                    .into_iter()
                    .map(|(start, end)| (start as i32, end as i32 - 1))
                    .collect();
            let offset = ranges.queries.offset(q_idx) as i32;
            crate::algorithm::blastx::query_setup::combine_masks(intervals)
                .into_iter()
                .filter(|&(from, to)| from <= to && !(from == 0 && to == letters.len() as i32 - 1))
                .map(|(from, to)| (from + offset, to + offset))
                .collect()
        })
        .collect()
}

/// The Seq-aligns of one batch: the elements of outfmt 8, one query at a time.
struct TblastnAligns<'a> {
    run: RunRecord<'a>,
    queries: &'a [FastaRecord],
    first: usize,
    hitlists: &'a [KappaResultHitList],
    subject_records: &'a [FastaRecord],
    subject_letters: SubjectLetters<'a>,
    shown: &'a ShownSubjects<'a>,
    ranges: &'a ReportRanges<'a>,
}

impl<'a> AlignSource for TblastnAligns<'a> {
    fn write<E: crate::report::serial::encoder::Encoder>(
        &mut self,
        out: &mut ArchiveOut<'_, E>,
    ) -> Result<()> {
        let mut segments = Segments::default();
        let mut subjects: Vec<SubjectRecord<'a>> = Vec::new();
        let mut hsps: Vec<HspRecord<'a>> = Vec::new();
        let hitlists: &'a [KappaResultHitList] = self.hitlists;
        for (offset, hitlist) in hitlists.iter().enumerate() {
            let q_idx = self.first + offset;
            query_records_of(
                hitlist,
                self.subject_records,
                self.subject_letters,
                self.shown,
                self.ranges,
                &mut subjects,
                &mut hsps,
            )?;
            if hsps.is_empty() {
                continue;
            }
            let query = &self.queries[q_idx];
            let query = QueryRecord {
                id: query.seq_id(),
                searched_length: i32::try_from(query.seq().len())?,
                offset: i32::try_from(self.ranges.queries.offset(q_idx))?,
            };
            out.write_with(|enc, flush| {
                write_seq_aligns(
                    enc,
                    &self.run,
                    &query,
                    &subjects,
                    &hsps,
                    &mut segments,
                    |enc| flush(enc),
                );
            });
        }
        Ok(())
    }
}

/// The report data of a run, as `stage_e_report::render` holds it.
pub(super) struct ArchiveData<'a> {
    pub queries: &'a [FastaRecord],
    pub subject_records: &'a [FastaRecord],
    pub subject_letters: SubjectLetters<'a>,
    pub shown: &'a ShownSubjects<'a>,
    pub hitlists: &'a [KappaResultHitList],
    pub parameters: &'a LocalSubjectParameters,
    pub ungapped_karlin: &'a [KarlinParams],
    pub query_validity: &'a [bool],
    pub query_batch_skipped: &'a [bool],
    pub gapped: KarlinParams,
    pub seg: Option<&'a SegParams>,
    pub mask_lowercase: bool,
    pub ranges: &'a ReportRanges<'a>,
}

/// Writes the archive of each batch of `render` to `writer`; `warnings` writes the warnings of a
/// batch's queries (to the diagnostics, which saves them as the batch's messages) before its
/// archive.
pub(super) fn write_tblastn_archives(
    render: &ArchiveRender<'_>,
    data: &ArchiveData<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    writer: &mut dyn Write,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Tblastn,
        db_title: data.shown.db_title(),
        gapped: true,
    };
    // The Bioseq inputs borrow their deflines: keep them alive for the whole call.
    let sequences = matches!(data.shown, ShownSubjects::Fasta { .. });
    let subject_deflines: Vec<Cow<'_, [u8]>> = if sequences {
        data.subject_records.iter().map(defline_of).collect()
    } else {
        Vec::new()
    };
    let subjects_input: Vec<BioseqInput<'_>> = if sequences {
        data.subject_records
            .iter()
            .zip(&subject_deflines)
            .map(|(record, defline)| BioseqInput {
                id: &record.local_id,
                ids: record.believed().map(|ids| ids.ids.as_slice()),
                seq_loc: None,
                defline,
                title: &record.title,
                letters: record.archive_letters(),
                protein: false,
            })
            .collect()
    } else {
        Vec::new()
    };
    for span in render.spans.iter().cloned() {
        if span.is_empty() {
            continue;
        }
        if let Some(warnings) = warnings.as_deref_mut() {
            for index in span.clone() {
                warnings.before_query(index, &mut *writer)?;
            }
        }
        let records = &data.queries[span.clone()];
        let deflines: Vec<Cow<'_, [u8]>> = records.iter().map(defline_of).collect();
        let queries: Vec<BioseqInput<'_>> = records
            .iter()
            .zip(&deflines)
            .map(|(record, defline)| BioseqInput {
                id: &record.local_id,
                ids: record.believed().map(|ids| ids.ids.as_slice()),
                seq_loc: record.archive_seq_loc(),
                defline,
                title: &record.title,
                letters: record.archive_letters(),
                protein: true,
            })
            .collect();
        let first = span.start;
        let searched = !data.query_batch_skipped[first];
        // Only a searched batch has results masks (the search made them): they are computed per batch.
        let masks = if searched {
            query_record_masks(records, first, data.seg, data.mask_lowercase, data.ranges)
        } else {
            Vec::new()
        };

        // `RequiredStart` / `RequiredEnd`: a `-query_loc` that is not the first query's whole
        // length.
        let placements = data.ranges.queries;
        let required = if placements.ranged() {
            let placement = &placements.0[first];
            let cut = data.queries[first].seq().len();
            (cut != placement.length)
                .then_some((placement.offset as i32, (placement.offset + cut) as i32 - 1))
        } else {
            None
        };

        // The results' masks are those of the search, so only a searched batch has them: one
        // `notset` mask per query with its SEG and lower-case intervals (a protein query).
        let mut results_masks: Vec<MaskEntry> = Vec::new();
        if searched {
            for (offset, intervals) in masks.into_iter().enumerate() {
                if !intervals.is_empty() {
                    results_masks.push(MaskEntry {
                        id: records[offset].local_id.clone(),
                        best: records[offset].believed().map(|ids| ids.best.clone()),
                        frame: Frame::NotSet,
                        intervals,
                    });
                }
            }
        }
        // `-lcase_masking`: the `LCaseMask` parameters of the request are the lower-case runs of
        // each whole query record, one `notset` mask per query that has any (`SEG` is never one).
        let mut lcase_masks: Vec<MaskEntry> = Vec::new();
        if data.mask_lowercase {
            for record in records {
                let runs = lowercase_runs(record.archive_letters());
                if !runs.is_empty() {
                    lcase_masks.push(MaskEntry {
                        id: record.local_id.clone(),
                        best: record.believed().map(|ids| ids.best.clone()),
                        frame: Frame::NotSet,
                        intervals: runs,
                    });
                }
            }
        }

        let context = if data.query_validity[first] {
            let length = &data.parameters.lengths[first];
            FirstContext {
                karlin: Some((data.ungapped_karlin[first].clone(), data.gapped.clone())),
                search_space: length.eff_searchsp,
                length_adjustment: length.length_adjustment,
                context0_space: length.eff_searchsp,
            }
        } else {
            FirstContext {
                karlin: None,
                search_space: 0,
                length_adjustment: 0,
                context0_space: data.parameters.lengths[first].eff_searchsp,
            }
        };
        let ka_blocks = match &context.karlin {
            Some((ungapped, gapped)) => [ka_block(ungapped, false), ka_block(gapped, true)],
            None => [KaBlock::missing(false), KaBlock::missing(true)],
        };
        let batch = BatchArchive {
            queries: &queries,
            subjects: &subjects_input,
            required,
            lcase_masks: &lcase_masks,
            masks: &results_masks,
            ka_blocks,
            search_space: context.search_space,
            length_adjustment: context.length_adjustment,
            searched_space: searched.then_some(context.context0_space),
        };
        let mut aligns = TblastnAligns {
            run: RunRecord {
                program: run.program,
                db_title: run.db_title,
                gapped: run.gapped,
            },
            queries: data.queries,
            first,
            hitlists: &data.hitlists[span.clone()],
            subject_records: data.subject_records,
            subject_letters: data.subject_letters,
            shown: data.shown,
            ranges: data.ranges,
        };
        write_archive(render.handle, &batch, &mut aligns, writer)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{try_parse_from, Cli, Commands};
    use crate::report::archive::options::{algorithm_options, ParamValue};
    use crate::report::serial::real::double_to_string_posix;

    // The `algorithm-options` that NCBI BLAST+ 2.17.0 wrote for tblastn (S/of8/b2/out, the OF8b2
    // probes; `check_archive_layout.py opts`), without the `EffectiveSearchSpace` that each
    // archive adds: the name, the options given besides the inputs, and the list as
    // `Name=variant:value`. LOSAT's TBLASTN refuses many of these options at the command line or
    // in its limits (word sizes 2 and 4, matrices, gap costs,
    // `-comp_based_stats 1` and `3`, `-ungapped`, `-max_intron_length`, `-subject_besthit`,
    // `-use_sw_tback`): the list that NCBI would record is what the archive's model knows
    // about them, pinned here.
    const CASES: &[(&str, &[&str], &str)] = &[
        ("def", &[], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbs0", &["-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("cbs0u", &["-comp_based_stats", "0u"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("cbs1", &["-comp_based_stats", "1"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:1 SmithWatermanMode=boolean:False"),
        ("cbs1u", &["-comp_based_stats", "1u"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:1 SmithWatermanMode=boolean:False"),
        ("cbs2", &["-comp_based_stats", "2"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbs2u", &["-comp_based_stats", "2u"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbs3", &["-comp_based_stats", "3"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:3 SmithWatermanMode=boolean:False"),
        ("cbs3sw", &["-comp_based_stats", "3", "-use_sw_tback"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:3 SmithWatermanMode=boolean:True"),
        ("cbs3u", &["-comp_based_stats", "3u"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:3 SmithWatermanMode=boolean:False"),
        ("cbsD", &["-comp_based_stats", "D"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbsF", &["-comp_based_stats", "F"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("cbst", &["-comp_based_stats", "t"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("fast", &["-task", "tblastn-fast"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("gc11", &["-db_gencode", "11"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:11 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("gc2", &["-db_gencode", "2"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:2 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("gc4", &["-db_gencode", "4"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:4 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("mts_warn", &["-max_target_seqs", "2"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:2 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_besthit", &["-subject_besthit"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 UseSubjectBestHit=boolean:True WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_cull0", &["-culling_limit", "0"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:0 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_cull5", &["-culling_limit", "5"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_cull_mts", &["-culling_limit", "3", "-max_target_seqs", "2"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:3 WindowSize=integer:40 HitlistSize=integer:2 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_dbsize", &["-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ev", &["-evalue", "0.01"], "EvalueThreshold=cutoff:e-value:0.01 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_evbig", &["-evalue", "1000"], "EvalueThreshold=cutoff:e-value:1000 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_fast_m80", &["-task", "tblastn-fast", "-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 GapOpeningCost=integer:10 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM80 WindowSize=integer:25 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_fast_thr", &["-task", "tblastn-fast", "-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:15 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_fast_ws", &["-task", "tblastn-fast", "-word_size", "3"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:3 WordThreshold=integer:20 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_gapext", &["-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapExtensionCost=integer:2 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_gapopen", &["-gapopen", "10"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_gaps", &["-gapopen", "9", "-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:2 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_hsps", &["-max_hsps", "1"], "EvalueThreshold=cutoff:e-value:10 MaxHspsPerSubject=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_intron0", &["-max_intron_length", "0"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_intron", &["-max_intron_length", "100"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:100 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_m45", &["-matrix", "BLOSUM45"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:14 GapExtensionCost=integer:2 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM45 WordThreshold=integer:16 WindowSize=integer:60 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_m45cbs0", &["-matrix", "BLOSUM45", "-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:14 GapExtensionCost=integer:2 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM45 WordThreshold=integer:16 WindowSize=integer:60 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("opt_m62", &["-matrix", "BLOSUM62"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM62 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_m62lc", &["-matrix", "blosum62"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:blosum62 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_m80", &["-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM80 WordThreshold=integer:14 WindowSize=integer:25 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_many", &["-evalue", "1e-5", "-threshold", "12", "-window_size", "30", "-gapopen", "10", "-gapextend", "1", "-max_hsps", "2", "-qcov_hsp_perc", "30", "-db_gencode", "4", "-seg", "12 2.2 2.5"], "EvalueThreshold=cutoff:e-value:1e-05 GapOpeningCost=integer:10 GapExtensionCost=integer:1 QueryCovHspPerc=real:30 MaxHspsPerSubject=integer:2 DbGeneticCode=integer:4 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:30 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_mts1", &["-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:1 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_mts4", &["-max_target_seqs", "4"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:4 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_naln", &["-num_alignments", "3"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:3 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ndesc", &["-num_descriptions", "3"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ndescaln", &["-num_descriptions", "3", "-num_alignments", "4"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:4 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_pam30", &["-matrix", "PAM30"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:PAM30 WordThreshold=integer:18 WindowSize=integer:15 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_pam70", &["-matrix", "PAM70"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:PAM70 WordThreshold=integer:16 WindowSize=integer:20 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_parse", &["-parse_deflines"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_qcov", &["-qcov_hsp_perc", "50"], "EvalueThreshold=cutoff:e-value:10 QueryCovHspPerc=real:50 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_soft", &["-soft_masking", "true"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:True SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_softf", &["-soft_masking", "false"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_sumstats", &["-sum_stats", "true"], "EvalueThreshold=cutoff:e-value:10 SumStatistics=boolean:True DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_sumstatsf", &["-sum_stats", "false"], "EvalueThreshold=cutoff:e-value:10 SumStatistics=boolean:False DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_thr13", &["-threshold", "13"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_thr_frac", &["-threshold", "11.7"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_thr", &["-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:15 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ungapped_cbs0", &["-ungapped", "-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:True LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("opt_ungapped_cbsF", &["-ungapped", "-comp_based_stats", "F"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:True LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("opt_win0", &["-window_size", "0"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:0 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_win", &["-window_size", "30"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:30 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws2", &["-word_size", "2"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:2 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws3", &["-word_size", "3"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:3 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws4", &["-word_size", "4"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:4 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws5", &["-word_size", "5"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:19 WordSize=integer:5 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws6", &["-word_size", "6"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:21 WordSize=integer:6 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws7", &["-word_size", "7", "-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:15 WordSize=integer:7 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_xdrop", &["-xdrop_gap", "10", "-xdrop_gap_final", "30"], "EvalueThreshold=cutoff:e-value:10 GapXDropoff=real:10 GapXDropoffFinal=real:30 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_xdropu", &["-xdrop_ungap", "20"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("qloc", &["-query_loc", "10-100"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("qloc_full", &["-query_loc", "1-140"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("qsloc", &["-query_loc", "10-100", "-subject_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("sloc", &["-subject_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("sw", &["-use_sw_tback"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:True"),
        // `-seg yes` and `-seg "10 1.8 2.1"` (frozen by the dropped rows `sa.arch.tbn.segyes` / `.segpar`, S1).
        ("segyes", &["-seg", "yes"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFiltering=boolean:True WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("segpar", &["-seg", "10 1.8 2.1"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:10 SegFilteringLocut=real:1.8 SegFilteringHicut=real:2.1 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
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

    /// The recorded inputs of a command line. The options that LOSAT's TBLASTN does not accept
    /// (they stop at the command line) are applied to the inputs by hand.
    fn inputs_of(extra: &[&str]) -> Option<OptionInputs> {
        let mut parsed: Vec<&str> = Vec::new();
        let mut by_hand: Vec<(&str, &str)> = Vec::new();
        let mut index = 0;
        while index < extra.len() {
            match extra[index] {
                "-subject_besthit" | "-use_sw_tback" => {
                    by_hand.push((extra[index], ""));
                    index += 1;
                }
                // Not recorded (the ungapped X-drop is never set by the application).
                "-xdrop_ungap" => index += 2,
                "-parse_deflines" => index += 1,
                other => {
                    parsed.push(other);
                    index += 1;
                }
            }
        }
        let mut argv = vec!["losat", "tblastn", "-query", "q.fa", "-subject", "s.fa"];
        argv.extend_from_slice(&parsed);
        argv.extend_from_slice(&["-outfmt", "11"]);
        let cli: Cli = try_parse_from(argv).ok()?;
        let Commands::Tblastn(args) = cli.command else {
            return None;
        };
        let mut inputs = option_inputs(&args);
        for (name, _) in by_hand {
            match name {
                "-subject_besthit" => inputs.subject_besthit = true,
                _ => inputs.use_sw_tback = true,
            }
        }
        Some(inputs)
    }

    #[test]
    fn the_option_lists_equal_those_of_ncbis_probes() {
        assert!(CASES.len() >= 70);
        for &(name, extra, expected) in CASES {
            let inputs =
                inputs_of(extra).unwrap_or_else(|| panic!("{name} {extra:?} does not parse"));
            assert_eq!(
                render(&algorithm_options(&inputs)),
                expected,
                "{name} {extra:?}"
            );
        }
    }

    #[test]
    fn the_service_is_plain() {
        assert_eq!(
            crate::report::archive::options::service(Program::Tblastn, "tblastn-fast"),
            "plain"
        );
    }

    #[test]
    fn lower_case_runs_are_inclusive_record_positions() {
        assert_eq!(lowercase_runs(b"MWmehSVlk"), [(2, 4), (7, 8)]);
        assert_eq!(lowercase_runs(b"mwME"), [(0, 1)]);
        assert_eq!(lowercase_runs(b"MWME"), []);
        assert_eq!(lowercase_runs(b""), []);
    }
}
