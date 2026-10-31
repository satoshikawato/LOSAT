//! BLASTP's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:281-283
//! ```c++
//!             if (isArchiveFormat) {
//!                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
//!                 m_Bah.ResetMessages();
//!             } else {
//! ```
//! One archive per query batch, written after the batch's search; `report::archive` writes
//! the document, this module gathers what BLASTP holds for its other formats: the hit lists
//! (the Seq-aligns of outfmt 8), the SEG masks, and the Karlin-Altschul blocks and search
//! spaces of the batch's first query. BLASTP has one context per query, so the first valid
//! context of a batch is its first query (a first query without one gives both blocks the
//! `-1` triple), and the batch's `EffectiveSearchSpace` option is that query's search space.

use std::borrow::Cow;
use std::io::Write;
use std::ops::Range;

use anyhow::Result;

use super::args::BlastpArgs;
use super::hsp::BlastpHitList;
use super::seqalign::{query_records, EngineCounts};
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::blastinput::value_parsers::SegSpec;
use crate::report::archive::bioseq::BioseqInput;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::{Filter, OptionInputs};
use crate::report::archive::{
    defline_of, ka_block, write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive,
    KaBlock,
};
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryNames, RunRecord, SubjectNames, SubjectRecord};
use crate::report::serial::seqalign::write_seq_aligns;
use crate::stats::tables::KarlinParams;

/// The recorded options of a BLASTP run (`options::algorithm_options`'s inputs), from the
/// arguments as given: a value is recorded only when it was given or the class of its argument
/// always sets it (`S/DESIGN_OF8.md` 2.2).
pub(crate) fn option_inputs(args: &BlastpArgs) -> OptionInputs {
    let mut inputs = OptionInputs::new(Program::Blastp, &args.task);
    inputs.dbsize = args.dbsize;
    inputs.searchsp = args.searchsp;
    inputs.evalue = args.evalue;
    inputs.gapopen = args.gap_open;
    inputs.gapextend = args.gap_extend;
    inputs.matrix = args.matrix.clone();
    inputs.word_size = args.word_size;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    inputs.max_hsps = args.max_hsps_per_subject.map(|count| count as i32);
    inputs.culling_limit = args.culling_limit;
    // `-threshold` is read as a double and recorded as an `int` (`static_cast<int>`).
    inputs.threshold = args.threshold.map(|threshold| threshold as i32);
    inputs.window_size = args.window_size;
    inputs.ungapped = args.ungapped;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
    // ```c++
    //     if (args[kArgLookupTableMaskingOnly]) {
    //         opt.SetMaskAtHash(args[kArgLookupTableMaskingOnly].AsBoolean());
    //     }
    // ```
    // `MaskAtHash` is recorded with the value given (the class always sets it: false when absent).
    inputs.soft_masking = args.soft_masking;
    inputs.seg = args
        .seg
        .as_deref()
        .and_then(|value| crate::blastinput::app::parse_seg_option(value, "BLASTP").ok())
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
    inputs.hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(500, |size| size as i32);
    if let Ok((mode, unified_p)) =
        crate::blastinput::app::parse_comp_based_stats(&args.comp_based_stats, true, false)
    {
        inputs.comp_based_stats = mode as i32;
        inputs.unified_p = unified_p;
    }
    inputs.use_sw_tback = args.use_sw_tback;
    inputs
}

/// The handle of a run that writes the archive, or `None` when no requested format is the
/// archive (`ArchiveHandle::for_program`, with BLASTP's recorded options).
pub(crate) fn handle_for(
    args: &BlastpArgs,
    archive_requested: bool,
) -> Result<Option<ArchiveHandle>> {
    if !archive_requested {
        return Ok(None);
    }
    // NCBI keeps the name without the white space at its ends (`NStr::TruncateSpaces`,
    // blast_args.cpp:2472-2476), and the request's `database` is that name.
    let database = args.db.as_deref().map(crate::seqdb::db_argument_name);
    ArchiveHandle::for_program(&option_inputs(args), database, "BLASTP").map(Some)
}

/// `report::archive::check_titles` for BLASTP's records.
pub(crate) fn check_titles(records: &[FastaRecord], role: &str) -> Result<()> {
    crate::report::archive::check_titles(records, role, "BLASTP")
}

pub(crate) use crate::report::archive::error_lines;

/// What the archive of a call needs: the call's queries (the batches of one block, or of the
/// whole search) and their data.
pub(crate) struct ArchiveCall<'a> {
    pub handle: &'a ArchiveHandle,
    /// The searched queries of the call (their records keep the whole letters of a `-query_loc`
    /// cut).
    pub queries: &'a [FastaRecord],
    pub query_lengths: &'a [usize],
    pub query_placements: &'a Placements,
    /// `-subject`: every subject record, whole (the Bioseqs of the request).
    pub subject_records: Option<&'a [FastaRecord]>,
    /// The subject names of the Seq-aligns.
    pub subjects: SubjectNames<'a>,
    pub subject_placements: &'a Placements,
    pub hit_lists: &'a [Option<BlastpHitList>],
    pub engine_counts: &'a EngineCounts,
    /// The SEG masks of each query, in record coordinates (`query_record_masks`).
    pub masks: &'a [Vec<(i32, i32)>],
    /// Whether the query's batch was not searched (every query of it invalid).
    pub unsearched: &'a [bool],
    /// The queries of each NCBI batch of the call, in order.
    pub spans: &'a [Range<usize>],
}

/// The values of a batch's first query that `search-stats`, `ka-blocks` and the
/// `EffectiveSearchSpace` option take.
pub(crate) struct FirstContext {
    /// The ungapped block and, for a gapped search, the gapped one; `None` when the first query
    /// has no valid context (both blocks are then missing).
    pub karlin: Option<(KarlinParams, Option<KarlinParams>)>,
    /// The effective search space and length adjustment of that context (0 without one).
    pub search_space: i64,
    pub length_adjustment: i64,
    /// The `EffectiveSearchSpace` option: that of the first query's context, `-searchsp` in an
    /// invalid one too.
    pub context0_space: i64,
}

/// The Seq-aligns of one batch of a call (the elements of outfmt 8, without the Seq-annot).
struct BlastpAligns<'a> {
    run: RunRecord<'a>,
    queries: QueryNames<'a>,
    subjects: SubjectNames<'a>,
    subject_placements: &'a Placements,
    hit_lists: &'a [Option<BlastpHitList>],
    engine_counts: &'a EngineCounts,
    /// The index of the batch's first query in the call.
    first: usize,
}

impl<'a> AlignSource for BlastpAligns<'a> {
    fn write<E: crate::report::serial::encoder::Encoder>(
        &mut self,
        out: &mut ArchiveOut<'_, E>,
    ) -> Result<()> {
        let mut segments = crate::report::serial::seqalign::Segments::default();
        let mut subject_records: Vec<SubjectRecord<'a>> = Vec::new();
        let mut hsp_records: Vec<HspRecord<'a>> = Vec::new();
        let hit_lists: &'a [Option<BlastpHitList>] = self.hit_lists;
        for (offset, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            let q_idx = self.first + offset;
            query_records(
                hit_list,
                q_idx,
                &self.subjects,
                self.subject_placements,
                self.engine_counts,
                self.queries.placements.offset(q_idx) as i32,
                &mut subject_records,
                &mut hsp_records,
            );
            let query = self.queries.record(q_idx);
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
pub(crate) fn write_blastp_archives(
    call: &ArchiveCall<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    first_context: &dyn Fn(usize) -> FirstContext,
    writer: &mut dyn Write,
) -> Result<()> {
    let whole = [0..call.queries.len()];
    let spans: &[Range<usize>] = if call.spans.is_empty() {
        &whole
    } else {
        call.spans
    };
    let run = RunRecord {
        program: Program::Blastp,
        db_title: call.subjects.db_title(),
        // The value cannot differ: `-ungapped` is refused before the search
        // (`validate_requested_blastp_support`), so every BLASTP run is gapped, as in `seqalign.rs`.
        gapped: true,
    };
    // The Bioseq inputs borrow their deflines: keep them alive for the whole call.
    let subject_deflines: Vec<Cow<'_, [u8]>> = call
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
                    ids: record.believed().map(|ids| ids.ids.as_slice()),
                    seq_loc: None,
                    defline,
                    title: &record.title,
                    letters: record.archive_letters(),
                    protein: true,
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
        let searched = !call.unsearched[first];
        let context = first_context(first);

        // `RequiredStart` / `RequiredEnd`: a `-query_loc` that is not the first query's whole
        // length.
        let required = if call.query_placements.ranged() {
            let placement = &call.query_placements.0[first];
            let cut = call.queries[first].seq().len();
            (cut != placement.length)
                .then_some((placement.offset as i32, (placement.offset + cut) as i32 - 1))
        } else {
            None
        };

        // The results' masks are those of the search, so only a searched batch has them: one
        // `notset` mask per query with the SEG intervals (a protein query). BLASTP has no
        // `LCaseMask` parameter (`-lcase_masking` is not supported); SEG is never one.
        let mut masks: Vec<MaskEntry> = Vec::new();
        if searched {
            for index in span.clone() {
                if !call.masks[index].is_empty() {
                    masks.push(MaskEntry {
                        id: call.queries[index].local_id.clone(),
                        best: call.queries[index].believed().map(|ids| ids.best.clone()),
                        frame: Frame::NotSet,
                        intervals: call.masks[index].clone(),
                    });
                }
            }
        }

        let ka_blocks = match &context.karlin {
            Some((ungapped, gapped)) => [
                ka_block(ungapped, false),
                gapped
                    .as_ref()
                    .map_or(KaBlock::missing(true), |gapped| ka_block(gapped, true)),
            ],
            None => [KaBlock::missing(false), KaBlock::missing(true)],
        };
        let batch = BatchArchive {
            queries: &queries,
            subjects: &subjects_input,
            required,
            lcase_masks: &[],
            masks: &masks,
            ka_blocks,
            search_space: context.search_space,
            length_adjustment: context.length_adjustment,
            searched_space: searched.then_some(context.context0_space),
        };
        let mut aligns = BlastpAligns {
            run,
            queries: QueryNames {
                records: call.queries,
                searched_lengths: call.query_lengths,
                placements: call.query_placements,
            },
            subjects: call.subjects,
            subject_placements: call.subject_placements,
            hit_lists: &call.hit_lists[span.clone()],
            engine_counts: call.engine_counts,
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
    use crate::report::archive::options::{algorithm_options, ParamValue};
    use crate::report::serial::real::double_to_string_posix;

    // The `algorithm-options` that NCBI BLAST+ 2.17.0 wrote for blastp (S/of8/a2/out, the OF8a2
    // probes; `check_archive_layout.py opts`), without the `EffectiveSearchSpace` that each
    // archive adds: the name, the options given besides the inputs, and the list as
    // `Name=variant:value`.
    const CASES: &[(&str, &[&str], &str)] = &[
        ("def", &[], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("fast", &["-task", "blastp-fast"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("short", &["-task", "blastp-short"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:16 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbs0", &["-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("cbs1", &["-comp_based_stats", "1"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:1 SmithWatermanMode=boolean:False"),
        ("cbs3", &["-comp_based_stats", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:3 SmithWatermanMode=boolean:False"),
        ("cbs3sw", &["-comp_based_stats", "3", "-use_sw_tback"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:3 SmithWatermanMode=boolean:True"),
        ("cbsF", &["-comp_based_stats", "F"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("cbst", &["-comp_based_stats", "t"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbsD", &["-comp_based_stats", "D"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("cbs2u", &["-comp_based_stats", "2u"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 UnifiedP=integer:1 SmithWatermanMode=boolean:False"),
        ("cbs1u", &["-comp_based_stats", "1U"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:1 UnifiedP=integer:1 SmithWatermanMode=boolean:False"),
        ("cbs0u", &["-comp_based_stats", "0u"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("cbs3u_sw", &["-comp_based_stats", "3u", "-use_sw_tback"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:3 UnifiedP=integer:1 SmithWatermanMode=boolean:True"),
        ("segyes", &["-seg", "yes"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:True WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("segno", &["-seg", "no"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("segpar", &["-seg", "10 1.8 2.1"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFilteringWindow=integer:10 SegFilteringLocut=real:1.8 SegFilteringHicut=real:2.1 WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        // SM-1 probes (ncbi-parity-20261019/SM/probe/arch): `-soft_masking` sets `MaskAtHash` after
        // the task's `ClearFilterOptions` (blast_args.cpp:387-391).
        ("opt_soft", &["-soft_masking", "true"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_softf", &["-soft_masking", "false"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_soft_seg", &["-soft_masking", "true", "-seg", "yes"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True SegFiltering=boolean:True WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_soft_fast", &["-soft_masking", "true", "-task", "blastp-fast"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 MaskAtHash=boolean:True SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_soft_short", &["-soft_masking", "true", "-task", "blastp-short"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:16 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:True WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("segyes_short", &["-seg", "yes", "-task", "blastp-short"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:16 DustFiltering=boolean:False SegFiltering=boolean:True RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt1", &["-evalue", "1e-5", "-threshold", "12", "-word_size", "2", "-max_hsps", "1", "-qcov_hsp_perc", "50", "-window_size", "30"], "EvalueThreshold=cutoff:e-value:1e-05 WordSize=integer:2 QueryCovHspPerc=real:50 MaxHspsPerSubject=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:12 WindowSize=integer:30 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws5", &["-word_size", "5"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:19 WordSize=integer:5 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws6", &["-word_size", "6"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:21 WordSize=integer:6 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws7", &["-word_size", "7", "-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:15 WordSize=integer:7 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_matrix80", &["-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:BLOSUM80 WordThreshold=integer:12 WindowSize=integer:25 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_matrix45", &["-matrix", "BLOSUM45"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:14 GapExtensionCost=integer:2 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:BLOSUM45 WordThreshold=integer:14 WindowSize=integer:60 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_pam30", &["-matrix", "PAM30"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:PAM30 WordThreshold=integer:16 WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_pam70", &["-matrix", "PAM70"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:PAM70 WordThreshold=integer:14 WindowSize=integer:20 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_gaps", &["-gapopen", "9", "-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:2 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_dbsize", &["-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_w0", &["-window_size", "0"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:0 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_mts1", &["-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:1 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_mts4", &["-max_target_seqs", "4"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:4 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_naln3", &["-num_alignments", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:3 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ndesc", &["-num_descriptions", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ndescaln", &["-num_descriptions", "3", "-num_alignments", "4"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:4 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_thr_frac", &["-threshold", "11.7"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ev", &["-evalue", "0.01"], "EvalueThreshold=cutoff:e-value:0.01 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_xdrop", &["-matrix", "BLOSUM62", "-gapopen", "11", "-gapextend", "1"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:BLOSUM62 WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_fast_ws3", &["-task", "blastp-fast", "-word_size", "3"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:3 WordThreshold=integer:20 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_fast_thr", &["-task", "blastp-fast", "-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:15 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_short_matrix", &["-task", "blastp-short", "-matrix", "BLOSUM62"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:BLOSUM62 GapOpeningCost=integer:11 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:11 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_short_ev", &["-task", "blastp-short", "-evalue", "100", "-word_size", "3"], "EvalueThreshold=cutoff:e-value:100 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:3 WordThreshold=integer:16 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_short_cbs", &["-task", "blastp-short", "-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:16 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("opt_short_ws6", &["-task", "blastp-short", "-word_size", "6"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:6 WordThreshold=integer:21 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_short_thr", &["-task", "blastp-short", "-threshold", "12"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:12 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_short_ws4", &["-task", "blastp-short", "-word_size", "4", "-matrix", "BLOSUM62"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:BLOSUM62 GapOpeningCost=integer:11 GapExtensionCost=integer:1 WordSize=integer:4 WordThreshold=integer:11 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws4", &["-word_size", "4"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:4 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ws2_matrix", &["-word_size", "2", "-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 WordSize=integer:2 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:BLOSUM80 WordThreshold=integer:12 WindowSize=integer:25 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_gapopen_only", &["-gapopen", "10"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_gapext_only", &["-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapExtensionCost=integer:2 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_fast_matrix", &["-task", "blastp-fast", "-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 GapOpeningCost=integer:10 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:BLOSUM80 WindowSize=integer:25 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_pam30_thr", &["-matrix", "PAM30", "-threshold", "11"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:PAM30 WordThreshold=integer:11 WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
        ("opt_ungapped_cbs0", &["-ungapped", "-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:True CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
        ("opt_ungapped_cbsF", &["-ungapped", "-comp_based_stats", "F"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:True CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
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

    /// `-task blastp-short` stops at the command line (LOSAT's BLASTP does not implement
    /// it); NCBI's probes of it are applied to the arguments by hand.
    fn parse(extra: &[&str]) -> Option<BlastpArgs> {
        let short = extra
            .windows(2)
            .any(|pair| pair == ["-task", "blastp-short"]);
        let mut argv = vec!["losat", "blastp", "-query", "q.fa", "-subject", "s.fa"];
        let mut index = 0;
        while index < extra.len() {
            if short && extra[index] == "-task" {
                index += 2;
            } else {
                argv.push(extra[index]);
                index += 1;
            }
        }
        argv.extend_from_slice(&["-outfmt", "11"]);
        let cli: Cli = try_parse_from(argv).ok()?;
        match cli.command {
            Commands::Blastp(mut args) => {
                if short {
                    args.task = "blastp-short".into();
                }
                Some(args)
            }
            _ => None,
        }
    }

    #[test]
    fn the_option_lists_equal_those_of_ncbis_probes() {
        assert!(CASES.len() >= 50);
        for &(name, extra, expected) in CASES {
            let args = parse(extra).unwrap_or_else(|| panic!("{name} {extra:?} does not parse"));
            assert_eq!(
                render(&algorithm_options(&option_inputs(&args))),
                expected,
                "{name} {extra:?}"
            );
        }
    }

    #[test]
    fn the_service_is_plain() {
        assert_eq!(
            crate::report::archive::options::service(Program::Blastp, "blastp-short"),
            "plain"
        );
    }
}
