//! BLASTP's BLAST archive (outfmt 11): the request options from the arguments, the archive of
//! each query batch from the batch's report data, and the error archive.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:280-283
//! One archive per query batch, written after the batch's search; `report::archive` writes
//! the document, this module gathers what BLASTP holds for its other formats: the hit lists
//! (the Seq-aligns of outfmt 8), the SEG masks, and the Karlin-Altschul blocks and search
//! spaces of the batch's first query. BLASTP has one context per query, so the first valid
//! context of a batch is its first query (a first query without one gives both blocks the
//! `-1` triple), and the batch's `EffectiveSearchSpace` option is that query's search space.

use std::io::Write;
use std::ops::Range;

use anyhow::Result;

use super::args::BlastpArgs;
use super::hsp::BlastpHitList;
use super::seqalign::{query_records, EngineCounts};
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::core::blast_util::ncbi_int4_from_double;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::options::OptionInputs;
use crate::report::archive::requests::{LcaseFrames, RequestPolicy, RunRequests};
use crate::report::archive::{
    ka_block, write_archive, AlignSource, ArchiveHandle, ArchiveOut, BatchArchive, KaBlock,
};
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryNames, RunRecord, SubjectNames, SubjectRecord};
use crate::report::serial::seqalign::write_seq_aligns;
use crate::stats::tables::KarlinParams;

/// The recorded options of a BLASTP run (`options::algorithm_options`'s inputs), from the
/// arguments as given: a value is recorded only when it was given or the class of its argument
/// always sets it (`S/DESIGN_OF8.md` 2.2). An application whose search is BLASTP's (psiblast,
/// rpsblast) gives its own (`BlastpArgs::request_inputs`), made from these and its arguments.
pub(crate) fn option_inputs(args: &BlastpArgs) -> OptionInputs {
    // The archive's `HitlistSize` is that of outfmt 11, whoever made the other options.
    let hitlist_size = args
        .max_target_seqs
        .or(args.num_alignments)
        .map_or(500, |size| size as i32);
    if let Some(inputs) = &args.request_inputs {
        let mut inputs = OptionInputs::clone(inputs);
        inputs.hitlist_size = hitlist_size;
        return inputs;
    }
    let mut inputs = OptionInputs::new(Program::Blastp, &args.task);
    inputs.dbsize = args.dbsize;
    inputs.searchsp = args.searchsp;
    inputs.evalue = args.evalue;
    inputs.gapopen = args.gap_open;
    inputs.gapextend = args.gap_extend;
    inputs.matrix = args.matrix.clone();
    // `-xdrop_gap` and `-xdrop_gap_final` as given (`SetGapXDropoff`, `SetGapXDropoffFinal`);
    // `-xdrop_ungap` has no entry (NCBI records none, `tblastn/archive.rs` `opt_xdropu`).
    inputs.xdrop_gap = args.xdrop_gap;
    inputs.xdrop_gap_final = args.xdrop_gap_final;
    inputs.word_size = args.word_size;
    inputs.qcov_hsp_perc = args.qcov_hsp_perc;
    inputs.max_hsps = args.max_hsps_per_subject.map(|count| count as i32);
    inputs.culling_limit = args.culling_limit;
    inputs.best_hit_overhang = args.best_hit_overhang;
    inputs.best_hit_score_edge = args.best_hit_score_edge;
    inputs.subject_besthit = args.subject_besthit;
    // `-threshold` is read as a double and recorded as an `int` (`static_cast<int>`), converted
    // as the x86-64 build does (`INT_MIN` past `Int4`, as the lookup table takes it:
    // `ncbi_int4_from_double`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:884-892
    inputs.threshold = args.threshold.map(ncbi_int4_from_double);
    inputs.window_size = args.window_size;
    inputs.ungapped = args.ungapped;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
    // `MaskAtHash` is recorded with the value given (the class always sets it: false when absent).
    inputs.soft_masking = args.soft_masking;
    inputs.seg = args
        .seg
        .as_deref()
        .and_then(crate::report::archive::options::seg_filter);
    inputs.hitlist_size = hitlist_size;
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
    options: &crate::report::archive::OptionsHandle,
) -> Option<ArchiveHandle> {
    if !archive_requested {
        return None;
    }
    // NCBI keeps the name without the white space at its ends (`NStr::TruncateSpaces`,
    // blast_args.cpp:2472-2476), and the request's `database` is that name.
    let database = args.db.as_deref().map(crate::seqdb::db_argument_name);
    Some(ArchiveHandle::for_program(
        &option_inputs(args),
        database,
        options,
    ))
}

/// The requests of a run's batches (`RunRequests`): protein queries and subjects (the `-subject`
/// records `subjects`, `None` for `-db`), and with `-lcase_masking` (`lcase_masking`) one `notset`
/// mask of each query's lower-case runs (SEG is never one); `pssm` is psiblast's `-in_pssm` as
/// read, which every request holds in place of the queries (`RunRequests::with_pssm_query`,
/// `PsiSearch::request`).
pub(crate) fn requests(
    subjects: Option<&[FastaRecord]>,
    lcase_masking: bool,
    pssm: Option<std::sync::Arc<crate::report::serial::pssm::PssmWithParameters>>,
) -> RunRequests<'_> {
    RunRequests::new(
        RequestPolicy {
            query_protein: true,
            subject_protein: true,
            lcase: lcase_masking.then_some(LcaseFrames::NotSet),
        },
        subjects,
    )
    .with_pssm_query(pssm)
}

pub(crate) use crate::report::archive::error_lines;

/// The queries of a call's batches and what their archives read (`for_each_batch`).
pub(crate) struct CallBatches<'a> {
    /// The searched queries of the call (their records keep the whole letters of a `-query_loc`
    /// cut).
    pub queries: &'a [FastaRecord],
    pub query_placements: &'a Placements,
    /// `-subject`: every subject record, whole (the Bioseqs of the request).
    pub subject_records: Option<&'a [FastaRecord]>,
    /// The SEG masks and, with `-lcase_masking`, the lower-case runs of each query, in record
    /// coordinates (`query_record_masks`).
    pub masks: &'a [Vec<(i32, i32)>],
    /// `-lcase_masking`.
    pub lcase_masking: bool,
    /// Whether the query's batch was not searched (every query of it invalid).
    pub unsearched: &'a [bool],
    /// The queries of each NCBI batch of the call, in order.
    pub spans: &'a [Range<usize>],
}

/// What the archive of a call needs: the call's batches (the batches of one block, or of the
/// whole search) and their results.
pub(crate) struct ArchiveCall<'a> {
    pub handle: &'a ArchiveHandle,
    /// `!args.ungapped`: the archive's Seq-aligns are those of outfmt 8 (`write_blastp_seqalign`).
    pub gapped: bool,
    /// The search is rpsblast's (`RunRecord::rps`).
    pub rps: bool,
    pub batches: CallBatches<'a>,
    pub query_lengths: &'a [usize],
    /// The subject names of the Seq-aligns.
    pub subjects: SubjectNames<'a>,
    pub subject_placements: &'a Placements,
    pub hit_lists: &'a [Option<BlastpHitList>],
    pub engine_counts: &'a EngineCounts,
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
    /// That context's query length, 0 for an invalid one (`archive::RequestLength::space`).
    pub context0_length: i64,
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

/// Calls `each` with the queries and the archive of each batch of `call`, in order
/// (`BatchArchive`: its request, `RunRequests::with_batch`, and its results).
/// `first_context(first_query)` gives the values of a batch's first query.
pub(crate) fn for_each_batch(
    call: &CallBatches<'_>,
    first_context: &dyn Fn(usize) -> FirstContext,
    mut each: impl FnMut(Range<usize>, &BatchArchive<'_>) -> Result<()>,
) -> Result<()> {
    let whole = [0..call.queries.len()];
    let spans: &[Range<usize>] = if call.spans.is_empty() {
        &whole
    } else {
        call.spans
    };
    // psiblast's `-in_pssm` writes no archive (`in_pssm::check_in_pssm_options`): the requests'
    // queries are the call's records.
    let requests = requests(call.subject_records, call.lcase_masking, None);
    for span in spans.iter().cloned() {
        if span.is_empty() {
            continue;
        }
        let first = span.start;
        let searched = !call.unsearched[first];
        let context = first_context(first);

        // The results' masks are those of the search, so only a searched batch has them: one
        // `notset` mask per query with the SEG intervals and lower-case runs (a protein query).
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
        requests.with_batch(
            &call.queries[span.clone()],
            call.query_placements.get(first),
            |request| {
                let batch = BatchArchive {
                    request: *request,
                    masks: &masks,
                    ka_blocks,
                    search_space: context.search_space,
                    length_adjustment: context.length_adjustment,
                    searched_space: searched.then_some(crate::report::archive::SearchedContext0 {
                        space: context.context0_space,
                        query_length: context.context0_length,
                    }),
                };
                each(span, &batch)
            },
        )?;
    }
    Ok(())
}

/// Writes the archive of each batch of `call` to `writer`; `warnings` writes the warnings
/// of a batch's queries (to the diagnostics, which saves them as the batch's messages)
/// before its archive, and its sink gets NCBI's posts of the bytes that the archive's stream
/// replaces (`write_archive`). `first_context(first_query)` gives the values of a batch's first
/// query.
pub(crate) fn write_blastp_archives(
    call: &ArchiveCall<'_>,
    warnings: &mut QueryWarnings<'_>,
    first_context: &dyn Fn(usize) -> FirstContext,
    writer: &mut dyn Write,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Blastp,
        rps: call.rps,
        db_title: call.subjects.db_title(),
        // An ungapped search's Seq-aligns are Dense-diags, as in outfmt 8 (`seqalign.rs`).
        gapped: call.gapped,
    };
    for_each_batch(&call.batches, first_context, |span, batch| {
        for index in span.clone() {
            warnings.before_query(index, &mut *writer)?;
        }
        let mut aligns = BlastpAligns {
            run,
            queries: QueryNames {
                records: call.batches.queries,
                searched_lengths: call.query_lengths,
                placements: call.batches.query_placements,
            },
            subjects: call.subjects,
            subject_placements: call.subject_placements,
            hit_lists: &call.hit_lists[span.clone()],
            engine_counts: call.engine_counts,
            first: span.start,
        };
        write_archive(
            call.handle,
            batch,
            &mut aligns,
            &mut *writer,
            &mut *warnings.sink,
        )
    })
}

#[cfg(test)]
#[path = "archive_tests.rs"]
mod tests;
