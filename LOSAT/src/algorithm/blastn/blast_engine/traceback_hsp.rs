//! The parts of the blastn traceback loop of one subject that work on one HSP. Those that depend on
//! the HSP alone: where its traceback starts (`prepare_traceback`), its traceback call as a
//! function of its index (`x_trace_compute`, for the look-ahead helpers and the ordered loop) and
//! the speculative batch, which makes the calls of the next HSPs on the threads of the pool ahead
//! of the ordered loop (`speculative_batch`, in the builds with the parallel feature outside serial
//! Wasm). Those that the ordered loop calls in turn for an HSP the tree does not contain: its
//! alignment (`trace_align_one`: the traceback call, or its value computed ahead) and the HSP it
//! makes, which goes into the tree and to the end of the list (`trace_build_hsp`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375-605 (the loop of
//! `Blast_TracebackFromHSPList`: `BlastGetOffsetsForGappedAlignment` or
//! `BlastGetStartForGappedAlignmentNucl`, `AdjustSubjectRange`, `BLAST_GreedyGappedAlignment` or
//! `BLAST_GappedAlignmentWithTraceback`, then `Blast_HSPUpdateWithTraceback`, `Blast_HSPTest`,
//! `Blast_HSPAdjustSubjectOffset` and `BlastIntervalTreeAddHSP`).
//!
//! This module is a child of `run` for the types it keeps private (`PrelimHit`, `QueryContext`,
//! `XAheadScratch`, `XSpecTraceback`, `BlastnTiming`). The ordered loop (the tree test, the use of
//! a value computed ahead and the calls of `trace_align_one` and `trace_build_hsp`, in NCBI order)
//! is `traceback_list::trace_hsp_list`. What the first three compute depends on the HSP alone, so
//! a value computed ahead on any thread is the one the ordered loop would compute.
//! `x_trace_compute` and the batch's jobs run on other threads of the pool, so what they read
//! (`TracebackCtx`) is `Sync`, and each pairs its fence mark and note (`x_fence`, counted per
//! thread) on the thread that computes, as `trace_align_one` does on the thread of the ordered
//! loop.

use super::traceback_finish::TracebackCtx;
use super::{
    adjust_blastn_offsets, adjust_subject_range, blast_get_offsets_for_gapped_alignment,
    blast_get_start_for_gapped_alignment_nucl, blastn_trace, calculate_blastn_context_statistics,
    context_query, extend_gapped_heuristic_with_traceback_with_scratch,
    format_gap_edit_ops_for_trace, greedy_gapped_alignment_with_traceback, stats_from_edit_ops,
    x_fence, AdjustedSubject, BlastIntervalTree, BlastnHsp, BlastnTiming, GapAlignScratch,
    GreedyAlignScratch, IndexMethod, PrelimHit, QueryContext, TreeHsp, XAheadScratch,
    XSpecTraceback,
};
use crate::common::GapEditOp;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use {crate::utils::threading::SearchPool, rayon::prelude::*};

/// Where the traceback of one HSP starts: `(trace_q_start, trace_s_start, start_shift,
/// adjusted_s_len, trace_s_start_adj)`, the start offsets on the query (NCBI's `Int4` one in the
/// context, outside it for a greedy HSP from a hit that passes an end of its context) and on the
/// subject, the subject range that `AdjustSubjectRange` keeps (its first letter and its length)
/// and the subject start in that range, NCBI's `Int4` values: on a subject of 90000 letters or
/// more a start more than 3000 letters past the end of its context gives a length below the start
/// in the range, and one more than 3000 letters before its context a start below 0 in the range
/// and a shift that can pass the subject's end, with a length below 0 then
/// (blast_gapalign.c:4163-4191).
pub(super) type TraceStart = (i32, usize, i32, i32, i32);

/// What the dynamic-programming traceback call returns
/// (`extend_gapped_heuristic_with_traceback_with_scratch`).
pub(super) type DpTraceback = (
    usize,
    usize,
    usize,
    usize,
    i32,
    usize,
    usize,
    usize,
    usize,
    Vec<GapEditOp>,
);

/// What the greedy traceback call returns (`greedy_gapped_alignment_with_traceback`; its offsets
/// are NCBI's `Int4` ones in the context and in the subject range).
pub(super) type GreedyTraceback = Option<(i32, i32, i32, i32, i32, usize, Vec<GapEditOp>)>;

/// The traceback of one HSP computed ahead of the ordered loop: where it starts and the value of
/// its traceback call.
pub(super) type XTraced = (TraceStart, XSpecTraceback<DpTraceback, GreedyTraceback>);

/// Where the traceback of `prelim` starts: its start offsets (`BlastGetOffsetsForGappedAlignment`
/// for an HSP without a gapped start, `BlastGetStartForGappedAlignmentNucl` otherwise) and the
/// subject range of its alignment (`AdjustSubjectRange`). `None` when no start is found (NCBI frees
/// the HSP and goes to the next one).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-447,457-460,463-473
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:403-405,436-472,583-612
/// N02 scheduling only: speculative values have no externally visible effects;
/// the original ordered contains/materialize/add sequence remains authoritative.
pub(super) fn prepare_traceback(
    traceback_ctx: &TracebackCtx<'_>,
    prelim: &PrelimHit,
) -> Option<TraceStart> {
    let &TracebackCtx {
        s_seq_blastna,
        encoded_query_concat_blastna,
        query_contexts,
        score_matrix,
        ..
    } = traceback_ctx;
    // NCBI's `query` (c++/src/algo/blast/core/blast_traceback.c:386-387): the context in the concatenated query, read at
    // the HSP's offsets, which lie outside it for a greedy HSP from a hit that passes an end of
    // its context. Both functions below give the same offsets on the context's letters alone for
    // an HSP in its context.
    let query = context_query(
        encoded_query_concat_blastna,
        &query_contexts[prelim.context_idx as usize],
    );
    let mut trace_q_start = prelim.seed_qs;
    let mut trace_s_start = prelim.seed_ss;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-460
    if trace_q_start == 0 && trace_s_start == 0 {
        let (q_start, s_start) = blast_get_offsets_for_gapped_alignment(
            query.buf,
            s_seq_blastna,
            query.at(prelim.prelim_qs),
            query.at(prelim.prelim_qe),
            prelim.prelim_ss,
            prelim.prelim_se,
            score_matrix,
        )?;
        trace_q_start = (q_start as isize - query.start as isize) as i32;
        trace_s_start = s_start;
    } else {
        let (q_start, s_start) = blast_get_start_for_gapped_alignment_nucl(
            query,
            s_seq_blastna,
            prelim.prelim_qs,
            prelim.prelim_qe,
            prelim.prelim_ss,
            prelim.prelim_se,
            trace_q_start,
            trace_s_start,
        );
        trace_q_start = q_start;
        trace_s_start = s_start;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:466-472
    let mut trace_s_start_adj = trace_s_start as i32;
    let mut adjusted_s_len = s_seq_blastna.len() as i32;
    let start_shift = adjust_subject_range(
        &mut trace_s_start_adj,
        &mut adjusted_s_len,
        trace_q_start,
        query.len,
    );

    Some((
        trace_q_start,
        trace_s_start,
        start_shift,
        adjusted_s_len,
        trace_s_start_adj,
    ))
}

/// NCBI's `adjusted_subject` and `adjusted_s_length` of the traceback that starts at `start`, on
/// the subject `s_seq` (`AdjustedSubject`): the whole subject, the pointer at `start_shift`, which
/// is 0 or `s_offset - max_extension_left > 0` (never below 0, past the subject's end for a start
/// far before its context), and the signed length.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:466-469
#[inline]
pub(super) fn adjusted_subject<'a>(s_seq: &'a [u8], start: &TraceStart) -> AdjustedSubject<'a> {
    let &(_, _, start_shift, adjusted_s_len, _) = start;
    AdjustedSubject {
        buf: s_seq,
        start: start_shift as usize,
        len: adjusted_s_len,
    }
}

/// The greedy traceback call of `prelim` from where `prepare_traceback` starts it (`start`), on
/// the scratch `scratch`: NCBI's `query` is the context in the concatenated query
/// (`ContextQuery`), which the alignment reads outside the context from a start outside it, and
/// NCBI's `adjusted_subject` the range of the whole subject (`adjusted_subject`), which it reads
/// outside the range from a start outside it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-507
#[inline]
pub(super) fn greedy_trace_call(
    traceback_ctx: &TracebackCtx<'_>,
    prelim: &PrelimHit,
    start: &TraceStart,
    scratch: &mut GreedyAlignScratch,
) -> GreedyTraceback {
    let &TracebackCtx {
        s_seq_blastna,
        encoded_query_concat_blastna,
        query_contexts,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop_final,
        ..
    } = traceback_ctx;
    let &(q_start, _, _, _, s_start) = start;
    greedy_gapped_alignment_with_traceback(
        context_query(
            encoded_query_concat_blastna,
            &query_contexts[prelim.context_idx as usize],
        ),
        adjusted_subject(s_seq_blastna, start),
        q_start,
        s_start,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop_final,
        scratch,
    )
}

/// The dynamic-programming traceback call of `prelim` from where `prepare_traceback` starts it
/// (`start`), on the scratch `scratch`, on the context's letters: a dynamic-programming HSP lies
/// in its context (its preliminary alignment runs on `query_tmp` from a word hit in it), so the
/// start is not negative, and neither is any value of its subject range: `max_extension_left`
/// and `max_extension_right` are 3000 or more, so the range lies in the subject and holds the
/// start (blast_gapalign.c:4163-4191).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:508-512
#[inline]
pub(super) fn dp_trace_call(
    traceback_ctx: &TracebackCtx<'_>,
    prelim: &PrelimHit,
    start: &TraceStart,
    scratch: &mut GapAlignScratch,
) -> DpTraceback {
    let &TracebackCtx {
        s_seq_blastna,
        encoded_queries_blastna,
        reward,
        penalty,
        gap_open,
        gap_extend,
        score_matrix,
        x_drop_final,
        ..
    } = traceback_ctx;
    let &(q_start, _, start_shift, adjusted_s_len, s_start) = start;
    let in_context =
        |v: i32| usize::try_from(v).expect("a dynamic-programming traceback starts in its context");
    let (q_start, start_shift, adjusted_s_len, s_start) = (
        in_context(q_start),
        in_context(start_shift),
        in_context(adjusted_s_len),
        in_context(s_start),
    );
    extend_gapped_heuristic_with_traceback_with_scratch(
        encoded_queries_blastna[prelim.context_idx as usize],
        &s_seq_blastna[start_shift..start_shift + adjusted_s_len],
        q_start,
        s_start,
        1,
        reward,
        penalty,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop_final,
        scratch,
    )
}

/// The traceback of `prelim_hits[x_index]`: where it starts and the value of its traceback call on
/// the scratch `x_scratch`, or `None` when it has no start or the subject met a fence
/// (LOSAT_X_BNTBRANGE). The look-ahead helpers call it on the other threads of the pool, the
/// ordered loop for an HSP no helper has taken.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-440,457,466-469,503-509
/// This function is the part of the loop body that depends on the HSP alone: the start offsets and
/// the subject range (`prepare_traceback`) and the traceback call, dynamic-programming or greedy.
/// Index in, value out. The tree test and the tree insertion stay in the ordered loop.
/// No NCBI counterpart: the traceback of the ordered loop (`trace_align_one`, called by
/// `traceback_list::trace_hsp_list`) as a function of the HSP's index, which the look-ahead
/// helpers run and the ordered loop runs through `Ahead::take` for an HSP that no helper has
/// taken.
pub(super) fn x_trace_compute(
    traceback_ctx: &TracebackCtx<'_>,
    prelim_hits: &[PrelimHit],
    x_index: usize,
    x_scratch: &mut XAheadScratch<'_>,
) -> Option<XTraced> {
    let &TracebackCtx {
        use_dp, x_fenced, ..
    } = traceback_ctx;
    let p = &prelim_hits[x_index];
    let prepared = prepare_traceback(traceback_ctx, p)?;
    if x_fence::x_fenced(x_fenced) {
        return None;
    }
    let x_mark = x_fence::x_fence_mark(x_fenced);
    if !use_dp {
        let result = greedy_trace_call(traceback_ctx, p, &prepared, &mut *x_scratch.greedy);
        x_fence::x_fence_note(x_fenced, x_mark);
        return Some((prepared, XSpecTraceback::Greedy(result)));
    }
    let result = dp_trace_call(traceback_ctx, p, &prepared, &mut *x_scratch.gap);
    x_fence::x_fence_note(x_fenced, x_mark);
    Some((prepared, XSpecTraceback::Dp(result)))
}

/// One speculative batch of the traceback loop: the loop's place, what the batch reads, and the
/// storage it reuses from batch to batch (the jobs, the results by their place in the batch, and
/// one scratch slot per thread with its results). The storage stays with
/// `traceback_list::trace_hsp_list`, whose ordered loop takes the results.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) struct SpeculativeBatch<'b, 's> {
    /// The HSPs of the list, the index the ordered loop is at and the batch's length
    /// (LOSAT_X_SPECBATCH).
    pub(super) prelim_hits: &'b [PrelimHit],
    pub(super) prelim_index: usize,
    pub(super) batch_size: usize,
    /// The tree as the ordered loop has it now (read only).
    pub(super) interval_tree: &'b BlastIntervalTree,
    /// The threads of the pool.
    pub(super) num_threads: usize,
    pub(super) parallel_pool: &'b SearchPool<'b>,
    pub(super) speculative_jobs: &'b mut Vec<usize>,
    pub(super) speculative_results: &'b mut Vec<Option<XTraced>>,
    pub(super) speculative_scratch: &'b mut Vec<(
        &'s mut GapAlignScratch,
        &'s mut GreedyAlignScratch,
        Vec<(usize, XTraced)>,
    )>,
}

/// The speculative batch of the traceback loop at `prelim_index`: of the next `batch_size` HSPs,
/// the ones the tree does not contain now get their traceback (`prepare_traceback` and the
/// traceback call) on the threads of the pool, each slot on its own scratch, and the results go to
/// `speculative_results` by their place in the batch, where the ordered loop takes them. The
/// ordered loop still tests each HSP against the tree when it reaches it; an HSP left out here
/// (contained now) that the tree no longer contains then is traced there.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,403-405,503-507,509-512
/// No NCBI counterpart: the batch itself, scheduling only (NCBI makes each traceback call in its
/// loop).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn speculative_batch(traceback_ctx: &TracebackCtx<'_>, batch: SpeculativeBatch<'_, '_>) {
    let &TracebackCtx {
        use_dp,
        min_diag_separation,
        x_fenced,
        query_contexts,
        ..
    } = traceback_ctx;
    let SpeculativeBatch {
        prelim_hits,
        prelim_index,
        batch_size,
        interval_tree,
        num_threads,
        parallel_pool,
        speculative_jobs,
        speculative_results,
        speculative_scratch,
    } = batch;
    let batch_end = (prelim_index + batch_size).min(prelim_hits.len());
    speculative_jobs.clear();
    speculative_jobs.extend((prelim_index..batch_end).filter(|&index| {
        let p = &prelim_hits[index];
        let h = TreeHsp {
            query_offset: p.prelim_qs,
            query_end: p.prelim_qe,
            subject_offset: p.prelim_ss as i32,
            subject_end: p.prelim_se as i32,
            score: p.prelim_score,
            query_frame: p.query_frame,
            query_length: query_contexts[p.context_idx as usize].len as i32,
            query_context_offset: p.query_context_offset,
            subject_frame_sign: 1,
        };
        interval_tree
            .containing_hsp(&h, p.query_context_offset, min_diag_separation)
            .is_none()
    }));
    let jobs = &*speculative_jobs;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:509-513
    crate::utils::threading::report_stage("blastn", "dp_traceback", jobs.len(), jobs.len() > 1);
    speculative_results.clear();
    speculative_results.resize_with(batch_end - prelim_index, || None);
    let pool = parallel_pool;
    pool.install(|| {
        speculative_scratch.par_iter_mut().enumerate().for_each(
            |(slot, (scratch, x_greedy_scratch, results))| {
                results.clear();
                results.extend(
                    jobs.iter()
                        .skip(slot)
                        .step_by(num_threads)
                        .filter_map(|&index| {
                            let p = &prelim_hits[index];
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-472,509-513
                            // Keep the exact preparation with its speculative DP.
                            let prepared = prepare_traceback(traceback_ctx, p)?;
                            if x_fence::x_fenced(x_fenced) {
                                return None;
                            }
                            let x_mark = x_fence::x_fence_mark(x_fenced);
                            if !use_dp {
                                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-507
                                // The greedy traceback call of the ordered loop, made ahead on this slot's
                                // scratch. The call depends on its arguments only, so the value is the one the ordered loop would
                                // compute.
                                // No NCBI counterpart: the call of the ordered loop
                                // (`trace_align_one`, called by `traceback_list::trace_hsp_list`), on this
                                // slot's scratch.
                                let result = greedy_trace_call(
                                    traceback_ctx,
                                    p,
                                    &prepared,
                                    x_greedy_scratch,
                                );
                                x_fence::x_fence_note(x_fenced, x_mark);
                                return Some((index, (prepared, XSpecTraceback::Greedy(result))));
                            }
                            let result = dp_trace_call(traceback_ctx, p, &prepared, scratch);
                            x_fence::x_fence_note(x_fenced, x_mark);
                            Some((index, (prepared, XSpecTraceback::Dp(result))))
                        }),
                );
            },
        )
    });
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-612
    // One owner restores each result to its original batch index.
    for (_, _, results) in speculative_scratch.iter_mut() {
        for (index, result) in results.drain(..) {
            speculative_results[index - prelim_index] = Some(result);
        }
    }
}

/// One HSP of the ordered traceback loop that the tree does not contain, with where its traceback
/// starts: what its alignment (`trace_align_one`) and the HSP it makes (`trace_build_hsp`) read
/// of it. Shared references and `Copy` values, made by the loop.
pub(super) struct TraceHsp<'h> {
    /// The HSP and its query context.
    pub(super) prelim: &'h PrelimHit,
    pub(super) ctx: &'h QueryContext,
    /// The context's letters in blastna without masks, which a dynamic-programming HSP's
    /// identities are counted on.
    pub(super) q_seq_nomask_blastna: &'h [u8],
    /// Whether the trace (`blastn_trace`) follows this HSP's seed, and the subject's frame sign
    /// in the tree.
    pub(super) trace_traceback_seed: bool,
    pub(super) subject_frame_sign: i32,
    /// Where its traceback starts (`TraceStart`, as `prepare_traceback` gives it).
    pub(super) trace_q_start: i32,
    pub(super) trace_s_start: usize,
    pub(super) start_shift: i32,
    pub(super) adjusted_s_len: i32,
    pub(super) trace_s_start_adj: i32,
}

/// The alignment of one HSP (`trace_align_one`): its ends on the query (NCBI's `Int4` offsets in
/// the context) and on the subject range of the alignment (NCBI's `Int4` offsets in the range,
/// below 0 for a greedy alignment from a start before the range), its score, the identities,
/// mismatches, gap opens and gap letters the traceback call counts (0 from the greedy call), its
/// length and its edit script.
pub(super) struct Aligned {
    pub(super) final_qs: i32,
    pub(super) final_qe: i32,
    pub(super) final_ss: i32,
    pub(super) final_se: i32,
    pub(super) score: i32,
    pub(super) matches: usize,
    pub(super) mismatches: usize,
    pub(super) gaps: usize,
    pub(super) gap_letters: usize,
    pub(super) aln_len: usize,
    pub(super) edit_ops: Vec<GapEditOp>,
}

/// The alignment of one HSP of the ordered traceback loop that the tree does not contain: the
/// traceback call, dynamic-programming or greedy, from where `prepare_traceback` starts it, on the
/// loop's scratch (`gap_scratch`, `greedy_align_scratch`), or the value of the same call made
/// ahead (`precomputed`: a look-ahead helper's or the speculative batch's). `None` when it makes
/// no HSP: the subject met a fence (LOSAT_X_BNTBRANGE) or the greedy call returns none. Times the
/// alignment (LOSAT_TIMING).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-520
/// The loop calls it in NCBI order, on its own thread, so the fence mark and note of a call made
/// here are paired on that thread.
pub(super) fn trace_align_one(
    traceback_ctx: &TracebackCtx<'_>,
    trace_hsp: &TraceHsp<'_>,
    precomputed: Option<XTraced>,
    gap_scratch: &mut GapAlignScratch,
    greedy_align_scratch: &mut GreedyAlignScratch,
) -> Option<Aligned> {
    let &TracebackCtx {
        timing_ref,
        timing_enabled,
        x_fenced,
        use_dp,
        ..
    } = traceback_ctx;
    let &TraceHsp {
        prelim,
        trace_q_start,
        trace_s_start,
        start_shift,
        adjusted_s_len,
        trace_s_start_adj,
        ..
    } = trace_hsp;
    let start: TraceStart = (
        trace_q_start,
        trace_s_start,
        start_shift,
        adjusted_s_len,
        trace_s_start_adj,
    );
    if let Some(timing) = timing_ref {
        BlastnTiming::record_count(&timing.traceback_full_traceback_hsps, 1);
    }
    let alignment_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    let (
        final_qs,
        final_qe,
        final_ss,
        final_se,
        score,
        matches,
        mismatches,
        gaps,
        gap_letters,
        aln_len,
        edit_ops,
    ) = if use_dp {
        let (
            final_qs,
            final_qe,
            final_ss,
            final_se,
            score,
            matches,
            mismatches,
            gaps,
            gap_letters,
            edit_ops,
        ) = {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:403-405,436-472,583-612
            // N02 scheduling only: speculative values have no externally visible effects;
            // the original ordered contains/materialize/add sequence remains authoritative.
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:273-297,558-585
            // Endpoint replacement can invalidate batch-start containment. A newly eligible
            // HSP therefore runs the unchanged DP here, at its original sequential position.
            if x_fence::x_fenced(x_fenced) {
                return None;
            }
            let x_mark = x_fence::x_fence_mark(x_fenced);
            let x_result = match precomputed {
                Some((_, XSpecTraceback::Dp(result))) => result,
                _ => dp_trace_call(traceback_ctx, prelim, &start, gap_scratch),
            };
            // EXPERIMENT (LOSAT_X_BNTBRANGE): an alignment that met a
            // fence (here or ahead) makes no HSP; the sequence is traced
            // again whole (NCBI breaks its loop, blast_traceback.c:516-520).
            x_fence::x_fence_note(x_fenced, x_mark);
            if x_fence::x_fenced(x_fenced) {
                return None;
            }
            x_result
        };
        let aln_len = matches + mismatches + gap_letters;
        (
            final_qs as i32,
            final_qe as i32,
            final_ss as i32,
            final_se as i32,
            score,
            matches,
            mismatches,
            gaps,
            gap_letters,
            aln_len,
            edit_ops,
        )
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-507
        // The `_` arm is the ported call. The other arm uses the
        // value of the same call made ahead on another scratch.
        // No NCBI counterpart: the same call, already made
        // for this HSP on another scratch, or made here.
        if x_fence::x_fenced(x_fenced) {
            return None;
        }
        let x_mark = x_fence::x_fence_mark(x_fenced);
        let x_greedy_result = match precomputed {
            Some((_, XSpecTraceback::Greedy(result))) => result,
            _ => greedy_trace_call(traceback_ctx, prelim, &start, &mut *greedy_align_scratch),
        };
        // EXPERIMENT (LOSAT_X_BNTBRANGE): as for the DP traceback above.
        x_fence::x_fence_note(x_fenced, x_mark);
        if x_fence::x_fenced(x_fenced) {
            return None;
        }
        match x_greedy_result {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-597
            // Greedy traceback keeps score/coordinates/edit-script only here.
            // Identity statistics are recomputed later during the NCBI-style
            // reevaluation and Blast_HSPTestIdentityAndLength pass.
            Some((final_qs, final_qe, final_ss, final_se, score, aln_len, edit_ops)) => (
                final_qs, final_qe, final_ss, final_se, score, 0usize, 0usize, 0usize, 0usize,
                aln_len, edit_ops,
            ),
            None => {
                if let (Some(timing), Some(alignment_start)) = (timing_ref, alignment_start) {
                    let alignment_elapsed_ns = alignment_start.elapsed().as_nanos() as u64;
                    BlastnTiming::record_ns(&timing.traceback_alignment_ns, alignment_elapsed_ns);
                    BlastnTiming::record_ns(
                        &timing.traceback_alignment_greedy_ns,
                        alignment_elapsed_ns,
                    );
                }
                return None;
            }
        }
    };
    if let (Some(timing), Some(alignment_start)) = (timing_ref, alignment_start) {
        let alignment_elapsed_ns = alignment_start.elapsed().as_nanos() as u64;
        BlastnTiming::record_ns(&timing.traceback_alignment_ns, alignment_elapsed_ns);
        if use_dp {
            BlastnTiming::record_ns(&timing.traceback_alignment_dp_ns, alignment_elapsed_ns);
        } else {
            BlastnTiming::record_ns(&timing.traceback_alignment_greedy_ns, alignment_elapsed_ns);
        }
    }
    Some(Aligned {
        final_qs,
        final_qe,
        final_ss,
        final_se,
        score,
        matches,
        mismatches,
        gaps,
        gap_letters,
        aln_len,
        edit_ops,
    })
}

/// Where the ordered traceback loop puts what one HSP makes: the tree of the HSPs made so far, the
/// list of HSPs, and the lengths of their edit scripts and alignments (LOSAT_TIMING).
pub(super) struct TraceOut<'o> {
    pub(super) interval_tree: &'o mut BlastIntervalTree,
    pub(super) hits: &'o mut Vec<BlastnHsp>,
    pub(super) traceback_edit_script_lengths: &'o mut Vec<u32>,
    pub(super) traceback_alignment_lengths: &'o mut Vec<u32>,
}

/// The HSP that the alignment of one HSP of the ordered traceback loop makes: its subject offsets
/// on the whole subject, its identities on the unmasked query (dynamic-programming traceback), the
/// identity and length test, then its insertion in the tree, its e-value and bit score and its
/// place at the end of `hits`. `false` when the test deletes it (the tree and `hits` are left as
/// they are). Records the lengths and times the step (LOSAT_TIMING).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-589,596-605
pub(super) fn trace_build_hsp(
    traceback_ctx: &TracebackCtx<'_>,
    trace_hsp: &TraceHsp<'_>,
    aligned: Aligned,
    out: TraceOut<'_>,
) -> bool {
    let &TracebackCtx {
        timing_ref,
        timing_enabled,
        blastn_trace_enabled,
        s_idx,
        s_id,
        s_seq_blastna,
        query_lengths,
        query_eff_searchsp,
        search_karlin_ref,
        round_down_evalue_score,
        x_drop_final: x_drop_trace,
        use_dp,
        percent_identity,
        ..
    } = traceback_ctx;
    let &TraceHsp {
        prelim,
        ctx,
        q_seq_nomask_blastna,
        trace_traceback_seed,
        subject_frame_sign,
        trace_q_start,
        trace_s_start,
        start_shift,
        adjusted_s_len,
        trace_s_start_adj,
        ..
    } = trace_hsp;
    let Aligned {
        final_qs,
        final_qe,
        final_ss,
        final_se,
        score,
        matches,
        mismatches,
        gaps,
        gap_letters,
        aln_len,
        edit_ops,
    } = aligned;
    let TraceOut {
        interval_tree,
        hits,
        traceback_edit_script_lengths,
        traceback_alignment_lengths,
    } = out;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:598-600
    let hsp_build_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    // The offsets in the whole subject are not negative: the alignment reads from the start's
    // place in the subject at most `s_off` letters back and `s_avail` forward
    // (`SeqPointer::left_of`, `right_of`).
    let in_subject =
        |v: i32| usize::try_from(v).expect("a traceback alignment lies in its subject");
    let final_ss = in_subject(final_ss + start_shift);
    let final_se = in_subject(final_se + start_shift);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-597
    // Traceback scoring uses query_blk->sequence; identity/length
    // reporting uses query_blk->sequence_nomask. Without masks written
    // into the query (`mask_query_letters`) these slices alias.
    // A dynamic-programming HSP lies in its context (`dp_trace_call`).
    let (matches, mismatches, gaps, gap_letters, aln_len) = if use_dp {
        let (matches, mismatches, gaps, gap_letters) = stats_from_edit_ops(
            q_seq_nomask_blastna,
            s_seq_blastna,
            final_qs as usize,
            final_ss,
            &edit_ops,
        );
        (
            matches,
            mismatches,
            gaps,
            gap_letters,
            matches + mismatches + gap_letters,
        )
    } else {
        (matches, mismatches, gaps, gap_letters, aln_len)
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-472
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-600
    if blastn_trace_enabled
        && (trace_traceback_seed
            || blastn_trace::should_trace_range(
                "traceback",
                prelim.context_idx,
                s_idx,
                s_id,
                final_qs,
                final_qe,
                final_ss,
                final_se,
                ctx.len,
                ctx.frame,
            ))
    {
        blastn_trace::log(
            "traceback",
            format!(
                "subject={}({}) context={} start=({}, {}) adjusted_start=({}, {}) start_shift={} adjusted_s_len={} final=q{}..{} s{}..{} raw_score={} x_drop={} aln_len={} identities={} mismatches={} gaps={} gap_letters={} edit_ops_len={} edit_ops={}",
                s_id,
                s_idx,
                prelim.context_idx,
                trace_q_start,
                trace_s_start,
                trace_q_start,
                trace_s_start_adj,
                start_shift,
                adjusted_s_len,
                final_qs,
                final_qe,
                final_ss,
                final_se,
                score,
                x_drop_trace,
                aln_len,
                matches,
                mismatches,
                gaps,
                gap_letters,
                edit_ops.len(),
                format_gap_edit_ops_for_trace(&edit_ops)
            ),
        );
    }

    if timing_enabled {
        traceback_edit_script_lengths.push(edit_ops.len().min(u32::MAX as usize) as u32);
        traceback_alignment_lengths.push(aln_len.min(u32::MAX as usize) as u32);
    }

    // `Blast_HSPTest`; NCBI's `min_hit_length` is always 0 (`core::blast_hits::hsp_test`).
    if use_dp
        && crate::core::blast_hits::hsp_test(matches as i64, aln_len as i64, percent_identity, 0)
    {
        if let Some(timing) = timing_ref {
            BlastnTiming::record_count(&timing.traceback_deleted_identity_length_hsps, 1);
        }
        if let (Some(timing), Some(hsp_build_start)) = (timing_ref, hsp_build_start) {
            BlastnTiming::record_duration(&timing.traceback_hsp_build_ns, hsp_build_start);
        }
        return false;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-612
    let final_tree_hsp = TreeHsp {
        query_offset: final_qs,
        query_end: final_qe,
        subject_offset: final_ss as i32,
        subject_end: final_se as i32,
        score,
        query_frame: prelim.query_frame,
        query_length: ctx.len as i32,
        query_context_offset: prelim.query_context_offset,
        subject_frame_sign,
    };
    interval_tree.add_hsp(
        final_tree_hsp,
        prelim.query_context_offset,
        IndexMethod::QueryAndSubject,
    );

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:234-250
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1887-1890
    let eff_searchsp = query_eff_searchsp[prelim.context_idx as usize];
    let (bit_score, eval) = calculate_blastn_context_statistics(
        score,
        &search_karlin_ref[prelim.context_idx as usize].gapped,
        eff_searchsp,
        round_down_evalue_score,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:234-250
    // NCBI reaps by E-value only after common-endpoint purging,
    // ambiguity re-evaluation, score resort, and final containment.
    // Keep this traceback HSP in the subject list so it can still
    // participate in those survivor-set decisions.

    let identity = if aln_len > 0 {
        ((matches as f64 / aln_len as f64) * 100.0).min(100.0)
    } else {
        0.0
    };

    let query_length = query_lengths[prelim.query_idx as usize];
    let (hit_q_start, hit_q_end, hit_s_start, hit_s_end) = adjust_blastn_offsets(
        final_qs,
        final_qe,
        final_ss,
        final_se,
        query_length,
        prelim.query_frame,
    );

    let gap_info = if edit_ops.is_empty() {
        None
    } else {
        Some(edit_ops)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4058-4077
    hits.push(BlastnHsp {
        identity,
        length: aln_len,
        mismatch: mismatches,
        gapopen: gaps,
        q_start: hit_q_start,
        q_end: hit_q_end,
        s_start: hit_s_start,
        s_end: hit_s_end,
        e_value: eval,
        bit_score,
        num_ident: matches,
        query_frame: prelim.query_frame,
        query_length,
        q_idx: prelim.query_idx,
        s_idx: s_idx as u32,
        raw_score: score,
        internal_q_offset_0: final_qs,
        internal_q_end_0: final_qe,
        internal_s_offset_0: final_ss,
        internal_s_end_0: final_se,
        internal_query_context_offset: prelim.query_context_offset,
        gap_info,
        num_positives: matches,
    });
    if let (Some(timing), Some(hsp_build_start)) = (timing_ref, hsp_build_start) {
        BlastnTiming::record_duration(&timing.traceback_hsp_build_ns, hsp_build_start);
    }
    true
}
