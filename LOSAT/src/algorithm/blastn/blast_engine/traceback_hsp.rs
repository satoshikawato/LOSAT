//! The parts of the blastn traceback loop of one subject that depend on one HSP alone: where its
//! traceback starts (`prepare_traceback`), its traceback call as a function of its index
//! (`x_trace_compute`, for the look-ahead helpers and the ordered loop) and the speculative batch,
//! which makes the calls of the next HSPs on the threads of the pool ahead of the ordered loop
//! (`speculative_batch`, in the builds with the parallel feature outside serial Wasm).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375-514 (the loop of
//! `Blast_TracebackFromHSPList` up to the traceback call: `BlastGetOffsetsForGappedAlignment` or
//! `BlastGetStartForGappedAlignmentNucl`, `AdjustSubjectRange`, then
//! `BLAST_GreedyGappedAlignment` or `BLAST_GappedAlignmentWithTraceback`).
//!
//! This module is a child of `run` for the types it keeps private (`PrelimHit`, `XAheadScratch`,
//! `XSpecTraceback`). The ordered loop (the tree test, the use of a value computed ahead and the
//! tree insertion, in NCBI order) stays in `x_trace_list` of `search_query_batch`'s
//! `process_subject`. What is computed here depends on the HSP alone, so a value computed ahead on
//! any thread is the one the ordered loop would compute. `x_trace_compute` and the batch's jobs run
//! on other threads of the pool, so what they read (`TracebackCtx`) is `Sync`, and each pairs its
//! fence mark and note (`x_fence`, counted per thread) on the thread that computes.

use super::traceback_finish::TracebackCtx;
use super::{
    adjust_subject_range, blast_get_offsets_for_gapped_alignment,
    blast_get_start_for_gapped_alignment_nucl, extend_gapped_heuristic_with_traceback_with_scratch,
    greedy_gapped_alignment_with_traceback, x_fence, PrelimHit, XAheadScratch, XSpecTraceback,
};
use crate::common::GapEditOp;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use {
    super::{BlastIntervalTree, GapAlignScratch, GreedyAlignScratch, QueryContext, TreeHsp},
    crate::utils::threading::SearchPool,
    rayon::prelude::*,
};

/// Where the traceback of one HSP starts: `(trace_q_start, trace_s_start, start_shift,
/// adjusted_s_len, trace_s_start_adj)`, the start offsets on the query and on the subject, the
/// subject range that `AdjustSubjectRange` keeps (its first letter and its length) and the
/// subject start in that range.
pub(super) type TraceStart = (usize, usize, usize, usize, usize);

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

/// What the greedy traceback call returns (`greedy_gapped_alignment_with_traceback`).
pub(super) type GreedyTraceback = Option<(usize, usize, usize, usize, i32, usize, Vec<GapEditOp>)>;

/// The traceback of one HSP computed ahead of the ordered loop: where it starts and the value of
/// its traceback call.
pub(super) type XTraced = (TraceStart, XSpecTraceback<DpTraceback, GreedyTraceback>);

/// Where the traceback of `prelim` starts: its start offsets (`BlastGetOffsetsForGappedAlignment`
/// for an HSP without a gapped start, `BlastGetStartForGappedAlignmentNucl` otherwise) and the
/// subject range of its alignment (`AdjustSubjectRange`). `None` when no start is found (NCBI frees
/// the HSP and goes to the next one).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-447,457,463-473
/// ```c
///          if (!kIsOutOfFrame && hsp->query.gapped_start == 0 &&
///                                hsp->subject.gapped_start == 0) {
///             Boolean retval =
///                BlastGetOffsetsForGappedAlignment(query, subject, sbp,
///                    hsp, &q_start, &s_start);
///             if (!retval)
///             {  /* Unable to find start for this HSP */
///                hsp_array[index] = Blast_HSPFree(hsp);
///                continue;
///             }
///             hsp->query.gapped_start = q_start;
///             hsp->subject.gapped_start = s_start;
/// ...
///                BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
/// ...
///          adjusted_s_length = subject_length;
///          adjusted_subject = subject;
///
///          if (!kTranslateSubject && !kSmithWaterman) {
///              AdjustSubjectRange(&s_start, &adjusted_s_length, q_start,
///                                 query_length, &start_shift);
///              adjusted_subject = subject + start_shift;
///              /* Shift the gapped start in HSP structure, to compensate for
///                 a shift in the other direction later. */
///              hsp->subject.gapped_start = s_start;
///          }
/// ```
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405,436-472,583-612
/// ```c
/// if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, hit_options->min_diag_separation)) {
///     BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
///     AdjustSubjectRange(&s_start, &adjusted_s_length, q_start, query_length, &start_shift);
///     /* traceback, identity test, then BlastIntervalTreeAddHSP */
/// }
/// ```
/// N02 scheduling only: speculative values have no externally visible effects;
/// the original ordered contains/materialize/add sequence remains authoritative.
pub(super) fn prepare_traceback(ctx: &TracebackCtx<'_>, prelim: &PrelimHit) -> Option<TraceStart> {
    let &TracebackCtx {
        s_seq_blastna,
        encoded_queries_blastna,
        score_matrix,
        ..
    } = ctx;
    let q_seq_blastna = encoded_queries_blastna[prelim.context_idx as usize];
    let mut trace_q_start = prelim.seed_qs;
    let mut trace_s_start = prelim.seed_ss;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:436-460
    // ```c
    // if (!kIsOutOfFrame && hsp->query.gapped_start == 0 &&
    //                       hsp->subject.gapped_start == 0) {
    //    Boolean retval =
    //       BlastGetOffsetsForGappedAlignment(query, subject, sbp,
    //           hsp, &q_start, &s_start);
    //    if (!retval) { ... }
    //    hsp->query.gapped_start = q_start;
    //    hsp->subject.gapped_start = s_start;
    // } else {
    //    ...
    //    BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
    //    q_start = hsp->query.gapped_start;
    //    s_start = hsp->subject.gapped_start;
    // }
    // ```
    if trace_q_start == 0 && trace_s_start == 0 {
        let (q_start, s_start) = blast_get_offsets_for_gapped_alignment(
            q_seq_blastna,
            s_seq_blastna,
            prelim.prelim_qs,
            prelim.prelim_qe,
            prelim.prelim_ss,
            prelim.prelim_se,
            score_matrix,
        )?;
        trace_q_start = q_start;
        trace_s_start = s_start;
    } else {
        let (q_start, s_start) = blast_get_start_for_gapped_alignment_nucl(
            q_seq_blastna,
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:466-472
    // ```c
    // AdjustSubjectRange(&s_start, &adjusted_s_length, q_start,
    //                    query_length, &start_shift);
    // adjusted_subject = subject + start_shift;
    // hsp->subject.gapped_start = s_start;
    // ```
    let mut adjusted_s_len = s_seq_blastna.len();
    let mut s_start_i32 = trace_s_start as i32;
    let mut s_len_i32 = adjusted_s_len as i32;
    let start_shift_i32 = adjust_subject_range(
        &mut s_start_i32,
        &mut s_len_i32,
        trace_q_start as i32,
        q_seq_blastna.len() as i32,
    );
    let start_shift = start_shift_i32 as usize;
    adjusted_s_len = s_len_i32 as usize;
    let trace_s_start_adj = s_start_i32 as usize;

    Some((
        trace_q_start,
        trace_s_start,
        start_shift,
        adjusted_s_len,
        trace_s_start_adj,
    ))
}

/// The traceback of `prelim_hits[x_index]`: where it starts and the value of its traceback call on
/// the scratch `x_scratch`, or `None` when it has no start or the subject met a fence
/// (LOSAT_X_BNTBRANGE). The look-ahead helpers call it on the other threads of the pool, the
/// ordered loop for an HSP no helper has taken.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-440,457,466-469,503-509
/// ```c
///          if (!kIsOutOfFrame && hsp->query.gapped_start == 0 &&
///                                hsp->subject.gapped_start == 0) {
///             Boolean retval =
///                BlastGetOffsetsForGappedAlignment(query, subject, sbp,
/// ...
///                BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
/// ...
///             AdjustSubjectRange(&s_start, &adjusted_s_length, q_start,
///                                query_length, &start_shift);
/// ...
///         } else if (kGreedyTraceback) {
///             BLAST_GreedyGappedAlignment(query, adjusted_subject,
/// ...
///         } else {
///           BLAST_GappedAlignmentWithTraceback(program_number, query,
/// ```
/// This function is the part of the loop body that depends on the HSP alone: the start offsets and
/// the subject range (`prepare_traceback`) and the traceback call, dynamic-programming or greedy.
/// Index in, value out. The tree test and the tree insertion stay in the ordered loop.
/// No NCBI counterpart: the traceback of the ordered loop (`x_trace_list` in `process_subject`)
/// as a function of the HSP's index (the job of the speculative batch).
pub(super) fn x_trace_compute(
    ctx: &TracebackCtx<'_>,
    prelim_hits: &[PrelimHit],
    x_index: usize,
    x_scratch: &mut XAheadScratch<'_>,
) -> Option<XTraced> {
    let &TracebackCtx {
        s_seq_blastna,
        encoded_queries_blastna,
        reward,
        penalty,
        gap_open,
        gap_extend,
        score_matrix,
        use_dp,
        x_fenced,
        x_drop_final,
        ..
    } = ctx;
    let p = &prelim_hits[x_index];
    let prepared = prepare_traceback(ctx, p)?;
    let (qs, _, shift, slen, ss) = prepared;
    if x_fence::x_fenced(x_fenced) {
        return None;
    }
    let x_mark = x_fence::x_fence_mark(x_fenced);
    if !use_dp {
        let result = greedy_gapped_alignment_with_traceback(
            encoded_queries_blastna[p.context_idx as usize],
            &s_seq_blastna[shift..shift + slen],
            slen,
            qs,
            ss,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop_final,
            &mut *x_scratch.greedy,
        );
        x_fence::x_fence_note(x_fenced, x_mark);
        return Some((prepared, XSpecTraceback::Greedy(result)));
    }
    let result = extend_gapped_heuristic_with_traceback_with_scratch(
        encoded_queries_blastna[p.context_idx as usize],
        &s_seq_blastna[shift..shift + slen],
        qs,
        ss,
        1,
        reward,
        penalty,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop_final,
        &mut *x_scratch.gap,
    );
    x_fence::x_fence_note(x_fenced, x_mark);
    Some((prepared, XSpecTraceback::Dp(result)))
}

/// One speculative batch of the traceback loop: the loop's place, what the batch reads, and the
/// storage it reuses from batch to batch (the jobs, the results by their place in the batch, and
/// one scratch slot per thread with its results). The storage stays with `x_trace_list`, whose
/// ordered loop takes the results.
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
    /// The tree as the ordered loop has it now (read only), and the queries' contexts.
    pub(super) interval_tree: &'b BlastIntervalTree,
    pub(super) query_contexts: &'b [QueryContext],
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,403-405,509-512
/// ```c
///    for (index=0; index < num_initial_hsps; index++) {
/// ...
///       if (program_number == eBlastTypeRpsBlast ||
///           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
///                                hit_options->min_diag_separation)) {
/// ...
///            BLAST_GappedAlignmentWithTraceback(program_number, query,
///                  adjusted_subject, gap_align, score_params, q_start, s_start,
///                  query_length, adjusted_s_length,
///                  fence_hit);
/// ```
/// No NCBI counterpart: the batch itself, scheduling only (NCBI makes each traceback call in its
/// loop).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn speculative_batch(ctx: &TracebackCtx<'_>, batch: SpeculativeBatch<'_, '_>) {
    let &TracebackCtx {
        s_seq_blastna,
        encoded_queries_blastna,
        reward,
        penalty,
        gap_open,
        gap_extend,
        score_matrix,
        use_dp,
        min_diag_separation,
        x_fenced,
        x_drop_final,
        ..
    } = ctx;
    let SpeculativeBatch {
        prelim_hits,
        prelim_index,
        batch_size: SPECULATIVE_TRACEBACK_BATCH_SIZE,
        interval_tree,
        query_contexts,
        num_threads,
        parallel_pool,
        speculative_jobs,
        speculative_results,
        speculative_scratch,
    } = batch;
    let batch_end = (prelim_index + SPECULATIVE_TRACEBACK_BATCH_SIZE).min(prelim_hits.len());
    speculative_jobs.clear();
    speculative_jobs.extend((prelim_index..batch_end).filter(|&index| {
        let p = &prelim_hits[index];
        let h = TreeHsp {
            query_offset: p.prelim_qs as i32,
            query_end: p.prelim_qe as i32,
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
    // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:509-513
    // BLAST_GappedAlignmentWithTraceback(...);
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
                            // NCBI reference: core/blast_traceback.c:436-472,509-513
                            // BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                            // AdjustSubjectRange(&s_start, &adjusted_s_length,
                            //                    q_start, query_length, &start_shift);
                            // BLAST_GappedAlignmentWithTraceback(...);
                            // Keep the exact preparation with its speculative DP.
                            let prepared = prepare_traceback(ctx, p)?;
                            let (qs, _, shift, slen, ss) = prepared;
                            if x_fence::x_fenced(x_fenced) {
                                return None;
                            }
                            let x_mark = x_fence::x_fence_mark(x_fenced);
                            if !use_dp {
                                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-507
                                // ```c
                                //         } else if (kGreedyTraceback) {
                                //             BLAST_GreedyGappedAlignment(query, adjusted_subject,
                                //                 query_length, adjusted_s_length, gap_align,
                                //                 score_params, q_start, s_start, FALSE, TRUE,
                                //                 fence_hit);
                                // ```
                                // The greedy traceback call of the ordered loop, made ahead on this slot's
                                // scratch. The call depends on its arguments only, so the value is the one the ordered loop would
                                // compute.
                                // No NCBI counterpart: the call of the ordered loop
                                // (`x_trace_list` in `process_subject`), on this slot's scratch.
                                let result = greedy_gapped_alignment_with_traceback(
                                    encoded_queries_blastna[p.context_idx as usize],
                                    &s_seq_blastna[shift..shift + slen],
                                    slen,
                                    qs,
                                    ss,
                                    reward,
                                    penalty,
                                    gap_open,
                                    gap_extend,
                                    x_drop_final,
                                    x_greedy_scratch,
                                );
                                x_fence::x_fence_note(x_fenced, x_mark);
                                return Some((index, (prepared, XSpecTraceback::Greedy(result))));
                            }
                            let result = extend_gapped_heuristic_with_traceback_with_scratch(
                                encoded_queries_blastna[p.context_idx as usize],
                                &s_seq_blastna[shift..shift + slen],
                                qs,
                                ss,
                                1,
                                reward,
                                penalty,
                                score_matrix,
                                gap_open,
                                gap_extend,
                                x_drop_final,
                                scratch,
                            );
                            x_fence::x_fence_note(x_fenced, x_mark);
                            Some((index, (prepared, XSpecTraceback::Dp(result))))
                        }),
                );
            },
        )
    });
    // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:583-612
    // Blast_HSPUpdateWithTraceback(gap_align, hsp);
    // BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
    // One owner restores each result to its original batch index.
    for (_, _, results) in speculative_scratch.iter_mut() {
        for (index, result) in results.drain(..) {
            speculative_results[index - prelim_index] = Some(result);
        }
    }
}
