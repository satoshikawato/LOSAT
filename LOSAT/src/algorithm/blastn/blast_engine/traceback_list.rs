//! The blastn traceback of one subject's kept HSPs (`traceback_subject_hits`): per query on the
//! pool (LOSAT_X_QGROUPPAR) or over the subject's list on the owner, the traceback loop, then the
//! purges and the final list (`traceback_finish`). The traceback loop takes the kept HSPs of one
//! list in turn, each tested against the interval tree, traced, and put into the tree and at the
//! end of the list of HSPs (`trace_hsp_list`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:259-723
//! (`Blast_TracebackFromHSPList`), of which :365-375,611 are the start of the traceback: the sort
//! assertion, the interval tree and the loop on the HSPs
//! ```c
//!    ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
//! ...
//!    tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
//!                                  0, (subject_length > 0 ? subject_length :
//!                                  subject_blk->length / CODON_LENGTH) + 1);
//!
//!    for (index=0; index < num_initial_hsps; index++) {
//! ...
//!    } /* End loop on HSPs */
//! ```
//!
//! This module is a child of `run` for the types it keeps private (`PrelimHit`, `XAheadScratch`,
//! `BlastnTiming`, `SubjectHits`). `trace_hsp_list` is the body of what was the closure
//! `x_trace_list` of `search_query_batch`'s `process_subject`, and `traceback_subject_hits` the
//! code of `process_subject` that called it. The values they read are fields of
//! `traceback_finish::TracebackCtx` (the search's and the subject's), of `TraceSchedule` (what
//! schedules the loop on the pool) and of `TraceDrive`; what the loop works on is
//! `TraceListState`. The steps of one HSP are in `traceback_hsp`. The query group tasks of
//! `x_qgroup_par` run `trace_hsp_list` on the threads of the pool, so what it reads is `Sync`.

use super::traceback_finish::{self, TracebackCtx};
use super::{
    blastn_trace, traceback_hsp, x_qgroup_par, BlastIntervalTree, BlastnHsp, BlastnTiming,
    FastaRecord, GapAlignScratch, GreedyAlignScratch, PrelimHit, SubjectHits, TreeHsp,
    XAheadScratch,
};
use crate::utils::threading::SearchPool;

/// What schedules the traceback loop of one list on the threads of the pool: the speculative
/// traceback batches (`speculative_traceback`; for the greedy traceback of megablast too,
/// `x_greedy_spec`, the rule of LOSAT_X_GREEDYSPEC kept as its `x_reach!` key), the window of the
/// look-ahead helpers (`x_ahead_window`, LOSAT_X_AHEAD's key), the threads and the pool. `Copy`
/// values and a shared reference, the same for every subject of the batch, so the struct is `Sync`
/// for the query group tasks.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375
/// ```c
///    for (index=0; index < num_initial_hsps; index++) {
/// ```
/// No NCBI counterpart: scheduling only. NCBI makes each traceback call in this loop on one
/// thread; here other threads of the pool may compute the call of a later HSP ahead of the loop,
/// which still tests the tree and adds to it in NCBI order.
pub(super) struct TraceSchedule<'a> {
    /// Whether the speculative batches may run (`search_query_batch`'s rule), and whether that
    /// rule holds for the greedy traceback (the `x_reach!` record of LOSAT_X_GREEDYSPEC).
    pub(super) speculative_traceback: bool,
    pub(super) x_greedy_spec: bool,
    /// The look-ahead window, `None` when the look-ahead helpers do not run.
    pub(super) x_ahead_window: Option<usize>,
    /// The threads of the pool, and the pool.
    pub(super) num_threads: usize,
    pub(super) parallel_pool: &'a SearchPool<'a>,
}

/// What the traceback loop of one list works on and fills: the interval tree, the scratch of the
/// loop's traceback calls, the scratch lists of the look-ahead helpers and of the speculative
/// batches (`None` for a query group task, which uses neither), the list the HSPs go to, and the
/// lengths of their edit scripts and alignments (LOSAT_TIMING). The owner's come from the
/// subject's scratch (`SubjectScratch`), a query group task's from its own state.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:509-512,600-601
/// ```c
///           BLAST_GappedAlignmentWithTraceback(program_number, query,
///                 adjusted_subject, gap_align, score_params, q_start, s_start,
///                 query_length, adjusted_s_length,
///                 fence_hit);
/// ...
///               status = BlastIntervalTreeAddHSP(hsp, tree, query_info,
///                                          eQueryAndSubject);
/// ```
/// NCBI's loop works on the `tree` and the `gap_align` of its call and traces the HSPs of the
/// list in place; here the HSPs it keeps are pushed to `hits`. The helpers' and the batches'
/// scratch has no NCBI counterpart (scheduling only).
pub(super) struct TraceListState<'s> {
    pub(super) interval_tree: &'s mut BlastIntervalTree,
    pub(super) gap_scratch: &'s mut GapAlignScratch,
    pub(super) greedy_align_scratch: &'s mut GreedyAlignScratch,
    pub(super) x_helper_scratch: Option<&'s mut Vec<(GapAlignScratch, GreedyAlignScratch)>>,
    pub(super) x_spec_scratch: Option<&'s mut Vec<(GapAlignScratch, GreedyAlignScratch)>>,
    pub(super) hits: &'s mut Vec<BlastnHsp>,
    pub(super) traceback_edit_script_lengths: &'s mut Vec<u32>,
    pub(super) traceback_alignment_lengths: &'s mut Vec<u32>,
}

/// The traceback loop over one list of kept HSPs: the subject's on the owner (with the look-ahead
/// helpers, whose scratch `x_helper_scratch` holds, or the speculative batches), or one query's
/// (LOSAT_X_QGROUPPAR, neither). The HSPs go to `hits`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,403-405,607-611
/// ```c
///    for (index=0; index < num_initial_hsps; index++) {
/// ...
///       if (program_number == eBlastTypeRpsBlast ||
///           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
///                                hit_options->min_diag_separation)) {
/// ...
///       } else {
///          /* Contained within another HSP, delete. */
///          hsp_array[index] = Blast_HSPFree(hsp);
///       }
///    } /* End loop on HSPs */
/// ```
/// The tree is the caller's, made for this subject's ranges (`Blast_IntervalTreeInit`,
/// blast_traceback.c:371-373); the loop empties it first. The steps of one HSP are
/// `traceback_hsp`'s.
pub(super) fn trace_hsp_list(
    traceback_ctx: &TracebackCtx<'_>,
    schedule: &TraceSchedule<'_>,
    prelim_hits: &[PrelimHit],
    state: TraceListState<'_>,
) {
    let &TracebackCtx {
        timing_ref,
        timing_enabled,
        blastn_trace_enabled,
        s_idx,
        s_id,
        encoded_queries_blastna,
        encoded_queries_nomask,
        query_contexts,
        min_diag_separation,
        ..
    } = traceback_ctx;
    let &TraceSchedule {
        speculative_traceback,
        x_greedy_spec,
        x_ahead_window,
        num_threads,
        #[cfg_attr(
            not(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            )),
            expect(
                unused_variables,
                reason = "read only by the speculative traceback batch, built with the parallel feature outside serial Wasm"
            )
        )]
        parallel_pool,
    } = schedule;
    let TraceListState {
        interval_tree,
        gap_scratch,
        greedy_align_scratch,
        x_helper_scratch,
        x_spec_scratch,
        hits,
        traceback_edit_script_lengths,
        traceback_alignment_lengths,
    } = state;
    let speculative_traceback = speculative_traceback && x_helper_scratch.is_some();
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let _ = &x_spec_scratch;
    if x_greedy_spec {
        crate::utils::xreach::x_reach!("LOSAT_X_GREEDYSPEC", speculative_traceback);
    }
    if !prelim_hits.is_empty() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:358-365
        // ```c
        // /* Make sure the HSPs in the HSP list are sorted by score, as they
        //    should be. */
        // #ifdef _DEBUG
        // {
        //     Blast_HSPListSortByScore(hsp_list);
        // }
        // #endif
        // ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
        // ```
        // The oracle is a release build: the HSPs are traced in the order of the
        // stored list (`kept`, each query's list for this subject in turn). That
        // is score order for a list kept as it came (sorted per subject chunk,
        // blast_engine.c:555, and merged by score), and e-value order for a list
        // that `Blast_HitListUpdate` sorted when the hit list became a heap
        // (blast_hits.c:3272-3284).
        if let Some(timing) = timing_ref {
            BlastnTiming::record_count(&timing.traceback_prelim_hsps, prelim_hits.len() as u64);
        }

        interval_tree.reset();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:66-70
        // ```c
        // if (tree->num_used == tree->num_alloc) {
        //     tree->num_alloc = 2 * tree->num_alloc;
        //     tree->nodes = (SIntervalNode *)realloc(tree->nodes, tree->num_alloc *
        //                                                  sizeof(SIntervalNode));
        // }
        // ```
        // Reserve from the preliminary HSP count before the hot add/contains loop.
        // This preserves NCBI node order while avoiding repeated Vec growth.
        interval_tree.reserve_nodes_for_hsps(prelim_hits.len());

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405,436-472,583-612
        // ```c
        // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, hit_options->min_diag_separation)) {
        //     BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
        //     AdjustSubjectRange(&s_start, &adjusted_s_length, q_start, query_length, &start_shift);
        //     /* traceback, identity test, then BlastIntervalTreeAddHSP */
        // }
        // ```
        // N02 scheduling only: speculative values have no externally visible effects;
        // the original ordered contains/materialize/add sequence remains authoritative.
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        let mut x_no_spec = Vec::new();
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        let x_spec_scratch = x_spec_scratch.unwrap_or(&mut x_no_spec);
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        let mut speculative_scratch: Vec<_> = if speculative_traceback {
            // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:509-513
            // BLAST_GappedAlignmentWithTraceback(...);
            // Each independent scratch slot retains its intermediate results.
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-507
            // ```c
            //         } else if (kGreedyTraceback) {
            //             BLAST_GreedyGappedAlignment(query, adjusted_subject,
            //                 query_length, adjusted_s_length, gap_align,
            // ```
            // Each slot of the batch also has the scratch of the greedy aligner, so that
            // the greedy traceback call can run on a slot instead of on the ordered loop.
            // Scratch only, made once per subject scratch (`x_spec_scratch`).
            while x_spec_scratch.len() < num_threads {
                x_spec_scratch.push((GapAlignScratch::new(), GreedyAlignScratch::new()));
            }
            x_spec_scratch
                .iter_mut()
                .take(num_threads)
                .map(|(gap, greedy)| (gap, greedy, Vec::new()))
                .collect()
        } else {
            Vec::new()
        };
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        let mut speculative_results = Vec::new();
        // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:403-405,598-601
        // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, ...)) {
        //     BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
        // }
        // Reuse only batch storage; live ordered containment remains below.
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        let mut speculative_jobs = Vec::new();
        // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:403-405,509-513,598-601
        // ```c
        // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, ...)) {
        //     BLAST_GappedAlignmentWithTraceback(...);
        //     ...
        //     BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
        // }
        // ```
        // Scheduling only; ordered containment and insertion stay unchanged.
        // Serial Wasm compiles out the speculative path and retains batch 8.
        const X_DEFAULT_SPECULATIVE_TRACEBACK_BATCH_SIZE: usize =
            if cfg!(all(target_arch = "wasm32", not(feature = "wasm-threads"))) {
                8
            } else {
                16
            };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375
        // ```c
        //    for (index=0; index < num_initial_hsps; index++) {
        // ```
        // The number of HSPs whose tracebacks are computed ahead of this loop in one batch. Scheduling
        // only: the loop still visits the HSPs in this order and decides each one itself.
        // EXPERIMENT (LOSAT_X_SPECBATCH=n): how many HSPs are traced ahead of
        // the ordered loop at a time. Scheduling only, like the batch itself.
        #[allow(non_snake_case)]
        #[cfg_attr(
            not(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            )),
            expect(
                unused_variables,
                reason = "used only by the speculative traceback batch, built with the parallel feature outside serial Wasm"
            )
        )]
        let SPECULATIVE_TRACEBACK_BATCH_SIZE: usize = {
            use std::sync::OnceLock;
            static SIZE: OnceLock<usize> = OnceLock::new();
            *SIZE.get_or_init(|| {
                std::env::var("LOSAT_X_SPECBATCH")
                    .ok()
                    .and_then(|v| v.parse::<usize>().ok())
                    .filter(|&n| (1..=4096).contains(&n))
                    .unwrap_or(X_DEFAULT_SPECULATIVE_TRACEBACK_BATCH_SIZE)
            })
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-472,503-512
        // (the part of the loop body that depends on the HSP alone). The traceback of one
        // HSP of this list by its index (`traceback_hsp::x_trace_compute`), for the
        // look-ahead helpers and the ordered loop.
        let x_trace_compute = |x_index: usize, x_scratch: &mut XAheadScratch<'_>| {
            traceback_hsp::x_trace_compute(traceback_ctx, prelim_hits, x_index, x_scratch)
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375
        // ```c
        //    for (index=0; index < num_initial_hsps; index++) {
        // ```
        // Look-ahead for the traceback loop. Without helpers `x_ahead` is `None` and
        // the loop runs as ported.
        let x_ahead = match x_ahead_window {
            Some(x_window)
                if x_helper_scratch.is_some()
                    && !blastn_trace_enabled
                    && prelim_hits.len() >= 2 * x_window =>
            {
                Some(crate::utils::xahead::Ahead::new(x_window))
            }
            _ => None,
        };
        if x_ahead_window.is_some() {
            crate::utils::xreach::x_reach!("LOSAT_X_AHEAD", x_ahead.is_some());
        }
        let mut x_no_helpers = Vec::new();
        let x_helper_scratch = x_helper_scratch.unwrap_or(&mut x_no_helpers);
        if x_ahead.is_some() {
            while x_helper_scratch.len() + 1 < num_threads {
                x_helper_scratch.push((GapAlignScratch::new(), GreedyAlignScratch::new()));
            }
        }
        let x_helpers = x_helper_scratch
            .iter_mut()
            .map(|(gap, greedy)| XAheadScratch { gap, greedy });
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,600-601
        // ```c
        //    for (index=0; index < num_initial_hsps; index++) {
        // ...
        //               status = BlastIntervalTreeAddHSP(hsp, tree, query_info,
        //                                          eQueryAndSubject);
        // ```
        // `with_helpers` runs the closure on this thread. The closure is the port of this NCBI loop, in NCBI
        // order; the other threads of the pool only evaluate tracebacks ahead.
        crate::utils::xahead::with_helpers(x_ahead.as_ref(), x_helpers, &x_trace_compute, || {
            for prelim_index in 0..prelim_hits.len() {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,404-405
                // ```c
                //    for (index=0; index < num_initial_hsps; index++) {
                // ...
                //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                //                                hit_options->min_diag_separation)) {
                // ```
                // The owner offers later HSPs to the helpers, leaving out the ones the tree already contains. The
                // tree is only read here.
                // No NCBI counterpart: offer the next HSPs to the helpers,
                // leaving out those the tree already contains.
                if let Some(x_ahead) = x_ahead.as_ref() {
                    x_ahead.announce(prelim_index, prelim_hits.len(), |x_index| {
                        let p = &prelim_hits[x_index];
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
                    });
                }
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405,436-472,583-612
                // ```c
                // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, hit_options->min_diag_separation)) {
                //     BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                //     AdjustSubjectRange(&s_start, &adjusted_s_length, q_start, query_length, &start_shift);
                //     /* traceback, identity test, then BlastIntervalTreeAddHSP */
                // }
                // ```
                // N02 scheduling only: speculative values have no externally visible effects;
                // the original ordered contains/materialize/add sequence remains authoritative.
                #[cfg(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                ))]
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,509-512
                // ```c
                //    for (index=0; index < num_initial_hsps; index++) {
                // ...
                //           BLAST_GappedAlignmentWithTraceback(program_number, query,
                //                 adjusted_subject, gap_align, score_params, q_start, s_start,
                //                 query_length, adjusted_s_length,
                //                 fence_hit);
                // ```
                // Look-ahead: with helpers the batch (`traceback_hsp::speculative_batch`) is
                // not used; the tracebacks come from `x_ahead.take`. Both ways only compute
                // the traceback call of this loop ahead of its turn.
                if speculative_traceback
                    && x_ahead.is_none()
                    && prelim_index % SPECULATIVE_TRACEBACK_BATCH_SIZE == 0
                {
                    traceback_hsp::speculative_batch(
                        traceback_ctx,
                        traceback_hsp::SpeculativeBatch {
                            prelim_hits,
                            prelim_index,
                            batch_size: SPECULATIVE_TRACEBACK_BATCH_SIZE,
                            interval_tree,
                            query_contexts,
                            num_threads,
                            parallel_pool,
                            speculative_jobs: &mut speculative_jobs,
                            speculative_results: &mut speculative_results,
                            speculative_scratch: &mut speculative_scratch,
                        },
                    );
                }

                let prelim = &prelim_hits[prelim_index];
                let ctx = &query_contexts[prelim.context_idx as usize];
                let q_seq_blastna = encoded_queries_blastna[prelim.context_idx as usize];
                let q_seq_nomask_blastna = encoded_queries_nomask[prelim.context_idx as usize];

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405
                // ```c
                // if (program_number == eBlastTypeRpsBlast ||
                //     !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                //                          hit_options->min_diag_separation)) {
                // ```
                let tree_precheck_start = if timing_enabled {
                    Some(std::time::Instant::now())
                } else {
                    None
                };
                let subject_frame_sign = 1i32;
                let prelim_tree_hsp = TreeHsp {
                    query_offset: prelim.prelim_qs as i32,
                    query_end: prelim.prelim_qe as i32,
                    subject_offset: prelim.prelim_ss as i32,
                    subject_end: prelim.prelim_se as i32,
                    score: prelim.prelim_score,
                    query_frame: prelim.query_frame,
                    query_length: ctx.len as i32,
                    query_context_offset: prelim.query_context_offset,
                    subject_frame_sign,
                };
                let traceback_containing_hsp = interval_tree.containing_hsp(
                    &prelim_tree_hsp,
                    prelim.query_context_offset,
                    min_diag_separation,
                );
                let prelim_traceback_contained = traceback_containing_hsp.is_some();
                let trace_traceback_seed = blastn_trace_enabled
                    && blastn_trace::should_trace_seed(
                        "traceback",
                        prelim.context_idx,
                        s_idx,
                        s_id,
                        prelim.seed_qs,
                        prelim.seed_ss,
                    );
                if let (Some(timing), Some(tree_precheck_start)) = (timing_ref, tree_precheck_start)
                {
                    BlastnTiming::record_duration(
                        &timing.traceback_tree_precheck_ns,
                        tree_precheck_start,
                    );
                }
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405
                // ```c
                // if (program_number == eBlastTypeRpsBlast ||
                //     !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                //                              hit_options->min_diag_separation)) {
                // ```
                if blastn_trace_enabled
                    && (trace_traceback_seed
                        || blastn_trace::should_trace_range(
                            "traceback",
                            prelim.context_idx,
                            s_idx,
                            s_id,
                            prelim.prelim_qs,
                            prelim.prelim_qe,
                            prelim.prelim_ss,
                            prelim.prelim_se,
                            ctx.len,
                            ctx.frame,
                        ))
                {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4071-4076
                    // ```c
                    // status = Blast_HSPInit(gap_align->query_start,
                    //               gap_align->query_stop, gap_align->subject_start,
                    //               gap_align->subject_stop,
                    //               init_hsp->offsets.qs_offsets.q_off,
                    //               init_hsp->offsets.qs_offsets.s_off, context,
                    //               query_frame, subject->frame, gap_align->score,
                    //               &(gap_align->edit_script), &new_hsp);
                    // ```
                    // The init_hsp q_off/s_off values become the traceback
                    // gapped_start seed that BlastGetStartForGappedAlignmentNucl
                    // may keep or move.
                    blastn_trace::log(
                        "traceback",
                        format!(
                            "subject={}({}) context={} prelim=q{}..{} s{}..{} seed=({}, {}) raw_score={} tree_contains={} containing_hsp={:?}",
                            s_id,
                            s_idx,
                            prelim.context_idx,
                            prelim.prelim_qs,
                            prelim.prelim_qe,
                            prelim.prelim_ss,
                            prelim.prelim_se,
                            prelim.seed_qs,
                            prelim.seed_ss,
                            prelim.prelim_score,
                            prelim_traceback_contained,
                            traceback_containing_hsp
                        ),
                    );
                }
                if prelim_traceback_contained {
                    if let Some(timing) = timing_ref {
                        BlastnTiming::record_count(&timing.traceback_tree_precheck_skipped_hsps, 1);
                    }
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:404-405,607-610
                    // ```c
                    //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                    //                                hit_options->min_diag_separation)) {
                    // ...
                    //       } else {
                    //          /* Contained within another HSP, delete. */
                    // ```
                    // When the tree contains the HSP, NCBI deletes it and makes no traceback call. The helpers are
                    // told that nobody needs its value.
                    if let Some(x_ahead) = x_ahead.as_ref() {
                        x_ahead.skip(prelim_index);
                    }
                    continue;
                }

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4017
                // ```c
                // if (init_hsp->ungapped_data) {
                //     init_hsp->offsets.qs_offsets.q_off =
                //         init_hsp->ungapped_data->q_start + init_hsp->ungapped_data->length/2;
                //     init_hsp->offsets.qs_offsets.s_off =
                //         init_hsp->ungapped_data->s_start + init_hsp->ungapped_data->length/2;
                // }
                // ```
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:451-460
                // ```c
                // BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                // q_start = hsp->query.gapped_start;
                // s_start = hsp->subject.gapped_start;
                // ```
                // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:403-405,436-472
                // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, ...)) {
                //     BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                //     AdjustSubjectRange(&s_start, &adjusted_s_length,
                //                        q_start, query_length, &start_shift);
                // }
                // Consume only after the live ordered containment check. An HSP
                // skipped at batch start can become eligible after endpoint
                // replacement, so missing results still prepare and run here.
                #[cfg(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                ))]
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:403-405,509-512
                // ```c
                //       if (program_number == eBlastTypeRpsBlast ||
                //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                //                                hit_options->min_diag_separation)) {
                // ...
                //           BLAST_GappedAlignmentWithTraceback(program_number, query,
                //                 adjusted_subject, gap_align, score_params, q_start, s_start,
                //                 query_length, adjusted_s_length,
                //                 fence_hit);
                // ```
                // Look-ahead. This is reached after the tree test said "not contained". The
                // value of the traceback call comes from a helper, or
                // `traceback_hsp::trace_align_one` makes the call as before.
                let precomputed = if let Some(x_ahead) = x_ahead.as_ref() {
                    // No NCBI counterpart: a helper's value, if one has
                    // started on this HSP; otherwise `traceback_hsp::trace_align_one`
                    // evaluates it.
                    let mut x_scratch = XAheadScratch {
                        gap: &mut *gap_scratch,
                        greedy: &mut *greedy_align_scratch,
                    };
                    let x_taken = x_ahead.take(prelim_index, &mut x_scratch, &x_trace_compute);
                    x_taken.flatten()
                } else if speculative_traceback {
                    speculative_results[prelim_index % SPECULATIVE_TRACEBACK_BATCH_SIZE].take()
                } else {
                    None
                };
                #[cfg(not(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                )))]
                let precomputed: Option<((usize, usize, usize, usize, usize), _)> = None;
                let start_offsets_start = timing_enabled.then(std::time::Instant::now);
                let prepared = precomputed
                    .as_ref()
                    .map(|(prepared, _)| *prepared)
                    .or_else(|| traceback_hsp::prepare_traceback(traceback_ctx, prelim));
                if let (Some(timing), Some(start)) = (timing_ref, start_offsets_start) {
                    BlastnTiming::record_duration(&timing.traceback_start_offsets_ns, start);
                }
                let Some((
                    trace_q_start,
                    trace_s_start,
                    start_shift,
                    adjusted_s_len,
                    trace_s_start_adj,
                )) = prepared
                else {
                    continue;
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-520
                // (`BLAST_GreedyGappedAlignment` or
                // `BLAST_GappedAlignmentWithTraceback`, and the break on a fence). The
                // alignment of this HSP (`traceback_hsp::trace_align_one`): the
                // traceback call, or its value made ahead in `precomputed`. An
                // alignment that met a fence, or a greedy call that returns none, makes
                // no HSP.
                let trace_hsp = traceback_hsp::TraceHsp {
                    prelim,
                    ctx,
                    q_seq_blastna,
                    q_seq_nomask_blastna,
                    trace_traceback_seed,
                    subject_frame_sign,
                    trace_q_start,
                    trace_s_start,
                    start_shift,
                    adjusted_s_len,
                    trace_s_start_adj,
                };
                let Some(aligned) = traceback_hsp::trace_align_one(
                    traceback_ctx,
                    &trace_hsp,
                    precomputed,
                    gap_scratch,
                    greedy_align_scratch,
                ) else {
                    continue;
                };

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-605
                // (`Blast_HSPUpdateWithTraceback`, `Blast_HSPTest`,
                // `Blast_HSPAdjustSubjectOffset`, `BlastIntervalTreeAddHSP`). The HSP
                // this alignment makes (`traceback_hsp::trace_build_hsp`) goes into the
                // tree and to the end of `hits`, unless the identity and length test
                // deletes it.
                if !traceback_hsp::trace_build_hsp(
                    traceback_ctx,
                    &trace_hsp,
                    aligned,
                    traceback_hsp::TraceOut {
                        interval_tree,
                        hits,
                        traceback_edit_script_lengths,
                        traceback_alignment_lengths,
                    },
                ) {
                    continue;
                }
            }
        });
        // No NCBI counterpart: LOSAT_X_STATS work counters (xstats feature); they do not change any
        // value NCBI computes.
        if let Some(x_ahead) = x_ahead.as_ref() {
            let (x_evaluated, x_taken) = x_ahead.counts();
            crate::utils::xstats::AHEAD_TRACE_EVALUATED
                .fetch_add(x_evaluated, std::sync::atomic::Ordering::Relaxed);
            crate::utils::xstats::AHEAD_TRACE_TAKEN
                .fetch_add(x_taken, std::sync::atomic::Ordering::Relaxed);
        }
    }
}

/// What the traceback of one subject's kept HSPs reads besides the search's values and the
/// schedule: whether the query groups may run for this batch (LOSAT_X_QGROUPPAR,
/// `x_qgroup_batch`), the batch's queries (their number) and the lengths that bound a query group
/// task's interval tree (the concatenated queries and the whole subject). `Copy` values and a
/// shared reference.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:371-373
/// ```c
///    tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
///                                  0, (subject_length > 0 ? subject_length :
///                                  subject_blk->length / CODON_LENGTH) + 1);
/// ```
/// The bounds of each query group task's tree are those of the subject's tree. No NCBI
/// counterpart for the query groups themselves (scheduling only).
pub(super) struct TraceDrive<'q> {
    pub(super) x_qgroup_batch: bool,
    pub(super) queries: &'q [FastaRecord],
    pub(super) query_concat_length: usize,
    pub(super) s_len_full: usize,
}

/// The traceback of one subject's kept HSPs, from the traceback loop to the final lists: per query
/// on the pool (LOSAT_X_QGROUPPAR: each query's HSPs with a tree and scratch of the task's own,
/// `SubjectHits::ByQuery`), or the loop over the subject's list on this thread with the
/// subject's tree and scratch (`state`, `SubjectHits::Mixed`); each list then goes through the
/// purges and the final list (`traceback_finish`). Returns the hits and the lengths after each
/// step (for the debug line).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:259-261,375,611,637-638,717-719
/// ```c
/// Blast_TracebackFromHSPList(EBlastProgramType program_number,
///    BlastHSPList* hsp_list, const BLAST_SequenceBlk* query_blk,
///    BLAST_SequenceBlk* subject_blk, const BlastQueryInfo* query_info_in,
/// ...
///    for (index=0; index < num_initial_hsps; index++) {
/// ...
///    } /* End loop on HSPs */
/// ...
///        Int4 extra_start =
///            Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
/// ...
///        s_HSPListPostTracebackUpdate(program_number, hsp_list, query_info_in,
///                                     score_params, hit_params, sbp,
///                                     stat_length);
/// ```
/// `Blast_TracebackFromHSPList` as one subject's traceback: the loop (`trace_hsp_list`), then the
/// steps after it (`traceback_finish::traceback_post_purge`, `traceback_finish::traceback_final_list`).
/// No NCBI counterpart: the query groups on the pool (LOSAT_X_QGROUPPAR), scheduling only; NCBI
/// traces the hit lists of one subject in turn on one thread (blast_traceback.c:1632-1633).
pub(super) fn traceback_subject_hits(
    traceback_ctx: &TracebackCtx<'_>,
    schedule: &TraceSchedule<'_>,
    drive: &TraceDrive<'_>,
    prelim_hits: &[PrelimHit],
    state: TraceListState<'_>,
) -> (SubjectHits, [usize; 5]) {
    let &TraceSchedule { parallel_pool, .. } = schedule;
    let &TraceDrive {
        x_qgroup_batch,
        queries,
        query_concat_length,
        s_len_full,
    } = drive;
    let TraceListState {
        interval_tree,
        gap_scratch,
        greedy_align_scratch,
        x_helper_scratch,
        x_spec_scratch,
        hits,
        traceback_edit_script_lengths,
        traceback_alignment_lengths,
    } = state;
    // Dispatch point of LOSAT_X_QGROUPPAR for the traceback: the kept HSPs are the queries'
    // hit lists one after another, so each query's are one slice; per query on the pool, the
    // final lists concatenated in query order and sorted by score (stable, so ties stay in
    // query order, as the serial sort leaves them: its input, the second purge's output, is
    // ordered by context first; see x_qgroup_par.rs), or the loop over the subject's list.
    let x_trace_counts = (x_qgroup_batch
        && x_qgroup_par::in_query_order(prelim_hits.iter().map(|hit| hit.query_idx)))
    .then(|| x_qgroup_par::counts(prelim_hits.iter().map(|hit| hit.query_idx), queries.len()));
    let x_trace_order = x_trace_counts
        .as_deref()
        .and_then(|counts| x_qgroup_par::plan(counts, x_qgroup_par::Stage::Traceback));
    if let Some(order) = x_trace_order.as_ref() {
        let x_starts = x_qgroup_par::starts(x_trace_counts.as_deref().unwrap_or(&[]));
        let x_kept: &[PrelimHit] = prelim_hits;
        let x_tree_q_max = (query_concat_length + 1) as i32;
        let x_tree_s_max = (s_len_full + 1) as i32;
        let x_results = x_qgroup_par::run(
            parallel_pool,
            x_qgroup_par::Stage::Traceback,
            order,
            queries.len(),
            || {
                (
                    BlastIntervalTree::new(0, x_tree_q_max, 0, x_tree_s_max),
                    GapAlignScratch::new(),
                    GreedyAlignScratch::new(),
                )
            },
            |(tree, gap, greedy), query| {
                let query = query as usize;
                let mut edit_lengths = Vec::new();
                let mut alignment_lengths = Vec::new();
                // (the list holds the query's preliminary count, as NCBI's traced list does:
                // the subject path below)
                let mut hits = Vec::with_capacity(x_starts[query + 1] - x_starts[query]);
                trace_hsp_list(
                    traceback_ctx,
                    schedule,
                    &x_kept[x_starts[query]..x_starts[query + 1]],
                    TraceListState {
                        interval_tree: tree,
                        gap_scratch: gap,
                        greedy_align_scratch: greedy,
                        x_helper_scratch: None,
                        x_spec_scratch: None,
                        hits: &mut hits,
                        traceback_edit_script_lengths: &mut edit_lengths,
                        traceback_alignment_lengths: &mut alignment_lengths,
                    },
                );
                let (list, s1, s2, s3, s4) =
                    traceback_finish::traceback_post_purge(traceback_ctx, hits);
                let (list, s5) = traceback_finish::traceback_final_list(traceback_ctx, list, tree);
                (list, [s1, s2, s3, s4, s5], edit_lengths, alignment_lengths)
            },
        );
        // Each query's list is sorted by score (`traceback_final_list`), so the lists go to
        // the hit lists as they are, in query order (`SubjectHits::ByQuery`), not joined
        // into the subject's list and split by query again.
        let mut lists: Vec<Vec<BlastnHsp>> = Vec::new();
        let mut steps = [0usize; 5];
        for (list, step_counts, edit_lengths, alignment_lengths) in x_results.into_iter().flatten()
        {
            if !list.is_empty() {
                lists.push(list);
            }
            for (sum, step) in steps.iter_mut().zip(step_counts) {
                *sum += step;
            }
            traceback_edit_script_lengths.extend(edit_lengths);
            traceback_alignment_lengths.extend(alignment_lengths);
        }
        (SubjectHits::ByQuery(lists), steps)
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:296,312,375-376
        // ```c
        //    Int4 num_initial_hsps = hsp_list->hspcnt;
        // ...
        //    hsp_array = hsp_list->hsp_array;
        // ...
        //    for (index=0; index < num_initial_hsps; index++) {
        //       hsp = hsp_array[index];
        // ```
        // NCBI traces the preliminary list in its own array, so the traced list never
        // holds more than the preliminary count. No NCBI counterpart (memory placement
        // only): the list here gets that count once, instead of growing by doubling (up to
        // twice the traced count) and being copied or shrunk afterwards.
        hits.reserve_exact(prelim_hits.len());
        trace_hsp_list(
            traceback_ctx,
            schedule,
            prelim_hits,
            TraceListState {
                interval_tree,
                gap_scratch,
                greedy_align_scratch,
                x_helper_scratch,
                x_spec_scratch,
                hits,
                traceback_edit_script_lengths,
                traceback_alignment_lengths,
            },
        );
        // No NCBI counterpart (memory placement only): the list goes to the purges as it is
        // (they work in its allocation, `x_in_handle_order`) and leaves as the subject's
        // result, instead of a copy while the scratch keeps the buffer at its largest size;
        // the next subject's list gets that subject's preliminary count (above).
        let (list, s1, s2, s3, s4) =
            traceback_finish::traceback_post_purge(traceback_ctx, std::mem::take(hits));
        let (list, s5) = traceback_finish::traceback_final_list(traceback_ctx, list, interval_tree);
        (SubjectHits::Mixed(list), [s1, s2, s3, s4, s5])
    }
}
