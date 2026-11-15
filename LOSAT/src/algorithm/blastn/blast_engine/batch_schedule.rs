//! The schedule of one blastn query batch's search (`batch_schedule`): the `LOSAT_TIMING` clocks, whether
//! the search asks for more than one thread (`requested_parallel`), the task's configuration
//! (`configure_task`) and the traceback's speculative scheduling (the speculative greedy traceback, the
//! look-ahead window of the helper threads). It reads the switches and the environment once per batch.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:82-88 (the thread count of a
//! preliminary search)
//! ```c
//!     if (num_threads > 1) {
//!         SetNumberOfThreads(num_threads);
//!     }
//! ```
//!
//! This module is a child of `run` for the types it keeps private (`BlastnTiming`).

use super::{configure_task, BlastnArgs, BlastnTiming};
use crate::algorithm::blastn::coordination::TaskConfig;
use std::sync::Arc;

/// What `search_query_batch` reads of the batch's schedule: the clocks, the thread request, the task's
/// configuration (finished later by `finalize_task_config`) and the speculative traceback's rules.
pub(super) struct BatchSchedule {
    pub(super) timing_enabled: bool,
    pub(super) t_total: Option<std::time::Instant>,
    pub(super) timing: Option<Arc<BlastnTiming>>,
    pub(super) requested_parallel: bool,
    pub(super) config: TaskConfig,
    pub(super) x_greedy_spec: bool,
    pub(super) speculative_traceback: bool,
    pub(super) x_ahead_window: Option<usize>,
}

/// The schedule of a batch searched with `num_threads` threads.
pub(super) fn batch_schedule(
    args: &BlastnArgs,
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        expect(
            unused_variables,
            reason = "read only by the parallel request, built with the parallel feature outside serial Wasm"
        )
    )]
    num_threads: usize,
) -> BatchSchedule {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:493-555
    // ```c
    // if (aux_struct->WordFinder) {
    //     aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
    //                            word_params, aux_struct->ewp,
    //                            aux_struct->offset_pairs,
    //                            kScanSubjectOffsetArraySize,
    //                            init_hitlist, ungapped_stats);
    //     if (init_hitlist->total == 0) continue;
    // }
    // ...
    // if (aux_struct->GetGappedScore) {
    //     status = aux_struct->GetGappedScore(program_number, query,
    //             query_info, subject, gap_align, score_params, ext_params,
    //             hit_params, word_params, init_hitlist, &hsp_list,
    //             gapped_stats, NULL);
    // }
    // ...
    // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    // Blast_HSPListSortByScore(hsp_list);
    // ```
    let timing_enabled = std::env::var_os("LOSAT_TIMING").is_some();
    let t_total = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    let timing = if timing_enabled {
        Some(Arc::new(BlastnTiming::new()))
    } else {
        None
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:478-536
    // ```c
    // while (TRUE) {
    //     status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
    //                                    dbseq_chunk_overlap);
    //     if (status == SUBJECT_SPLIT_DONE) break;
    //     if (status == SUBJECT_SPLIT_NO_RANGE) continue;
    //     ...
    //     if (aux_struct->WordFinder) {
    //         aux_struct->WordFinder(...);
    //         if (init_hitlist->total == 0) continue;
    //     }
    //     ...
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let requested_parallel = num_threads > 1;
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let requested_parallel = false;

    // Configure task-specific parameters (initial configuration)
    let config = configure_task(args);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:509-513
    // ```c
    // BLAST_GappedAlignmentWithTraceback(program_number, query,
    //     adjusted_subject, gap_align, score_params, q_start, s_start,
    //     query_length, adjusted_s_length, fence_hit);
    // ```
    // Preserve ordered debug-coordinate logging; speculative DP is otherwise pure.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,503-507
    // ```c
    //    for (index=0; index < num_initial_hsps; index++) {
    // ...
    //         } else if (kGreedyTraceback) {
    //             BLAST_GreedyGappedAlignment(query, adjusted_subject,
    //                 query_length, adjusted_s_length, gap_align,
    //                 score_params, q_start, s_start, FALSE, TRUE,
    //                 fence_hit);
    // ```
    // The dynamic-programming traceback and, in a parallel search, the greedy traceback of megablast
    // (the call above) are scheduled ahead of the ordered loop (the rule of `LOSAT_X_GREEDYSPEC`,
    // kept as its `x_reach!` key). The ordered loop still tests the tree and adds to it in the NCBI
    // order. This is scheduling only.
    // No NCBI counterpart: the same scheduling for the greedy
    // traceback of megablast. A greedy alignment is a function of its
    // arguments only: every cell it reads from its scratch was written
    // earlier in the same call, and a larger `max_dist` left by an earlier
    // call only saves retries.
    let x_greedy_spec = requested_parallel && !config.use_dp;
    let speculative_traceback = requested_parallel
        && (config.use_dp || x_greedy_spec)
        && std::env::var_os("LOSAT_DEBUG_COORDS").is_none()
        && std::env::var_os("LOSAT_DEBUG_COORDS_START").is_none();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3918-3919
    // ```c
    //    for (index=0; index<init_hitlist->total; index++)
    // ...
    //       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
    //                                         hit_options->min_diag_separation))
    // ```
    // Look-ahead for the preliminary loop (`prelim_hsp_list::gapped_score_list`) and the traceback
    // loop (`Blast_TracebackFromHSPList`, next block). The helper threads only evaluate the
    // alignment call of an index ahead of the loop. The loops, their tree tests and their tree
    // insertions stay on the calling thread in the NCBI order.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,404-405
    // ```c
    //    for (index=0; index < num_initial_hsps; index++) {
    // ...
    //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
    //                                hit_options->min_diag_separation)) {
    // ```
    // The same for the loop of `Blast_TracebackFromHSPList`.
    // No NCBI counterpart: the other threads of the pool
    // evaluate the gapped alignments of the two ordered loops (preliminary
    // and traceback) ahead of the loop; see utils/xahead.rs. Same conditions
    // as the speculative traceback above.
    let x_ahead_window = if speculative_traceback {
        Some(crate::utils::xahead::WINDOW)
    } else {
        None
    };
    if x_ahead_window.is_none() {
        crate::utils::xreach::x_reach!("LOSAT_X_AHEAD", false);
    }
    BatchSchedule {
        timing_enabled,
        t_total,
        timing,
        requested_parallel,
        config,
        x_greedy_spec,
        speculative_traceback,
        x_ahead_window,
    }
}
