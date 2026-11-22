//! The composition-based redo of a BLASTP call's matches (`run_call`, after the choice of the
//! ordinary traceback): the parameters of each query's redo (`RedoParamsBuilder`, which
//! `run_call`'s `build_redo_align_params` calls), each query's batch (`kappa_query_batches`), the
//! plan of the redo (`plan_redo_and_expand_lists`), which chooses its path (`RedoPath`), and the
//! four paths, which read `RedoInputs`: the windows of `LOSAT_X_BPREDOPAR` (`x_windows_redo`), the
//! matches of one query on the pool (`match_parallel_redo`), the queries on the pool
//! (`query_parallel_redo`) and the one-thread loop (`one_thread_redo`). Each function is the
//! statements of `run_call` at its place; `run_call` calls the plan, then the path that it chose,
//! and its `?` passes the path's errors on.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2977-2981
//!
//! This module is a child of `blast_engine` for the private items of that file it reads: the six
//! `blastp_*` functions that `build` calls (`blastp_local_scaling_factor`,
//! `blastp_compo_adjust_mode`, `blastp_redo_matrix`, `blastp_redo_near_identical_cutoff`,
//! `blastp_redo_x_dropoff`, `blastp_redo_cutoff_score`), `blastp_use_matrix_workspace`,
//! `preliminary_hit_from_local_match_hsp`, `query_nomask_sequence`, `X_NEWTON_LANES_MATCHES`, and
//! for the parallel paths' choice `blastp_should_parallelize_stage` and `x_prelim_hsplist_count`.

use crate::utils::glibc_fma_libm::GlibcLibm;
use std::ops::Range;
use std::sync::OnceLock;

use anyhow::Result;

use super::{
    blastp_compo_adjust_mode, blastp_local_scaling_factor, blastp_redo_cutoff_score,
    blastp_redo_matrix, blastp_redo_near_identical_cutoff, blastp_redo_x_dropoff,
    blastp_use_matrix_workspace, preliminary_hit_from_local_match_hsp, query_nomask_sequence,
    X_NEWTON_LANES_MATCHES,
};
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::blastp::gapalign::GapAlignScratch;
use crate::algorithm::blastp::hsp::{blast_compo_early_termination, BlastCompoHeap, BlastpHitList};
use crate::algorithm::blastp::kappa::{
    self, build_query_workspace, postprocess_preliminary_hits, x_newton_lanes_prefetch_matches,
    BlastRedoAlignParams, BlastpKappaSubjectRangeCache,
};
use crate::algorithm::blastp::x_prelim_compact::{
    x_expand_hit_list, x_finish_compact_lists, XPrelimHsp,
};
use crate::algorithm::blastp::x_query_blocks;
use crate::algorithm::blastp::x_redo_subjects::XRedoSubjects;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::blastinput::seq_range::RangedBatch;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::{
    build_matrix_info, BlastCompositionWorkspace, BlastMatrixInfo,
};
use crate::core::composition_adjustment::redo_alignment::BlastCompoGappingParams;
use crate::report::record::ReportNeeds;
use crate::stats::search_space::SearchSpace;
use crate::stats::{BlastGumbelBlk, KarlinParams};
use crate::utils::matrix::BLASTAA_SIZE;
use crate::utils::threading::SearchPool;

// The imports of the parallel code (the plan's choice of a parallel path and the parallel paths).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use {
    super::{blastp_should_parallelize_stage, x_prelim_hsplist_count},
    crate::algorithm::blastp::kappa::BlastpPostprocessResult,
    crate::algorithm::blastp::x_redo_par::{self, XRedoParSearch},
    rayon::prelude::*,
};

/// What the parameters of a query's redo are made of (`build`): each field is `run_call`'s value
/// of the same name (`Copy` scalars by value, the rest by reference, so that the moved body reads
/// them as it did in `run_call`), and the start matrix that the
/// parameters of every query share, made at the first `build`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2401-2407
pub(super) struct RedoParamsBuilder<'a> {
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) contexts: &'a [QueryContext],
    pub(super) gapped_params: &'a KarlinParams,
    pub(super) hit_cutoff_score_min: i32,
    pub(super) ideal_ungapped_lambda: f64,
    pub(super) x_drop_gapped_final: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2228-2236
    //                              self->startFreqRatios, self->ungappedLambda);
    // This non-position-based initial matrix is invariant within this search.
    // Initialize at the existing first params request, preserving no-hit/error
    // boundaries, and clone completed values into independently mutable params.
    pub(super) initial_matrix_info: OnceLock<Result<BlastMatrixInfo>>,
}

impl RedoParamsBuilder<'_> {
    /// The parameters of the redo of query `q_idx` (`run_call`'s `build_redo_align_params`, which
    /// the redo's paths call for each query, or each match, that they redo); the start matrix is
    /// made at the first call, and each call clones it (or returns its error).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2471-2479
    pub(super) fn build(&self, q_idx: usize) -> Result<BlastRedoAlignParams> {
        let RedoParamsBuilder {
            args,
            contexts,
            gapped_params,
            hit_cutoff_score_min,
            ideal_ungapped_lambda,
            x_drop_gapped_final,
            ref initial_matrix_info,
        } = *self;
        let local_scaling_factor =
            blastp_local_scaling_factor(args.comp_based_stats.mode, args.scoring.matrix);
        let scaled_gapped_lambda = gapped_params.lambda / local_scaling_factor;
        let do_link_hsps = args.sum_stats;
        let matrix_info = initial_matrix_info
            .get_or_init(|| {
                build_matrix_info(
                    blastp_redo_matrix(args.scoring.matrix),
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2216-2232
                    ideal_ungapped_lambda / local_scaling_factor,
                )
            })
            .as_ref()
            .map_err(|error| anyhow::anyhow!("{error}"))?
            .clone();
        Ok(BlastRedoAlignParams {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2216-2232
            matrix_info,
            // NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:97-112
            //
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2117-2133
            gapping_params: BlastCompoGappingParams {
                gap_open: ((args.scoring.gap_open as f64) * local_scaling_factor).round() as i32,
                gap_extend: ((args.scoring.gap_extend as f64) * local_scaling_factor).round()
                    as i32,
                x_dropoff: blastp_redo_x_dropoff(
                    args.x_drop_final_bits,
                    scaled_gapped_lambda,
                    x_drop_gapped_final,
                ),
                context: std::cell::Cell::new(None),
            },
            compo_adjust_mode: blastp_compo_adjust_mode(args.comp_based_stats.mode),
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1107-1118
            alphsize: BLASTAA_SIZE as i32,
            composition_test_index: i32::from(args.comp_based_stats.unified_p),
            log_k: gapped_params.k.glibc_ln(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3584-3690
            score_divisor: local_scaling_factor,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039
            // `-use_sw_tback` (item X4): the one-thread loop redoes each match with NCBI's
            // Smith-Waterman redo (`one_thread_redo`).
            smith_waterman: args.use_sw_tback,
            near_identical_cutoff: blastp_redo_near_identical_cutoff(scaled_gapped_lambda),
            position_based: false,
            ccat_query_length: contexts[q_idx].aa_len as i32,
            query_is_translated: false,
            subject_is_translated: false,
            cutoff_score: blastp_redo_cutoff_score(
                do_link_hsps,
                hit_cutoff_score_min,
                local_scaling_factor,
            ),
            do_link_hsps,
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:382-406
        })
    }
}

/// Each query's context of the E-value adjustment with a composition P-value, in query order
/// (`kappa::evalue_adjust_contexts`): the next query's of its batch (`kappa_query_batches`), or
/// its own. With unified P-values a batch whose adjustment NCBI makes with a search space of 0
/// is refused here, by the ids the reports show (`query_ids`).
pub(super) fn evalue_adjust_contexts(
    args: &ResolvedBlastpArgs,
    contexts: &[QueryContext],
    search_spaces: &[SearchSpace],
    kappa_query_batches: &[Range<usize>],
    query_ids: &[std::sync::Arc<[u8]>],
) -> Result<Vec<kappa::EvalueAdjustContext>> {
    let refuse = args
        .comp_based_stats
        .unified_p
        .then(|| kappa::UnifiedPRefusal {
            program: args.app.name(),
            query_ids,
        });
    kappa::evalue_adjust_contexts(
        contexts,
        search_spaces,
        kappa_query_batches,
        args.searchsp.unwrap_or(0),
        refuse,
    )
}

/// Each query's batch as a range of the queries, one range per query in query order
/// (`query_count`): the heaps that the early termination test of the query's redo reads.
pub(super) fn kappa_query_batches(
    query_count: usize,
    query_batches: &[RangedBatch],
) -> Vec<Range<usize>> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-260,274-276
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3013,3119-3128
    // NCBI searches one query batch at a time, so the heaps that the redo of a query reads
    // (BlastCompo_EarlyTermination) are those of the queries of its batch. `redone_matches`
    // holds the heaps of all batches; this is each query's batch as a range of it. A query
    // outside every batch is not searched and stands alone.
    let mut kappa_query_batches: Vec<Range<usize>> =
        (0..query_count).map(|q_idx| q_idx..q_idx + 1).collect();
    for batch in query_batches {
        for q_idx in batch.searched.clone() {
            if let Some(range) = kappa_query_batches.get_mut(q_idx) {
                *range = batch.searched.clone();
            }
        }
    }
    kappa_query_batches
}

/// What the redo's plan and its four paths read of `run_call`: each field is `run_call`'s value of
/// the same name (the `Copy` scalars by value, the rest by reference, so that the moved bodies read
/// them as they did in `run_call`; `subjects` are the redo's subjects, `kappa_query_batches` holds
/// one range per query, and `build_redo_align_params` makes the parameters of a query's redo with
/// `RedoParamsBuilder`). The last three fields only the parallel code reads (the plan's choice of
/// a parallel path and the windows' `XRedoParSearch`). The paths' threads read these values, as
/// NCBI's threads share those of its parallel region.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3432-3445
pub(super) struct RedoInputs<'a, P> {
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) blastp_parallel_pool: &'a SearchPool<'a>,
    pub(super) needs: ReportNeeds,
    pub(super) subjects: XRedoSubjects<'a>,
    pub(super) contexts: &'a [QueryContext],
    pub(super) search_spaces: &'a [SearchSpace],
    pub(super) gapped_params: &'a KarlinParams,
    pub(super) gapped_gumbel: &'a Option<BlastGumbelBlk>,
    pub(super) kappa_query_batches: &'a [Range<usize>],
    /// Each query's context of the E-value adjustment with a composition P-value
    /// (`kappa::evalue_adjust_contexts`).
    pub(super) evalue_adjust_contexts: &'a [kappa::EvalueAdjustContext],
    pub(super) build_redo_align_params: &'a P,
    /// The marks of the lists that the redo reaches and its fetch excludes, by query and place
    /// (`database::exclude_fetched`; empty without a database or a negative seqid list): only the
    /// status of the Smith-Waterman redo reads them (`one_thread_redo`); the X-drop redo has no
    /// status, so its paths read the lists without them.
    pub(super) fetch_excluded: &'a [super::database::FetchExcludedMatch],
    pub(super) x_win_expand: bool,
    pub(super) ordinary_traceback: bool,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) blastp_num_threads: usize,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) timing_enabled: bool,
}

/// The path of the redo that its plan chose, with what only that path reads: the windows of
/// `LOSAT_X_BPREDOPAR` (`x_windows_redo`), the matches of one query on the pool
/// (`match_parallel_redo`, that query), the queries on the pool (`query_parallel_redo`, the queries
/// with matches) or the one-thread loop (`one_thread_redo`), the only path of a build without
/// threads. NCBI has one loop over the matches, which its threads split
/// (`plan_redo_and_expand_lists`).
pub(super) enum RedoPath {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    Windows,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    MatchParallel(usize),
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    QueryParallel(Vec<usize>),
    OneThread,
}

/// The plan of the redo: the path that runs (`RedoPath`), chosen in the order in which the
/// paths are tried (the windows, the matches of one query, the queries, else the one-thread
/// loop), and the reach records of the sides that the chosen path takes, in their order
/// (`LOSAT_X_BPREDOPAR`'s read switch, the graduated SEGSHARE and NEWTONLANES,
/// `LOSAT_X_BPSEGWIN`). Between the first record and the others it expands the compact lists that
/// the windows would have expanded, where they could run (`x_win_expand`) but do not.
///
/// NCBI chooses only the number of threads of its one loop (`actual_num_threads`); it has no
/// counterpart of the choice among LOSAT's paths, each of which gives the output of the one loop.
///
/// With `-use_sw_tback` (item X4) the one-thread loop runs: the Smith-Waterman redo of a match
/// asks the query's heap before its first alignment (`BlastCompo_HeapWouldInsert`), so a match's
/// redo reads the heap that the matches before it in the loop's order leave.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3079-3082
pub(super) fn plan_redo_and_expand_lists<P>(
    redo: &RedoInputs<'_, P>,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) -> RedoPath {
    let RedoInputs {
        x_win_expand,
        ordinary_traceback,
        ..
    } = *redo;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let RedoInputs {
        args,
        subjects,
        blastp_num_threads,
        ..
    } = *redo;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let kappa_parallel_query_indices: Vec<usize> = (0..preliminary_hit_lists.len())
        .filter(|&q_idx| {
            x_prelim_hsplist_count(preliminary_hit_lists, x_compact_hit_lists, q_idx) > 0
        })
        .collect();
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let use_parallel_kappa_redo = {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3429-3449
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:405-409,417-418
        // The query-parallel redo below relies on a heap that no match has reached not being
        // filled to its cutoff (hitlist_size >= 1); otherwise the one-thread loop is used.
        kappa_parallel_query_indices.len() > 1
            && !args.use_sw_tback
            && !args.redone_match_heap().filled_to_cutoff()
            && blastp_should_parallelize_stage(
                "kappa-query-redo",
                args.num_threads,
                blastp_num_threads,
                subjects.len(),
                kappa_parallel_query_indices.len(),
            )
    };
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let use_parallel_kappa_redo = false;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let kappa_match_parallel_query_index = {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3429-3449
        if use_parallel_kappa_redo || kappa_parallel_query_indices.len() != 1 || args.use_sw_tback {
            None
        } else {
            let q_idx = kappa_parallel_query_indices[0];
            let match_jobs =
                x_prelim_hsplist_count(preliminary_hit_lists, x_compact_hit_lists, q_idx);
            if blastp_should_parallelize_stage(
                "kappa-match-redo",
                args.num_threads,
                blastp_num_threads,
                subjects.len(),
                match_jobs,
            ) {
                Some(q_idx)
            } else {
                None
            }
        }
    };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let use_parallel_kappa_match_redo = kappa_match_parallel_query_index.is_some();
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let use_parallel_kappa_match_redo = false;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3398-3399,3447-3449
    // EXPERIMENT (LOSAT_X_BPREDOPAR, a strict-set switch): NCBI splits one loop over the matches
    // of all queries between its threads. With the switch, more than one thread and at least two
    // matches, the matches run on the pool in jobs of up to four matches of one query, and the
    // early termination test, the heap insertions and the error that ends the redo are replayed
    // in the order of the one-thread loop below; it replaces both parallel paths
    // (`x_redo_par.rs`).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let x_redo_par_switch = x_redo_par::x_bp_redo_par();
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let x_redo_par_on = x_redo_par_switch
        && !args.use_sw_tback
        && blastp_should_parallelize_stage(
            "kappa-redo-par",
            args.num_threads,
            blastp_num_threads,
            subjects.len(),
            // The compact lists count while not expanded.
            x_redo_par::x_matches(
                preliminary_hit_lists,
                if x_win_expand {
                    x_compact_hit_lists
                } else {
                    &[]
                },
            ),
        );
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_redo_par_switch && !ordinary_traceback {
        crate::utils::xreach::x_reach!("LOSAT_X_BPREDOPAR", x_redo_par_on);
    }
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let x_redo_par_run = x_redo_par_on;
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let x_redo_par_run = false;
    // The other paths read the expanded lists of all queries.
    if x_win_expand && !x_redo_par_run {
        x_finish_compact_lists(preliminary_hit_lists, x_compact_hit_lists);
    }

    // Graduated 2026-10-22 (was `LOSAT_X_SEGSHARE`): the query-parallel redo and the windows of
    // `x_redo_par` share the masked ranges; the match-parallel redo shares none, and the
    // one-thread loop shares none.
    // The one-thread (`one_thread_redo`) and the query-parallel redo loops solve the first Newton
    // problems of the next `X_NEWTON_LANES_MATCHES` matches ahead (Newton lanes); the
    // match-parallel redo solves none.
    // With `-comp_based_stats 0` no redo runs, so neither side is reached.
    if !ordinary_traceback {
        crate::utils::xreach::x_reach!(
            "LOSAT_X_SEGSHARE",
            x_redo_par_run || (!use_parallel_kappa_match_redo && use_parallel_kappa_redo)
        );
        // `-use_sw_tback` solves no Newton problem ahead (`x_newton_lanes_prefetch_matches`
        // returns at once for the Smith-Waterman redo).
        crate::utils::xreach::x_reach!(
            "LOSAT_X_NEWTONLANES",
            (x_redo_par_run || !use_parallel_kappa_match_redo) && !redo.args.use_sw_tback
        );
    }
    // EXPERIMENT (LOSAT_X_BPSEGWIN): a redo not in windows empties no store; the windows record
    // their side (`x_redo_par::x_redo`).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if kappa::x_bp_seg_win().is_some() && !x_redo_par_on && !ordinary_traceback {
        crate::utils::xreach::x_reach!("LOSAT_X_BPSEGWIN", false);
    }
    // The path, in the order of the choice: the windows, the matches of one query, the queries.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_redo_par_run {
        return RedoPath::Windows;
    } else if let Some(q_idx) = kappa_match_parallel_query_index {
        return RedoPath::MatchParallel(q_idx);
    } else if use_parallel_kappa_redo {
        return RedoPath::QueryParallel(kappa_parallel_query_indices);
    }
    RedoPath::OneThread
}

/// The redo in windows of `LOSAT_X_BPREDOPAR` (`x_redo_par::x_redo_par_run`, which reads
/// `XRedoParSearch`, made here of `RedoInputs`): the matches of all queries on the pool in jobs,
/// replayed in stream order, the compact lists expanded as the windows reach them (`x_win_expand`);
/// the lists it has redone are released, and counted in `x_prelim_freed`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3414-3418
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn x_windows_redo<P>(
    redo: &RedoInputs<'_, P>,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    redone_matches: &mut [BlastCompoHeap],
    x_prelim_freed: &mut usize,
) -> Result<()>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let RedoInputs {
        args,
        blastp_parallel_pool,
        needs,
        subjects,
        contexts,
        search_spaces,
        gapped_params,
        gapped_gumbel,
        kappa_query_batches,
        evalue_adjust_contexts,
        build_redo_align_params,
        x_win_expand,
        blastp_num_threads,
        timing_enabled,
        ..
    } = *redo;
    // What `x_redo_par` reads of the redo.
    let x_redo_par_search = XRedoParSearch {
        pool: blastp_parallel_pool,
        threads: blastp_num_threads,
        contexts,
        query_nomask: query_nomask_sequence,
        search_spaces,
        subjects,
        query_batches: kappa_query_batches,
        evalue_adjust_contexts,
        matrix: args.scoring.matrix,
        gapped_params,
        gapped_gumbel: gapped_gumbel.as_ref(),
        evalue: args.evalue,
        build_redo_align_params,
        timing: timing_enabled,
        engine_positives: needs.seqalign,
    };
    // The parallel paths read the expanded lists of all queries (`plan_redo_and_expand_lists`
    // expands them before the paths), or the windows expand the compact lists as they reach them.
    debug_assert!(x_win_expand || x_compact_hit_lists.iter().all(Option::is_none));
    let released = x_redo_par::x_redo_par_run(
        &x_redo_par_search,
        preliminary_hit_lists,
        x_win_expand.then_some(&mut x_compact_hit_lists[..]),
        redone_matches,
    )?;
    *x_prelim_freed += released;
    Ok(())
}

/// The redo of the matches of one query on the pool (`q_idx`, `RedoPath::MatchParallel`): the matches
/// are redone on the pool, each work unit with its own scratch, and the early termination test and
/// the heap insertions are replayed after them in stream order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3511-3511
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn match_parallel_redo<P>(
    redo: &RedoInputs<'_, P>,
    q_idx: usize,
    preliminary_hit_lists: &[Option<BlastpHitList>],
    redone_matches: &mut [BlastCompoHeap],
) -> Result<()>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let RedoInputs {
        args,
        blastp_parallel_pool,
        needs,
        subjects,
        contexts,
        search_spaces,
        gapped_params,
        gapped_gumbel,
        kappa_query_batches,
        evalue_adjust_contexts,
        build_redo_align_params,
        ..
    } = *redo;
    let preliminary_hit_list = preliminary_hit_lists[q_idx]
        .as_ref()
        .expect("single-query Kappa redo parallel path requires preliminary hits");
    let local_matches = &preliminary_hit_list.hsplist_array[..preliminary_hit_list.hsplist_count];
    let ctx = &contexts[q_idx];
    let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
    let query_nomask = query_nomask_sequence(ctx);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3429-3449
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3493-3503
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3429-3459
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2309-2336,3493-3504
    // Only immutable query preparation is shared within this search.
    let query_workspace = std::sync::OnceLock::new();
    let pool = blastp_parallel_pool;
    let precomputed: Vec<Result<Option<BlastpPostprocessResult>>> = pool.install(|| {
        local_matches
            .par_iter()
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334,3493-3503
            // Each Rayon work unit owns its mutable scratch. Initialization
            // count is scheduler-dependent, never assumed to equal threads.
            .map_init(
                || {
                    let mut composition_workspace = BlastCompositionWorkspace::new_blosum62();
                    if args.scoring.matrix != ScoringMatrix::Blosum62 {
                        blastp_use_matrix_workspace(
                            &mut composition_workspace,
                            args.scoring.matrix,
                        );
                    }
                    (composition_workspace, GapAlignScratch::new(), Vec::new())
                },
                |(composition_workspace, kappa_gap_scratch, kappa_preliminary_hits),
                 local_match|
                 -> Result<Option<BlastpPostprocessResult>> {
                    kappa_preliminary_hits.clear();
                    kappa_preliminary_hits.extend(
                        local_match
                            .hsps
                            .iter()
                            .map(preliminary_hit_from_local_match_hsp),
                    );
                    if kappa_preliminary_hits.is_empty() {
                        return Ok(None);
                    }

                    let s_idx = local_match.oid;
                    let redo_align_params = build_redo_align_params(q_idx)?;
                    // Read from the database here.
                    let subject = subjects.read(s_idx)?;
                    let subject_raw = subject.residues();
                    let query_workspace = query_workspace.get_or_init(|| {
                        build_query_workspace(
                            query_raw,
                            ctx.aa_len as i32,
                            search_spaces[q_idx].length_adjustment as i32,
                            search_spaces[q_idx].effective_space,
                            evalue_adjust_contexts[q_idx],
                        )
                    });
                    let mut kappa_subject_range_cache = BlastpKappaSubjectRangeCache::new(None);

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3577-3585
                    postprocess_preliminary_hits(
                        kappa_preliminary_hits,
                        query_workspace,
                        composition_workspace,
                        query_nomask,
                        subject_raw,
                        args.scoring.matrix,
                        gapped_params,
                        gapped_gumbel.as_ref(),
                        &redo_align_params,
                        args.evalue,
                        kappa_gap_scratch,
                        &mut kappa_subject_range_cache,
                        needs.seqalign,
                    )
                    .map(Some)
                },
            )
            .collect()
    });

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3817-3849
    //
    // The expensive redo calculation above runs in NCBI's static-loop
    // spirit, but heap admission is replayed in the original HSP stream
    // order so LOSAT keeps thread-count-independent output.
    let batch = kappa_query_batches[q_idx].clone();
    for (local_match, postprocessed) in local_matches.iter().zip(precomputed.into_iter()) {
        if local_match.hsps.is_empty() {
            continue;
        }
        if blast_compo_early_termination(local_match.best_evalue, &redone_matches[batch.clone()]) {
            continue;
        }
        let Some(postprocessed) = postprocessed? else {
            continue;
        };
        let best_score = postprocessed.best_score;
        let best_evalue = postprocessed.best_evalue;
        let Some(hsp_list) = postprocessed.hsp_list else {
            continue;
        };
        let redone_match = &mut redone_matches[q_idx];
        if redone_match.would_insert(best_evalue, best_score, local_match.oid) {
            redone_match.insert(hsp_list, best_evalue, best_score, local_match.oid);
        }
    }
    Ok(())
}

/// The redo of the queries on the pool (`kappa_parallel_query_indices`, `RedoPath::QueryParallel`),
/// each query with its own heap; the last query of a batch of several keeps its results, whose
/// early termination test and heap insertions are replayed after the queries in stream order,
/// with the heaps of its batch.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3282-3283,3288-3291
/// NCBI's threads split the matches, each thread with its own heaps; LOSAT runs the queries on the
/// pool, each query with its own heap.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn query_parallel_redo<P>(
    redo: &RedoInputs<'_, P>,
    kappa_parallel_query_indices: &[usize],
    preliminary_hit_lists: &[Option<BlastpHitList>],
    redone_matches: &mut [BlastCompoHeap],
) -> Result<()>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let RedoInputs {
        args,
        blastp_parallel_pool,
        needs,
        subjects,
        contexts,
        search_spaces,
        gapped_params,
        gapped_gumbel,
        kappa_query_batches,
        evalue_adjust_contexts,
        build_redo_align_params,
        ..
    } = *redo;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3119-3128
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3436-3459
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3429-3459
    let pool = blastp_parallel_pool;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636,3448-3449
    // Dispatch point: NCBI SEGs the subject range inside the per-match loop (omp for over
    // matches). Here the masked ranges are kept in one store for all queries, as the serial
    // redo already does. Reuse of an identical result only.
    // One store of masked subject ranges for all queries, as in `one_thread_redo`.
    let x_shared_ranges = Some(std::sync::Arc::new(
        crate::algorithm::blastp::kappa::XSharedSubjectRanges::new(),
    ));
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1580
    // The test reads the heaps of the queries of the batch (`kappa_query_batches`).
    // `one_thread_redo` redoes the queries in order, so when it tests a match of one
    // query, the heaps of the queries of its batch before it are complete and those of
    // the queries after it are as BlastCompo_HeapInitialize left them: not filled to the
    // cutoff (`use_parallel_kappa_redo`), which makes the test FALSE. The queries are
    // redone here at the same time, with the same outcome:
    // - a query that is not the last one of its batch redoes every match;
    // - a query that is alone in its batch tests its own heap, the only one of the batch;
    // - the last query of a batch of several redoes every match and keeps the results;
    //   the test and the heap insertions are replayed in stream order after the parallel
    //   loop, when the heaps of the other queries of the batch are complete.
    //
    // A heap that no match has reached, for the Newton lanes prefetch: the
    // prefetch leaves out the matches that the test ends for the heaps it is given, and
    // a query that does not test its own heap redoes every match.
    let x_untouched_heap = args.redone_match_heap();
    // (query, its heap, the kept results of a batch's last query by match position)
    type KappaQueryRedo = (
        usize,
        BlastCompoHeap,
        Vec<(usize, Result<BlastpPostprocessResult>)>,
    );
    let query_heaps: Result<Vec<KappaQueryRedo>> = pool.install(|| {
        kappa_parallel_query_indices
            .par_iter()
            .map(|&q_idx| -> Result<KappaQueryRedo> {
                let mut composition_workspace = BlastCompositionWorkspace::new_blosum62();
                if args.scoring.matrix != ScoringMatrix::Blosum62 {
                    blastp_use_matrix_workspace(&mut composition_workspace, args.scoring.matrix);
                }
                let mut kappa_preliminary_hits = Vec::new();
                let mut kappa_gap_scratch = GapAlignScratch::new();
                // See XSharedSubjectRanges.
                let mut kappa_subject_range_cache =
                    BlastpKappaSubjectRangeCache::new(x_shared_ranges.clone());
                let mut redone_match = args.redone_match_heap();
                let mut kept_for_replay = Vec::new();
                let batch = &kappa_query_batches[q_idx];
                let alone_in_batch = batch.len() == 1;
                let replayed = !alone_in_batch && q_idx + 1 == batch.end;
                let redo_align_params = build_redo_align_params(q_idx)?;

                let ctx = &contexts[q_idx];
                let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
                let query_nomask = query_nomask_sequence(ctx);
                let query_workspace = build_query_workspace(
                    query_raw,
                    ctx.aa_len as i32,
                    search_spaces[q_idx].length_adjustment as i32,
                    search_spaces[q_idx].effective_space,
                    evalue_adjust_contexts[q_idx],
                );
                let Some(preliminary_hit_list) = preliminary_hit_lists[q_idx].as_ref() else {
                    return Ok((q_idx, redone_match, kept_for_replay));
                };

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:289-318
                // No NCBI counterpart: position of the next Newton lanes prefetch.
                let mut x_next_prefetch = 0usize;
                for (match_position, local_match) in preliminary_hit_list
                    .hsplist_array
                    .iter()
                    .take(preliminary_hit_list.hsplist_count)
                    .enumerate()
                {
                    let s_idx = local_match.oid;
                    kappa_preliminary_hits.clear();
                    kappa_preliminary_hits.extend(
                        local_match
                            .hsps
                            .iter()
                            .map(preliminary_hit_from_local_match_hsp),
                    );
                    if kappa_preliminary_hits.is_empty() {
                        continue;
                    }

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529
                    // Tested here only for a query that is alone in its batch: FALSE for
                    // a query that is not the last one of its batch, and replayed after
                    // the parallel loop for the last one of several (see above).
                    if alone_in_batch
                        && blast_compo_early_termination(
                            local_match.best_evalue,
                            std::slice::from_ref(&redone_match),
                        )
                    {
                        continue;
                    }
                    // No NCBI counterpart: the Newton lanes solve the first Newton problems of this
                    // and the next matches ahead, four per vector (kappa::x_newton_lanes_prefetch_matches);
                    // a result is used only by a call with the same input bits. It does not change any value
                    // NCBI computes.
                    if match_position >= x_next_prefetch {
                        x_next_prefetch = match_position + X_NEWTON_LANES_MATCHES;
                        x_newton_lanes_prefetch_matches(
                            &preliminary_hit_list.hsplist_array[match_position
                                ..x_next_prefetch.min(preliminary_hit_list.hsplist_count)],
                            subjects,
                            std::slice::from_ref(if alone_in_batch {
                                &redone_match
                            } else {
                                &x_untouched_heap
                            }),
                            &query_workspace,
                            &composition_workspace,
                            &redo_align_params,
                            &mut kappa_subject_range_cache,
                        );
                    }
                    // Read from the database here, after
                    // the test (a query alone in its batch; the others ahead of the replay),
                    // as NCBI reads it, and after the first Newton problems, which read it
                    // ahead and drop it; an error is kept for the replay as the redo's own
                    // errors are (below).
                    let subject = match subjects.read(s_idx) {
                        Ok(subject) => subject,
                        Err(error) if replayed => {
                            kept_for_replay.push((match_position, Err(error)));
                            continue;
                        }
                        Err(error) => return Err(error),
                    };
                    let subject_raw = subject.residues();
                    let postprocessed = postprocess_preliminary_hits(
                        &kappa_preliminary_hits,
                        &query_workspace,
                        &mut composition_workspace,
                        query_nomask,
                        subject_raw,
                        args.scoring.matrix,
                        gapped_params,
                        gapped_gumbel.as_ref(),
                        &redo_align_params,
                        args.evalue,
                        &mut kappa_gap_scratch,
                        &mut kappa_subject_range_cache,
                        needs.seqalign,
                    );
                    if replayed {
                        // An error is kept as well: the one-thread loop reports it only
                        // for a match that the test does not end.
                        kept_for_replay.push((match_position, postprocessed));
                        continue;
                    }
                    let postprocessed = postprocessed?;
                    let best_score = postprocessed.best_score;
                    let best_evalue = postprocessed.best_evalue;
                    let Some(hsp_list) = postprocessed.hsp_list else {
                        continue;
                    };
                    if redone_match.would_insert(best_evalue, best_score, s_idx) {
                        redone_match.insert(hsp_list, best_evalue, best_score, s_idx);
                    }
                }

                Ok((q_idx, redone_match, kept_for_replay))
            })
            .collect()
    });
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3817-3850
    let mut kept_queries = Vec::new();
    for (q_idx, heap, kept_for_replay) in query_heaps? {
        redone_matches[q_idx] = heap;
        if !kept_for_replay.is_empty() {
            kept_queries.push((q_idx, kept_for_replay));
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529,3705-3719
    // The last query of each batch of several: the steps of `one_thread_redo`
    // in stream order, with the heaps of the batch and the results redone above.
    for (q_idx, kept_for_replay) in kept_queries {
        let batch = kappa_query_batches[q_idx].clone();
        let preliminary_hit_list = preliminary_hit_lists[q_idx]
            .as_ref()
            .expect("kept redo results belong to a query with preliminary hits");
        for (match_position, postprocessed) in kept_for_replay {
            let local_match = &preliminary_hit_list.hsplist_array[match_position];
            if blast_compo_early_termination(
                local_match.best_evalue,
                &redone_matches[batch.clone()],
            ) {
                continue;
            }
            let postprocessed = postprocessed?;
            let best_score = postprocessed.best_score;
            let best_evalue = postprocessed.best_evalue;
            let Some(hsp_list) = postprocessed.hsp_list else {
                continue;
            };
            let redone_match = &mut redone_matches[q_idx];
            if redone_match.would_insert(best_evalue, best_score, local_match.oid) {
                redone_match.insert(hsp_list, best_evalue, best_score, local_match.oid);
            }
        }
    }
    Ok(())
}

/// NCBI's error of a query batch whose Smith-Waterman redo (`kappa::smith_waterman::SwRedoFailed`)
/// fails at the batch's last redone match: NCBI's one-thread redo keeps the status of its last
/// match (`status_code_tld[0]`; a failed earlier match is dropped and the next redone match sets
/// the status again), a status other than 0 at the end fails the batch's traceback, and the
/// application writes `BLAST engine error: Traceback failed` (exit status 3) after the reports of
/// the batches before (without the epilog) and the failing batch's reader messages and search
/// strategy request. LOSAT gives NCBI's one-thread output at every `-num_threads`. The error names
/// the batch's first query (`batch_start`, a query of this call: the search of all queries, the only
/// call that runs `-use_sw_tback`, `x_query_block_plan`), as the error of a header file that a
/// batch reaches does (`database::HeaderReadError`), whose batches before `search` searches again
/// and reports alone.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3775-3779
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:251-253
/// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:274-276 (the exception leaves the
/// batch loop before the batch's formatting and `PrintEpilog`)
///
/// NCBI ends such a batch by `Traceback failed` only when the batch's first query holds no redone
/// match (`first_query`, the heap of `batch_start`). With a status other than 0 the clean-up frees
/// the arrays of every query's heap, leaving each heap's count, and then empties the heap of the
/// first query alone (`s_ClearHeap` of `redoneMatches_tld[i]`, the address of its element 0):
/// `BlastCompo_HeapPop` converts nothing (`array` is NULL) and, when the count is above 0, reads
/// `heapArray[1]` through the NULL `heapArray` (SIGSEGV, at address 0x28 on ws, 2026-11-18). A
/// heap's count is above 0 once `BlastCompo_HeapInsert` (blast_kappa.c:3705-3719) has put a match
/// in it, and nothing pops it before the clean-up, so NCBI crashes exactly when a match of the
/// batch's first query was redone and kept. The heaps of the other queries are freed without a
/// pop; the heaps of the function's own `redoneMatches` (blast_kappa.c:3006, 3119-3130), which its
/// clean-up empties at 3807-3813, receive no match (the loop inserts into its block's
/// `redoneMatches`, :3472 and :3496, that thread's heaps), so they pop nothing. LOSAT refuses that
/// batch (rule 3: NCBI crashes; PD-LOSAT-NCBI-DEFECTS 1.53), in the place of NCBI's error and
/// after the same reports (the S1 form of 2026-11-18: what NCBI formats before its failing step,
/// not the part of its stdio buffer that its crash loses). Exit status 1, its message.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3795-3805
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2518-2524
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:444-457
#[cold]
#[inline(never)]
fn sw_batch_failed(batch_start: usize, first_query: &BlastCompoHeap) -> anyhow::Error {
    let error = if first_query.n > 0 {
        anyhow::anyhow!(
            "a query batch whose last match fails its Smith-Waterman redo (-use_sw_tback) after a match of the batch's first query was redone and kept, with which NCBI BLAST+ 2.17.0's blastp crashes, is not supported by LOSAT's BLASTP"
        )
    } else {
        crate::blastinput::app::engine_error("Traceback failed")
    };
    super::database::HeaderReadError {
        query: batch_start,
        error,
    }
    .into()
}

/// LOSAT's refusal of a match whose Smith-Waterman redo does not end in NCBI
/// (`kappa::smith_waterman::SwRedoEndless`; PD-LOSAT-NCBI-DEFECTS 1.53, rule 3, within the class of
/// RPS-BLAST's refusal of 1.30). NCBI's search of the match's batch never ends after NCBI has
/// printed the reports of the batches before (`blastp_app.cpp:257-276`, cited at
/// [`sw_batch_failed`]): the refusal fails that batch as NCBI's own error of a batch's search
/// does, after those reports (without the epilog), as TBLASTX's and TBLASTN's refusals of a step
/// where NCBI's search of a batch does not finish do (`blastinput::app::ErrorAfterBatchReports`,
/// the Owner's decision of 2026-11-16 on line RG's Q-X2c: "the results of the query batches before
/// are written"). Exit status 1, its message.
#[cold]
#[inline(never)]
fn sw_endless(batch_start: usize) -> anyhow::Error {
    super::database::HeaderReadError {
        query: batch_start,
        error: anyhow::anyhow!(
            "a match whose Smith-Waterman redo (-use_sw_tback) does not end in NCBI BLAST+ 2.17.0 (the e-value of a score of 0 passes the threshold, the match has several preliminary HSPs and the first query letter scores the subject's first letter above 0) is not supported by LOSAT's BLASTP"
        ),
    }
    .into()
}

/// NCBI's Smith-Waterman redo at a list whose subject its fetch excludes (`mark`, a list that
/// `database::exclude_fetched` removed; `batch_heaps`, the heaps of the list's query batch): the
/// list is freed without a change of the status when the early termination test holds, and
/// otherwise the excluded fetch sets the status back to 0, which clears an earlier failure of
/// the batch (`sw_last_failed`). Nothing is fetched, redone or inserted.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3540
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3560-3575
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:306-314
#[inline]
fn fetch_excluded_reached(
    mark: &super::database::FetchExcludedMatch,
    batch_heaps: &[BlastCompoHeap],
    sw_last_failed: &mut bool,
) {
    if !blast_compo_early_termination(mark.best_evalue, batch_heaps) {
        *sw_last_failed = false;
    }
}

/// NCBI's Smith-Waterman redo at the lists of one query whose subjects its fetch excludes
/// (`marks`, in stream order), read after the query's loop. `status_position` is the place of
/// the query's last kept list whose redo set the status (`None` where no redo of the query did).
/// A list before that place changes nothing, since that redo sets the status again. After it no
/// list of the query enters a heap (each one ends early), so each later excluded list meets the
/// heaps as they are after the loop (`batch_heaps`) and goes through `fetch_excluded_reached`.
/// With no such redo the query changed no heap, and every excluded list goes through it.
///
/// No NCBI counterpart: NCBI tests each list in stream order in its loop. The status this gives
/// is NCBI's one-thread status of the batch.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3774-3779
fn query_fetch_excluded_reached(
    marks: &[super::database::FetchExcludedMatch],
    status_position: Option<usize>,
    batch_heaps: &[BlastCompoHeap],
    sw_last_failed: &mut bool,
) {
    let after_status = status_position.map_or(0, |at| at + 1);
    for mark in marks.iter().filter(|mark| mark.position >= after_status) {
        fetch_excluded_reached(mark, batch_heaps, sw_last_failed);
    }
}

/// The one-thread loop: the queries redone in order on the calling thread, each query's matches
/// in stream order, each query's compact list expanded when the loop reaches it and its preliminary
/// list freed after its redo (counted in `x_prelim_freed`). It is the only path of a build without
/// threads, and the path of a search where no parallel path is chosen.
///
/// NCBI's blastp traceback calls `Blast_RedoAlignmentCore_MT` with the search's thread count (at
/// least 1):
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1488,1492-1494
///
/// With one thread its parallel region runs on the calling thread, as this loop does:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3429-3431
pub(super) fn one_thread_redo<P>(
    redo: &RedoInputs<'_, P>,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    redone_matches: &mut [BlastCompoHeap],
    x_prelim_freed: &mut usize,
) -> Result<()>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams>,
{
    let RedoInputs {
        args,
        blastp_parallel_pool,
        needs,
        subjects,
        contexts,
        search_spaces,
        gapped_params,
        gapped_gumbel,
        kappa_query_batches,
        evalue_adjust_contexts,
        build_redo_align_params,
        fetch_excluded,
        ..
    } = *redo;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334
    let mut composition_workspace = BlastCompositionWorkspace::new_blosum62();
    if args.scoring.matrix != ScoringMatrix::Blosum62 {
        blastp_use_matrix_workspace(&mut composition_workspace, args.scoring.matrix);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3656-3667
    let mut kappa_preliminary_hits = Vec::new();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3017-3019
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1934-1942
    let mut kappa_gap_scratch = GapAlignScratch::new();
    // In a block search of the subjects encoded for all blocks, the store of SEG-masked subject
    // ranges of all blocks (`x_query_blocks::x_block_seg_store`): an entry depends on its key
    // only.
    let x_block_seg_store = x_query_blocks::x_block_seg_store(blastp_parallel_pool, subjects);
    let mut kappa_subject_range_cache =
        BlastpKappaSubjectRangeCache::new(x_block_seg_store.clone());
    // `-use_sw_tback` (item X4): the Smith-Waterman redo's memory, and whether the last match
    // that the redo of the query batch has reached failed (NCBI's `status_code_tld`).
    let mut sw_scratch = args
        .use_sw_tback
        .then(kappa::smith_waterman::SwMatchScratch::new);
    let mut sw_last_failed = false;
    // The marks of the lists that the Smith-Waterman redo reaches and whose fetch excludes the
    // subject (`fetch_excluded_reached`), by query; the X-drop redo reads none.
    let mut fetch_excluded = if sw_scratch.is_some() {
        fetch_excluded
    } else {
        &[][..]
    };

    for q_idx in 0..kappa_query_batches.len() {
        let batch = kappa_query_batches[q_idx].clone();
        if batch.start == q_idx {
            sw_last_failed = false;
        }
        let excluded_here = fetch_excluded
            .iter()
            .take_while(|mark| mark.query == q_idx)
            .count();
        let (excluded_here, excluded_after) = fetch_excluded.split_at(excluded_here);
        fetch_excluded = excluded_after;
        // The place of the query's last list whose Smith-Waterman redo set the status.
        let mut x_status_position = None;
        let redo_align_params = build_redo_align_params(q_idx)?;

        let ctx = &contexts[q_idx];
        let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
        let query_nomask = query_nomask_sequence(ctx);
        // This query's list of full HSPs (`x_prelim_compact.rs`).
        if let Some(compact) = x_compact_hit_lists[q_idx].take() {
            preliminary_hit_lists[q_idx] = Some(x_expand_hit_list(compact));
        }
        let Some(preliminary_hit_list) = preliminary_hit_lists[q_idx].as_ref() else {
            if sw_last_failed && q_idx + 1 == batch.end {
                return Err(sw_batch_failed(batch.start, &redone_matches[batch.start]));
            }
            continue;
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:821-847
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3244-3248
        let query_workspace = build_query_workspace(
            query_raw,
            ctx.aa_len as i32,
            search_spaces[q_idx].length_adjustment as i32,
            search_spaces[q_idx].effective_space,
            evalue_adjust_contexts[q_idx],
        );

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:133-151
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:289-318
        //
        // For composition-based searches the preliminary stream is read in
        // query order, and within each query from best to worst preliminary
        // HSPList. Use the sorted preliminary hitlist directly instead of
        // draining an unordered map, so `localMatch->best_evalue` and the
        // admission order match NCBI.
        // No NCBI counterpart: position of the next Newton lanes prefetch.
        let mut x_next_prefetch = 0usize;
        for (x_match_position, local_match) in preliminary_hit_list
            .hsplist_array
            .iter()
            .take(preliminary_hit_list.hsplist_count)
            .enumerate()
        {
            let s_idx = local_match.oid;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3580-3581
            kappa_preliminary_hits.clear();
            kappa_preliminary_hits.extend(
                local_match
                    .hsps
                    .iter()
                    .map(preliminary_hit_from_local_match_hsp),
            );
            // An empty match (NCBI's `hsp_array` not NULL with `hspcnt` 0) would reach NCBI's redo
            // and set `status_code_tld` to 0 (`-use_sw_tback`: the batch's failure status); LOSAT's
            // preliminary lists hold no match without HSPs, so this skip is never taken there.
            if kappa_preliminary_hits.is_empty() {
                continue;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3528
            // `redoneMatches` and `numQueries` are those of the query's batch
            // (`kappa_query_batches`).
            if blast_compo_early_termination(
                local_match.best_evalue,
                &redone_matches[batch.clone()],
            ) {
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529
            // No NCBI counterpart: the Newton lanes solve the first Newton problems of this
            // and the next matches ahead, four per vector (kappa::x_newton_lanes_prefetch_matches);
            // a result is used only by a call with the same input bits. It does not change any value
            // NCBI computes.
            if x_match_position >= x_next_prefetch {
                x_next_prefetch = x_match_position + X_NEWTON_LANES_MATCHES;
                x_newton_lanes_prefetch_matches(
                    &preliminary_hit_list.hsplist_array
                        [x_match_position..x_next_prefetch.min(preliminary_hit_list.hsplist_count)],
                    subjects,
                    &redone_matches[batch.clone()],
                    &query_workspace,
                    &composition_workspace,
                    &redo_align_params,
                    &mut kappa_subject_range_cache,
                );
            }
            // Read from the database here, after the test, as
            // NCBI reads it, and after the first Newton problems, which read it ahead and drop it.
            let subject = subjects.read(s_idx)?;
            let subject_raw = subject.residues();
            let postprocessed = if let Some(sw_scratch) = sw_scratch.as_mut() {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3644-3646,3735-3737
                // The Smith-Waterman redo (item X4). A failed redo leaves the match without HSPs
                // and its status until the next match that the redo reaches; a redo that does not
                // end in NCBI is refused.
                let redone_match = &mut redone_matches[q_idx];
                match kappa::postprocess_preliminary_hits_smith_waterman(
                    &kappa_preliminary_hits,
                    &query_workspace,
                    &mut composition_workspace,
                    query_nomask,
                    subject_raw,
                    args.scoring.matrix,
                    gapped_params,
                    gapped_gumbel.as_ref(),
                    &redo_align_params,
                    args.evalue,
                    &mut kappa_gap_scratch,
                    &mut kappa_subject_range_cache,
                    needs.seqalign,
                    sw_scratch,
                    &mut |evalue, score| redone_match.would_insert(evalue, score, s_idx),
                ) {
                    Ok(postprocessed) => {
                        sw_last_failed = false;
                        x_status_position = Some(x_match_position);
                        postprocessed
                    }
                    Err(error) if error.is::<kappa::smith_waterman::SwRedoFailed>() => {
                        sw_last_failed = true;
                        x_status_position = Some(x_match_position);
                        continue;
                    }
                    Err(error) if error.is::<kappa::smith_waterman::SwRedoEndless>() => {
                        return Err(sw_endless(batch.start));
                    }
                    Err(error) => return Err(error),
                }
            } else {
                postprocess_preliminary_hits(
                    &kappa_preliminary_hits,
                    &query_workspace,
                    &mut composition_workspace,
                    query_nomask,
                    subject_raw,
                    args.scoring.matrix,
                    gapped_params,
                    gapped_gumbel.as_ref(),
                    &redo_align_params,
                    args.evalue,
                    &mut kappa_gap_scratch,
                    &mut kappa_subject_range_cache,
                    needs.seqalign,
                )?
            };
            let best_score = postprocessed.best_score;
            let best_evalue = postprocessed.best_evalue;
            let Some(hsp_list) = postprocessed.hsp_list else {
                continue;
            };
            let redone_match = &mut redone_matches[q_idx];
            if redone_match.would_insert(best_evalue, best_score, s_idx) {
                redone_match.insert(hsp_list, best_evalue, best_score, s_idx);
            }
        }
        // The excluded lists of the query, read after its loop (none without `-use_sw_tback`, so
        // the loop above has no test of them).
        query_fetch_excluded_reached(
            excluded_here,
            x_status_position,
            &redone_matches[batch.clone()],
            &mut sw_last_failed,
        );
        // Nothing reads this query's preliminary list after its redo; it is freed here as NCBI
        // frees each list it has redone (`x_prelim_free.rs`).
        preliminary_hit_lists[q_idx] = None;
        *x_prelim_freed += 1;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3774-3781
        // A failed Smith-Waterman redo of the batch's last redone match fails NCBI's traceback of
        // the batch (`Traceback failed`, traceback_stage.cpp:251-253), or crashes NCBI's clean-up
        // when the batch's first query holds a redone match (`sw_batch_failed`).
        if sw_last_failed && q_idx + 1 == batch.end {
            return Err(sw_batch_failed(batch.start, &redone_matches[batch.start]));
        }
    }
    // EXPERIMENT (LOSAT_X_BPQBLOCK): the block search keeps the store for the next block only as
    // LOSAT_X_BPSEGWIN keeps it for the next window (`x_redo_par.rs` `x_redo`).
    if let Some(store) = &x_block_seg_store {
        let seg_win = kappa::x_bp_seg_win();
        let emptied = seg_win.is_none_or(|min| store.x_len() >= min);
        if emptied {
            store.x_clear();
        }
        // The switch's own side here is a store kept for the next block (`=0` empties it after
        // every block, as each block of the search before the store).
        if seg_win.is_some() {
            crate::utils::xreach::x_reach!("LOSAT_X_BPSEGWIN", !emptied);
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "redo_run_tests.rs"]
mod tests;
