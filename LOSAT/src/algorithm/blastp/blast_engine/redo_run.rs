//! The composition-based redo of a BLASTP call's matches (`run_call`, after the choice of the
//! ordinary traceback): the parameters of each query's redo (`RedoParamsBuilder`, which
//! `run_call`'s `build_redo_align_params` calls), and the redo's three parallel paths, which read
//! `ParallelRedo`: the windows of `LOSAT_X_BPREDOPAR` (`x_windows_redo`), the matches of one query
//! on the pool (`match_parallel_redo`) and the queries on the pool (`query_parallel_redo`). Each path
//! is the statements of `run_call` at its place; `run_call` calls the one that its plan of the
//! redo chose, and its `?` passes the path's errors on. The plan and the one-thread loop stay in
//! `run_call`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2977-2981
//! ```c
//!  *  Recompute alignments for each match found by the gapped BLAST
//!  *  algorithm.
//!  */
//! Int2
//! Blast_RedoAlignmentCore_MT(EBlastProgramType program_number,
//! ```
//!
//! This module is a child of `blast_engine` for the private items of that file it reads: the six
//! `blastp_*` functions that `build` calls (`blastp_local_scaling_factor`,
//! `blastp_compo_adjust_mode`, `blastp_redo_matrix`, `blastp_redo_near_identical_cutoff`,
//! `blastp_redo_x_dropoff`, `blastp_redo_cutoff_score`), `blastp_use_matrix_workspace`,
//! `preliminary_hit_from_local_match_hsp`, `query_nomask_sequence`, `X_NEWTON_LANES_MATCHES`, and
//! five of its imports that only this module reads (`build_matrix_info`,
//! `BlastCompoGappingParams`, `BlastRedoAlignParams`, `BlastpPostprocessResult`, `BLASTAA_SIZE`).

use std::sync::OnceLock;

use anyhow::Result;

use super::{
    blastp_compo_adjust_mode, blastp_local_scaling_factor, blastp_redo_cutoff_score,
    blastp_redo_matrix, blastp_redo_near_identical_cutoff, blastp_redo_x_dropoff,
    build_matrix_info, BlastCompoGappingParams, BlastRedoAlignParams, BLASTAA_SIZE,
};
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::core::composition_adjustment::adjust_scores::BlastMatrixInfo;
use crate::stats::KarlinParams;

// The parallel paths' imports.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use {
    super::{
        blastp_use_matrix_workspace, preliminary_hit_from_local_match_hsp, query_nomask_sequence,
        BlastpPostprocessResult, X_NEWTON_LANES_MATCHES,
    },
    crate::algorithm::blastp::gapalign::GapAlignScratch,
    crate::algorithm::blastp::hsp::{blast_compo_early_termination, BlastCompoHeap, BlastpHitList},
    crate::algorithm::blastp::kappa::{
        build_query_workspace, postprocess_preliminary_hits, x_newton_lanes_prefetch_matches,
        BlastpKappaSubjectRangeCache,
    },
    crate::algorithm::blastp::x_prelim_compact::XPrelimHsp,
    crate::algorithm::blastp::x_redo_par::{self, XRedoParSearch},
    crate::algorithm::blastp::x_redo_subjects::XRedoSubjects,
    crate::config::ScoringMatrix,
    crate::core::composition_adjustment::adjust_scores::BlastCompositionWorkspace,
    crate::report::record::ReportNeeds,
    crate::stats::search_space::SearchSpace,
    crate::stats::BlastGumbelBlk,
    crate::utils::threading::SearchPool,
    rayon::prelude::*,
    std::ops::Range,
};

/// What the parameters of a query's redo are made of (`build`): each field is `run_call`'s value
/// of the same name (`Copy` scalars by value, the rest by reference, so that the moved body reads
/// them as it did in `run_call`), and the start matrix that the
/// parameters of every query share, made at the first `build`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2401-2407
/// ```c
/// /**
///  * Read the parameters required for the Blast_RedoOneMatch* functions from
///  * the corresponding parameters in standard BLAST datatypes.  Return a new
///  * object representing these parameters.
///  */
/// static Blast_RedoAlignParams *
/// s_GetAlignParams(BlastKappa_GappingParamsContext * context,
/// ```
pub(super) struct RedoParamsBuilder<'a> {
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) contexts: &'a [QueryContext],
    pub(super) gapped_params: &'a KarlinParams,
    pub(super) hit_cutoff_score_min: i32,
    pub(super) ideal_ungapped_lambda: f64,
    pub(super) x_drop_gapped_final: i32,
    // NCBI reference: c++/src/algo/blast/core/blast_kappa.c:2228-2236
    // self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
    // status = s_GetStartFreqRatios(self->startFreqRatios, matrixName);
    // if (status == 0) {
    //     Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
    //                              self->startFreqRatios, self->ungappedLambda);
    // }
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
    /// ```c
    ///     } else {
    ///         return
    ///             Blast_RedoAlignParamsNew(&scaledMatrixInfo, &gapping_params,
    ///                                      compo_adjust_mode, positionBased,
    ///                                      query_is_translated,
    ///                                      subject_is_translated,
    ///                                      queryInfo->max_length, cutoff_s, cutoff_e,
    ///                                      do_link_hsps, &redo_align_callbacks,
    ///                                      near_identical_cutoff);
    /// ```
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
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2216-2231
                    // ```c
                    // self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
                    // status = s_GetStartFreqRatios(self->startFreqRatios, matrixName);
                    // Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
                    //                          self->startFreqRatios, self->ungappedLambda);
                    // ```
                    ideal_ungapped_lambda / local_scaling_factor,
                )
            })
            .as_ref()
            .map_err(|error| anyhow::anyhow!("{error}"))?
            .clone();
        Ok(BlastRedoAlignParams {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2216-2231
            // ```c
            // self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
            // status = s_GetStartFreqRatios(self->startFreqRatios, matrixName);
            // Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
            //                          self->startFreqRatios, self->ungappedLambda);
            // ```
            matrix_info,
            // NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:97-106
            // ```c
            // typedef struct BlastCompo_GappingParams {
            //     int gap_open;
            //     int gap_extend;
            //     int decline_align;
            //     int x_dropoff;
            //     void * context;
            // } BlastCompo_GappingParams;
            // ```
            //
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2117-2133
            // ```c
            // kbp->Lambda /= scale_factor;
            // sp->gap_open = (Int4)BLAST_Nint(sp->gap_open  * scale_factor);
            // sp->gap_extend = (Int4)BLAST_Nint(sp->gap_extend * scale_factor);
            // ```
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
            // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1107-1118
            // ```c
            // Blast_RedoOneMatch(..., int ** matrix, int alphsize,
            //                    ..., int compositionTestIndex, ...)
            // ```
            alphsize: BLASTAA_SIZE as i32,
            composition_test_index: i32::from(args.comp_based_stats.unified_p),
            log_k: gapped_params.k.ln(),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3584-3690
            // ```c
            // localScalingFactor      /* i */
            // ...
            // s_HSPListNormalizeScores(hsp_list, kbp->Lambda, kbp->logK,
            //                          localScalingFactor);
            // ```
            score_divisor: local_scaling_factor,
            smith_waterman: false,
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
            // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:382-406
            // ```c
            // Nearly identical alignments are computed with exact subjects,
            // and others with segged subjects; this makes comparing end
            // points more difficult.
            // ```
        })
    }
}

/// What the redo's three parallel paths read of `run_call`: each field is `run_call`'s value of
/// the same name (the `Copy` scalars by value, the rest by reference, so that the moved bodies read
/// them as they did in `run_call`; `subjects` are the redo's subjects,
/// `build_redo_align_params` makes the parameters of a query's redo with `RedoParamsBuilder`, and
/// `x_redo_par_search` is what the windows read). The paths' threads read these values, as NCBI's
/// threads share those of its parallel region.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3432-3445
/// ```c
///     shared(interrupt, seqsrc_tld, score_params_tld, hit_params_tld, \
///     gap_align_tld, results_tld, \
///     redoneMatches_tld, \
///     STDERR_COMMA \
///     numMatches, theseMatches, \
///     numFrames, program_number, subjectBlk_tld, positionBased, \
///     default_db_genetic_code, localScalingFactor, queryInfo, \
///     sbp, smithWaterman, numQueries_tld, compositionTestIndex_tld, forbidden_tld, \
///     NRrecord_tld, actual_num_threads, sbp_tld, \
///     matrix_tld, query_info_tld, numContexts_tld, \
///     genetic_code_string, queryBlk, compo_adjust_mode, \
///     alignments_tld, incoming_align_set_tld, savedParams_tld, \
///     scoringParams, redo_align_params_tld, \
///     status_code_tld)
/// ```
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) struct ParallelRedo<'a, P> {
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) blastp_parallel_pool: &'a SearchPool<'a>,
    pub(super) needs: ReportNeeds,
    pub(super) subjects: XRedoSubjects<'a>,
    pub(super) contexts: &'a [QueryContext],
    pub(super) search_spaces: &'a [SearchSpace],
    pub(super) gapped_params: &'a KarlinParams,
    pub(super) gapped_gumbel: &'a Option<BlastGumbelBlk>,
    pub(super) kappa_query_batches: &'a [Range<usize>],
    pub(super) kappa_parallel_query_indices: &'a [usize],
    pub(super) kappa_match_parallel_query_index: Option<usize>,
    pub(super) x_redo_par_search: &'a XRedoParSearch<'a, P>,
    pub(super) build_redo_align_params: &'a P,
    pub(super) x_win_expand: bool,
}

/// The redo in windows of `LOSAT_X_BPREDOPAR` (`x_redo_par::x_redo_par_run`): the matches of all
/// queries on the pool in jobs, replayed in stream order, the compact lists expanded as the
/// windows reach them (`x_win_expand`); the lists it has redone are released, and counted in
/// `x_prelim_freed`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3414-3418
/// ```c
///         /*
///          * Convert linked list of matches into array.
///          */
///         theseMatches =
///                 (BlastHSPList**) calloc(numMatches, sizeof(BlastHSPList*));
/// ```
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn x_windows_redo<P>(
    redo: &ParallelRedo<'_, P>,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    redone_matches: &mut [BlastCompoHeap],
    x_prelim_freed: &mut usize,
) -> Result<()>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let ParallelRedo {
        x_redo_par_search,
        x_win_expand,
        ..
    } = *redo;
    // The parallel paths read the expanded lists of all queries (`run_call` expands them before
    // the paths), or the windows expand the compact lists as they reach them.
    debug_assert!(x_win_expand || x_compact_hit_lists.iter().all(Option::is_none));
    let released = x_redo_par::x_redo_par_run(
        x_redo_par_search,
        preliminary_hit_lists,
        x_win_expand.then_some(&mut x_compact_hit_lists[..]),
        redone_matches,
    )?;
    *x_prelim_freed += released;
    Ok(())
}

/// The redo of the matches of one query on the pool (`kappa_match_parallel_query_index`): the matches
/// are redone on the pool, each work unit with its own scratch, and the early termination test and
/// the heap insertions are replayed after them in stream order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3511-3511
/// ```c
///                 BlastHSPList* localMatch = theseMatches[b];
/// ```
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn match_parallel_redo<P>(
    redo: &ParallelRedo<'_, P>,
    preliminary_hit_lists: &[Option<BlastpHitList>],
    redone_matches: &mut [BlastCompoHeap],
) -> Result<()>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let ParallelRedo {
        args,
        blastp_parallel_pool,
        needs,
        subjects,
        contexts,
        search_spaces,
        gapped_params,
        gapped_gumbel,
        kappa_query_batches,
        kappa_match_parallel_query_index,
        build_redo_align_params,
        ..
    } = *redo;
    let q_idx = kappa_match_parallel_query_index
        .expect("single-query Kappa redo parallel path requires a query index");
    let preliminary_hit_list = preliminary_hit_lists[q_idx]
        .as_ref()
        .expect("single-query Kappa redo parallel path requires preliminary hits");
    let local_matches = &preliminary_hit_list.hsplist_array[..preliminary_hit_list.hsplist_count];
    let ctx = &contexts[q_idx];
    let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
    let query_nomask = query_nomask_sequence(ctx);

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3429-3449
    // ```c
    // #pragma omp parallel ... num_threads(actual_num_threads)
    // #pragma omp for schedule(static)
    // for (b = 0; b < numMatches; ++b) {
    // ```
    //
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3493-3503
    // ```c
    // scoringParams = score_params_tld[tid];
    // redoneMatches = redoneMatches_tld[tid];
    // NRrecord = NRrecord_tld[tid];
    // redo_align_params = redo_align_params_tld[tid];
    // ```
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3429-3459
    // ```c
    // #pragma omp parallel default(none) num_threads(actual_num_threads)
    // #pragma omp for schedule(static)
    // for (b = 0; b < numMatches; ++b) {
    // ```
    // NCBI reference: c++/src/algo/blast/core/blast_kappa.c:2309-2336,3493-3504
    // query_info->seq.data = &query_data[query_info->origin];
    // query_info = query_info_tld[tid];
    // Only immutable query preparation is shared within this search.
    let query_workspace = std::sync::OnceLock::new();
    let pool = blastp_parallel_pool;
    let precomputed: Vec<Result<Option<BlastpPostprocessResult>>> = pool.install(|| {
        local_matches
            .par_iter()
            // NCBI reference: c++/src/algo/blast/core/blast_kappa.c:3329-3334,3493-3503
            // NRrecord_tld[i] = Blast_CompositionWorkspaceNew();
            // NRrecord = NRrecord_tld[tid];
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
                        )
                    });
                    let mut kappa_subject_range_cache = BlastpKappaSubjectRangeCache::new(None);

                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3577-3585
                    // ```c
                    // *pStatusCode = s_ResultHspToDistinctAlign(...);
                    // ```
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

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3525-3529
    // ```c
    // if (BlastCompo_EarlyTermination(localMatch->best_evalue,
    //                                 redoneMatches,
    //                                 numQueries)) {
    // ```
    //
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3817-3849
    // ```c
    // /* Reduce results from all threads and continue with business as usual */
    // local_results = SThreadLocalDataArrayConsolidateResults(thread_data);
    // ```
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

/// The redo of the queries on the pool (`kappa_parallel_query_indices`), each query with its own
/// heap; the last query of a batch of several keeps its results, whose early termination test and
/// heap insertions are replayed after the queries in stream order, with the heaps of its batch.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3282-3283,3288-3291
/// ```c
///         redoneMatches_tld[i] =
///                 (BlastCompo_Heap*) calloc(numQueries, sizeof(BlastCompo_Heap));
/// ...
///         for (query_index = 0; query_index < numQueries; query_index++) {
///             status_code =
///                 BlastCompo_HeapInitialize(&redoneMatches_tld[i][query_index],
///                                           hitParams->options->hitlist_size,
/// ```
/// NCBI's threads split the matches, each thread with its own heaps; LOSAT runs the queries on the
/// pool, each query with its own heap.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn query_parallel_redo<P>(
    redo: &ParallelRedo<'_, P>,
    preliminary_hit_lists: &[Option<BlastpHitList>],
    redone_matches: &mut [BlastCompoHeap],
) -> Result<()>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let ParallelRedo {
        args,
        blastp_parallel_pool,
        needs,
        subjects,
        contexts,
        search_spaces,
        gapped_params,
        gapped_gumbel,
        kappa_query_batches,
        kappa_parallel_query_indices,
        build_redo_align_params,
        ..
    } = *redo;
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3119-3128
    // ```c
    // redoneMatches = calloc(numQueries, sizeof(BlastCompo_Heap));
    // for (query_index = 0;  query_index < numQueries;  query_index++) {
    //     BlastCompo_HeapInitialize(&redoneMatches[query_index], ...);
    // }
    // ```
    //
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3329-3334
    // ```c
    // NRrecord_tld[i] = Blast_CompositionWorkspaceNew();
    // status_code = Blast_CompositionWorkspaceInit(...);
    // ```
    //
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3436-3459
    // ```c
    // #pragma omp for schedule(static)
    // for (b = 0; b < numMatches; ++b) {
    //     ... redoneMatches = redoneMatches_tld[tid];
    // ```
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3429-3459
    // ```c
    // #pragma omp parallel default(none) num_threads(actual_num_threads)
    // #pragma omp for schedule(static)
    // for (b = 0; b < numMatches; ++b) {
    // ```
    let pool = blastp_parallel_pool;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636,3448-3449
    // ```c
    //                 status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
    //                                              subject_maybe_biased);
    // ...
    // #pragma omp for schedule(static)
    //         for (b = 0; b < numMatches; ++b) {
    // ```
    // Dispatch point: NCBI SEGs the subject range inside the per-match loop (omp for over
    // matches). Here the masked ranges are kept in one store for all queries, as the serial
    // redo already does. Reuse of an identical result only.
    // One store of masked subject ranges for all queries, as in `run_call`'s one-thread loop.
    let x_shared_ranges = Some(std::sync::Arc::new(
        crate::algorithm::blastp::kappa::XSharedSubjectRanges::new(),
    ));
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529
    // ```c
    //                 if (BlastCompo_EarlyTermination(
    //                         localMatch->best_evalue,
    //                         redoneMatches,
    //                         numQueries
    //                 )) {
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1580
    // ```c
    // BlastCompo_EarlyTermination(double evalue,
    //                             BlastCompo_Heap significantMatches[],
    //                             int numQueries)
    // {
    //     int i;
    //     for (i = 0;  i < numQueries;  i++) {
    //         if (BlastCompo_HeapFilledToCutoff(&significantMatches[i])) {
    //             double ecutoff = significantMatches[i].ecutoff;
    //             /* Only matches with evalue <= ethresh will be saved. */
    //             if (evalue <= EVALUE_STRETCH * ecutoff) {
    //                 ...
    //                 return FALSE;
    //             }
    //         } else {
    //             return FALSE;
    //         }
    //     }
    //     return TRUE;
    // }
    // ```
    // The test reads the heaps of the queries of the batch (`kappa_query_batches`). The
    // `run_call`'s one-thread loop redoes the queries in order, so when it tests a match of one
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
    let x_untouched_heap =
        BlastCompoHeap::new(args.traceback_hitlist_size(), args.inclusion_ethresh);
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
                let mut redone_match =
                    BlastCompoHeap::new(args.traceback_hitlist_size(), args.inclusion_ethresh);
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
                );
                let Some(preliminary_hit_list) = preliminary_hit_lists[q_idx].as_ref() else {
                    return Ok((q_idx, redone_match, kept_for_replay));
                };

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:289-318
                // ```c
                // hit_list = results->hitlist_array[index];
                // last_hsplist_index = hit_list->hsplist_count - 1;
                // *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
                // ```
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
                    // ```c
                    //                 if (BlastCompo_EarlyTermination(
                    //                         localMatch->best_evalue,
                    //                         redoneMatches,
                    //                         numQueries
                    //                 )) {
                    // ```
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
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3817-3850
    // ```c
    // /* Reduce results from all threads and continue with business as usual */
    // local_results = SThreadLocalDataArrayConsolidateResults(thread_data);
    // ```
    let mut kept_queries = Vec::new();
    for (q_idx, heap, kept_for_replay) in query_heaps? {
        redone_matches[q_idx] = heap;
        if !kept_for_replay.is_empty() {
            kept_queries.push((q_idx, kept_for_replay));
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529,3705-3719
    // ```c
    //                 if (BlastCompo_EarlyTermination(
    //                         localMatch->best_evalue,
    //                         redoneMatches,
    //                         numQueries
    //                 )) {
    // ...
    //                     if (BlastCompo_HeapWouldInsert(
    //                             &redoneMatches[query_index],
    //                             best_evalue,
    //                             best_score,
    //                             localMatch->oid
    //                     )) {
    //                         *pStatusCode =
    //                                 BlastCompo_HeapInsert(
    //                                         &redoneMatches[query_index],
    //                                         hsp_list,
    //                                         best_evalue,
    //                                         best_score,
    //                                         localMatch->oid,
    //                                         &discarded_aligns
    //                                 );
    // ```
    // The last query of each batch of several: the steps of `run_call`'s one-thread loop
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
