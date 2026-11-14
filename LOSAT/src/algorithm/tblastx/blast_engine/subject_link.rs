//! The end of a TBLASTX subject's search after its frame loop, up to its linked HSP list: the end
//! of its seed stage (`end_subject_seed_stage`), the frame bases of its HSP list
//! (`subject_frame_bases`), the linking call (`link_subject_hits`, twice), the sort and the reap by
//! the preliminary E-value after the first one (`sort_and_reap_after_first_link`), the
//! re-evaluation with the subject's letters (`reevaluate_subject_hits`) and the second linking
//! call (`relink_subject_hits`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-899 (in
//! `s_BlastSearchEngineCore`: `BLAST_LinkHsps`, then `s_Blast_HSPListReapByPrelimEvalue`) and
//! :1492-1520 (in `BLAST_PreliminarySearchEngine`: `Blast_HSPListReevaluateUngapped`, then
//! `BLAST_LinkHsps` again); c++/src/algo/blast/core/link_hsps.c:1761-1813 (`BLAST_LinkHsps`). The
//! end of the seed stage has no NCBI counterpart: NCBI frees its tables after the subject loop.
//!
//! This module is a child of `run_impl` for the modules it keeps private (`x_seed_release`,
//! `x_batch_pipe`, `x_search_timing`). Its functions are statements of `process_subject` moved
//! here in their order. The subject's guards stay in `process_subject`, which lends them
//! (`x_pipe_seed`, `x_lap`) and keeps the helped mode's guard that `end_subject_seed_stage`
//! returns until the subject ends. The HSP lists go in and come out by value; what the functions
//! read besides them is a struct of references and `Copy` values.

use super::x_batch_pipe::XSeedGuard;
use super::x_search_timing::XLap;
use super::x_seed_release::{self, XSeedStage};
use super::{
    apply_sum_stats_even_gap_linking_with_parallel, link_one_hit_with_int4_ends,
    reevaluate_ungapped_hsp_list, sort_ungapped_hits_by_score_ncbi, stage_dump, trace_hsp_target,
    trace_ungapped_hit_if_match, ungapped_hits_is_sorted_by_score_ncbi, AtomicOrdering, AtomicU64,
    DiagStruct, DiagnosticCounters, Instant, KarlinParams, LinkCall, LinkingParams, QueryContext,
    QueryFrame, TbxScore, UngappedHit, XStage,
};
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
use crate::algorithm::tblastx::sum_stats_linking::XHelped;
use crate::algorithm::tblastx::x_link_timing::XLinkTimer;

// No NCBI counterpart: the seed release and the batch pipe decide when tables are freed and where
// the linking calls run, not a value NCBI computes (see the note in the body).
/// Dispatch point 5 of the seed release: ends the subject's lookup stage (`x_seed_stage`; when its
/// end frees the batch's lookup table, the worker's diagonal table `diag_array` is dropped too)
/// and the subject's part of the batch pipe's seed-stage end with its `hsps` HSPs (when
/// `use_parallel`; otherwise `x_pipe_seed` ends with none when it drops), then sets the
/// thread's helped mode for the linking calls (`x_on_board`: the batch is on the search's board).
/// Returns whether the subject's linking calls may use the pool, and the helped mode's guard.
pub(super) fn end_subject_seed_stage(
    x_seed_stage: XSeedStage<'_, BlastAaLookupTable>,
    diag_array: &mut Vec<DiagStruct>,
    x_pipe_seed: &mut XSeedGuard<'_>,
    hsps: usize,
    use_parallel: bool,
    x_on_board: bool,
) -> (bool, XHelped) {
    // Dispatch point 5 of the seed release (per subject): the seed stage of this subject ends
    // here; nothing after it reads the batch's lookup table (the closure moves the stage that
    // its `lookup_ref` borrows from into this fn, and `end` consumes it; the unit path's stage
    // ended at dispatch point 3) or the worker's diagonal table. The batch's last lookup stage
    // frees the lookup table and returns true, and this worker, which gets no further subject,
    // drops its diagonal table before the subject's first linking call, where the heap of the
    // search peaks.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873
    // ```c
    // if (hit_params->link_hsp_params) {
    //     status = BLAST_LinkHsps(program_number, hsp_list_out, query_info,
    //               subject->length, gap_align->sbp, hit_params->link_hsp_params,
    //               score_options->gapped_calculation);
    // ```
    // NCBI links with its lookup table and diagonal table still allocated (freed after the
    // subject loop, blast_engine.c:1605); `BLAST_LinkHsps` reads neither.
    if x_seed_stage.end() {
        x_seed_release::drop_worker_table(diag_array);
    }
    // EXPERIMENT (LOSAT_X_TBXPIPE): the seed stage of this subject has ended with the HSP list
    // that the first linking call takes; the batch's last one may start the next batch. While a
    // later batch may run, the subject's two linking calls do not use the pool (`x_link_on_pool`,
    // tblastx/x_batch_pipe.rs item 6); the lists are the same with or without it.
    let x_link_on_pool = use_parallel && x_pipe_seed.end(hsps);
    // EXPERIMENT (LOSAT_X_TBXBOARD): a chain that links without the pool puts its groups and
    // sorts on the search's board instead of running them on this thread
    // (tblastx/x_link_slices.rs `x_helped`); the guard ends with the subject.
    let x_helped = crate::algorithm::tblastx::sum_stats_linking::x_link_helped(
        use_parallel && !x_link_on_pool && x_on_board,
    );
    (x_link_on_pool, x_helped)
}

/// The start of each subject frame of `s_frames_listed` in the subject's translation as NCBI
/// concatenates it, which the linking calls read.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1098-1101 (in
/// `BLAST_GetAllTranslations`, :1045).
pub(super) fn subject_frame_bases(s_frames_listed: &[QueryFrame]) -> Vec<i32> {
    // Build NCBI-style subject frame base offsets for sum-statistics linking.
    // In NCBI, HSP coords live in a concatenated translation buffer with sentinels.
    // LOSAT uses per-frame sequences; for linking we emulate absolute offsets by
    // concatenating frames in the same order as `generate_frames()`.
    let mut subject_frame_bases: Vec<i32> = Vec::with_capacity(s_frames_listed.len());
    let mut base: i32 = 0;
    for f in s_frames_listed {
        subject_frame_bases.push(base);
        // NCBI concatenation shares the trailing NULLB sentinel between frames:
        //   offset += length + 1;
        // where `length` is the number of residues (excluding sentinels).
        // Source: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
        base += f.aa_seq.len() as i32 - 1;
    }
    subject_frame_bases
}

/// What both linking calls of a subject read besides its HSP list, made per subject by
/// `process_subject`: the linking cutoffs' Karlin block, the subject's linking parameters, the
/// contexts with their length adjustments and search spaces, the frame bases, and whether the
/// calls may use the pool.
pub(super) struct SubjectLinkArgs<'a> {
    pub(super) linking_params_for_cutoff: &'a KarlinParams,
    pub(super) linking_params: &'a LinkingParams,
    pub(super) contexts_ref: &'a [QueryContext],
    pub(super) subject_frame_bases: &'a [i32],
    /// NCBI Parity: Pass pre-computed length_adjustment and eff_searchsp per context.
    /// These values are stored in `query_info->contexts[ctx]` in NCBI and referenced
    /// by link_hsps.c for BLAST_SmallGapSumE/BLAST_LargeGapSumE calculations.
    pub(super) length_adj_per_context: &'a [i64],
    /// See `length_adj_per_context`.
    pub(super) eff_searchsp_per_context: &'a [i64],
    pub(super) x_link_on_pool: bool,
}

/// The subject's HSP list `hits` linked with sum statistics (even gaps).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1761-1813 (`BLAST_LinkHsps`),
/// called at c++/src/algo/blast/core/blast_engine.c:871-873 and :1516-1520.
pub(super) fn link_subject_hits(
    link: &SubjectLinkArgs<'_>,
    hits: Vec<UngappedHit>,
) -> Vec<UngappedHit> {
    let &SubjectLinkArgs {
        linking_params_for_cutoff,
        linking_params,
        contexts_ref,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        x_link_on_pool,
    } = link;
    apply_sum_stats_even_gap_linking_with_parallel(
        hits,
        linking_params_for_cutoff,
        linking_params,
        contexts_ref,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        x_link_on_pool,
    )
}

// No NCBI counterpart: which form of the first linking call a subject of the one-hit word finder
// (`-window_size 0`) takes. NCBI links every list the same way; LOSAT's default kernel takes HSPs
// whose end is after their start (`link_one_hit_with_int4_ends` in `linking.rs`).
/// Whether the subject's list (of the one-hit word finder) holds an HSP of length 0 or less.
#[inline(never)]
pub(super) fn holds_int4_inverted(hits: &[UngappedHit]) -> bool {
    hits.iter().any(UngappedHit::is_int4_inverted)
}

/// The first linking call of a subject of the one-hit word finder whose list holds an HSP of
/// length 0 or less (`holds_int4_inverted`), with NCBI's `Int4` ends; those HSPs are dropped after
/// it, where NCBI's re-evaluation deletes them (`link_one_hit_with_int4_ends`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873 (the first
/// `BLAST_LinkHsps`) and c++/src/algo/blast/core/blast_hits.c:2703-2731
/// (`Blast_HSPListReevaluateUngapped`, called at blast_engine.c:1492-1497).
#[cold]
#[inline(never)]
pub(super) fn link_subject_hits_int4(
    link: &SubjectLinkArgs<'_>,
    hits: Vec<UngappedHit>,
) -> Vec<UngappedHit> {
    let &SubjectLinkArgs {
        linking_params_for_cutoff,
        linking_params,
        contexts_ref,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        x_link_on_pool,
    } = link;
    link_one_hit_with_int4_ends(
        hits,
        linking_params_for_cutoff,
        linking_params,
        contexts_ref,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        x_link_on_pool,
    )
}

/// What the stages of a subject after its first linking call read of the batch, made per subject
/// by `process_subject` from the names it captures: the contexts and their cutoff scores, the
/// E-value threshold, and the flags and counters of LOSAT_TIMING, LOSAT_DIAGNOSTICS and
/// LOSAT_DEBUG_OUTPUT_FILTER.
pub(super) struct SubjectStageArgs<'a> {
    pub(super) contexts_ref: &'a [QueryContext],
    pub(super) cutoff_scores: &'a [i32],
    pub(super) evalue_threshold: f64,
    pub(super) timing_enabled: bool,
    pub(super) reeval_ns: &'a AtomicU64,
    pub(super) reeval_calls: &'a AtomicU64,
    pub(super) linking_ns: &'a AtomicU64,
    pub(super) linking_calls: &'a AtomicU64,
    pub(super) diag_enabled: bool,
    pub(super) diagnostics: &'a DiagnosticCounters,
    pub(super) debug_output_filter: bool,
}

/// The first linking call's list sorted by score (when it is not), then reaped by the preliminary
/// E-value, with its stage dumps and the `[DEBUG PRELIM_REAP]` line of LOSAT_DIAGNOSTICS.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1802-1803 (the sort that ends
/// `BLAST_LinkHsps`) and c++/src/algo/blast/core/blast_engine.c:643-674
/// (`s_Blast_HSPListReapByPrelimEvalue`), called at :899.
pub(super) fn sort_and_reap_after_first_link(
    mut prelinked_ungapped_hits: Vec<UngappedHit>,
    stage: &SubjectStageArgs<'_>,
    x_lap: &mut XLap<'_>,
) -> Vec<UngappedHit> {
    let &SubjectStageArgs {
        evalue_threshold,
        diag_enabled,
        ..
    } = stage;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c:1802-1803
    // ```c
    //     /* Sort the HSP array by score */
    //     Blast_HSPListSortByScore(hsp_list);
    // ```
    // BLAST_LinkHsps ends by sorting the list it linked, here with the scores of the
    // preliminary (compressed-subject) search. The re-evaluation that follows may change or
    // trim those scores (random bases for ambiguous subject letters, two-hit
    // extension tails), after which Blast_HSPListReevaluateUngapped sorts again
    // (blast_hits.c:2733-2734): HSPs that tie after re-evaluation keep the order
    // of this first sort, not the order of the linking chains. That order is the
    // input order of the second BLAST_LinkHsps (blast_engine.c:1515-1520), whose
    // `>=` chain choices (link_hsps.c:613-622, 759, 887) give the better chain to
    // the later of two comparator-equal HSPs of different frames.
    // (The sort at the end of the second BLAST_LinkHsps is the score sort of each
    // subject's list in `hit_collector::QueryLists::add`.)
    if !ungapped_hits_is_sorted_by_score_ncbi(&prelinked_ungapped_hits) {
        sort_ungapped_hits_by_score_ncbi(&mut prelinked_ungapped_hits);
    }
    if stage_dump::enabled() {
        stage_dump::dump_ungapped_hits(
            "after_prelim_link_hsps_before_reevaluate",
            &prelinked_ungapped_hits,
        );
    }
    let prelim_linked_count = prelinked_ungapped_hits.len();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:654-662
    // ```c
    //    cutoff = hit_params->prelim_evalue;
    // ...
    //       if (hsp->evalue > cutoff) {
    // ```
    // In `s_Blast_HSPListReapByPrelimEvalue` (blast_engine.c:643, called at :899), whose
    // cutoff is `-evalue` without composition-based statistics (blast_parameters.c:868).
    // An HSP is removed when its e-value is greater than the cutoff, so a NaN cutoff
    // (`-evalue -nan`) keeps every HSP.
    #[expect(
        clippy::neg_cmp_op_on_partial_ord,
        reason = "NCBI removes an HSP only when evalue > cutoff (blast_engine.c:662); `<=` would drop NaN values"
    )]
    prelinked_ungapped_hits.retain(|h| !(h.e_value > evalue_threshold));
    x_lap.lap(XStage::Reap);
    if diag_enabled && prelim_linked_count != prelinked_ungapped_hits.len() {
        eprintln!(
            "[DEBUG PRELIM_REAP] linked_before={} kept={} filtered_by_prelim_evalue={} threshold={}",
            prelim_linked_count,
            prelinked_ungapped_hits.len(),
            prelim_linked_count - prelinked_ungapped_hits.len(),
            evalue_threshold
        );
    }
    if stage_dump::enabled() {
        stage_dump::dump_ungapped_hits(
            "after_prelim_evalue_reap_before_reevaluate",
            &prelinked_ungapped_hits,
        );
    }
    prelinked_ungapped_hits
}

/// The subject's HSPs re-evaluated with the subject's letters (`s_frames` to score them,
/// `s_frames_report` for the identities) and sorted by score, with the stage dumps of the result.
/// A long subject's statistics count the HSPs the re-evaluation removes
/// (`collect_hsp_saving_stats`, `stats_hsp_filtered_by_reeval`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2609-2737
/// (`Blast_HSPListReevaluateUngapped`, whose sort is at :2733-2734), called at
/// c++/src/algo/blast/core/blast_engine.c:1492-1497.
pub(super) fn reevaluate_subject_hits<S: TbxScore>(
    sc: S,
    prelinked_ungapped_hits: Vec<UngappedHit>,
    (s_frames, s_frames_report): (&[QueryFrame], &[QueryFrame]),
    stage: &SubjectStageArgs<'_>,
    collect_hsp_saving_stats: bool,
    stats_hsp_filtered_by_reeval: &mut usize,
    x_lap: &mut XLap<'_>,
) -> Vec<UngappedHit> {
    let &SubjectStageArgs {
        contexts_ref,
        cutoff_scores,
        timing_enabled,
        reeval_ns,
        reeval_calls,
        ..
    } = stage;
    // NCBI: Blast_HSPListReevaluateUngapped equivalent
    // Reference: blast_engine.c:1492-1497, blast_hits.c:2609-2737
    // Perform batch reevaluation on all HSPs after merging all frames
    let mut ungapped_hits = reevaluate_ungapped_hsp_list(
        sc,
        prelinked_ungapped_hits,
        contexts_ref,
        s_frames,
        s_frames_report,
        cutoff_scores,
        timing_enabled,
        reeval_ns,
        reeval_calls,
        collect_hsp_saving_stats,
        stats_hsp_filtered_by_reeval,
    );
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2733-2734
    // ```c
    // /* Sort the HSP array by score (scores may have changed!) */
    // Blast_HSPListSortByScore(hsp_list);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1374-1381
    // ```c
    // if (!Blast_HSPListIsSortedByScore(hsp_list)) {
    //     qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //           ScoreCompareHSPs);
    // }
    // ```
    // This sort happens before the second BLAST_LinkHsps call in
    // blast_engine.c:1515-1520 and controls tie order for link_hsps.c
    // comparator-equal short HSPs.
    if !ungapped_hits_is_sorted_by_score_ncbi(&ungapped_hits) {
        sort_ungapped_hits_by_score_ncbi(&mut ungapped_hits);
    }
    x_lap.lap(XStage::Reevaluate);
    // Record the post-reevaluation list and the exact input to link_hsps.
    if stage_dump::enabled() {
        stage_dump::dump_ungapped_hits("after_reevaluate", &ungapped_hits);
        stage_dump::dump_ungapped_hits("before_link_hsps", &ungapped_hits);
    }
    ungapped_hits
}

/// The second linking call on the re-evaluated list (`link_subject_hits`, armed on
/// `x_link_timer`), timed for LOSAT_TIMING's `sum_stats_linking` line together with its stage dump
/// and trace, then the list's `hsps_after_chain` count (LOSAT_DIAGNOSTICS).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1510-1520 (the re-link of
/// `BLAST_PreliminarySearchEngine` after `Blast_HSPListReevaluateUngapped`).
pub(super) fn relink_subject_hits(
    link: &SubjectLinkArgs<'_>,
    ungapped_hits: Vec<UngappedHit>,
    stage: &SubjectStageArgs<'_>,
    x_link_timer: &XLinkTimer,
    x_lap: &mut XLap<'_>,
) -> Vec<UngappedHit> {
    let &SubjectStageArgs {
        timing_enabled,
        linking_ns,
        linking_calls,
        diag_enabled,
        diagnostics,
        ..
    } = stage;
    let t_linking = if timing_enabled {
        Some(Instant::now())
    } else {
        None
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c:553-558
    // ```c
    // for (frame_index=0; frame_index<num_query_frames; frame_index++)
    // {
    //    hp_start->next = hp_frame_start[frame_index];
    //    number_of_hsps = hp_frame_number[frame_index];
    // }
    // ```
    x_link_timer.arm(LinkCall::Second);
    let linked = link_subject_hits(link, ungapped_hits);
    x_link_timer.disarm();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:955-974
    // ```c
    //          best[ordering_method]->start_of_chain = TRUE;
    //          best[ordering_method]->hsp->evalue    = prob[ordering_method];
    // ...
    //             H->linked_set = linked_set;
    //             H->ordering_method = ordering_method;
    //             H->hsp->evalue = prob[ordering_method];
    // ```
    // Record linked_set/start_of_chain/order/e-value before output
    // coordinate conversion can obscure frame-relative HSP identity.
    if stage_dump::enabled() {
        stage_dump::dump_ungapped_hits("after_link_hsps_before_output_conversion", &linked);
    }
    if trace_hsp_target().is_some() {
        for h in &linked {
            trace_ungapped_hit_if_match("after_linking", h);
        }
    }
    if let Some(t) = t_linking {
        let elapsed = t.elapsed();
        linking_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
        linking_calls.fetch_add(1, AtomicOrdering::Relaxed);
    }
    x_lap.lap(XStage::Link2);
    if diag_enabled {
        diagnostics
            .base
            .hsps_after_chain
            .fetch_add(linked.len(), AtomicOrdering::Relaxed);
    }
    linked
}
