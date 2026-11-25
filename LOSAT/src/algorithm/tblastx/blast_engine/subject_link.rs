//! The linking parameters of a TBLASTX subject (`subject_linking_params`, made before its frame
//! loop; none without sum statistics) and the end of its search after its frame loop, up to its
//! linked HSP list: the end of its seed stage (`end_subject_seed_stage`), the frame bases of its
//! HSP list (`subject_frame_bases`), the linking call (`first_link_subject_hits`, then
//! `link_subject_hits`; without sum statistics, the HSPs' own E-values), the sort and the reap by
//! the preliminary E-value after the first one (`sort_and_reap_after_first_link`), the
//! re-evaluation with the subject's letters (`reevaluate_subject_hits`) and the second linking
//! call (`relink_subject_hits`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1448-1456 (in
//! `BLAST_PreliminarySearchEngine`: `CalculateLinkHSPCutoffs`), :870-899 (in
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
    apply_sum_stats_even_gap_linking_with_parallel, calculate_link_hsp_cutoffs_ncbi,
    link_one_hit_with_int4_ends, reevaluate_ungapped_hsp_list, sort_ungapped_hits_by_score_int4,
    sort_ungapped_hits_by_score_ncbi, stage_dump, trace_hsp_target, trace_ungapped_hit_if_match,
    ungapped_hits_is_sorted_by_score_ncbi, AtomicOrdering, AtomicU64, DiagStruct,
    DiagnosticCounters, Instant, KarlinParams, LinkCall, LinkingParams, QueryContext, QueryFrame,
    TbxScore, UngappedHit, XStage, BLAST_GAP_DECAY_RATE,
};
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
use crate::algorithm::tblastx::sum_stats_linking::XHelped;
use crate::algorithm::tblastx::x_link_timing::XLinkTimer;
use crate::stats::karlin::evalue_from_raw_score;

/// The subject's linking parameters: its length, the gap decay rate and the linking cutoffs of a
/// subject of `subject_len_nucl` nucleotides, computed once per subject before its scan; both
/// linking calls read them. A search without sum statistics (`-sum_stats false`) has none.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:774-783 (in
/// `BlastHitSavingParametersNew`)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1448-1456
pub(super) fn subject_linking_params(
    sum_stats: bool,
    avg_query_length: i32,
    subject_len_nucl: i64,
    link_db_length_nucl: i64,
    cutoff_score_min: i32,
    linking_params_for_cutoff: &KarlinParams,
) -> Option<LinkingParams> {
    sum_stats.then(|| LinkingParams {
        subject_len_nucl,
        gap_decay_rate: BLAST_GAP_DECAY_RATE,
        cutoffs: calculate_link_hsp_cutoffs_ncbi(
            avg_query_length,
            subject_len_nucl,
            link_db_length_nucl,
            cutoff_score_min,
            linking_params_for_cutoff,
        ),
    })
}

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
        // where `length` is the number of residues (excluding sentinels).
        // Source: c++/src/algo/blast/core/blast_util.c:1098-1101
        base += f.aa_seq.len() as i32 - 1;
    }
    subject_frame_bases
}

/// What both linking calls of a subject read besides its HSP list, made per subject by
/// `process_subject`: the linking cutoffs' Karlin block, the subject's linking parameters (`None`
/// without sum statistics), the contexts with their length adjustments and search spaces, the
/// frame bases, and whether the calls may use the pool.
pub(super) struct SubjectLinkArgs<'a> {
    pub(super) linking_params_for_cutoff: &'a KarlinParams,
    pub(super) linking_params: Option<&'a LinkingParams>,
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

/// The subject's HSP list `hits` linked with sum statistics (even gaps), or, without them
/// (`linking_params` is `None`), with each HSP's own E-value (`evalues_without_linking`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1761-1813 (`BLAST_LinkHsps`),
/// called at c++/src/algo/blast/core/blast_engine.c:871-873 and :1516-1520, whose `else`
/// branches (:874-890 and :1521-1526) take the E-values without linking.
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
    let Some(linking_params) = linking_params else {
        return evalues_without_linking(link, hits);
    };
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

/// The E-values of the subject's HSPs in a search without sum statistics: each HSP's own, from
/// its score with its context's Karlin block and search space. The list keeps its order and the
/// HSPs their linked-set size `num` (0, which every reader takes as 1).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:874-890 (the first call,
/// in `s_BlastSearchEngineCore`) and :1521-1526 (the second, after
/// `Blast_HSPListReevaluateUngapped`)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1887-1893 (in
/// `Blast_HSPListGetEvalues`; `sbp->gbp` is unset for an ungapped search)
/// The gap decay rate passed is 0, so the divisor stays 1, and the scaling factor is 1 (not RPS):
/// both leave the value as `BLAST_KarlinStoE_simple` gives it (`evalue_from_raw_score`). The HSP's
/// `num` stays as `Blast_HSPNew`'s `calloc` left it (blast_hits.c:143, called by `Blast_HSPInit` at
/// :162), as no linking sets it.
#[cold]
#[inline(never)]
fn evalues_without_linking(
    link: &SubjectLinkArgs<'_>,
    mut hits: Vec<UngappedHit>,
) -> Vec<UngappedHit> {
    for hit in &mut hits {
        let context = hit.ctx_idx as usize;
        hit.e_value = evalue_from_raw_score(
            hit.raw_score,
            &link.contexts_ref[context].karlin_params,
            link.eff_searchsp_per_context[context] as f64,
        );
    }
    hits
}

/// The first linking call of the subject's HSP list from its frame loop (`hits`): an empty list
/// stays as it is; a list of the one-hit word finder (`window` 0 or less) that holds an HSP of
/// length 0 or less links with NCBI's `Int4` ends (`link_subject_hits_int4`), which keeps those
/// HSPs; any other, and every list of a search without sum statistics, goes to `link_subject_hits`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-890 (cited at the call in
/// `process_subject`).
pub(super) fn first_link_subject_hits(
    link: &SubjectLinkArgs<'_>,
    hits: Vec<UngappedHit>,
    window: i32,
) -> Vec<UngappedHit> {
    match link.linking_params {
        _ if hits.is_empty() => hits,
        Some(linking_params) if window <= 0 && holds_int4_inverted(&hits) => {
            link_subject_hits_int4(link, linking_params, hits)
        }
        _ => link_subject_hits(link, hits),
    }
}

// No NCBI counterpart: which form of the first linking call a subject of the one-hit word finder
// (`-window_size 0`) takes. NCBI links every list the same way; LOSAT's default kernel takes HSPs
// whose end is after their start (`link_one_hit_with_int4_ends` in `linking.rs`).
/// Whether the subject's list (of the one-hit word finder) holds an HSP of length 0 or less.
#[inline(never)]
fn holds_int4_inverted(hits: &[UngappedHit]) -> bool {
    hits.iter().any(UngappedHit::is_int4_inverted)
}

/// The first linking call of a subject of the one-hit word finder whose list holds an HSP of
/// length 0 or less (`holds_int4_inverted`), with NCBI's `Int4` ends; those HSPs stay in the list,
/// with their linked E-values, up to the re-evaluation's fetch of the subject, which NCBI makes
/// for every list that still holds an HSP after the reap, and leave after it, where NCBI's
/// re-evaluation deletes them (`drop_unread_int4_inverted`, `link_one_hit_with_int4_ends`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873 (the first
/// `BLAST_LinkHsps`) and c++/src/algo/blast/core/blast_hits.c:2703-2731
/// (`Blast_HSPListReevaluateUngapped`, called at blast_engine.c:1492-1497).
#[cold]
#[inline(never)]
fn link_subject_hits_int4(
    link: &SubjectLinkArgs<'_>,
    linking_params: &LinkingParams,
    hits: Vec<UngappedHit>,
) -> Vec<UngappedHit> {
    let &SubjectLinkArgs {
        linking_params_for_cutoff,
        linking_params: _,
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

/// The first linking call's list sorted by score (when it is not; a list that was not linked,
/// without sum statistics, keeps its order; a list of the one-hit word finder, `window` 0, with
/// its ends read as NCBI's `Int4`), then reaped by the preliminary E-value, with its stage dumps
/// and the `[DEBUG PRELIM_REAP]` line of LOSAT_DIAGNOSTICS.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1802-1803 (the sort that ends
/// `BLAST_LinkHsps`) and c++/src/algo/blast/core/blast_engine.c:643-674
/// (`s_Blast_HSPListReapByPrelimEvalue`), called at :899.
pub(super) fn sort_and_reap_after_first_link(
    mut prelinked_ungapped_hits: Vec<UngappedHit>,
    link: &SubjectLinkArgs<'_>,
    stage: &SubjectStageArgs<'_>,
    window: i32,
    x_lap: &mut XLap<'_>,
) -> Vec<UngappedHit> {
    let &SubjectStageArgs {
        evalue_threshold,
        diag_enabled,
        ..
    } = stage;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1802-1803
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
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-874
    // Without sum statistics nothing sorts the list here: `Blast_HSPListGetEvalues` keeps the
    // order of the frame loop's list, and the re-evaluation's sort (blast_hits.c:2733-2734) is
    // the next.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1353
    // A list of the one-hit word finder may still hold HSPs of length 0 or less, whose ends
    // `ScoreCompareHSPs` reads as `Int4` (`score_compare_ungapped_hits_int4`, the comparator of
    // that list's sort in the frame loop too); on ends of 0 or more the two comparators agree.
    if link.linking_params.is_some() {
        if window == 0 {
            sort_ungapped_hits_by_score_int4(&mut prelinked_ungapped_hits);
        } else if !ungapped_hits_is_sorted_by_score_ncbi(&prelinked_ungapped_hits) {
            sort_ungapped_hits_by_score_ncbi(&mut prelinked_ungapped_hits);
        }
    }
    if stage_dump::enabled() {
        stage_dump::dump_ungapped_hits(
            "after_prelim_link_hsps_before_reevaluate",
            &prelinked_ungapped_hits,
        );
    }
    let prelim_linked_count = prelinked_ungapped_hits.len();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:654-662
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

/// The HSPs of length 0 or less of the one-hit word finder (`window` 0 or less) removed from the
/// subject's list (`hits`) before its re-evaluation, which deletes them in NCBI, with sum
/// statistics or without.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:687,699 (in
/// `Blast_HSPReevaluateWithAmbiguitiesUngapped`) and :452-453 (in `s_UpdateReevaluatedHSP`)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1480-1497 and
/// c++/src/algo/blast/core/blast_hits.c:2651-2665 (the fetch of the subject, first in
/// `Blast_HSPListReevaluateUngapped`)
/// NCBI's re-evaluation reads no letter of such an HSP and scores it 0, below every context's
/// cutoff (at least 1). Up to here such an HSP is in NCBI's list, at NCBI's place: the first
/// linking call (`link_subject_hits_int4`) and the sort after it
/// (`sort_and_reap_after_first_link`) read its ends as NCBI's `Int4`, and the reap reads its
/// E-value, the linking call's (without sum statistics its own, `evalues_without_linking`;
/// nothing sorts that list before the re-evaluation, and the `u32` and `Int4` comparators of the
/// frame loop's sorts agree on HSPs of length 1 or more). So the HSPs kept are NCBI's, in the order
/// of NCBI's list. They are removed after the fetch of the subject, which NCBI makes for a list
/// that still holds them (`SubjectEmit::reevaluation_fetch_keeps`, the frames of its letters),
/// and before the refusal of a subject of fewer than 12 bases, whose frame NCBI reads only for an
/// HSP of length 1 or more.
#[inline]
pub(super) fn drop_unread_int4_inverted(window: i32, hits: &mut Vec<UngappedHit>) {
    if window <= 0 {
        hits.retain(|hit| !hit.is_int4_inverted());
    }
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1492-1497, c++/src/algo/blast/core/blast_hits.c:2609-2737
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2733-2734
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1381
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
    x_link_timer.arm(LinkCall::Second);
    let linked = link_subject_hits(link, ungapped_hits);
    x_link_timer.disarm();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:955-974
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

#[cfg(test)]
#[path = "subject_link_tests.rs"]
mod tests;
