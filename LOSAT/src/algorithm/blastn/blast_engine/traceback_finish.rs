//! The last two steps of the blastn traceback of one subject (or of one query's HSPs of it,
//! LOSAT_X_QGROUPPAR): the purges and re-evaluation after the traceback loop
//! (`traceback_post_purge`), then the score sort, the final interval tree and the e-values
//! (`traceback_final_list`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:635-692 (the end of
//! `Blast_TracebackFromHSPList`) and c++/src/algo/blast/core/blast_traceback.c:199-253
//! (`s_HSPListPostTracebackUpdate`, which `BLAST_ComputeTraceback_MT` runs on the list next).
//!
//! This module is a child of `run` for the types it keeps private (`BlastnTiming`, the e-value
//! helpers). The two functions are the bodies of what were two closures of
//! `SubjectPasses::search_subject`; the values those closures captured are among the
//! fields of `TracebackCtx`, which `traceback_hsp`, `traceback_list` and `traceback_subject` share,
//! and which pool threads read (the query group tasks of `x_qgroup_par`, the look-ahead helpers,
//! and the speculative batch's jobs).

use super::{
    blast_hsp_test_identity_and_length, blastn_trace, calculate_blastn_context_statistics,
    context_query, format_gap_edit_ops_for_trace, hsp_survives_evalue_reap,
    purge_hsps_with_common_endpoints, purge_hsps_with_common_endpoints_ex,
    reevaluate_hsp_with_ambiguities_gapped_ex, sort_hsps_by_score, BlastIntervalTree, BlastnHsp,
    BlastnTiming, ContextKarlin, IndexMethod, QueryContext, ReevalParams, TreeHsp,
};
use crate::algorithm::blastn::alignment::gapped::BlastnaMatrix;

/// The values of one subject's search that its traceback reads: these two steps, the traceback
/// loop (`traceback_list`), the traceback of one HSP (`traceback_hsp`) and the subject's side
/// (`traceback_subject`). Built once per subject by `traceback_subject::traceback_subject`, once
/// the subject's blastna letters are made; every field is a shared reference or a `Copy` value
/// (mostly the types the closures captured), so the struct is `Sync` for the query group tasks,
/// the look-ahead helpers and the speculative batch's jobs on the pool.
pub(super) struct TracebackCtx<'a> {
    /// `LOSAT_TIMING` counters, `None` when timing is off.
    pub(super) timing_ref: Option<&'a BlastnTiming>,
    pub(super) timing_enabled: bool,
    pub(super) blastn_trace_enabled: bool,
    /// The subject's index in the batch and its trace id.
    pub(super) s_idx: usize,
    pub(super) s_id: &'a str,
    /// The subject's letters in blastna.
    pub(super) s_seq_blastna: &'a [u8],
    /// EXPERIMENT (LOSAT_X_BNTBRANGE): the subject's fence flag (`alignment::x_fence`), `None`
    /// unless its letters are decoded in ranges.
    pub(super) x_fenced: Option<&'a std::sync::atomic::AtomicBool>,
    /// The queries' contexts: their letters, their lengths and frames (`QueryContext`, which the
    /// traceback loop reads), the cutoffs, the effective search spaces and the Karlin blocks, by
    /// context.
    pub(super) encoded_queries_blastna: &'a Vec<&'a [u8]>,
    pub(super) encoded_queries_nomask: &'a [&'a [u8]],
    /// The same letters concatenated, with masks and without
    /// (`TraceBatch::encoded_query_concat_blastna`, `TraceBatch::encoded_query_concat_nomask`).
    pub(super) encoded_query_concat_blastna: &'a [u8],
    pub(super) encoded_query_concat_nomask: &'a [u8],
    pub(super) query_contexts: &'a [QueryContext],
    /// The queries' lengths, by query.
    pub(super) query_lengths: &'a [usize],
    pub(super) cutoff_scores: &'a Vec<i32>,
    pub(super) hit_saving_cutoff_scores: &'a Vec<i32>,
    pub(super) query_eff_searchsp: &'a Vec<i64>,
    pub(super) search_karlin_ref: &'a Vec<ContextKarlin>,
    pub(super) round_down_evalue_score: bool,
    /// The scoring of the alignment.
    pub(super) reward: i32,
    pub(super) penalty: i32,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    /// The X-dropoff of the traceback alignment.
    pub(super) x_drop_final: i32,
    pub(super) score_matrix: &'a BlastnaMatrix,
    pub(super) use_dp: bool,
    /// The hit-saving options.
    pub(super) percent_identity: f64,
    pub(super) min_diag_separation: i32,
    pub(super) evalue_threshold: f64,
}

/// The passes after the traceback loop up to the second purge, on one list: the subject's, or
/// one query's (LOSAT_X_QGROUPPAR). Returns the list and its lengths after each step (for the
/// debug line).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:635-669
/// (`Blast_TracebackFromHSPList`, from the first purge to the second), which calls
/// `Blast_HSPListPurgeHSPsWithCommonEndpoints` (c++/src/algo/blast/core/blast_hits.c:2455-2535),
/// `Blast_HSPReevaluateWithAmbiguitiesGapped` (blast_hits.c:479-647),
/// `Blast_HSPTestIdentityAndLength` (blast_hits.c:1004-1025) and
/// `Blast_HSPListPurgeNullHSPs` (blast_hits.c:2225-2256).
pub(super) fn traceback_post_purge(
    ctx: &TracebackCtx<'_>,
    local_hits: Vec<BlastnHsp>,
) -> (Vec<BlastnHsp>, usize, usize, usize, usize) {
    let &TracebackCtx {
        timing_ref,
        timing_enabled,
        blastn_trace_enabled,
        s_idx,
        s_id,
        s_seq_blastna,
        encoded_query_concat_blastna,
        encoded_query_concat_nomask,
        query_contexts,
        cutoff_scores,
        hit_saving_cutoff_scores,
        query_eff_searchsp,
        search_karlin_ref,
        round_down_evalue_score,
        reward,
        penalty,
        gap_open,
        gap_extend,
        score_matrix,
        use_dp,
        percent_identity,
        ..
    } = ctx;

    // Step 2: Endpoint purging pass 1 (trim, purge=false)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:637-638
    let hits_step1 = local_hits.len();
    let purge_pass1_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    let (mut local_hits, mut extra_start) = purge_hsps_with_common_endpoints_ex(local_hits, false);
    if let (Some(timing), Some(purge_pass1_start)) = (timing_ref, purge_pass1_start) {
        BlastnTiming::record_duration(&timing.purge_common_endpoint_pass1_ns, purge_pass1_start);
        BlastnTiming::record_count(
            &timing.traceback_removed_endpoint_pass1_hsps,
            hits_step1.saturating_sub(local_hits.len()) as u64,
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:633-647
    if blastn_trace_enabled && blastn_trace::should_trace_subject("purge", None, s_idx, s_id) {
        blastn_trace::log(
            "purge",
            format!(
                "subject={}({}) pass=common_endpoint_trim before={} after={} extra_start={}",
                s_id,
                s_idx,
                hits_step1,
                local_hits.len(),
                extra_start
            ),
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:640-644
    if !use_dp {
        extra_start = 0;
    }
    let hits_step2 = local_hits.len();

    // Step 3: Re-evaluate trimmed HSPs
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:647-665
    // A trimmed HSP can still be extendable, so it is scored again

    for hit in local_hits.iter_mut().skip(extra_start) {
        // Get sequences for re-evaluation
        let context_idx = (hit.q_idx as usize) * 2 + if hit.query_frame < 0 { 1 } else { 0 };
        let cutoff = hit_saving_cutoff_scores
            .get(context_idx)
            .copied()
            .unwrap_or_else(|| cutoff_scores.get(context_idx).copied().unwrap_or(0));

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1887-1891
        let eff_searchsp = query_eff_searchsp.get(context_idx).copied().unwrap_or(0);
        let reeval_params = ReevalParams {
            lambda: search_karlin_ref[context_idx].gapped.lambda,
            k: search_karlin_ref[context_idx].gapped.k,
            eff_searchsp,
            round_down_evalue_score,
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:653-665 (reevaluate with blastna sequences)
        // NCBI's `query` and `query_nomask` (c++/src/algo/blast/core/blast_traceback.c:386-389): the context in the
        // concatenated query, with masks and without, read at the HSP's offsets, which lie
        // outside the context for an HSP outside it.
        let ctx = &query_contexts[context_idx];
        let query = context_query(encoded_query_concat_blastna, ctx);
        let query_nomask = context_query(encoded_query_concat_nomask, ctx);
        let s_seq_eval = s_seq_blastna;
        let reeval_start = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let delete = reevaluate_hsp_with_ambiguities_gapped_ex(
            hit,
            query,
            s_seq_eval,
            reward,
            penalty,
            gap_open,
            gap_extend,
            cutoff,
            score_matrix,
            Some(&reeval_params),
        );
        if let (Some(timing), Some(reeval_start)) = (timing_ref, reeval_start) {
            BlastnTiming::record_duration(&timing.reevaluate_ambiguities_ns, reeval_start);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:647-665
        let trace_reeval = blastn_trace_enabled
            && blastn_trace::should_trace_range(
                "purge",
                context_idx as u32,
                s_idx,
                s_id,
                hit.internal_q_offset_0,
                hit.internal_q_end_0,
                hit.internal_s_offset_0,
                hit.internal_s_end_0,
                hit.query_length,
                hit.query_frame,
            );
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:656-663
        if delete {
            if let Some(timing) = timing_ref {
                BlastnTiming::record_count(&timing.traceback_deleted_reevaluation_hsps, 1);
            }
            if trace_reeval {
                blastn_trace::log(
                    "purge",
                    format!(
                        "subject={}({}) context={} pass=reeval q{}..{} s{}..{} deleted=true reason=score_or_gap_info cutoff={}",
                        s_id,
                        s_idx,
                        context_idx,
                        hit.internal_q_offset_0,
                        hit.internal_q_end_0,
                        hit.internal_s_offset_0,
                        hit.internal_s_end_0,
                        cutoff
                    ),
                );
            }
            hit.raw_score = i32::MIN; // Mark for removal
            continue;
        }
        let identity_length_start = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let delete_identity_length =
            blast_hsp_test_identity_and_length(hit, query_nomask, s_seq_eval, percent_identity);
        if let (Some(timing), Some(identity_length_start)) = (timing_ref, identity_length_start) {
            BlastnTiming::record_duration(&timing.identity_length_test_ns, identity_length_start);
        }
        if delete_identity_length {
            if let Some(timing) = timing_ref {
                BlastnTiming::record_count(&timing.traceback_deleted_identity_length_hsps, 1);
            }
            if trace_reeval {
                blastn_trace::log(
                    "purge",
                    format!(
                        "subject={}({}) context={} pass=identity_length q{}..{} s{}..{} deleted=true raw_score={} identities={} aln_len={} identity={:.6}",
                        s_id,
                        s_idx,
                        context_idx,
                        hit.internal_q_offset_0,
                        hit.internal_q_end_0,
                        hit.internal_s_offset_0,
                        hit.internal_s_end_0,
                        hit.raw_score,
                        hit.num_ident,
                        hit.length,
                        hit.identity
                    ),
                );
            }
            hit.raw_score = i32::MIN; // Mark for removal
        } else if trace_reeval {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:647-665
            blastn_trace::log(
                "purge",
                format!(
                    "subject={}({}) context={} pass=reeval q{}..{} s{}..{} deleted=false raw_score={} evalue={:.12e} bit_score={:.12} identities={} aln_len={} identity={:.6} gap_info={}",
                    s_id,
                    s_idx,
                    context_idx,
                    hit.internal_q_offset_0,
                    hit.internal_q_end_0,
                    hit.internal_s_offset_0,
                    hit.internal_s_end_0,
                    hit.raw_score,
                    hit.e_value,
                    hit.bit_score,
                    hit.num_ident,
                    hit.length,
                    hit.identity,
                    hit.gap_info
                        .as_deref()
                        .map(format_gap_edit_ops_for_trace)
                        .unwrap_or_else(|| "[]".to_string())
                ),
            );
        }
    }
    local_hits.retain(|h| h.raw_score != i32::MIN);
    let hits_step3 = local_hits.len();

    // Step 4: Endpoint purging pass 2 (delete, purge=true) - BLASTN only
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:667-668
    let purge_pass2_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    let local_hits = purge_hsps_with_common_endpoints(local_hits);
    let hits_step4 = local_hits.len();
    if let (Some(timing), Some(purge_pass2_start)) = (timing_ref, purge_pass2_start) {
        BlastnTiming::record_duration(&timing.purge_common_endpoint_pass2_ns, purge_pass2_start);
        BlastnTiming::record_count(
            &timing.traceback_removed_endpoint_pass2_hsps,
            hits_step3.saturating_sub(hits_step4) as u64,
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:666-668
    if blastn_trace_enabled && blastn_trace::should_trace_subject("purge", None, s_idx, s_id) {
        blastn_trace::log(
            "purge",
            format!(
                "subject={}({}) pass=post_reeval before={} after_reeval={} after_delete_purge={}",
                s_id, s_idx, hits_step2, hits_step3, hits_step4
            ),
        );
    }

    (local_hits, hits_step1, hits_step2, hits_step3, hits_step4)
}

/// The passes after the second purge on one list: the score sort, the final tree and the
/// e-values. Returns the HSPs kept and their number after the tree (for the debug line).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:671-692
/// (`Blast_TracebackFromHSPList`, from the score sort to the end of the final tree pass), which
/// calls `Blast_HSPListSortByScore` (c++/src/algo/blast/core/blast_hits.c:1374-1383) and
/// `Blast_IntervalTreeReset` (c++/src/algo/blast/core/blast_itree.c:193-202); then
/// c++/src/algo/blast/core/blast_traceback.c:199-253 (`s_HSPListPostTracebackUpdate`) for the
/// e-values and the reap: `Blast_HSPListGetEvalues` (blast_hits.c:1811-1905),
/// `Blast_HSPListReapByEvalue` (blast_hits.c:1976-2008).
pub(super) fn traceback_final_list(
    ctx: &TracebackCtx<'_>,
    mut local_hits: Vec<BlastnHsp>,
    interval_tree: &mut BlastIntervalTree,
) -> (Vec<BlastnHsp>, usize) {
    let &TracebackCtx {
        timing_ref,
        timing_enabled,
        blastn_trace_enabled,
        s_idx,
        s_id,
        query_eff_searchsp,
        search_karlin_ref,
        round_down_evalue_score,
        min_diag_separation,
        evalue_threshold,
        ..
    } = ctx;
    // Step 5: Re-sort by gapped score
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:671-672
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353 (ScoreCompareHSPs)
    let score_resort_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:671-672
    sort_hsps_by_score(&mut local_hits);
    if let (Some(timing), Some(score_resort_start)) = (timing_ref, score_resort_start) {
        BlastnTiming::record_duration(&timing.traceback_score_resort_ns, score_resort_start);
    }

    // Step 6: Phase 2 - Tree reset and second containment pass
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:678
    interval_tree.reset();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:66-73
    // The final containment pass adds at most one tree HSP for each surviving
    // HSP. Reserving here only removes allocator churn; containment and order
    // are still controlled by BlastIntervalTreeContainsHSP/AddHSP.
    interval_tree.reserve_nodes_for_hsps(local_hits.len());

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:679-692
    // Purge HSPs that lie inside an earlier one. The list is in score order, so a contained
    // HSP can never outscore its container.
    // No NCBI counterpart (memory placement only): the kept HSPs move to the front of the
    // same list, in order (`kept`), instead of into a second list of its length.
    let mut kept = 0usize;
    for read in 0..local_hits.len() {
        let hit = &local_hits[read];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:537-550
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3908-3914
        // NCBI uses canonical coordinates: subject.offset < subject.end always
        let query_context_offset = hit.internal_query_context_offset;
        let tree_hsp = TreeHsp {
            query_offset: hit.internal_q_offset_0,
            query_end: hit.internal_q_end_0,
            subject_offset: hit.internal_s_offset_0 as i32,
            subject_end: hit.internal_s_end_0 as i32,
            score: hit.raw_score,
            query_frame: hit.query_frame,
            query_length: hit.query_length as i32,
            query_context_offset,
            subject_frame_sign: 1,
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:678-711
        let final_tree_contains_start = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let final_tree_containing_hsp =
            interval_tree.containing_hsp(&tree_hsp, query_context_offset, min_diag_separation);
        let final_tree_contains = final_tree_containing_hsp.is_some();
        if let (Some(timing), Some(final_tree_contains_start)) =
            (timing_ref, final_tree_contains_start)
        {
            BlastnTiming::record_duration(
                &timing.traceback_final_tree_contains_ns,
                final_tree_contains_start,
            );
        }
        if blastn_trace_enabled
            && blastn_trace::should_trace_range(
                "purge",
                hit.q_idx * 2 + if hit.query_frame < 0 { 1 } else { 0 },
                s_idx,
                s_id,
                hit.internal_q_offset_0,
                hit.internal_q_end_0,
                hit.internal_s_offset_0,
                hit.internal_s_end_0,
                hit.query_length,
                hit.query_frame,
            )
        {
            blastn_trace::log(
                "purge",
                format!(
                    "subject={}({}) context={} pass=final_interval_tree q{}..{} s{}..{} raw_score={} contained={} containing_hsp={:?} min_diag_separation={}",
                    s_id,
                    s_idx,
                    hit.q_idx * 2 + if hit.query_frame < 0 { 1 } else { 0 },
                    hit.internal_q_offset_0,
                    hit.internal_q_end_0,
                    hit.internal_s_offset_0,
                    hit.internal_s_end_0,
                    hit.raw_score,
                    final_tree_contains,
                    final_tree_containing_hsp,
                    min_diag_separation
                ),
            );
        }
        if !final_tree_contains {
            let final_tree_add_start = if timing_enabled {
                Some(std::time::Instant::now())
            } else {
                None
            };
            interval_tree.add_hsp(tree_hsp, query_context_offset, IndexMethod::QueryAndSubject);
            if let (Some(timing), Some(final_tree_add_start)) = (timing_ref, final_tree_add_start) {
                BlastnTiming::record_duration(
                    &timing.traceback_final_tree_add_ns,
                    final_tree_add_start,
                );
            }
            local_hits.swap(kept, read);
            kept += 1;
        } else if let Some(timing) = timing_ref {
            BlastnTiming::record_count(&timing.traceback_removed_final_tree_hsps, 1);
        }
        // else: HSP is contained within another, skip (implicit delete)
    }
    local_hits.truncate(kept);
    let mut final_hits = local_hits;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:678-711
    // Keep the final-containment count separate from the later E-value reap.
    let final_tree_hit_count = final_hits.len();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:234-250
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1983-2003
    // Apply final gapped E-value assignment and reaping after all
    // post-traceback survivor-set operations, matching
    // s_HSPListPostTracebackUpdate.
    let before_evalue_reap = final_hits.len();
    for hit in final_hits.iter_mut() {
        let context_idx = (hit.q_idx as usize) * 2 + if hit.query_frame < 0 { 1 } else { 0 };
        let eff_searchsp = query_eff_searchsp.get(context_idx).copied().unwrap_or(0);
        let (bit_score, eval) = calculate_blastn_context_statistics(
            hit.raw_score,
            &search_karlin_ref[context_idx].gapped,
            eff_searchsp,
            round_down_evalue_score,
        );
        hit.bit_score = bit_score;
        hit.e_value = eval;
    }
    final_hits.retain(|hit| hsp_survives_evalue_reap(hit.e_value, evalue_threshold));
    if let Some(timing) = timing_ref {
        BlastnTiming::record_count(
            &timing.traceback_deleted_evalue_cutoff_hsps,
            before_evalue_reap.saturating_sub(final_hits.len()) as u64,
        );
    }
    (final_hits, final_tree_hit_count)
}
