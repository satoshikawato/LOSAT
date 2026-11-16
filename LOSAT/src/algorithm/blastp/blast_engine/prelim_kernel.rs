//! The gapped half of a BLASTP subject's preliminary search (`run_call`'s kernel `search_subject`,
//! after its scan and ungapped extensions): the subject's ungapped hits sorted by score and, with
//! `-chaining` and BLOSUM62, chained; the queries whose gapped extensions start restricted; the
//! gapped extension of each hit outside the interval tree of the HSPs kept so far, with the exact
//! retry of a query whose restricted score is inconclusive; the purge of the HSPs with common
//! endpoints; and the subject's HSP lists per query with their E-values, reaped by the preliminary
//! E-value (`BlastpPrelimKernel::gapped_subject`). The method is the kernel's statements at their
//! place, in their order: the kernel calls it once per subject with letters, after the scan.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:522-526
//! ```c
//!         if (aux_struct->GetGappedScore) {
//!             status = aux_struct->GetGappedScore(program_number, query,
//!                     query_info,
//!                     subject, gap_align, score_params, ext_params, hit_params,
//!                     word_params, init_hitlist, &hsp_list, gapped_stats, NULL);
//! ```
//!
//! This module is a child of `blast_engine` for the private items of that file it reads (the
//! subject's scratch, reference and result, the preliminary values, the initial hits, the clocks and
//! the helpers of the gapped stage) and its sibling modules `x_hits_flat` and `x_ih_pool`.

use std::sync::Arc;

use super::{
    blastp_adjust_hsp_offsets_and_get_query_data, blastp_advance_preliminary_hit_index,
    blastp_get_ungapped_hsp_context_with_bounds, blastp_restricted_align_enabled,
    blastp_timing_start, blastp_trace_init_hsp_overlaps_target, build_restricted_align_array_into,
    chain_blastp_init_hsps, extend_preliminary_blastp_hit, ncbi_qsort_init_hsps_by_score,
    preliminary_tree_hsp, prepare_preliminary_hits_for_kappa_in_place,
    rebuild_preliminary_interval_tree_from_list, remove_preliminary_hits_for_query,
    should_retry_exact_preliminary_gapped_alignment, split_preliminary_hits_by_query_into,
    x_build_compact_preliminary_hsp_list, x_hits_flat, x_ih_pool, BlastpContextLookupBounds,
    BlastpPrelimValues, BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch, BlastpTiming,
    InitHSP, BLASTP_RESTRICTED_MULT,
};
use crate::algorithm::blastn::interval_tree::{IndexMethod, TreeHsp};
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::blastp::gapalign::BlastpGappedAlignmentMode;
use crate::algorithm::blastp::hsp::reap_hsplist_by_evalue;
use crate::algorithm::blastp::kappa::{blastp_trace_hsp_target, blastp_trace_matches_pair};
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::stats::{BlastGumbelBlk, KarlinParams};

/// What the gapped half of the kernel reads of `run_call` (`gapped_subject`): each field is
/// `run_call`'s value of the same name, which the kernel closure captured before (references of the
/// captures' types, `Copy` values by value), so that the moved statements read them as they did in
/// the closure. `run_call` makes it once, before the kernel.
pub(super) struct BlastpPrelimKernel<'a> {
    // The call's options and its queries.
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) query_ids: &'a Vec<Arc<[u8]>>,
    pub(super) contexts: &'a Vec<QueryContext>,
    pub(super) context_lookup_bounds: BlastpContextLookupBounds,
    pub(super) concatenated_query: &'a Vec<u8>,
    // The gapped extension's X-drop and the statistics of the preliminary E-values.
    pub(super) x_drop_gapped: i32,
    pub(super) gapped_params: &'a KarlinParams,
    pub(super) gapped_gumbel: &'a Option<BlastGumbelBlk>,
    // The trace's first query (LOSAT_DEBUG_BLASTP_TRACE_HSP) and the clocks (LOSAT_TIMING).
    pub(super) x_trace_query_base: u32,
    pub(super) timing_enabled: bool,
    pub(super) timing: &'a Option<Arc<BlastpTiming>>,
}

impl BlastpPrelimKernel<'_> {
    /// The gapped half of the kernel for one subject with letters (`search_subject`, after the
    /// scan): `init_hsps` are the subject's ungapped hits (the scratch's list, which goes back to the
    /// scratch here), and the subject's HSP lists go into `subject_result`. The arrays that the method
    /// takes of the scratch (the restricted alignment's flags, the HSP list and its counts per query,
    /// the lists per query) go back to it before the method returns; no statement leaves it early.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3664-3666
    /// ```c
    /// Int2 BLAST_GetGappedScore (EBlastProgramType program_number,
    ///         BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
    ///         BLAST_SequenceBlk* subject,
    /// ```
    pub(super) fn gapped_subject(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject: BlastpSubjectRef<'_>,
        subject_result: &mut BlastpSubjectResult,
        values: &BlastpPrelimValues<'_>,
        mut init_hsps: Vec<InitHSP>,
    ) {
        // `run_call`'s values as the closure read them: a reference or a `Copy` value each.
        let Self {
            args,
            query_ids,
            contexts,
            context_lookup_bounds,
            concatenated_query,
            x_drop_gapped,
            gapped_params,
            gapped_gumbel,
            x_trace_query_base,
            timing_enabled,
            timing,
        } = *self;
        // The kernel's locals that this half reads: the cutoffs, the subject's letters after the
        // sentinel byte, and the interval tree that the scan reset for the subject.
        let (word_cutoff_scores, hit_cutoff_scores) =
            (values.word_cutoff_scores, values.hit_cutoff_scores);
        let subject_sequence = &subject.aa_seq[1..];
        let interval_tree = &mut scratch.preliminary_interval_tree;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:228-234
        // ```c
        // Blast_InitHitListSortByScore(init_hitlist);
        // return status;
        // ```
        ncbi_qsort_init_hsps_by_score(&mut init_hsps);
        let init_hsps = if args.chaining && args.matrix_name == "BLOSUM62" {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3716-3724
            // ```c
            // if (ext_params->options->chaining &&
            //     is_prot && strcmp(score_params->options->matrix, "BLOSUM62") == 0) {
            //     status = s_ChainingAlignment(..., init_hitlist);
            // }
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:569-576
            // ```c
            // CBlastOptionsLocal::SetMatrixName(const char* matrix)
            // {
            //     if (!matrix)
            //         return;
            //
            //     sfree(m_ScoringOpts->matrix);
            //     m_ScoringOpts->matrix = strdup(matrix);
            // }
            // ```
            // The name is the typed one (`-matrix blosum62` does not chain; the C setter
            // `BlastScoringOptionsSetMatrix` upper-cases it, the C++ one does not).
            chain_blastp_init_hsps(
                init_hsps,
                contexts,
                word_cutoff_scores,
                hit_cutoff_scores,
                args.scoring.gap_open,
                args.scoring.gap_extend,
                &mut scratch.chaining_keep,
                &mut scratch.chaining_nodes,
            )
        } else {
            init_hsps
        };
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3735-3789
        // ```c
        // restricted_align_array = calloc(...);
        // for (i = 0;i < init_hitlist->total;i++) {
        //     init_hsp = &init_hitlist->init_hsp_array[i];
        //     ...
        // }
        // ```
        let mut restricted_align_array = std::mem::take(&mut scratch.restricted_align_array);
        if blastp_restricted_align_enabled(args.app, args.evalue) {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3735-3789
            // ```c
            // if (hit_params->restricted_align && !score_params->options->is_ooframe) {
            //     restricted_align_array = calloc(...);
            //     for (i = 0;i < init_hitlist->total;i++) {
            //         init_hsp = &init_hitlist->init_hsp_array[i];
            //         ...
            //     }
            // }
            // ```
            build_restricted_align_array_into(
                &init_hsps,
                contexts,
                context_lookup_bounds,
                hit_cutoff_scores,
                query_ids.len(),
                &mut scratch.restricted_align_found,
                &mut restricted_align_array,
            );
        } else {
            restricted_align_array.clear();
            restricted_align_array.resize(query_ids.len(), false);
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3971-4008
        // ```c
        // BlastHSPList* new_hsp_list = Blast_HSPListNew(kHspNumMax);
        // ...
        // Blast_HSPListFree(hsp_list);
        // hsp_list = new_hsp_list;
        // ```
        //
        // Keep the accepted preliminary HSPs in one ordered list for the
        // current subject so exact-retry rebuilds preserve NCBI's HSPList
        // order instead of HashMap iteration order.
        let mut preliminary_hsp_list = std::mem::take(&mut scratch.preliminary_hsp_list);
        preliminary_hsp_list.clear();
        let mut preliminary_hsp_counts_by_query =
            std::mem::take(&mut scratch.preliminary_hsp_counts_by_query);
        preliminary_hsp_counts_by_query.clear();
        preliminary_hsp_counts_by_query.resize(query_ids.len(), 0);
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1481-1497
        // ```c
        // if (!gapped_calculation) {
        //     /* The following must be performed for any ungapped
        //        search with a nucleotide database. */
        //     status = Blast_HSPListReevaluateUngapped(...);
        // }
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3697-3698
        // ```c
        // Int4 redo_index = -1; /* these are used for recomputing alignmnets if the */
        // Int4 redo_query = -1; /* approximate alignment score is inconclusive */
        // ```
        let mut redo_index: Option<usize> = None;
        let mut redo_query: Option<usize> = None;
        let mut first_init_hsp_index_by_query: Vec<Option<usize>> = Vec::new();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3870-3878
        // ```c
        // if (index < redo_index && query_index != redo_query) {
        //     continue;
        // }
        // ```
        let mut next_same_query_index: Vec<Option<usize>> = Vec::new();
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80
        // ```c
        // BlastGapDP* dp_mem; /**< scratch structures for dynamic programming */
        // Int4 dp_mem_alloc;  /**< current number of structures allocated */
        // ```
        //
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319
        // ```c
        // s_BlastProtGappedAlignment(..., BlastGapAlignStruct* gap_align, ...)
        // ```
        let preliminary_gapped_start = blastp_timing_start(timing_enabled);
        let mut hit_index = 0usize;
        while hit_index < init_hsps.len() {
            let mut init_hsp = init_hsps[hit_index];
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3852-3868
            // ```c
            // tmp_init_hsp = init_hsp_array[index];
            // ...
            // s_AdjustHspOffsetsAndGetQueryData(query, query_info, init_hsp,
            //                                   &query_tmp, &context);
            // query_index = Blast_GetQueryIndexFromContext(context, program_number);
            // ```
            let (context_idx, ctx, query_sequence) = blastp_adjust_hsp_offsets_and_get_query_data(
                concatenated_query,
                contexts,
                context_lookup_bounds,
                &mut init_hsp,
            );
            let query_index = usize::try_from(ctx.q_idx)
                .expect("NCBI BLAST requires blastp query indices to fit in usize");
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3870-3878
            // ```c
            // /* If redo_index > -1 and redo_query > -1, the main loop is recomputing
            //    gappaed alignments for a single query, becuase the approximate
            //    alignment score was inconclusive. This recomputing was triggered when
            //    index was equal to redo_index. Until index reaches redo_index again,
            //    skip all concatenated queries with index different from redo_query */
            // if (index < redo_index && query_index != redo_query) {
            //     continue;
            // }
            // ```
            if let (Some(redo_index), Some(redo_query)) = (redo_index, redo_query) {
                if hit_index < redo_index && query_index != redo_query {
                    hit_index += 1;
                    continue;
                }
            }
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3918
            // ```c
            // tmp_hsp.query.offset = q_start;
            // tmp_hsp.query.end = q_end;
            // tmp_hsp.subject.offset = s_start;
            // tmp_hsp.subject.end = s_end;
            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                   hit_options->min_diag_separation))
            // ```
            let ungapped_tree_hsp = TreeHsp {
                query_offset: init_hsp.ungapped_data.q_start,
                query_end: init_hsp.query_end_absolute(),
                subject_offset: init_hsp.ungapped_data.s_start,
                subject_end: init_hsp.subject_end(),
                score: init_hsp.ungapped_data.score,
                query_frame: 1,
                query_length: ctx.aa_len as i32,
                query_context_offset: ctx.frame_base,
                subject_frame_sign: 0,
            };
            let tree_contains = interval_tree.contains_hsp(&ungapped_tree_hsp, ctx.frame_base, 0);
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3919
            // ```c
            // tmp_hsp.query.offset = q_start;
            // tmp_hsp.query.end = q_end;
            // tmp_hsp.subject.offset = s_start;
            // tmp_hsp.subject.end = s_end;
            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                   hit_options->min_diag_separation))
            // ```
            if let Some(target) = blastp_trace_hsp_target() {
                if blastp_trace_matches_pair(target, ctx.q_idx + x_trace_query_base, s_idx as u32)
                    && blastp_trace_init_hsp_overlaps_target(target, &init_hsp)
                {
                    eprintln!(
                        "[BLASTP_TRACE_HSP] init hit_index={} q_idx={} s_idx={} ungapped=q{}..{} s{}..{} seed=({}, {}) score={} cutoff={} mode={:?} tree_contains={}",
                        hit_index,
                        ctx.q_idx + x_trace_query_base,
                        s_idx,
                        init_hsp.ungapped_data.q_start + 1,
                        init_hsp.query_end_absolute(),
                        init_hsp.ungapped_data.s_start + 1,
                        init_hsp.subject_end(),
                        init_hsp.q_seed_absolute,
                        init_hsp.s_seed,
                        init_hsp.ungapped_data.score,
                        hit_cutoff_scores[context_idx],
                        if restricted_align_array
                            .get(query_index)
                            .copied()
                            .unwrap_or(false)
                        {
                            BlastpGappedAlignmentMode::Restricted
                        } else {
                            BlastpGappedAlignmentMode::Exact
                        },
                        tree_contains
                    );
                }
            }
            if tree_contains {
                hit_index = blastp_advance_preliminary_hit_index(
                    hit_index,
                    redo_index,
                    &next_same_query_index,
                );
                continue;
            }
            let restricted_alignment = restricted_align_array
                .get(query_index)
                .copied()
                .unwrap_or(false);
            let cutoff = hit_cutoff_scores[context_idx];
            let restricted_cutoff = (BLASTP_RESTRICTED_MULT * cutoff as f64) as i32;
            let mode = if restricted_alignment {
                BlastpGappedAlignmentMode::Restricted
            } else {
                BlastpGappedAlignmentMode::Exact
            };
            let preliminary_hsp = extend_preliminary_blastp_hit(
                &mut init_hsp,
                s_idx as u32,
                context_idx,
                ctx,
                query_sequence,
                subject_sequence,
                args.scoring.matrix,
                args.scoring.gap_open,
                args.scoring.gap_extend,
                x_drop_gapped,
                mode,
                &mut scratch.preliminary_gap_scratch,
            );
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3935-3953
            // ```c
            // max_offset = BlastGetStartForGappedAlignment(...);
            // init_hsp->offsets.qs_offsets.s_off +=
            //     max_offset - init_hsp->offsets.qs_offsets.q_off;
            // init_hsp->offsets.qs_offsets.q_off = max_offset;
            // status = s_BlastProtGappedAlignment(..., init_hsp, ...);
            // ```
            if let Some(target) = blastp_trace_hsp_target() {
                if blastp_trace_matches_pair(target, ctx.q_idx + x_trace_query_base, s_idx as u32)
                    && usize::try_from(preliminary_hsp.query_start + 1).ok() == Some(target.q_start)
                    && usize::try_from(preliminary_hsp.query_end).ok() == Some(target.q_end)
                    && usize::try_from(preliminary_hsp.subject_start + 1).ok()
                        == Some(target.s_start)
                    && usize::try_from(preliminary_hsp.subject_end).ok() == Some(target.s_end)
                {
                    eprintln!(
                        "[BLASTP_TRACE_HSP] preliminary_gapped hit_index={} q={}..{} s={}..{} gapped_start=({}, {}) raw_score={} cutoff={} mode={:?}",
                        hit_index,
                        preliminary_hsp.query_start + 1,
                        preliminary_hsp.query_end,
                        preliminary_hsp.subject_start + 1,
                        preliminary_hsp.subject_end,
                        preliminary_hsp.gapped_query_start,
                        preliminary_hsp.gapped_subject_start,
                        preliminary_hsp.raw_score,
                        cutoff,
                        mode
                    );
                }
            }
            if should_retry_exact_preliminary_gapped_alignment(
                restricted_alignment,
                preliminary_hsp.raw_score,
                cutoff,
                restricted_cutoff,
            ) {
                if let Some(restricted_slot) = restricted_align_array.get_mut(query_index) {
                    *restricted_slot = false;
                }
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3975-4000
                // ```c
                // Blast_IntervalTreeReset(tree);
                // /* remove all HSPs computed for the current query */
                // for (index2 = 0; index2 < hsp_list->hspcnt; index2++) { ... }
                // ```
                let accepted_for_query = preliminary_hsp_counts_by_query
                    .get(query_index)
                    .copied()
                    .unwrap_or(0);
                let (removed_hsps, rebuild_nodes) = if accepted_for_query == 0 {
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3975-3997
                    // ```c
                    // Blast_IntervalTreeReset(tree);
                    // /* remove all HSPs computed for the current query */
                    // for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
                    //     ...
                    // }
                    // ```
                    //
                    // With no accepted HSP for this query, NCBI's rebuild
                    // would copy every existing non-query HSP back into the
                    // tree. The current tree already contains exactly that set.
                    (0, 0)
                } else {
                    let removed_hsps = remove_preliminary_hits_for_query(
                        &mut preliminary_hsp_list,
                        query_index,
                        contexts,
                    );
                    debug_assert_eq!(
                        removed_hsps, accepted_for_query,
                        "NCBI redo removal should drop all accepted HSPs for the retried query"
                    );
                    if let Some(count) = preliminary_hsp_counts_by_query.get_mut(query_index) {
                        *count = 0;
                    }
                    let rebuild_nodes = rebuild_preliminary_interval_tree_from_list(
                        interval_tree,
                        &preliminary_hsp_list,
                        contexts,
                    );
                    (removed_hsps, rebuild_nodes)
                };
                if let Some(timing) = timing.as_ref() {
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3975-4007
                    // ```c
                    // Blast_IntervalTreeReset(tree);
                    // /* remove all HSPs computed for the current query */
                    // ...
                    // redo_index = index;
                    // redo_query = query_index;
                    // ```
                    BlastpTiming::record_count(&timing.preliminary_exact_retry_count, 1);
                    BlastpTiming::record_count(
                        &timing.preliminary_retry_removed_hsps,
                        removed_hsps,
                    );
                    if accepted_for_query != 0 {
                        BlastpTiming::record_count(&timing.preliminary_interval_tree_rebuilds, 1);
                        BlastpTiming::record_count(
                            &timing.preliminary_interval_tree_rebuild_nodes,
                            rebuild_nodes,
                        );
                    }
                }
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4001-4007
                // ```c
                // redo_index = index;
                // redo_query = query_index;
                // index = -1;
                // continue;
                // ```
                redo_index = Some(hit_index);
                redo_query = Some(query_index);
                if first_init_hsp_index_by_query.is_empty() {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3870-3878
                    // ```c
                    // if (index < redo_index && query_index != redo_query) {
                    //     continue;
                    // }
                    // ```
                    //
                    // NCBI restarts at the first initial hit and skips until
                    // it finds the retried query. Rust precomputes that first
                    // per-query hit lazily and jumps to the same first
                    // processable index, preserving the subsequent hit order.
                    first_init_hsp_index_by_query.resize(query_ids.len(), None);
                    next_same_query_index.clear();
                    next_same_query_index.resize(init_hsps.len(), None);
                    let mut last_init_hsp_index_by_query = vec![None; query_ids.len()];
                    for (index, init_hsp) in init_hsps.iter().enumerate() {
                        let context_idx = blastp_get_ungapped_hsp_context_with_bounds(
                            contexts,
                            context_lookup_bounds,
                            init_hsp,
                        );
                        let init_query_index = usize::try_from(contexts[context_idx].q_idx)
                            .expect("NCBI BLAST requires blastp query indices to fit in usize");
                        if let Some(slot) = first_init_hsp_index_by_query.get_mut(init_query_index)
                        {
                            if slot.is_none() {
                                *slot = Some(index);
                            }
                        }
                        if let Some(last_index) =
                            last_init_hsp_index_by_query.get_mut(init_query_index)
                        {
                            if let Some(previous_index) = *last_index {
                                next_same_query_index[previous_index] = Some(index);
                            }
                            *last_index = Some(index);
                        }
                    }
                }
                hit_index = first_init_hsp_index_by_query
                    .get(query_index)
                    .and_then(|index| *index)
                    .unwrap_or(0);
                continue;
            }
            if preliminary_hsp.raw_score < cutoff {
                hit_index = blastp_advance_preliminary_hit_index(
                    hit_index,
                    redo_index,
                    &next_same_query_index,
                );
                continue;
            }
            interval_tree.add_hsp(
                preliminary_tree_hsp(&preliminary_hsp, ctx.frame_base),
                ctx.frame_base,
                IndexMethod::QueryAndSubject,
            );
            preliminary_hsp_list.push(preliminary_hsp);
            if let Some(count) = preliminary_hsp_counts_by_query.get_mut(query_index) {
                *count += 1;
            }
            hit_index =
                blastp_advance_preliminary_hit_index(hit_index, redo_index, &next_same_query_index);
        }
        if let Some(timing) = timing.as_ref() {
            BlastpTiming::record_call(
                &timing.preliminary_gapped_ns,
                &timing.preliminary_gapped_calls,
                preliminary_gapped_start,
            );
        }
        scratch.init_hsps = init_hsps;
        x_ih_pool::at_subject_end(&mut scratch.init_hsps);
        scratch.init_hsps.clear();

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:540-555
        // ```c
        // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
        // ...
        // Blast_HSPListSortByScore(hsp_list);
        // ```
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2525
        // ```c
        // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
        // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
        // ```
        let endpoint_duplicates_purged =
            prepare_preliminary_hits_for_kappa_in_place(&mut preliminary_hsp_list);
        if let Some(timing) = timing.as_ref() {
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2525
            // ```c
            // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
            // ...
            // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
            // ```
            BlastpTiming::record_count(
                &timing.preliminary_endpoint_duplicates_purged,
                endpoint_duplicates_purged,
            );
        }
        if x_hits_flat::x_bphitsflat() {
            x_hits_flat::x_lists_by_query(
                &mut preliminary_hsp_list,
                contexts,
                query_ids.len(),
                &mut scratch.x_hits_order,
                &mut scratch.x_hits_starts,
                s_idx as u32,
                subject.aa_len,
                gapped_params,
                gapped_gumbel.as_ref(),
                values.search_spaces,
                values.prelim_evalue,
                values.sum_stats.map(|sum_stats| &sum_stats.preliminary),
                subject_result,
            );
        } else {
            let mut preliminary_hits_by_query =
                std::mem::take(&mut scratch.preliminary_hits_by_query);
            split_preliminary_hits_by_query_into(
                &mut preliminary_hsp_list,
                contexts,
                query_ids.len(),
                &mut preliminary_hits_by_query,
            );
            for (q_idx, preliminary_hits) in preliminary_hits_by_query
                .iter_mut()
                .take(query_ids.len())
                .enumerate()
            {
                if preliminary_hits.is_empty() {
                    continue;
                }
                let ctx = &contexts[q_idx];
                // The compact list, with the E-values of `Blast_HSPListGetEvalues` and the reap by the
                // preliminary E-value (`x_prelim_compact.rs`); a list that keeps HSPs gives back the
                // room of the reaped ones (`x_prelim_free.rs`).
                let mut compact_hsp_list = x_build_compact_preliminary_hsp_list(
                    preliminary_hits.iter(),
                    s_idx as u32,
                    q_idx as u32,
                    ctx.aa_len,
                    subject.aa_len,
                    gapped_params,
                    gapped_gumbel.as_ref(),
                    &values.search_spaces[q_idx],
                    values.sum_stats.map(|sum_stats| &sum_stats.preliminary),
                );
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1523-1539
                // ```c
                // Blast_HSPListGetEvalues(program_number, query_info,
                //                         stat_length, hsp_list,
                //                         gapped_calculation, FALSE,
                //                         sbp, 0, 1.0);
                // status = s_Blast_HSPListReapByPrelimEvalue(hsp_list, hit_params);
                // ```
                //
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:640-670
                // ```c
                // cutoff = hit_params->prelim_evalue;
                // ...
                // if (hsp->evalue > cutoff) {
                //     hsp_array[index] = Blast_HSPFree(hsp_array[index]);
                // }
                // ```
                reap_hsplist_by_evalue(&mut compact_hsp_list, values.prelim_evalue);
                if !compact_hsp_list.hsps.is_empty() {
                    compact_hsp_list.hsps.shrink_to_fit();
                    subject_result
                        .x_compact_hsp_lists
                        .push((q_idx, compact_hsp_list));
                }
                preliminary_hits.clear();
            }
            scratch.preliminary_hits_by_query = preliminary_hits_by_query;
        }
        scratch.restricted_align_array = restricted_align_array;
        scratch.restricted_align_array.clear();
        scratch.preliminary_hsp_counts_by_query = preliminary_hsp_counts_by_query;
        scratch.preliminary_hsp_counts_by_query.clear();
        scratch.preliminary_hsp_list = preliminary_hsp_list;
    }
}
