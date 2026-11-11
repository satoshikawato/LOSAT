//! The preliminary HSP list of one blastn subject chunk, made from the chunk's initial (ungapped)
//! hits: the loop of the gapped preliminary search over one list of them (`gapped_score_list`),
//! each hit tested against the interval tree, extended score-only and, at or above its context's
//! hit-saving cutoff, put into the tree and the list; and that extension as a function of the hit
//! alone (`x_prelim_extend`), which the look-ahead helpers compute ahead of the loop.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3664-4114
//! (`BLAST_GetGappedScore`), of which :3850 is the loop on the initial hits
//! ```c
//! Int2 BLAST_GetGappedScore (EBlastProgramType program_number,
//! ...
//!    for (index=0; index<init_hitlist->total; index++)
//! ```
//!
//! This module is a child of `run` for the types it keeps private (`UngappedHit`, `PrelimHit`,
//! `SubjectChunk`, `QueryContext`, `XAheadScratch`). `gapped_score_list` is the body of what was
//! the closure `x_prelim_list` of `search_query_batch`'s `collect_prelim_hits_for_chunk`, and
//! `x_prelim_extend` that of the closure of the same name. The values they read are the fields of
//! `GappedPrelim`; what the loop works on is `GappedListState`. The query group tasks of
//! `x_qgroup_par` run `gapped_score_list`, and the look-ahead helpers `x_prelim_extend`, on the
//! threads of the pool, so what they read is `Sync`.

use super::x_greedy_prelim::XPrelimSubject;
use super::{
    blastn_trace, extend_gapped_heuristic_with_scratch, BlastIntervalTree, GapAlignScratch,
    GreedyAlignScratch, IndexMethod, PrelimHit, QueryContext, SubjectChunk, TreeHsp, UngappedHit,
    XAheadScratch, COMPRESSION_RATIO,
};
use crate::algorithm::blastn::alignment::gapped::BlastnaMatrix;

/// What the gapped preliminary loop of one subject chunk (`gapped_score_list`) and one hit's
/// extension (`x_prelim_extend`) read: the trace and verbose switches, the subject (its index, its
/// trace id, its letters for the greedy aligner), the chunk (its offset, packed letters and
/// length), the queries' contexts, the scoring, the call (the owner's or a query group task's) and
/// the chunk's number of initial hits, and the look-ahead's schedule. Built once per chunk by
/// `collect_prelim_hits_for_chunk`, where the closures were; every field is a shared reference or
/// a `Copy` value, so the struct is `Sync` for the query group tasks and the look-ahead helpers on
/// the pool.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3664-3675
/// ```c
/// Int2 BLAST_GetGappedScore (EBlastProgramType program_number,
///         BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
///         BLAST_SequenceBlk* subject,
///         BlastGapAlignStruct* gap_align,
///         const BlastScoringParameters* score_params,
///         const BlastExtensionParameters* ext_params,
///         const BlastHitSavingParameters* hit_params,
///         const BlastInitialWordParameters* word_params,
///         BlastInitHitList* init_hitlist,
///         BlastHSPList** hsp_list_ptr, BlastGappedStats* gapped_stats,
///         Boolean * fence_hit)
/// ```
/// NCBI's function takes the query, the subject and the parameters; the look-ahead's fields
/// have no NCBI counterpart (scheduling only).
pub(super) struct GappedPrelim<'a> {
    /// The trace lines (`blastn_trace`) and the verbose lines.
    pub(super) blastn_trace_enabled: bool,
    pub(super) verbose: bool,
    /// The subject's index in the batch, its trace id, and its letters one per byte for the
    /// greedy aligner (`XPrelimSubject`).
    pub(super) s_idx: usize,
    pub(super) s_id: &'a str,
    pub(super) x_prelim_subject: &'a XPrelimSubject<'a>,
    /// The chunk (its offset in the subject), its packed letters and its length.
    pub(super) chunk: &'a SubjectChunk,
    pub(super) s_seq_packed: &'a [u8],
    pub(super) s_len: usize,
    /// The queries' contexts: their letters in blastna, their lengths and frames, and their
    /// hit-saving cutoffs, by context.
    pub(super) encoded_queries_blastna: &'a [&'a [u8]],
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) hit_saving_cutoff_scores: &'a [i32],
    /// The scoring, the aligner (dynamic programming or greedy) and its X-dropoff, and the
    /// interval tree's diagonal separation.
    pub(super) reward: i32,
    pub(super) penalty: i32,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    pub(super) score_matrix: &'a BlastnaMatrix,
    pub(super) use_dp: bool,
    pub(super) x_drop_gapped: i32,
    pub(super) min_diag_separation: i32,
    /// Whether this is the owner's serial call for the chunk (the look-ahead's condition), and the
    /// chunk's number of initial hits (the verbose line).
    pub(super) reuse_prelim_hits: bool,
    pub(super) total_ungapped: usize,
    /// The look-ahead window (`None` when the helpers do not run), the rule of LOSAT_X_GREEDYSPEC
    /// (kept as its `x_reach!` key) and the threads of the pool.
    pub(super) x_ahead_window: Option<usize>,
    pub(super) x_greedy_spec: bool,
    pub(super) num_threads: usize,
}

/// What the gapped preliminary loop over one list works on and fills: the interval tree, the
/// scratch of its score-only calls, the scratch list of the look-ahead helpers (`None` for a query
/// group task, which runs without helpers) and the list that the hits it keeps go to. The owner's
/// tree, greedy scratch, helper list and hit list come from the subject's scratch
/// (`SubjectScratch`), its gap scratch from `collect_prelim_hits_for_chunk`; a query group task's
/// tree and scratch come from its state, its list from the task.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3830-3831,4084-4088
/// ```c
///       tree = Blast_IntervalTreeInit(0, query->length+1,
///                                     0, subject->length+1);
/// ...
///              status = Blast_HSPListSaveHSP(hsp_list, new_hsp);
/// ...
///              status = BlastIntervalTreeAddHSP(new_hsp, tree, query_info,
///                                      eQueryAndSubject);
/// ```
/// NCBI's loop works on its `tree`, the `gap_align` of its call and its `hsp_list`; the helpers'
/// scratch has no NCBI counterpart (scheduling only).
pub(super) struct GappedListState<'s> {
    pub(super) interval_tree: &'s mut BlastIntervalTree,
    pub(super) gap_scratch: &'s mut GapAlignScratch,
    pub(super) greedy_align_scratch: &'s mut GreedyAlignScratch,
    pub(super) x_helper_scratch: Option<&'s mut Vec<(GapAlignScratch, GreedyAlignScratch)>>,
    pub(super) prelim_hits: &'s mut Vec<PrelimHit>,
}

/// One ungapped hit's score-only gapped extension on the scratch `x_scratch`: the value of the
/// aligner call of the ordered loop (`gapped_score_list`) for that hit, which the look-ahead
/// helpers compute ahead of the loop on the other threads of the pool, and the loop takes
/// (`Ahead::take`) or computes itself.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2959-2963,4012-4025,4041-4046
/// ```c
///    Int4 x_dropoff =  gap_align->gap_x_dropoff;
///
///    if(init_hsp->ungapped_data->score < x_dropoff) {
///           x_dropoff = init_hsp->ungapped_data->score;
///    }
/// ...
///          } else if (is_greedy) {
///             if (init_hsp->ungapped_data) {
///                 init_hsp->offsets.qs_offsets.q_off =
///                     init_hsp->ungapped_data->q_start + init_hsp->ungapped_data->length/2;
/// ...
///             status = BLAST_GreedyGappedAlignment(
/// ...
///             if (s_end >= (Int4)init_hsp->offsets.qs_offsets.s_off + 8) {
///                init_hsp->offsets.qs_offsets.s_off += 3;
/// ...
///             status = s_BlastDynProgNtGappedAlignment(&query_tmp, subject,
/// ```
/// This function is what each branch of the ordered loop (`gapped_score_list`) computes for hit
/// `x_index`, as a function of the index alone: the seed shift by 3, the X-drop limited by the
/// ungapped score, and the aligner call, for the dynamic-programming branch and the greedy branch.
/// The branches of that loop make the same calls on the same arguments, so the function gives the
/// same value. The tree test and the tree insertion are not part of it.
/// No NCBI counterpart: the gapped extension of `gapped_score_list`'s loop
/// as a function of the hit's index (the two branches of the loop,
/// without their trace logging).
pub(super) fn x_prelim_extend(
    gapped: &GappedPrelim<'_>,
    uh: &UngappedHit,
    x_scratch: &mut XAheadScratch<'_>,
) -> Option<(usize, usize, usize, usize, i32, usize, usize)> {
    let &GappedPrelim {
        x_prelim_subject,
        chunk,
        s_seq_packed,
        s_len,
        encoded_queries_blastna,
        reward,
        penalty,
        gap_open,
        gap_extend,
        score_matrix,
        use_dp,
        x_drop_gapped,
        ..
    } = gapped;
    let q_seq_blastna = encoded_queries_blastna[uh.context_idx as usize];
    if use_dp {
        let mut seed_qs = uh.seed_q_off as usize;
        let mut seed_ss = uh.seed_s_off as usize;
        if uh.se as usize >= (uh.seed_s_off as usize).saturating_add(8) {
            seed_qs = seed_qs.saturating_add(3);
            seed_ss = seed_ss.saturating_add(3);
        }
        let x_drop_score_only = x_drop_gapped.min(uh.score);
        let (p_qs, p_qe, p_ss, p_se, p_score, _, _, _, _, _) = extend_gapped_heuristic_with_scratch(
            q_seq_blastna,
            s_seq_packed,
            s_len,
            seed_qs,
            seed_ss,
            1,
            reward,
            penalty,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop_score_only,
            &mut *x_scratch.gap,
            true,
        );
        Some((p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss))
    } else {
        let ungapped_len = uh.qe.saturating_sub(uh.qs) as usize;
        let seed_qs = uh.qs as usize + ungapped_len / 2;
        let seed_ss = uh.ss as usize + ungapped_len / 2;
        x_prelim_subject.score_only(
            chunk.offset,
            s_len,
            (uh.ss as usize, uh.se as usize),
            q_seq_blastna,
            seed_qs,
            seed_ss,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop_gapped,
            &mut *x_scratch.greedy,
        )
    }
}

/// The preliminary gapped loop over one list of ungapped hits in score order: the
/// chunk's on the owner (with the look-ahead helpers, whose scratch
/// `x_helper_scratch` holds), or one query's (LOSAT_X_QGROUPPAR, no helpers). Returns
/// the alignments run and the hits the tree held.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850
/// ```c
///    for (index=0; index<init_hitlist->total; index++)
/// ```
/// Look-ahead for the preliminary loop. When the conditions below do not hold,
/// `x_ahead` is `None` and the loop runs as ported.
/// Only the caller's own (single-chunk) loop gets helpers; chunk
/// batches already occupy the pool.
pub(super) fn gapped_score_list(
    gapped: &GappedPrelim<'_>,
    ungapped_hits: &[UngappedHit],
    state: GappedListState<'_>,
) -> (usize, usize) {
    let &GappedPrelim {
        blastn_trace_enabled,
        verbose,
        s_idx,
        s_id,
        x_prelim_subject,
        chunk,
        s_seq_packed,
        s_len,
        encoded_queries_blastna,
        query_contexts,
        hit_saving_cutoff_scores,
        reward,
        penalty,
        gap_open,
        gap_extend,
        score_matrix,
        use_dp,
        x_drop_gapped,
        min_diag_separation,
        reuse_prelim_hits,
        total_ungapped,
        x_ahead_window,
        x_greedy_spec,
        num_threads,
    } = gapped;
    let GappedListState {
        interval_tree,
        gap_scratch,
        greedy_align_scratch,
        x_helper_scratch,
        prelim_hits,
    } = state;
    let mut dbg_containment_skipped = 0usize;
    let mut dbg_gapped_calls = 0usize;
    let x_ungapped: &[UngappedHit] = ungapped_hits;
    let x_prelim_compute = |x_index: usize, x_scratch: &mut XAheadScratch<'_>| {
        x_prelim_extend(gapped, &x_ungapped[x_index], x_scratch)
    };
    let x_ahead = match x_ahead_window {
        Some(x_window)
            if x_helper_scratch.is_some()
                && reuse_prelim_hits
                && !blastn_trace_enabled
                && x_ungapped.len() >= 2 * x_window =>
        {
            Some(crate::utils::xahead::Ahead::new(x_window))
        }
        _ => None,
    };
    if x_ahead_window.is_some() {
        crate::utils::xreach::x_reach!("LOSAT_X_AHEAD", x_ahead.is_some());
    }
    if x_greedy_spec && x_ahead_window.is_some() {
        crate::utils::xreach::x_reach!("LOSAT_X_GREEDYSPEC", x_ahead.is_some());
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,4087-4088
    // ```c
    //    for (index=0; index<init_hitlist->total; index++)
    // ...
    //             status = BlastIntervalTreeAddHSP(new_hsp, tree, query_info,
    //                                     eQueryAndSubject);
    // ```
    // `with_helpers` runs the closure on this thread. The closure is the port of this NCBI loop, in
    // NCBI order; the other threads of the pool only evaluate values ahead (`utils/xahead.rs`).
    // Process each ungapped hit in score order
    crate::utils::xahead::with_helpers(x_ahead.as_ref(), x_helpers, &x_prelim_compute, || {
        for (idx, uh) in ungapped_hits.iter().enumerate() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3918-3919
            // ```c
            //    for (index=0; index<init_hitlist->total; index++)
            // ...
            //       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                         hit_options->min_diag_separation))
            // ```
            // The owner offers later indexes to the helpers, leaving out the ones the tree already contains.
            // The tree is only read here. A wrong guess costs time only, because the loop decides each index
            // itself.
            // No NCBI counterpart: offer the next hits to the helpers,
            // leaving out those the tree already contains.
            if let Some(x_ahead) = x_ahead.as_ref() {
                x_ahead.announce(idx, x_ungapped.len(), |x_index| {
                    let u = &x_ungapped[x_index];
                    let h = TreeHsp {
                        query_offset: u.qs as i32,
                        query_end: u.qe as i32,
                        subject_offset: u.ss as i32,
                        subject_end: u.se as i32,
                        score: u.score,
                        query_frame: u.query_frame,
                        query_length: query_contexts[u.context_idx as usize].len as i32,
                        query_context_offset: u.query_context_offset,
                        subject_frame_sign: 1,
                    };
                    interval_tree
                        .containing_hsp(&h, u.query_context_offset, min_diag_separation)
                        .is_none()
                });
            }
            // NCBI reference: blast_gapalign.c:3886-4090 (no per-HSP logging)
            if verbose && idx % 100 == 0 {
                eprintln!("[INFO] Gapped extension: {}/{}", idx, total_ungapped);
            }
            // Get query sequence for this hit (context-specific)
            let ctx = &query_contexts[uh.context_idx as usize];
            let q_seq_blastna = encoded_queries_blastna[uh.context_idx as usize];
            let subject_len = s_len;

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2793
            // ```c
            // if (!compressed_subject) {
            //    s = subject + s_off;
            //    rem = 4;
            // } else {
            //    s = subject + s_off/4;
            //    rem = s_off % 4;
            // }
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:503-507 (traceback uses uncompressed subject)
            // The dynamic-programming score-only call below reads the packed subject
            // (`s_seq_score`); the greedy one runs on the bases one per byte
            // (`x_prelim_subject`), NCBI's `compressed_subject = FALSE` form.
            let s_seq_score = s_seq_packed;

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3913
            // ```c
            // tmp_hsp.query.offset = q_start;
            // tmp_hsp.query.end = q_end;
            // tmp_hsp.query.frame = query_info->contexts[context].frame;
            // tmp_hsp.subject.offset = s_start;
            // tmp_hsp.subject.end = s_end;
            // ```
            // NCBI uses 0-based coordinates internally for the tree
            let subject_frame_sign = 1i32;
            let ungapped_tree_hsp = TreeHsp {
                query_offset: uh.qs as i32,
                query_end: uh.qe as i32,
                subject_offset: uh.ss as i32,
                subject_end: uh.se as i32,
                score: uh.score,
                query_frame: uh.query_frame,
                query_length: ctx.len as i32,
                query_context_offset: uh.query_context_offset,
                subject_frame_sign,
            };

            // NCBI reference: blast_gapalign.c:3918 BlastIntervalTreeContainsHSP
            // Check if UNGAPPED HSP is contained in existing GAPPED HSPs
            let containing_hsp = interval_tree.containing_hsp(
                &ungapped_tree_hsp,
                uh.query_context_offset,
                min_diag_separation,
            );
            let is_contained = containing_hsp.is_some();

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3919
            // ```c
            // tmp_hsp.query.offset = q_start;
            // tmp_hsp.query.end = q_end;
            // tmp_hsp.subject.offset = s_start;
            // tmp_hsp.subject.end = s_end;
            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                   hit_options->min_diag_separation))
            // ```
            if blastn_trace_enabled
                && blastn_trace::should_trace_range(
                    "prelim",
                    uh.context_idx,
                    s_idx,
                    s_id,
                    uh.qs as usize,
                    uh.qe as usize,
                    (uh.ss as usize).saturating_add(chunk.offset),
                    (uh.se as usize).saturating_add(chunk.offset),
                    ctx.len,
                    ctx.frame,
                )
            {
                blastn_trace::log(
                "prelim",
                format!(
                    "subject={}({}) context={} ungapped=q{}..{} s{}..{} raw_score={} tree_contains={} containing_hsp={:?} min_diag_separation={}",
                    s_id,
                    s_idx,
                    uh.context_idx,
                    uh.qs as usize,
                    uh.qe as usize,
                    (uh.ss as usize).saturating_add(chunk.offset),
                    (uh.se as usize).saturating_add(chunk.offset),
                    uh.score,
                    is_contained,
                    containing_hsp,
                    min_diag_separation
                ),
            );
            }

            if is_contained {
                // NCBI: Skip gapped extension if ungapped HSP is contained
                dbg_containment_skipped += 1;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3918-3921
                // ```c
                //       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
                //                                         hit_options->min_diag_separation))
                //       {
                // ```
                // When the tree contains the hit, NCBI skips the aligner call for this index. The helpers are
                // told that nobody needs its value.
                if let Some(x_ahead) = x_ahead.as_ref() {
                    x_ahead.skip(idx);
                }
                continue;
            }

            // Select gapped-start seed within the ungapped HSP.
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4046
            //
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3924-3927,4058
            // ```c
            // if (is_rpsblast)
            //    cutoff = hit_params->cutoffs[rps_cutoff_index].cutoff_score;
            // else
            //    cutoff = hit_params->cutoffs[context].cutoff_score;
            // ...
            // if (gap_align->score >= cutoff) {
            // ```
            // The score-only gapped HSP enters the preliminary HSP list and
            // interval tree only if it reaches the hit-saving cutoff, not
            // the lower BlastInitialWordParameters cutoff used for saving
            // ungapped HSPs.
            let cutoff_score = hit_saving_cutoff_scores[uh.context_idx as usize];
            let (prelim_qs, prelim_qe, prelim_ss, prelim_se, prelim_score, seed_qs, seed_ss) =
                if let Some(x_ahead) = x_ahead.as_ref() {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4025,4041-4046
                    // ```c
                    //          } else if (is_greedy) {
                    //             if (init_hsp->ungapped_data) {
                    //                 init_hsp->offsets.qs_offsets.q_off =
                    //                     init_hsp->ungapped_data->q_start + init_hsp->ungapped_data->length/2;
                    // ...
                    //             status = BLAST_GreedyGappedAlignment(
                    // ...
                    //             if (s_end >= (Int4)init_hsp->offsets.qs_offsets.s_off + 8) {
                    //                init_hsp->offsets.qs_offsets.s_off += 3;
                    // ...
                    //             status = s_BlastDynProgNtGappedAlignment(&query_tmp, subject,
                    // ```
                    // Look-ahead: the value of the aligner call comes from a helper or is
                    // computed here with `x_prelim_compute`. Without helpers the `else if use_dp` branch and the greedy
                    // branch below port these calls.
                    // No NCBI counterpart: the value of the branches
                    // below, taken from a helper or evaluated here.
                    let mut x_scratch = XAheadScratch {
                        gap: &mut *gap_scratch,
                        greedy: &mut *greedy_align_scratch,
                    };
                    let x_value = match x_ahead.take(idx, &mut x_scratch, &x_prelim_compute) {
                        Some(x_value) => x_value,
                        None => x_prelim_compute(idx, &mut x_scratch),
                    };
                    match x_value {
                        Some(x_value) => x_value,
                        None => continue,
                    }
                } else if use_dp {
                    // DP seed selection (blastn)
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4033-4045
                    // ```c
                    // if (s_end >= (Int4)init_hsp->offsets.qs_offsets.s_off + 8) {
                    //    init_hsp->offsets.qs_offsets.s_off += 3;
                    //    init_hsp->offsets.qs_offsets.q_off += 3;
                    // }
                    // status = s_BlastDynProgNtGappedAlignment(...);
                    // ```
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4071-4076
                    // ```c
                    // status = Blast_HSPInit(...,
                    //        init_hsp->offsets.qs_offsets.q_off,
                    //        init_hsp->offsets.qs_offsets.s_off, ...);
                    // ```
                    let mut seed_qs = uh.seed_q_off as usize;
                    let mut seed_ss = uh.seed_s_off as usize;
                    if uh.se as usize >= (uh.seed_s_off as usize).saturating_add(8) {
                        seed_qs = seed_qs.saturating_add(3);
                        seed_ss = seed_ss.saturating_add(3);
                    }

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2959-2963 (x_dropoff limited by ungapped score)
                    let x_drop_score_only = x_drop_gapped.min(uh.score);
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2973-3003
                    // ```c
                    // offset_adjustment = COMPRESSION_RATIO -
                    //      (init_hsp->offsets.qs_offsets.s_off % COMPRESSION_RATIO);
                    // q_length = init_hsp->offsets.qs_offsets.q_off + offset_adjustment;
                    // s_length = init_hsp->offsets.qs_offsets.s_off + offset_adjustment;
                    // score_left = s_BlastAlignPackedNucl(query, subject, q_length, s_length, ...);
                    // score_right = s_BlastAlignPackedNucl(query+q_length-1,
                    //    subject+(s_length+3)/COMPRESSION_RATIO - 1, ...);
                    // ```
                    if blastn_trace_enabled
                        && blastn_trace::should_trace_seed(
                            "prelim",
                            uh.context_idx,
                            s_idx,
                            s_id,
                            seed_qs,
                            seed_ss.saturating_add(chunk.offset),
                        )
                    {
                        let offset_adjustment = COMPRESSION_RATIO - (seed_ss % COMPRESSION_RATIO);
                        let mut score_q_anchor = seed_qs.saturating_add(offset_adjustment);
                        let mut score_s_anchor = seed_ss.saturating_add(offset_adjustment);
                        if score_q_anchor > q_seq_blastna.len() || score_s_anchor > subject_len {
                            score_q_anchor = score_q_anchor.saturating_sub(COMPRESSION_RATIO);
                            score_s_anchor = score_s_anchor.saturating_sub(COMPRESSION_RATIO);
                        }
                        blastn_trace::log(
                        "prelim",
                        format!(
                            "subject={}({}) context={} score_only_seed=({}, {}) local_seed_s={} offset_adjustment={} score_anchor=({}, {}) global_score_anchor_s={} x_drop_score_only={} ungapped_score={}",
                            s_id,
                            s_idx,
                            uh.context_idx,
                            seed_qs,
                            seed_ss.saturating_add(chunk.offset),
                            seed_ss,
                            offset_adjustment,
                            score_q_anchor,
                            score_s_anchor,
                            score_s_anchor.saturating_add(chunk.offset),
                            x_drop_score_only,
                            uh.score
                        ),
                    );
                    }

                    // Preliminary DP gapped extension (score-only)
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:176-183
                    let (p_qs, p_qe, p_ss, p_se, p_score, _, _, _, _, _) =
                        extend_gapped_heuristic_with_scratch(
                            q_seq_blastna,
                            s_seq_score,
                            subject_len,
                            seed_qs,
                            seed_ss,
                            1,
                            reward,
                            penalty,
                            score_matrix,
                            gap_open,
                            gap_extend,
                            x_drop_score_only,
                            gap_scratch,
                            true,
                        );
                    (p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss)
                } else {
                    // Greedy seed selection (megablast)
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4017
                    let ungapped_len = uh.qe.saturating_sub(uh.qs) as usize;
                    let seed_qs = uh.qs as usize + ungapped_len / 2;
                    let seed_ss = uh.ss as usize + ungapped_len / 2;

                    let prelim = match x_prelim_subject.score_only(
                        chunk.offset,
                        subject_len,
                        (uh.ss as usize, uh.se as usize),
                        q_seq_blastna,
                        seed_qs,
                        seed_ss,
                        reward,
                        penalty,
                        gap_open,
                        gap_extend,
                        x_drop_gapped,
                        greedy_align_scratch,
                    ) {
                        Some(value) => value,
                        None => {
                            continue;
                        }
                    };

                    let (p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss) = prelim;
                    (p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss)
                };

            dbg_gapped_calls += 1;
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4058-4091
            // ```c
            // if (gap_align->score >= cutoff) {
            //     status = Blast_HSPInit(gap_align->query_start,
            //                           gap_align->query_stop,
            //                           gap_align->subject_start,
            //                           gap_align->subject_stop, ...);
            //     status = BlastIntervalTreeAddHSP(new_hsp, tree, query_info,
            //                                      eQueryAndSubject);
            // }
            // ```
            let trace_prelim_seed = blastn_trace_enabled
                && blastn_trace::should_trace_seed(
                    "prelim",
                    uh.context_idx,
                    s_idx,
                    s_id,
                    seed_qs,
                    seed_ss.saturating_add(chunk.offset),
                );
            if blastn_trace_enabled
                && (trace_prelim_seed
                    || blastn_trace::should_trace_range(
                        "prelim",
                        uh.context_idx,
                        s_idx,
                        s_id,
                        prelim_qs,
                        prelim_qe,
                        prelim_ss.saturating_add(chunk.offset),
                        prelim_se.saturating_add(chunk.offset),
                        ctx.len,
                        ctx.frame,
                    ))
            {
                blastn_trace::log(
                "prelim",
                format!(
                    "subject={}({}) context={} seed=({}, {}) prelim=q{}..{} s{}..{} raw_score={} cutoff={} x_drop_score_only={} accepted={}",
                    s_id,
                    s_idx,
                    uh.context_idx,
                    seed_qs,
                    seed_ss.saturating_add(chunk.offset),
                    prelim_qs,
                    prelim_qe,
                    prelim_ss.saturating_add(chunk.offset),
                    prelim_se.saturating_add(chunk.offset),
                    prelim_score,
                    cutoff_score,
                    if use_dp { x_drop_gapped.min(uh.score) } else { x_drop_gapped },
                    prelim_score >= cutoff_score
                ),
            );
            }
            if prelim_score < cutoff_score {
                continue;
            }

            // Add preliminary GAPPED HSP to interval tree for containment checks
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3913
            // ```c
            // tmp_hsp.query.offset = q_start;
            // tmp_hsp.query.end = q_end;
            // tmp_hsp.query.frame = query_info->contexts[context].frame;
            // tmp_hsp.subject.offset = s_start;
            // tmp_hsp.subject.end = s_end;
            // ```
            let gapped_tree_hsp = TreeHsp {
                query_offset: prelim_qs as i32,
                query_end: prelim_qe as i32,
                subject_offset: prelim_ss as i32,
                subject_end: prelim_se as i32,
                score: prelim_score,
                query_frame: uh.query_frame,
                query_length: ctx.len as i32,
                query_context_offset: uh.query_context_offset,
                subject_frame_sign: 1,
            };
            interval_tree.add_hsp(
                gapped_tree_hsp,
                uh.query_context_offset,
                IndexMethod::QueryAndSubject,
            );

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4031
            // ```c
            // if (init_hsp->ungapped_data) {
            //    init_hsp->offsets.qs_offsets.q_off =
            //        init_hsp->ungapped_data->q_start + init_hsp->ungapped_data->length/2;
            //    init_hsp->offsets.qs_offsets.s_off =
            //        init_hsp->ungapped_data->s_start + init_hsp->ungapped_data->length/2;
            // }
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4058-4076
            // ```c
            // if (gap_align->score >= cutoff) {
            //    status = Blast_HSPInit(gap_align->query_start,
            //              gap_align->query_stop, gap_align->subject_start,
            //              gap_align->subject_stop,
            //              init_hsp->offsets.qs_offsets.q_off,
            //              init_hsp->offsets.qs_offsets.s_off, context,
            //              query_frame, subject->frame, gap_align->score,
            //              &(gap_align->edit_script), &new_hsp);
            // }
            // ```
            prelim_hits.push(PrelimHit {
                context_idx: uh.context_idx,
                query_idx: uh.query_idx,
                query_frame: uh.query_frame,
                query_context_offset: uh.query_context_offset,
                prelim_qs,
                prelim_qe,
                prelim_ss,
                prelim_se,
                prelim_score,
                seed_qs,
                seed_ss,
                prelim_evalue: 0.0,
            });
            continue;
        }
    });
    // No NCBI counterpart: LOSAT_X_STATS work counters (xstats feature); they do not change any
    // value NCBI computes.
    if let Some(x_ahead) = x_ahead.as_ref() {
        let (x_evaluated, x_taken) = x_ahead.counts();
        crate::utils::xstats::AHEAD_PRELIM_EVALUATED
            .fetch_add(x_evaluated, std::sync::atomic::Ordering::Relaxed);
        crate::utils::xstats::AHEAD_PRELIM_TAKEN
            .fetch_add(x_taken, std::sync::atomic::Ordering::Relaxed);
    }

    (dbg_gapped_calls, dbg_containment_skipped)
}
