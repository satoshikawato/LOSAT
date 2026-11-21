//! The preliminary HSP list of one blastn subject chunk, made from the chunk's initial (ungapped)
//! hits once its word finder is done. A gapped search sorts the hits by score, or groups them by
//! query (`sort_ungapped_hits_for_gapped`), then (`gapped_hsp_list`) runs the loop of the gapped
//! preliminary search over one list of them (`gapped_score_list`: each hit tested against the
//! interval tree, extended score-only and, at or above its context's hit-saving cutoff, put into
//! the tree and the list) per query on the pool (LOSAT_X_QGROUPPAR) or over the chunk's list, then
//! purges and sorts the list and moves it to the subject's coordinates. That extension as a
//! function of the hit alone (`x_prelim_extend`) is what the look-ahead helpers compute ahead of
//! the loop. An ungapped search (`ungapped_hsp_list`) makes the list from the sorted hits.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:503-578 (the chunk loop of
//! `s_BlastSearchEngineOneContext` after the word finder: `GetGappedScore`, which is
//! `BLAST_GetGappedScore`, c++/src/algo/blast/core/blast_gapalign.c:3664-4114, with its loop on
//! the initial hits at :3850; the purge and the sort; or `BLAST_GetUngappedHSPList`,
//! blast_gapalign.c:4719-4775; the offsets)
//!
//! This module is a child of `run` for the types it keeps private (`UngappedHit`, `PrelimHit`,
//! `SubjectChunk`, `XAheadScratch`, the sorts). Its functions are the code of
//! `prelim_subject::prelim_chunk_hits` after the scan loop: `gapped_score_list` and
//! `x_prelim_extend` the bodies of what were the closures `x_prelim_list` and `x_prelim_extend`,
//! and `ungapped_hsp_list`, `sort_ungapped_hits_for_gapped` and `gapped_hsp_list` the code around
//! them. The values they read are the fields of `GappedPrelim` (the batch's
//! `prelim_subject::PrelimBatch` and the chunk's own) and `GappedDrive`; what the loop works on is
//! `GappedListState`. The query group tasks of `x_qgroup_par` run `gapped_score_list`, and the
//! look-ahead helpers `x_prelim_extend`, on the threads of the pool, so what they read is `Sync`.

use super::prelim_subject::PrelimBatch;
use super::x_greedy_prelim::{GreedyPrelim, XPrelimSubject};
use super::{
    blastn_trace, context_query, extend_gapped_heuristic_with_scratch,
    purge_prelim_hits_with_common_endpoints, score_compare_prelim_hits,
    score_compare_ungapped_hits, sort_prelim_hits_by_score, x_bn_mapped, x_bn_mapped_reach,
    x_qgroup_par, BlastIntervalTree, FastaRecord, GapAlignScratch, GreedyAlignScratch, IndexMethod,
    PrelimHit, SubjectChunk, TreeHsp, UngappedHit, XAheadScratch, COMPRESSION_RATIO, X_SORT,
};
use std::sync::atomic::AtomicU64;

/// The HSP list of one subject chunk of an ungapped (`-ungapped`) search, made from its initial
/// hits: the hits sorted by score and counted, then the list in the subject's coordinates, sorted
/// by score.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:503,522-526,560-563
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4719,4742,4760-4772
/// An ungapped (`-ungapped`) search sorts the initial hits by score first, as a gapped
/// search does (`BlastNaWordFinder` ends with `if (word_params->ungapped_extension)
/// Blast_InitHitListSortByScore(init_hitlist)`, na_ungapped.c:1691-1692, and
/// `ungapped_extension` is TRUE for every program but PHI: blast_parameters.c:189-192),
/// with the same comparator (`score_compare_ungapped_hits` = `score_compare_match`:
/// score descending, s_start, length descending, the concatenated query's q_start;
/// stable, as the oracle's glibc qsort). The order of ties (a plus and a minus hit
/// with the same score and context-relative ranges) then comes out as NCBI's, and
/// `Blast_HSPListSortByScore` below keeps it. Then the HSP list of the chunk is made
/// from those hits (no gapped extension), the subject chunk's offset is taken
/// (blast_engine.c:564-580), and there is no purge of common end points (a gapped
/// search's step).
pub(super) fn ungapped_hsp_list(
    ungapped_hits: &mut [UngappedHit],
    chunk: &SubjectChunk,
    good_init_extends: &AtomicU64,
) -> Vec<PrelimHit> {
    // The same sort as the gapped search's (`sort_ungapped_hits_for_gapped`, review 2 S5): with
    // LOSAT_X_BNMAPPED (default ON) the scratch of the stable sort is a mapped buffer, not a plain
    // `sort_by` scratch of up to 8 MB left in the thread's arena. The comparator and the stability
    // are the same, so the order is `sort_by`'s. LOSAT_X_QGROUPPAR is not used for an ungapped
    // search.
    if x_bn_mapped() {
        x_bn_mapped_reach(
            ungapped_hits
                .len()
                .saturating_mul(std::mem::size_of::<UngappedHit>()),
        );
        crate::utils::x_mapped::x_sort_by_mapped(
            ungapped_hits,
            score_compare_ungapped_hits,
            &X_SORT,
        )
    } else {
        ungapped_hits.sort_by(score_compare_ungapped_hits)
    }
    good_init_extends.fetch_add(
        ungapped_hits.len() as u64,
        std::sync::atomic::Ordering::Relaxed,
    );
    let mut chunk_prelim_hits: Vec<PrelimHit> = ungapped_hits
        .iter()
        .map(|uh| {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4760-4767
            // NCBI's HSP of a hit that passes an end of its context (`UngappedHit::saved`) has
            // the hit's `Int4` offsets, outside the context.
            PrelimHit {
                context_idx: uh.context_idx,
                query_idx: uh.query_idx,
                query_frame: uh.query_frame,
                query_context_offset: uh.query_context_offset,
                prelim_qs: uh.qs,
                prelim_qe: uh.qe,
                prelim_ss: (uh.ss as usize).saturating_add(chunk.offset),
                prelim_se: (uh.se as usize).saturating_add(chunk.offset),
                prelim_score: uh.score,
                seed_qs: uh.seed_q_off as i32,
                seed_ss: (uh.seed_s_off as usize).saturating_add(chunk.offset),
                prelim_evalue: 0.0,
            }
        })
        .collect();
    sort_prelim_hits_by_score(&mut chunk_prelim_hits);
    chunk_prelim_hits
}

/// The initial hits of one subject chunk of a gapped search, sorted by score for the gapped loop
/// or grouped by query for its query group tasks, and counted: returns `(x_counts, x_order,
/// total_ungapped)`, where `x_counts` (the hits of each query) is `Some` when the hits are grouped
/// by query (LOSAT_X_QGROUPPAR on, the owner's serial call), and `x_order` (the order of the tasks)
/// is `Some` only when the loop then runs per query (`x_qgroup_par::plan`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1688-1692
/// The end of NCBI's word finder (`BlastNaWordFinder`): the hits counted, then sorted. No NCBI
/// counterpart for the query groups (scheduling only).
pub(super) fn sort_ungapped_hits_for_gapped(
    ungapped_hits: &mut [UngappedHit],
    x_qgroup_batch: bool,
    reuse_prelim_hits: bool,
    queries: &[FastaRecord],
    good_init_extends: &AtomicU64,
) -> (Option<Vec<u32>>, Option<Vec<u32>>, usize) {
    // ============================================================
    // NCBI ARCHITECTURE: Batch process ungapped hits in score order
    // ============================================================
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3824 - ASSERT(Blast_InitHitListIsSortedByScore(init_hitlist))
    // NCBI processes init_hsp_array sorted by score DESCENDING
    // High-score HSPs are processed first, their GAPPED HSPs added to interval tree
    // Lower-score HSPs are then checked for containment against these GAPPED HSPs

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-310
    // The oracle's qsort (glibc 2.39) is a stable merge sort: hits that
    // compare equal keep the order in which they were saved. `sort_by`
    // is stable too.
    // EXPERIMENT (LOSAT_X_QGROUPPAR): with the stage run per query, each query's hits are
    // taken in the order the seed stage saved them and sorted by the task (the comparator and
    // the sort are those below, so each query's hits come in the relative order the batch's
    // sort gives them).
    let x_counts = (x_qgroup_batch && reuse_prelim_hits).then(|| {
        x_qgroup_par::counts(ungapped_hits.iter().map(|hit| hit.query_idx), queries.len())
    });
    let x_order = x_counts
        .as_deref()
        .and_then(|counts| x_qgroup_par::plan(counts, x_qgroup_par::Stage::Prelim));
    match (&x_order, &x_counts) {
        (Some(_), Some(counts)) => {
            x_qgroup_par::group_in_place(ungapped_hits, |hit| hit.query_idx, counts)
        }
        // EXPERIMENT (LOSAT_X_BNMAPPED): `x_mapped::x_sort_by_mapped` is a stable sort
        // with the same comparator, so the order is that of `sort_by`.
        _ if x_bn_mapped() => {
            x_bn_mapped_reach(
                ungapped_hits
                    .len()
                    .saturating_mul(std::mem::size_of::<UngappedHit>()),
            );
            crate::utils::x_mapped::x_sort_by_mapped(
                ungapped_hits,
                score_compare_ungapped_hits,
                &X_SORT,
            )
        }
        _ => ungapped_hits.sort_by(score_compare_ungapped_hits),
    }

    let total_ungapped = ungapped_hits.len();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1688-1689
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_diagnostics.c:106-112
    good_init_extends.fetch_add(total_ungapped as u64, std::sync::atomic::Ordering::Relaxed);
    (x_counts, x_order, total_ungapped)
}

/// What the gapped preliminary loop of one subject chunk (`gapped_score_list`) and one hit's
/// extension (`x_prelim_extend`) read: the batch's (`PrelimBatch`, or its copy per subject: the
/// trace and verbose switches, the queries' contexts, the scoring, the look-ahead's schedule), and
/// the chunk's own: the subject (its index, its trace id, its letters for the greedy aligner), the
/// chunk (its offset, packed letters and length), whether the chunk is the owner's serial call (not
/// one of a chunk batch on the pool) and the chunk's number of initial hits. Built once per chunk
/// by `prelim_subject::prelim_chunk_hits`, where the closures were; every field is a shared
/// reference or a `Copy` value, so the struct is `Sync` for the query group tasks and the
/// look-ahead helpers on the pool.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3664-3675
/// NCBI's function takes the query, the subject and the parameters; the look-ahead's fields
/// have no NCBI counterpart (scheduling only).
pub(super) struct GappedPrelim<'a> {
    /// The batch's values: the trace lines (`blastn_trace`) and the verbose lines; the queries'
    /// contexts (their letters in blastna, their lengths and frames, and their hit-saving cutoffs,
    /// by context); the scoring, the aligner (dynamic programming or greedy) and its X-dropoff,
    /// and the interval tree's diagonal separation; the look-ahead window (`None` when the helpers
    /// do not run), the rule of LOSAT_X_GREEDYSPEC (kept as its `x_reach!` key) and the threads of
    /// the pool.
    pub(super) batch: &'a PrelimBatch<'a>,
    /// The subject's index in the batch, its trace id, and its letters one per byte for the
    /// greedy aligner (`XPrelimSubject`).
    pub(super) s_idx: usize,
    pub(super) s_id: &'a str,
    pub(super) x_prelim_subject: &'a XPrelimSubject<'a>,
    /// The chunk (its offset in the subject), its packed letters and its length.
    pub(super) chunk: &'a SubjectChunk,
    pub(super) s_seq_packed: &'a [u8],
    pub(super) s_len: usize,
    /// Whether this is the owner's serial call for the chunk (the look-ahead's condition), and the
    /// chunk's number of initial hits (the verbose line).
    pub(super) reuse_prelim_hits: bool,
    pub(super) total_ungapped: usize,
}

/// What the gapped preliminary loop over one list works on and fills: the interval tree, the
/// scratch of its score-only calls, the scratch list of the look-ahead helpers (`None` for a query
/// group task, which runs without helpers) and the list that the hits it keeps go to. The owner's
/// tree, greedy scratch, helper list and hit list come from the subject's scratch
/// (`SubjectScratch`), its gap scratch from `prelim_subject::prelim_chunk_hits`; a query group
/// task's tree and scratch come from its state, its list from the task.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3830-3831,4084-4088
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
) -> Option<GreedyPrelim> {
    let &GappedPrelim {
        batch,
        x_prelim_subject,
        chunk,
        s_seq_packed,
        s_len,
        ..
    } = gapped;
    let &PrelimBatch {
        encoded_query_concat_blastna,
        encoded_queries_blastna,
        query_contexts,
        reward,
        penalty,
        gap_open,
        gap_extend,
        score_matrix,
        use_dp,
        x_drop_gapped,
        ..
    } = batch;
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
        // The dynamic-programming alignment runs on the context (`query_tmp`), from a seed in it.
        Some((
            p_qs as i32,
            p_qe as i32,
            p_ss,
            p_se,
            p_score,
            seed_qs as i32,
            seed_ss,
        ))
    } else {
        // NCBI's `Int4` seed, outside the context for a hit that passes an end of it
        // (`UngappedHit::greedy_seed`); the alignment reads the concatenated query there
        // (`greedy_gapped_alignment_internal`).
        let (seed_qs, seed_ss) = uh.greedy_seed();
        x_prelim_subject.score_only(
            chunk.offset,
            s_len,
            (uh.ss as usize, uh.se as usize),
            context_query(
                encoded_query_concat_blastna,
                &query_contexts[uh.context_idx as usize],
            ),
            seed_qs,
            seed_ss as usize,
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
        batch,
        s_idx,
        s_id,
        x_prelim_subject,
        chunk,
        s_seq_packed,
        s_len,
        reuse_prelim_hits,
        total_ungapped,
    } = gapped;
    let &PrelimBatch {
        blastn_trace_enabled,
        verbose,
        encoded_query_concat_blastna,
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
        x_ahead_window,
        x_greedy_spec,
        num_threads,
        ..
    } = batch;
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
    // `with_helpers` runs the closure on this thread. The closure is the port of this NCBI loop, in
    // NCBI order; the other threads of the pool only evaluate values ahead (`utils/xahead.rs`).
    // Process each ungapped hit in score order
    crate::utils::xahead::with_helpers(x_ahead.as_ref(), x_helpers, &x_prelim_compute, || {
        for (idx, uh) in ungapped_hits.iter().enumerate() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3918-3919
            // The owner offers later indexes to the helpers, leaving out the ones the tree already contains.
            // The tree is only read here. A wrong guess costs time only, because the loop decides each index
            // itself.
            // No NCBI counterpart: offer the next hits to the helpers,
            // leaving out those the tree already contains.
            if let Some(x_ahead) = x_ahead.as_ref() {
                x_ahead.announce(idx, x_ungapped.len(), |x_index| {
                    let u = &x_ungapped[x_index];
                    let h = TreeHsp {
                        query_offset: u.qs,
                        query_end: u.qe,
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3886-4090 (no per-HSP logging)
            if verbose && idx % 100 == 0 {
                eprintln!("[INFO] Gapped extension: {}/{}", idx, total_ungapped);
            }
            // Get query sequence for this hit (context-specific)
            let ctx = &query_contexts[uh.context_idx as usize];
            let q_seq_blastna = encoded_queries_blastna[uh.context_idx as usize];
            let subject_len = s_len;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2762-2793
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-507 (traceback uses uncompressed subject)
            // The dynamic-programming score-only call below reads the packed subject
            // (`s_seq_score`); the greedy one runs on the bases one per byte
            // (`x_prelim_subject`), NCBI's `compressed_subject = FALSE` form.
            let s_seq_score = s_seq_packed;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3908-3914
            // NCBI uses 0-based coordinates internally for the tree: the `Int4` offsets in the
            // context, outside it for a hit that passes an end of it (`UngappedHit::saved`), and
            // so for a greedy alignment from such a hit (`x_prelim_extend`); the tree holds them
            // at the context's offset (`BlastIntervalTreeAddHSP`), inside the concatenated query.
            let subject_frame_sign = 1i32;
            let ungapped_tree_hsp = TreeHsp {
                query_offset: uh.qs,
                query_end: uh.qe,
                subject_offset: uh.ss as i32,
                subject_end: uh.se as i32,
                score: uh.score,
                query_frame: uh.query_frame,
                query_length: ctx.len as i32,
                query_context_offset: uh.query_context_offset,
                subject_frame_sign,
            };

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3918 BlastIntervalTreeContainsHSP
            // Check if UNGAPPED HSP is contained in existing GAPPED HSPs
            let containing_hsp = interval_tree.containing_hsp(
                &ungapped_tree_hsp,
                uh.query_context_offset,
                min_diag_separation,
            );
            let is_contained = containing_hsp.is_some();

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3908-3919
            if blastn_trace_enabled
                && blastn_trace::should_trace_range(
                    "prelim",
                    uh.context_idx,
                    s_idx,
                    s_id,
                    uh.qs,
                    uh.qe,
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
                    uh.qs,
                    uh.qe,
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
                // When the tree contains the hit, NCBI skips the aligner call for this index. The helpers are
                // told that nobody needs its value.
                if let Some(x_ahead) = x_ahead.as_ref() {
                    x_ahead.skip(idx);
                }
                continue;
            }

            // Select gapped-start seed within the ungapped HSP.
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4046
            //
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3924-3927,4058
            // The score-only gapped HSP enters the preliminary HSP list and
            // interval tree only if it reaches the hit-saving cutoff, not
            // the lower BlastInitialWordParameters cutoff used for saving
            // ungapped HSPs.
            let cutoff_score = hit_saving_cutoff_scores[uh.context_idx as usize];
            let (prelim_qs, prelim_qe, prelim_ss, prelim_se, prelim_score, seed_qs, seed_ss) =
                if let Some(x_ahead) = x_ahead.as_ref() {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4025,4041-4046
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
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4033-4045
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4071-4076
                    let mut seed_qs = uh.seed_q_off as usize;
                    let mut seed_ss = uh.seed_s_off as usize;
                    if uh.se as usize >= (uh.seed_s_off as usize).saturating_add(8) {
                        seed_qs = seed_qs.saturating_add(3);
                        seed_ss = seed_ss.saturating_add(3);
                    }

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2959-2963 (x_dropoff limited by ungapped score)
                    let x_drop_score_only = x_drop_gapped.min(uh.score);
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2973-3003
                    if blastn_trace_enabled
                        && blastn_trace::should_trace_seed(
                            "prelim",
                            uh.context_idx,
                            s_idx,
                            s_id,
                            seed_qs as i32,
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
                    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:176-183
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
                    // The dynamic-programming alignment runs on the context (`query_tmp`), from a
                    // seed in it.
                    (
                        p_qs as i32,
                        p_qe as i32,
                        p_ss,
                        p_se,
                        p_score,
                        seed_qs as i32,
                        seed_ss,
                    )
                } else {
                    // Greedy seed selection (megablast)
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4017
                    // (`UngappedHit::greedy_seed`: NCBI's `Int4` seed, outside the context for a
                    // hit that passes an end of it).
                    let (seed_qs, seed_ss) = uh.greedy_seed();

                    let prelim = match x_prelim_subject.score_only(
                        chunk.offset,
                        subject_len,
                        (uh.ss as usize, uh.se as usize),
                        context_query(encoded_query_concat_blastna, ctx),
                        seed_qs,
                        seed_ss as usize,
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4058-4091
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3908-3914
            let gapped_tree_hsp = TreeHsp {
                query_offset: prelim_qs,
                query_end: prelim_qe,
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

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4031
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4058-4077
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

/// No NCBI counterpart: LOSAT's `-verbose` line of a chunk's gapped stage (its time, the
/// alignments run, the hits the tree held and the chunk's initial hits); a chunk without initial
/// hits writes it from `prelim_subject::prelim_chunk_hits` (no stage ran: zeros). Out of line,
/// off the chunk's path.
#[cold]
#[inline(never)]
pub(super) fn verbose_gapped_line(
    gapped_elapsed: std::time::Duration,
    dbg_gapped_calls: usize,
    dbg_containment_skipped: usize,
    total_ungapped: usize,
) {
    eprintln!(
        "[INFO] Gapped extension took {:?}: calls={}, skipped={}, total_ungapped={}",
        gapped_elapsed, dbg_gapped_calls, dbg_containment_skipped, total_ungapped
    );
}

/// What the gapped HSP list of one subject chunk (`gapped_hsp_list`) reads besides
/// `GappedPrelim`: the query groups' counts and order (`sort_ungapped_hits_for_gapped`; the order
/// is `None` unless the loop runs per query, LOSAT_X_QGROUPPAR) and the start of the gapped stage's
/// clock (`None` unless LOSAT_TIMING or `-verbose` reads it). The batch's queries (their number),
/// the length that bounds a query group task's interval tree (the concatenated queries; the
/// chunk's is `GappedPrelim::s_len`), the pool and the LOSAT_TIMING counters are the batch's
/// (`GappedPrelim::batch`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3830-3831
/// The bounds of each query group task's tree are those of the chunk's tree. No NCBI counterpart
/// for the query groups themselves and the clock (scheduling and diagnostics only).
pub(super) struct GappedDrive<'d> {
    /// Borrowed from the caller, which owns them (they are dropped where they were before).
    pub(super) x_counts: &'d Option<Vec<u32>>,
    pub(super) x_order: &'d Option<Vec<u32>>,
    pub(super) gapped_start_time: Option<std::time::Instant>,
}

/// The gapped HSP list of one subject chunk, from its initial hits in score order: the gapped loop
/// (`gapped_score_list`) per query on the pool (LOSAT_X_QGROUPPAR: each query's hits sorted, looped
/// with a tree and scratch of the task's own and purged; the lists concatenated in query order) or
/// over the chunk's list on this thread with the subject's tree and scratch (`state`); then the
/// LOSAT_TIMING counters and the verbose line, the list taken from `prelim_hits`, purged and sorted
/// by score (only sorted, on the pool, when the tasks purged it), and moved to the subject's
/// coordinates.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:522-526,545,555,578
/// No NCBI counterpart: the query groups on the pool (LOSAT_X_QGROUPPAR), scheduling only; NCBI
/// runs the loop of one chunk on one thread (blast_gapalign.c:3850).
pub(super) fn gapped_hsp_list(
    gapped: &GappedPrelim<'_>,
    drive: &GappedDrive<'_>,
    ungapped_hits: &mut [UngappedHit],
    state: GappedListState<'_>,
) -> Vec<PrelimHit> {
    let &GappedPrelim {
        batch,
        chunk,
        s_len,
        reuse_prelim_hits,
        total_ungapped,
        ..
    } = gapped;
    let &PrelimBatch {
        verbose,
        queries,
        query_concat_length,
        parallel_pool,
        timing_ref,
        ..
    } = batch;
    let &GappedDrive {
        x_counts,
        x_order,
        gapped_start_time,
    } = drive;
    let GappedListState {
        interval_tree,
        gap_scratch,
        greedy_align_scratch,
        x_helper_scratch,
        prelim_hits,
    } = state;
    // Debug counters for containment analysis
    let mut dbg_containment_skipped = 0usize;
    let mut dbg_gapped_calls = 0usize;
    // Dispatch point of LOSAT_X_QGROUPPAR for the preliminary loop: per query on the
    // pool, each query's purged list concatenated in query order (the serial purge's
    // output, see x_qgroup_par.rs), or the loop over the chunk's list.
    let mut x_purged = false;
    if let Some(order) = x_order.as_ref() {
        let x_starts = x_qgroup_par::starts(x_counts.as_deref().unwrap_or(&[]));
        let mut x_rest: &mut [UngappedHit] = ungapped_hits;
        let mut x_slices = Vec::with_capacity(queries.len());
        for query in 0..queries.len() {
            let (head, tail) = x_rest.split_at_mut(x_starts[query + 1] - x_starts[query]);
            x_slices.push(std::sync::Mutex::new(head));
            x_rest = tail;
        }
        let x_tree_q_max = (query_concat_length + 1) as i32;
        let x_tree_s_max = (s_len + 1) as i32;
        let x_results = x_qgroup_par::run(
            parallel_pool,
            x_qgroup_par::Stage::Prelim,
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
                let mut hits = x_slices[query as usize]
                    .lock()
                    .expect("each query's hits are taken by one task");
                hits.sort_by(score_compare_ungapped_hits);
                tree.reset();
                let mut out = Vec::new();
                let dbg = gapped_score_list(
                    gapped,
                    &hits,
                    GappedListState {
                        interval_tree: tree,
                        gap_scratch: gap,
                        greedy_align_scratch: greedy,
                        x_helper_scratch: None,
                        prelim_hits: &mut out,
                    },
                );
                (purge_prelim_hits_with_common_endpoints(out), dbg)
            },
        );
        let mut x_hits: Vec<PrelimHit> = Vec::new();
        for (hits, (calls, skipped)) in x_results.into_iter().flatten() {
            dbg_gapped_calls += calls;
            dbg_containment_skipped += skipped;
            x_hits.extend(hits);
        }
        *prelim_hits = x_hits;
        x_purged = true;
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3706-3707,3826-3831
        // query->length is the concatenated length for both strands (2*len + 1). The tree is
        // the subject scratch's (`x_tree`), given this chunk's ranges here, inside the gapped
        // stage of a chunk with initial hits (a chunk without one returns before it,
        // `prelim_subject::prelim_chunk_hits`); the query group tasks above set up their own.
        interval_tree.reset_with_bounds(
            0,                                // q_min
            (query_concat_length + 1) as i32, // q_max
            0,                                // s_min
            (s_len + 1) as i32,               // s_max
        );
        let (calls, skipped) = gapped_score_list(
            gapped,
            ungapped_hits,
            GappedListState {
                interval_tree,
                gap_scratch,
                greedy_align_scratch,
                x_helper_scratch,
                prelim_hits,
            },
        );
        dbg_gapped_calls = calls;
        dbg_containment_skipped = skipped;
    }

    // DEBUG: Log gapped extension stats (the clock is read only when LOSAT_TIMING or `-verbose`
    // is on, `prelim_subject::prelim_chunk_hits`).
    let gapped_elapsed = gapped_start_time
        .map(|start| start.elapsed())
        .unwrap_or_default();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3917-3922
    if let Some(timing) = timing_ref {
        timing.gapped_ns.fetch_add(
            gapped_elapsed.as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
        timing.gapped_calls.fetch_add(
            dbg_gapped_calls as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
    }
    if verbose {
        verbose_gapped_line(
            gapped_elapsed,
            dbg_gapped_calls,
            dbg_containment_skipped,
            total_ungapped,
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:488-491
    #[expect(
        clippy::drain_collect,
        reason = "without reuse the scratch keeps its buffer for the next chunk; mem::take would hand the buffer away and allocate again"
    )]
    let mut chunk_prelim_hits: Vec<PrelimHit> = if reuse_prelim_hits {
        std::mem::take(prelim_hits)
    } else {
        prelim_hits.drain(..).collect()
    };
    if !chunk_prelim_hits.is_empty() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4058-4088
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:540-557
        // NCBI purges redundant preliminary score-only gapped HSPs
        // before subject chunk offset adjustment and final traceback.
        // With LOSAT_X_QGROUPPAR the queries' lists were purged by their tasks, and the
        // same stable sort runs on the pool.
        if x_purged {
            x_qgroup_par::sort_stable(
                parallel_pool,
                &mut chunk_prelim_hits,
                score_compare_prelim_hits,
            );
        } else {
            chunk_prelim_hits = purge_prelim_hits_with_common_endpoints(chunk_prelim_hits);
            sort_prelim_hits_by_score(&mut chunk_prelim_hits);
        }
        for hit in chunk_prelim_hits.iter_mut() {
            hit.prelim_ss = hit.prelim_ss.saturating_add(chunk.offset);
            hit.prelim_se = hit.prelim_se.saturating_add(chunk.offset);
            hit.seed_ss = hit.seed_ss.saturating_add(chunk.offset);
        }
    }
    chunk_prelim_hits
}
