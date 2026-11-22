//! The blastn traceback of one subject (`traceback_subject`, `search_subject`'s traceback stage):
//! the HSPs that the preliminary stage kept, the subject's letters in blastna, which the traceback
//! aligns on (`traceback_subject_blastna`), then the traceback of the kept HSPs
//! (`traceback_list::traceback_subject_hits`), or that of an ungapped (`-ungapped`) search, which
//! aligns nothing (`ungapped_traceback_hits`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1521-1713 (one subject of
//! the loop of `BLAST_ComputeTraceback_MT`), of which :1526-1528,1556-1573 fetch the subject for
//! the traceback and :1691-1696 keep the preliminary HSPs of a search without a traceback.
//!
//! This module is a child of `run` for the types it keeps private (`BlastnSubjectRef`,
//! `SubjectCtx`, `PrelimHit`, `SubjectScratch`, `SubjectHits`, the e-value helpers).
//! `traceback_subject` is the statements of `search_subject`'s traceback side, in their order;
//! what it reads of a batch is `TraceBatch` (or `special_subject`'s copy with the statistics of
//! each subject), and what `search_subject` read of the call is `SubjectCtx`. The values the
//! other two functions read are the fields of `traceback_finish::TracebackCtx`, which
//! `traceback_subject` builds.

use super::traceback_finish::TracebackCtx;
use super::{
    adjust_blastn_offsets, calculate_blastn_context_statistics, encode_iupac_to_blastna,
    traceback_list, ungapped, x_blastna_into_pool, BlastnHsp, BlastnSubjectRef, BlastnTiming,
    ContextKarlin, FastaRecord, GapAlignScratch, PrelimHit, QueryContext, SubjectCtx, SubjectHits,
    SubjectScratch,
};
use crate::algorithm::blastn::alignment::gapped::BlastnaMatrix;
use crate::utils::threading::SearchPool;
use indicatif::ProgressBar;

/// What the traceback of a subject reads of its batch, in groups below. Built once per batch by
/// `search_query_batch`, beside `prelim_subject::PrelimBatch`, with the names there (the types
/// that `TracebackCtx`, `traceback_list::TraceSchedule` and `traceback_list::TraceDrive` take);
/// with the statistics of each subject also once per subject, by `SubjectPasses::special_subject`:
/// a copy with the subject's own cutoffs (and, in an `-ungapped` search, spaces). Every field is a
/// shared reference or a `Copy` value, so the struct is `Sync`, as `search_subject`, which reads
/// it, must be.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1436-1440
/// (`BLAST_ComputeTraceback_MT`'s parameters)
/// The schedule on the pool, the progress bar and the diagnostics have no NCBI counterpart.
pub(super) struct TraceBatch<'a> {
    /// The shortest subject that the preliminary stage searches (`prelim_subject::
    /// shortest_subject`): a shorter one has no HSPs to trace.
    pub(super) shortest_subject: usize,
    /// The queries: each context's letters with and without the query's masks, each query's
    /// length, the length of the concatenated query (the interval tree's bounds), the Karlin
    /// blocks by context and whether an e-value is rounded down from its score.
    pub(super) encoded_queries_blastna: &'a Vec<&'a [u8]>,
    pub(super) encoded_queries_nomask: &'a [&'a [u8]],
    /// The same letters concatenated, with masks (`PrelimBatch::encoded_query_concat_blastna`,
    /// NCBI's `query_blk->sequence`) and without (`PrelimBatch::encoded_query_concat_nomask`,
    /// `query_blk->sequence_nomask`): what the greedy traceback, the reevaluation and the
    /// identities read at an HSP's offsets, which can lie outside its context
    /// (`ContextQuery`).
    pub(super) encoded_query_concat_blastna: &'a [u8],
    pub(super) encoded_query_concat_nomask: &'a [u8],
    pub(super) query_lengths: &'a [usize],
    pub(super) query_concat_length: usize,
    pub(super) search_karlin_ref: &'a Vec<ContextKarlin>,
    pub(super) round_down_evalue_score: bool,
    /// The queries and their contexts (as `PrelimBatch` holds them), each context's effective
    /// search space, cutoff and hit-saving cutoff (`batch_setup::batch_cutoffs`: the setup's; in
    /// an `-ungapped` search, which has no traceback, the preliminary stage's, `PrelimBatch`'s; in
    /// `special_subject`'s copy the subject's own cutoffs), and the LOSAT_TIMING counters (`None`
    /// when timing is off).
    pub(super) queries: &'a [FastaRecord],
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) query_eff_searchsp: &'a Vec<i64>,
    pub(super) cutoff_scores: &'a Vec<i32>,
    pub(super) hit_saving_cutoff_scores: &'a Vec<i32>,
    pub(super) timing_ref: Option<&'a BlastnTiming>,
    /// The alignment: the score matrix, the reward and penalty, the gap costs, the X-dropoff of the
    /// traceback and the aligner (dynamic programming or greedy).
    pub(super) score_matrix: &'a BlastnaMatrix,
    pub(super) reward: i32,
    pub(super) penalty: i32,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    pub(super) x_drop_final: i32,
    pub(super) use_dp: bool,
    /// The hit-saving options: the minimum percent identity, the interval tree's diagonal
    /// separation and the e-value threshold.
    pub(super) percent_identity: f64,
    pub(super) min_diag_separation: i32,
    pub(super) evalue_threshold: f64,
    /// The schedule on the pool: whether the speculative batches may run, the rule of
    /// LOSAT_X_GREEDYSPEC (kept as its `x_reach!` key), the look-ahead window (`None` when the
    /// helpers do not run), the threads and the pool, and the query group tasks
    /// (LOSAT_X_QGROUPPAR).
    pub(super) speculative_traceback: bool,
    pub(super) x_greedy_spec: bool,
    pub(super) x_ahead_window: Option<usize>,
    pub(super) num_threads: usize,
    pub(super) parallel_pool: &'a SearchPool<'a>,
    pub(super) x_qgroup_batch: bool,
    /// The progress bar (one step per subject traced) and the diagnostics: BLEMIR_DEBUG,
    /// LOSAT_DEBUG_BLASTN and the trace lines (`blastn_trace`).
    pub(super) progress_bar: &'a Option<ProgressBar>,
    pub(super) debug_mode: bool,
    pub(super) blastn_debug: bool,
    pub(super) blastn_trace_enabled: bool,
}

/// The traceback of one subject, `search_subject`'s traceback stage: the early return of a subject
/// shorter than a word and the debug line of its cutoffs (each stage has its own: the preliminary
/// stage reads a database sequence's masks before them); the HSPs that the preliminary stage kept
/// (`kept`, copied into the list this function owns); if there are any, the subject's letters in
/// blastna (`traceback_subject_blastna`) and the values the traceback reads (`TracebackCtx`); then
/// an ungapped search's hits (`ungapped_traceback_hits`), or the interval tree set to the
/// subject's bounds and the traceback of the kept HSPs with the purges and the final list
/// (`traceback_list::traceback_subject_hits`), its debug line, and the lengths of its edit scripts
/// and alignments and its time (LOSAT_TIMING). The hits go to `subject_hits`; the progress bar
/// moves on by one. `subject_ctx` is what `search_subject` read of the call after its early
/// return without a valid context.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1526,1632-1633,1644,1691-1695,1701-1706
/// (one subject of the loop of `BLAST_ComputeTraceback_MT`, :1521-1713)
/// The progress bar has no NCBI counterpart.
pub(super) fn traceback_subject<const UNGAPPED: bool>(
    batch: &TraceBatch<'_>,
    subject_ctx: &SubjectCtx<'_>,
    subject: BlastnSubjectRef<'_>,
    kept: &Vec<PrelimHit>,
    gap_scratch: &mut GapAlignScratch,
    subject_scratch: &mut SubjectScratch,
    subject_hits: &mut Option<SubjectHits>,
) {
    let &TraceBatch {
        shortest_subject,
        encoded_queries_blastna,
        encoded_queries_nomask,
        encoded_query_concat_blastna,
        encoded_query_concat_nomask,
        query_lengths,
        query_concat_length,
        search_karlin_ref,
        round_down_evalue_score,
        score_matrix,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop_final,
        use_dp,
        percent_identity,
        min_diag_separation,
        evalue_threshold,
        speculative_traceback,
        x_greedy_spec,
        x_ahead_window,
        num_threads,
        parallel_pool,
        x_qgroup_batch,
        progress_bar,
        debug_mode,
        blastn_debug,
        blastn_trace_enabled,
        queries,
        query_contexts,
        query_eff_searchsp,
        cutoff_scores,
        hit_saving_cutoff_scores,
        timing_ref,
    } = batch;
    let &SubjectCtx {
        s_idx,
        s_id,
        s_len_full,
    } = subject_ctx;
    let timing_enabled = timing_ref.is_some();
    // EXPERIMENT (LOSAT_X_BNTBRANGE): set when a traceback alignment met a fence.
    let x_fenced = subject.x_fenced;
    // A subject that the preliminary stage did not search kept nothing; returning here keeps its
    // debug line off, as before.
    if s_len_full < shortest_subject {
        return;
    }

    // BLASTN debug: Log cutoff scores for this query-subject pair
    if blastn_debug {
        subject_ctx.debug_cutoffs(cutoff_scores);
    }

    let mut combined_prelim_hits: Vec<PrelimHit> = Vec::new();
    // The traceback of the HSPs that the preliminary stage kept (`search_subjects`).
    combined_prelim_hits.clone_from(kept);
    if combined_prelim_hits.is_empty() {
        return;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:358-373
    let traceback_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };

    let hits = &mut subject_scratch.hits;
    hits.clear();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1567-1573 (the
    // subject of the traceback, `traceback_subject_blastna`). The letters made for it are this
    // buffer's, which lives as long as `s_seq_blastna` is read.
    let mut s_seq_blastna_vec: Vec<u8> = Vec::new();
    let s_seq_blastna: &[u8] =
        traceback_subject_blastna(subject, parallel_pool, &mut s_seq_blastna_vec);

    // The values the traceback reads (`traceback_subject.rs`, `traceback_hsp.rs` and
    // `traceback_finish.rs`), for the query group tasks, the look-ahead helpers and the
    // speculative batch's jobs on the pool as for the owner thread.
    let traceback_ctx = TracebackCtx {
        timing_ref,
        timing_enabled,
        blastn_trace_enabled,
        s_idx,
        s_id,
        s_seq_blastna,
        x_fenced,
        encoded_queries_blastna,
        encoded_queries_nomask,
        encoded_query_concat_blastna,
        encoded_query_concat_nomask,
        query_contexts,
        query_lengths,
        cutoff_scores,
        hit_saving_cutoff_scores,
        query_eff_searchsp,
        search_karlin_ref,
        round_down_evalue_score,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop_final,
        score_matrix,
        use_dp,
        percent_identity,
        min_diag_separation,
        evalue_threshold,
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1526,1701-1706
    // The traceback of an ungapped (`-ungapped`) search aligns nothing and purges nothing
    // (`ungapped_traceback_hits`); a subject with hits gives them to the results.
    if UNGAPPED {
        let final_hits = ungapped_traceback_hits(&traceback_ctx, &combined_prelim_hits);
        if !final_hits.is_empty() {
            *subject_hits = Some(SubjectHits::Mixed(final_hits));
        }
        if let Some(bar) = progress_bar.as_ref() {
            bar.inc(1);
        }
        return;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:371-373
    // The tree is the scratch's (`x_tree`), given this subject's ranges.
    let interval_tree = &mut subject_scratch.x_tree;
    interval_tree.reset_with_bounds(
        0,
        (query_concat_length + 1) as i32,
        0,
        (s_len_full + 1) as i32,
    );
    let mut traceback_edit_script_lengths: Vec<u32> = Vec::new();
    let mut traceback_alignment_lengths: Vec<u32> = Vec::new();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375
    // What schedules the traceback loop of a list (`traceback_list::trace_hsp_list`) on the
    // pool: the speculative batches and the look-ahead helpers of the owner's loop.
    // No NCBI counterpart: scheduling only (NCBI runs the loop on one thread).
    let trace_schedule = traceback_list::TraceSchedule {
        speculative_traceback,
        x_greedy_spec,
        x_ahead_window,
        num_threads,
        parallel_pool,
    };

    // =================================================================
    // NCBI POST-GAPPED PROCESSING
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:633-692
    // =================================================================

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1632-1633,1644
    // The traceback of the subject's kept HSPs (`traceback_list::traceback_subject_hits`): the
    // traceback loop, the purges and the final list, per query on the pool (LOSAT_X_QGROUPPAR)
    // or over the subject's list on this thread with the subject's scratch.
    let prelim_hits = &mut combined_prelim_hits;
    let (final_hits, [hits_step1, hits_step2, hits_step3, hits_step4, final_tree_hit_count]) =
        traceback_list::traceback_subject_hits(
            &traceback_ctx,
            &trace_schedule,
            &traceback_list::TraceDrive {
                x_qgroup_batch,
                queries,
                query_concat_length,
                s_len_full,
            },
            prelim_hits,
            traceback_list::TraceListState {
                interval_tree,
                gap_scratch,
                greedy_align_scratch: &mut subject_scratch.greedy_align_scratch,
                x_helper_scratch: Some(&mut subject_scratch.x_ahead_scratch),
                x_spec_scratch: Some(&mut subject_scratch.x_spec_scratch),
                hits,
                traceback_edit_script_lengths: &mut traceback_edit_script_lengths,
                traceback_alignment_lengths: &mut traceback_alignment_lengths,
            },
        );
    prelim_hits.clear();

    // Phase 2 debug output
    if debug_mode || blastn_debug {
        let phase2_tree_removed = hits_step4 - final_tree_hit_count;
        let evalue_removed = final_tree_hit_count.saturating_sub(final_hits.len());
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:234-246
        // Report tree containment and the later E-value reap as separate
        // phases, preserving the same ordering as NCBI traceback cleanup.
        eprintln!(
                "[DEBUG] Phase 2 breakdown: step1(extract)={}, step2(purge1)={} (-{}), step3(reeval)={} (-{}), step4(purge2)={} (-{}), step6(tree)={} (-{}), evalue_reap={} (-{})",
                hits_step1,
                hits_step2, hits_step1 - hits_step2,
                hits_step3, hits_step2 - hits_step3,
                hits_step4, hits_step3 - hits_step4,
                final_tree_hit_count, phase2_tree_removed,
                final_hits.len(), evalue_removed
            );
    }

    if !final_hits.is_empty() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:633-692 (post-gapped processing is per-subject)
        *subject_hits = Some(final_hits);
    }
    if let Some(timing) = timing_ref {
        timing.record_traceback_lengths(
            std::mem::take(&mut traceback_edit_script_lengths),
            std::mem::take(&mut traceback_alignment_lengths),
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:358-373
    if let Some(timing) = timing_ref {
        if let Some(traceback_start) = traceback_start {
            let elapsed_ns = traceback_start.elapsed().as_nanos() as u64;
            timing
                .traceback_ns
                .fetch_add(elapsed_ns, std::sync::atomic::Ordering::Relaxed);
            timing
                .traceback_calls
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }
    if let Some(bar) = progress_bar.as_ref() {
        bar.inc(1);
    }
}

/// The subject's letters in blastna for its traceback: a database sequence's own, or the letters
/// made into `s_seq_blastna_vec`, which the caller keeps while it reads them.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1528, 1567-1573
/// The traceback fetches the uncompressed (blastna) subject of each subject with
/// preliminary HSPs; the preliminary stage reads only the packed one. The blastna letters
/// are made here, for this subject only, instead of being kept for every subject of the
/// search. A database sequence comes with its blastna letters, read with its
/// ambiguity codes (`database::traceback_listed`).
pub(super) fn traceback_subject_blastna<'v>(
    subject: BlastnSubjectRef<'v>,
    parallel_pool: &SearchPool<'_>,
    s_seq_blastna_vec: &'v mut Vec<u8>,
) -> &'v [u8] {
    match (subject.blastna, subject.store) {
        (Some(blastna), _) => blastna,
        // EXPERIMENT (LOSAT_X_BNSTORE): made from the packed subject, on the threads of the
        // pool (one thread decoding a long subject left the others waiting at every batch).
        (None, Some(store)) => {
            *s_seq_blastna_vec = Vec::new();
            x_blastna_into_pool(store, s_seq_blastna_vec, parallel_pool);
            s_seq_blastna_vec.as_slice()
        }
        (None, None) => {
            *s_seq_blastna_vec = encode_iupac_to_blastna(
                subject
                    .letters
                    .expect("a subject reaches the traceback with its letters"),
            );
            s_seq_blastna_vec.as_slice()
        }
    }
}

/// The hits of the traceback of an ungapped (`-ungapped`) search on one subject: it aligns
/// nothing and purges nothing. The HSPs the preliminary stage kept (re-evaluated, with their
/// e-values) become the subject's hits in the kept order, with the identities counted again on
/// the subject's BLASTNA letters (the ranges the traceback decoded include the HSPs').
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1526,1691-1696
pub(super) fn ungapped_traceback_hits(
    traceback_ctx: &TracebackCtx<'_>,
    combined_prelim_hits: &[PrelimHit],
) -> Vec<BlastnHsp> {
    let &TracebackCtx {
        s_idx,
        s_seq_blastna,
        encoded_query_concat_nomask,
        query_eff_searchsp,
        search_karlin_ref,
        round_down_evalue_score,
        query_lengths,
        ..
    } = traceback_ctx;

    let mut final_hits: Vec<BlastnHsp> = Vec::with_capacity(combined_prelim_hits.len());
    for prelim in combined_prelim_hits.iter() {
        let context_idx = prelim.context_idx as usize;
        let length = (prelim.prelim_qe - prelim.prelim_qs) as usize;
        // In the concatenated letters from NCBI's start, which an HSP that NCBI keeps outside its
        // context has outside it (`ungapped_subject_stage`).
        let matches = ungapped::count_identities_ungapped(
            encoded_query_concat_nomask,
            s_seq_blastna,
            (prelim.query_context_offset as isize + prelim.prelim_qs as isize) as usize,
            prelim.prelim_ss,
            length,
        );
        let eff_searchsp = query_eff_searchsp[context_idx];
        let (bit_score, eval) = calculate_blastn_context_statistics(
            prelim.prelim_score,
            &search_karlin_ref[context_idx].gapped,
            eff_searchsp,
            round_down_evalue_score,
        );
        let identity = if length > 0 {
            ((matches as f64 / length as f64) * 100.0).min(100.0)
        } else {
            0.0
        };
        let query_length = query_lengths[prelim.query_idx as usize];
        let (hit_q_start, hit_q_end, hit_s_start, hit_s_end) = adjust_blastn_offsets(
            prelim.prelim_qs,
            prelim.prelim_qe,
            prelim.prelim_ss,
            prelim.prelim_se,
            query_length,
            prelim.query_frame,
        );
        final_hits.push(BlastnHsp {
            identity,
            length,
            mismatch: length - matches,
            gapopen: 0,
            q_start: hit_q_start,
            q_end: hit_q_end,
            s_start: hit_s_start,
            s_end: hit_s_end,
            e_value: eval,
            bit_score,
            num_ident: matches,
            query_frame: prelim.query_frame,
            query_length,
            q_idx: prelim.query_idx,
            s_idx: s_idx as u32,
            raw_score: prelim.prelim_score,
            internal_q_offset_0: prelim.prelim_qs,
            internal_q_end_0: prelim.prelim_qe,
            internal_s_offset_0: prelim.prelim_ss,
            internal_s_end_0: prelim.prelim_se,
            internal_query_context_offset: prelim.query_context_offset,
            gap_info: None,
            num_positives: matches,
        });
    }
    final_hits
}
