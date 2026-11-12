//! The diagnostics at the end of a TBLASTX subject's search, after its frame loop: the
//! HSP-saving statistics of a long subject (`print_hsp_saving_stats`), the chain roles and counts
//! of LOSAT_DEBUG_OUTPUT_FILTER (`count_chain_roles`, `print_output_filter`) and the mask trace of
//! LOSAT_TRACE_HSP_MASKS (`trace_hsp_masks`).
//!
//! No NCBI counterpart: they print what the search computed and change none of it. NCBI prints
//! output of this kind only in its debug builds (`-verbose`,
//! c++/src/algo/blast/blastinput/blast_args.cpp:3265-3274).
//!
//! This module is a child of `run_impl`, as subject_link.rs is. Its functions are statements of
//! `process_subject` moved here in their order (`print_hsp_saving_stats` is called from it, the
//! others from subject_report.rs); what they read is a struct of references and `Copy` values.

use super::{
    trace_hsp_target, trace_match_target, AtomicOrdering, DiagnosticCounters, QueryContext,
    QueryFrame, TblastxHsp, TbxScore, UngappedHit,
};

/// The HSP-saving statistics of a subject, which its frame loop and its re-evaluation counted.
pub(super) struct HspSavingStats<'a> {
    /// The subject is longer than 600,000 nt and LOSAT_DEBUG_HSP_SAVING or LOSAT_DIAGNOSTICS is
    /// set.
    pub(super) collect_hsp_saving_stats: bool,
    pub(super) stats_hsp_saved: usize,
    pub(super) stats_hsp_filtered_by_cutoff: usize,
    pub(super) stats_hsp_filtered_by_reeval: usize,
    pub(super) stats_hsp_filtered_by_hsp_test: usize,
    /// The scores of the saved HSPs (sorted when they are printed).
    pub(super) stats_score_distribution: &'a mut [i32],
    pub(super) subject_len_nucl: i64,
    pub(super) cutoff_score_min: i32,
}

// No NCBI counterpart: debugging output and a counter of LOSAT_DIAGNOSTICS (see the note in the
// body).
/// The `[DEBUG HSP_STATS]` lines of a long subject's HSPs after the re-evaluation
/// (`ungapped_hits`), its `hsps_before_chain` count (LOSAT_DIAGNOSTICS), then its
/// `[DEBUG HSP_SAVING]` lines.
pub(super) fn print_hsp_saving_stats(
    stats: HspSavingStats<'_>,
    ungapped_hits: &[UngappedHit],
    diag_enabled: bool,
    diagnostics: &DiagnosticCounters,
) {
    let HspSavingStats {
        collect_hsp_saving_stats,
        stats_hsp_saved,
        stats_hsp_filtered_by_cutoff,
        stats_hsp_filtered_by_reeval,
        stats_hsp_filtered_by_hsp_test,
        stats_score_distribution,
        subject_len_nucl,
        cutoff_score_min,
    } = stats;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
    // ```c
    // #if _BLAST_DEBUG
    // arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
    //                  true);
    // #endif /* _BLAST_DEBUG */
    // ```
    // DEBUG: Print HSP statistics for long sequences only when requested.
    if collect_hsp_saving_stats {
        eprintln!(
            "[DEBUG HSP_STATS] After reevaluation: {} HSPs",
            ungapped_hits.len()
        );
        eprintln!("[DEBUG HSP_STATS] Saved by cutoff: {}", stats_hsp_saved);
        eprintln!(
            "[DEBUG HSP_STATS] Filtered by cutoff: {}",
            stats_hsp_filtered_by_cutoff
        );
        eprintln!(
            "[DEBUG HSP_STATS] Filtered by reeval: {}",
            stats_hsp_filtered_by_reeval
        );
        if !stats_score_distribution.is_empty() {
            let min_score = stats_score_distribution.iter().min().unwrap();
            let max_score = stats_score_distribution.iter().max().unwrap();
            let avg_score = stats_score_distribution.iter().sum::<i32>() as f64
                / stats_score_distribution.len() as f64;
            eprintln!(
                "[DEBUG HSP_STATS] Score range: {} - {} (avg: {:.2})",
                min_score, max_score, avg_score
            );
        }
    }
    if diag_enabled {
        diagnostics
            .base
            .hsps_before_chain
            .fetch_add(ungapped_hits.len(), AtomicOrdering::Relaxed);
    }
    // Output HSP saving statistics for long sequences
    if collect_hsp_saving_stats
        && (stats_hsp_saved > 0
            || stats_hsp_filtered_by_cutoff > 0
            || stats_hsp_filtered_by_reeval > 0
            || stats_hsp_filtered_by_hsp_test > 0)
    {
        let total_attempted = stats_hsp_saved
            + stats_hsp_filtered_by_cutoff
            + stats_hsp_filtered_by_reeval
            + stats_hsp_filtered_by_hsp_test;
        eprintln!(
            "[DEBUG HSP_SAVING] subject_len_nucl={}, cutoff={}",
            subject_len_nucl, cutoff_score_min
        );
        eprintln!("[DEBUG HSP_SAVING] total_attempted={}, saved={}, filtered_by_cutoff={}, filtered_by_reeval={}, filtered_by_hsp_test={}",
                total_attempted, stats_hsp_saved, stats_hsp_filtered_by_cutoff, stats_hsp_filtered_by_reeval, stats_hsp_filtered_by_hsp_test);
        if !stats_score_distribution.is_empty() {
            stats_score_distribution.sort();
            let min_score = stats_score_distribution[0];
            let max_score = stats_score_distribution[stats_score_distribution.len() - 1];
            let median_score = stats_score_distribution[stats_score_distribution.len() / 2];
            let low_score_count = stats_score_distribution.iter().filter(|&&s| s < 30).count();
            eprintln!("[DEBUG HSP_SAVING] score_distribution: min={}, max={}, median={}, low_score(<30)={}/{} ({:.2}%)",
                    min_score, max_score, median_score, low_score_count, stats_score_distribution.len(),
                    if !stats_score_distribution.is_empty() { (low_score_count as f64 / stats_score_distribution.len() as f64) * 100.0 } else { 0.0 });
        }
    }
}

/// The roles of a subject's linked HSPs for LOSAT_DEBUG_OUTPUT_FILTER (`count_chain_roles`).
pub(super) struct ChainRoles {
    total_linked: usize,
    stats_single_hsps: usize,
    stats_chain_heads: usize,
    stats_chain_members: usize,
}

// No NCBI counterpart: debugging output (see the note in the body).
/// The number of linked HSPs and, with `debug_output_filter`, how many are single, chain heads and
/// chain members.
pub(super) fn count_chain_roles(linked: &[UngappedHit], debug_output_filter: bool) -> ChainRoles {
    let total_linked = linked.len();
    let mut stats_single_hsps = 0usize;
    let mut stats_chain_heads = 0usize;
    let mut stats_chain_members = 0usize;
    if debug_output_filter {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
        // ```c
        // #if _BLAST_DEBUG
        // arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
        //                  true);
        // #endif /* _BLAST_DEBUG */
        // ```
        // DEBUG: Collect statistics before filtering only when requested.
        for h in linked {
            if h.linked_set && !h.start_of_chain {
                stats_chain_members += 1;
            } else if !h.linked_set {
                stats_single_hsps += 1;
            } else if h.start_of_chain {
                stats_chain_heads += 1;
            }
        }
    }
    ChainRoles {
        total_linked,
        stats_single_hsps,
        stats_chain_heads,
        stats_chain_members,
    }
}

// No NCBI counterpart: debugging output (see the note in the body).
/// The `[DEBUG OUTPUT_FILTER]` lines of a subject: its chain roles, the HSPs the E-value reap
/// removed (`filtered_by_evalue`) and its report hits (`final_hits`).
pub(super) fn print_output_filter(
    roles: &ChainRoles,
    evalue_threshold: f64,
    filtered_by_evalue: usize,
    final_hits: &[TblastxHsp],
) {
    let &ChainRoles {
        total_linked,
        stats_single_hsps,
        stats_chain_heads,
        stats_chain_members,
    } = roles;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
    // ```c
    // #if _BLAST_DEBUG
    // arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
    //                  true);
    // #endif /* _BLAST_DEBUG */
    // ```
    // DEBUG: Print output filtering statistics.
    // NCBI reference: link_hsps.c:1018-1020 - continue is NOT an output filter
    // Chain members are included in output via link pointer traversal.
    eprintln!("[DEBUG OUTPUT_FILTER] Total linked HSPs: {}", total_linked);
    eprintln!(
        "[DEBUG OUTPUT_FILTER] Single HSPs (linked_set=false): {}",
        stats_single_hsps
    );
    eprintln!(
        "[DEBUG OUTPUT_FILTER] Chain heads (linked_set=true, start_of_chain=true): {}",
        stats_chain_heads
    );
    eprintln!("[DEBUG OUTPUT_FILTER] Chain members (linked_set=true, start_of_chain=false): {} (included in output)", stats_chain_members);
    eprintln!(
        "[DEBUG OUTPUT_FILTER] Filtered by E-value (threshold={}): {}",
        evalue_threshold, filtered_by_evalue
    );
    eprintln!(
        "[DEBUG OUTPUT_FILTER] Expected after E-value filter: {} (all HSPs - E-value filtered)",
        total_linked - filtered_by_evalue
    );
    eprintln!(
        "[DEBUG OUTPUT_FILTER] Final hits after filtering: {}",
        final_hits.len()
    );
}

/// The traced HSP (`h`) for `trace_hsp_masks`: its query context and subject frames, its length
/// and frame offsets, and its report coordinates.
pub(super) struct MaskTrace<'a> {
    pub(super) h: &'a UngappedHit,
    pub(super) ctx: &'a QueryContext,
    pub(super) s_frame: &'a QueryFrame,
    pub(super) s_score_frame: &'a QueryFrame,
    pub(super) len: usize,
    pub(super) q0: usize,
    pub(super) s0: usize,
    pub(super) q_start: usize,
    pub(super) q_end: usize,
    pub(super) s_start: usize,
    pub(super) s_end: usize,
}

// No NCBI counterpart: a trace on stderr (see the note in the body).
/// With LOSAT_TRACE_HSP_MASKS (the caller reads it), the `[TRACE_HSP_MASKS]` line of a final HSP
/// that LOSAT_TRACE_HSP names: its score on the masked and on the unmasked query, and the masked
/// runs.
pub(super) fn trace_hsp_masks<S: TbxScore>(sc: S, trace: &MaskTrace<'_>) {
    let &MaskTrace {
        h,
        ctx,
        s_frame,
        s_score_frame,
        len,
        q0,
        s0,
        q_start,
        q_end,
        s_start,
        s_end,
    } = trace;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1379-1404
    // ```c
    // query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
    // query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
    // Blast_MaskTheResidues(buffer, query_length, kIsNucl, ...);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:699-700
    // ```c
    // sum += matrix[*query & kResidueMask][*subject];
    // query++;
    // ```
    // Diagnostic only: compare the traced HSP's masked working-query
    // score against the preserved unmasked query copy used for identity
    // reporting. Normal runtime behavior is unchanged unless both
    // LOSAT_TRACE_HSP and LOSAT_TRACE_HSP_MASKS are set.
    if let Some(target) = trace_hsp_target() {
        if trace_match_target(target, q_start, q_end, s_start, s_end) {
            // NCBI uses query_nomask = query_blk->sequence_nomask + query_offset,
            // with sequence_nomask pointing past the leading NULLB.
            // Reference: blast_filter.c:1381-1382, blast_util.c:112-116.
            let q_seq_nomask_full: &[u8] = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
            let q_seq_nomask = &q_seq_nomask_full[1..q_seq_nomask_full.len() - 1];
            let s_seq = &s_frame.aa_seq[1..s_frame.aa_seq.len() - 1];
            let q_seq_masked = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
            let s_seq_scoring = &s_score_frame.aa_seq[1..s_score_frame.aa_seq.len() - 1];
            let mut masked_residues = 0usize;
            let mut masked_score = 0i32;
            let mut unmasked_score = 0i32;
            let mut report_unmasked_score = 0i32;
            let mut mask_runs = Vec::new();
            let mut current_mask_start: Option<usize> = None;
            for rel in 0..len {
                let q_pos = q0 + rel;
                let s_pos = s0 + rel;
                if q_pos >= q_seq_masked.len()
                    || q_pos >= q_seq_nomask.len()
                    || s_pos >= s_seq.len()
                    || s_pos >= s_seq_scoring.len()
                {
                    break;
                }
                let q_masked = q_seq_masked[q_pos];
                let q_unmasked = q_seq_nomask[q_pos];
                let subject_scoring = s_seq_scoring[s_pos];
                let subject_report = s_seq[s_pos];
                if q_masked != q_unmasked {
                    masked_residues += 1;
                    if current_mask_start.is_none() {
                        current_mask_start = Some(rel);
                    }
                } else if let Some(start) = current_mask_start.take() {
                    mask_runs.push(format!("{}..{}", start, rel));
                }
                masked_score += sc.score(q_masked, subject_scoring);
                unmasked_score += sc.score(q_unmasked, subject_scoring);
                report_unmasked_score += sc.score(q_unmasked, subject_report);
            }
            if let Some(start) = current_mask_start.take() {
                mask_runs.push(format!("{}..{}", start, len));
            }
            eprintln!(
                "[TRACE_HSP_MASKS] q={}-{} s={}-{} ctx_idx={} q_frame={} s_frame={} len={} raw_score={} masked_score={} unmasked_score={} report_unmasked_score={} masked_residues={} mask_runs={}",
                q_start,
                q_end,
                s_start,
                s_end,
                h.ctx_idx,
                ctx.frame,
                s_frame.frame,
                len,
                h.raw_score,
                masked_score,
                unmasked_score,
                report_unmasked_score,
                masked_residues,
                if mask_runs.is_empty() {
                    "none".to_string()
                } else {
                    mask_runs.join(",")
                }
            );
        }
    }
}
