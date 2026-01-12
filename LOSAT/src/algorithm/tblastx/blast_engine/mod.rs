//! TBLASTX search coordination and execution
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:438-1497
//!
//! This module contains the main `run()` and `run_with_neighbor_map()` functions
//! that coordinate the TBLASTX search process:
//! - Query and subject sequence reading
//! - Lookup table construction
//! - Parallel subject scanning with two-hit extension
//! - HSP reevaluation, linking, and culling
//! - Output generation
//!
//! # Module Structure
//!
//! - `mod.rs` - Shared structures and helpers (WorkerState, atomic ops, reevaluation)
//! - `run.rs` - Main run() function for backbone mode
//! - `run_neighbor_map.rs` - run_with_neighbor_map() for neighbor map mode

mod run_impl;
mod run_neighbor_map;

// Re-export the main run function
pub use run_impl::run;

// Re-export from run_neighbor_map for internal use
pub(crate) use run_neighbor_map::run_with_neighbor_map;

// Imports used by submodules
pub(crate) use anyhow::{Context, Result};
pub(crate) use bio::io::fasta;
pub(crate) use indicatif::{ProgressBar, ProgressStyle};
pub(crate) use rayon::prelude::*;
pub(crate) use std::collections::HashSet;
pub(crate) use std::sync::atomic::{AtomicI32, AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
pub(crate) use std::sync::mpsc::channel;
pub(crate) use std::time::{Duration, Instant};

pub(crate) use crate::common::{write_output_ncbi_order, Hit};
pub(crate) use crate::config::ScoringMatrix;
pub(crate) use crate::stats::{lookup_protein_params_ungapped, KarlinParams};
pub(crate) use crate::utils::genetic_code::GeneticCode;
pub(crate) use crate::utils::seg::SegMasker;

pub(crate) use super::args::TblastxArgs;
pub(crate) use super::chaining::UngappedHit;
pub(crate) use super::constants::{GAP_TRIGGER_BIT_SCORE, X_DROP_UNGAPPED_BITS};
pub(crate) use super::ncbi_cutoffs::{
    compute_eff_lengths_subject_mode_tblastx, cutoff_score_for_update_tblastx,
    cutoff_score_max_for_tblastx, gap_trigger_raw_score, x_drop_raw_score, BLAST_GAP_DECAY_RATE,
};
pub(crate) use super::diagnostics::{
    diagnostics_enabled, print_summary as print_diagnostics_summary, DiagnosticCounters,
};
pub(crate) use super::extension::{convert_coords, extend_hit_two_hit};
pub(crate) use super::lookup::{build_ncbi_lookup, encode_kmer, NeighborLookup, QueryContext};
pub(crate) use super::reevaluate::{
    get_num_identities_and_positives_ungapped, hsp_test,
    reevaluate_ungapped_hit_ncbi_translated,
};
pub(crate) use super::sum_stats_linking::{
    apply_sum_stats_even_gap_linking, compute_avg_query_length_ncbi,
    find_smallest_lambda_params, LinkingParams,
};
pub(crate) use super::hsp_culling;
pub(crate) use super::translation::{generate_frames, QueryFrame};
pub(crate) use super::tracing::{
    trace_final_hit_if_match, trace_hsp_target, trace_match_target, trace_ungapped_hit_if_match,
};
pub(crate) use crate::stats::karlin::bit_score as calc_bit_score;

// Import OffsetPair from scan submodule
pub(crate) use super::scan::OffsetPair;

// Import DiagStruct from blast_extend module (NCBI blast_extend.c equivalent)
pub(crate) use super::blast_extend::DiagStruct;

// Import InitHSP and related functions from blast_gapalign module (NCBI blast_gapalign.c equivalent)
pub(crate) use super::blast_gapalign::{get_ungapped_hsp_list, trace_init_hsp_if_match, InitHSP};

// Import subject scanning functions from blast_aascan module (NCBI blast_aascan.c equivalent)
pub(crate) use super::blast_aascan::s_blast_aa_scan_subject;

// ---------------------------------------------------------------------------
// Shared structures
// ---------------------------------------------------------------------------

pub(crate) struct WorkerState {
    pub tx: std::sync::mpsc::Sender<Vec<Hit>>,
    pub offset_pairs: Vec<OffsetPair>,
    pub diag_array: Vec<DiagStruct>,
    pub diag_offset: i32,
}

// ---------------------------------------------------------------------------
// Atomic helpers
// ---------------------------------------------------------------------------

#[inline]
pub(crate) fn atomic_max_usize(dst: &AtomicUsize, val: usize) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val > cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

#[inline]
pub(crate) fn atomic_min_i32(dst: &AtomicI32, val: i32) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val < cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

#[inline]
pub(crate) fn atomic_max_i32(dst: &AtomicI32, val: i32) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val > cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

// ---------------------------------------------------------------------------
// HSP reevaluation
// ---------------------------------------------------------------------------

/// NCBI Blast_HSPListReevaluateUngapped equivalent
/// Reference: blast_hits.c:2609-2737, blast_engine.c:1492-1497
///
/// NCBI code flow:
/// 1. For each HSP in the list
/// 2. Get context and compute query_start (context start position)
/// 3. Call Blast_HSPReevaluateWithAmbiguitiesUngapped
/// 4. Call Blast_HSPGetNumIdentitiesAndPositives and Blast_HSPTest
/// 5. Delete HSPs that fail tests
/// 6. Sort by score (scores may have changed)
///
/// NCBI reference (verbatim, blast_hits.c:2694-2706):
/// ```c
/// context = hsp->context;
/// query_start = query_blk->sequence + query_info->contexts[context].query_offset;
/// if (kTranslateSubject)
///     subject_start = Blast_HSPGetTargetTranslation(target_t, hsp, NULL);
/// if (kNucleotideSubject) {
///     delete_hsp = Blast_HSPReevaluateWithAmbiguitiesUngapped(hsp, query_start,
///         subject_start, word_params, sbp, kTranslateSubject);
/// }
/// ```
///
/// This function performs batch reevaluation on all HSPs in the list,
/// updating scores and trimming boundaries. HSPs that fail reevaluation
/// or subsequent tests are removed from the list.
pub(crate) fn reevaluate_ungapped_hsp_list(
    ungapped_hits: Vec<UngappedHit>,
    contexts: &[QueryContext],
    s_frames: &[QueryFrame],
    cutoff_scores: &[Vec<i32>],
    args: &TblastxArgs,
    timing_enabled: bool,
    reeval_ns: &AtomicU64,
    reeval_calls: &AtomicU64,
    is_long_sequence: bool,
    stats_hsp_filtered_by_reeval: &mut usize,
) -> Vec<UngappedHit> {
    let mut kept_hits = Vec::new();

    for mut hit in ungapped_hits {
        let ctx = &contexts[hit.ctx_idx];
        let s_frame = &s_frames[hit.s_f_idx];
        let cutoff = cutoff_scores[hit.ctx_idx][hit.s_f_idx];

        // NCBI: query_start = query_blk->sequence + query_info->contexts[context].query_offset;
        // query->sequence points past the leading NULLB sentinel.
        // Reference: blast_query_info.c:311-315, blast_util.c:112-116.
        let query_full = &ctx.aa_seq;
        let query = &query_full[1..query_full.len() - 1];
        let subject_full = &s_frame.aa_seq;
        let subject = &subject_full[1..subject_full.len() - 1];

        // NCBI: query = query_start + hsp->query.offset (0-based offsets)
        let qs = hit.q_aa_start;
        let ss = hit.s_aa_start;
        let len_u = hit.q_aa_end.saturating_sub(hit.q_aa_start);

        if len_u == 0 {
            continue;
        }

        // NCBI: Blast_HSPReevaluateWithAmbiguitiesUngapped (blast_hits.c:2702-2705)
        // Reference: blast_hits.c:675-733
        // NCBI reference (verbatim, blast_hits.c:688):
        //   Int4 cutoff_score = word_params->cutoffs[hsp->context].cutoff_score;
        //   return s_UpdateReevaluatedHSPUngapped(hsp, cutoff_score, score, ...);
        // If score < cutoff_score, the HSP is deleted (s_UpdateReevaluatedHSPUngapped returns FALSE)
        let t0 = if timing_enabled { Some(std::time::Instant::now()) } else { None };
        let (new_qs, new_ss, new_len, new_score) = if let Some(result) =
            reevaluate_ungapped_hit_ncbi_translated(query, subject, qs, ss, len_u, cutoff)
        {
            result
        } else {
            // Deleted by NCBI reevaluation (score < cutoff_score)
            // NCBI reference: blast_hits.c:730-732 (s_UpdateReevaluatedHSPUngapped returns FALSE)
            if is_long_sequence {
                *stats_hsp_filtered_by_reeval += 1;
            }
            if let Some(t) = t0 {
                let elapsed = t.elapsed();
                reeval_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
                reeval_calls.fetch_add(1, AtomicOrdering::Relaxed);
            }
            continue;
        };
        if let Some(t) = t0 {
            let elapsed = t.elapsed();
            reeval_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
            reeval_calls.fetch_add(1, AtomicOrdering::Relaxed);
        }

        // NCBI: Blast_HSPGetNumIdentitiesAndPositives and Blast_HSPTest
        // Reference: blast_hits.c:2708-2720
        // NCBI uses sequence_nomask that is also offset past the leading NULLB.
        // Reference: blast_filter.c:1381-1382, blast_util.c:112-116.
        let query_nomask_full = ctx.aa_seq_nomask.as_ref().unwrap_or(&ctx.aa_seq);
        let query_nomask = &query_nomask_full[1..query_nomask_full.len() - 1];
        let num_ident = if let Some((num_ident, _align_length)) =
            get_num_identities_and_positives_ungapped(
                query_nomask,
                subject,
                new_qs,
                new_ss,
                new_len,
            ) {
            num_ident
        } else {
            continue;
        };

        // NCBI: Blast_HSPTest (blast_hits.c:2719)
        let percent_identity = args.percent_identity;
        let min_hit_length = args.min_hit_length;
        if hsp_test(num_ident, new_len, percent_identity, min_hit_length) {
            continue;
        }

        // Update hit with reevaluated coordinates and score (0-based offsets)
        let new_qe = new_qs + new_len;
        let new_se = new_ss + new_len;
        hit.q_aa_start = new_qs;
        hit.q_aa_end = new_qe;
        hit.s_aa_start = new_ss;
        hit.s_aa_end = new_se;
        hit.raw_score = new_score;
        hit.num_ident = num_ident;

        trace_ungapped_hit_if_match("after_reevaluate", &hit);
        kept_hits.push(hit);
    }

    // NCBI: Sort the HSP array by score (scores may have changed!)
    // Reference: blast_hits.c:2734
    kept_hits.sort_by(|a, b| b.raw_score.cmp(&a.raw_score));

    kept_hits
}
