//! Common utilities for HSP chaining and filtering
//!
//! This module provides shared utilities for HSP chaining used by both BLASTN and TBLASTX.
//! The actual chaining algorithms remain algorithm-specific due to differences in:
//! - Coordinate systems (nucleotide vs amino acid)
//! - Data structures (Hit vs ExtendedHit)
//! - Frame handling (TBLASTX-specific)
//!
//! However, common utilities like overlap calculation and filtering can be shared.

use crate::common::Hit;

/// Calculate overlap between two intervals
///
/// Returns the number of positions that overlap between the two intervals.
/// If there's no overlap, returns 0.
pub fn calculate_overlap(start1: usize, end1: usize, start2: usize, end2: usize) -> usize {
    let overlap_start = start1.max(start2);
    let overlap_end = end1.min(end2);

    if overlap_start <= overlap_end {
        overlap_end - overlap_start + 1
    } else {
        0
    }
}

/// Filter overlapping HSPs, keeping higher-scoring ones
///
/// This implements a simple overlap-based filtering where HSPs that overlap
/// significantly with higher-scoring HSPs are removed.
///
/// # Arguments
/// * `hits` - Vector of hits to filter (will be sorted by bit score)
/// * `overlap_threshold` - Fraction of overlap required for filtering (0.0 to 1.0)
///
/// # Returns
/// Filtered vector of hits with redundant overlaps removed
pub fn filter_overlapping_hsps(hits: Vec<Hit>, overlap_threshold: f64) -> Vec<Hit> {
    if hits.is_empty() || overlap_threshold <= 0.0 {
        return hits;
    }

    // Sort by bit score descending
    let mut sorted = hits;
    sorted.sort_by(|a, b| {
        b.bit_score
            .partial_cmp(&a.bit_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut kept: Vec<Hit> = Vec::new();

    for hit in sorted {
        let dominated = kept.iter().any(|kept_hit| {
            // Must be same query-subject pair
            // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
            // ```c
            // typedef struct BlastHSPList {
            //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
            //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
            //                       Set to 0 if not applicable */
            // } BlastHSPList;
            // ```
            if hit.q_idx != kept_hit.q_idx || hit.s_idx != kept_hit.s_idx {
                return false;
            }

            // Calculate overlap fractions
            let q_overlap = calculate_overlap(hit.q_start, hit.q_end, kept_hit.q_start, kept_hit.q_end);
            let s_overlap = calculate_overlap(hit.s_start, hit.s_end, kept_hit.s_start, kept_hit.s_end);

            let hit_q_len = (hit.q_end - hit.q_start + 1) as f64;
            let hit_s_len = (hit.s_end - hit.s_start + 1) as f64;

            let q_frac = q_overlap as f64 / hit_q_len;
            let s_frac = s_overlap as f64 / hit_s_len;

            q_frac > overlap_threshold && s_frac > overlap_threshold
        });

        if !dominated {
            kept.push(hit);
        }
    }

    kept
}

/// Check if two HSPs can be chained based on gap and diagonal drift
///
/// This is a common utility for checking if two HSPs are close enough to be chained.
/// The actual chaining logic remains algorithm-specific.
///
/// # Arguments
/// * `q_start1`, `q_end1` - Query coordinates of first HSP
/// * `s_start1`, `s_end1` - Subject coordinates of first HSP
/// * `q_start2`, `q_end2` - Query coordinates of second HSP
/// * `s_start2`, `s_end2` - Subject coordinates of second HSP
/// * `max_gap` - Maximum gap allowed between HSPs
/// * `max_diag_drift` - Maximum diagonal drift allowed
///
/// # Returns
/// `true` if the HSPs can be chained, `false` otherwise
pub fn can_chain_hsps(
    q_start1: usize,
    q_end1: usize,
    s_start1: usize,
    s_end1: usize,
    q_start2: usize,
    q_end2: usize,
    s_start2: usize,
    _s_end2: usize,
    max_gap: usize,
    max_diag_drift: isize,
) -> bool {
    let diag1 = s_start1 as isize - q_start1 as isize;
    let diag2 = s_start2 as isize - q_start2 as isize;
    let diag_drift = (diag2 - diag1).abs();

    if diag_drift > max_diag_drift {
        return false;
    }

    // Calculate gap or overlap
    let q_distance = q_start2 as isize - q_end1 as isize;
    let s_distance = s_start2 as isize - s_end1 as isize;

    // Allow some overlap (negative distance) but limit it
    let hit_q_len = (q_end1 - q_start1) as isize;
    let hit_q_len2 = (q_end2 - q_start2) as isize;
    let max_overlap = -(hit_q_len.min(hit_q_len2) / 2);

    let q_ok = q_distance >= max_overlap && q_distance <= max_gap as isize;
    let s_ok = s_distance >= max_overlap && s_distance <= max_gap as isize;

    q_ok && s_ok
}

