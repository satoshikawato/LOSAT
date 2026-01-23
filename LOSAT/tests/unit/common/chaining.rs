//! Unit tests for common/chaining.rs

use LOSAT::algorithm::common::chaining::{calculate_overlap, can_chain_hsps, filter_overlapping_hsps};
use LOSAT::common::Hit;
use super::super::helpers::make_hit;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
// ```c
// typedef struct BlastHSPList {
//    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
//    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
//                       Set to 0 if not applicable */
// } BlastHSPList;
// ```
const Q1_IDX: u32 = 0;
const Q2_IDX: u32 = 1;
const S1_IDX: u32 = 0;
const S2_IDX: u32 = 1;

#[test]
fn test_calculate_overlap_no_overlap() {
    // No overlap case
    assert_eq!(calculate_overlap(1, 100, 101, 200), 0);
    assert_eq!(calculate_overlap(101, 200, 1, 100), 0);
}

#[test]
fn test_calculate_overlap_full_overlap() {
    // Full overlap
    assert_eq!(calculate_overlap(1, 100, 1, 100), 100);
    assert_eq!(calculate_overlap(50, 150, 50, 150), 101);
}

#[test]
fn test_calculate_overlap_partial_overlap() {
    // Partial overlap
    assert_eq!(calculate_overlap(1, 100, 50, 150), 51); // Overlap: 50-100
    assert_eq!(calculate_overlap(50, 150, 1, 100), 51);
    assert_eq!(calculate_overlap(10, 50, 30, 70), 21); // Overlap: 30-50
}

#[test]
fn test_calculate_overlap_adjacent() {
    // Adjacent intervals (no overlap)
    assert_eq!(calculate_overlap(1, 100, 101, 200), 0);
}

#[test]
fn test_calculate_overlap_one_inside_other() {
    // One interval completely inside the other
    assert_eq!(calculate_overlap(1, 100, 10, 90), 81);
    assert_eq!(calculate_overlap(10, 90, 1, 100), 81);
}

#[test]
fn test_filter_overlapping_hsps_empty() {
    let hits = vec![];
    let filtered = filter_overlapping_hsps(hits, 0.5);
    assert_eq!(filtered.len(), 0);
}

#[test]
fn test_filter_overlapping_hsps_single_hit() {
    let hits = vec![make_hit(Q1_IDX, S1_IDX, 1, 100, 1, 100, 100.0)];
    let filtered = filter_overlapping_hsps(hits, 0.5);
    assert_eq!(filtered.len(), 1);
}

#[test]
fn test_filter_overlapping_hsps_no_overlap() {
    let hits = vec![
        make_hit(Q1_IDX, S1_IDX, 1, 100, 1, 100, 100.0),
        make_hit(Q1_IDX, S1_IDX, 200, 300, 200, 300, 90.0),
    ];
    let filtered = filter_overlapping_hsps(hits, 0.5);
    assert_eq!(filtered.len(), 2); // Both should be kept
}

#[test]
fn test_filter_overlapping_hsps_high_overlap() {
    let hits = vec![
        make_hit(Q1_IDX, S1_IDX, 1, 100, 1, 100, 100.0), // Higher score
        make_hit(Q1_IDX, S1_IDX, 10, 90, 10, 90, 80.0),  // Overlaps significantly
    ];
    let filtered = filter_overlapping_hsps(hits, 0.5);
    assert_eq!(filtered.len(), 1); // Lower scoring hit should be filtered
    assert_eq!(filtered[0].bit_score, 100.0);
}

#[test]
fn test_filter_overlapping_hsps_low_overlap() {
    let hits = vec![
        make_hit(Q1_IDX, S1_IDX, 1, 100, 1, 100, 100.0),
        make_hit(Q1_IDX, S1_IDX, 80, 120, 80, 120, 90.0), // Overlap: (80, 100) = 21 positions
    ];
    let filtered = filter_overlapping_hsps(hits, 0.5);
    // Overlap is 21/41 = 0.512 > 0.5, so lower-scoring hit should be filtered
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].bit_score, 100.0);
}

#[test]
fn test_filter_overlapping_hsps_different_subjects() {
    let hits = vec![
        make_hit(Q1_IDX, S1_IDX, 1, 100, 1, 100, 100.0),
        make_hit(Q1_IDX, S2_IDX, 1, 100, 1, 100, 90.0), // Different subject
    ];
    let filtered = filter_overlapping_hsps(hits, 0.5);
    // Different subjects should not filter each other
    assert_eq!(filtered.len(), 2);
}

#[test]
fn test_filter_overlapping_hsps_different_queries() {
    let hits = vec![
        make_hit(Q1_IDX, S1_IDX, 1, 100, 1, 100, 100.0),
        make_hit(Q2_IDX, S1_IDX, 1, 100, 1, 100, 90.0), // Different query
    ];
    let filtered = filter_overlapping_hsps(hits, 0.5);
    // Different queries should not filter each other
    assert_eq!(filtered.len(), 2);
}

#[test]
fn test_filter_overlapping_hsps_sorted_by_score() {
    let hits = vec![
        make_hit(Q1_IDX, S1_IDX, 1, 100, 1, 100, 50.0),  // Lower score
        make_hit(Q1_IDX, S1_IDX, 10, 90, 10, 90, 100.0), // Higher score, overlaps
    ];
    let filtered = filter_overlapping_hsps(hits, 0.5);
    // Should keep higher scoring hit
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].bit_score, 100.0);
}

#[test]
fn test_filter_overlapping_hsps_zero_threshold() {
    let hits = vec![
        make_hit(Q1_IDX, S1_IDX, 1, 100, 1, 100, 100.0),
        make_hit(Q1_IDX, S1_IDX, 10, 90, 10, 90, 80.0),
    ];
    let filtered = filter_overlapping_hsps(hits, 0.0);
    // Zero threshold should return all hits
    assert_eq!(filtered.len(), 2);
}

#[test]
fn test_can_chain_hsps_close_hsps() {
    // HSPs that are close and on same diagonal
    let q_start1 = 100;
    let q_end1 = 200;
    let s_start1 = 100;
    let s_end1 = 200;
    let q_start2 = 210;
    let q_end2 = 300;
    let s_start2 = 210;
    let s_end2 = 300;
    let max_gap = 20;
    let max_diag_drift = 5;

    assert!(can_chain_hsps(
        q_start1, q_end1, s_start1, s_end1,
        q_start2, q_end2, s_start2, s_end2,
        max_gap, max_diag_drift
    ));
}

#[test]
fn test_can_chain_hsps_too_far() {
    // HSPs that are too far apart
    let q_start1 = 100;
    let q_end1 = 200;
    let s_start1 = 100;
    let s_end1 = 200;
    let q_start2 = 500;
    let q_end2 = 600;
    let s_start2 = 500;
    let s_end2 = 600;
    let max_gap = 20;
    let max_diag_drift = 5;

    assert!(!can_chain_hsps(
        q_start1, q_end1, s_start1, s_end1,
        q_start2, q_end2, s_start2, s_end2,
        max_gap, max_diag_drift
    ));
}

#[test]
fn test_can_chain_hsps_large_diagonal_drift() {
    // HSPs with large diagonal drift
    let q_start1 = 100;
    let q_end1 = 200;
    let s_start1 = 100;
    let s_end1 = 200;
    let q_start2 = 210;
    let q_end2 = 300;
    let s_start2 = 250; // Large diagonal drift
    let s_end2 = 350;
    let max_gap = 20;
    let max_diag_drift = 5;

    assert!(!can_chain_hsps(
        q_start1, q_end1, s_start1, s_end1,
        q_start2, q_end2, s_start2, s_end2,
        max_gap, max_diag_drift
    ));
}

#[test]
fn test_can_chain_hsps_small_overlap() {
    // HSPs with small overlap (should be allowed)
    let q_start1 = 100;
    let q_end1 = 200;
    let s_start1 = 100;
    let s_end1 = 200;
    let q_start2 = 190; // Small overlap
    let q_end2 = 290;
    let s_start2 = 190;
    let s_end2 = 290;
    let max_gap = 20;
    let max_diag_drift = 5;

    // Small overlap should be allowed
    assert!(can_chain_hsps(
        q_start1, q_end1, s_start1, s_end1,
        q_start2, q_end2, s_start2, s_end2,
        max_gap, max_diag_drift
    ));
}

#[test]
fn test_can_chain_hsps_large_overlap() {
    // HSPs with large overlap (should not be allowed)
    let q_start1 = 100;
    let q_end1 = 200;
    let s_start1 = 100;
    let s_end1 = 200;
    let q_start2 = 110; // Large overlap (> 50%)
    let q_end2 = 190;
    let s_start2 = 110;
    let s_end2 = 190;
    let max_gap = 20;
    let max_diag_drift = 5;

    // Large overlap should not be allowed
    assert!(!can_chain_hsps(
        q_start1, q_end1, s_start1, s_end1,
        q_start2, q_end2, s_start2, s_end2,
        max_gap, max_diag_drift
    ));
}

