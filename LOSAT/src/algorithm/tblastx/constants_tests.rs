//! Unit tests for tblastx/constants.rs

use crate::algorithm::tblastx::constants::*;

#[test]
fn test_x_drop_constants() {
    // X-drop values should match NCBI BLAST protein defaults
    assert_eq!(X_DROP_GAPPED_PRELIM, 15);
    assert_eq!(X_DROP_GAPPED_FINAL, 25);
}

#[test]
fn test_min_ungapped_score() {
    // MIN_UNGAPPED_SCORE should be positive
    // Set to 14 to capture more low-scoring hits
    const { assert!(MIN_UNGAPPED_SCORE > 0) };
    assert_eq!(MIN_UNGAPPED_SCORE, 14);
}

#[test]
fn test_sentinel_constants() {
    // NCBI BLAST sentinel (NULLB) is 0.
    assert_eq!(SENTINEL_BYTE, 0);
}
