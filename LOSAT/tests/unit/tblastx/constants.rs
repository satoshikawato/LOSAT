//! Unit tests for tblastx/constants.rs

use LOSAT::algorithm::tblastx::constants::*;

#[test]
fn test_x_drop_constants() {
    // X-drop values should match NCBI BLAST protein defaults
    // NCBI uses 7 bits and converts to raw score at runtime:
    //   ceil(7 * ln(2) / lambda) = ceil(15.27...) = 16  (BLOSUM62 ungapped lambda≈0.3176)
    assert_eq!(X_DROP_UNGAPPED, 16);
    assert_eq!(X_DROP_GAPPED_PRELIM, 15);
    assert_eq!(X_DROP_GAPPED_FINAL, 25);
}

#[test]
fn test_two_hit_window() {
    // TWO_HIT_WINDOW should be reasonable for protein searches
    assert!(TWO_HIT_WINDOW > 0);
    assert_eq!(TWO_HIT_WINDOW, 40);
}

#[test]
fn test_stop_codon_encoding() {
    // STOP_CODON should be 24 (index in NCBI matrix order: ARNDCQEGHILKMFPSTWYVBJZX*)
    assert_eq!(STOP_CODON, 24);
}

#[test]
fn test_gap_penalties() {
    // Gap penalties should match NCBI BLAST protein defaults
    assert_eq!(GAP_OPEN, -11);
    assert_eq!(GAP_EXTEND, -1);
}

#[test]
fn test_high_score_threshold() {
    // HIGH_SCORE_THRESHOLD should be positive
    assert!(HIGH_SCORE_THRESHOLD > 0);
    assert_eq!(HIGH_SCORE_THRESHOLD, 60);
}

#[test]
fn test_min_ungapped_score() {
    // MIN_UNGAPPED_SCORE should be positive
    // Set to 14 to capture more low-scoring hits
    assert!(MIN_UNGAPPED_SCORE > 0);
    assert_eq!(MIN_UNGAPPED_SCORE, 14);
}

#[test]
fn test_chaining_constants() {
    // Chaining constants should be reasonable for amino acid coordinates
    assert!(MAX_GAP_AA > 0);
    assert_eq!(MAX_GAP_AA, 333);
    
    assert!(MAX_DIAG_DRIFT_AA > 0);
    assert_eq!(MAX_DIAG_DRIFT_AA, 33);
}

#[test]
fn test_sentinel_constants() {
    // NCBI BLAST sentinel (NULLB) is 0.
    assert_eq!(SENTINEL_BYTE, 0);
    // Sentinel score must match NCBI BLOSUM62 defscore (-4).
    assert_eq!(SENTINEL_PENALTY, -4);
}
