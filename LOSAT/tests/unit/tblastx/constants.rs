//! Unit tests for tblastx/constants.rs

use LOSAT::algorithm::tblastx::constants::*;

#[test]
fn test_x_drop_constants() {
    // X-drop values should match NCBI BLAST protein defaults
    assert_eq!(X_DROP_UNGAPPED, 11);
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
fn test_max_hits_per_kmer() {
    // MAX_HITS_PER_KMER should limit memory usage
    assert!(MAX_HITS_PER_KMER > 0);
    assert_eq!(MAX_HITS_PER_KMER, 200);
}

#[test]
fn test_stop_codon_encoding() {
    // STOP_CODON should be 25 (after Z which is 24)
    assert_eq!(STOP_CODON, 25);
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
    assert!(MIN_UNGAPPED_SCORE > 0);
    assert_eq!(MIN_UNGAPPED_SCORE, 22);
}

#[test]
fn test_chaining_constants() {
    // Chaining constants should be reasonable for amino acid coordinates
    assert!(MAX_GAP_AA > 0);
    assert_eq!(MAX_GAP_AA, 333);
    
    assert!(MAX_DIAG_DRIFT_AA > 0);
    assert_eq!(MAX_DIAG_DRIFT_AA, 33);
}

