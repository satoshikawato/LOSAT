//! Unit tests for blastn/constants.rs

use LOSAT::algorithm::blastn::constants::*;

#[test]
fn test_x_drop_constants() {
    // X-drop values should match NCBI BLAST defaults
    assert_eq!(X_DROP_UNGAPPED, 20);
    assert_eq!(X_DROP_GAPPED_FINAL, 100);
}

#[test]
fn test_two_hit_window() {
    // TWO_HIT_WINDOW should be a reasonable value for two-hit extension
    assert!(TWO_HIT_WINDOW > 0);
    assert_eq!(TWO_HIT_WINDOW, 64);
}

#[test]
fn test_max_hits_per_kmer() {
    // MAX_HITS_PER_KMER should limit memory usage
    assert!(MAX_HITS_PER_KMER > 0);
    assert_eq!(MAX_HITS_PER_KMER, 200);
}

#[test]
fn test_min_ungapped_score_thresholds() {
    // Thresholds should be positive and reasonable
    assert!(MIN_UNGAPPED_SCORE_MEGABLAST > 0);
    assert!(MIN_UNGAPPED_SCORE_BLASTN > 0);
    
    // blastn threshold should be higher (more filtering needed for smaller word size)
    assert!(MIN_UNGAPPED_SCORE_BLASTN >= MIN_UNGAPPED_SCORE_MEGABLAST);
    
    assert_eq!(MIN_UNGAPPED_SCORE_MEGABLAST, 20);
    assert_eq!(MIN_UNGAPPED_SCORE_BLASTN, 50);
}

#[test]
fn test_max_direct_lookup_word_size() {
    // MAX_DIRECT_LOOKUP_WORD_SIZE should be reasonable to avoid excessive memory
    // 4^13 = 67M entries, which is manageable
    assert!(MAX_DIRECT_LOOKUP_WORD_SIZE <= 13);
    assert_eq!(MAX_DIRECT_LOOKUP_WORD_SIZE, 13);
}

#[test]
fn test_greedy_alignment_constants() {
    // GREEDY_MAX_COST should be positive
    assert!(GREEDY_MAX_COST > 0);
    assert_eq!(GREEDY_MAX_COST, 1000);
    
    // GREEDY_MAX_COST_FRACTION should be positive
    assert!(GREEDY_MAX_COST_FRACTION > 0);
    assert_eq!(GREEDY_MAX_COST_FRACTION, 2);
}

#[test]
fn test_invalid_offset() {
    // INVALID_OFFSET should be a negative value that's unlikely to occur naturally
    assert!(INVALID_OFFSET < 0);
    assert_eq!(INVALID_OFFSET, -2);
}

#[test]
fn test_invalid_diag() {
    // INVALID_DIAG should be a very large value
    assert!(INVALID_DIAG > 0);
    assert_eq!(INVALID_DIAG, 100_000_000);
}

