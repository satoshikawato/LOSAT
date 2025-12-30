//! Unit tests for blastn/extension.rs

use LOSAT::algorithm::blastn::extension::extend_hit_ungapped;
// Note: make_nucleotide_sequence is not actually used in these tests

#[test]
fn test_extend_hit_ungapped_perfect_match() {
    let q_seq = b"ACGTACGT";
    let s_seq = b"ACGTACGT";
    
    // Start at position 2 (G)
    let (q_start, q_end, s_start, s_end, score) = extend_hit_ungapped(q_seq, s_seq, 2, 2, 1, -2, None);
    
    // Should extend to cover the entire matching region
    assert!(score > 0);
    assert!(q_start <= 2);
    assert!(q_end > 2);
    assert!(s_start <= 2);
    assert!(s_end > 2);
}

#[test]
fn test_extend_hit_ungapped_with_mismatches() {
    let q_seq = b"ACGTACGT";
    let s_seq = b"ACGTCCGT"; // One mismatch at position 4
    
    let (q_start, q_end, s_start, s_end, score) = extend_hit_ungapped(q_seq, s_seq, 2, 2, 1, -2, None);
    
    // Should still extend but score will be lower due to mismatch
    assert!(score > 0);
    assert!(q_start <= 2);
    assert!(q_end > 2);
}

#[test]
fn test_extend_hit_ungapped_x_drop_termination() {
    // Create sequences with many mismatches after initial match
    let q_seq = b"ACGTACGTACGTACGT";
    let s_seq = b"ACGTNNNNNNNNNNNN"; // Many mismatches after initial match
    
    let (q_start, q_end, s_start, s_end, score) = extend_hit_ungapped(q_seq, s_seq, 2, 2, 1, -2, None);
    
    // X-drop should terminate extension when score drops too much
    assert!(score >= 0);
    // Extension should be limited by X-drop
}

#[test]
fn test_extend_hit_ungapped_at_sequence_start() {
    let q_seq = b"ACGT";
    let s_seq = b"ACGT";
    
    // Start at position 0 (beginning of sequence)
    let (q_start, q_end, s_start, s_end, score) = extend_hit_ungapped(q_seq, s_seq, 0, 0, 1, -2, None);
    
    // Should extend rightward only (can't extend left from position 0)
    assert!(q_start == 0);
    assert!(q_end > 0);
    assert!(score > 0);
}

#[test]
fn test_extend_hit_ungapped_at_sequence_end() {
    let q_seq = b"ACGT";
    let s_seq = b"ACGT";
    
    // Start at last position
    let last_pos = q_seq.len() - 1;
    let (q_start, q_end, s_start, s_end, score) = extend_hit_ungapped(q_seq, s_seq, last_pos, last_pos, 1, -2, None);
    
    // Should extend leftward only (can't extend right beyond end)
    assert!(q_end <= q_seq.len());
    assert!(score > 0);
}

#[test]
fn test_extend_hit_ungapped_different_reward_penalty() {
    let q_seq = b"ACGTACGT";
    let s_seq = b"ACGTACGT";
    
    // Test with higher reward
    let (_, _, _, _, score1) = extend_hit_ungapped(q_seq, s_seq, 2, 2, 2, -3, None);
    
    // Test with lower reward
    let (_, _, _, _, score2) = extend_hit_ungapped(q_seq, s_seq, 2, 2, 1, -2, None);
    
    // Higher reward should give higher score for same match
    assert!(score1 > score2);
}

#[test]
fn test_extend_hit_ungapped_coordinates() {
    let q_seq = b"ACGTACGTACGT";
    let s_seq = b"ACGTACGTACGT";
    
    // Start at middle position
    let start_pos = 4;
    let (q_start, q_end, s_start, s_end, _score) = extend_hit_ungapped(q_seq, s_seq, start_pos, start_pos, 1, -2, None);
    
    // Coordinates should be valid
    assert!(q_start <= start_pos);
    assert!(q_end > start_pos);
    assert!(s_start <= start_pos);
    assert!(s_end > start_pos);
    assert!(q_end <= q_seq.len());
    assert!(s_end <= s_seq.len());
}

#[test]
fn test_extend_hit_ungapped_short_sequences() {
    let q_seq = b"AC";
    let s_seq = b"AC";
    
    let (q_start, q_end, s_start, s_end, score) = extend_hit_ungapped(q_seq, s_seq, 0, 0, 1, -2, None);
    
    // Should handle short sequences correctly
    assert!(q_start <= 1);
    assert!(q_end <= q_seq.len());
    assert!(score >= 0);
}

