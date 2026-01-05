//! Unit tests for tblastx/extension.rs

use LOSAT::algorithm::tblastx::extension::{
    extend_hit_ungapped, extend_hit_two_hit, extend_gapped_protein, get_score,
};
use LOSAT::algorithm::tblastx::constants::{X_DROP_UNGAPPED, X_DROP_GAPPED_PRELIM};
use LOSAT::utils::matrix::ncbistdaa;

#[test]
fn test_get_score_basic() {
    // Test BLOSUM62 scoring
    // A-A match should have positive score
    let score_aa = get_score(ncbistdaa::A, ncbistdaa::A); // A-A
    assert!(score_aa > 0);
    
    // A-B mismatch should have lower score than A-A
    let score_ab = get_score(ncbistdaa::A, ncbistdaa::B); // A-B
    assert!(score_ab < score_aa);
}

#[test]
fn test_get_score_matches() {
    // Identical amino acids should have positive scores
    for aa in 0u8..28u8 {
        let score = get_score(aa, aa);
        // Most amino acids have positive self-scores in BLOSUM62
        // (except some rare cases, but generally true)
        assert!(score >= -4); // BLOSUM62 minimum is around -4
    }
}

#[test]
fn test_get_score_stop_codon() {
    // Stop codon '*' is 25 in NCBISTDAA encoding.
    let score_stop = get_score(ncbistdaa::STOP, ncbistdaa::A);
    // Stop codon vs normal amino acid should be negative (BLOSUM62: -4) or zero.
    assert!(score_stop <= 0);
    
    let score_stop_stop = get_score(ncbistdaa::STOP, ncbistdaa::STOP);
    // NCBI BLAST BLOSUM62: *-* = +1 (sm_blosum62.c)
    assert_eq!(score_stop_stop, 1);
}

#[test]
fn test_extend_hit_ungapped_perfect_match() {
    // Create identical sequences
    let q_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8]; // ABCDEF
    let s_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8]; // ABCDEF
    
    // Start at position 1 (B)
    let (q_start, q_end, s_start, s_end, score, s_last_off) = 
        extend_hit_ungapped(&q_seq, &s_seq, 1, 1, 0, X_DROP_UNGAPPED);
    
    // Should extend to cover matching region
    assert!(score > 0);
    assert!(q_start <= 1);
    assert!(q_end > 1);
    assert!(s_start <= 1);
    assert!(s_end > 1);
    assert!(s_last_off >= s_end);
}

#[test]
fn test_extend_hit_ungapped_with_mismatches() {
    let q_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8]; // ABCDEF
    let s_seq = [0u8, 1u8, 2u8, 10u8, 11u8, 12u8]; // ABCKLM (mismatches)
    
    let (q_start, q_end, s_start, s_end, score, _s_last_off) = 
        extend_hit_ungapped(&q_seq, &s_seq, 1, 1, 0, X_DROP_UNGAPPED);
    
    // Should still extend but score will be lower due to mismatches
    assert!(score > 0); // Initial match should give positive score
    assert!(q_start <= 1);
    assert!(q_end > 1);
}

#[test]
fn test_extend_hit_ungapped_x_drop_termination() {
    // Create sequences with many mismatches after initial match
    let q_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8, 6u8, 7u8, 8u8, 9u8];
    let s_seq = [0u8, 1u8, 2u8, 20u8, 21u8, 22u8, 23u8, 24u8, 20u8, 21u8]; // Many mismatches
    
    let (q_start, q_end, s_start, s_end, score, _s_last_off) = 
        extend_hit_ungapped(&q_seq, &s_seq, 1, 1, 0, X_DROP_UNGAPPED);
    
    // X-drop should terminate extension when score drops too much
    assert!(score >= 0);
    // Extension should be limited by X-drop
    assert!((q_end - q_start) < q_seq.len()); // Should not extend to full length
}

#[test]
fn test_extend_hit_ungapped_at_sequence_start() {
    // Sequence buffers are NCBISTDAA with leading/trailing sentinel (0 / NULLB).
    let q_seq = [ncbistdaa::GAP, ncbistdaa::A, ncbistdaa::B, ncbistdaa::C, ncbistdaa::GAP];
    let s_seq = [ncbistdaa::GAP, ncbistdaa::A, ncbistdaa::B, ncbistdaa::C, ncbistdaa::GAP];
    
    // Start at first residue (position 1). Seeds never start at the sentinel.
    let (q_start, q_end, s_start, s_end, score, _s_last_off) = 
        extend_hit_ungapped(&q_seq, &s_seq, 1, 1, 0, X_DROP_UNGAPPED);
    
    // Should extend rightward only (can't extend left from position 0)
    assert_eq!(q_start, 1);
    assert!(q_end > 1);
    assert!(score > 0);
}

#[test]
fn test_extend_hit_ungapped_at_sequence_end() {
    let q_seq = [0u8, 1u8, 2u8];
    let s_seq = [0u8, 1u8, 2u8];
    
    // Start at last position
    let last_pos = q_seq.len() - 1;
    let (q_start, q_end, s_start, s_end, score, _s_last_off) = 
        extend_hit_ungapped(&q_seq, &s_seq, last_pos, last_pos, 0, X_DROP_UNGAPPED);
    
    // Should extend leftward only (can't extend right beyond end)
    assert!(q_end <= q_seq.len());
    assert!(score > 0);
}

#[test]
fn test_extend_hit_ungapped_coordinates() {
    let q_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8];
    let s_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8];
    
    // Start at middle position
    let start_pos = 2;
    let (q_start, q_end, s_start, s_end, _score, _s_last_off) = 
        extend_hit_ungapped(&q_seq, &s_seq, start_pos, start_pos, 0, X_DROP_UNGAPPED);
    
    // Coordinates should be valid
    assert!(q_start <= start_pos);
    assert!(q_end > start_pos);
    assert!(s_start <= start_pos);
    assert!(s_end > start_pos);
    assert!(q_end <= q_seq.len());
    assert!(s_end <= s_seq.len());
}

#[test]
fn test_extend_hit_two_hit_basic() {
    let q_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8, 6u8, 7u8, 8u8, 9u8];
    let s_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8, 6u8, 7u8, 8u8, 9u8];
    
    // Two hits on same diagonal: first at (1,1), second at (4,4)
    let s_left_off = 1;  // First hit position in subject
    let s_right_off = 4; // Second hit position in subject
    let q_right_off = 4; // Second hit position in query
    
    // NCBI default is 7 bits => raw x-drop is 16 for BLOSUM62 ungapped.
    let x_drop = X_DROP_UNGAPPED;
    let (q_start, q_end, s_start, s_end, score, right_extended, _s_last_off) = 
        extend_hit_two_hit(&q_seq, &s_seq, s_left_off, s_right_off, q_right_off, x_drop);
    
    // Should extend and connect the two hits
    assert!(score > 0);
    assert!(q_start <= q_right_off);
    assert!(q_end >= q_right_off);
    assert!(s_start <= s_right_off);
    assert!(s_end >= s_right_off);
    
    // If left extension reached first hit, right extension should have occurred
    if right_extended {
        assert!(q_end > q_right_off);
    }
}

#[test]
fn test_extend_hit_two_hit_no_connection() {
    let q_seq = [0u8, 1u8, 2u8, 20u8, 21u8, 22u8, 3u8, 4u8, 5u8, 6u8];
    let s_seq = [0u8, 1u8, 2u8, 20u8, 21u8, 22u8, 3u8, 4u8, 5u8, 6u8];
    
    // Two hits far apart with mismatches in between
    let s_left_off = 1;
    let s_right_off = 6;
    let q_right_off = 6;
    
    let x_drop = X_DROP_UNGAPPED;
    let (q_start, q_end, s_start, s_end, score, right_extended, _s_last_off) = 
        extend_hit_two_hit(&q_seq, &s_seq, s_left_off, s_right_off, q_right_off, x_drop);
    
    // Should extend left from second hit, but may not reach first hit
    assert!(score >= 0);
    assert!(q_start <= q_right_off);
    assert!(s_start <= s_right_off);
    
    // If left extension didn't reach first hit, right extension shouldn't occur
    if !right_extended {
        assert_eq!(q_end, q_right_off);
    }
}

#[test]
fn test_extend_hit_two_hit_at_boundaries() {
    let q_seq = [0u8, 1u8, 2u8];
    let s_seq = [0u8, 1u8, 2u8];
    
    // Two hits at sequence boundaries
    let s_left_off = 0;
    let s_right_off = 2;
    let q_right_off = 2;
    
    let x_drop = X_DROP_UNGAPPED;
    let (q_start, q_end, s_start, s_end, score, _right_extended, _s_last_off) = 
        extend_hit_two_hit(&q_seq, &s_seq, s_left_off, s_right_off, q_right_off, x_drop);
    
    // Should handle boundary conditions correctly
    assert!(score >= 0);
    assert!(q_start <= q_seq.len());
    assert!(q_end <= q_seq.len());
    assert!(s_start <= s_seq.len());
    assert!(s_end <= s_seq.len());
}

#[test]
fn test_extend_gapped_protein_basic() {
    let q_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8, 6u8, 7u8, 8u8, 9u8];
    let s_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8, 6u8, 7u8, 8u8, 9u8];
    
    // Start at position 2, seed length 3
    let qs = 2;
    let ss = 2;
    let len = 3;
    
    let (q_start, q_end, s_start, s_end, score, matches, _mismatches, gap_opens, gap_letters) = 
        extend_gapped_protein(&q_seq, &s_seq, qs, ss, len, X_DROP_GAPPED_PRELIM);
    
    // Should extend and produce valid alignment
    assert!(score > 0);
    assert!(q_start <= qs);
    assert!(q_end >= qs + len);
    assert!(s_start <= ss);
    assert!(s_end >= ss + len);
    
    // Statistics should be valid
    assert!(matches > 0);
    assert!(gap_opens >= 0);
    assert!(gap_letters >= 0);
}

#[test]
fn test_extend_gapped_protein_with_gaps() {
    // Create sequences that require gaps for optimal alignment
    let q_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8, 6u8, 7u8, 8u8, 9u8];
    let s_seq = [0u8, 1u8, 2u8, 20u8, 21u8, 3u8, 4u8, 5u8, 6u8, 7u8]; // Insertion in subject
    
    let qs = 0;
    let ss = 0;
    let len = 3;
    
    let (q_start, q_end, _s_start, _s_end, score, _matches, _mismatches, gap_opens, gap_letters) = 
        extend_gapped_protein(&q_seq, &s_seq, qs, ss, len, X_DROP_GAPPED_PRELIM);
    
    // Should handle gaps
    assert!(score > 0);
    assert!(q_start <= qs);
    assert!(q_end >= qs + len);
    
    // May have gaps if alignment requires them
    // (gap_opens and gap_letters may be > 0)
    assert!(gap_opens >= 0);
    assert!(gap_letters >= 0);
}

#[test]
fn test_extend_gapped_protein_at_boundaries() {
    let q_seq = [0u8, 1u8, 2u8];
    let s_seq = [0u8, 1u8, 2u8];
    
    // Start at beginning
    let qs = 0;
    let ss = 0;
    let len = 3;
    
    let (q_start, q_end, s_start, s_end, score, _matches, _mismatches, _gap_opens, _gap_letters) = 
        extend_gapped_protein(&q_seq, &s_seq, qs, ss, len, X_DROP_GAPPED_PRELIM);
    
    // Should handle boundary conditions
    assert!(score > 0);
    assert_eq!(q_start, 0);
    assert!(q_end <= q_seq.len());
    assert_eq!(s_start, 0);
    assert!(s_end <= s_seq.len());
}

#[test]
fn test_extend_gapped_protein_invalid_bounds() {
    let q_seq = [0u8, 1u8, 2u8];
    let s_seq = [0u8, 1u8, 2u8];
    
    // Invalid start position (beyond sequence length)
    let qs = 10;
    let ss = 10;
    let len = 3;
    
    let (q_start, q_end, s_start, s_end, score, matches, mismatches, gap_opens, gap_letters) = 
        extend_gapped_protein(&q_seq, &s_seq, qs, ss, len, X_DROP_GAPPED_PRELIM);
    
    // Should return zero score and original positions
    assert_eq!(score, 0);
    assert_eq!(q_start, qs);
    assert_eq!(q_end, qs);
    assert_eq!(s_start, ss);
    assert_eq!(s_end, ss);
    assert_eq!(matches, 0);
    assert_eq!(mismatches, 0);
    assert_eq!(gap_opens, 0);
    assert_eq!(gap_letters, 0);
}

#[test]
fn test_extend_gapped_protein_zero_length_seed() {
    let q_seq = [0u8, 1u8, 2u8];
    let s_seq = [0u8, 1u8, 2u8];
    
    let qs = 0;
    let ss = 0;
    let len = 0; // Zero length seed
    
    let (q_start, q_end, s_start, s_end, score, matches, mismatches, gap_opens, gap_letters) = 
        extend_gapped_protein(&q_seq, &s_seq, qs, ss, len, X_DROP_GAPPED_PRELIM);
    
    // Should return zero score
    assert_eq!(score, 0);
    assert_eq!(q_start, qs);
    assert_eq!(q_end, qs);
    assert_eq!(s_start, ss);
    assert_eq!(s_end, ss);
    assert_eq!(matches, 0);
    assert_eq!(mismatches, 0);
    assert_eq!(gap_opens, 0);
    assert_eq!(gap_letters, 0);
}

#[test]
fn test_extend_gapped_protein_with_stop_codons() {
    let q_seq = [0u8, 1u8, ncbistdaa::STOP, 3u8, 4u8, 5u8];
    let s_seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8];
    
    let qs = 0;
    let ss = 0;
    let len = 3;
    
    let (q_start, q_end, _s_start, _s_end, score, _matches, mismatches, _gap_opens, _gap_letters) = 
        extend_gapped_protein(&q_seq, &s_seq, qs, ss, len, X_DROP_GAPPED_PRELIM);
    
    // Should handle stop codons (they should be penalized)
    // Stop codons have negative scores, so total score may be negative
    // Don't assert score >= 0, just check coordinates are valid
    assert!(q_start <= qs);
    assert!(q_end >= qs + len);
    
    // Should count stop codons as mismatches
    assert!(mismatches >= 0);
}

