//! Unit tests for tblastx/extension.rs

use crate::algorithm::tblastx::extension::get_score;
use crate::utils::matrix::ncbistdaa;

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
        if aa == ncbistdaa::GAP {
            // NCBI blast core keeps NULLB/GAP rows at BLAST_SCORE_MIN.
            assert_eq!(score, i16::MIN as i32);
        } else {
            assert!(score >= -4);
        }
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
