//! Unit tests for common/evalue.rs

use LOSAT::algorithm::common::evalue::{
    bit_score, calculate_evalue_alignment_length, calculate_evalue_database_search,
};
use LOSAT::stats::tables::KarlinParams;
use super::super::helpers::{assert_bit_score_close, assert_evalue_close, default_nucl_params, default_protein_params};

#[test]
fn test_bit_score_calculation() {
    let params = default_protein_params();
    let score = 100;

    let bs = bit_score(score, &params);

    // Bit score should be positive for positive raw scores
    assert!(bs > 0.0);

    // Verify the formula: S' = (lambda * S - ln(K)) / ln(2)
    let expected = (params.lambda * (score as f64) - params.k.ln()) / 2.0_f64.ln();
    assert_bit_score_close(bs, expected, 0.001);
}

#[test]
fn test_bit_score_with_nucleotide_params() {
    let params = default_nucl_params();
    let score = 50;

    let bs = bit_score(score, &params);

    // Bit score should be positive
    assert!(bs > 0.0);

    // Verify the formula
    let expected = (params.lambda * (score as f64) - params.k.ln()) / 2.0_f64.ln();
    assert_bit_score_close(bs, expected, 0.001);
}

#[test]
fn test_bit_score_zero_score() {
    let params = default_protein_params();
    let score = 0;

    let bs = bit_score(score, &params);

    // Bit score for zero should be negative (since ln(K) > 0)
    assert!(bs < 0.0);
}

#[test]
fn test_evalue_database_search() {
    let params = default_nucl_params();
    let score = 100;
    let q_len = 1000;
    let db_len = 10000;
    let db_num_seqs = 5;

    let (bit_score, e_value) = calculate_evalue_database_search(score, q_len, db_len, db_num_seqs, &params);

    // Bit score should be positive
    assert!(bit_score > 0.0);

    // E-value should be positive
    assert!(e_value > 0.0);

    // For a high score, e-value should be small
    assert!(e_value < 1.0);
}

#[test]
fn test_evalue_database_search_large_database() {
    let params = default_nucl_params();
    let score = 50;
    let q_len = 1000;
    let db_len = 1_000_000;
    let db_num_seqs = 100;

    let (bit_score, e_value) = calculate_evalue_database_search(score, q_len, db_len, db_num_seqs, &params);

    // E-value should increase with larger database
    assert!(e_value > 0.0);
}

#[test]
fn test_evalue_alignment_length() {
    let params = default_protein_params();
    let score = 100;
    let aln_len = 200;

    let (bit_score, e_value) = calculate_evalue_alignment_length(score, aln_len, &params);

    // Bit score should be positive
    assert!(bit_score > 0.0);

    // E-value should be positive
    assert!(e_value > 0.0);

    // For a high score, e-value should be small
    assert!(e_value < 1.0);
}

#[test]
fn test_evalue_alignment_length_minimum_space() {
    let params = default_protein_params();
    let score = 50;
    let aln_len = 100; // Small alignment

    let (bit_score1, e_value1) = calculate_evalue_alignment_length(score, aln_len, &params);

    // Even with small alignment, minimum space (60M) should be used
    // So e-value should be based on minimum space, not aln_len^2
    let (bit_score2, e_value2) = calculate_evalue_alignment_length(score, 50, &params);

    // Bit scores should be the same (same score, same params)
    assert_bit_score_close(bit_score1, bit_score2, 0.001);

    // E-values should be similar (both using minimum space)
    assert_evalue_close(e_value1, e_value2, 0.1);
}

#[test]
fn test_evalue_alignment_length_large_alignment() {
    let params = default_protein_params();
    let score = 100;
    let aln_len = 1000; // Large alignment

    let (bit_score, e_value) = calculate_evalue_alignment_length(score, aln_len, &params);

    // For large alignments, search space should be aln_len^2 (if > 60M)
    // aln_len^2 = 1,000,000 which is < 60M, so still uses minimum
    let aln_len2 = 10000; // aln_len^2 = 100M > 60M
    let (bit_score2, e_value2) = calculate_evalue_alignment_length(score, aln_len2, &params);

    // Bit scores should be the same
    assert_bit_score_close(bit_score, bit_score2, 0.001);

    // E-value for larger alignment should be larger (larger search space)
    assert!(e_value2 > e_value);
}

#[test]
fn test_evalue_consistency() {
    // Test that e-value calculation is consistent
    let params = default_protein_params();
    let score = 80;

    let (bs1, ev1) = calculate_evalue_alignment_length(score, 200, &params);
    let (bs2, ev2) = calculate_evalue_alignment_length(score, 200, &params);

    // Same inputs should give same outputs
    assert_bit_score_close(bs1, bs2, 0.001);
    assert_evalue_close(ev1, ev2, 0.0001);
}

#[test]
fn test_evalue_monotonicity() {
    // Higher scores should give lower e-values
    let params = default_protein_params();
    let aln_len = 200;

    let (_, ev1) = calculate_evalue_alignment_length(50, aln_len, &params);
    let (_, ev2) = calculate_evalue_alignment_length(100, aln_len, &params);
    let (_, ev3) = calculate_evalue_alignment_length(150, aln_len, &params);

    // E-values should decrease as score increases
    assert!(ev1 > ev2);
    assert!(ev2 > ev3);
}

