//! E-value calculation numerical precision tests
//!
//! This module tests the numerical precision of E-value calculations,
//! comparing LOSAT's implementation with NCBI's expected behavior.
//!
//! Reference: NCBI BLAST source code
//!   - blast_stat.c: BLAST_SmallGapSumE, s_BlastSumP
//!   - ncbi_math.c: BLAST_LnFactorial, BLAST_LnGammaInt

use LOSAT::stats::sum_statistics::{
    e_to_p, large_gap_sum_e, ln_factorial, ln_gamma_int, normalize_score, p_to_e, small_gap_sum_e,
};

/// Test ln_factorial_int precision for small n (typical case: 2-10 HSPs)
#[test]
fn test_ln_factorial_int_small() {
    // For small n, direct calculation should be very accurate
    let test_cases = vec![
        (1, 0.0),
        (2, 2.0_f64.ln()),
        (3, 6.0_f64.ln()),
        (5, 120.0_f64.ln()),
        (10, 3628800.0_f64.ln()),
    ];

    for (n, expected) in test_cases {
        let result = ln_factorial(n as f64);
        if expected == 0.0 {
            // Special case: both result and expected are 0
            assert_eq!(result, expected, "ln_factorial({}) should be 0", n);
        } else {
            let relative_error = ((result - expected) / expected).abs();
            // For small n, expect very high precision
            assert!(
                relative_error < 1e-10,
                "ln_factorial({}) = {}, expected {}, relative error = {}",
                n,
                result,
                expected,
                relative_error
            );
        }
    }
}

/// Test ln_factorial_int precision for medium n (50-100 HSPs)
#[test]
fn test_ln_factorial_int_medium() {
    // For medium n, compare with Stirling's approximation
    let test_cases = vec![50, 75, 100];

    for n in test_cases {
        let direct = ln_factorial(n as f64);
        // Stirling's approximation: ln(n!) ≈ n*ln(n) - n + 0.5*ln(2*pi*n)
        let stirling = (n as f64) * (n as f64).ln() - (n as f64)
            + 0.5 * (2.0 * std::f64::consts::PI * n as f64).ln();
        let relative_error = ((direct - stirling) / stirling).abs();
        // For medium n, direct calculation should still be accurate
        // but may have some accumulated error
        // Note: Stirling's approximation itself has error, so we use a more lenient threshold
        assert!(
            relative_error < 1e-4,
            "ln_factorial({}) = {}, Stirling = {}, relative error = {}",
            n,
            direct,
            stirling,
            relative_error
        );
    }
}

/// Test ln_factorial_int precision for large n (400+ HSPs)
#[test]
fn test_ln_factorial_int_large() {
    // For large n, compare with Stirling's approximation
    // This is where NCBI uses lgamma which is more numerically stable
    let test_cases = vec![200, 300, 400, 500];

    for n in test_cases {
        let direct = ln_factorial(n as f64);
        // Stirling's approximation
        let stirling = (n as f64) * (n as f64).ln() - (n as f64)
            + 0.5 * (2.0 * std::f64::consts::PI * n as f64).ln();
        let relative_error = ((direct - stirling) / stirling).abs();
        // For large n, expect some accumulated error but still reasonable
        // NCBI uses lgamma which is more stable, but direct calculation
        // should still be within acceptable bounds for E-value calculation
        assert!(
            relative_error < 1e-4,
            "ln_factorial({}) = {}, Stirling = {}, relative error = {}",
            n,
            direct,
            stirling,
            relative_error
        );
    }
}

/// Test small_gap_sum_e precision for single HSP
#[test]
fn test_small_gap_sum_e_single_hsp() {
    let xsum = 10.0;
    let searchsp = 1_000_000_i64;
    let result = small_gap_sum_e(50, 1, xsum, 100, 1000, searchsp, 1.0);

    let expected = (searchsp as f64) * (-xsum).exp();
    let relative_error = ((result - expected) / expected).abs();
    assert!(
        relative_error < 1e-10,
        "small_gap_sum_e(single) = {}, expected {}, relative error = {}",
        result,
        expected,
        relative_error
    );
}

/// Test small_gap_sum_e precision for small chains (2-10 HSPs)
#[test]
fn test_small_gap_sum_e_small_chain() {
    // Test with 2, 5, 10 HSPs
    let test_cases = vec![2, 5, 10];

    for num_hsps in test_cases {
        // Use a reasonable xsum value
        let xsum = (num_hsps as f64) * 5.0; // 5.0 per HSP
        let searchsp = 1_000_000_i64;
        let query_len = 1000;
        let subject_len = 10000;
        let starting_points = 50;

        let result = small_gap_sum_e(
            starting_points,
            num_hsps as i16,
            xsum,
            query_len,
            subject_len,
            searchsp,
            1.0,
        );

        // Result should be finite and positive
        assert!(
            result.is_finite() && result > 0.0,
            "small_gap_sum_e({} HSPs) = {}, should be finite and positive",
            num_hsps,
            result
        );
    }
}

/// Test small_gap_sum_e precision for medium chains (50-100 HSPs)
#[test]
fn test_small_gap_sum_e_medium_chain() {
    let test_cases = vec![50, 75, 100];

    for num_hsps in test_cases {
        let xsum = (num_hsps as f64) * 5.0;
        let searchsp = 1_000_000_i64;
        let query_len = 1000;
        let subject_len = 10000;
        let starting_points = 50;

        let result = small_gap_sum_e(
            starting_points,
            num_hsps as i16,
            xsum,
            query_len,
            subject_len,
            searchsp,
            1.0,
        );

        assert!(
            result.is_finite() && result > 0.0,
            "small_gap_sum_e({} HSPs) = {}, should be finite and positive",
            num_hsps,
            result
        );
    }
}

/// Test small_gap_sum_e precision for large chains (400+ HSPs)
/// This is the critical case mentioned in the investigation
#[test]
fn test_small_gap_sum_e_large_chain() {
    let test_cases = vec![200, 300, 400, 500];

    for num_hsps in test_cases {
        let xsum = (num_hsps as f64) * 5.0;
        let searchsp = 1_000_000_i64;
        let query_len = 1000;
        let subject_len = 10000;
        let starting_points = 50;

        let result = small_gap_sum_e(
            starting_points,
            num_hsps as i16,
            xsum,
            query_len,
            subject_len,
            searchsp,
            1.0,
        );

        assert!(
            result.is_finite() && result > 0.0,
            "small_gap_sum_e({} HSPs) = {}, should be finite and positive",
            num_hsps,
            result
        );
    }
}

/// Test xsum accumulation precision
/// Simulates the accumulation of xsum across multiple HSPs
#[test]
fn test_xsum_accumulation_precision() {

    let lambda: f64 = 0.267;
    let k: f64 = 0.041;
    let log_k = k.ln();

    // Simulate accumulating xsum for 400 HSPs
    let num_hsps = 400;
    let mut xsum = 0.0;
    let scores: Vec<i32> = (0..num_hsps).map(|i| 50 + (i % 20) as i32).collect();

    for score in &scores {
        xsum += normalize_score(*score, lambda, log_k);
    }

    // Calculate expected xsum
    let sum_scores: i32 = scores.iter().sum();
    let expected_xsum = lambda * (sum_scores as f64) - (num_hsps as f64) * log_k;

    // Compare accumulated vs direct calculation
    let relative_error = ((xsum - expected_xsum) / expected_xsum.abs()).abs();
    assert!(
        relative_error < 1e-10,
        "xsum accumulation: accumulated = {}, expected = {}, relative error = {}",
        xsum,
        expected_xsum,
        relative_error
    );
}

/// Test p_to_e and e_to_p round-trip precision
#[test]
fn test_p_to_e_round_trip() {
    let test_cases = vec![0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999];

    for p in test_cases {
        let e = p_to_e(p);
        let p_back = e_to_p(e);
        let relative_error = ((p - p_back) / p).abs();
        assert!(
            relative_error < 1e-10,
            "p_to_e round-trip: p = {}, e = {}, p_back = {}, relative error = {}",
            p,
            e,
            p_back,
            relative_error
        );
    }
}

/// Test edge cases for E-value calculation
#[test]
fn test_evalue_edge_cases() {
    // Very small E-value
    let result1 = small_gap_sum_e(50, 1, 100.0, 100, 1000, 1_000_000, 1.0);
    assert!(result1 > 0.0 && result1.is_finite());

    // Very large xsum (should give very small E-value)
    let result2 = small_gap_sum_e(50, 1, 200.0, 100, 1000, 1_000_000, 1.0);
    assert!(result2 > 0.0 && result2.is_finite());
    assert!(result2 < result1, "Larger xsum should give smaller E-value");

    // Very small xsum (should give large E-value)
    let result3 = small_gap_sum_e(50, 1, 1.0, 100, 1000, 1_000_000, 1.0);
    assert!(result3 > 0.0 && result3.is_finite());
    assert!(result3 > result1, "Smaller xsum should give larger E-value");
}

