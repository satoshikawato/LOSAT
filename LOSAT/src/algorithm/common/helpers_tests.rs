//! Test utilities and helpers for unit tests
//!
//! This module provides reusable test utilities such as:
//! - Mock sequence generators
//! - Test data fixtures
//! - Assertion helpers for alignment results
//! - NCBI BLAST output comparison utilities

use crate::stats::tables::KarlinParams;

/// Create default Karlin-Altschul parameters for nucleotide alignments (BLASTN)
pub fn default_nucl_params() -> KarlinParams {
    KarlinParams {
        lambda: 0.625,
        k: 0.041,
        h: 0.85,
        alpha: 1.5,
        beta: -2.0,
    }
}

/// Create default Karlin-Altschul parameters for protein alignments (TBLASTX)
pub fn default_protein_params() -> KarlinParams {
    KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    }
}

/// Assert that two floating point values are approximately equal
pub fn assert_approx_eq(a: f64, b: f64, epsilon: f64) {
    assert!(
        (a - b).abs() < epsilon,
        "Values not approximately equal: {} vs {} (epsilon: {})",
        a,
        b,
        epsilon
    );
}

/// Assert that two e-values are within acceptable tolerance
pub fn assert_evalue_close(actual: f64, expected: f64, tolerance: f64) {
    // E-values can vary significantly, so we use a relative tolerance
    let relative_diff = (actual - expected).abs() / expected.max(1e-10);
    assert!(
        relative_diff < tolerance,
        "E-values not close: {} vs {} (relative diff: {})",
        actual,
        expected,
        relative_diff
    );
}

/// Assert that two bit scores are close
pub fn assert_bit_score_close(actual: f64, expected: f64, epsilon: f64) {
    assert_approx_eq(actual, expected, epsilon);
}
