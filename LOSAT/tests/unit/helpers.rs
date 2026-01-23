//! Test utilities and helpers for unit tests
//!
//! This module provides reusable test utilities such as:
//! - Mock sequence generators
//! - Test data fixtures
//! - Assertion helpers for alignment results
//! - NCBI BLAST output comparison utilities

use LOSAT::common::Hit;
use LOSAT::stats::tables::KarlinParams;

/// Create a test Hit with default values
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
// ```c
// typedef struct BlastHSPList {
//    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
//    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
//                       Set to 0 if not applicable */
// } BlastHSPList;
// ```
pub fn make_hit(
    q_idx: u32,
    s_idx: u32,
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    bit_score: f64,
) -> Hit {
    Hit {
        identity: 90.0,
        length: q_end - q_start + 1,
        mismatch: 0,
        gapopen: 0,
        q_start,
        q_end,
        s_start,
        s_end,
        e_value: 1e-10,
        bit_score,
        q_idx,
        s_idx,
        raw_score: 100,
        gap_info: None,
    }
}

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

/// Generate a simple nucleotide sequence for testing
pub fn make_nucleotide_sequence(length: usize) -> Vec<u8> {
    let bases = b"ACGT";
    (0..length)
        .map(|i| bases[i % bases.len()])
        .collect()
}

/// Generate a simple protein sequence for testing
pub fn make_protein_sequence(length: usize) -> Vec<u8> {
    let amino_acids = b"ACDEFGHIKLMNPQRSTVWY";
    (0..length)
        .map(|i| amino_acids[i % amino_acids.len()])
        .collect()
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

/// NCBI BLAST reference data and comparison utilities
pub mod ncbi_reference;

