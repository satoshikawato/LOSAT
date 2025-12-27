//! NCBI-compatible length adjustment calculation for Karlin-Altschul statistics.
//!
//! This module implements the length adjustment algorithm from NCBI BLAST's
//! `BLAST_ComputeLengthAdjustment` function. The length adjustment accounts for
//! the fact that alignments cannot extend beyond sequence boundaries.
//!
//! **Implementation Status**: This is a faithful port of NCBI BLAST's implementation.
//! The code has been directly translated from C to Rust, line by line, to ensure
//! exact compatibility with NCBI BLAST's behavior.
//!
//! Reference: NCBI BLAST source code (blast_stat.c)
//!   - Function: `BLAST_ComputeLengthAdjustment`
//!   - Location: `c++/src/algo/blast/core/blast_stat.c:5041-5126`
//!   - The implementation matches NCBI BLAST exactly, including the comparison
//!     `ell_min == ell_max` (line 5093) and direct cast `(Int4) ell_min` (lines 5110, 5122)

use super::tables::KarlinParams;

/// Result of length adjustment computation
#[derive(Debug, Clone, Copy)]
pub struct LengthAdjustmentResult {
    /// The computed length adjustment value (integer, as in NCBI BLAST)
    pub length_adjustment: i64,
    /// Whether the iteration converged
    pub converged: bool,
}

/// Compute the length adjustment using NCBI BLAST's algorithm.
///
/// This is a faithful port of `BLAST_ComputeLengthAdjustment` from NCBI BLAST.
///
/// The algorithm finds the fixed point of:
///   ell = alpha_d_lambda * (logK + log(ss)) + beta
/// where:
///   ss = (m - ell) * (n - N * ell)
///   m = query_length
///   n = db_length (total)
///   N = db_num_seqs
///
/// The iteration uses bisection-like bounds to ensure convergence.
///
/// # Arguments
/// * `query_length` - Length of the query sequence
/// * `db_length` - Total length of the database
/// * `db_num_seqs` - Number of sequences in the database
/// * `params` - Karlin-Altschul parameters (lambda, K, alpha, beta)
///
/// # Returns
/// A `LengthAdjustmentResult` containing the adjustment value and convergence status.
pub fn compute_length_adjustment_ncbi(
    query_length: i64,
    db_length: i64,
    db_num_seqs: i64,
    params: &KarlinParams,
) -> LengthAdjustmentResult {
    const MAX_ITERATIONS: i32 = 20;

    let m = query_length as f64;
    let n = db_length as f64;
    let n_seqs = db_num_seqs as f64;

    let k = params.k;
    let log_k = k.ln();
    let alpha_d_lambda = params.alpha / params.lambda;
    let beta = params.beta;

    // Handle edge cases
    if m <= 0.0 || n <= 0.0 || k <= 0.0 || params.lambda <= 0.0 {
        return LengthAdjustmentResult {
            length_adjustment: 0,
            converged: false,
        };
    }

    // Choose ell_max to be the largest nonnegative value that satisfies:
    //   K * (m - ell) * (n - N * ell) > MAX(m, n)
    //
    // Use quadratic formula: 2c / (-b + sqrt(b*b - 4*a*c))
    // where a = N, b = -(m*N + n), c = n*m - MAX(m,n)/K
    let a = n_seqs;
    let mb = m * n_seqs + n; // This is -b (negated)
    let c = n * m - m.max(n) / k;

    if c < 0.0 {
        // No valid length adjustment possible
        return LengthAdjustmentResult {
            length_adjustment: 0,
            converged: true,
        };
    }

    let discriminant = mb * mb - 4.0 * a * c;
    if discriminant < 0.0 {
        return LengthAdjustmentResult {
            length_adjustment: 0,
            converged: false,
        };
    }

    let ell_max = 2.0 * c / (mb + discriminant.sqrt());
    let mut ell_min = 0.0_f64;
    let mut ell_max = ell_max;
    let mut ell_next = 0.0_f64;
    let mut converged = false;

    for i in 1..=MAX_ITERATIONS {
        let ell = ell_next;
        let ss = (m - ell) * (n - n_seqs * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.ln()) + beta;

        if ell_bar >= ell {
            // ell is no bigger than the true fixed point
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            // Match NCBI BLAST's exact comparison: ell_min == ell_max
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:5093
            if ell_min == ell_max {
                // There are no more points to check
                break;
            }
        } else {
            // ell is greater than the true fixed point
            ell_max = ell;
        }

        if ell_min <= ell_bar && ell_bar <= ell_max {
            // ell_bar is in range. Accept it
            ell_next = ell_bar;
        } else {
            // ell_bar is not in range. Reject it
            ell_next = if i == 1 {
                ell_max
            } else {
                (ell_min + ell_max) / 2.0
            };
        }
    }

    // Determine the final length adjustment
    // Match NCBI BLAST's exact behavior: (Int4) ell_min (direct cast, not floor)
    // Note: f64 as i64 truncates toward zero, which is equivalent to floor for positive values
    let mut length_adjustment = ell_min as i64;

    if converged {
        // If ell_fixed is the (unknown) true fixed point, then we
        // wish to set length_adjustment to floor(ell_fixed).
        // We assume that floor(ell_min) = floor(ell_fixed)
        // But verify that ceil(ell_min) != floor(ell_fixed)
        // This matches NCBI BLAST's logic exactly
        let ell_ceil = ell_min.ceil();
        if ell_ceil <= ell_max {
            let ss = (m - ell_ceil) * (n - n_seqs * ell_ceil);
            if alpha_d_lambda * (log_k + ss.ln()) + beta >= ell_ceil {
                // ceil(ell_min) == floor(ell_fixed)
                length_adjustment = ell_ceil as i64;
            }
        }
    }
    // Note: NCBI BLAST sets length_adjustment even if not converged (uses ell_min)
    // Our implementation already does this above

    LengthAdjustmentResult {
        length_adjustment,
        converged,
    }
}

/// Simplified length adjustment for single-sequence comparisons.
///
/// This is a convenience wrapper for when comparing a query against a single subject
/// (db_num_seqs = 1).
pub fn compute_length_adjustment_simple(
    query_length: i64,
    subject_length: i64,
    params: &KarlinParams,
) -> LengthAdjustmentResult {
    compute_length_adjustment_ncbi(query_length, subject_length, 1, params)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_blosum62_params() -> KarlinParams {
        KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        }
    }

    fn default_blastn_params() -> KarlinParams {
        KarlinParams {
            lambda: 1.28,
            k: 0.46,
            h: 0.85,
            alpha: 1.5,
            beta: -2.0,
        }
    }

    #[test]
    fn test_length_adjustment_basic() {
        let params = default_blosum62_params();
        let result = compute_length_adjustment_ncbi(100, 10000, 10, &params);

        // Should converge
        assert!(result.converged);
        // Length adjustment should be positive
        assert!(result.length_adjustment >= 0);
        // Length adjustment should be less than query length
        assert!(result.length_adjustment < 100);
    }

    #[test]
    fn test_length_adjustment_single_sequence() {
        let params = default_blosum62_params();
        let result = compute_length_adjustment_simple(100, 1000, &params);

        assert!(result.converged);
        assert!(result.length_adjustment >= 0);
        assert!(result.length_adjustment < 100);
    }

    #[test]
    fn test_length_adjustment_nucleotide() {
        let params = default_blastn_params();
        let result = compute_length_adjustment_ncbi(1000, 100000, 100, &params);

        assert!(result.converged);
        assert!(result.length_adjustment >= 0);
    }

    #[test]
    fn test_length_adjustment_edge_cases() {
        let params = default_blosum62_params();

        // Very short query
        let result = compute_length_adjustment_ncbi(10, 10000, 10, &params);
        assert!(result.length_adjustment >= 0);
        assert!(result.length_adjustment <= 10);

        // Single sequence database
        let result = compute_length_adjustment_ncbi(100, 100, 1, &params);
        assert!(result.length_adjustment >= 0);
    }

    #[test]
    fn test_length_adjustment_zero_inputs() {
        let params = default_blosum62_params();

        // Zero query length
        let result = compute_length_adjustment_ncbi(0, 10000, 10, &params);
        assert_eq!(result.length_adjustment, 0);
        assert!(!result.converged);

        // Zero database length
        let result = compute_length_adjustment_ncbi(100, 0, 10, &params);
        assert_eq!(result.length_adjustment, 0);
        assert!(!result.converged);
    }

    #[test]
    fn test_length_adjustment_includes_beta() {
        // Test that beta is properly included in the calculation
        let params_with_beta = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0, // Negative beta
        };

        let params_zero_beta = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: 0.0, // Zero beta
        };

        let result_with_beta = compute_length_adjustment_ncbi(100, 10000, 10, &params_with_beta);
        let result_zero_beta = compute_length_adjustment_ncbi(100, 10000, 10, &params_zero_beta);

        // With negative beta, the length adjustment should be smaller
        assert!(result_with_beta.length_adjustment <= result_zero_beta.length_adjustment);
    }

    #[test]
    fn test_length_adjustment_db_num_seqs_effect() {
        let params = default_blosum62_params();

        // Same total database length, different number of sequences
        let result_few_seqs = compute_length_adjustment_ncbi(100, 10000, 10, &params);
        let result_many_seqs = compute_length_adjustment_ncbi(100, 10000, 100, &params);

        // More sequences should generally result in different (often smaller) adjustment
        // because the constraint K*(m-A)*(n-N*A) > MAX(m,n) is tighter
        assert!(result_few_seqs.converged);
        assert!(result_many_seqs.converged);
    }
}
