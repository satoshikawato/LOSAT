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
//! Reference: NCBI BLAST source code (c++/src/algo/blast/core/blast_stat.c)
//!   - Function: `BLAST_ComputeLengthAdjustment`
//!   - Location: `c++/src/algo/blast/core/blast_stat.c`
//!   - The implementation matches NCBI BLAST exactly, with only minor differences
//!     due to floating-point precision between C and Rust.

use super::tables::KarlinParams;
use crate::utils::glibc_fma_libm::GlibcLibm;

/// Result of length adjustment computation
///
/// NCBI's `BLAST_ComputeLengthAdjustment` returns 0 when `converged`, else 1
/// (blast_stat.c:5125) and its only caller, `BLAST_CalcEffLengths`
/// (blast_setup.c:821), drops it; no LOSAT caller read it either.
#[derive(Debug, Clone, Copy)]
pub struct LengthAdjustmentResult {
    /// The computed length adjustment value (integer, as in NCBI BLAST)
    pub length_adjustment: i64,
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
/// A `LengthAdjustmentResult` containing the adjustment value.
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
    let log_k = k.glibc_ln();
    let alpha_d_lambda = params.alpha / params.lambda;
    let beta = params.beta;

    // Handle edge cases
    if m <= 0.0 || n <= 0.0 || k <= 0.0 || params.lambda <= 0.0 {
        return LengthAdjustmentResult {
            length_adjustment: 0,
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
        };
    }

    let discriminant = mb * mb - 4.0 * a * c;
    if discriminant < 0.0 {
        return LengthAdjustmentResult {
            length_adjustment: 0,
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
        let ell_bar = alpha_d_lambda * (log_k + ss.glibc_ln()) + beta;

        if ell_bar >= ell {
            // ell is no bigger than the true fixed point
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_max {
                // Stop iterating when the lower and upper bounds meet (no more points to check).
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
            if alpha_d_lambda * (log_k + ss.glibc_ln()) + beta >= ell_ceil {
                // ceil(ell_min) == floor(ell_fixed)
                length_adjustment = ell_ceil as i64;
            }
        }
    }
    // Note: NCBI BLAST sets length_adjustment even if not converged (uses ell_min)
    // Our implementation already does this above

    LengthAdjustmentResult { length_adjustment }
}

#[cfg(test)]
#[path = "length_adjustment_tests.rs"]
mod tests;
