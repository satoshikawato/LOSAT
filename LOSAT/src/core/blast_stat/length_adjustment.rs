//! Length adjustment algorithms for BLAST.
//!
//! Reference: blast_stat.c BLAST_ComputeLengthAdjustment

use super::karlin_params::KarlinParams;

/// Result of length adjustment computation
#[derive(Debug, Clone, Copy)]
pub struct LengthAdjustmentResult {
    /// The computed length adjustment value
    pub length_adjustment: i64,
    /// Whether the iteration converged
    pub converged: bool,
}

/// Compute the length adjustment using NCBI BLAST's algorithm.
///
/// This is a faithful port of `BLAST_ComputeLengthAdjustment` from NCBI BLAST.
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

    if m <= 0.0 || n <= 0.0 || k <= 0.0 || params.lambda <= 0.0 {
        return LengthAdjustmentResult {
            length_adjustment: 0,
            converged: false,
        };
    }

    let a = n_seqs;
    let mb = m * n_seqs + n;
    let c = n * m - m.max(n) / k;

    if c < 0.0 {
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

    let ell_max_initial = 2.0 * c / (mb + discriminant.sqrt());
    let mut ell_min = 0.0_f64;
    let mut ell_max = ell_max_initial;
    let mut ell_next = 0.0_f64;
    let mut converged = false;

    for i in 1..=MAX_ITERATIONS {
        let ell = ell_next;
        let ss = (m - ell) * (n - n_seqs * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.ln()) + beta;

        if ell_bar >= ell {
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_max {
                break;
            }
        } else {
            ell_max = ell;
        }

        if ell_min <= ell_bar && ell_bar <= ell_max {
            ell_next = ell_bar;
        } else {
            ell_next = if i == 1 {
                ell_max
            } else {
                (ell_min + ell_max) / 2.0
            };
        }
    }

    let mut length_adjustment = ell_min as i64;

    if converged {
        let ell_ceil = ell_min.ceil();
        if ell_ceil <= ell_max {
            let ss = (m - ell_ceil) * (n - n_seqs * ell_ceil);
            if alpha_d_lambda * (log_k + ss.ln()) + beta >= ell_ceil {
                length_adjustment = ell_ceil as i64;
            }
        }
    }

    LengthAdjustmentResult {
        length_adjustment,
        converged,
    }
}

/// Simplified length adjustment for single-sequence comparisons.
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

    #[test]
    fn test_length_adjustment_basic() {
        let params = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };
        let result = compute_length_adjustment_ncbi(100, 10000, 10, &params);
        assert!(result.converged);
        assert!(result.length_adjustment >= 0);
        assert!(result.length_adjustment < 100);
    }
}
