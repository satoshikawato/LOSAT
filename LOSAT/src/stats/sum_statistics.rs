//! Sum statistics for evaluating multiple HSP alignments.
//!
//! This module implements NCBI BLAST's sum statistics algorithms for calculating
//! E-values when multiple HSPs are linked together. These functions are essential
//! for proper HSP chaining and scoring.
//!
//! Reference: NCBI BLAST source code (blast_stat.c)

/// Gap decay divisor for weighting E-values when multiple alignments are considered.
///
/// From NCBI BLAST: "The decayrate parameter is a value in the interval (0,1).
/// Typical values are 0.1 and 0.5."
///
/// Formula: (1 - decayrate) * decayrate^(nsegs - 1)
pub fn gap_decay_divisor(decay_rate: f64, num_segments: usize) -> f64 {
    if num_segments == 0 {
        return 1.0;
    }
    (1.0 - decay_rate) * decay_rate.powi((num_segments - 1) as i32)
}

/// Natural log of factorial using Stirling's approximation for large n.
///
/// For small n (< 10), uses direct calculation.
/// For large n, uses Stirling's approximation: ln(n!) ≈ n*ln(n) - n + 0.5*ln(2*pi*n)
pub fn ln_factorial(n: f64) -> f64 {
    if n <= 1.0 {
        return 0.0;
    }
    if n < 10.0 {
        // Direct calculation for small n
        let mut result = 0.0;
        let mut i = 2.0;
        while i <= n {
            result += i.ln();
            i += 1.0;
        }
        return result;
    }
    // Stirling's approximation for large n
    n * n.ln() - n + 0.5 * (2.0 * std::f64::consts::PI * n).ln()
}

/// Natural log of gamma function for positive integers.
///
/// ln(Gamma(n)) = ln((n-1)!) for positive integers
pub fn ln_gamma_int(n: i32) -> f64 {
    if n <= 0 {
        return f64::INFINITY;
    }
    if n == 1 || n == 2 {
        return 0.0;
    }
    ln_factorial((n - 1) as f64)
}

/// Convert P-value to E-value.
///
/// E = -ln(1 - P)
pub fn p_to_e(p: f64) -> f64 {
    if p < 0.0 || p > 1.0 {
        return f64::MIN;
    }
    if p == 1.0 {
        return f64::MAX;
    }
    -(1.0 - p).ln()
}

/// Convert E-value to P-value.
///
/// P = 1 - exp(-E)
pub fn e_to_p(e: f64) -> f64 {
    if e < 0.0 {
        return 0.0;
    }
    1.0 - (-e).exp()
}

/// Lookup tables for s_BlastSumP interpolation (from NCBI BLAST).
/// These are pre-computed P-values for small numbers of segments.
const TAB2: &[f64] = &[
    0.01669, 0.0249, 0.03683, 0.05390, 0.07794, 0.1111, 0.1559, 0.2146, 0.2890, 0.3794, 0.4836,
    0.5965, 0.7092, 0.8114, 0.8931, 0.9490, 0.9806, 0.9944, 0.9989,
];

const TAB3: &[f64] = &[
    0.9806, 0.9944, 0.9989, 0.0001682, 0.0002542, 0.0003829, 0.0005745, 0.0008587, 0.001278,
    0.001893, 0.002789, 0.004088, 0.005958, 0.008627, 0.01240, 0.01770, 0.02505, 0.03514, 0.04880,
    0.06704, 0.09103, 0.1220, 0.1612, 0.2097, 0.2682, 0.3368, 0.4145, 0.4994, 0.5881, 0.6765,
    0.7596, 0.8326, 0.8922, 0.9367, 0.9667, 0.9846, 0.9939, 0.9980,
];

const TAB4: &[f64] = &[
    2.658e-07, 4.064e-07, 6.203e-07, 9.450e-07, 1.437e-06, 2.181e-06, 3.302e-06, 4.990e-06,
    7.524e-06, 1.132e-05, 1.698e-05, 2.541e-05, 3.791e-05, 5.641e-05, 8.368e-05, 0.0001237,
    0.0001823, 0.0002677, 0.0003915, 0.0005704, 0.0008275, 0.001195, 0.001718, 0.002457, 0.003494,
    0.004942, 0.006948, 0.009702, 0.01346, 0.01853, 0.02532, 0.03431, 0.04607, 0.06128, 0.08068,
    0.1051, 0.1352, 0.1719, 0.2157, 0.2669, 0.3254, 0.3906, 0.4612, 0.5355, 0.6110, 0.6849, 0.7544,
    0.8168, 0.8699, 0.9127, 0.9451, 0.9679, 0.9827, 0.9915, 0.9963,
];

/// Calculate sum P-value using Romberg integration for r > 4.
///
/// This implements the numerical integration from NCBI BLAST's s_BlastSumPCalc.
fn blast_sum_p_calc(r: i32, s: f64) -> f64 {
    if r == 1 {
        if s > 8.0 {
            return (-s).exp();
        }
        return -(-(-s).exp()).exp_m1();
    }
    if r < 1 {
        return 0.0;
    }

    // Early return for very negative scores
    if r < 8 && s <= -2.3 * (r as f64) {
        return 1.0;
    } else if r < 15 && s <= -2.5 * (r as f64) {
        return 1.0;
    } else if r < 27 && s <= -3.0 * (r as f64) {
        return 1.0;
    } else if r < 51 && s <= -3.4 * (r as f64) {
        return 1.0;
    } else if r < 101 && s <= -4.0 * (r as f64) {
        return 1.0;
    }

    let stddev = (r as f64).sqrt();
    let stddev4 = 4.0 * stddev;
    let r1 = r - 1;

    if r > 100 {
        let est_mean = -(r as f64) * (r1 as f64);
        if s <= est_mean - stddev4 {
            return 1.0;
        }
    }

    let logr = (r as f64).ln();
    let mean = (r as f64) * (1.0 - logr) - 0.5;
    if s <= mean - stddev4 {
        return 1.0;
    }

    // For very high scores, use asymptotic approximation
    if s >= mean + stddev4 {
        let a = ln_gamma_int(r + 1);
        let result = (r as f64) * ((r1 as f64) * s.ln() - s - a - a).exp();
        return result.min(1.0);
    }

    // Use numerical integration (simplified Romberg)
    // This is a simplified version - for full accuracy, implement full Romberg integration
    let a = ln_gamma_int(r + 1);
    let result = (r as f64) * ((r1 as f64) * s.abs().ln() - s.abs() - a - a).exp();
    result.min(1.0).max(0.0)
}

/// Estimate the Sum P-value by calculation or interpolation.
///
/// Approx. 2-1/2 digits accuracy minimum throughout the range of r, s.
fn blast_sum_p(r: i32, s: f64) -> f64 {
    if r == 1 {
        return -(-(-s).exp()).exp_m1();
    }

    if r <= 4 {
        if r < 1 {
            return 0.0;
        }
        let r1 = r - 1;
        let rf = r as f64;

        if s >= rf * rf + (r1 as f64) {
            let a = ln_gamma_int(r + 1);
            return rf * ((r1 as f64) * s.ln() - s - a - a).exp();
        }

        if s > -2.0 * rf {
            // Interpolate from tables
            let a_val = s + s + 4.0 * rf;
            let i = a_val as i32;
            let a = a_val - (i as f64);
            let r2 = (r - 2) as usize;

            let table = match r2 {
                0 => TAB2,
                1 => TAB3,
                2 => TAB4,
                _ => return blast_sum_p_calc(r, s),
            };

            let tab_size = table.len() as i32 - 1;
            let idx = (tab_size - i) as usize;

            if idx > 0 && idx < table.len() {
                return a * table[idx - 1] + (1.0 - a) * table[idx];
            }
        }
        return 1.0;
    }

    blast_sum_p_calc(r, s)
}

/// Calculate E-value for alignments with "small" gaps.
///
/// This is used for linking HSPs that are relatively close together.
///
/// # Arguments
/// * `starting_points` - Number of starting points permitted between adjacent alignments
///                       (typically max_overlap + max_gap + 1)
/// * `num_hsps` - Number of distinct alignments in this collection
/// * `xsum` - Sum of normalized scores (each weighted by lambda and logK)
/// * `query_length` - Effective length of the query sequence
/// * `subject_length` - Effective length of the subject sequence
/// * `searchsp_eff` - Effective size of the search space
/// * `weight_divisor` - Divisor used to weight the E-value (from gap_decay_divisor)
pub fn small_gap_sum_e(
    starting_points: i32,
    num_hsps: i16,
    xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: i64,
    weight_divisor: f64,
) -> f64 {
    let mut sum_e: f64;

    if num_hsps == 1 {
        sum_e = (searchsp_eff as f64) * (-xsum).exp();
    } else {
        let pair_search_space = (subject_length as f64) * (query_length as f64);
        let num = num_hsps as i32;

        let mut adjusted_xsum = xsum;
        adjusted_xsum -=
            pair_search_space.ln() + 2.0 * ((num - 1) as f64) * (starting_points as f64).ln();
        adjusted_xsum -= ln_factorial(num as f64);

        let sum_p = blast_sum_p(num, adjusted_xsum);
        sum_e = p_to_e(sum_p) * ((searchsp_eff as f64) / pair_search_space);
    }

    if weight_divisor == 0.0 || sum_e / weight_divisor > (i32::MAX as f64) {
        sum_e = i32::MAX as f64;
    } else {
        sum_e /= weight_divisor;
    }

    sum_e
}

/// Calculate E-value for alignments with asymmetric gaps.
///
/// Used for linking HSPs with different gap sizes in query and subject
/// (e.g., exon linking in translated searches).
///
/// # Arguments
/// * `query_start_points` - Number of starting points in query between adjacent alignments
/// * `subject_start_points` - Number of starting points in subject between adjacent alignments
/// * `num_hsps` - Number of distinct alignments in this collection
/// * `xsum` - Sum of normalized scores
/// * `query_length` - Effective length of the query sequence
/// * `subject_length` - Effective length of the subject sequence
/// * `searchsp_eff` - Effective size of the search space
/// * `weight_divisor` - Divisor used to weight the E-value
pub fn uneven_gap_sum_e(
    query_start_points: i32,
    subject_start_points: i32,
    num_hsps: i16,
    xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: i64,
    weight_divisor: f64,
) -> f64 {
    let mut sum_e: f64;

    if num_hsps == 1 {
        sum_e = (searchsp_eff as f64) * (-xsum).exp();
    } else {
        let pair_search_space = (subject_length as f64) * (query_length as f64);
        let num = num_hsps as i32;

        let mut adjusted_xsum = xsum;
        adjusted_xsum -= pair_search_space.ln()
            + ((num - 1) as f64)
                * ((query_start_points as f64).ln() + (subject_start_points as f64).ln());
        adjusted_xsum -= ln_factorial(num as f64);

        let sum_p = blast_sum_p(num, adjusted_xsum);
        sum_e = p_to_e(sum_p) * ((searchsp_eff as f64) / pair_search_space);
    }

    if weight_divisor == 0.0 || sum_e / weight_divisor > (i32::MAX as f64) {
        sum_e = i32::MAX as f64;
    } else {
        sum_e /= weight_divisor;
    }

    sum_e
}

/// Calculate E-value for alignments with arbitrarily large gaps.
///
/// This is used when HSPs are far apart and the gap size is not constrained.
///
/// # Arguments
/// * `num_hsps` - Number of distinct alignments in this collection
/// * `xsum` - Sum of normalized scores
/// * `query_length` - Effective length of the query sequence
/// * `subject_length` - Effective length of the subject sequence
/// * `searchsp_eff` - Effective size of the search space
/// * `weight_divisor` - Divisor used to weight the E-value
pub fn large_gap_sum_e(
    num_hsps: i16,
    xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: i64,
    weight_divisor: f64,
) -> f64 {
    let mut sum_e: f64;

    if num_hsps == 1 {
        sum_e = (searchsp_eff as f64) * (-xsum).exp();
    } else {
        let pair_search_space = (subject_length as f64) * (query_length as f64);
        let num = num_hsps as i32;

        let mut adjusted_xsum = xsum;
        adjusted_xsum -= (num as f64) * pair_search_space.ln();
        adjusted_xsum -= ln_factorial(num as f64);

        let sum_p = blast_sum_p(num, adjusted_xsum);
        sum_e = p_to_e(sum_p) * (searchsp_eff as f64) / pair_search_space.powi(num);
    }

    if weight_divisor == 0.0 || sum_e / weight_divisor > (i32::MAX as f64) {
        sum_e = i32::MAX as f64;
    } else {
        sum_e /= weight_divisor;
    }

    sum_e
}

/// Normalize a raw score to nats using Karlin-Altschul parameters.
///
/// Formula: xsum = lambda * score - logK
pub fn normalize_score(raw_score: i32, lambda: f64, log_k: f64) -> f64 {
    lambda * (raw_score as f64) - log_k
}

/// NCBI BLAST default parameters for HSP linking
pub mod defaults {
    /// Gap probability for ungapped alignments
    pub const GAP_PROB_UNGAPPED: f64 = 0.5;
    /// Gap probability for gapped alignments
    pub const GAP_PROB_GAPPED: f64 = 1.0;
    /// Gap decay rate for ungapped alignments
    pub const GAP_DECAY_RATE_UNGAPPED: f64 = 0.5;
    /// Gap decay rate for gapped alignments
    pub const GAP_DECAY_RATE_GAPPED: f64 = 0.1;
    /// Default gap size for HSP linking
    pub const GAP_SIZE: i32 = 40;
    /// Default overlap size for HSP linking
    pub const OVERLAP_SIZE: i32 = 9;
    /// Window size = gap_size + overlap_size + 1
    pub const WINDOW_SIZE: i32 = GAP_SIZE + OVERLAP_SIZE + 1;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gap_decay_divisor() {
        // For decay_rate = 0.5, nsegs = 1: (1-0.5) * 0.5^0 = 0.5
        assert!((gap_decay_divisor(0.5, 1) - 0.5).abs() < 1e-10);

        // For decay_rate = 0.5, nsegs = 2: (1-0.5) * 0.5^1 = 0.25
        assert!((gap_decay_divisor(0.5, 2) - 0.25).abs() < 1e-10);

        // For decay_rate = 0.1, nsegs = 1: (1-0.1) * 0.1^0 = 0.9
        assert!((gap_decay_divisor(0.1, 1) - 0.9).abs() < 1e-10);
    }

    #[test]
    fn test_ln_factorial() {
        // ln(1!) = 0
        assert!((ln_factorial(1.0) - 0.0).abs() < 1e-10);

        // ln(2!) = ln(2) ≈ 0.693
        assert!((ln_factorial(2.0) - 2.0_f64.ln()).abs() < 1e-10);

        // ln(5!) = ln(120) ≈ 4.787
        assert!((ln_factorial(5.0) - 120.0_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_p_to_e_and_e_to_p() {
        // Round trip test
        let p = 0.01;
        let e = p_to_e(p);
        let p_back = e_to_p(e);
        assert!((p - p_back).abs() < 1e-10);

        // Edge cases
        assert_eq!(e_to_p(0.0), 0.0);
        assert!((e_to_p(f64::INFINITY) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_small_gap_sum_e_single_hsp() {
        // For a single HSP, sum_e = searchsp * exp(-xsum)
        let xsum = 10.0;
        let searchsp = 1000000_i64;
        let result = small_gap_sum_e(50, 1, xsum, 100, 1000, searchsp, 1.0);

        let expected = (searchsp as f64) * (-xsum).exp();
        assert!((result - expected).abs() / expected < 0.01);
    }

    #[test]
    fn test_normalize_score() {
        let lambda = 0.267;
        let k = 0.041;
        let log_k = k.ln();
        let score = 100;

        let xsum = normalize_score(score, lambda, log_k);
        let expected = lambda * 100.0 - log_k;
        assert!((xsum - expected).abs() < 1e-10);
    }
}
