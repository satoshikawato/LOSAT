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

    // In our usage, n is almost always an integer (e.g. number of HSPs),
    // so prefer an exact log-factorial for stability and NCBI-compatibility.
    if (n.fract()).abs() < f64::EPSILON && n <= (i32::MAX as f64) {
        return ln_factorial_int(n as i32);
    }

    // Fallback: Stirling's approximation for non-integers/large values.
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
    // Exact for integers: ln(Gamma(n)) = ln((n-1)!)
    ln_factorial_int(n - 1)
}

/// Convert P-value to E-value.
///
/// E = -ln(1 - P)
pub fn p_to_e(p: f64) -> f64 {
    if p < 0.0 || p > 1.0 {
        return i32::MIN as f64;
    }
    if p == 1.0 {
        return i32::MAX as f64;
    }
    // NCBI: -BLAST_Log1p(-p)
    -(-p).ln_1p()
}

/// Convert E-value to P-value.
///
/// P = 1 - exp(-E)
pub fn e_to_p(e: f64) -> f64 {
    if e < 0.0 {
        return 0.0;
    }
    // NCBI: -BLAST_Expm1(-e)
    -(-e).exp_m1()
}

/// Exact ln(n!) for integers (n >= 0).
///
/// This is used heavily in sum-statistics and we want it to be stable and
/// as close as possible to NCBI's `BLAST_LnFactorial` when called with integer n.
///
/// NCBI reference: ncbi_math.c:473-480
/// ```c
/// double BLAST_LnFactorial (double x) {
///     if( x <= 0.0)
///         return 0.0;
///     else
///         return s_LnGamma(x + 1.0);  // s_LnGamma uses lgamma from C standard library
/// }
/// ```
///
/// Note: NCBI uses lgamma(n+1) for all n, but direct calculation is exact for integers
/// and matches NCBI's precision for typical HSP chain sizes (2-10 HSPs, rarely up to 400+).
/// For very large n (400+), direct calculation may have accumulated error, but tests
/// show it remains within acceptable bounds for E-value calculation.
///
fn ln_factorial_int(n: i32) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    
    // Direct sum calculation: exact for integers and matches NCBI's precision
    // for typical use cases. NCBI uses lgamma(n+1) which is more numerically
    // stable for very large n, but direct calculation is sufficient for our
    // use case (HSP chain sizes typically 2-10, rarely up to 400+).
    let mut sum = 0.0;
    for i in 2..=n {
        sum += (i as f64).ln();
    }
    sum
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
    // Faithful port of NCBI BLAST's s_BlastSumPCalc + BLAST_RombergIntegrate.
    // Reference:
    // - ncbi-blast/c++/src/algo/blast/core/blast_stat.c:s_BlastSumPCalc
    // - ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:BLAST_RombergIntegrate

    const SUMP_EPSILON: f64 = 0.002;

    if r == 1 {
        if s > 8.0 {
            return (-s).exp();
        }
        // -BLAST_Expm1(-exp(-s))
        return -(-(-s).exp()).exp_m1();
    }
    if r < 1 {
        return 0.0;
    }

    // Early return for very negative scores
    let rf = r as f64;
    if r < 8 {
        if s <= -2.3 * rf {
            return 1.0;
        }
    } else if r < 15 {
        if s <= -2.5 * rf {
            return 1.0;
        }
    } else if r < 27 {
        if s <= -3.0 * rf {
            return 1.0;
        }
    } else if r < 51 {
        if s <= -3.4 * rf {
            return 1.0;
        }
    } else if r < 101 {
        if s <= -4.0 * rf {
            return 1.0;
        }
    }

    let stddev = rf.sqrt();
    let stddev4 = 4.0 * stddev;
    let r1 = r - 1;

    if r > 100 {
        // Calculate lower bound on the mean using inequality log(r) <= r
        let est_mean = -(rf) * (r1 as f64);
        if s <= est_mean - stddev4 {
            return 1.0;
        }
    }

    // mean is close to the mode and is readily calculated
    let logr = rf.ln();
    let mean = rf * (1.0 - logr) - 0.5;
    if s <= mean - stddev4 {
        return 1.0;
    }

    // Choose integration upper bound t and minimum iterations itmin
    let (t, mut itmin) = if s >= mean {
        (s + 6.0 * stddev, 1_i32)
    } else {
        (mean + 6.0 * stddev, 2_i32)
    };

    // adj1 = (r-2)*log(r) - lnGamma(r-1) - lnGamma(r)
    let adj1 = (r - 2) as f64 * logr - ln_gamma_int(r1) - ln_gamma_int(r);
    let num_hsps_minus_2 = r - 2;

    // Callback for outer integral
    let mut outer_integral = |x: f64, adj2: f64, sdvir: f64| -> f64 {
        let y = (x - sdvir).exp();
        if !y.is_finite() {
            return 0.0;
        }
        if num_hsps_minus_2 == 0 {
            return (adj2 - y).exp();
        }
        if x == 0.0 {
            return 0.0;
        }
        ((num_hsps_minus_2 as f64) * x.ln() + adj2 - y).exp()
    };

    // Callback for inner integral (calls outer via Romberg)
    let mut inner_integral = |s_var: f64| -> f64 {
        let adj2 = adj1 - s_var;
        let sdvir = s_var / rf;
        let mx = if s_var > 0.0 { sdvir + 3.0 } else { 3.0 };
        let mut outer = |x: f64| outer_integral(x, adj2, sdvir);
        romberg_integrate(&mut outer, 0.0, mx, SUMP_EPSILON, 0, 1)
    };

    // Integrate inner integral from s..t, with potential retries for small s<mean
    loop {
        let d = romberg_integrate(&mut inner_integral, s, t, SUMP_EPSILON, 0, itmin);
        if !d.is_finite() {
            return d;
        }
        if !(s < mean && d < 0.4 && itmin < 4) {
            return d.min(1.0);
        }
        itmin += 1;
    }
}

/// Estimate the Sum P-value by calculation or interpolation.
///
/// Approx. 2-1/2 digits accuracy minimum throughout the range of r, s.
fn blast_sum_p(r: i32, s: f64) -> f64 {
    // Faithful port of NCBI BLAST's s_BlastSumP.
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
            // Interpolate from tables (NCBI indexing)
            let mut a = s + s + (4.0 * rf);
            let mut i = a as i32;
            a -= i as f64;
            let r2 = (r - 2) as usize;

            let (table, tab_size) = match r2 {
                0 => (TAB2, (TAB2.len() as i32) - 1),
                1 => (TAB3, (TAB3.len() as i32) - 1),
                2 => (TAB4, (TAB4.len() as i32) - 1),
                _ => return blast_sum_p_calc(r, s),
            };

            i = tab_size - i;
            let idx = i as usize;
            // NCBI assumes idx in-range for s > -2*r; clamp defensively.
            if idx > 0 && idx < table.len() {
                return a * table[idx - 1] + (1.0 - a) * table[idx];
            }
        }
        return 1.0;
    }

    blast_sum_p_calc(r, s)
}

/// Romberg numerical integrator (NCBI BLAST compatible).
///
/// Reference: `ncbi_math.c:BLAST_RombergIntegrate`.
fn romberg_integrate<F>(f: &mut F, p: f64, q: f64, eps: f64, epsit: i32, itmin: i32) -> f64
where
    F: FnMut(f64) -> f64,
{
    const MAX_DIAGS: usize = 20;

    // itmin = min. no. of iterations to perform
    let mut itmin = itmin.max(1);
    itmin = itmin.min((MAX_DIAGS - 1) as i32);

    // epsit = min. no. of consecutive iterations that must satisfy epsilon
    let mut epsit = epsit.max(1);
    epsit = epsit.min(3);

    let epsck = itmin - epsit;

    let mut romb = [0.0_f64; MAX_DIAGS];
    let mut npts: i32 = 1;
    let mut h = q - p;

    let x0 = f(p);
    if !x0.is_finite() {
        return x0;
    }
    let y0 = f(q);
    if !y0.is_finite() {
        return y0;
    }
    romb[0] = 0.5 * h * (x0 + y0); // trapezoidal rule

    let mut epsit_cnt: i32 = 0;
    for i in 1..MAX_DIAGS {
        let mut sum = 0.0;
        // sum of ordinates for x = p+0.5*h, p+1.5*h, ..., q-0.5*h
        let mut x = p + 0.5 * h;
        for _k in 0..npts {
            let y = f(x);
            if !y.is_finite() {
                return y;
            }
            sum += y;
            x += h;
        }
        romb[i] = 0.5 * (romb[i - 1] + h * sum); // new trapezoidal estimate

        // Update Romberg array with new column
        let mut n: f64 = 4.0;
        for j in (0..i).rev() {
            romb[j] = (n * romb[j + 1] - romb[j]) / (n - 1.0);
            n *= 4.0;
        }

        if (i as i32) > epsck {
            if (romb[1] - romb[0]).abs() > eps * romb[0].abs() {
                epsit_cnt = 0;
            } else {
                epsit_cnt += 1;
                if (i as i32) >= itmin && epsit_cnt >= epsit {
                    return romb[0];
                }
            }
        }

        npts *= 2;
        h *= 0.5;
    }

    f64::INFINITY
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
        adjusted_xsum -= pair_search_space.ln()
            + 2.0 * ((num - 1) as f64) * (starting_points as f64).ln();
        adjusted_xsum -= ln_factorial_int(num);

        let sum_p = blast_sum_p(num, adjusted_xsum);
        sum_e = p_to_e(sum_p) * ((searchsp_eff as f64) / pair_search_space);
    }

    // NCBI: if( weight_divisor == 0.0 || (sum_e /= weight_divisor) > INT4_MAX )
    if weight_divisor == 0.0 {
        sum_e = i32::MAX as f64;
    } else {
        sum_e /= weight_divisor;
        if sum_e > (i32::MAX as f64) {
            sum_e = i32::MAX as f64;
        }
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
        adjusted_xsum -= ln_factorial_int(num);

        let sum_p = blast_sum_p(num, adjusted_xsum);
        sum_e = p_to_e(sum_p) * ((searchsp_eff as f64) / pair_search_space);
    }

    // NCBI: if( weight_divisor == 0.0 || (sum_e /= weight_divisor) > INT4_MAX )
    if weight_divisor == 0.0 {
        sum_e = i32::MAX as f64;
    } else {
        sum_e /= weight_divisor;
        if sum_e > (i32::MAX as f64) {
            sum_e = i32::MAX as f64;
        }
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
        let num = num_hsps as i32;

        // NCBI: xsum -= num*log(subject_length*query_length) - ln_factorial(num)
        // i.e. xsum = xsum - num*log(prod) + ln_factorial(num)
        let prod = (subject_length as f64) * (query_length as f64);
        let adjusted_xsum = xsum - (num as f64) * prod.ln() + ln_factorial_int(num);

        let sum_p = blast_sum_p(num, adjusted_xsum);

        // NCBI: sum_e = PtoE(sum_p) * (searchsp_eff / (query_length*subject_length))
        sum_e = p_to_e(sum_p) * ((searchsp_eff as f64) / prod);
    }

    // NCBI: if( weight_divisor == 0.0 || (sum_e /= weight_divisor) > INT4_MAX )
    if weight_divisor == 0.0 {
        sum_e = i32::MAX as f64;
    } else {
        sum_e /= weight_divisor;
        if sum_e > (i32::MAX as f64) {
            sum_e = i32::MAX as f64;
        }
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
        let lambda: f64 = 0.267;
        let k: f64 = 0.041;
        let log_k = k.ln();
        let score = 100;

        let xsum = normalize_score(score, lambda, log_k);
        let expected = lambda * 100.0 - log_k;
        assert!((xsum - expected).abs() < 1e-10);
    }
}
