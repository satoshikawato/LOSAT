//! BLAST Statistical Functions
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c
//!
//! This module provides Karlin-Altschul statistical calculations for BLAST:
//! - Bit score and E-value calculation
//! - Length adjustment for effective search space
//! - Sum statistics for HSP linking
//! - Statistical parameter lookup tables
//!
//! # Organization
//!
//! The module is organized into sections matching NCBI BLAST's blast_stat.c:
//! - Parameter definitions (KarlinParams, SearchSpace)
//! - Score/E-value calculations
//! - Length adjustment algorithm
//! - Sum statistics for multi-HSP E-values
//! - Parameter lookup tables

use crate::config::{NuclScoringSpec, ProteinScoringSpec, ScoringMatrix};
use crate::utils::matrix::{blosum62_score, BLASTAA_SIZE};

// ============================================================================
// KARLIN PARAMETERS
// ============================================================================

/// Karlin-Altschul statistical parameters
#[derive(Debug, Clone, Copy)]
pub struct KarlinParams {
    /// Lambda parameter for bit score calculation
    pub lambda: f64,
    /// K parameter for E-value calculation
    pub k: f64,
    /// H parameter (entropy) for length adjustment
    pub h: f64,
    /// Alpha parameter for length correction mean
    pub alpha: f64,
    /// Beta parameter for length correction
    pub beta: f64,
}

impl Default for KarlinParams {
    fn default() -> Self {
        Self {
            lambda: 0.625,
            k: 0.041,
            h: 0.85,
            alpha: 1.5,
            beta: -2.0,
        }
    }
}

// ============================================================================
// SEARCH SPACE
// ============================================================================

/// Effective search space with length-adjusted values
#[derive(Debug, Clone, Copy)]
pub struct SearchSpace {
    /// Effective query length after adjustment
    pub effective_query_len: f64,
    /// Effective database length after adjustment
    pub effective_db_len: f64,
    /// Effective search space (product of effective lengths)
    pub effective_space: f64,
    /// The length adjustment value used (for debugging/verification)
    pub length_adjustment: i64,
}

impl SearchSpace {
    /// Create a simple search space without length adjustment (BLEMIR default)
    pub fn simple(query_len: usize, db_len: usize) -> Self {
        let q = query_len as f64;
        let d = db_len as f64;
        Self {
            effective_query_len: q,
            effective_db_len: d,
            effective_space: q * d,
            length_adjustment: 0,
        }
    }

    /// Create search space with NCBI-style length adjustment for single-sequence comparison.
    pub fn with_length_adjustment(query_len: usize, db_len: usize, params: &KarlinParams) -> Self {
        let m = query_len as f64;
        let n = db_len as f64;

        let result = compute_length_adjustment_simple(query_len as i64, db_len as i64, params);
        let length_adj = result.length_adjustment as f64;

        let effective_m = (m - length_adj).max(1.0);
        let effective_n = (n - length_adj).max(1.0);

        Self {
            effective_query_len: effective_m,
            effective_db_len: effective_n,
            effective_space: effective_m * effective_n,
            length_adjustment: result.length_adjustment,
        }
    }

    /// Create search space for database search (multiple subjects) using NCBI BLAST algorithm.
    pub fn for_database_search(
        query_len: usize,
        total_db_len: usize,
        num_sequences: usize,
        params: &KarlinParams,
        use_length_adjustment: bool,
    ) -> Self {
        if !use_length_adjustment {
            return Self::simple(query_len, total_db_len);
        }

        let m = query_len as f64;
        let n = total_db_len as f64;
        let n_seqs = num_sequences as f64;

        let result = compute_length_adjustment_ncbi(
            query_len as i64,
            total_db_len as i64,
            num_sequences as i64,
            params,
        );
        let length_adj = result.length_adjustment as f64;

        let effective_m = (m - length_adj).max(1.0);
        let effective_n = (n - (length_adj * n_seqs)).max(1.0);

        Self {
            effective_query_len: effective_m,
            effective_db_len: effective_n,
            effective_space: effective_m * effective_n,
            length_adjustment: result.length_adjustment,
        }
    }
}

/// Compute effective search space for translated searches (TBLASTX)
pub fn compute_tblastx_search_space(
    query_nucl_len: usize,
    subject_nucl_len: usize,
    params: &KarlinParams,
    use_length_adjustment: bool,
) -> SearchSpace {
    let query_aa_len = query_nucl_len / 3;
    let subject_aa_len = subject_nucl_len / 3;

    if use_length_adjustment {
        SearchSpace::with_length_adjustment(query_aa_len, subject_aa_len, params)
    } else {
        SearchSpace::simple(query_aa_len, subject_aa_len)
    }
}

// ============================================================================
// BIT SCORE AND E-VALUE CALCULATIONS
// Reference: blast_stat.c BLAST_KarlinStoE_simple, BlastKarlinEtoS_simple
// ============================================================================

// NCBI blast_stat.c uses kSmallFloat = 1e-297 in BlastKarlinEtoS_simple
const K_SMALL_FLOAT: f64 = 1.0e-297;

/// Calculate bit score from raw score using Karlin-Altschul statistics
///
/// Formula: S' = (lambda * S - ln(K)) / ln(2)
pub fn bit_score(raw_score: i32, params: &KarlinParams) -> f64 {
    let ln2 = 2.0_f64.ln();
    (params.lambda * (raw_score as f64) - params.k.ln()) / ln2
}

/// Calculate E-value from bit score and search space
///
/// Formula: E = m * n * 2^(-S')
pub fn evalue(bit_score: f64, search_space: &SearchSpace) -> f64 {
    search_space.effective_space * 2.0_f64.powf(-bit_score)
}

/// Calculate both bit score and E-value from raw score
pub fn calculate_statistics(
    raw_score: i32,
    params: &KarlinParams,
    search_space: &SearchSpace,
) -> (f64, f64) {
    let bs = bit_score(raw_score, params);
    let ev = evalue(bs, search_space);
    (bs, ev)
}

/// Calculate raw score from E-value (inverse calculation)
///
/// Formula: S = (ln(K) + ln(m*n) - ln(E)) / lambda
pub fn raw_score_from_evalue(e_value: f64, params: &KarlinParams, search_space: &SearchSpace) -> i32 {
    if e_value <= 0.0 {
        return i32::MAX;
    }
    let e = e_value.max(K_SMALL_FLOAT);
    let score = (params.k.ln() + search_space.effective_space.ln() - e.ln()) / params.lambda;
    score.ceil() as i32
}

/// Calculate raw score from E-value with gap decay adjustment (NCBI BLAST_Cutoffs compatible)
///
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4112-4118
pub fn raw_score_from_evalue_with_decay(
    e_value: f64,
    params: &KarlinParams,
    search_space: &SearchSpace,
    dodecay: bool,
    gap_decay_rate: f64,
) -> i32 {
    if e_value <= 0.0 {
        return i32::MAX;
    }

    let adjusted_e = if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
        e_value * gap_decay_divisor(gap_decay_rate, 1)
    } else {
        e_value
    };

    let adjusted_e = adjusted_e.max(K_SMALL_FLOAT);
    let score = (params.k.ln() + search_space.effective_space.ln() - adjusted_e.ln()) / params.lambda;
    score.ceil() as i32
}

/// Calculate raw score from bit score (inverse calculation)
///
/// Formula: S = (S' * ln(2) + ln(K)) / lambda
pub fn raw_score_from_bit_score(bit_score: f64, params: &KarlinParams) -> i32 {
    let ln2 = 2.0_f64.ln();
    let score = (bit_score * ln2 + params.k.ln()) / params.lambda;
    score.ceil() as i32
}

/// Calculate E-value from raw score and search space (NCBI BLAST_KarlinStoE_simple compatible)
///
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4157-4171
pub fn evalue_from_raw_score(raw_score: i32, params: &KarlinParams, search_space: f64) -> f64 {
    if params.lambda < 0.0 || params.k < 0.0 {
        return -1.0;
    }
    search_space * ((-params.lambda * (raw_score as f64)) + params.k.ln()).exp()
}

/// Simple E-value calculation without length adjustment (BLEMIR default)
pub fn simple_evalue(raw_score: i32, q_len: usize, db_len: usize, params: &KarlinParams) -> (f64, f64) {
    let search_space = (q_len as f64) * (db_len as f64);
    let bs = bit_score(raw_score, params);
    let ev = search_space * 2.0_f64.powf(-bs);
    (bs, ev)
}

// ============================================================================
// LENGTH ADJUSTMENT
// Reference: blast_stat.c BLAST_ComputeLengthAdjustment
// ============================================================================

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

// ============================================================================
// SUM STATISTICS
// Reference: blast_stat.c s_BlastSumP, BLAST_SmallGapSumE, etc.
// ============================================================================

/// Gap decay divisor for weighting E-values when multiple alignments are considered.
///
/// Formula: (1 - decayrate) * decayrate^(nsegs - 1)
pub fn gap_decay_divisor(decay_rate: f64, num_segments: usize) -> f64 {
    if num_segments == 0 {
        return 1.0;
    }
    (1.0 - decay_rate) * decay_rate.powi((num_segments - 1) as i32)
}

/// Natural log of factorial using direct calculation for integers.
pub fn ln_factorial(n: f64) -> f64 {
    if n <= 1.0 {
        return 0.0;
    }
    if (n.fract()).abs() < f64::EPSILON && n <= (i32::MAX as f64) {
        return ln_factorial_int(n as i32);
    }
    n * n.ln() - n + 0.5 * (2.0 * std::f64::consts::PI * n).ln()
}

/// Natural log of gamma function for positive integers.
pub fn ln_gamma_int(n: i32) -> f64 {
    if n <= 0 {
        return f64::INFINITY;
    }
    if n == 1 || n == 2 {
        return 0.0;
    }
    ln_factorial_int(n - 1)
}

/// Convert P-value to E-value.
pub fn p_to_e(p: f64) -> f64 {
    if p < 0.0 || p > 1.0 {
        return i32::MIN as f64;
    }
    if p == 1.0 {
        return i32::MAX as f64;
    }
    -(-p).ln_1p()
}

/// Convert E-value to P-value.
pub fn e_to_p(e: f64) -> f64 {
    if e < 0.0 {
        return 0.0;
    }
    -(-e).exp_m1()
}

fn ln_factorial_int(n: i32) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 2..=n {
        sum += (i as f64).ln();
    }
    sum
}

// Lookup tables for s_BlastSumP interpolation
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

fn blast_sum_p_calc(r: i32, s: f64) -> f64 {
    const SUMP_EPSILON: f64 = 0.002;

    if r == 1 {
        if s > 8.0 {
            return (-s).exp();
        }
        return -(-(-s).exp()).exp_m1();
    }
    if r < 1 {
        return 0.0;
    }

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
        let est_mean = -(rf) * (r1 as f64);
        if s <= est_mean - stddev4 {
            return 1.0;
        }
    }

    let logr = rf.ln();
    let mean = rf * (1.0 - logr) - 0.5;
    if s <= mean - stddev4 {
        return 1.0;
    }

    let (t, mut itmin) = if s >= mean {
        (s + 6.0 * stddev, 1_i32)
    } else {
        (mean + 6.0 * stddev, 2_i32)
    };

    let adj1 = (r - 2) as f64 * logr - ln_gamma_int(r1) - ln_gamma_int(r);
    let num_hsps_minus_2 = r - 2;

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

    let mut inner_integral = |s_var: f64| -> f64 {
        let adj2 = adj1 - s_var;
        let sdvir = s_var / rf;
        let mx = if s_var > 0.0 { sdvir + 3.0 } else { 3.0 };
        let mut outer = |x: f64| outer_integral(x, adj2, sdvir);
        romberg_integrate(&mut outer, 0.0, mx, SUMP_EPSILON, 0, 1)
    };

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
            if idx > 0 && idx < table.len() {
                return a * table[idx - 1] + (1.0 - a) * table[idx];
            }
        }
        return 1.0;
    }

    blast_sum_p_calc(r, s)
}

fn romberg_integrate<F>(f: &mut F, p: f64, q: f64, eps: f64, epsit: i32, itmin: i32) -> f64
where
    F: FnMut(f64) -> f64,
{
    const MAX_DIAGS: usize = 20;

    let mut itmin = itmin.max(1);
    itmin = itmin.min((MAX_DIAGS - 1) as i32);

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
    romb[0] = 0.5 * h * (x0 + y0);

    let mut epsit_cnt: i32 = 0;
    for i in 1..MAX_DIAGS {
        let mut sum = 0.0;
        let mut x = p + 0.5 * h;
        for _k in 0..npts {
            let y = f(x);
            if !y.is_finite() {
                return y;
            }
            sum += y;
            x += h;
        }
        romb[i] = 0.5 * (romb[i - 1] + h * sum);

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
        let prod = (subject_length as f64) * (query_length as f64);
        let adjusted_xsum = xsum - (num as f64) * prod.ln() + ln_factorial_int(num);

        let sum_p = blast_sum_p(num, adjusted_xsum);
        sum_e = p_to_e(sum_p) * ((searchsp_eff as f64) / prod);
    }

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
pub fn normalize_score(raw_score: i32, lambda: f64, log_k: f64) -> f64 {
    lambda * (raw_score as f64) - log_k
}

/// NCBI BLAST default parameters for HSP linking
pub mod defaults {
    pub const GAP_PROB_UNGAPPED: f64 = 0.5;
    pub const GAP_PROB_GAPPED: f64 = 1.0;
    pub const GAP_DECAY_RATE_UNGAPPED: f64 = 0.5;
    pub const GAP_DECAY_RATE_GAPPED: f64 = 0.1;
    pub const GAP_SIZE: i32 = 40;
    pub const OVERLAP_SIZE: i32 = 9;
    pub const WINDOW_SIZE: i32 = GAP_SIZE + OVERLAP_SIZE + 1;
}

// ============================================================================
// KARLIN PARAMETER CALCULATION FROM COMPOSITION
// Reference: blast_stat.c Blast_KarlinBlkUngappedCalc
// ============================================================================

/// Standard amino acid frequencies (Robinson probabilities)
const STD_AA_FREQS: [f64; 20] = [
    0.07805, 0.01926, 0.05364, 0.06295, 0.01487, 0.03374, 0.06661, 0.07129, 0.02105, 0.05142,
    0.05744, 0.05068, 0.01471, 0.03965, 0.04728, 0.06141, 0.05506, 0.01330, 0.03216, 0.06891,
];

const NCBISTDAA_TO_STD_AA: [usize; 28] = [
    20, 0, 2, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 20, 18, 5, 20, 20, 20, 10,
];

/// Compute amino acid composition from sequence
pub fn compute_aa_composition(seq: &[u8], length: usize) -> [f64; BLASTAA_SIZE] {
    let mut comp: [u64; BLASTAA_SIZE] = [0; BLASTAA_SIZE];

    let start = 1;
    let end = seq.len().saturating_sub(1);
    let actual_len = end.saturating_sub(start).min(length);

    for i in start..end.min(start + actual_len) {
        let residue = seq[i] as usize;
        if residue < BLASTAA_SIZE {
            comp[residue] += 1;
        }
    }

    let sum: u64 = comp.iter().sum();
    let mut freq = [0.0; BLASTAA_SIZE];

    if sum > 0 {
        for i in 0..BLASTAA_SIZE {
            freq[i] = comp[i] as f64 / sum as f64;
        }
    }

    freq
}

/// Compute standard amino acid composition
pub fn compute_std_aa_composition() -> [f64; BLASTAA_SIZE] {
    let mut freq = [0.0; BLASTAA_SIZE];

    for (ncbi_idx, &std_idx) in NCBISTDAA_TO_STD_AA.iter().enumerate() {
        if std_idx < 20 {
            freq[ncbi_idx] = STD_AA_FREQS[std_idx];
        }
    }

    let sum: f64 = freq.iter().sum();
    if sum > 0.0 {
        for f in freq.iter_mut() {
            *f /= sum;
        }
    }

    freq
}

/// Score frequency profile
#[derive(Debug, Clone)]
pub struct ScoreFreqProfile {
    sprob: Vec<f64>,
    obs_min: i32,
    obs_max: i32,
    score_avg: f64,
    score_min: i32,
}

impl ScoreFreqProfile {
    pub fn new(score_min: i32, score_max: i32) -> Self {
        let range = (score_max - score_min + 1) as usize;
        Self {
            sprob: vec![0.0; range],
            obs_min: 0,
            obs_max: 0,
            score_avg: 0.0,
            score_min,
        }
    }

    pub fn get_prob(&self, score: i32) -> f64 {
        let idx = (score - self.score_min) as usize;
        if idx < self.sprob.len() { self.sprob[idx] } else { 0.0 }
    }

    pub fn set_prob(&mut self, score: i32, prob: f64) {
        let idx = (score - self.score_min) as usize;
        if idx < self.sprob.len() { self.sprob[idx] = prob; }
    }

    pub fn obs_min(&self) -> i32 { self.obs_min }
    pub fn obs_max(&self) -> i32 { self.obs_max }
    pub fn score_avg(&self) -> f64 { self.score_avg }
    pub fn sprob(&self) -> &[f64] { &self.sprob }
}

/// Compute score frequency profile from two amino acid compositions
pub fn compute_score_freq_profile(
    comp1: &[f64; BLASTAA_SIZE],
    comp2: &[f64; BLASTAA_SIZE],
    score_min: i32,
    score_max: i32,
) -> ScoreFreqProfile {
    let mut sfp = ScoreFreqProfile::new(score_min, score_max);

    for score in score_min..=score_max {
        sfp.set_prob(score, 0.0);
    }

    for i in 0..BLASTAA_SIZE {
        for j in 0..BLASTAA_SIZE {
            let score = blosum62_score(i as u8, j as u8);
            if score >= score_min {
                let prob = comp1[i] * comp2[j];
                let current = sfp.get_prob(score);
                sfp.set_prob(score, current + prob);
            }
        }
    }

    let mut score_sum = 0.0;
    let mut obs_min = score_max;
    let mut obs_max = score_min;

    for score in score_min..=score_max {
        let prob = sfp.get_prob(score);
        if prob > 0.0 {
            score_sum += prob;
            obs_max = obs_max.max(score);
            if obs_min == score_max { obs_min = score; }
        }
    }

    sfp.obs_min = if obs_min <= obs_max { obs_min } else { 0 };
    sfp.obs_max = if obs_min <= obs_max { obs_max } else { 0 };

    let mut score_avg = 0.0;
    if score_sum > 0.0001 || score_sum < -0.0001 {
        for score in score_min..=score_max {
            let prob = sfp.get_prob(score);
            if prob > 0.0 { score_avg += (score as f64) * prob; }
        }
        score_avg /= score_sum;
    }

    sfp.score_avg = score_avg;
    sfp
}

fn gcd(a: i32, b: i32) -> i32 {
    let mut a = a.abs();
    let mut b = b.abs();
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

fn compute_lambda_nr(sfp: &ScoreFreqProfile, initial_guess: f64) -> Result<f64, String> {
    const LAMBDA_ACCURACY: f64 = 1e-5;
    const LAMBDA_ITER_MAX: i32 = 17;

    let low = sfp.obs_min();
    let high = sfp.obs_max();

    if sfp.score_avg() >= 0.0 {
        return Err("Expected score must be negative".to_string());
    }

    if low >= high {
        return Err("Invalid score range".to_string());
    }

    let sprob = sfp.sprob();
    let score_min = sfp.score_min;

    let mut d = (-low).abs();
    for i in 1..=(high - low) {
        let score = low + i;
        let idx = (score - score_min) as usize;
        if idx < sprob.len() && sprob[idx] > 0.0 {
            d = gcd(d, i);
            if d == 1 { break; }
        }
    }

    let mut lambda = initial_guess;

    for _iter in 0..LAMBDA_ITER_MAX {
        let mut sum = 0.0;
        let mut sum_deriv = 0.0;

        for score in low..=high {
            let idx = (score - score_min) as usize;
            if idx < sprob.len() {
                let prob = sprob[idx];
                if prob > 0.0 {
                    let exp_term = (lambda * score as f64).exp();
                    sum += prob * exp_term;
                    sum_deriv += prob * exp_term * (score as f64);
                }
            }
        }

        if sum_deriv.abs() < 1e-10 { break; }

        let f = sum - 1.0;
        let f_prime = sum_deriv;
        let lambda_new = lambda - f / f_prime;

        if (lambda_new - lambda).abs() < LAMBDA_ACCURACY {
            lambda = lambda_new;
            break;
        }

        lambda = lambda_new;
    }

    let mut sum = 0.0;
    for score in low..=high {
        let idx = (score - score_min) as usize;
        if idx < sprob.len() {
            let prob = sprob[idx];
            if prob > 0.0 { sum += prob * (lambda * score as f64).exp(); }
        }
    }

    if (sum - 1.0).abs() > LAMBDA_ACCURACY {
        return Err(format!("Lambda convergence failed: sum={}", sum));
    }

    Ok(lambda)
}

fn compute_h_from_lambda(sfp: &ScoreFreqProfile, lambda: f64) -> Result<f64, String> {
    if lambda < 0.0 {
        return Err("Lambda must be non-negative".to_string());
    }

    let low = sfp.obs_min();
    let high = sfp.obs_max();

    if low >= high {
        return Err("Invalid score range".to_string());
    }

    let sprob = sfp.sprob();
    let score_min = sfp.score_min;
    let etonlam = (-lambda).exp();

    let mut sum = (low as f64) * sprob[(low - score_min) as usize];
    for score in (low + 1)..=high {
        let idx = (score - score_min) as usize;
        if idx < sprob.len() {
            sum = (score as f64) * sprob[idx] + etonlam * sum;
        }
    }

    let scale = etonlam.powi(high);
    let h = if scale > 0.0 {
        lambda * sum / scale
    } else {
        lambda * (lambda * high as f64 + sum.ln()).exp()
    };

    Ok(h)
}

fn compute_k_from_lambda_h(sfp: &ScoreFreqProfile, lambda: f64, _h: f64) -> Result<f64, String> {
    let low = sfp.obs_min();
    let high = sfp.obs_max();

    if low >= 0 || high <= 0 {
        return Err("Scores must span negative and positive values".to_string());
    }

    let k = _h / lambda;

    if k <= 0.0 {
        return Err("Computed K is non-positive".to_string());
    }

    Ok(k)
}

/// Compute Karlin-Altschul parameters from score frequency profile
pub fn compute_karlin_params_ungapped(sfp: &ScoreFreqProfile) -> Result<KarlinParams, String> {
    const LAMBDA0_DEFAULT: f64 = 0.5;

    let lambda = compute_lambda_nr(sfp, LAMBDA0_DEFAULT)?;
    let h = compute_h_from_lambda(sfp, lambda)?;
    let k = compute_k_from_lambda_h(sfp, lambda, h)?;

    Ok(KarlinParams {
        lambda,
        k,
        h,
        alpha: 0.7916,
        beta: -3.2,
    })
}

/// Apply check_ideal logic: use ideal params if computed Lambda >= ideal Lambda
pub fn apply_check_ideal(computed: KarlinParams, ideal: KarlinParams) -> KarlinParams {
    if computed.lambda >= ideal.lambda { ideal } else { computed }
}

// ============================================================================
// PARAMETER LOOKUP TABLES
// Reference: blast_stat.c BLAST_GetNuclValuesArray, BLAST_GetProteinGapParams
// ============================================================================

#[derive(Debug, Clone, Copy)]
struct ParamEntry {
    gap_open: i32,
    gap_extend: i32,
    lambda: f64,
    k: f64,
    h: f64,
    alpha: f64,
    beta: f64,
}

impl ParamEntry {
    const fn new(gap_open: i32, gap_extend: i32, lambda: f64, k: f64, h: f64, alpha: f64, beta: f64) -> Self {
        Self { gap_open, gap_extend, lambda, k, h, alpha, beta }
    }

    fn to_karlin_params(&self) -> KarlinParams {
        KarlinParams { lambda: self.lambda, k: self.k, h: self.h, alpha: self.alpha, beta: self.beta }
    }
}

// Nucleotide parameters
const BLASTN_1_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.39, 0.747, 1.38, 1.00, 0.0),
    ParamEntry::new(3, 3, 1.39, 0.747, 1.38, 1.00, 0.0),
];

const BLASTN_1_4: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.383, 0.738, 1.36, 1.02, 0.0),
    ParamEntry::new(1, 2, 1.36, 0.67, 1.2, 1.1, 0.0),
    ParamEntry::new(0, 2, 1.26, 0.43, 0.90, 1.4, -1.0),
    ParamEntry::new(2, 1, 1.35, 0.61, 1.1, 1.2, -1.0),
    ParamEntry::new(1, 1, 1.22, 0.35, 0.72, 1.7, -3.0),
];

const BLASTN_2_7: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.69, 0.73, 1.34, 0.515, 0.0),
    ParamEntry::new(2, 4, 0.68, 0.67, 1.2, 0.55, 0.0),
    ParamEntry::new(0, 4, 0.63, 0.43, 0.90, 0.7, -1.0),
    ParamEntry::new(4, 2, 0.675, 0.62, 1.1, 0.6, -1.0),
    ParamEntry::new(2, 2, 0.61, 0.35, 0.72, 1.7, -3.0),
];

const BLASTN_1_3: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.374, 0.711, 1.31, 1.05, 0.0),
    ParamEntry::new(2, 2, 1.37, 0.70, 1.2, 1.1, 0.0),
    ParamEntry::new(1, 2, 1.35, 0.64, 1.1, 1.2, -1.0),
    ParamEntry::new(0, 2, 1.25, 0.42, 0.83, 1.5, -2.0),
    ParamEntry::new(2, 1, 1.34, 0.60, 1.1, 1.2, -1.0),
    ParamEntry::new(1, 1, 1.21, 0.34, 0.71, 1.7, -2.0),
];

const BLASTN_2_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.675, 0.65, 1.1, 0.6, -1.0),
    ParamEntry::new(2, 4, 0.67, 0.59, 1.1, 0.6, -1.0),
    ParamEntry::new(0, 4, 0.62, 0.39, 0.78, 0.8, -2.0),
    ParamEntry::new(4, 2, 0.67, 0.61, 1.0, 0.65, -2.0),
    ParamEntry::new(2, 2, 0.56, 0.32, 0.59, 0.95, -4.0),
];

const BLASTN_1_2: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.28, 0.46, 0.85, 1.5, -2.0),
    ParamEntry::new(2, 2, 1.33, 0.62, 1.1, 1.2, 0.0),
    ParamEntry::new(1, 2, 1.30, 0.52, 0.93, 1.4, -2.0),
    ParamEntry::new(0, 2, 1.19, 0.34, 0.66, 1.8, -3.0),
    ParamEntry::new(3, 1, 1.32, 0.57, 1.0, 1.3, -1.0),
    ParamEntry::new(2, 1, 1.29, 0.49, 0.92, 1.4, -1.0),
    ParamEntry::new(1, 1, 1.14, 0.26, 0.52, 2.2, -5.0),
];

const BLASTN_2_3: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.55, 0.21, 0.46, 1.2, -5.0),
    ParamEntry::new(4, 4, 0.63, 0.42, 0.84, 0.75, -2.0),
    ParamEntry::new(2, 4, 0.615, 0.37, 0.72, 0.85, -3.0),
    ParamEntry::new(0, 4, 0.55, 0.21, 0.46, 1.2, -5.0),
    ParamEntry::new(3, 3, 0.615, 0.37, 0.68, 0.9, -3.0),
    ParamEntry::new(6, 2, 0.63, 0.42, 0.84, 0.75, -2.0),
    ParamEntry::new(5, 2, 0.625, 0.41, 0.78, 0.8, -2.0),
    ParamEntry::new(4, 2, 0.61, 0.35, 0.68, 0.9, -3.0),
    ParamEntry::new(2, 2, 0.515, 0.14, 0.33, 1.55, -9.0),
];

const BLASTN_3_4: &[ParamEntry] = &[
    ParamEntry::new(6, 3, 0.389, 0.25, 0.56, 0.7, -5.0),
    ParamEntry::new(5, 3, 0.375, 0.21, 0.47, 0.8, -6.0),
    ParamEntry::new(4, 3, 0.351, 0.14, 0.35, 1.0, -9.0),
    ParamEntry::new(6, 2, 0.362, 0.16, 0.45, 0.8, -4.0),
    ParamEntry::new(5, 2, 0.330, 0.092, 0.28, 1.2, -13.0),
    ParamEntry::new(4, 2, 0.281, 0.046, 0.16, 1.8, -23.0),
];

const BLASTN_4_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.22, 0.061, 0.22, 1.0, -15.0),
    ParamEntry::new(6, 5, 0.28, 0.21, 0.47, 0.6, -7.0),
    ParamEntry::new(5, 5, 0.27, 0.17, 0.39, 0.7, -9.0),
    ParamEntry::new(4, 5, 0.25, 0.10, 0.31, 0.8, -10.0),
    ParamEntry::new(3, 5, 0.23, 0.065, 0.25, 0.9, -11.0),
];

const BLASTN_1_1: &[ParamEntry] = &[
    ParamEntry::new(3, 2, 1.09, 0.31, 0.55, 2.0, -2.0),
    ParamEntry::new(2, 2, 1.07, 0.27, 0.49, 2.2, -3.0),
    ParamEntry::new(1, 2, 1.02, 0.21, 0.36, 2.8, -6.0),
    ParamEntry::new(0, 2, 0.80, 0.064, 0.17, 4.8, -16.0),
    ParamEntry::new(4, 1, 1.08, 0.28, 0.54, 2.0, -2.0),
    ParamEntry::new(3, 1, 1.06, 0.25, 0.46, 2.3, -4.0),
    ParamEntry::new(2, 1, 0.99, 0.17, 0.30, 3.3, -10.0),
];

const BLASTN_3_2: &[ParamEntry] = &[ParamEntry::new(5, 5, 0.208, 0.030, 0.072, 2.9, -47.0)];

const BLASTN_5_4: &[ParamEntry] = &[
    ParamEntry::new(10, 6, 0.163, 0.068, 0.16, 1.0, -19.0),
    ParamEntry::new(8, 6, 0.146, 0.039, 0.11, 1.3, -29.0),
];

// Protein parameters
const BLOSUM45: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.2291, 0.0924, 0.2514, 0.9113, -5.7),
    ParamEntry::new(13, 3, 0.207, 0.049, 0.14, 1.5, -22.0),
    ParamEntry::new(12, 3, 0.199, 0.039, 0.11, 1.8, -34.0),
    ParamEntry::new(11, 3, 0.190, 0.031, 0.095, 2.0, -38.0),
    ParamEntry::new(10, 3, 0.179, 0.023, 0.075, 2.4, -51.0),
    ParamEntry::new(16, 2, 0.210, 0.051, 0.14, 1.5, -24.0),
    ParamEntry::new(15, 2, 0.203, 0.041, 0.12, 1.7, -31.0),
    ParamEntry::new(14, 2, 0.195, 0.032, 0.10, 1.9, -36.0),
    ParamEntry::new(13, 2, 0.185, 0.024, 0.084, 2.2, -45.0),
    ParamEntry::new(12, 2, 0.171, 0.016, 0.061, 2.8, -65.0),
    ParamEntry::new(19, 1, 0.205, 0.040, 0.11, 1.9, -43.0),
    ParamEntry::new(18, 1, 0.198, 0.032, 0.10, 2.0, -43.0),
    ParamEntry::new(17, 1, 0.189, 0.024, 0.079, 2.4, -57.0),
    ParamEntry::new(16, 1, 0.176, 0.016, 0.063, 2.8, -67.0),
];

const BLOSUM50: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.2318, 0.112, 0.3362, 0.6895, -4.0),
    ParamEntry::new(13, 3, 0.212, 0.063, 0.19, 1.1, -16.0),
    ParamEntry::new(12, 3, 0.206, 0.055, 0.17, 1.2, -18.0),
    ParamEntry::new(11, 3, 0.197, 0.042, 0.14, 1.4, -25.0),
    ParamEntry::new(10, 3, 0.186, 0.031, 0.11, 1.7, -34.0),
    ParamEntry::new(9, 3, 0.172, 0.022, 0.082, 2.1, -48.0),
    ParamEntry::new(16, 2, 0.215, 0.066, 0.20, 1.05, -15.0),
    ParamEntry::new(15, 2, 0.210, 0.058, 0.17, 1.2, -20.0),
    ParamEntry::new(14, 2, 0.202, 0.045, 0.14, 1.4, -27.0),
    ParamEntry::new(13, 2, 0.193, 0.035, 0.12, 1.6, -32.0),
    ParamEntry::new(12, 2, 0.181, 0.025, 0.095, 1.9, -41.0),
    ParamEntry::new(19, 1, 0.212, 0.057, 0.18, 1.2, -21.0),
    ParamEntry::new(18, 1, 0.207, 0.050, 0.15, 1.4, -28.0),
    ParamEntry::new(17, 1, 0.198, 0.037, 0.12, 1.6, -33.0),
    ParamEntry::new(16, 1, 0.186, 0.025, 0.10, 1.9, -42.0),
    ParamEntry::new(15, 1, 0.171, 0.015, 0.063, 2.7, -76.0),
];

const BLOSUM62: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2),
    ParamEntry::new(11, 2, 0.297, 0.082, 0.27, 1.1, -10.0),
    ParamEntry::new(10, 2, 0.291, 0.075, 0.23, 1.3, -15.0),
    ParamEntry::new(9, 2, 0.279, 0.058, 0.19, 1.5, -19.0),
    ParamEntry::new(8, 2, 0.264, 0.045, 0.15, 1.8, -26.0),
    ParamEntry::new(7, 2, 0.239, 0.027, 0.10, 2.5, -46.0),
    ParamEntry::new(6, 2, 0.201, 0.012, 0.061, 3.3, -58.0),
    ParamEntry::new(13, 1, 0.292, 0.071, 0.23, 1.2, -11.0),
    ParamEntry::new(12, 1, 0.283, 0.059, 0.19, 1.5, -19.0),
    ParamEntry::new(11, 1, 0.267, 0.041, 0.14, 1.9, -30.0),
    ParamEntry::new(10, 1, 0.243, 0.024, 0.10, 2.5, -44.0),
    ParamEntry::new(9, 1, 0.206, 0.010, 0.052, 4.0, -87.0),
];

const BLOSUM80: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3430, 0.177, 0.6568, 0.5222, -1.6),
    ParamEntry::new(25, 2, 0.342, 0.17, 0.66, 0.52, -1.6),
    ParamEntry::new(13, 2, 0.336, 0.15, 0.57, 0.59, -3.0),
    ParamEntry::new(9, 2, 0.319, 0.11, 0.42, 0.76, -6.0),
    ParamEntry::new(8, 2, 0.308, 0.090, 0.35, 0.89, -9.0),
    ParamEntry::new(7, 2, 0.293, 0.070, 0.27, 1.1, -14.0),
    ParamEntry::new(6, 2, 0.268, 0.045, 0.19, 1.4, -19.0),
    ParamEntry::new(11, 1, 0.314, 0.095, 0.35, 0.90, -9.0),
    ParamEntry::new(10, 1, 0.299, 0.071, 0.27, 1.1, -14.0),
    ParamEntry::new(9, 1, 0.279, 0.048, 0.20, 1.4, -19.0),
];

const BLOSUM90: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3346, 0.190, 0.7547, 0.4434, -1.4),
    ParamEntry::new(9, 2, 0.310, 0.12, 0.46, 0.67, -6.0),
    ParamEntry::new(8, 2, 0.300, 0.099, 0.39, 0.76, -7.0),
    ParamEntry::new(7, 2, 0.283, 0.072, 0.30, 0.93, -11.0),
    ParamEntry::new(6, 2, 0.259, 0.048, 0.22, 1.2, -16.0),
    ParamEntry::new(11, 1, 0.302, 0.093, 0.39, 0.78, -8.0),
    ParamEntry::new(10, 1, 0.290, 0.075, 0.28, 1.04, -15.0),
    ParamEntry::new(9, 1, 0.265, 0.044, 0.20, 1.3, -19.0),
];

const PAM30: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3400, 0.283, 1.754, 0.1938, -0.3),
    ParamEntry::new(7, 2, 0.305, 0.15, 0.87, 0.35, -3.0),
    ParamEntry::new(6, 2, 0.287, 0.11, 0.68, 0.42, -4.0),
    ParamEntry::new(5, 2, 0.264, 0.079, 0.45, 0.59, -7.0),
    ParamEntry::new(10, 1, 0.309, 0.15, 0.88, 0.34, -3.0),
    ParamEntry::new(9, 1, 0.294, 0.11, 0.61, 0.48, -6.0),
    ParamEntry::new(8, 1, 0.270, 0.072, 0.40, 0.68, -10.0),
];

const PAM70: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3345, 0.229, 1.029, 0.3250, -0.9),
    ParamEntry::new(8, 2, 0.301, 0.12, 0.54, 0.56, -5.0),
    ParamEntry::new(7, 2, 0.286, 0.093, 0.43, 0.67, -7.0),
    ParamEntry::new(6, 2, 0.264, 0.064, 0.29, 0.90, -12.0),
    ParamEntry::new(11, 1, 0.305, 0.12, 0.52, 0.59, -6.0),
    ParamEntry::new(10, 1, 0.291, 0.091, 0.41, 0.71, -9.0),
    ParamEntry::new(9, 1, 0.270, 0.060, 0.28, 0.97, -14.0),
];

const PAM250: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.2252, 0.0868, 0.2223, 0.98, -5.0),
    ParamEntry::new(15, 3, 0.205, 0.049, 0.13, 1.6, -23.0),
    ParamEntry::new(14, 3, 0.200, 0.043, 0.12, 1.7, -26.0),
    ParamEntry::new(13, 3, 0.194, 0.036, 0.10, 1.9, -31.0),
    ParamEntry::new(12, 3, 0.186, 0.029, 0.085, 2.2, -41.0),
    ParamEntry::new(11, 3, 0.174, 0.020, 0.070, 2.5, -48.0),
    ParamEntry::new(17, 2, 0.204, 0.047, 0.12, 1.7, -28.0),
    ParamEntry::new(16, 2, 0.198, 0.038, 0.11, 1.8, -29.0),
    ParamEntry::new(15, 2, 0.191, 0.031, 0.087, 2.2, -44.0),
    ParamEntry::new(14, 2, 0.182, 0.024, 0.073, 2.5, -53.0),
    ParamEntry::new(13, 2, 0.171, 0.017, 0.059, 2.9, -64.0),
    ParamEntry::new(21, 1, 0.205, 0.045, 0.11, 1.8, -34.0),
    ParamEntry::new(20, 1, 0.199, 0.037, 0.10, 1.9, -35.0),
    ParamEntry::new(19, 1, 0.192, 0.029, 0.083, 2.3, -52.0),
    ParamEntry::new(18, 1, 0.183, 0.021, 0.070, 2.6, -60.0),
    ParamEntry::new(17, 1, 0.171, 0.014, 0.052, 3.3, -86.0),
];

/// Look up Karlin-Altschul parameters for nucleotide scoring scheme
pub fn lookup_nucl_params(spec: &NuclScoringSpec) -> KarlinParams {
    let reward = spec.reward;
    let penalty = spec.penalty.abs();
    let gap_open = spec.gap_open.abs();
    let gap_extend = spec.gap_extend.abs();

    let table: &[ParamEntry] = match (reward, penalty) {
        (1, 5) => BLASTN_1_5,
        (1, 4) => BLASTN_1_4,
        (2, 7) => BLASTN_2_7,
        (1, 3) => BLASTN_1_3,
        (2, 5) => BLASTN_2_5,
        (1, 2) => BLASTN_1_2,
        (2, 3) => BLASTN_2_3,
        (3, 4) => BLASTN_3_4,
        (4, 5) => BLASTN_4_5,
        (1, 1) => BLASTN_1_1,
        (3, 2) => BLASTN_3_2,
        (5, 4) => BLASTN_5_4,
        _ => {
            return KarlinParams { lambda: 1.28, k: 0.46, h: 0.85, alpha: 1.5, beta: -2.0 };
        }
    };

    for entry in table {
        if entry.gap_open == gap_open && entry.gap_extend == gap_extend {
            return entry.to_karlin_params();
        }
    }

    if !table.is_empty() {
        return table[0].to_karlin_params();
    }

    KarlinParams::default()
}

/// Look up Karlin-Altschul parameters for protein scoring scheme
pub fn lookup_protein_params(spec: &ProteinScoringSpec) -> KarlinParams {
    let gap_open = spec.gap_open;
    let gap_extend = spec.gap_extend;

    let table: &[ParamEntry] = match spec.matrix {
        ScoringMatrix::Blosum45 => BLOSUM45,
        ScoringMatrix::Blosum50 => BLOSUM50,
        ScoringMatrix::Blosum62 => BLOSUM62,
        ScoringMatrix::Blosum80 => BLOSUM80,
        ScoringMatrix::Blosum90 => BLOSUM90,
        ScoringMatrix::Pam30 => PAM30,
        ScoringMatrix::Pam70 => PAM70,
        ScoringMatrix::Pam250 => PAM250,
    };

    for entry in table {
        if entry.gap_open == gap_open && entry.gap_extend == gap_extend {
            return entry.to_karlin_params();
        }
    }

    if !table.is_empty() {
        for entry in table {
            if entry.gap_open != i32::MAX {
                return entry.to_karlin_params();
            }
        }
        return table[0].to_karlin_params();
    }

    KarlinParams { lambda: 0.267, k: 0.041, h: 0.14, alpha: 1.9, beta: -30.0 }
}

/// Check if scores need to be rounded down for even-score-only matrices.
pub fn requires_even_scores(reward: i32, penalty: i32) -> bool {
    let penalty = penalty.abs();
    matches!((reward, penalty), (2, 7) | (2, 5) | (2, 3) | (3, 4))
}

/// Look up UNGAPPED Karlin-Altschul parameters for protein scoring scheme.
pub fn lookup_protein_params_ungapped(matrix: ScoringMatrix) -> KarlinParams {
    let table: &[ParamEntry] = match matrix {
        ScoringMatrix::Blosum45 => BLOSUM45,
        ScoringMatrix::Blosum50 => BLOSUM50,
        ScoringMatrix::Blosum62 => BLOSUM62,
        ScoringMatrix::Blosum80 => BLOSUM80,
        ScoringMatrix::Blosum90 => BLOSUM90,
        ScoringMatrix::Pam30 => PAM30,
        ScoringMatrix::Pam70 => PAM70,
        ScoringMatrix::Pam250 => PAM250,
    };

    for entry in table {
        if entry.gap_open == i32::MAX && entry.gap_extend == i32::MAX {
            return entry.to_karlin_params();
        }
    }

    if !table.is_empty() {
        return table[0].to_karlin_params();
    }

    KarlinParams { lambda: 0.3176, k: 0.134, h: 0.4012, alpha: 0.7916, beta: -3.2 }
}

/// Look up GAPPED Karlin-Altschul parameters for protein scoring scheme.
pub fn lookup_protein_params_gapped(matrix: ScoringMatrix) -> KarlinParams {
    let table: &[ParamEntry] = match matrix {
        ScoringMatrix::Blosum45 => BLOSUM45,
        ScoringMatrix::Blosum50 => BLOSUM50,
        ScoringMatrix::Blosum62 => BLOSUM62,
        ScoringMatrix::Blosum80 => BLOSUM80,
        ScoringMatrix::Blosum90 => BLOSUM90,
        ScoringMatrix::Pam30 => PAM30,
        ScoringMatrix::Pam70 => PAM70,
        ScoringMatrix::Pam250 => PAM250,
    };

    for entry in table {
        if entry.gap_open == 11 && entry.gap_extend == 1 {
            return entry.to_karlin_params();
        }
    }

    for entry in table {
        if entry.gap_open != i32::MAX {
            return entry.to_karlin_params();
        }
    }

    KarlinParams { lambda: 0.267, k: 0.041, h: 0.14, alpha: 1.9, beta: -30.0 }
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_score() {
        let params = KarlinParams { lambda: 0.267, k: 0.041, h: 0.14, alpha: 1.9, beta: -30.0 };
        let score = 100;
        let bs = bit_score(score, &params);
        assert!(bs > 0.0);
        let expected = (0.267 * 100.0 - 0.041_f64.ln()) / 2.0_f64.ln();
        assert!((bs - expected).abs() < 0.001);
    }

    #[test]
    fn test_gap_decay_divisor() {
        assert!((gap_decay_divisor(0.5, 1) - 0.5).abs() < 1e-10);
        assert!((gap_decay_divisor(0.5, 2) - 0.25).abs() < 1e-10);
        assert!((gap_decay_divisor(0.1, 1) - 0.9).abs() < 1e-10);
    }

    #[test]
    fn test_length_adjustment_basic() {
        let params = KarlinParams { lambda: 0.267, k: 0.041, h: 0.14, alpha: 1.9, beta: -30.0 };
        let result = compute_length_adjustment_ncbi(100, 10000, 10, &params);
        assert!(result.converged);
        assert!(result.length_adjustment >= 0);
        assert!(result.length_adjustment < 100);
    }

    #[test]
    fn test_lookup_nucl_params_megablast() {
        let spec = NuclScoringSpec { reward: 1, penalty: -2, gap_open: 0, gap_extend: 0 };
        let params = lookup_nucl_params(&spec);
        assert!((params.lambda - 1.28).abs() < 0.0001);
        assert!((params.k - 0.46).abs() < 0.0001);
    }

    #[test]
    fn test_lookup_protein_params_blosum62() {
        let spec = ProteinScoringSpec { matrix: ScoringMatrix::Blosum62, gap_open: 11, gap_extend: 1 };
        let params = lookup_protein_params(&spec);
        assert!((params.lambda - 0.267).abs() < 0.01);
    }
}
