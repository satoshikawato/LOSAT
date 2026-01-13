//! Sum statistics for HSP linking.
//!
//! Reference: blast_stat.c s_BlastSumP, BLAST_SmallGapSumE, etc.

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

pub(crate) fn ln_factorial_int(n: i32) -> f64 {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gap_decay_divisor() {
        assert!((gap_decay_divisor(0.5, 1) - 0.5).abs() < 1e-10);
        assert!((gap_decay_divisor(0.5, 2) - 0.25).abs() < 1e-10);
        assert!((gap_decay_divisor(0.1, 1) - 0.9).abs() < 1e-10);
    }
}
