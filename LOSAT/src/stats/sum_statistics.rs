//! Sum statistics for evaluating multiple HSP alignments.
//!
//! This module implements NCBI BLAST's sum statistics algorithms for calculating
//! E-values when multiple HSPs are linked together. These functions are essential
//! for proper HSP chaining and scoring.
//!
//! Reference: NCBI BLAST source code (blast_stat.c)

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/ncbi_math.h:152-163
const LOGDERIV_ORDER_MAX: usize = 4;
const POLYGAMMA_ORDER_MAX: usize = LOGDERIV_ORDER_MAX;
const NCBIMATH_PI: f64 = 3.1415926535897932384626433832795;
const NCBIMATH_LN2: f64 = 0.69314718055994530941723212145818;
const NCBIMATH_LNPI: f64 = 1.1447298858494001741434273513531;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:60-62
const DBL_EPSILON: f64 = 2.2204460492503131e-16;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:140-151
const DEFAULT_GAMMA_COEF: [f64; 11] = [
    4.694580336184385e+04,
    -1.560605207784446e+05,
    2.065049568014106e+05,
    -1.388934775095388e+05,
    5.031796415085709e+04,
    -9.601592329182778e+03,
    8.785855930895250e+02,
    -3.155153906098611e+01,
    2.908143421162229e-01,
    -2.319827630494973e-04,
    1.251639670050933e-10,
];

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:296-309
const PRECOMPUTED_FACTORIAL: [f64; 35] = [
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0,
    3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0,
    1307674368000.0, 20922789888000.0, 355687428096000.0,
    6402373705728000.0, 121645100408832000.0, 2432902008176640000.0,
    51090942171709440000.0, 1124000727777607680000.0,
    25852016738884976640000.0, 620448401733239439360000.0,
    15511210043330985984000000.0, 403291461126605635584000000.0,
    10888869450418352160768000000.0, 304888344611713860501504000000.0,
    8841761993739701954543616000000.0, 265252859812191058636308480000000.0,
    8222838654177922817725562880000000.0, 263130836933693530167218012160000000.0,
    8683317618811886495518194401280000000.0, 295232799039604140847618609643520000000.0,
];

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:33-55
fn blast_expm1(x: f64) -> f64 {
    let absx = x.abs();
    if absx > 0.33 {
        return x.exp() - 1.0;
    }
    if absx < 1.0e-16 {
        return x;
    }
    x * (1.0
        + x * (1.0 / 2.0
            + x * (1.0 / 6.0
                + x * (1.0 / 24.0
                    + x * (1.0 / 120.0
                        + x * (1.0 / 720.0
                            + x * (1.0 / 5040.0
                                + x * (1.0 / 40320.0
                                    + x * (1.0 / 362880.0
                                        + x * (1.0 / 3628800.0
                                            + x * (1.0 / 39916800.0
                                                + x * (1.0 / 479001600.0
                                                    + x / 6227020800.0))))))))))))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:64-83
fn blast_log1p(x: f64) -> f64 {
    if x.abs() >= 0.2 {
        return (x + 1.0).ln();
    }
    let mut sum = 0.0;
    let mut y = x;
    let mut i = 0;
    while i < 500 {
        i += 1;
        sum += y / (i as f64);
        if y.abs() < DBL_EPSILON {
            break;
        }
        y *= x;
        i += 1;
        sum -= y / (i as f64);
        if y < DBL_EPSILON {
            break;
        }
        y *= x;
    }
    sum
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:444-470
fn blast_powi(mut x: f64, mut n: i32) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if x == 0.0 {
        if n < 0 {
            return f64::INFINITY;
        }
        return 0.0;
    }
    if n < 0 {
        x = 1.0 / x;
        n = -n;
    }
    let mut y = 1.0;
    while n > 0 {
        if (n & 1) != 0 {
            y *= x;
        }
        n /= 2;
        x *= x;
    }
    y
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:94-137
fn s_log_derivative(order: i32, u: &[f64]) -> f64 {
    if order < 0 || order as usize > LOGDERIV_ORDER_MAX {
        return f64::INFINITY;
    }
    if order > 0 && u[0] == 0.0 {
        return f64::INFINITY;
    }

    let mut y = [0.0; LOGDERIV_ORDER_MAX + 1];
    for i in 1..=order as usize {
        y[i] = u[i] / u[0];
    }

    match order {
        0 => {
            if u[0] > 0.0 {
                u[0].ln()
            } else {
                f64::INFINITY
            }
        }
        1 => y[1],
        2 => y[2] - y[1] * y[1],
        3 => y[3] - 3.0 * y[2] * y[1] + 2.0 * y[1] * y[1] * y[1],
        4 => {
            let tmp = y[1] * y[1];
            let mut value = y[4] - 4.0 * y[3] * y[1] - 3.0 * y[2] * y[2] + 12.0 * y[2] * tmp;
            value -= 6.0 * tmp * tmp;
            value
        }
        _ => f64::INFINITY,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:161-221
fn s_general_ln_gamma(x: f64, order: i32) -> f64 {
    let xx = x - 1.0;
    let xgamma_dim = DEFAULT_GAMMA_COEF.len() as f64;
    let tx = xx + xgamma_dim;
    let mut y = [0.0; POLYGAMMA_ORDER_MAX + 1];

    for i in 0..=order {
        let mut tmp = tx;
        let mut idx = DEFAULT_GAMMA_COEF.len();
        let mut value;
        if i == 0 {
            idx -= 1;
            value = DEFAULT_GAMMA_COEF[idx] / tmp;
            while idx > 0 {
                idx -= 1;
                tmp -= 1.0;
                value += DEFAULT_GAMMA_COEF[idx] / tmp;
            }
        } else {
            idx -= 1;
            value = DEFAULT_GAMMA_COEF[idx] / blast_powi(tmp, i + 1);
            while idx > 0 {
                idx -= 1;
                tmp -= 1.0;
                value += DEFAULT_GAMMA_COEF[idx] / blast_powi(tmp, i + 1);
            }
            let tmp_factorial = blast_factorial(i);
            value *= if i % 2 == 0 { tmp_factorial } else { -tmp_factorial };
        }
        y[i as usize] = value;
    }
    y[0] += 1.0;

    let mut value = s_log_derivative(order, &y);
    let mut tmp = tx + 0.5;
    match order {
        0 => {
            value += (NCBIMATH_LNPI + NCBIMATH_LN2) / 2.0 + (xx + 0.5) * tmp.ln() - tmp;
        }
        1 => {
            value += tmp.ln() - xgamma_dim / tmp;
        }
        2 => {
            value += (tmp + xgamma_dim) / (tmp * tmp);
        }
        3 => {
            value -= (1.0 + 2.0 * xgamma_dim / tmp) / (tmp * tmp);
        }
        4 => {
            value += 2.0 * (1.0 + 3.0 * xgamma_dim / tmp) / (tmp * tmp * tmp);
        }
        _ => {
            tmp = blast_factorial(order - 2)
                * blast_powi(tmp, 1 - order)
                * (1.0 + (order as f64 - 1.0) * xgamma_dim / tmp);
            if order % 2 == 0 {
                value += tmp;
            } else {
                value -= tmp;
            }
        }
    }
    value
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:235-284
fn s_poly_gamma(x: f64, order: i32) -> f64 {
    if order < 0 || order as usize > POLYGAMMA_ORDER_MAX {
        return f64::INFINITY;
    }

    if x >= 1.0 {
        return s_general_ln_gamma(x, order);
    }

    if x < 0.0 {
        let mut value = s_general_ln_gamma(1.0 - x, order);
        if (order - 1) % 2 != 0 {
            value = -value;
        }
        if order == 0 {
            let mut sx = (NCBIMATH_PI * x).sin().abs();
            if (x < -0.1 && ((x.ceil() == x) || sx < 2.0 * DBL_EPSILON)) || sx == 0.0 {
                return f64::INFINITY;
            }
            value += NCBIMATH_LNPI - sx.ln();
        } else {
            let mut y = [0.0; POLYGAMMA_ORDER_MAX + 1];
            let mut tmp = 1.0;
            let mut angle = x * NCBIMATH_PI;
            y[0] = angle.sin();
            for k in 1..=order as usize {
                tmp *= NCBIMATH_PI;
                angle += NCBIMATH_PI / 2.0;
                y[k] = tmp * angle.sin();
            }
            value -= s_log_derivative(order, &y);
        }
        value
    } else {
        let mut value = s_general_ln_gamma(1.0 + x, order);
        if order == 0 {
            if x == 0.0 {
                return f64::INFINITY;
            }
            value -= x.ln();
        } else {
            let tmp = blast_factorial(order - 1) * blast_powi(x, -order);
            value += if order % 2 == 0 { tmp } else { -tmp };
        }
        value
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:292-295
fn s_ln_gamma(x: f64) -> f64 {
    s_poly_gamma(x, 0)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:312-321
fn blast_factorial(n: i32) -> f64 {
    if n < 0 {
        return 0.0;
    }
    if (n as usize) < PRECOMPUTED_FACTORIAL.len() {
        return PRECOMPUTED_FACTORIAL[n as usize];
    }
    (s_ln_gamma(n as f64 + 1.0)).exp()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:323-329
fn blast_ln_gamma_int(n: i32) -> f64 {
    if n > 1 && (n as usize) < PRECOMPUTED_FACTORIAL.len() {
        return PRECOMPUTED_FACTORIAL[(n - 1) as usize].ln();
    }
    s_ln_gamma(n as f64)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:473-480
fn blast_ln_factorial(x: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else {
        s_ln_gamma(x + 1.0)
    }
}

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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4081
    (1.0 - decay_rate) * blast_powi(decay_rate, (num_segments - 1) as i32)
}

/// Natural log of factorial.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:473-480
pub fn ln_factorial(n: f64) -> f64 {
    blast_ln_factorial(n)
}

/// Natural log of gamma function for positive integers.
///
/// ln(Gamma(n)) = ln((n-1)!) for positive integers
pub fn ln_gamma_int(n: i32) -> f64 {
    if n <= 0 {
        return f64::INFINITY;
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:323-329
    blast_ln_gamma_int(n)
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
    -blast_log1p(-p)
}

/// Convert E-value to P-value.
///
/// P = 1 - exp(-E)
pub fn e_to_p(e: f64) -> f64 {
    if e < 0.0 {
        return 0.0;
    }
    // NCBI: -BLAST_Expm1(-e)
    -blast_expm1(-e)
}

/// Exact ln(n!) for integers (n >= 0).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:473-480
fn ln_factorial_int(n: i32) -> f64 {
    blast_ln_factorial(n as f64)
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
        return -blast_expm1(-(-s).exp());
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
      return -blast_expm1(-(-s).exp());
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
