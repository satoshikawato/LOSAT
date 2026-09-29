use crate::config::{NuclScoringSpec, ProteinScoringSpec, ScoringMatrix};

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

/// Entry in the statistical parameter table
/// Format: (gap_open, gap_extend, lambda, k, h, alpha, beta)
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
    const fn new(
        gap_open: i32,
        gap_extend: i32,
        lambda: f64,
        k: f64,
        h: f64,
        alpha: f64,
        beta: f64,
    ) -> Self {
        Self {
            gap_open,
            gap_extend,
            lambda,
            k,
            h,
            alpha,
            beta,
        }
    }

    fn to_karlin_params(&self) -> KarlinParams {
        KarlinParams {
            lambda: self.lambda,
            k: self.k,
            h: self.h,
            alpha: self.alpha,
            beta: self.beta,
        }
    }
}

// ============================================================================
// NUCLEOTIDE STATISTICAL PARAMETERS (from NCBI blast_stat.c)
// ============================================================================
// NCBI reference: blast_stat.c:611-614
// Format: { gap_open, gap_extend, lambda, k, h, alpha, beta, ... }
// LOSAT uses: ParamEntry::new(gap_open, gap_extend, lambda, k, h, alpha, beta)

/// Parameters for reward=1, penalty=-5
/// NCBI reference: blast_stat.c:611-614 (blastn_values_1_5)
const BLASTN_1_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.39, 0.747, 1.38, 1.00, 0.0),
    ParamEntry::new(3, 3, 1.39, 0.747, 1.38, 1.00, 0.0),
];

/// Parameters for reward=1, penalty=-4
/// NCBI reference: blast_stat.c:617-623 (blastn_values_1_4)
const BLASTN_1_4: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.383, 0.738, 1.36, 1.02, 0.0),
    ParamEntry::new(1, 2, 1.36, 0.67, 1.2, 1.1, 0.0),
    ParamEntry::new(0, 2, 1.26, 0.43, 0.90, 1.4, -1.0),
    ParamEntry::new(2, 1, 1.35, 0.61, 1.1, 1.2, -1.0),
    ParamEntry::new(1, 1, 1.22, 0.35, 0.72, 1.7, -3.0),
];

/// Parameters for reward=2, penalty=-7 (even scores only)
/// NCBI reference: blast_stat.c:629-635 (blastn_values_2_7)
/// Note: These parameters can only be applied to even scores. Any odd score must be
/// rounded down to the nearest even number before calculating the e-value.
const BLASTN_2_7: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.69, 0.73, 1.34, 0.515, 0.0),
    ParamEntry::new(2, 4, 0.68, 0.67, 1.2, 0.55, 0.0),
    ParamEntry::new(0, 4, 0.63, 0.43, 0.90, 0.7, -1.0),
    ParamEntry::new(4, 2, 0.675, 0.62, 1.1, 0.6, -1.0),
    ParamEntry::new(2, 2, 0.61, 0.35, 0.72, 1.7, -3.0),
];

/// Parameters for reward=1, penalty=-3
/// NCBI reference: blast_stat.c:638-645 (blastn_values_1_3)
const BLASTN_1_3: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.374, 0.711, 1.31, 1.05, 0.0),
    ParamEntry::new(2, 2, 1.37, 0.70, 1.2, 1.1, 0.0),
    ParamEntry::new(1, 2, 1.35, 0.64, 1.1, 1.2, -1.0),
    ParamEntry::new(0, 2, 1.25, 0.42, 0.83, 1.5, -2.0),
    ParamEntry::new(2, 1, 1.34, 0.60, 1.1, 1.2, -1.0),
    ParamEntry::new(1, 1, 1.21, 0.34, 0.71, 1.7, -2.0),
];

/// Parameters for reward=2, penalty=-5 (even scores only)
/// NCBI reference: blast_stat.c:651-657 (blastn_values_2_5)
/// Note: These parameters can only be applied to even scores. Any odd score must be
/// rounded down to the nearest even number before calculating the e-value.
const BLASTN_2_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.675, 0.65, 1.1, 0.6, -1.0),
    ParamEntry::new(2, 4, 0.67, 0.59, 1.1, 0.6, -1.0),
    ParamEntry::new(0, 4, 0.62, 0.39, 0.78, 0.8, -2.0),
    ParamEntry::new(4, 2, 0.67, 0.61, 1.0, 0.65, -2.0),
    ParamEntry::new(2, 2, 0.56, 0.32, 0.59, 0.95, -4.0),
];

/// Parameters for reward=1, penalty=-2 (megablast task default)
/// NCBI reference: blast_stat.c:660-668 (blastn_values_1_2)
const BLASTN_1_2: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.28, 0.46, 0.85, 1.5, -2.0),
    ParamEntry::new(2, 2, 1.33, 0.62, 1.1, 1.2, 0.0),
    ParamEntry::new(1, 2, 1.30, 0.52, 0.93, 1.4, -2.0),
    ParamEntry::new(0, 2, 1.19, 0.34, 0.66, 1.8, -3.0),
    ParamEntry::new(3, 1, 1.32, 0.57, 1.0, 1.3, -1.0),
    ParamEntry::new(2, 1, 1.29, 0.49, 0.92, 1.4, -1.0),
    ParamEntry::new(1, 1, 1.14, 0.26, 0.52, 2.2, -5.0),
];

/// Parameters for reward=2, penalty=-3 (blastn task default)
/// NCBI reference: blast_stat.c:674-684 (blastn_values_2_3)
/// Note: These parameters can only be applied to even scores. Any odd score must be
/// rounded down to the nearest even number before calculating the e-value.
/// NCBI blastn task default: gap_open=5, gap_extend=2 → ParamEntry::new(5, 2, 0.625, 0.41, 0.78, 0.8, -2.0)
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

/// Parameters for reward=3, penalty=-4
/// NCBI reference: blast_stat.c:687-694 (blastn_values_3_4)
/// Note: These parameters can only be applied to even scores. Any odd score must be
/// rounded down to the nearest even number before calculating the e-value.
const BLASTN_3_4: &[ParamEntry] = &[
    ParamEntry::new(6, 3, 0.389, 0.25, 0.56, 0.7, -5.0),
    ParamEntry::new(5, 3, 0.375, 0.21, 0.47, 0.8, -6.0),
    ParamEntry::new(4, 3, 0.351, 0.14, 0.35, 1.0, -9.0),
    ParamEntry::new(6, 2, 0.362, 0.16, 0.45, 0.8, -4.0),
    ParamEntry::new(5, 2, 0.330, 0.092, 0.28, 1.2, -13.0),
    ParamEntry::new(4, 2, 0.281, 0.046, 0.16, 1.8, -23.0),
];

/// Parameters for reward=4, penalty=-5
/// NCBI reference: blast_stat.c:697-703 (blastn_values_4_5)
const BLASTN_4_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.22, 0.061, 0.22, 1.0, -15.0),
    ParamEntry::new(6, 5, 0.28, 0.21, 0.47, 0.6, -7.0),
    ParamEntry::new(5, 5, 0.27, 0.17, 0.39, 0.7, -9.0),
    ParamEntry::new(4, 5, 0.25, 0.10, 0.31, 0.8, -10.0),
    ParamEntry::new(3, 5, 0.23, 0.065, 0.25, 0.9, -11.0),
];

/// Parameters for reward=1, penalty=-1
/// NCBI reference: blast_stat.c:706-714 (blastn_values_1_1)
const BLASTN_1_1: &[ParamEntry] = &[
    ParamEntry::new(3, 2, 1.09, 0.31, 0.55, 2.0, -2.0),
    ParamEntry::new(2, 2, 1.07, 0.27, 0.49, 2.2, -3.0),
    ParamEntry::new(1, 2, 1.02, 0.21, 0.36, 2.8, -6.0),
    ParamEntry::new(0, 2, 0.80, 0.064, 0.17, 4.8, -16.0),
    ParamEntry::new(4, 1, 1.08, 0.28, 0.54, 2.0, -2.0),
    ParamEntry::new(3, 1, 1.06, 0.25, 0.46, 2.3, -4.0),
    ParamEntry::new(2, 1, 0.99, 0.17, 0.30, 3.3, -10.0),
];

/// Parameters for reward=3, penalty=-2
/// NCBI reference: blast_stat.c:717-719 (blastn_values_3_2)
const BLASTN_3_2: &[ParamEntry] = &[ParamEntry::new(5, 5, 0.208, 0.030, 0.072, 2.9, -47.0)];

/// Parameters for reward=5, penalty=-4
/// NCBI reference: blast_stat.c:722-724 (blastn_values_5_4)
const BLASTN_5_4: &[ParamEntry] = &[
    ParamEntry::new(10, 6, 0.163, 0.068, 0.16, 1.0, -19.0),
    ParamEntry::new(8, 6, 0.146, 0.039, 0.11, 1.3, -29.0),
];

// ============================================================================
// PROTEIN STATISTICAL PARAMETERS (from NCBI blast_stat.c)
// ============================================================================

/// BLOSUM45 parameters
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

/// BLOSUM50 parameters
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

/// BLOSUM62 parameters (default for protein)
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

/// BLOSUM80 parameters
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

/// BLOSUM90 parameters
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

/// PAM30 parameters
const PAM30: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3400, 0.283, 1.754, 0.1938, -0.3),
    ParamEntry::new(7, 2, 0.305, 0.15, 0.87, 0.35, -3.0),
    ParamEntry::new(6, 2, 0.287, 0.11, 0.68, 0.42, -4.0),
    ParamEntry::new(5, 2, 0.264, 0.079, 0.45, 0.59, -7.0),
    ParamEntry::new(10, 1, 0.309, 0.15, 0.88, 0.34, -3.0),
    ParamEntry::new(9, 1, 0.294, 0.11, 0.61, 0.48, -6.0),
    ParamEntry::new(8, 1, 0.270, 0.072, 0.40, 0.68, -10.0),
];

/// PAM70 parameters
const PAM70: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3345, 0.229, 1.029, 0.3250, -0.9),
    ParamEntry::new(8, 2, 0.301, 0.12, 0.54, 0.56, -5.0),
    ParamEntry::new(7, 2, 0.286, 0.093, 0.43, 0.67, -7.0),
    ParamEntry::new(6, 2, 0.264, 0.064, 0.29, 0.90, -12.0),
    ParamEntry::new(11, 1, 0.305, 0.12, 0.52, 0.59, -6.0),
    ParamEntry::new(10, 1, 0.291, 0.091, 0.41, 0.71, -9.0),
    ParamEntry::new(9, 1, 0.270, 0.060, 0.28, 0.97, -14.0),
];

/// PAM250 parameters
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

/// NCBI reference: c++/src/algo/blast/core/ncbi_math.c:405-419
/// ```c
/// Int4 BLAST_Gcd(Int4 a, Int4 b)
/// {
///    Int4   c;
///
///    b = ABS(b);
///    if (b > a)
///       c=a, a=b, b=c;
///
///    while (b != 0) {
///       c = a%b;
///       a = b;
///       b = c;
///    }
///    return a;
/// }
/// ```
fn blast_gcd(a: i32, b: i32) -> i32 {
    let (mut a, mut b) = (a, b.abs());
    if b > a {
        std::mem::swap(&mut a, &mut b);
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a
}

/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3955-3963
/// ```c
/// static double s_GetUngappedBeta(Int4 reward, Int4 penalty)
/// {
///     double beta = 0;
///     if ((reward == 1 && penalty == -1) ||
///         (reward == 2 && penalty == -3))
///         beta = -2;
///
///     return beta;
/// }
/// ```
fn ungapped_beta(reward: i32, penalty: i32) -> f64 {
    if (reward, penalty) == (1, -1) || (reward, penalty) == (2, -3) {
        -2.0
    } else {
        0.0
    }
}

/// The Karlin-Altschul values of a BLASTN reward/penalty pair (NCBI
/// `s_GetNuclValuesArray`): the table of the pair divided by its greatest common divisor,
/// with the gap costs multiplied and Lambda and alpha divided by the divisor.
///
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3237-3372
/// ```c
///     int divisor = BLAST_Gcd(reward, penalty);
///
///     *round_down = FALSE;
///     ...
///     if (divisor != 1)
///     {
///        reward /= divisor;
///        penalty /= divisor;
///     }
///
///     if (reward == 1 && penalty == -5) {
///         if ((status=s_SplitArrayOf8(blastn_values_1_5, &kValues, &kValues_non_affine, &split)))
///            return status;
///
///         *array_size = sizeof(blastn_values_1_5)/sizeof(array_of_8);
///         *gap_open_max = 3;
///         *gap_extend_max = 3;
///     ...
///     } else if (reward == 2 && penalty == -3) {
///         ...
///         *round_down = TRUE;
///         *array_size = sizeof(blastn_values_2_3)/sizeof(array_of_8);
///         *gap_open_max = 6;
///         *gap_extend_max = 4;
///     ...
///     } else  { /* Unsupported reward-penalty */
///         status = -1;
///         if (error_return) {
///             char buffer[256];
///             snprintf(buffer, sizeof(buffer), "Substitution scores %d and %d are not supported",
///                 reward, penalty);
///             Blast_MessageWrite(error_return, eBlastSevError, kBlastMessageNoContext, buffer);
///         }
///     }
///     if (split)
///         (*array_size)--;
///     ...
///         status = s_AdjustGapParametersByGcd(*normal, *non_affine, *array_size, gap_open_max, gap_extend_max, divisor);
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3159-3170 (s_SplitArrayOf8)
/// ```c
///     if (input[0][0] == 0 && input[0][1] == 0)
///     {
///             *normal = input+1;
///             *non_affine = input;
///             *split = TRUE;
///     }
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3192-3215 (s_AdjustGapParametersByGcd)
/// ```c
///     (*gap_existence_max) *= divisor;
///     (*gap_extend_max) *= divisor;
///     ...
///                 normal[i][0] *= divisor;
///                 normal[i][1] *= divisor;
///                 normal[i][2] /= divisor;
///                 normal[i][5] /= divisor;
///     ...
///        linear[0][0] *= divisor;
///        linear[0][1] *= divisor;
///        linear[0][2] /= divisor;
///        linear[0][5] /= divisor;
/// ```
/// The command line accepts no reward of 0, so the divisor is at least 1.
pub struct NuclValues {
    /// The affine gap costs of the table, with their values.
    normal: Vec<ParamEntry>,
    /// The values of gap costs 0 and 0 (megablast's linear gaps), when the table has them.
    linear: Option<ParamEntry>,
    gap_open_max: i32,
    gap_extend_max: i32,
    /// Scores are rounded down to even numbers for e-values and bit scores.
    pub round_down: bool,
}

impl NuclValues {
    pub fn new(reward: i32, penalty: i32) -> Result<Self, String> {
        let divisor = blast_gcd(reward, penalty).max(1);
        let (reward, penalty) = (reward / divisor, penalty / divisor);
        let (table, gap_open_max, gap_extend_max, round_down) = match (reward, penalty) {
            (1, -5) => (BLASTN_1_5, 3, 3, false),
            (1, -4) => (BLASTN_1_4, 2, 2, false),
            (2, -7) => (BLASTN_2_7, 4, 4, true),
            (1, -3) => (BLASTN_1_3, 2, 2, false),
            (2, -5) => (BLASTN_2_5, 4, 4, true),
            (1, -2) => (BLASTN_1_2, 2, 2, false),
            (2, -3) => (BLASTN_2_3, 6, 4, true),
            (3, -4) => (BLASTN_3_4, 6, 3, true),
            (1, -1) => (BLASTN_1_1, 4, 2, false),
            (3, -2) => (BLASTN_3_2, 5, 5, false),
            (4, -5) => (BLASTN_4_5, 12, 8, false),
            (5, -4) => (BLASTN_5_4, 25, 10, false),
            _ => {
                return Err(format!(
                    "Substitution scores {reward} and {penalty} are not supported"
                ))
            }
        };
        let adjust = |entry: &ParamEntry| ParamEntry {
            gap_open: entry.gap_open * divisor,
            gap_extend: entry.gap_extend * divisor,
            lambda: entry.lambda / divisor as f64,
            alpha: entry.alpha / divisor as f64,
            ..*entry
        };
        let split = table[0].gap_open == 0 && table[0].gap_extend == 0;
        Ok(Self {
            normal: table[usize::from(split)..].iter().map(adjust).collect(),
            linear: split.then(|| adjust(&table[0])),
            gap_open_max: gap_open_max * divisor,
            gap_extend_max: gap_extend_max * divisor,
            round_down,
        })
    }

    /// The table row of the gap costs: `Some` row, `None` for gap costs beyond the table
    /// (which use the ungapped block), or NCBI's message for unsupported gap costs
    /// (`Blast_KarlinBlkNuclGappedCalc`). The message names the given (not the divided)
    /// scores.
    ///
    /// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3875-3935
    /// ```c
    ///     if (gap_open == 0 && gap_extend == 0 && linear)
    ///     {
    ///         kbp->Lambda = linear[0][kLambdaIndex];
    ///     ...
    ///         for (index = 0; index < num_combinations; ++index) {
    ///             if (normal[index][kGapOpenIndex] == gap_open &&
    ///                 normal[index][kGapExtIndex] == gap_extend) {
    ///     ...
    ///         if (index == num_combinations) {
    ///         /* If gap costs are larger than maximal provided in tables, copy
    ///            the values from the ungapped Karlin block. */
    ///             if (gap_open >= gap_open_max && gap_extend >= gap_extend_max) {
    ///                 Blast_KarlinBlkCopy(kbp, kbp_ungap);
    ///             } else if (error_return) {
    ///     ...
    ///                 out_sz = snprintf(buffer, (size_t)buffer_sz,"Gap existence and extension values %ld and %ld "
    ///                         "are not supported for substitution scores %ld and %ld\n",
    ///                         (long) gap_open, (long) gap_extend, (long) reward, (long) penalty);
    ///     ...
    ///                      out_sz = snprintf(buffer+len, (size_t)buffer_sz, "%ld and %ld are supported existence and extension values\n",
    ///                         (long) normal[i][kGapOpenIndex],  (long) normal[i][kGapExtIndex]);
    ///     ...
    ///                 out_sz = snprintf(buffer+len, (size_t)buffer_sz, "%ld and %ld are supported existence and extension values\n",
    ///                      (long) gap_open_max, (long) gap_extend_max);
    ///     ...
    ///                 out_sz = snprintf(buffer+len, (size_t)buffer_sz, "Any values more stringent than %ld and %ld are supported\n",
    ///                      (long) gap_open_max, (long) gap_extend_max);
    /// ```
    fn row(&self, spec: &NuclScoringSpec) -> Result<Option<ParamEntry>, String> {
        let row = if spec.gap_open == 0 && spec.gap_extend == 0 && self.linear.is_some() {
            self.linear
        } else {
            self.normal
                .iter()
                .find(|entry| {
                    entry.gap_open == spec.gap_open && entry.gap_extend == spec.gap_extend
                })
                .copied()
        };
        if row.is_some()
            || (spec.gap_open >= self.gap_open_max && spec.gap_extend >= self.gap_extend_max)
        {
            return Ok(row);
        }
        let mut message = format!(
            "Gap existence and extension values {} and {} are not supported for substitution scores {} and {}\n",
            spec.gap_open, spec.gap_extend, spec.reward, spec.penalty
        );
        let supported = self
            .normal
            .iter()
            .map(|entry| (entry.gap_open, entry.gap_extend))
            .chain([(self.gap_open_max, self.gap_extend_max)]);
        for (gap_open, gap_extend) in supported {
            message += &format!(
                "{gap_open} and {gap_extend} are supported existence and extension values\n"
            );
        }
        message += &format!(
            "Any values more stringent than {} and {} are supported\n",
            self.gap_open_max, self.gap_extend_max
        );
        Err(message)
    }

    /// NCBI's message when the gap costs are not supported for these scores.
    pub fn check_gaps(&self, spec: &NuclScoringSpec) -> Result<(), String> {
        self.row(spec).map(|_| ())
    }

    /// The gapped Karlin block of a query context, with the alpha and beta of its length
    /// adjustment (NCBI `Blast_KarlinBlkNuclGappedCalc` and `Blast_GetNuclAlphaBeta`).
    /// `ungapped` is the context's ungapped block, which gap costs beyond the table use.
    ///
    /// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3995-4026
    /// ```c
    ///     if (gapped_calculation && normal) {
    ///         if (gap_open == 0 && gap_extend == 0 && linear)
    ///         {
    ///             *alpha = linear[0][kAlphaIndex];
    ///             *beta = linear[0][kBetaIndex];
    ///     ...
    ///     if (!found)
    ///     {
    ///         *alpha = kbp->Lambda/kbp->H;
    ///         *beta = s_GetUngappedBeta(reward, penalty);
    ///     }
    /// ```
    pub fn gapped(
        &self,
        spec: &NuclScoringSpec,
        ungapped: &KarlinParams,
    ) -> Result<KarlinParams, String> {
        Ok(match self.row(spec)? {
            Some(row) => row.to_karlin_params(),
            None => KarlinParams {
                alpha: ungapped.lambda / ungapped.h,
                beta: ungapped_beta(spec.reward, spec.penalty),
                ..*ungapped
            },
        })
    }
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

    // Find matching gap penalties
    for entry in table {
        if entry.gap_open == gap_open && entry.gap_extend == gap_extend {
            return entry.to_karlin_params();
        }
    }

    // If no exact match, try to find closest match
    // First, try ungapped (first entry with MAX values)
    if !table.is_empty() {
        // Look for the "best" entry (typically gap_open=11, gap_extend=1 for BLOSUM62)
        for entry in table {
            if entry.gap_open != i32::MAX {
                return entry.to_karlin_params();
            }
        }
        return table[0].to_karlin_params();
    }

    // Default BLOSUM62 with gap_open=11, gap_extend=1
    KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_options.c:908-936
// ```c
// if ((status=Blast_KarlinBlkGappedLoadFromTables(NULL, options->gap_open,
//       options->gap_extend, options->matrix, std_matrix_only)) != 0)
// {
//     ...
//     return BLASTERR_OPTION_VALUE_INVALID;
// }
// ```
pub fn protein_scoring_supported(spec: &ProteinScoringSpec) -> bool {
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

    table
        .iter()
        .any(|entry| entry.gap_open == spec.gap_open && entry.gap_extend == spec.gap_extend)
}

/// Look up UNGAPPED Karlin-Altschul parameters for protein scoring scheme.
///
/// NCBI BLAST uses ungapped params (kbp_std) for gap_trigger calculation.
/// These are stored as the first entry in each matrix table with gap_open=MAX, gap_extend=MAX.
///
/// Reference: ncbi-blast blast_parameters.c:340-345
/// ```c
/// if (sbp->kbp_std) {
///    kbp = sbp->kbp_std[context];
///    gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda);
/// }
/// ```
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

    // Ungapped params are stored with gap_open=MAX, gap_extend=MAX (first entry)
    for entry in table {
        if entry.gap_open == i32::MAX && entry.gap_extend == i32::MAX {
            return entry.to_karlin_params();
        }
    }

    // Fallback: first entry should be ungapped
    if !table.is_empty() {
        return table[0].to_karlin_params();
    }

    // Default BLOSUM62 ungapped
    KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.4012,
        alpha: 0.7916,
        beta: -3.2,
    }
}

/// Look up GAPPED Karlin-Altschul parameters for protein scoring scheme.
///
/// NCBI BLAST uses gapped params (kbp_gap_std with gap_open=11, gap_extend=1 for BLOSUM62)
/// for length_adjustment calculation in BLAST_CalcEffLengths.
///
/// Reference: ncbi-blast blast_setup.c:814-819
/// ```c
/// BLAST_GetAlphaBeta(sbp->name, &alpha, &beta,
///                    scoring_options->gapped_calculation,  // TRUE
///                    gap_open, gap_extend, sbp->kbp_std[index]);
/// ```
///
/// For tblastx, scoring_options->gapped_calculation = TRUE (default),
/// so BLAST_GetAlphaBeta returns gapped alpha/beta values.
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

    // Gapped params: gap_open=11, gap_extend=1 (BLOSUM62 default)
    for entry in table {
        if entry.gap_open == 11 && entry.gap_extend == 1 {
            return entry.to_karlin_params();
        }
    }

    // Fallback: first non-ungapped entry
    for entry in table {
        if entry.gap_open != i32::MAX {
            return entry.to_karlin_params();
        }
    }

    // Default BLOSUM62 gapped (11/1)
    KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(reward: i32, penalty: i32, gap_open: i32, gap_extend: i32) -> NuclScoringSpec {
        NuclScoringSpec {
            reward,
            penalty,
            gap_open,
            gap_extend,
        }
    }

    const UNGAPPED: KarlinParams = KarlinParams {
        lambda: 1.1,
        k: 0.3,
        h: 0.5,
        alpha: 0.0,
        beta: 0.0,
    };

    fn values(spec: &NuclScoringSpec) -> Result<(KarlinParams, bool), String> {
        let values = NuclValues::new(spec.reward, spec.penalty)?;
        Ok((values.gapped(spec, &UNGAPPED)?, values.round_down))
    }

    #[test]
    fn task_defaults_use_the_table_rows() {
        // blastn_values_1_2[0] (linear) and blastn_values_2_3 row { 5, 2, ... }.
        let (megablast, round_down) = values(&spec(1, -2, 0, 0)).unwrap();
        assert_eq!(
            (
                megablast.lambda,
                megablast.k,
                megablast.h,
                megablast.alpha,
                megablast.beta
            ),
            (1.28, 0.46, 0.85, 1.5, -2.0)
        );
        assert!(!round_down);
        let (blastn, round_down) = values(&spec(2, -3, 5, 2)).unwrap();
        assert_eq!(
            (blastn.lambda, blastn.k, blastn.h, blastn.alpha, blastn.beta),
            (0.625, 0.41, 0.78, 0.8, -2.0)
        );
        assert!(round_down);
    }

    #[test]
    fn scores_with_a_common_divisor_use_the_divided_table() {
        // 2/-4 is 1/-2: the gap costs double, Lambda and alpha halve.
        let (params, round_down) = values(&spec(2, -4, 4, 4)).unwrap();
        assert_eq!(
            (params.lambda, params.k, params.h, params.alpha, params.beta),
            (1.33 / 2.0, 0.62, 1.1, 1.2 / 2.0, 0.0)
        );
        assert!(!round_down);
        // 4/-6 is 2/-3, whose scores are rounded down.
        assert!(values(&spec(4, -6, 10, 4)).unwrap().1);
        assert!(values(&spec(2, -4, 3, 3))
            .unwrap_err()
            .contains("\n4 and 4 are supported"));
    }

    #[test]
    fn gap_costs_beyond_the_table_use_the_ungapped_block() {
        let (params, _) = values(&spec(1, -2, 5, 2)).unwrap();
        assert_eq!(
            (params.lambda, params.k, params.h, params.alpha, params.beta),
            (1.1, 0.3, 0.5, 1.1 / 0.5, 0.0)
        );
        assert_eq!(values(&spec(2, -3, 7, 4)).unwrap().0.beta, -2.0);
    }

    #[test]
    fn unsupported_scores_have_ncbi_messages() {
        assert_eq!(
            values(&spec(2, -5, 5, 2)).unwrap_err(),
            "Gap existence and extension values 5 and 2 are not supported for substitution scores 2 and -5\n\
             2 and 4 are supported existence and extension values\n\
             0 and 4 are supported existence and extension values\n\
             4 and 2 are supported existence and extension values\n\
             2 and 2 are supported existence and extension values\n\
             4 and 4 are supported existence and extension values\n\
             Any values more stringent than 4 and 4 are supported\n"
        );
        assert_eq!(
            values(&spec(1, -6, 0, 0)).unwrap_err(),
            "Substitution scores 1 and -6 are not supported"
        );
        // The message names the divided scores; 5/-4 has no linear row.
        assert_eq!(
            values(&spec(2, -12, 0, 0)).unwrap_err(),
            "Substitution scores 1 and -6 are not supported"
        );
        assert!(values(&spec(5, -4, 0, 0))
            .unwrap_err()
            .starts_with("Gap existence and extension values 0 and 0 are not supported"));
    }

    #[test]
    fn test_lookup_protein_params_blosum62() {
        let spec = ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        };
        let params = lookup_protein_params(&spec);
        assert!((params.lambda - 0.267).abs() < 0.01);
    }
}
