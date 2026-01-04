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

/// Parameters for reward=1, penalty=-5
const BLASTN_1_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.39, 0.747, 1.38, 1.00, 0.0),
    ParamEntry::new(3, 3, 1.39, 0.747, 1.38, 1.00, 0.0),
];

/// Parameters for reward=1, penalty=-4
const BLASTN_1_4: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.383, 0.738, 1.36, 1.02, 0.0),
    ParamEntry::new(1, 2, 1.36, 0.67, 1.2, 1.1, 0.0),
    ParamEntry::new(0, 2, 1.26, 0.43, 0.90, 1.4, -1.0),
    ParamEntry::new(2, 1, 1.35, 0.61, 1.1, 1.2, -1.0),
    ParamEntry::new(1, 1, 1.22, 0.35, 0.72, 1.7, -3.0),
];

/// Parameters for reward=2, penalty=-7 (even scores only)
const BLASTN_2_7: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.69, 0.73, 1.34, 0.515, 0.0),
    ParamEntry::new(2, 4, 0.68, 0.67, 1.2, 0.55, 0.0),
    ParamEntry::new(0, 4, 0.63, 0.43, 0.90, 0.7, -1.0),
    ParamEntry::new(4, 2, 0.675, 0.62, 1.1, 0.6, -1.0),
    ParamEntry::new(2, 2, 0.61, 0.35, 0.72, 1.7, -3.0),
];

/// Parameters for reward=1, penalty=-3
const BLASTN_1_3: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.374, 0.711, 1.31, 1.05, 0.0),
    ParamEntry::new(2, 2, 1.37, 0.70, 1.2, 1.1, 0.0),
    ParamEntry::new(1, 2, 1.35, 0.64, 1.1, 1.2, -1.0),
    ParamEntry::new(0, 2, 1.25, 0.42, 0.83, 1.5, -2.0),
    ParamEntry::new(2, 1, 1.34, 0.60, 1.1, 1.2, -1.0),
    ParamEntry::new(1, 1, 1.21, 0.34, 0.71, 1.7, -2.0),
];

/// Parameters for reward=2, penalty=-5 (even scores only)
const BLASTN_2_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.675, 0.65, 1.1, 0.6, -1.0),
    ParamEntry::new(2, 4, 0.67, 0.59, 1.1, 0.6, -1.0),
    ParamEntry::new(0, 4, 0.62, 0.39, 0.78, 0.8, -2.0),
    ParamEntry::new(4, 2, 0.67, 0.61, 1.0, 0.65, -2.0),
    ParamEntry::new(2, 2, 0.56, 0.32, 0.59, 0.95, -4.0),
];

/// Parameters for reward=1, penalty=-2 (BLASTN default)
const BLASTN_1_2: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.28, 0.46, 0.85, 1.5, -2.0),
    ParamEntry::new(2, 2, 1.33, 0.62, 1.1, 1.2, 0.0),
    ParamEntry::new(1, 2, 1.30, 0.52, 0.93, 1.4, -2.0),
    ParamEntry::new(0, 2, 1.19, 0.34, 0.66, 1.8, -3.0),
    ParamEntry::new(3, 1, 1.32, 0.57, 1.0, 1.3, -1.0),
    ParamEntry::new(2, 1, 1.29, 0.49, 0.92, 1.4, -1.0),
    ParamEntry::new(1, 1, 1.14, 0.26, 0.52, 2.2, -5.0),
];

/// Parameters for reward=2, penalty=-3 (even scores only)
/// Note: NCBI blastn task default is gap_open=5, gap_extend=2
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
const BLASTN_3_4: &[ParamEntry] = &[
    ParamEntry::new(6, 3, 0.389, 0.25, 0.56, 0.7, -5.0),
    ParamEntry::new(5, 3, 0.375, 0.21, 0.47, 0.8, -6.0),
    ParamEntry::new(4, 3, 0.351, 0.14, 0.35, 1.0, -9.0),
    ParamEntry::new(6, 2, 0.362, 0.16, 0.45, 0.8, -4.0),
    ParamEntry::new(5, 2, 0.330, 0.092, 0.28, 1.2, -13.0),
    ParamEntry::new(4, 2, 0.281, 0.046, 0.16, 1.8, -23.0),
];

/// Parameters for reward=4, penalty=-5
const BLASTN_4_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.22, 0.061, 0.22, 1.0, -15.0),
    ParamEntry::new(6, 5, 0.28, 0.21, 0.47, 0.6, -7.0),
    ParamEntry::new(5, 5, 0.27, 0.17, 0.39, 0.7, -9.0),
    ParamEntry::new(4, 5, 0.25, 0.10, 0.31, 0.8, -10.0),
    ParamEntry::new(3, 5, 0.23, 0.065, 0.25, 0.9, -11.0),
];

/// Parameters for reward=1, penalty=-1
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
const BLASTN_3_2: &[ParamEntry] = &[ParamEntry::new(5, 5, 0.208, 0.030, 0.072, 2.9, -47.0)];

/// Parameters for reward=5, penalty=-4
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
            // Default to reward=1, penalty=-2 parameters
            return KarlinParams {
                lambda: 1.28,
                k: 0.46,
                h: 0.85,
                alpha: 1.5,
                beta: -2.0,
            };
        }
    };

    // Find matching gap penalties or use ungapped (first entry)
    for entry in table {
        if entry.gap_open == gap_open && entry.gap_extend == gap_extend {
            return entry.to_karlin_params();
        }
    }

    // If no exact match, use ungapped parameters (first entry)
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

/// Check if scores need to be rounded down for even-score-only matrices.
///
/// For certain reward/penalty combinations, NCBI BLAST requires that odd scores
/// be rounded down to the nearest even number before calculating E-values.
/// This matches the behavior in NCBI BLAST's `s_GetNuclValuesArray` function.
pub fn requires_even_scores(reward: i32, penalty: i32) -> bool {
    let penalty = penalty.abs();
    // NCBI BLAST sets round_down=TRUE for these combinations
    matches!((reward, penalty), (2, 7) | (2, 5) | (2, 3) | (3, 4))
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

    #[test]
    fn test_lookup_nucl_params_default() {
        let spec = NuclScoringSpec {
            reward: 1,
            penalty: -2,
            gap_open: 0,
            gap_extend: 0,
        };
        let params = lookup_nucl_params(&spec);
        assert!((params.lambda - 1.28).abs() < 0.01);
        assert!((params.k - 0.46).abs() < 0.01);
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
