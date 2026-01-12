//! NCBI-style Karlin-Altschul parameter calculation from amino acid composition
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c
//!   - Blast_ResFreqString: compute amino acid composition
//!   - BlastScoreFreqCalc: compute score frequency profile
//!   - Blast_KarlinBlkUngappedCalc: compute Lambda, K, H
//!   - check_ideal logic: use kbp_ideal if computed Lambda >= ideal Lambda

use crate::stats::KarlinParams;
use crate::utils::matrix::{blosum62_score, ncbistdaa, BLASTAA_SIZE};

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_stat.h:121-122
const BLAST_SCORE_MIN: i32 = i16::MIN as i32;
const BLAST_SCORE_MAX: i32 = i16::MAX as i32;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:56-72
const BLAST_SCORE_RANGE_MAX: i32 = BLAST_SCORE_MAX - BLAST_SCORE_MIN;
const BLAST_KARLIN_K_SUMLIMIT_DEFAULT: f64 = 0.0001;
const BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT: f64 = 1.0e-5;
const BLAST_KARLIN_LAMBDA_ITER_DEFAULT: i32 = 17;
const BLAST_KARLIN_LAMBDA0_DEFAULT: f64 = 0.5;
const BLAST_KARLIN_K_ITER_MAX: i32 = 100;

/// Standard amino acid frequencies (Robinson probabilities)
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1818
///   STD_AMINO_ACID_FREQS Robinson_prob
/// These are the standard frequencies used for database composition
const STD_AA_FREQS: [f64; 20] = [
    0.07805, // A
    0.01926, // R
    0.05364, // N
    0.06295, // D
    0.01487, // C
    0.03374, // Q
    0.06661, // E
    0.07129, // G
    0.02105, // H
    0.05142, // I
    0.05744, // L
    0.05068, // K
    0.01471, // M
    0.03965, // F
    0.04728, // P
    0.06141, // S
    0.05506, // T
    0.01330, // W
    0.03216, // Y
    0.06891, // V
];

/// NCBISTDAA to standard amino acid mapping (for standard composition)
/// Standard composition uses 20 standard amino acids (A-R-N-D-C-Q-E-G-H-I-L-K-M-F-P-S-T-W-Y-V)
/// NCBISTDAA: 0='-', 1=A, 2=B, 3=C, 4=D, 5=E, 6=F, 7=G, 8=H, 9=I, 10=K, 11=L, 12=M, 13=N, 14=P, 15=Q, 16=R, 17=S, 18=T, 19=V, 20=W, 21=X, 22=Y, 23=Z, 24=U, 25='*', 26=O, 27=J
/// Standard AA order: A(0), R(1), N(2), D(3), C(4), Q(5), E(6), G(7), H(8), I(9), L(10), K(11), M(12), F(13), P(14), S(15), T(16), W(17), Y(18), V(19)
const NCBISTDAA_TO_STD_AA: [usize; 28] = [
    20, // 0: '-' -> invalid
    0,  // 1: A -> 0
    2,  // 2: B (N or D) -> use N (2)
    4,  // 3: C -> 4
    3,  // 4: D -> 3
    6,  // 5: E -> 6
    13, // 6: F -> 13
    7,  // 7: G -> 7
    8,  // 8: H -> 8
    9,  // 9: I -> 9
    11, // 10: K -> 11
    10, // 11: L -> 10
    12, // 12: M -> 12
    2,  // 13: N -> 2
    14, // 14: P -> 14
    5,  // 15: Q -> 5
    1,  // 16: R -> 1
    15, // 17: S -> 15
    16, // 18: T -> 16
    19, // 19: V -> 19
    17, // 20: W -> 17
    20, // 21: X (unknown) -> invalid
    18, // 22: Y -> 18
    5,  // 23: Z (E or Q) -> use Q (5)
    20, // 24: U -> invalid
    20, // 25: '*' -> invalid
    20, // 26: O -> invalid
    10, // 27: J (L or I) -> use L (10)
];

/// Compute amino acid composition from sequence
/// Reference: NCBI Blast_ResFreqString (blast_stat.c:2078-2091)
///   - BlastResCompStr: count residues
///   - Blast_ResFreqResComp: normalize to frequencies
///
/// # Arguments
/// * `seq` - Amino acid sequence in NCBISTDAA encoding (0-27)
/// * `length` - Length of sequence (excluding sentinels)
///
/// # Returns
/// Array of 28 frequencies (one per NCBISTDAA residue)
pub fn compute_aa_composition(seq: &[u8], length: usize) -> [f64; BLASTAA_SIZE] {
    let mut comp: [u64; BLASTAA_SIZE] = [0; BLASTAA_SIZE];
    
    // Count residues (NCBI BlastResCompStr)
    // Skip sentinels (first and last byte)
    let start = 1;
    let end = seq.len().saturating_sub(1);
    let actual_len = end.saturating_sub(start).min(length);
    
    for i in start..end.min(start + actual_len) {
        let residue = seq[i] as usize;
        if residue < BLASTAA_SIZE {
            comp[residue] += 1;
        }
    }
    
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1984-2014
    // Ambiguous residues (proteins) are zeroed via Blast_ResCompStr; for BLASTAA
    // this is just 'X' per blast_setup.c:382-385 (BLAST_ScoreSetAmbigRes).
    comp[ncbistdaa::X as usize] = 0;

    // Normalize to frequencies (NCBI Blast_ResFreqResComp)
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
/// Reference: NCBI Blast_ResFreqStdComp (blast_stat.c:1887-1917)
/// Uses Robinson standard frequencies
pub fn compute_std_aa_composition() -> [f64; BLASTAA_SIZE] {
    let mut freq = [0.0; BLASTAA_SIZE];
    
    // Map standard frequencies to NCBISTDAA positions
    for (ncbi_idx, &std_idx) in NCBISTDAA_TO_STD_AA.iter().enumerate() {
        if std_idx < 20 {
            freq[ncbi_idx] = STD_AA_FREQS[std_idx];
        }
    }
    
    // Normalize (should already be normalized, but ensure)
    let sum: f64 = freq.iter().sum();
    if sum > 0.0 {
        for f in freq.iter_mut() {
            *f /= sum;
        }
    }
    
    freq
}

/// Score frequency profile
/// Stores probability distribution of alignment scores
#[derive(Debug, Clone)]
pub struct ScoreFreqProfile {
    /// Score probabilities: sprob[score] = probability of score
    /// Indexed by score value (negative scores are valid)
    sprob: Vec<f64>,
    /// Minimum score with non-zero probability
    obs_min: i32,
    /// Maximum score with non-zero probability
    obs_max: i32,
    /// Average score (must be negative for valid Karlin params)
    score_avg: f64,
    /// Score range offset (sprob is centered around 0)
    score_min: i32,
}

impl ScoreFreqProfile {
    pub fn new(score_min: i32, score_max: i32) -> Self {
        let range = (score_max - score_min + 1) as usize;
        let mut sprob = vec![0.0; range];
        
        Self {
            sprob,
            obs_min: 0,
            obs_max: 0,
            score_avg: 0.0,
            score_min,
        }
    }
    
    pub fn get_prob(&self, score: i32) -> f64 {
        let idx = (score - self.score_min) as usize;
        if idx < self.sprob.len() {
            self.sprob[idx]
        } else {
            0.0
        }
    }
    
    pub fn set_prob(&mut self, score: i32, prob: f64) {
        let idx = (score - self.score_min) as usize;
        if idx < self.sprob.len() {
            self.sprob[idx] = prob;
        }
    }
    
    pub fn obs_min(&self) -> i32 {
        self.obs_min
    }
    
    pub fn obs_max(&self) -> i32 {
        self.obs_max
    }
    
    pub fn score_avg(&self) -> f64 {
        self.score_avg
    }
    
    pub fn sprob(&self) -> &[f64] {
        &self.sprob
    }
}

/// Compute score frequency profile from two amino acid compositions
/// Reference: NCBI BlastScoreFreqCalc (blast_stat.c:2151-2205)
///
/// # Arguments
/// * `comp1` - Query composition (28 frequencies)
/// * `comp2` - Subject/database composition (28 frequencies, typically standard)
/// * `score_min` - Minimum score in matrix
/// * `score_max` - Maximum score in matrix
///
/// # Returns
/// Score frequency profile
pub fn compute_score_freq_profile(
    comp1: &[f64; BLASTAA_SIZE],
    comp2: &[f64; BLASTAA_SIZE],
    score_min: i32,
    score_max: i32,
) -> ScoreFreqProfile {
    let mut sfp = ScoreFreqProfile::new(score_min, score_max);
    
    // Initialize all probabilities to zero
    for score in score_min..=score_max {
        sfp.set_prob(score, 0.0);
    }
    
    // Compute score probabilities: P(score) = sum_{i,j: matrix[i][j]=score} comp1[i] * comp2[j]
    // Reference: blast_stat.c:2171-2181
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
    
    // Find observed min/max and compute average
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2183-2205
    let mut score_sum = 0.0;
    let mut obs_min = BLAST_SCORE_MIN;
    let mut obs_max = BLAST_SCORE_MIN;

    for score in score_min..=score_max {
        let prob = sfp.get_prob(score);
        if prob > 0.0 {
            score_sum += prob;
            obs_max = score;
            if obs_min == BLAST_SCORE_MIN {
                obs_min = score;
            }
        }
    }

    sfp.obs_min = obs_min;
    sfp.obs_max = obs_max;

    // Compute average score (also normalize probabilities)
    let mut score_avg = 0.0;
    if score_sum > 0.0001 || score_sum < -0.0001 {
        for score in obs_min..=obs_max {
            let prob = sfp.get_prob(score) / score_sum;
            sfp.set_prob(score, prob);
            score_avg += (score as f64) * prob;
        }
    }

    sfp.score_avg = score_avg;
    
    sfp
}

/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:405-418
fn blast_gcd(mut a: i32, mut b: i32) -> i32 {
    let mut c: i32;

    b = b.abs();
    if b > a {
        c = a;
        a = b;
        b = c;
    }

    while b != 0 {
        c = a % b;
        a = b;
        b = c;
    }
    a
}

/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2099-2109
fn blast_score_chk(lo: i32, hi: i32) -> Result<(), String> {
    if lo >= 0 || hi <= 0 || lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX {
        return Err("Invalid score range".to_string());
    }
    if hi - lo > BLAST_SCORE_RANGE_MAX {
        return Err("Score range exceeds BLAST_SCORE_RANGE_MAX".to_string());
    }
    Ok(())
}

/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:33-55
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

/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:444-470
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

/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2491-2563
fn nlm_karlin_lambda_nr(
    sfp: &ScoreFreqProfile,
    d: i32,
    low: i32,
    high: i32,
    lambda0: f64,
    tolx: f64,
    itmax: i32,
    max_newton: i32,
) -> Result<f64, String> {
    if d <= 0 {
        return Err("GCD must be positive".to_string());
    }

    let x0 = (-lambda0).exp();
    let mut x = if x0 > 0.0 && x0 < 1.0 { x0 } else { 0.5 };
    let mut a = 0.0;
    let mut b = 1.0;
    let mut f = 4.0;
    let mut is_newton = false;

    let mut k = 0;
    while k < itmax {
        let mut g = 0.0;
        let fold = f;
        let was_newton = is_newton;
        is_newton = false;

        // Horner's rule for evaluating polynomial and derivative
        f = sfp.get_prob(low);
        let mut i = low + d;
        while i < 0 {
            g = x * g + f;
            f = f * x + sfp.get_prob(i);
            i += d;
        }
        g = x * g + f;
        f = f * x + sfp.get_prob(0) - 1.0;
        i = d;
        while i <= high {
            g = x * g + f;
            f = f * x + sfp.get_prob(i);
            i += d;
        }

        if f > 0.0 {
            a = x;
        } else if f < 0.0 {
            b = x;
        } else {
            break;
        }
        if b - a < 2.0 * a * (1.0 - b) * tolx {
            x = (a + b) / 2.0;
            break;
        }

        if k >= max_newton
            || (was_newton && f.abs() > 0.9 * fold.abs())
            || g >= 0.0
        {
            x = (a + b) / 2.0;
        } else {
            let p = -f / g;
            let y = x + p;
            if y <= a || y >= b {
                x = (a + b) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() < tolx * x * (1.0 - x) {
                    break;
                }
            }
        }

        k += 1;
    }

    Ok(-x.ln() / d as f64)
}

/// Compute Lambda using Newton-Raphson method
/// Reference: NCBI Blast_KarlinLambdaNR (blast_stat.c:2567-2598)
/// Uses NlmKarlinLambdaNR internally (blast_stat.c:2491-2563)
fn compute_lambda_nr(sfp: &ScoreFreqProfile, initial_guess: f64) -> Result<f64, String> {
    let low = sfp.obs_min();
    let high = sfp.obs_max();
    
    if sfp.score_avg() >= 0.0 {
        return Err("Expected score must be negative".to_string());
    }

    blast_score_chk(low, high)?;

    // Find greatest common divisor of all scores with non-zero probability
    let mut d = -low;
    let range = high - low;
    for i in 1..=range {
        if d <= 1 {
            break;
        }
        if sfp.get_prob(low + i) != 0.0 {
            d = blast_gcd(d, i);
        }
    }

    nlm_karlin_lambda_nr(
        sfp,
        d,
        low,
        high,
        initial_guess,
        BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
        20,
        20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT,
    )
}

/// Compute H from Lambda
/// Reference: NCBI BlastKarlinLtoH (blast_stat.c:2607-2633)
fn compute_h_from_lambda(sfp: &ScoreFreqProfile, lambda: f64) -> Result<f64, String> {
    if lambda < 0.0 {
        return Err("Lambda must be non-negative".to_string());
    }
    
    let low = sfp.obs_min();
    let high = sfp.obs_max();

    blast_score_chk(low, high)?;

    let etonlam = (-lambda).exp();

    let mut sum = (low as f64) * sfp.get_prob(low);
    for score in (low + 1)..=high {
        sum = (score as f64) * sfp.get_prob(score) + etonlam * sum;
    }
    
    let scale = blast_powi(etonlam, high);
    let h = if scale > 0.0 {
        lambda * sum / scale
    } else {
        // Underflow: use log form
        lambda * (lambda * high as f64 + sum.ln()).exp()
    };
    
    Ok(h)
}

/// Compute K from Lambda and H
/// Reference: NCBI BlastKarlinLHtoK (blast_stat.c:2247-2418)
fn compute_k_from_lambda_h(sfp: &ScoreFreqProfile, lambda: f64, h: f64) -> Result<f64, String> {
    if lambda <= 0.0 || h <= 0.0 {
        return Err("Lambda and H must be positive".to_string());
    }

    if sfp.score_avg() >= 0.0 {
        return Err("Expected score must be negative".to_string());
    }

    let mut low = sfp.obs_min();
    let mut high = sfp.obs_max();
    blast_score_chk(low, high)?;

    let mut range = high - low;
    let mut prob_array_start_low = Vec::with_capacity((range + 1) as usize);
    for i in 0..=range {
        prob_array_start_low.push(sfp.get_prob(low + i));
    }

    let mut divisor = -low;
    for i in 1..=range {
        if divisor <= 1 {
            break;
        }
        if prob_array_start_low[i as usize] != 0.0 {
            divisor = blast_gcd(divisor, i);
        }
    }

    high /= divisor;
    low /= divisor;
    let lambda = lambda * divisor as f64;
    range = high - low;

    let mut first_term_closed_form = h / lambda;
    let exp_minus_lambda = (-lambda).exp();

    if low == -1 && high == 1 {
        let low_prob = sfp.get_prob(low * divisor);
        let high_prob = sfp.get_prob(high * divisor);
        let diff = low_prob - high_prob;
        let k = diff * diff / low_prob;
        return Ok(k);
    }

    if low == -1 || high == 1 {
        if high != 1 {
            let score_avg = sfp.score_avg() / divisor as f64;
            first_term_closed_form = (score_avg * score_avg) / first_term_closed_form;
        }
        return Ok(first_term_closed_form * (1.0 - exp_minus_lambda));
    }

    let sumlimit = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
    let iterlimit = BLAST_KARLIN_K_ITER_MAX;
    let array_len = (iterlimit as usize) * (range as usize) + 1;
    let mut alignment_score_probabilities = vec![0.0; array_len];
    let mut outer_sum = 0.0;
    let mut low_alignment_score = 0;
    let mut high_alignment_score = 0;
    let mut inner_sum = 1.0;
    let mut oldsum = 1.0;
    let mut oldsum2 = 1.0;
    alignment_score_probabilities[0] = 1.0;

    let mut iter_counter = 0;
    while iter_counter < iterlimit && inner_sum > sumlimit {
        let mut first = range;
        let mut last = range;
        low_alignment_score += low;
        high_alignment_score += high;

        let mut ptr_p_idx = (high_alignment_score - low_alignment_score) as isize;
        while ptr_p_idx >= 0 {
            let mut ptr1_idx = ptr_p_idx - first as isize;
            let ptr1e_idx = ptr_p_idx - last as isize;
            let mut ptr2_idx = first as isize;

            inner_sum = 0.0;
            while ptr1_idx >= ptr1e_idx {
                inner_sum += alignment_score_probabilities[ptr1_idx as usize]
                    * prob_array_start_low[ptr2_idx as usize];
                ptr1_idx -= 1;
                ptr2_idx += 1;
            }
            if first > 0 {
                first -= 1;
            }
            if ptr_p_idx <= range as isize {
                last -= 1;
            }
            alignment_score_probabilities[ptr_p_idx as usize] = inner_sum;
            ptr_p_idx -= 1;
        }

        let mut ptr_p_idx = 0isize;
        inner_sum = alignment_score_probabilities[ptr_p_idx as usize];
        let mut i = low_alignment_score + 1;
        while i < 0 {
            ptr_p_idx += 1;
            inner_sum = alignment_score_probabilities[ptr_p_idx as usize]
                + inner_sum * exp_minus_lambda;
            i += 1;
        }
        inner_sum *= exp_minus_lambda;

        while i <= high_alignment_score {
            ptr_p_idx += 1;
            inner_sum += alignment_score_probabilities[ptr_p_idx as usize];
            i += 1;
        }

        oldsum2 = oldsum;
        oldsum = inner_sum;

        iter_counter += 1;
        inner_sum /= iter_counter as f64;
        outer_sum += inner_sum;
    }

    let k = -(-2.0 * outer_sum).exp()
        / (first_term_closed_form * blast_expm1(-lambda));

    if k <= 0.0 {
        return Err("Computed K is non-positive".to_string());
    }

    Ok(k)
}

/// Compute Karlin-Altschul parameters from score frequency profile
/// Reference: NCBI Blast_KarlinBlkUngappedCalc (blast_stat.c:2699-2734)
///
/// # Arguments
/// * `sfp` - Score frequency profile
///
/// # Returns
/// Karlin parameters (Lambda, K, H) or error
pub fn compute_karlin_params_ungapped(
    sfp: &ScoreFreqProfile,
) -> Result<KarlinParams, String> {
    // Calculate Lambda
    let lambda = compute_lambda_nr(sfp, BLAST_KARLIN_LAMBDA0_DEFAULT)?;
    
    // Calculate H
    let h = compute_h_from_lambda(sfp, lambda)?;
    
    // Calculate K
    let k = compute_k_from_lambda_h(sfp, lambda, h)?;
    
    // Alpha and beta are not computed for ungapped (use defaults or lookup)
    // For ungapped, these are typically not used in E-value calculation
    Ok(KarlinParams {
        lambda,
        k,
        h,
        alpha: 0.7916, // BLOSUM62 ungapped default
        beta: -3.2,     // BLOSUM62 ungapped default
    })
}

/// Apply check_ideal logic: use ideal params if computed Lambda >= ideal Lambda
/// Reference: NCBI blast_stat.c:2796-2797
///
/// # Arguments
/// * `computed` - Computed Karlin parameters from sequence composition
/// * `ideal` - Ideal Karlin parameters (kbp_ideal)
///
/// # Returns
/// Final parameters (either computed or ideal)
pub fn apply_check_ideal(computed: KarlinParams, ideal: KarlinParams) -> KarlinParams {
    if computed.lambda >= ideal.lambda {
        // Use ideal (more conservative, smaller Lambda)
        ideal
    } else {
        // Use computed (smaller Lambda, less conservative)
        computed
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_compute_aa_composition() {
        // Test sequence: "ACDEFGHIKLMNPQRSTVWY" (20 standard AAs in NCBISTDAA encoding)
        // NCBISTDAA: A=1, C=3, D=4, E=5, F=6, G=7, H=8, I=9, K=10, L=11, M=12, N=13, P=14, Q=15, R=16, S=17, T=18, V=19, W=20, Y=22
        // Sequence with sentinels: [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 0]
        let seq = b"\x00\x01\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f\x10\x11\x12\x13\x14\x16\x00";
        let comp = compute_aa_composition(seq, 20);
        
        // Each AA should appear once, so frequency = 1/20 = 0.05
        let expected_freq = 1.0 / 20.0;
        // Check standard AAs (skip B=2, X=21, Z=23, U=24, *=25, O=26, J=27)
        let standard_aas = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22];
        for &i in &standard_aas {
            assert!((comp[i] - expected_freq).abs() < 0.01, "AA {} frequency: expected {}, got {}", i, expected_freq, comp[i]);
        }
        
        // Non-standard AAs should have frequency 0
        assert!(comp[2] < 0.01, "B should not appear");
        assert!(comp[21] < 0.01, "X should not appear");
        assert!(comp[23] < 0.01, "Z should not appear");
    }
    
    #[test]
    fn test_compute_std_aa_composition() {
        let comp = compute_std_aa_composition();
        let sum: f64 = comp.iter().sum();
        assert!((sum - 1.0).abs() < 0.01, "Standard composition should sum to 1.0");
    }
    
    #[test]
    fn test_compute_score_freq_profile() {
        let comp1 = compute_std_aa_composition();
        let comp2 = compute_std_aa_composition();
        let sfp = compute_score_freq_profile(&comp1, &comp2, -4, 11);
        
        assert!(sfp.obs_min() <= sfp.obs_max());
        assert!(sfp.score_avg() < 0.0, "Average score should be negative");
    }
    
    #[test]
    fn test_apply_check_ideal() {
        let ideal = KarlinParams {
            lambda: 0.3176,
            k: 0.134,
            h: 0.4012,
            alpha: 0.7916,
            beta: -3.2,
        };
        
        // Computed Lambda >= ideal Lambda -> use ideal
        let computed_high = KarlinParams {
            lambda: 0.35,
            k: 0.15,
            h: 0.45,
            alpha: 0.8,
            beta: -3.0,
        };
        let result = apply_check_ideal(computed_high, ideal);
        assert_eq!(result.lambda, ideal.lambda);
        
        // Computed Lambda < ideal Lambda -> use computed
        let computed_low = KarlinParams {
            lambda: 0.25,
            k: 0.10,
            h: 0.30,
            alpha: 0.7,
            beta: -2.5,
        };
        let result = apply_check_ideal(computed_low, ideal);
        assert_eq!(result.lambda, computed_low.lambda);
    }
}
