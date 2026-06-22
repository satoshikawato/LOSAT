//! Karlin parameter calculation from composition.
//!
//! Reference: blast_stat.c Blast_KarlinBlkUngappedCalc

use super::karlin_params::KarlinParams;
use crate::utils::matrix::{blosum62_score, ncbistdaa, BLASTAA_SIZE};

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_stat.h:121-122
// ```c
// #define BLAST_SCORE_MIN INT2_MIN
// #define BLAST_SCORE_MAX INT2_MAX
// ```
const BLAST_SCORE_MIN: i32 = i16::MIN as i32;
const BLAST_SCORE_MAX: i32 = i16::MAX as i32;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:56-72
// ```c
// #define BLAST_SCORE_RANGE_MAX (BLAST_SCORE_MAX - BLAST_SCORE_MIN)
// #define BLAST_KARLIN_K_SUMLIMIT_DEFAULT 0.0001
// #define BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT 1.0e-5
// #define BLAST_KARLIN_LAMBDA_ITER_DEFAULT 17
// #define BLAST_KARLIN_LAMBDA0_DEFAULT 0.5
// #define BLAST_KARLIN_K_ITER_MAX 100
// ```
const BLAST_SCORE_RANGE_MAX: i32 = BLAST_SCORE_MAX - BLAST_SCORE_MIN;
const BLAST_KARLIN_K_SUMLIMIT_DEFAULT: f64 = 0.0001;
const BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT: f64 = 1.0e-5;
const BLAST_KARLIN_LAMBDA_ITER_DEFAULT: i32 = 17;
const BLAST_KARLIN_LAMBDA0_DEFAULT: f64 = 0.5;
const BLAST_KARLIN_K_ITER_MAX: i32 = 100;

/// Standard amino acid frequencies (Robinson probabilities).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1794-1816
/// ```c
/// static BLAST_LetterProb Robinson_prob[] = {
///       { 'A', 78.05 },
///       { 'C', 19.25 },
///       { 'D', 53.64 },
///       ...
///       { 'Y', 32.16 }
///    };
/// ```
const STD_AA_FREQS: [f64; 20] = [
    78.05, 19.25, 53.64, 62.95, 38.56, 73.77, 21.99, 51.42, 57.44, 90.19, 22.43, 44.87, 52.03,
    42.64, 51.29, 71.20, 58.41, 64.41, 13.30, 32.16,
];

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1863-1879
// ```c
// for (index=0; index<(int)DIM(STD_AMINO_ACID_FREQS); index++)
// {
//    if (alphabet_code == BLASTAA_SEQ_CODE)
//    {
//       residues[index] =
//          AMINOACID_TO_NCBISTDAA[toupper((unsigned char) STD_AMINO_ACID_FREQS[index].ch)];
//    }
// }
// ```
const STD_AA_NCBISTDAA: [u8; 20] = [
    ncbistdaa::A,
    ncbistdaa::C,
    ncbistdaa::D,
    ncbistdaa::E,
    ncbistdaa::F,
    ncbistdaa::G,
    ncbistdaa::H,
    ncbistdaa::I,
    ncbistdaa::K,
    ncbistdaa::L,
    ncbistdaa::M,
    ncbistdaa::N,
    ncbistdaa::P,
    ncbistdaa::Q,
    ncbistdaa::R,
    ncbistdaa::S,
    ncbistdaa::T,
    ncbistdaa::V,
    ncbistdaa::W,
    ncbistdaa::Y,
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1984-2014
    // ```c
    // for (lp = str, lpmax = lp+length; lp < lpmax; lp++)
    // {
    //    ++rcp->comp[(int)(*lp & mask)];
    // }
    // ...
    // for (index=0; index<sbp->ambig_occupy; index++)
    // {
    //    rcp->comp[sbp->ambiguous_res[index]] = 0;
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:382-385
    // ```c
    // sbp->read_in_matrix = TRUE;
    // BLAST_ScoreSetAmbigRes(sbp, 'X');
    // sbp->name = BLAST_StrToUpper(scoring_options->matrix);
    // ```
    comp[ncbistdaa::X as usize] = 0;

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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1887-1901
    // ```c
    // retval = Blast_GetStdAlphabet(sbp->alphabet_code, residues, DIM(STD_AMINO_ACID_FREQS));
    // for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++)
    // {
    //    rfp->prob[residues[index]] = STD_AMINO_ACID_FREQS[index].p;
    // }
    // ```
    for (index, &residue) in STD_AA_NCBISTDAA.iter().enumerate() {
        freq[residue as usize] = STD_AA_FREQS[index];
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1909-1910
    // ```c
    // Blast_ResFreqNormalize(sbp, rfp, 1.0);
    // ```
    let sum: f64 = freq.iter().sum();
    if sum > 0.0 {
        for f in freq.iter_mut() {
            *f /= sum;
        }
    }

    freq
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2833-2852
// ```c
// Blast_ScoreBlkKbpIdealCalc(BlastScoreBlk* sbp)
// {
//    stdrfp = Blast_ResFreqNew(sbp);
//    Blast_ResFreqStdComp(sbp, stdrfp);
//    sfp = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
//    BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
//    sbp->kbp_ideal = Blast_KarlinBlkNew();
//    Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);
// }
// ```
pub fn compute_blosum62_ideal_karlin_params() -> Result<KarlinParams, String> {
    let std_comp = compute_std_aa_composition();
    let sfp = compute_score_freq_profile(&std_comp, &std_comp, -4, 11);
    compute_karlin_params_ungapped(&sfp)
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2183-2192
    // ```c
    // obs_min = BLAST_SCORE_MIN;
    // obs_max = BLAST_SCORE_MIN;
    // for (score = score_min; score <= score_max; score++) {
    //     if (sprob[score] != 0.0) {
    //         score_sum += sprob[score];
    //         obs_max = score;
    //         if (obs_min == BLAST_SCORE_MIN)
    //             obs_min = score;
    //     }
    // }
    // ```
    let mut obs_min = BLAST_SCORE_MIN;
    let mut obs_max = BLAST_SCORE_MIN;

    for score in score_min..=score_max {
        let prob = sfp.get_prob(score);
        if prob > 0.0 {
            score_sum += prob;
            obs_max = obs_max.max(score);
            if obs_min == BLAST_SCORE_MIN {
                obs_min = score;
            }
        }
    }

    sfp.obs_min = obs_min;
    sfp.obs_max = obs_max;

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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:405-418
// ```c
// Int4 BLAST_Gcd(Int4 a, Int4 b)
// {
//     ...
// }
// ```
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2099-2109
// ```c
// static Int2 BlastScoreChk(Int4 lo, Int4 hi)
// {
//     if (lo >= 0 || hi <= 0 || lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX)
//         return 1;
//     if (hi - lo > BLAST_SCORE_RANGE_MAX)
//         return 1;
//     return 0;
// }
// ```
fn blast_score_chk(lo: i32, hi: i32) -> Result<(), String> {
    if lo >= 0 || hi <= 0 || lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX {
        return Err("Invalid score range".to_string());
    }
    if hi - lo > BLAST_SCORE_RANGE_MAX {
        return Err("Score range exceeds BLAST_SCORE_RANGE_MAX".to_string());
    }
    Ok(())
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:33-55
// ```c
// double BLAST_Expm1(double x)
// {
//     ...
// }
// ```
#[rustfmt::skip]
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:444-470
// ```c
// double BLAST_Powi(double x, Int4 n)
// {
//     ...
// }
// ```
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2491-2563
// ```c
// static double
// NlmKarlinLambdaNR(double* probs, Int4 d, Int4 low, Int4 high, double lambda0,
//                   double tolx, Int4 itmax, Int4 maxNewton, Int4 * itn)
// {
//     ...
// }
// ```
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

        if k >= max_newton || (was_newton && f.abs() > 0.9 * fold.abs()) || g >= 0.0 {
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2567-2598
// ```c
// double Blast_KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess)
// {
//     ...
// }
// ```
fn compute_lambda_nr(sfp: &ScoreFreqProfile, initial_guess: f64) -> Result<f64, String> {
    let low = sfp.obs_min();
    let high = sfp.obs_max();

    if sfp.score_avg() >= 0.0 {
        return Err("Expected score must be negative".to_string());
    }

    blast_score_chk(low, high)?;

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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:551-570
// ```c
// static double
// s_CalcLambda(double probs[], int min_score, int max_score, double lambda0)
// {
//     ...
//     return Blast_KarlinLambdaNR(&freq, lambda0);
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2567-2598
// ```c
// double Blast_KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess)
// {
//     ...
// }
// ```
pub fn compute_lambda_from_score_probs(
    score_probs: &[f64],
    obs_min: i32,
    obs_max: i32,
    initial_lambda: f64,
) -> Result<f64, String> {
    if obs_min > obs_max {
        return Err("Invalid observed score range".to_string());
    }

    let mut sfp = ScoreFreqProfile::new(obs_min, obs_max);
    let mut score_sum = 0.0;
    let mut score_avg = 0.0;
    let mut first_obs = None;
    let mut last_obs = None;

    for (index, &prob) in score_probs.iter().enumerate() {
        let score = obs_min + index as i32;
        sfp.set_prob(score, prob);
        if prob > 0.0 {
            score_sum += prob;
            score_avg += prob * score as f64;
            first_obs.get_or_insert(score);
            last_obs = Some(score);
        }
    }

    if score_sum <= 0.0 {
        return Err("Score probabilities sum to zero".to_string());
    }

    for (index, &prob) in score_probs.iter().enumerate() {
        let score = obs_min + index as i32;
        sfp.set_prob(score, prob / score_sum);
    }

    sfp.obs_min = first_obs.unwrap_or(obs_min);
    sfp.obs_max = last_obs.unwrap_or(obs_max);
    sfp.score_avg = score_avg / score_sum;

    compute_lambda_nr(&sfp, initial_lambda)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2607-2633
// ```c
// static double BlastKarlinLtoH(Blast_ScoreFreq* sfp, double lambda)
// {
//     ...
// }
// ```
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
        lambda * (lambda * high as f64 + sum.ln()).exp()
    };

    Ok(h)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2247-2418
// ```c
// static double BlastKarlinLHtoK(Blast_ScoreFreq* sfp, double lambda, double H)
// {
//     ...
// }
// ```
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

    let array_len = (BLAST_KARLIN_K_ITER_MAX as usize) * (range as usize) + 1;
    let mut alignment_score_probabilities = vec![0.0; array_len];
    let mut outer_sum = 0.0;
    let mut low_alignment_score = 0;
    let mut high_alignment_score = 0;
    let mut inner_sum = 1.0;
    alignment_score_probabilities[0] = 1.0;

    let mut iter_counter = 0;
    while iter_counter < BLAST_KARLIN_K_ITER_MAX && inner_sum > BLAST_KARLIN_K_SUMLIMIT_DEFAULT {
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
            inner_sum =
                alignment_score_probabilities[ptr_p_idx as usize] + inner_sum * exp_minus_lambda;
            i += 1;
        }
        inner_sum *= exp_minus_lambda;

        while i <= high_alignment_score {
            ptr_p_idx += 1;
            inner_sum += alignment_score_probabilities[ptr_p_idx as usize];
            i += 1;
        }

        iter_counter += 1;
        inner_sum /= iter_counter as f64;
        outer_sum += inner_sum;
    }

    let k = -(-2.0 * outer_sum).exp() / (first_term_closed_form * blast_expm1(-lambda));
    if k <= 0.0 {
        return Err("Computed K is non-positive".to_string());
    }
    Ok(k)
}

/// Compute Karlin-Altschul parameters from score frequency profile
pub fn compute_karlin_params_ungapped(sfp: &ScoreFreqProfile) -> Result<KarlinParams, String> {
    let lambda = compute_lambda_nr(sfp, BLAST_KARLIN_LAMBDA0_DEFAULT)?;
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1105-1144
// ```c
// Blast_CompositionBasedStats(..., double (*calc_lambda)(double*,int,int,double), ...)
// {
//     ...
//     correctUngappedLambda =
//         calc_lambda(scoreArray, obs_min, obs_max, ss->ungappedLambda);
//     ...
// }
// ```
pub fn compute_karlin_params_from_score_probs(
    score_probs: &[f64],
    obs_min: i32,
    obs_max: i32,
    initial_lambda: f64,
) -> Result<KarlinParams, String> {
    if obs_min > obs_max {
        return Err("Invalid observed score range".to_string());
    }

    let mut sfp = ScoreFreqProfile::new(obs_min, obs_max);
    let mut score_sum = 0.0;
    let mut score_avg = 0.0;
    let mut first_obs = None;
    let mut last_obs = None;

    for (index, &prob) in score_probs.iter().enumerate() {
        let score = obs_min + index as i32;
        sfp.set_prob(score, prob);
        if prob > 0.0 {
            score_sum += prob;
            score_avg += prob * score as f64;
            first_obs.get_or_insert(score);
            last_obs = Some(score);
        }
    }

    if score_sum <= 0.0 {
        return Err("Score probabilities sum to zero".to_string());
    }

    for (index, &prob) in score_probs.iter().enumerate() {
        let score = obs_min + index as i32;
        sfp.set_prob(score, prob / score_sum);
    }

    sfp.obs_min = first_obs.unwrap_or(obs_min);
    sfp.obs_max = last_obs.unwrap_or(obs_max);
    sfp.score_avg = score_avg / score_sum;

    let lambda = compute_lambda_nr(&sfp, initial_lambda)?;
    let h = compute_h_from_lambda(&sfp, lambda)?;
    let k = compute_k_from_lambda_h(&sfp, lambda, h)?;

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
    if computed.lambda >= ideal.lambda {
        ideal
    } else {
        computed
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_aa_composition_excludes_x_from_denominator() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:382-385
        // ```c
        // sbp->read_in_matrix = TRUE;
        // BLAST_ScoreSetAmbigRes(sbp, 'X');
        // ```
        let seq = b"\x00\x01\x15\x03\x00";
        let comp = compute_aa_composition(seq, 3);

        assert!((comp[ncbistdaa::A as usize] - 0.5).abs() < 1.0e-12);
        assert!((comp[ncbistdaa::C as usize] - 0.5).abs() < 1.0e-12);
        assert_eq!(comp[ncbistdaa::X as usize], 0.0);
    }

    #[test]
    fn test_compute_std_aa_composition_matches_robinson_order() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1794-1816
        // ```c
        // static BLAST_LetterProb Robinson_prob[] = {
        //       { 'A', 78.05 },
        //       { 'C', 19.25 },
        //       ...
        //       { 'R', 51.29 },
        //       ...
        // };
        // ```
        let comp = compute_std_aa_composition();
        let robinson_sum: f64 = STD_AA_FREQS.iter().sum();

        assert!((comp.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
        assert!((comp[ncbistdaa::A as usize] - 78.05 / robinson_sum).abs() < 1.0e-12);
        assert!((comp[ncbistdaa::C as usize] - 19.25 / robinson_sum).abs() < 1.0e-12);
        assert!((comp[ncbistdaa::R as usize] - 51.29 / robinson_sum).abs() < 1.0e-12);
        assert_eq!(comp[ncbistdaa::X as usize], 0.0);
    }
}
