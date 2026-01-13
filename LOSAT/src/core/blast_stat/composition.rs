//! Karlin parameter calculation from composition.
//!
//! Reference: blast_stat.c Blast_KarlinBlkUngappedCalc

use super::karlin_params::KarlinParams;
use crate::utils::matrix::{blosum62_score, BLASTAA_SIZE};

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
