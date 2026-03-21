//! Karlin parameter calculation from composition.
//!
//! Reference: blast_stat.c Blast_KarlinBlkUngappedCalc

use super::karlin_params::KarlinParams;
use crate::utils::matrix::{blosum62_score, ncbistdaa, BLASTAA_SIZE};

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
    let mut obs_min = score_max;
    let mut obs_max = score_min;

    for score in score_min..=score_max {
        let prob = sfp.get_prob(score);
        if prob > 0.0 {
            score_sum += prob;
            obs_max = obs_max.max(score);
            if obs_min == score_max {
                obs_min = score;
            }
        }
    }

    sfp.obs_min = if obs_min <= obs_max { obs_min } else { 0 };
    sfp.obs_max = if obs_min <= obs_max { obs_max } else { 0 };

    let mut score_avg = 0.0;
    if score_sum > 0.0001 || score_sum < -0.0001 {
        for score in score_min..=score_max {
            let prob = sfp.get_prob(score);
            if prob > 0.0 {
                score_avg += (score as f64) * prob;
            }
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
            if d == 1 {
                break;
            }
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

        if sum_deriv.abs() < 1e-10 {
            break;
        }

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
            if prob > 0.0 {
                sum += prob * (lambda * score as f64).exp();
            }
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

    let lambda = compute_lambda_nr(&sfp, initial_lambda.max(1.0e-6))?;
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
