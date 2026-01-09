//! Bit score and E-value calculations.
//!
//! Reference: blast_stat.c BLAST_KarlinStoE_simple, BlastKarlinEtoS_simple

use super::karlin_params::KarlinParams;
use super::search_space::SearchSpace;
use super::sum_statistics::gap_decay_divisor;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_score() {
        let params = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };
        let score = 100;
        let bs = bit_score(score, &params);
        assert!(bs > 0.0);
        let expected = (0.267 * 100.0 - 0.041_f64.ln()) / 2.0_f64.ln();
        assert!((bs - expected).abs() < 0.001);
    }
}
