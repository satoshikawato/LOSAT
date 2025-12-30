use super::search_space::SearchSpace;
use super::tables::KarlinParams;

/// Calculate bit score from raw score using Karlin-Altschul statistics
///
/// Formula: S' = (lambda * S - ln(K)) / ln(2)
/// where S is the raw score, S' is the bit score
pub fn bit_score(raw_score: i32, params: &KarlinParams) -> f64 {
    let ln2 = 2.0_f64.ln();
    (params.lambda * (raw_score as f64) - params.k.ln()) / ln2
}

/// Calculate E-value from bit score and search space
///
/// Formula: E = m * n * 2^(-S')
/// where m*n is the effective search space, S' is the bit score
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
pub fn raw_score_from_evalue(
    e_value: f64,
    params: &KarlinParams,
    search_space: &SearchSpace,
) -> i32 {
    if e_value <= 0.0 {
        return i32::MAX;
    }
    let score = (params.k.ln() + search_space.effective_space.ln() - e_value.ln()) / params.lambda;
    score.ceil() as i32
}

/// Calculate raw score from E-value with gap decay adjustment (NCBI BLAST_Cutoffs compatible)
///
/// When dodecay is true, the E-value is adjusted by gap_decay_divisor to get a more
/// lenient cutoff score suitable for sum statistics linking.
///
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4112-4118
///
/// Formula (with dodecay):
///   adjusted_e = e_value * gap_decay_divisor(gap_decay_rate, 1)
///   S = (ln(K) + ln(m*n) - ln(adjusted_e)) / lambda
pub fn raw_score_from_evalue_with_decay(
    e_value: f64,
    params: &KarlinParams,
    search_space: &SearchSpace,
    dodecay: bool,
    gap_decay_rate: f64,
) -> i32 {
    use super::sum_statistics::gap_decay_divisor;
    
    if e_value <= 0.0 {
        return i32::MAX;
    }
    
    let adjusted_e = if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
        // Apply gap decay divisor to get a more lenient cutoff
        // NCBI: e *= BLAST_GapDecayDivisor(gap_decay_rate, 1)
        e_value * gap_decay_divisor(gap_decay_rate, 1)
    } else {
        e_value
    };
    
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
///
/// Formula: E = searchsp * exp(-Lambda * S + logK)
/// where S is the raw score, searchsp is the effective search space
pub fn evalue_from_raw_score(
    raw_score: i32,
    params: &KarlinParams,
    search_space: f64,
) -> f64 {
    if params.lambda < 0.0 || params.k < 0.0 {
        return -1.0;
    }
    search_space * ((-params.lambda * (raw_score as f64)) + params.k.ln()).exp()
}

/// Simple E-value calculation without length adjustment (BLEMIR default)
///
/// This is the original BLEMIR calculation for backward compatibility
pub fn simple_evalue(
    raw_score: i32,
    q_len: usize,
    db_len: usize,
    params: &KarlinParams,
) -> (f64, f64) {
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

        // Test with a typical protein alignment score
        let score = 100;
        let bs = bit_score(score, &params);

        // Bit score should be positive for positive raw scores
        assert!(bs > 0.0);

        // Verify the formula: S' = (lambda * S - ln(K)) / ln(2)
        let expected = (0.267 * 100.0 - 0.041_f64.ln()) / 2.0_f64.ln();
        assert!((bs - expected).abs() < 0.001);
    }

    #[test]
    fn test_evalue() {
        let _params = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };

        let search_space = SearchSpace {
            effective_query_len: 100.0,
            effective_db_len: 1000.0,
            effective_space: 100000.0,
            length_adjustment: 0,
        };

        let bs = 50.0; // bit score
        let ev = evalue(bs, &search_space);

        // E-value should be very small for high bit scores
        assert!(ev < 1.0);

        // Verify the formula: E = m * n * 2^(-S')
        let expected = 100000.0 * 2.0_f64.powf(-50.0);
        assert!((ev - expected).abs() < 1e-10);
    }

    #[test]
    fn test_round_trip() {
        let params = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };

        let original_score = 100;
        let bs = bit_score(original_score, &params);
        let recovered_score = raw_score_from_bit_score(bs, &params);

        // Should recover the original score (within rounding)
        assert!((original_score - recovered_score).abs() <= 1);
    }
}
