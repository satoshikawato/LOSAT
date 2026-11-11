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

/// Calculate E-value from raw score and search space (NCBI BLAST_KarlinStoE_simple compatible)
///
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4157-4171
///
/// Formula: E = searchsp * exp(-Lambda * S + logK)
/// where S is the raw score, searchsp is the effective search space
pub fn evalue_from_raw_score(raw_score: i32, params: &KarlinParams, search_space: f64) -> f64 {
    if params.lambda < 0.0 || params.k < 0.0 {
        return -1.0;
    }
    search_space * ((-params.lambda * (raw_score as f64)) + params.k.ln()).exp()
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
}
