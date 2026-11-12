use super::tables::KarlinParams;

/// Calculate bit score from raw score using Karlin-Altschul statistics
///
/// Formula: S' = (lambda * S - ln(K)) / ln(2)
/// where S is the raw score, S' is the bit score
pub fn bit_score(raw_score: i32, params: &KarlinParams) -> f64 {
    let ln2 = 2.0_f64.ln();
    (params.lambda * (raw_score as f64) - params.k.ln()) / ln2
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
}
