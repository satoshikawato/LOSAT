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
    let expected = (0.267 * 100.0 - 0.041_f64.glibc_ln()) / 2.0_f64.glibc_ln();
    assert!((bs - expected).abs() < 0.001);
}
