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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4157-4171
#[test]
fn test_evalue_from_raw_score_formula_and_guards() {
    let kbp = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };

    let e = evalue_from_raw_score(48, &kbp, 10_000.0);
    let expected = 10_000.0_f64 * (-0.267 * 48.0_f64 + 0.041_f64.glibc_ln()).glibc_exp();
    assert_eq!(e.to_bits(), expected.to_bits());
    assert_eq!(
        evalue_from_raw_score_with_log_k(48, &kbp, kbp.k.glibc_ln(), 10_000.0).to_bits(),
        expected.to_bits()
    );

    let e_max = evalue_from_raw_score(48, &kbp, i64::MAX as f64);
    let expected_max = (i64::MAX as f64) * (-0.267 * 48.0_f64 + 0.041_f64.glibc_ln()).glibc_exp();
    assert_eq!(e_max.to_bits(), expected_max.to_bits());

    assert_eq!(
        evalue_from_raw_score(48, &kbp, 0.0).to_bits(),
        0.0_f64.to_bits()
    );

    // Negative Lambda, K, or H each return -1.0 even when search space is 0.
    for invalid in [
        KarlinParams {
            lambda: -0.267,
            ..kbp
        },
        KarlinParams { k: -0.041, ..kbp },
        KarlinParams { h: -0.14, ..kbp },
    ] {
        assert_eq!(
            evalue_from_raw_score(48, &invalid, 10_000.0).to_bits(),
            (-1.0_f64).to_bits()
        );
        assert_eq!(
            evalue_from_raw_score(48, &invalid, 0.0).to_bits(),
            (-1.0_f64).to_bits()
        );
    }

    // Zero H or zero Lambda is not negative (`< 0.0`), matching C's strict `< 0.` check.
    let zero_h = KarlinParams { h: 0.0, ..kbp };
    assert_eq!(
        evalue_from_raw_score(48, &zero_h, 10_000.0).to_bits(),
        expected.to_bits()
    );
    let zero_lambda = KarlinParams { lambda: 0.0, ..kbp };
    assert_eq!(
        evalue_from_raw_score(48, &zero_lambda, 10_000.0).to_bits(),
        (10_000.0_f64 * 0.041_f64.glibc_ln().glibc_exp()).to_bits()
    );
}
