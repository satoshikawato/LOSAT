use super::*;

#[test]
fn test_compute_aa_composition() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2009-2012
    // Test sequence: "ACDEFGHIKLMNPQRSTVWY" (20 standard AAs in NCBISTDAA encoding)
    // NCBISTDAA: A=1, C=3, D=4, E=5, F=6, G=7, H=8, I=9, K=10, L=11, M=12, N=13, P=14, Q=15, R=16, S=17, T=18, V=19, W=20, Y=22
    // Sequence with sentinels: [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 0]
    let seq =
        b"\x00\x01\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f\x10\x11\x12\x13\x14\x16\x00";
    let comp = compute_aa_composition(seq, 20);

    // Each AA should appear once, so frequency = 1/20 = 0.05
    let expected_freq = 1.0 / 20.0;
    // Check standard AAs (skip B=2, X=21, Z=23, U=24, *=25, O=26, J=27)
    let standard_aas = [
        1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
    ];
    for &i in &standard_aas {
        assert!(
            (comp[i] - expected_freq).abs() < 0.01,
            "AA {} frequency: expected {}, got {}",
            i,
            expected_freq,
            comp[i]
        );
    }

    // Non-standard AAs should have frequency 0
    assert!(comp[2] < 0.01, "B should not appear");
    assert!(comp[21] < 0.01, "X should not appear");
    assert!(comp[23] < 0.01, "Z should not appear");
}

#[test]
fn test_compute_aa_composition_excludes_x_from_denominator() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:382-385
    let seq = b"\x00\x01\x15\x03\x00";
    let comp = compute_aa_composition(seq, 3);

    assert!((comp[ncbistdaa::A as usize] - 0.5).abs() < 1.0e-12);
    assert!((comp[ncbistdaa::C as usize] - 0.5).abs() < 1.0e-12);
    assert_eq!(comp[ncbistdaa::X as usize], 0.0);
}

#[test]
fn test_compute_std_aa_composition() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1794-1816
    let comp = compute_std_aa_composition();
    let sum: f64 = comp.iter().sum();
    assert!((sum - 1.0).abs() < 1.0e-12);
    let robinson_sum: f64 = STD_AA_FREQS.iter().sum();
    assert!((comp[ncbistdaa::A as usize] - 78.05 / robinson_sum).abs() < 1.0e-12);
    assert!((comp[ncbistdaa::C as usize] - 19.25 / robinson_sum).abs() < 1.0e-12);
    assert!((comp[ncbistdaa::R as usize] - 51.29 / robinson_sum).abs() < 1.0e-12);
    assert_eq!(comp[ncbistdaa::X as usize], 0.0);
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
