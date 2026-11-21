use super::*;

#[test]
fn test_gap_decay_divisor() {
    // For decay_rate = 0.5, nsegs = 1: (1-0.5) * 0.5^0 = 0.5
    assert!((gap_decay_divisor(0.5, 1) - 0.5).abs() < 1e-10);

    // For decay_rate = 0.5, nsegs = 2: (1-0.5) * 0.5^1 = 0.25
    assert!((gap_decay_divisor(0.5, 2) - 0.25).abs() < 1e-10);

    // For decay_rate = 0.1, nsegs = 1: (1-0.1) * 0.1^0 = 0.9
    assert!((gap_decay_divisor(0.1, 1) - 0.9).abs() < 1e-10);
}

#[test]
fn test_ln_factorial() {
    // ln(1!) = 0
    assert!((blast_ln_factorial(1.0) - 0.0).abs() < 1e-10);

    // ln(2!) = ln(2) ≈ 0.693
    assert!((blast_ln_factorial(2.0) - 2.0_f64.glibc_ln()).abs() < 1e-10);

    // ln(5!) = ln(120) ≈ 4.787
    assert!((blast_ln_factorial(5.0) - 120.0_f64.glibc_ln()).abs() < 1e-10);
}

#[test]
fn test_p_to_e_and_e_to_p() {
    // Round trip test
    let p = 0.01;
    let e = p_to_e(p);
    let p_back = e_to_p(e);
    assert!((p - p_back).abs() < 1e-10);

    // Edge cases
    assert_eq!(e_to_p(0.0), 0.0);
    assert!((e_to_p(f64::INFINITY) - 1.0).abs() < 1e-10);
}

#[test]
fn test_small_gap_sum_e_single_hsp() {
    // For a single HSP, sum_e = searchsp * exp(-xsum)
    let xsum = 10.0;
    let searchsp = 1000000_i64;
    let result = small_gap_sum_e(50, 1, xsum, 100, 1000, searchsp, 1.0);

    let expected = (searchsp as f64) * (-xsum).glibc_exp();
    assert!((result - expected).abs() / expected < 0.01);
}

#[test]
fn test_normalize_score() {
    let lambda: f64 = 0.267;
    let k: f64 = 0.041;
    let log_k = k.glibc_ln();
    let score = 100;

    let xsum = normalize_score(score, lambda, log_k);
    let expected = lambda * 100.0 - log_k;
    assert!((xsum - expected).abs() < 1e-10);
}

#[test]
fn sum_p_calc_reused_gives_the_computed_value() {
    // Repeats after one and after two other calls, a key that left the two kept, and a
    // value of `s` that differs from a kept one in its last bit.
    let s: f64 = -366.25;
    let args = [
        (157, s),
        (157, s),
        (12, -20.0),
        (157, s),
        (30, -50.5),
        (12, -20.0),
        (157, s),
        (157, f64::from_bits(s.to_bits() + 1)),
        (157, s),
    ];
    for (r, s) in args {
        assert_eq!(
            blast_sum_p_calc_reused(r, s).to_bits(),
            blast_sum_p_calc(r, s).to_bits(),
            "r {r} s {s}"
        );
    }
}
