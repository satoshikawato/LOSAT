// EXPERIMENT: Dekker's product error equals the fused one inside the range
// in which the Wasm build uses it.
// Boost (NCBI 598d8ae6 copy): c++/src/algo/blast/core/boost_erf.c:143-143
// The test checks, for 2 million random pairs, that Dekker's error term has the bits of
// the fused multiply-add error term.
#[test]
fn x_dekker_product_error_matches_fused_multiply_add() {
    let mut state = 0x1234_5678_9ABC_DEF1u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut checked = 0u32;
    for _ in 0..2_000_000 {
        let mantissa_a = (next() >> 12) | (1u64 << 52);
        let mantissa_b = (next() >> 12) | (1u64 << 52);
        let exp_a = (next() % 1700) as i32 - 850 - 52;
        let exp_b = (next() % 1700) as i32 - 850 - 52;
        let mut a = (mantissa_a as f64) * 2f64.powi(exp_a.clamp(-1000, 900));
        let mut b = (mantissa_b as f64) * 2f64.powi(exp_b.clamp(-1000, 900));
        if next() & 1 == 1 {
            a = -a;
        }
        if next() & 2 == 2 {
            b = -b;
        }
        // Also values with few significant bits, as the series divisors have.
        if next() % 7 == 0 {
            a = (next() % 41) as f64;
        }
        let product = a * b;
        const TINY: f64 = super::X_TWO_POW_M900;
        const HUGE: f64 = super::X_TWO_POW_900;
        assert_eq!(TINY, 2f64.powi(-900));
        assert_eq!(HUGE, 2f64.powi(900));
        if product.abs() > TINY && product.abs() < HUGE && a.abs() < HUGE && b.abs() < HUGE {
            let fused = a.mul_add(b, -product);
            let dekker = super::x_dekker_product_error(a, b, product);
            assert_eq!(fused.to_bits(), dekker.to_bits(), "a={a:e} b={b:e}");
            checked += 1;
        }
    }
    assert!(checked > 500_000, "only {checked} products in range");
}

use super::*;
use crate::config::ScoringMatrix;

#[test]
fn x_erfc_memo_keeps_the_reference_bits_as_it_grows() {
    // No NCBI counterpart (reuse only): `X_ERFC_MEMO` grows from 64 to 4,096 slots as this
    // thread fills it; every value, computed or read back after the growth, is ErfC's.
    let zs: Vec<f64> = (0..6_000).map(|i| -5.9 + f64::from(i) * 1.9e-3).collect();
    for _ in 0..2 {
        for &z in &zs {
            assert_eq!(
                erfc_ncbi(z).to_bits(),
                erf_impl_ncbi(z, true).to_bits(),
                "z={z:e}"
            );
        }
    }
}

// LOSAT_X_ERFABSORB (graduated 2026-10-26): the reference ErfC is 2.0 for every z < -6, and
// `erfc_ncbi`, which does not call it there, has the reference's bits everywhere.
// Boost (NCBI 598d8ae6 copy): c++/src/algo/blast/core/boost_erf.c:58-59
// The argument of `erfc_ncbi` bounds `ErfImpl` at -z with `invert` set below 2.2e-17 for z <= -6; this
// checks the reference on a grid of [-28.5, -6], on random arguments below -6 (dense next to
// -6 and 28) and on the doubles just below -6, and that `erfc_ncbi`
// equals `erf_impl_ncbi` on random arguments in [-40, 40]. A reference call costs about
// 5,000 cycles, so the sets are kept to about 100,000 calls.
#[test]
fn x_erfc_below_minus_six_is_two_and_the_switched_function_has_the_reference_bits() {
    let two = 2.0f64.to_bits();
    let check = |z: f64| {
        assert_eq!(erf_impl_ncbi(z, true).to_bits(), two, "ErfC({z:e})");
    };
    for i in 0..=22_500u32 {
        check(-28.5 + i as f64 * 1e-3);
    }
    let mut below_six = -6.0f64;
    for _ in 0..20_000 {
        below_six = f64::from_bits(below_six.to_bits() + 1);
        check(below_six);
    }
    let mut state = 0x0DDB_A11C_AFE5_EED5u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..15_000 {
        let u = next();
        check(-6.0 - u * u * 24.0);
        check(-28.0 + (u - 0.5) * 1e-3);
        check(-6.0 - u * 1e300);
    }
    check(f64::NEG_INFINITY);
    for _ in 0..30_000 {
        let z = (next() - 0.5) * 80.0;
        assert_eq!(
            erfc_ncbi(z).to_bits(),
            erf_impl_ncbi(z, true).to_bits(),
            "z={z:e}"
        );
    }
    assert_eq!(
        erfc_ncbi(f64::NAN).to_bits(),
        erf_impl_ncbi(f64::NAN, true).to_bits()
    );
}

// Boost (NCBI 598d8ae6 copy): c++/src/algo/blast/core/boost_erf.c:252-255
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:5180-5231:
// Compare every captured preliminary and scaled post-redo StoE call.
#[test]
fn pinned_tblastn_spouge_calls_match_ncbi_double_bits() {
    let cases = [
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/multi_query_20260924_default.spouge.tsv")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/multi_query_20260924_control.spouge.tsv")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/seg_hard_query_20260924_default.spouge.tsv")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/seg_hard_query_20260924_control.spouge.tsv")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/run_20260923_default.spouge.tsv")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/run_20260923_control.spouge.tsv")),
    ];
    for case in cases {
        for row in case.lines() {
            let f: Vec<_> = row.split('\t').collect();
            let mut karlin =
                crate::stats::tables::lookup_protein_params_gapped(ScoringMatrix::Blosum62);
            karlin.lambda = f[4].parse().unwrap();
            karlin.k = f[5].parse().unwrap();
            let gumbel = BlastGumbelBlk {
                lambda: f[7].parse().unwrap(),
                c: f[8].parse().unwrap(),
                g: f[9].parse().unwrap(),
                a: f[10].parse().unwrap(),
                alpha: f[11].parse().unwrap(),
                sigma: f[12].parse().unwrap(),
                a_un: f[13].parse().unwrap(),
                alpha_un: f[14].parse().unwrap(),
                b: f[15].parse().unwrap(),
                beta: f[16].parse().unwrap(),
                tau: f[17].parse().unwrap(),
                db_length: f[18].parse().unwrap(),
            };
            let actual = blast_spouge_stoe(
                f[1].parse().unwrap(),
                &karlin,
                &gumbel,
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
            );
            let expected: f64 = f[20].parse().unwrap();
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "{row} actual={actual:.17e}"
            );
        }
    }
}

#[test]
fn pinned_tblastn_erfc_calls_match_ncbi_double_bits() {
    let cases = [
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/erfc_20260924/multi_query_20260924_default.tsv"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/erfc_20260924/multi_query_20260924_control.tsv"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/erfc_20260924/seg_hard_query_20260924_default.tsv"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/erfc_20260924/seg_hard_query_20260924_control.tsv"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/erfc_20260924/run_20260923_default.tsv"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/erfc_20260924/run_20260923_control.tsv"
        )),
    ];
    for case in cases {
        for row in case.lines() {
            let fields: Vec<_> = row.split('\t').collect();
            let z: f64 = fields[1].parse().unwrap();
            let expected: f64 = fields[2].parse().unwrap();
            let actual = erfc_ncbi(z);
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "ErfC({z:.17e}) actual {actual:.17e} expected {expected:.17e}"
            );
        }
    }
}

// Boost (NCBI 598d8ae6 copy): c++/src/algo/blast/core/boost_erf.c:143,171,199,229
// x86-64 NCBI rounds the long-double product to 64 bits and then stores it to the double:
// two roundings. For this argument (P_n_F of swissprot blastp Query_17 x P01717.2) the
// 64-bit product 0x9.6644bd26f0b6cp-11 lies halfway between two doubles and the store rounds
// to even, up, while the exact product lies just below halfway (one rounding gives ...16d).
#[test]
fn erfc_rounds_the_long_double_product_twice_as_ncbi() {
    let z = f64::from_bits(0x4000_08CF_53B9_A16E);
    assert_eq!(erfc_ncbi(z).to_bits(), 0x3F72_CC89_7A4D_E16E);
}

// Boost (NCBI 598d8ae6 copy): c++/src/algo/blast/core/boost_erf.c:207-229
// From z near 25.9 the low part of the double-double exponential is subnormal, from 26.6 its
// high part, and from about 27.3 it is zero (NCBI's 80-bit expl is still normal there, and the
// stored product is zero all the same). The x87 model gets subnormal and zero operands there:
// no panic under the test profile's overflow checks, results at most 1e-290 that do not
// increase with z, and zero from 27.3 on. NCBI's bits there are not reproduced (the e-value is
// 0.0 either way, docs/evidence/spouge_x87_20261027/manifest.txt).
#[test]
fn erfc_in_the_underflow_range_of_the_double_double_exponential() {
    let mut previous = f64::INFINITY;
    for z in [26.0, 26.5, 26.62, 27.0, 27.2257, 27.3, 27.99] {
        let r = erfc_ncbi(z);
        assert!(
            (0.0..=1e-290).contains(&r) && r <= previous,
            "ErfC({z}) = {r:e}"
        );
        if z >= 27.3 {
            assert_eq!(r.to_bits(), 0, "ErfC({z})");
        }
        previous = r;
    }
    for i in 0..20_000 {
        let z = 25.9 + i as f64 * (2.1 / 20_000.0);
        let r = erf_impl_ncbi(z, true);
        assert!((0.0..=1e-290).contains(&r), "ErfC({z}) = {r:e}");
    }
}

fn scaled_blosum62_11_1() -> (ProteinScoringSpec, KarlinParams) {
    let spec = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: 11,
        gap_extend: 1,
    };
    let kbp = KarlinParams {
        lambda: 0.267 / 32.0,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    (spec, kbp)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:5223-5229
// swissprot blastp Query_17 (120 letters) x P01717.2 (112 letters), score 6558 of the
// matrix scaled by 32: NCBI's e-value 7.11685934590109512e-21 (outfmt 8
// `{ 71168593459011, 10, -34 }`) follows from the doubly rounded ErfC above.
#[test]
fn spouge_stoe_of_swissprot_query_17_x_p01717_has_the_ncbi_bits() {
    let (spec, kbp) = scaled_blosum62_11_1();
    let gbp = lookup_protein_gumbel_params(&spec, 186_063_409).unwrap();
    assert_eq!(
        blast_spouge_stoe(6558, &kbp, &gbp, 120, 112).to_bits(),
        0x3BC0_CDE0_A46C_923E
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:5216-5229
// The rows of docs/evidence/spouge_x87_20261027 (its manifest.txt): ErfC arguments and
// Spouge calls of a random sample where one rounding of the long-double product differs
// from NCBI's two.
#[test]
fn pinned_double_rounding_erfc_and_spouge_calls_match_ncbi_double_bits() {
    let hex = |field: &str| u64::from_str_radix(field, 16).unwrap();
    let erfc = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/spouge_x87_20261027/erfc_double_rounding.tsv"
    ));
    for row in erfc.lines() {
        let f: Vec<_> = row.split('\t').collect();
        let z = f64::from_bits(hex(f[0]));
        assert_eq!(erfc_ncbi(z).to_bits(), hex(f[1]), "ErfC({z:.17e})");
    }
    let (spec, kbp) = scaled_blosum62_11_1();
    let spouge = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/spouge_x87_20261027/spouge_double_rounding.tsv"
    ));
    for row in spouge.lines() {
        let f: Vec<_> = row.split('\t').collect();
        let gbp = lookup_protein_gumbel_params(&spec, f[3].parse().unwrap()).unwrap();
        let actual = blast_spouge_stoe(
            f[0].parse().unwrap(),
            &kbp,
            &gbp,
            f[1].parse().unwrap(),
            f[2].parse().unwrap(),
        );
        let want = hex(f[4]);
        // The ErfC port above is plain double and double-double arithmetic, so its rows
        // are exact on every platform. `blast_spouge_stoe` also calls the platform `exp`
        // three times (the two `exp(-0.5 * f * f)` terms and `exp(-lambda * y)`), and the
        // E value multiplies them into a product. Apple's libm differs from glibc's by
        // an ULP per call (PD-LOSAT-PLATFORM-LIBM): platform-checks run 41062695644 read
        // 2 ULP on macOS arm64 (row `54521 99 256`) and 4 ULP on macOS x64 (row
        // `24337 600 112`). macOS accepts up to 8 ULP, every other OS stays exact.
        let max_ulp = if cfg!(target_os = "macos") { 8 } else { 0 };
        assert!(
            actual.to_bits().abs_diff(want) <= max_ulp,
            "{row}: bits {} want {want} (at most {max_ulp} ULP apart)",
            actual.to_bits()
        );
    }
    assert_eq!((erfc.lines().count(), spouge.lines().count()), (337, 177));
}

#[test]
fn test_lookup_protein_gumbel_params_blosum62_11_1() {
    let spec = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: 11,
        gap_extend: 1,
    };
    let gbp = lookup_protein_gumbel_params(&spec, 66738).expect("BLOSUM62 11/1 must exist");
    assert!((gbp.lambda - 0.267).abs() < 1e-12);
    assert!((gbp.c - 0.669720).abs() < 1e-12);
    assert!((gbp.a - 1.9).abs() < 1e-12);
    assert!((gbp.alpha - 42.602800).abs() < 1e-12);
    assert!((gbp.sigma - 43.636200).abs() < 1e-12);
    assert_eq!(gbp.db_length, 66738);
}

// The two Gumbel blocks this file built from literals before it read the table (dev
// fb8b10ee6), with the same arithmetic.
fn literal_gumbel(
    gap_open: i32,
    gap_extend: i32,
    [lambda, c, a, alpha, sigma, a_un, alpha_un]: [f64; 7],
) -> BlastGumbelBlk {
    let g = (gap_open + gap_extend) as f64;
    BlastGumbelBlk {
        lambda,
        c,
        g,
        a,
        alpha,
        sigma,
        a_un,
        alpha_un,
        b: 2.0 * g * (a_un - a),
        beta: 2.0 * g * (alpha_un - alpha),
        tau: 2.0 * g * (alpha_un - sigma),
        db_length: 1234,
    }
}

fn gumbel_bits(g: &BlastGumbelBlk) -> [u64; 11] {
    [
        g.lambda, g.c, g.g, g.a, g.alpha, g.sigma, g.a_un, g.alpha_un, g.b, g.beta, g.tau,
    ]
    .map(f64::to_bits)
}

#[test]
fn gumbel_blocks_from_the_table_equal_the_former_literals_bit_for_bit() {
    for (matrix, gap_open, gap_extend, literals) in [
        (
            ScoringMatrix::Blosum45,
            14,
            2,
            [0.195, 0.685753, 1.9, 60.736200, 61.102300, 0.9113, 9.611060],
        ),
        (
            ScoringMatrix::Blosum62,
            11,
            1,
            [0.267, 0.669720, 1.9, 42.602800, 43.636200, 0.7916, 4.964660],
        ),
    ] {
        let spec = ProteinScoringSpec {
            matrix,
            gap_open,
            gap_extend,
        };
        let table = lookup_protein_gumbel_params(&spec, 1234).expect("in the table");
        let literal = literal_gumbel(gap_open, gap_extend, literals);
        assert_eq!(gumbel_bits(&table), gumbel_bits(&literal), "{matrix:?}");
        assert_eq!(table.db_length, 1234);
    }
}

#[test]
fn every_gap_cost_row_has_a_gumbel_block() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3724-3749 fills the block for any row of the table, e.g. PAM30
    // 9/1 {9, 1, INT2_MAX, 0.294, 0.11, 0.61, 0.48, -6, 0.492716, 2.463920, 3.186150}
    // and the ungapped row's a 0.1938 and Alpha 0.161818.
    let pam30 = lookup_protein_gumbel_params(
        &ProteinScoringSpec {
            matrix: ScoringMatrix::Pam30,
            gap_open: 9,
            gap_extend: 1,
        },
        0,
    )
    .expect("PAM30 9/1");
    assert_eq!(
        (
            pam30.lambda,
            pam30.c,
            pam30.g,
            pam30.a,
            pam30.alpha,
            pam30.sigma
        ),
        (0.294, 0.492716, 10.0, 0.48, 2.463920, 3.186150)
    );
    assert_eq!((pam30.a_un, pam30.alpha_un), (0.1938, 0.161818));
    let mut rows = 0;
    for matrix in ScoringMatrix::ALL {
        for row in &crate::stats::protein_options::scoring_matrix_table(matrix).values[1..] {
            let spec = ProteinScoringSpec {
                matrix,
                gap_open: row.0[0] as i32,
                gap_extend: row.0[1] as i32,
            };
            assert!(lookup_protein_gumbel_params(&spec, 0).is_some(), "{spec:?}");
            rows += 1;
        }
    }
    assert_eq!(rows, 89);
    assert!(lookup_protein_gumbel_params(
        &ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 5,
            gap_extend: 5,
        },
        0,
    )
    .is_none());
}

#[test]
fn test_blast_spouge_stoe_matches_ncbi_formula_for_scaled_lambda() {
    let gbp = BlastGumbelBlk {
        lambda: 0.267,
        c: 0.669720,
        g: 12.0,
        a: 1.9,
        alpha: 42.602800,
        sigma: 43.636200,
        a_un: 0.7916,
        alpha_un: 4.964660,
        b: 2.0 * 12.0 * (0.7916 - 1.9),
        beta: 2.0 * 12.0 * (4.964660 - 42.602800),
        tau: 2.0 * 12.0 * (4.964660 - 43.636200),
        db_length: 5000,
    };
    let kbp = KarlinParams {
        lambda: 0.267 / 32.0,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };

    let evalue = blast_spouge_stoe(1062, &kbp, &gbp, 100, 200);
    assert!((evalue - 1.5863473416482987).abs() < 1e-12);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:5176-5231
#[test]
fn test_blast_spouge_stoe_matches_ncbi_formula_for_representative_tail_inputs() {
    let gbp = BlastGumbelBlk {
        lambda: 0.267,
        c: 0.669720,
        g: 12.0,
        a: 1.9,
        alpha: 42.602800,
        sigma: 43.636200,
        a_un: 0.7916,
        alpha_un: 4.964660,
        b: 2.0 * 12.0 * (0.7916 - 1.9),
        beta: 2.0 * 12.0 * (4.964660 - 42.602800),
        tau: 2.0 * 12.0 * (4.964660 - 43.636200),
        db_length: 5000,
    };

    let moderate_kbp = KarlinParams {
        lambda: 0.267 / 16.0,
        k: 0.082,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    let moderate = blast_spouge_stoe(128, &moderate_kbp, &gbp, 64, 96);
    assert!((moderate - 4250.4559852614275).abs() < 1e-9);

    let long_tail_kbp = KarlinParams {
        lambda: 0.267 / 64.0,
        k: 0.021,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    let long_tail = blast_spouge_stoe(2048, &long_tail_kbp, &gbp, 250, 400);
    assert!((long_tail - 4.058690669998566).abs() < 1e-12);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:5236-5268
#[test]
fn test_blast_spouge_etos_matches_ncbi_formula_for_representative_evalues() {
    let gbp = BlastGumbelBlk {
        lambda: 0.267,
        c: 0.669720,
        g: 12.0,
        a: 1.9,
        alpha: 42.602800,
        sigma: 43.636200,
        a_un: 0.7916,
        alpha_un: 4.964660,
        b: 2.0 * 12.0 * (0.7916 - 1.9),
        beta: 2.0 * 12.0 * (4.964660 - 42.602800),
        tau: 2.0 * 12.0 * (4.964660 - 43.636200),
        db_length: 5000,
    };

    let default_kbp = KarlinParams {
        lambda: 0.267 / 32.0,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    assert_eq!(blast_spouge_etos(1e-6, &default_kbp, &gbp, 100, 200), 2537);

    let tight_kbp = KarlinParams {
        lambda: 0.267 / 16.0,
        k: 0.082,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    assert_eq!(blast_spouge_etos(1e-30, &tight_kbp, &gbp, 64, 96), 3749);
}
