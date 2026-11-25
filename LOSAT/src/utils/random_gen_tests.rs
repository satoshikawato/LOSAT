// NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:217-280
// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:215-264,287-309
use super::*;

// Pinned draws from NCBI's `CRandom` (`eGetRand_LFG`, `Reset` startup state), covering both the
// `m_RJ` wrap at draw 13 and the `m_RK` wrap at draw 33.
const DEFAULT_RAW_32: [u32; 40] = [
    0x13bc_f39b,
    0xdbc8_a221,
    0x5998_680b,
    0x7e21_67bb,
    0x54f0_c660,
    0xeda2_b322,
    0x0244_57df,
    0xcb0c_f05d,
    0xc3e2_41a9,
    0x3477_8fcb,
    0xeca0_1abc,
    0x4711_0c0a,
    0xacc7_efb7,
    0xd978_0b83,
    0x4b17_d62f,
    0xa4e3_688d,
    0xdc78_5742,
    0x058d_ad71,
    0xfb50_e401,
    0xe959_cec1,
    0x7092_26a3,
    0x33a4_8732,
    0xcb53_9ca9,
    0x99fa_9971,
    0x7618_bf3e,
    0xd954_a601,
    0xaf62_0b2d,
    0x1291_2429,
    0x27a6_c24c,
    0xeadc_9539,
    0x54b3_0e07,
    0xdb18_6f84,
    0xbe98_e713,
    0x844f_1a3e,
    0x0f6d_2953,
    0x24ec_04b4,
    0x181c_012c,
    0xcb09_859e,
    0xc6f7_5923,
    0xb1a6_630c,
];

const DEFAULT_GET_RAND: [u32; 40] = [
    0x09de_79cd,
    0x6de4_5110,
    0x2ccc_3405,
    0x3f10_b3dd,
    0x2a78_6330,
    0x76d1_5991,
    0x0122_2bef,
    0x6586_782e,
    0x61f1_20d4,
    0x1a3b_c7e5,
    0x7650_0d5e,
    0x2388_8605,
    0x5663_f7db,
    0x6cbc_05c1,
    0x258b_eb17,
    0x5271_b446,
    0x6e3c_2ba1,
    0x02c6_d6b8,
    0x7da8_7200,
    0x74ac_e760,
    0x3849_1351,
    0x19d2_4399,
    0x65a9_ce54,
    0x4cfd_4cb8,
    0x3b0c_5f9f,
    0x6caa_5300,
    0x57b1_0596,
    0x0948_9214,
    0x13d3_6126,
    0x756e_4a9c,
    0x2a59_8703,
    0x6d8c_37c2,
    0x5f4c_7389,
    0x4227_8d1f,
    0x07b6_94a9,
    0x1276_025a,
    0x0c0e_0096,
    0x6584_c2cf,
    0x637b_ac91,
    0x58d3_3186,
];

// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:241-264
#[test]
fn default_matches_ncbi_reset_sequence() {
    let mut raw = NcbiRandom::default();
    let mut draws = NcbiRandom::default();
    for (step, (&expected_raw, &expected_draw)) in
        DEFAULT_RAW_32.iter().zip(&DEFAULT_GET_RAND).enumerate()
    {
        let got_raw = raw.rand_32_bits();
        let got_draw = draws.get_rand();
        assert_eq!(got_raw, expected_raw, "default raw mismatch at step {step}");
        assert_eq!(
            got_draw, expected_draw,
            "default draw mismatch at step {step}"
        );
        assert_eq!(got_draw, got_raw >> 1);
        assert!(got_draw <= 0x7fff_ffff, "draw exceeded CRandom::GetMax()");
    }
}

// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:227-231,287-309
#[test]
fn seeded_generator_matches_ncbi_set_seed_draws() {
    // Pinned draws from NCBI's seeded `CRandom` (first 8 draws, plus draw 33 after one full
    // state cycle and draw 100).
    let cases: &[(u32, [u32; 8], u32, u32)] = &[
        (
            0,
            [
                0x1006_04f3,
                0x5171_7c3d,
                0x0fb6_c7ae,
                0x41a5_a561,
                0x0db4_2601,
                0x390c_c2e8,
                0x598b_ba2c,
                0x24a8_f31c,
            ],
            0x598f_86ec,
            0x7a28_00b8,
        ),
        (
            1,
            [
                0x0467_aa76,
                0x3bee_ff9f,
                0x1923_cd41,
                0x2d15_b818,
                0x6375_f40d,
                0x09fd_e49a,
                0x3826_ca7c,
                0x799f_d2ac,
            ],
            0x5b5b_d788,
            0x38ab_6ace,
        ),
        (
            42,
            [
                0x280b_2c86,
                0x4a09_0a3c,
                0x1b99_b1e1,
                0x6208_b752,
                0x1f7f_f40e,
                0x009c_4a1c,
                0x5efc_674c,
                0x1529_a0bc,
            ],
            0x2514_c084,
            0x3bb7_683f,
        ),
        (
            150,
            [
                0x413d_0000,
                0x36fc_775e,
                0x1598_0c1b,
                0x3550_9c50,
                0x4d42_e154,
                0x2656_8134,
                0x4867_490c,
                0x6d4f_f17c,
            ],
            0x6746_c254,
            0x1b28_2951,
        ),
        (
            1_000,
            [
                0x2d74_8e9f,
                0x4bba_b119,
                0x6194_8fda,
                0x6f6e_be45,
                0x0ac1_06d5,
                0x66f8_6238,
                0x6743_72ac,
                0x0902_3d9c,
            ],
            0x5faa_684c,
            0x2b7e_64b4,
        ),
        (
            0x7fff_ffff,
            [
                0x5ba4_5f6f,
                0x26f3_f8dc,
                0x4649_c21a,
                0x1635_92ab,
                0x77f2_57f4,
                0x681b_a136,
                0x7af0_a9dc,
                0x4fb2_138c,
            ],
            0x57c3_3650,
            0x7ba4_96a3,
        ),
        (
            0x8000_0000,
            [
                0x5006_04f3,
                0x1171_7c3d,
                0x4fb6_c7ae,
                0x01a5_a561,
                0x4db4_2601,
                0x390c_c2e8,
                0x598b_ba2c,
                0x24a8_f31c,
            ],
            0x598f_86ec,
            0x3a28_00b8,
        ),
        (
            u32::MAX,
            [
                0x1ba4_5f6f,
                0x66f3_f8dc,
                0x0649_c21a,
                0x5635_92ab,
                0x37f2_57f4,
                0x681b_a136,
                0x7af0_a9dc,
                0x4fb2_138c,
            ],
            0x57c3_3650,
            0x3ba4_96a3,
        ),
    ];

    for &(seed, expected_head, expected_33, expected_100) in cases {
        let mut rng = NcbiRandom::new(seed);
        for (step, &expected) in expected_head.iter().enumerate() {
            let draw = rng.get_rand();
            assert_eq!(draw, expected, "seed {seed} mismatch at step {step}");
            assert!(draw <= 0x7fff_ffff);
        }
        for _ in 8..33 {
            assert!(rng.get_rand() <= 0x7fff_ffff);
        }
        assert_eq!(
            rng.get_rand(),
            expected_33,
            "seed {seed} mismatch at step 33"
        );
        for _ in 34..100 {
            assert!(rng.get_rand() <= 0x7fff_ffff);
        }
        assert_eq!(
            rng.get_rand(),
            expected_100,
            "seed {seed} mismatch at step 100"
        );
    }
}

// NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:251-280
#[test]
fn get_rand_index_and_range_match_ncbi() {
    let mut pow2 = NcbiRandom::new(150);
    let got_pow2: Vec<u32> = (0..8).map(|_| pow2.get_rand_index(16)).collect();
    assert_eq!(got_pow2, [8, 6, 2, 6, 9, 4, 9, 13]);

    let mut non_pow2 = NcbiRandom::new(150);
    let got_non_pow2: Vec<u32> = (0..8).map(|_| non_pow2.get_rand_index(7)).collect();
    assert_eq!(got_non_pow2, [3, 3, 2, 1, 6, 2, 2, 1]);

    // Large non-power-of-two size exercises the rejection loop of `GetRandIndex`.
    let mut reject = NcbiRandom::new(150);
    let got_reject: Vec<u32> = (0..8)
        .map(|_| reject.get_rand_index(3_000_000_000))
        .collect();
    assert_eq!(
        got_reject,
        [
            2_189_033_472,
            1_845_030_589,
            724_572_214,
            1_788_950_689,
            2_592_457_384,
            1_286_406_760,
            2_429_456_920,
            2_629_917_463,
        ]
    );

    let mut range = NcbiRandom::new(150);
    let got_range: Vec<u32> = (0..8).map(|_| range.get_rand_range(10, 25)).collect();
    assert_eq!(got_range, [18, 16, 12, 16, 19, 14, 19, 23]);
}

#[test]
fn cloned_generator_produces_identical_continuation() {
    let mut first = NcbiRandom::new(42);
    for _ in 0..37 {
        first.get_rand();
    }
    let mut second = first.clone();
    for step in 0..100 {
        assert_eq!(
            first.get_rand(),
            second.get_rand(),
            "clone diverged at step {step}"
        );
    }
}
