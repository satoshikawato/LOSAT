// NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:217-241
use super::*;

/// Independent reference port of NCBI's signed-index `CRandom` (`eGetRand_LFG`).
struct ReferenceLfg {
    state: [u32; STATE_SIZE],
    rj: i32,
    rk: i32,
}

impl ReferenceLfg {
    fn from_default_table() -> Self {
        Self {
            state: [
                0xd53f_1852,
                0xdfc7_8b83,
                0x4f25_6096,
                0x0e64_3df7,
                0x82c3_59bf,
                0xc779_4dfa,
                0xd5e9_ffaa,
                0x2c8c_b64a,
                0x2f07_b334,
                0xad5a_7eb5,
                0x96dc_0cde,
                0x6fc2_4589,
                0xa585_3646,
                0xe715_76e2,
                0x0dae_30df,
                0xb09c_e711,
                0x5e56_ef87,
                0x4b4b_0082,
                0x6f4f_340e,
                0xc5bb_17e8,
                0xd788_d765,
                0x6749_8087,
                0x9d7a_ba26,
                0x2613_51d4,
                0x411e_e7ea,
                0x0393_a263,
                0x2c5a_5835,
                0xc115_fcd8,
                0x25e9_132c,
                0xd0c6_e906,
                0xc2bc_5b2d,
                0x6c06_5c98,
                0x6e37_bd55,
            ],
            rj: 12,
            rk: 32,
        }
    }

    fn seeded(seed: u32) -> Self {
        let mut state = [0u32; STATE_SIZE];
        state[0] = seed;
        for i in 1..STATE_SIZE {
            state[i] = state[i - 1]
                .wrapping_mul(1_103_515_245)
                .wrapping_add(12_345);
        }
        let mut rng = Self {
            state,
            rj: 12,
            rk: 32,
        };
        for _ in 0..330 {
            rng.get_rand();
        }
        rng
    }

    fn rand_32_bits(&mut self) -> u32 {
        let r = self.state[self.rk as usize].wrapping_add(self.state[self.rj as usize]);
        self.state[self.rk as usize] = r;
        self.rj -= 1;
        self.rk -= 1;
        if self.rk < 0 {
            self.rk = 32;
        } else if self.rj < 0 {
            self.rj = 32;
        }
        r
    }

    fn get_rand(&mut self) -> u32 {
        self.rand_32_bits() >> 1
    }
}

// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:241-264
#[test]
fn default_state_matches_ncbi_reset_table_and_draws() {
    let mut unified = NcbiRandom::default();
    let mut reference = ReferenceLfg::from_default_table();
    assert_eq!(unified.state, reference.state);
    assert_eq!(unified.rj, reference.rj);
    assert_eq!(unified.rk, reference.rk);
    for step in 0..256 {
        let draw = unified.get_rand();
        assert_eq!(
            draw,
            reference.get_rand(),
            "default draw mismatch at {step}"
        );
        assert!(draw <= 0x7fff_ffff, "draw exceeded CRandom::GetMax()");
    }
    let mut raw_unified = NcbiRandom::default();
    let mut raw_ref = ReferenceLfg::from_default_table();
    for step in 0..128 {
        assert_eq!(
            raw_unified.rand_32_bits(),
            raw_ref.rand_32_bits(),
            "raw 32-bit draw mismatch at {step}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:241-264,287-309
// Reversing the 330-step LFG warmup is a bijection on the 33-word state, so any state created by
// `SetSeed` must invert back to the linear congruential recurrence; `DEFAULT_STATE` does not,
// proving NCBI's hard-coded `Reset` table cannot be produced by `SetSeed` for any `u32` seed.
#[test]
fn default_state_is_distinct_from_every_seeded_state() {
    let mut pre_warmup = DEFAULT_STATE;
    let mut rj = STATE_OFFSET as usize;
    let mut rk = LAST_INDEX as usize;
    for _ in 0..10 * STATE_SIZE {
        rj = (rj + 1) % STATE_SIZE;
        rk = (rk + 1) % STATE_SIZE;
        pre_warmup[rk] = pre_warmup[rk].wrapping_sub(pre_warmup[rj]);
    }
    let expected_second = pre_warmup[0]
        .wrapping_mul(1_103_515_245)
        .wrapping_add(12_345);
    assert_ne!(pre_warmup[1], expected_second);
    assert_ne!(
        NcbiRandom::default().get_rand(),
        NcbiRandom::new(0).get_rand()
    );
    assert_ne!(
        NcbiRandom::default().get_rand(),
        NcbiRandom::new(1).get_rand()
    );
}

// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:227-231,287-309
#[test]
fn seeded_generator_matches_reference_across_boundary_seeds() {
    for seed in [0u32, 1, 42, 1_000, 0x7fff_ffff, 0x8000_0000, u32::MAX] {
        let mut unified = NcbiRandom::new(seed);
        let mut reference = ReferenceLfg::seeded(seed);
        assert_eq!(
            unified.state, reference.state,
            "state mismatch for seed {seed}"
        );
        for step in 0..200 {
            let draw = unified.get_rand();
            assert_eq!(
                draw,
                reference.get_rand(),
                "draw mismatch for seed {seed} at step {step}"
            );
            assert!(draw <= 0x7fff_ffff);
        }
    }
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
