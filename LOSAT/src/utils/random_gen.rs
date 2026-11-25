//! NCBI's `CRandom` lagged Fibonacci random number generator (`eGetRand_LFG`).
//!
//! NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:51-175,217-241
//! NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:98,215-264,287-309

// NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:159-165
// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:98
const STATE_SIZE: usize = 33;
const STATE_OFFSET: i32 = 12;
const LAST_INDEX: i32 = (STATE_SIZE - 1) as i32;

// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:241-252
// Hard-coded startup state used by `CRandom::Reset` (the default constructor). This 33-word table
// is distinct from any `SetSeed` output, and copying it directly avoids warmup overhead.
const DEFAULT_STATE: [u32; STATE_SIZE] = [
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
];

/// NCBI's `CRandom` lagged Fibonacci generator (`eGetRand_LFG`, lags 33 and 13).
///
/// Use [`NcbiRandom::default`] for NCBI's default constructor (`Reset` startup state) and
/// [`NcbiRandom::new`] for a seeded generator (`SetSeed` linear congruential initialization and
/// 330-step warmup).
#[derive(Debug, Clone)]
pub(crate) struct NcbiRandom {
    state: [u32; STATE_SIZE],
    rj: i32,
    rk: i32,
}

impl Default for NcbiRandom {
    // NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:215-224,234-264
    #[inline]
    fn default() -> Self {
        Self {
            state: DEFAULT_STATE,
            rj: STATE_OFFSET,
            rk: LAST_INDEX,
        }
    }
}

impl NcbiRandom {
    // NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:227-231,287-309
    /// Constructs a generator seeded with `seed` (NCBI `CRandom` seeded constructor and
    /// `SetSeed`).
    pub(crate) fn new(seed: u32) -> Self {
        let mut state = [0u32; STATE_SIZE];
        state[0] = seed;
        for i in 1..STATE_SIZE {
            state[i] = 1_103_515_245u32
                .wrapping_mul(state[i - 1])
                .wrapping_add(12_345);
        }
        let mut random = Self {
            state,
            rj: STATE_OFFSET,
            rk: LAST_INDEX,
        };
        for _ in 0..10 * STATE_SIZE {
            random.get_rand();
        }
        random
    }

    // NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:238-241
    /// `GetRand()`: the next number in `[0, 0x7fffffff]`, discarding the lowest bit of
    /// `rand_32_bits`.
    #[inline]
    pub(crate) fn get_rand(&mut self) -> u32 {
        self.rand_32_bits() >> 1
    }

    // NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:217-235
    // `m_RJ` and `m_RK` start at 12 and 32 (20 apart modulo 33) and decrement together, wrapping
    // to 32 when negative, matching NCBI's signed `int` indices and avoiding extra branchless
    // register pressure when inlined into the DUST inner loop.
    /// `x_GetRand32Bits()`: the next 32-bit sum of the lagged Fibonacci generator.
    #[inline]
    pub(crate) fn rand_32_bits(&mut self) -> u32 {
        let r = self.state[self.rk as usize].wrapping_add(self.state[self.rj as usize]);
        self.state[self.rk as usize] = r;
        self.rj -= 1;
        self.rk -= 1;
        if self.rk < 0 {
            self.rk = LAST_INDEX;
        } else if self.rj < 0 {
            self.rj = LAST_INDEX;
        }
        r
    }
}

#[cfg(test)]
#[path = "random_gen_tests.rs"]
mod tests;
