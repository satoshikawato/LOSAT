//! NCBI's `CRandom` lagged Fibonacci random number generator (`eGetRand_LFG`).
//!
//! NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:51-175,217-280
//! NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:98,215-264,287-309

// NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:159-165
// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:98
const STATE_SIZE: usize = 33;
const STATE_OFFSET: i32 = 12;
const LAST_INDEX: i32 = (STATE_SIZE - 1) as i32;

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
    // Build the 33-word `CRandom::Reset` startup table inline as an array literal rather than
    // copying from a module-level `const` so x86-64 codegen emits inline stores instead of a
    // `memcpy` call on per-sequence masker construction.
    #[inline]
    fn default() -> Self {
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
        // NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:303-308
        // The 330 draws of `SetSeed` on the local array: each draw adds the word at `m_RJ` to the
        // word at `m_RK`, and both indices step down by one with `m_RJ` 20 below `m_RK` (mod 33).
        // One pass of `rk` from 32 down to 0 is 33 draws (`m_RJ` from 12 down to 0, then from 32
        // down to 13), and after 10 passes the indices are back at 12 and 32.
        for _ in 0..10 {
            for rk in (20..STATE_SIZE).rev() {
                state[rk] = state[rk].wrapping_add(state[rk - 20]);
            }
            for rk in (0..20).rev() {
                state[rk] = state[rk].wrapping_add(state[rk + 13]);
            }
        }
        Self {
            state,
            rj: STATE_OFFSET,
            rk: LAST_INDEX,
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:238-241
    /// `GetRand()`: the next number in `[0, 0x7fffffff]`, discarding the lowest bit of
    /// `rand_32_bits`.
    #[inline]
    pub(crate) fn get_rand(&mut self) -> u32 {
        self.rand_32_bits() >> 1
    }

    // NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:251-264
    // The multiplication of a power of two keeps the high bits; another size rejects the draws
    // past the last whole multiple of it (the subtraction wraps in 32 bits, as NCBI's).
    /// `GetRandIndex` of `size`: a number in `[0, size)`.
    fn get_rand_index(&mut self, size: u32) -> u32 {
        if size & size.wrapping_sub(1) == 0 {
            return ((u64::from(self.rand_32_bits()) * u64::from(size)) >> 32) as u32;
        }
        loop {
            let bits = self.rand_32_bits();
            let r = bits % size;
            if bits <= r.wrapping_sub(size) {
                return r;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:277-280
    /// `GetRand` of a range: a number in `[min_value, max_value]`.
    pub(crate) fn get_rand_range(&mut self, min_value: u32, max_value: u32) -> u32 {
        min_value
            .wrapping_add(self.get_rand_index(max_value.wrapping_sub(min_value).wrapping_add(1)))
    }

    // NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:217-235
    // Uses NCBI's signed `int` index representation (`m_RJ` starts at 12, `m_RK` at 32; they are
    // 20 apart modulo 33 and wrap to 32 when negative). Matches dev's DUST inner-loop instruction
    // count in Callgrind.
    /// `x_GetRand32Bits()`: the next 32-bit sum of the lagged Fibonacci generator.
    #[inline]
    fn rand_32_bits(&mut self) -> u32 {
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
