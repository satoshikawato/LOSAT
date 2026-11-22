//! `std::shuffle` with a `std::mt19937`, as NCBI's makeblastdb runs it for `-randomize`: the
//! order in which `-input_type blastdb` copies the sequences of its input.
//!
//! The engine is the C++ standard's (ISO/IEC 14882, [rand.eng.mers] and [rand.predef]: `mt19937`
//! is the Mersenne Twister of 624 words of 32 bits, whose 10000th output after the default seed
//! 5489 is 4123659995). The standard leaves `std::shuffle` and `uniform_int_distribution` to the
//! C++ library. The walk here is the one observed of GCC 13's library, with which NCBI builds
//! makeblastdb 2.17.0 for Linux (the official x64 build: GCC 13.2; bioconda's: GCC 13.4), from
//! black-box records: a generator that writes each value it gives and an element whose swap
//! writes the two positions recorded every draw and every swap of the library's shuffle (line RG
//! item D3d, 2026-11-21: 3,567 shuffles of 0 to 1,000,000 elements under 14 seeds, 11.1 million
//! records, no difference), and makeblastdb itself, under a fixed clock, copied databases of 1 to
//! 70,001 sequences in the order this walk gives (39 copies with each build, no difference). It is
//! written from those records, not from the library's source.
//!
//! The walk swaps each position from the second on with a position at or before it, Fisher and
//! Yates' shuffle, and takes the positions from the engine's 32-bit values as below:
//! * A value below `s` is the high 32 bits of `draw * s`; while the low 32 bits are below
//!   `(2^32 - s) mod s` the next draw is taken instead (Lemire's method: D. Lemire, "Fast Random
//!   Integer Generation in an Interval", ACM TOMACS 29(1), 2019).
//! * While `n * n` fits in 32 bits (`n` elements, at most 65,535), one value gives two positions:
//!   when `n - 1` is odd, the second element first swaps with a value below 2; then each pair of
//!   positions `i`, `i + 1` takes a value `r` below `(i + 1) * (i + 2)`, and `i` swaps with
//!   `r / (i + 2)`, then `i + 1` with `r % (i + 2)`.
//! * Otherwise each position `i` swaps with a value below `i + 1`.
//!
//! NCBI's builds on Windows (MSVC's library) and macOS (libc++) shuffle with their own libraries,
//! whose order is not measured: LOSAT gives the Linux builds' order on every platform.

/// `std::mt19937`.
pub(crate) struct Mt19937 {
    state: [u32; 624],
    /// The word of `state` that gives the next output; 624 when the words are used up.
    next: usize,
}

impl Mt19937 {
    /// The engine seeded with `seed` (`mt19937` constructed with `seed`; a wider seed is taken
    /// modulo 2^32, as the standard's seeding does).
    pub(crate) fn new(seed: u32) -> Self {
        let mut state = [0u32; 624];
        state[0] = seed;
        for i in 1..624 {
            let previous = state[i - 1];
            state[i] = 1_812_433_253u32
                .wrapping_mul(previous ^ (previous >> 30))
                .wrapping_add(i as u32);
        }
        Self { state, next: 624 }
    }

    /// The next output.
    pub(crate) fn next_u32(&mut self) -> u32 {
        if self.next == 624 {
            self.twist();
        }
        let mut y = self.state[self.next];
        self.next += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// The next 624 words of the state.
    fn twist(&mut self) {
        for k in 0..624 {
            let y = (self.state[k] & 0x8000_0000) | (self.state[(k + 1) % 624] & 0x7fff_ffff);
            let mut word = self.state[(k + 397) % 624] ^ (y >> 1);
            if y & 1 == 1 {
                word ^= 0x9908_b0df;
            }
            self.state[k] = word;
        }
        self.next = 0;
    }

    /// A value below `s` (1 to 2^32).
    fn below(&mut self, s: u64) -> u64 {
        let mut product = u64::from(self.next_u32()) * s;
        if (product & 0xffff_ffff) < s {
            let threshold = ((1u64 << 32) - s) % s;
            while (product & 0xffff_ffff) < threshold {
                product = u64::from(self.next_u32()) * s;
            }
        }
        product >> 32
    }
}

/// `std::shuffle(items.begin(), items.end(), engine)`.
pub(crate) fn shuffle<T>(items: &mut [T], engine: &mut Mt19937) {
    let n = items.len() as u64;
    if n < 2 {
        return;
    }
    if n * n <= u64::from(u32::MAX) {
        let mut i = 1u64;
        if (n - 1) % 2 == 1 {
            items.swap(1, engine.below(2) as usize);
            i = 2;
        }
        while i < n {
            let r = engine.below((i + 1) * (i + 2));
            items.swap(i as usize, (r / (i + 2)) as usize);
            items.swap(i as usize + 1, (r % (i + 2)) as usize);
            i += 2;
        }
    } else {
        for i in 1..n {
            items.swap(i as usize, engine.below(i + 1) as usize);
        }
    }
}

#[cfg(test)]
#[path = "std_shuffle_tests.rs"]
mod tests;
