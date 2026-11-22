//! DUST filter implementation for masking low-complexity regions in nucleotide sequences.
//!
//! This implements the symmetric DUST algorithm as described in NCBI BLAST.
//! The algorithm identifies low-complexity regions by calculating triplet-based
//! complexity scores within sliding windows.
//!
//! NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp, c++/include/algo/dustmask/symdust.hpp

use std::collections::VecDeque;

// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:241-264
// NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:217-240
/// `CRandom` of `CIupac2Ncbi2na_converter` (the random numbers that resolve `N`): one per
/// `CSymDustMasker`, so the calls of one masker draw from one sequence of numbers.
pub(crate) struct NcbiLfgRandom {
    state: [u32; 33],
    rj: i32,
    rk: i32,
}

impl NcbiLfgRandom {
    pub(crate) fn new() -> Self {
        Self {
            state: [
                0xd53f1852, 0xdfc78b83, 0x4f256096, 0x0e643df7, 0x82c359bf, 0xc7794dfa, 0xd5e9ffaa,
                0x2c8cb64a, 0x2f07b334, 0xad5a7eb5, 0x96dc0cde, 0x6fc24589, 0xa5853646, 0xe71576e2,
                0x0dae30df, 0xb09ce711, 0x5e56ef87, 0x4b4b0082, 0x6f4f340e, 0xc5bb17e8, 0xd788d765,
                0x67498087, 0x9d7aba26, 0x261351d4, 0x411ee7ea, 0x0393a263, 0x2c5a5835, 0xc115fcd8,
                0x25e9132c, 0xd0c6e906, 0xc2bc5b2d, 0x6c065c98, 0x6e37bd55,
            ],
            rj: 12,
            rk: 32,
        }
    }

    #[inline]
    fn next(&mut self) -> u32 {
        let r = self.state[self.rk as usize].wrapping_add(self.state[self.rj as usize]);
        self.state[self.rk as usize] = r;
        self.rj -= 1;
        self.rk -= 1;
        if self.rk < 0 {
            self.rk = 32;
        } else if self.rj < 0 {
            self.rj = 32;
        }
        r >> 1
    }
}

/// DUST filter parameters
#[derive(Debug, Clone)]
pub struct DustParams {
    /// Score threshold (default: 20, valid range: 2-64)
    pub level: u32,
    /// Maximum window size (default: 64, valid range: 8-64)
    pub window: usize,
    /// Maximum distance to merge consecutive masked intervals (default: 1, valid range: 1-32)
    pub linker: usize,
}

impl Default for DustParams {
    fn default() -> Self {
        Self {
            level: 20,
            window: 64,
            linker: 1,
        }
    }
}

impl DustParams {
    pub fn new(level: u32, window: usize, linker: usize) -> Self {
        // Validate and clamp parameters to NCBI BLAST ranges
        let level = if (2..=64).contains(&level) { level } else { 20 };
        let window = if (8..=64).contains(&window) {
            window
        } else {
            64
        };
        let linker = if (1..=32).contains(&linker) {
            linker
        } else {
            1
        };
        Self {
            level,
            window,
            linker,
        }
    }
}

/// Represents a masked interval in a sequence (0-based, inclusive start, exclusive end)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaskedInterval {
    pub start: usize,
    pub end: usize,
}

impl MaskedInterval {
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }
}

/// Perfect interval - represents a region that exceeds the complexity threshold
#[derive(Debug, Clone, Copy)]
struct PerfectInterval {
    start: usize,
    end: usize,
    score: u32,
    len: usize,
}

impl PerfectInterval {
    fn new(start: usize, end: usize, score: u32, len: usize) -> Self {
        Self {
            start,
            end,
            score,
            len,
        }
    }
}

/// DUST masker implementation following NCBI BLAST's symmetric DUST algorithm
pub struct DustMasker {
    window: usize,
    linker: usize,
    low_k: u8,
    thresholds: Vec<u32>,
}

impl DustMasker {
    /// Create a new DUST masker with the given parameters
    pub fn new(level: u32, window: usize, linker: usize) -> Self {
        let params = DustParams::new(level, window, linker);

        // low_k: max triplet multiplicity that guarantees window score is not above threshold
        let low_k = (params.level / 5) as u8;

        // Build threshold table: thresholds[i] = i * level for i in 1..window-2
        // thresholds[0] = 1 (special case)
        let mut thresholds = Vec::with_capacity(params.window - 2);
        thresholds.push(1);
        for i in 1..(params.window - 2) {
            thresholds.push(i as u32 * params.level);
        }

        Self {
            window: params.window,
            linker: params.linker,
            low_k,
            thresholds,
        }
    }

    /// Create a DUST masker with default parameters
    #[cfg(test)]
    pub fn with_defaults() -> Self {
        Self::new(20, 64, 1)
    }

    #[inline]
    fn convert_iupac_to_ncbi2na(base: u8, rng: &mut NcbiLfgRandom) -> u8 {
        // NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:75-84
        match base {
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' | b'U' | b'u' => 3,
            b'N' | b'n' => (rng.next() & 0x3) as u8,
            _ => 0,
        }
    }

    /// `CSymDustMasker::operator()` of `seq`: the masked intervals of the whole sequence, with a
    /// new masker's random numbers (blastn's `Blast_FindDustFilterLoc` makes one masker per query).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:286-290
    pub fn mask_sequence(&self, seq: &[u8]) -> Vec<MaskedInterval> {
        if seq.is_empty() {
            return Vec::new();
        }
        self.mask_range(seq, 0, seq.len() - 1, &mut NcbiLfgRandom::new())
    }

    /// `CSymDustMasker::operator()` of `seq`, `start` and `stop`: the masked intervals of the
    /// positions `start..=stop` (NCBI's `stop` is the last position), drawing the random numbers of
    /// `N` from `rng`, the masker's converter (dustmasker calls one masker over the pieces of a
    /// sequence between runs of `N`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:214-229
    /// `mask_range_reference` ports this function. The ring-buffer version (was
    /// LOSAT_X_DUSTRING) runs the same steps and only needs a window of at most `X_RING` (64,
    /// the largest NCBI window); a larger window uses the reference.
    ///
    /// Inlined into each caller with `x_mask_range`, so blastn's `mask_sequence` gets a copy
    /// with `start` 0, `stop` the last position and its own random state. One shared copy for
    /// blastn and dustmasker kept those values live in the inner loop, which then moved
    /// registers on every letter: blastn's DUST ran 7.6 % more instructions than dev's
    /// (callgrind of ba17fdf17, megablast, 12 MB of cDNA).
    #[inline(always)]
    pub(crate) fn mask_range(
        &self,
        seq: &[u8],
        start: usize,
        stop: usize,
        rng: &mut NcbiLfgRandom,
    ) -> Vec<MaskedInterval> {
        if seq.is_empty() {
            return Vec::new();
        }
        let stop = stop.min(seq.len() - 1);
        let start = start.min(stop);
        if self.window <= X_RING {
            self.x_mask_range(seq, start, stop, rng)
        } else {
            self.mask_range_reference(seq, start, stop, rng)
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:229-282
    // This is `CSymDustMasker::operator()` with the window in `XTripletWindow` and the base
    // codes from `X_BASE_CODE`. The bases are converted in the same order, so a random
    // number for an `N` is drawn at the same point. It gives the same intervals.
    // `mask_range_reference` step for step, with the window in `XTripletWindow` and the
    // base codes from a table.
    #[inline(always)]
    fn x_mask_range(
        &self,
        seq: &[u8],
        start: usize,
        stop: usize,
        rng: &mut NcbiLfgRandom,
    ) -> Vec<MaskedInterval> {
        // NCBI's loop test (`GetPos` of the iterator is at most `stop`) as `pos < end` over the letters up to `stop`: the
        // compiler then knows that each `seq[pos]` is in bounds.
        let seq = &seq[..=stop];
        let end = seq.len();
        let mut result = Vec::new();
        let mut current_start = start;
        // NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:75-84
        // 4 marks the bases that draw a random number.
        let code = |base: u8, rng: &mut NcbiLfgRandom| -> u8 {
            let c = X_BASE_CODE[base as usize];
            if c < 4 {
                c
            } else {
                (rng.next() & 0x3) as u8
            }
        };

        while stop > current_start + 2 {
            let mut perfect_list: VecDeque<PerfectInterval> = VecDeque::new();
            let mut window = XTripletWindow::new(self.window, self.low_k, &self.thresholds);

            let mut current_triplet =
                (code(seq[current_start], rng) << 2) + code(seq[current_start + 1], rng);
            let mut pos = current_start + 2;
            let mut done = false;

            while !done && pos < end {
                if !perfect_list.is_empty() {
                    self.save_masked_regions(
                        &mut result,
                        window.start,
                        current_start,
                        &mut perfect_list,
                    );
                }

                let new_triplet = ((current_triplet << 2) & 0x3F) + code(seq[pos], rng);
                current_triplet = new_triplet;
                pos += 1;

                if window.shift_window(new_triplet, &mut perfect_list) {
                    if window.needs_processing() {
                        window.find_perfect(&mut perfect_list);
                    }
                } else {
                    while pos < end {
                        if !perfect_list.is_empty() {
                            self.save_masked_regions(
                                &mut result,
                                window.start,
                                current_start,
                                &mut perfect_list,
                            );
                        }

                        let new_triplet = ((current_triplet << 2) & 0x3F) + code(seq[pos], rng);
                        current_triplet = new_triplet;

                        if window.shift_window(new_triplet, &mut perfect_list) {
                            done = true;
                            break;
                        }
                        pos += 1;
                    }
                }
            }

            let mut wstart = window.start;
            while !perfect_list.is_empty() {
                self.save_masked_regions(&mut result, wstart, current_start, &mut perfect_list);
                wstart += 1;
            }

            if window.start > 0 {
                current_start += window.start;
            } else {
                break;
            }
        }

        result
    }

    fn mask_range_reference(
        &self,
        seq: &[u8],
        start: usize,
        stop: usize,
        rng: &mut NcbiLfgRandom,
    ) -> Vec<MaskedInterval> {
        let mut result = Vec::new();
        let mut current_start = start;

        while stop > current_start + 2 {
            // Initialize perfect intervals list for this window
            let mut perfect_list: VecDeque<PerfectInterval> = VecDeque::new();

            // Create triplet window tracker
            let mut window = TripletWindow::new(self.window, self.low_k, &self.thresholds);

            // NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:235-249
            let mut current_triplet = (Self::convert_iupac_to_ncbi2na(seq[current_start], rng)
                << 2)
                + Self::convert_iupac_to_ncbi2na(seq[current_start + 1], rng);
            let mut pos = current_start + 2;
            let mut done = false;

            while !done && pos <= stop {
                // Save masked regions from previous window position
                self.save_masked_regions(
                    &mut result,
                    window.start(),
                    current_start,
                    &mut perfect_list,
                );

                // Shift window by adding new triplet
                let new_triplet =
                    ((current_triplet << 2) & 0x3F) + Self::convert_iupac_to_ncbi2na(seq[pos], rng);
                current_triplet = new_triplet;
                pos += 1;

                if window.shift_window(new_triplet, &mut perfect_list) {
                    if window.needs_processing() {
                        window.find_perfect(&mut perfect_list);
                    }
                } else {
                    // Window contains only one triplet value - fast path
                    while pos <= stop {
                        self.save_masked_regions(
                            &mut result,
                            window.start(),
                            current_start,
                            &mut perfect_list,
                        );

                        let new_triplet = ((current_triplet << 2) & 0x3F)
                            + Self::convert_iupac_to_ncbi2na(seq[pos], rng);
                        current_triplet = new_triplet;

                        if window.shift_window(new_triplet, &mut perfect_list) {
                            done = true;
                            break;
                        }
                        pos += 1;
                    }
                }
            }

            // Append remaining perfect intervals to result
            let mut wstart = window.start();
            while !perfect_list.is_empty() {
                self.save_masked_regions(&mut result, wstart, current_start, &mut perfect_list);
                wstart += 1;
            }

            // Move to next segment
            if window.start() > 0 {
                current_start += window.start();
            } else {
                break;
            }
        }

        result
    }

    /// Save masked regions from perfect intervals
    fn save_masked_regions(
        &self,
        result: &mut Vec<MaskedInterval>,
        wstart: usize,
        offset: usize,
        perfect_list: &mut VecDeque<PerfectInterval>,
    ) {
        if perfect_list.is_empty() {
            return;
        }

        // Get the last (oldest) perfect interval
        if let Some(p) = perfect_list.back() {
            if p.start < wstart {
                let interval_start = p.start + offset;
                // NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:185-208
                // NCBI stores b1.second as an inclusive CSeq_loc endpoint
                // (dust_filter.cpp:101-112 passes GetTo() through). LOSAT's
                // `MaskedInterval` is 0-based half-open, so add one here.
                let interval_end = p.end + offset + 1;

                // Try to merge with previous interval if within linker distance
                if let Some(last) = result.last_mut() {
                    if last.end.saturating_sub(1) + self.linker >= interval_start {
                        last.end = last.end.max(interval_end);
                    } else {
                        result.push(MaskedInterval::new(interval_start, interval_end));
                    }
                } else {
                    result.push(MaskedInterval::new(interval_start, interval_end));
                }

                // Remove processed perfect intervals
                while let Some(p) = perfect_list.back() {
                    if p.start < wstart {
                        perfect_list.pop_back();
                    } else {
                        break;
                    }
                }
            }
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:43,173
/// NCBI accepts windows up to 64, so a window holds at most 62 triplets. 64 slots with
/// index arithmetic modulo 64 are enough.
/// Capacity of the ring buffer (the largest DUST window).
const X_RING: usize = 64;

/// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:75-84
/// The same mapping as `convert_iupac_to_ncbi2na` in a 256-entry table (it also lists the
/// lower case letters and `U` that function accepts). Entry 4 means that the caller
/// draws a random number, as NCBI does for `N`.
/// `convert_iupac_to_ncbi2na` as a table; 4 = draws a random number.
static X_BASE_CODE: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t[b'N' as usize] = 4;
    t[b'n' as usize] = 4;
    t
};

/// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:40-48
/// The same state as NCBI's `triplets` class. Only the container for `triplet_list_`
/// differs (a fixed ring instead of a `std::deque`). The perfect list `P` is passed to
/// each method instead of being held by reference.
/// `TripletWindow` (was LOSAT_X_DUSTRING) with the triplet deque in a
/// fixed ring buffer. Deque index `i` (0 = newest, as after `push_front`) is
/// `ring[(head + i) % X_RING]`; every method below is the `TripletWindow`
/// method of the same name with only that substitution.
struct XTripletWindow<'a> {
    ring: [u8; X_RING],
    head: usize,
    len: usize,
    start: usize,
    stop: usize,
    max_size: usize,
    low_k: u8,
    l: usize,
    c_w: [u8; 64],
    c_v: [u8; 64],
    r_w: u32,
    r_v: u32,
    num_diff: u32,
    thresholds: &'a [u32],
    /// No NCBI counterpart: the vector `find_perfect` merges the new perfect intervals into,
    /// kept with the window between calls (emptied each time) instead of made for each call.
    x_merged: Vec<PerfectInterval>,
}

impl<'a> XTripletWindow<'a> {
    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:43-48
    /// Starts with the same values, with an empty ring in place of the empty deque.
    fn new(window: usize, low_k: u8, thresholds: &'a [u32]) -> Self {
        Self {
            ring: [0; X_RING],
            head: 0,
            len: 0,
            start: 0,
            stop: 0,
            max_size: window - 2,
            low_k,
            l: 0,
            c_w: [0; 64],
            c_v: [0; 64],
            r_w: 0,
            r_v: 0,
            num_diff: 0,
            thresholds,
            x_merged: Vec::new(),
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:103
    /// `triplet_list_[index]` of the deque: index 0 is the newest triplet.
    #[inline(always)]
    fn at(&self, index: usize) -> u8 {
        self.ring[(self.head + index) & (X_RING - 1)]
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:54-55
    /// Returns the oldest triplet and removes it.
    #[inline(always)]
    fn pop_back(&mut self) -> u8 {
        self.len -= 1;
        self.ring[(self.head + self.len) & (X_RING - 1)]
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:60
    /// Makes `triplet` the newest element (index 0).
    #[inline(always)]
    fn push_front(&mut self, triplet: u8) {
        self.head = (self.head + X_RING - 1) & (X_RING - 1);
        self.ring[self.head] = triplet;
        self.len += 1;
    }

    /// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:275-277
    /// The same two steps in the same order.
    #[inline(always)]
    fn add_triplet(sum: &mut u32, counts: &mut [u8; 64], triplet: u8) {
        let idx = (triplet & 63) as usize;
        *sum += counts[idx] as u32;
        counts[idx] += 1;
    }

    /// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:287-289
    /// The same two steps in the same order.
    #[inline(always)]
    fn rem_triplet(sum: &mut u32, counts: &mut [u8; 64], triplet: u8) {
        let idx = (triplet & 63) as usize;
        counts[idx] -= 1;
        *sum -= counts[idx] as u32;
    }

    /// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:251-256
    /// The same test. The three conditions are combined without short-circuit, which is
    /// safe because the table is read at a clamped index and that value is then ignored.
    /// The three tests of `TripletWindow::needs_processing` as one
    /// condition (the table is read at a clamped index when the second test
    /// fails, and that read is then ignored).
    #[inline(always)]
    fn needs_processing(&self) -> bool {
        let count = self.stop - self.l;
        let last = self.thresholds.len() - 1;
        (count < self.len)
            & (count < self.thresholds.len())
            & (10 * self.r_w > self.thresholds[count.min(last)])
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:76-114
    /// The same steps in the same order on the same integers. The two `if`s that depend on the
    /// sequence (whether `c_w[s]` is 0 and whether `L` is `start_`) are written as arithmetic on
    /// the 0/1 value of the condition, which gives the same values.
    #[inline(always)]
    fn shift_window(&mut self, triplet: u8, perfect_list: &mut VecDeque<PerfectInterval>) -> bool {
        if self.len >= self.max_size {
            if self.num_diff <= 1 {
                return self.shift_high(triplet, perfect_list);
            }

            // The two `if`s of the reference depend on the sequence in a way
            // a branch predictor cannot learn; they are written as arithmetic
            // on the 0/1 value of their conditions.
            let old = (self.pop_back() & 63) as usize;
            self.c_w[old] -= 1;
            self.r_w -= self.c_w[old] as u32;
            self.num_diff -= (self.c_w[old] == 0) as u32;

            // When the suffix start equals the window start, advance it and drop the expired
            // triplet from the suffix counts (`rem_triplet_info`).
            let at_start = (self.l == self.start) as usize;
            self.l += at_start;
            let suffix_count = self.c_v[old] - at_start as u8;
            self.c_v[old] = suffix_count;
            self.r_v -= suffix_count as u32 * at_start as u32;

            self.start += 1;
        }

        self.push_front(triplet);
        self.num_diff += (self.c_w[(triplet & 63) as usize] == 0) as u32;
        Self::add_triplet(&mut self.r_w, &mut self.c_w, triplet);
        Self::add_triplet(&mut self.r_v, &mut self.c_v, triplet);

        if self.c_v[(triplet & 63) as usize] > self.low_k {
            let mut off = self.len - (self.l - self.start) - 1;
            loop {
                let t = self.at(off);
                Self::rem_triplet(&mut self.r_v, &mut self.c_v, t);
                self.l += 1;
                if t == triplet {
                    break;
                }
                if off == 0 {
                    break;
                }
                off -= 1;
            }
        }

        self.stop += 1;

        if self.len >= self.max_size && self.num_diff <= 1 {
            perfect_list.clear();
            perfect_list.push_front(PerfectInterval::new(self.start, self.stop + 1, 0, 0));
            return false;
        }

        true
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:52-69
    /// The same steps in the same order.
    #[inline(never)]
    fn shift_high(&mut self, triplet: u8, perfect_list: &mut VecDeque<PerfectInterval>) -> bool {
        let old_triplet = self.pop_back();
        Self::rem_triplet(&mut self.r_w, &mut self.c_w, old_triplet);
        if self.c_w[(old_triplet & 63) as usize] == 0 {
            self.num_diff -= 1;
        }
        self.start += 1;

        self.push_front(triplet);
        if self.c_w[(triplet & 63) as usize] == 0 {
            self.num_diff += 1;
        }
        Self::add_triplet(&mut self.r_w, &mut self.c_w, triplet);
        self.stop += 1;

        if self.num_diff <= 1 {
            perfect_list.push_front(PerfectInterval::new(self.start, self.stop + 1, 0, 0));
            return false;
        }

        true
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:123-141,161-163
    /// The one-pass merge of `TripletWindow::find_perfect` (see there), reading the
    /// triplets from the ring.
    #[inline(never)]
    fn find_perfect(&mut self, perfect_list: &mut VecDeque<PerfectInterval>) {
        let suffix_len = self.stop - self.l;
        if suffix_len >= self.len {
            return;
        }
        let mut counts = self.c_v;
        let mut score = self.r_v;
        let mut max_perfect_score = 0u32;
        let mut max_len = 0usize;
        let mut pos = self.l.saturating_sub(1);
        let mut count = suffix_len;
        let mut read = 0usize;
        let mut merged = std::mem::take(&mut self.x_merged);
        merged.clear();
        let mut inserted = false;
        let old_len = perfect_list.len();
        for idx in suffix_len..self.len {
            let triplet = self.at(idx);
            let cnt = counts[(triplet & 63) as usize];
            Self::add_triplet(&mut score, &mut counts, triplet);
            if cnt > 0 && count < self.thresholds.len() && score * 10 > self.thresholds[count] {
                while read < old_len && pos <= perfect_list[read].start {
                    let p = perfect_list[read];
                    if max_perfect_score == 0
                        || max_len * p.score as usize > max_perfect_score as usize * p.len
                    {
                        max_perfect_score = p.score;
                        max_len = p.len;
                    }
                    if inserted {
                        merged.push(p);
                    }
                    read += 1;
                }
                if max_perfect_score == 0
                    || score as usize * max_len >= max_perfect_score as usize * count
                {
                    max_perfect_score = score;
                    max_len = count;
                    if !inserted {
                        inserted = true;
                        merged.reserve(old_len + 8);
                        merged.extend(perfect_list.iter().take(read).copied());
                    }
                    merged.push(PerfectInterval::new(
                        pos,
                        self.stop + 1,
                        max_perfect_score,
                        count,
                    ));
                }
            }
            count += 1;
            pos = pos.saturating_sub(1);
        }
        if inserted {
            merged.extend(perfect_list.iter().skip(read).copied());
            perfect_list.clear();
            perfect_list.extend(merged.drain(..));
        }
        self.x_merged = merged;
    }
}

/// Triplet window tracker for DUST algorithm
/// Uses VecDeque to match NCBI BLAST's deque with push_front/pop_back semantics
// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:40-45
struct TripletWindow<'a> {
    triplet_list: VecDeque<u8>,
    start: usize,
    stop: usize,
    max_size: usize,
    low_k: u8,
    l: usize,      // suffix start position (L in NCBI code)
    c_w: [u8; 64], // triplet counts for whole window
    c_v: [u8; 64], // triplet counts for suffix
    r_w: u32,      // running sum for whole window
    r_v: u32,      // running sum for suffix
    num_diff: u32,
    thresholds: &'a [u32],
}

impl<'a> TripletWindow<'a> {
    fn new(window: usize, low_k: u8, thresholds: &'a [u32]) -> Self {
        Self {
            triplet_list: VecDeque::with_capacity(window),
            start: 0,
            stop: 0,
            max_size: window - 2,
            low_k,
            l: 0, // suffix start position
            c_w: [0; 64],
            c_v: [0; 64],
            r_w: 0,
            r_v: 0,
            num_diff: 0,
            thresholds,
        }
    }

    fn start(&self) -> usize {
        self.start
    }

    #[inline]
    fn add_triplet(sum: &mut u32, counts: &mut [u8; 64], triplet: u8) {
        let idx = triplet as usize;
        *sum += counts[idx] as u32;
        counts[idx] += 1;
    }

    #[inline]
    fn rem_triplet(sum: &mut u32, counts: &mut [u8; 64], triplet: u8) {
        let idx = triplet as usize;
        counts[idx] -= 1;
        *sum -= counts[idx] as u32;
    }

    fn needs_processing(&self) -> bool {
        let count = self.stop - self.l;
        if count >= self.triplet_list.len() {
            return false;
        }
        if count >= self.thresholds.len() {
            return false;
        }
        10 * self.r_w > self.thresholds[count]
    }

    fn shift_window(&mut self, triplet: u8, perfect_list: &mut VecDeque<PerfectInterval>) -> bool {
        if self.triplet_list.len() >= self.max_size {
            if self.num_diff <= 1 {
                return self.shift_high(triplet, perfect_list);
            }

            // Remove oldest triplet from back (NCBI: pop_back)
            let old_triplet = self.triplet_list.pop_back().unwrap();
            Self::rem_triplet(&mut self.r_w, &mut self.c_w, old_triplet);
            if self.c_w[old_triplet as usize] == 0 {
                self.num_diff -= 1;
            }

            if self.l == self.start {
                self.l += 1;
                Self::rem_triplet(&mut self.r_v, &mut self.c_v, old_triplet);
            }

            self.start += 1;
        }

        // Add new triplet at front (NCBI: push_front)
        self.triplet_list.push_front(triplet);
        if self.c_w[triplet as usize] == 0 {
            self.num_diff += 1;
        }
        Self::add_triplet(&mut self.r_w, &mut self.c_w, triplet);
        Self::add_triplet(&mut self.r_v, &mut self.c_v, triplet);

        // Update suffix start if triplet count exceeds low_k
        // Compute the deque offset of suffix start L from `triplet_list_` length and `start_`.
        // With push_front/pop_back, index 0 is newest, back() is oldest
        // So off maps position L to deque index near the back
        if self.c_v[triplet as usize] > self.low_k {
            let mut off = self.triplet_list.len() - (self.l - self.start) - 1;
            loop {
                let t = self.triplet_list[off];
                Self::rem_triplet(&mut self.r_v, &mut self.c_v, t);
                self.l += 1;
                if t == triplet {
                    break;
                }
                if off == 0 {
                    break;
                }
                off -= 1;
            }
        }

        self.stop += 1;

        if self.triplet_list.len() >= self.max_size && self.num_diff <= 1 {
            perfect_list.clear();
            perfect_list.push_front(PerfectInterval::new(self.start, self.stop + 1, 0, 0));
            return false;
        }

        true
    }

    fn shift_high(&mut self, triplet: u8, perfect_list: &mut VecDeque<PerfectInterval>) -> bool {
        // Remove oldest triplet from back (NCBI: pop_back)
        let old_triplet = self.triplet_list.pop_back().unwrap();
        Self::rem_triplet(&mut self.r_w, &mut self.c_w, old_triplet);
        if self.c_w[old_triplet as usize] == 0 {
            self.num_diff -= 1;
        }
        self.start += 1;

        // Add new triplet at front (NCBI: push_front)
        self.triplet_list.push_front(triplet);
        if self.c_w[triplet as usize] == 0 {
            self.num_diff += 1;
        }
        Self::add_triplet(&mut self.r_w, &mut self.c_w, triplet);
        self.stop += 1;

        if self.num_diff <= 1 {
            perfect_list.push_front(PerfectInterval::new(self.start, self.stop + 1, 0, 0));
            return false;
        }

        true
    }

    // NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:123-166
    // Same walk over the triplets with the same tests and the same integer updates. Only
    // the list update differs: NCBI inserts into a `std::list`; here the old list and the
    // new elements are merged into a new vector in the order the inserts would give. The
    // argument is in the comment below. (Was LOSAT_X_DUSTFAST; the literal std::list port
    // that it replaced was deleted when the merge graduated.)
    //
    // NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:137-166
    // NCBI's P is a std::list, so the insert is O(1); a VecDeque insert moves
    // every later element. Insert positions never move backwards during one
    // call, and stepping over an element this call inserted leaves
    // max_perfect_score/max_len unchanged (it compares the element with
    // itself), so the result is the old list with the new elements merged in
    // at the positions where the walk stood.
    fn find_perfect(&mut self, perfect_list: &mut VecDeque<PerfectInterval>) {
        let suffix_len = self.stop - self.l;
        if suffix_len >= self.triplet_list.len() {
            return;
        }
        let mut counts = self.c_v;
        let mut score = self.r_v;
        let mut max_perfect_score = 0u32;
        let mut max_len = 0usize;
        let mut pos = self.l.saturating_sub(1);
        let mut count = suffix_len;
        // `read` walks the list as it was on entry; `merged` is only filled
        // once something has been inserted.
        let mut read = 0usize;
        let mut merged: Vec<PerfectInterval> = Vec::new();
        let mut inserted = false;
        let old_len = perfect_list.len();
        for idx in suffix_len..self.triplet_list.len() {
            let triplet = self.triplet_list[idx];
            let cnt = counts[triplet as usize];
            Self::add_triplet(&mut score, &mut counts, triplet);
            if cnt > 0 && count < self.thresholds.len() && score * 10 > self.thresholds[count] {
                while read < old_len && pos <= perfect_list[read].start {
                    let p = perfect_list[read];
                    if max_perfect_score == 0
                        || max_len * p.score as usize > max_perfect_score as usize * p.len
                    {
                        max_perfect_score = p.score;
                        max_len = p.len;
                    }
                    if inserted {
                        merged.push(p);
                    }
                    read += 1;
                }
                if max_perfect_score == 0
                    || score as usize * max_len >= max_perfect_score as usize * count
                {
                    max_perfect_score = score;
                    max_len = count;
                    if !inserted {
                        inserted = true;
                        merged.reserve(old_len + 8);
                        merged.extend(perfect_list.iter().take(read).copied());
                    }
                    merged.push(PerfectInterval::new(
                        pos,
                        self.stop + 1,
                        max_perfect_score,
                        count,
                    ));
                }
            }
            count += 1;
            pos = pos.saturating_sub(1);
        }
        if inserted {
            merged.extend(perfect_list.iter().skip(read).copied());
            perfect_list.clear();
            perfect_list.extend(merged);
        }
    }
}

#[cfg(test)]
#[path = "dust_tests.rs"]
mod tests;
