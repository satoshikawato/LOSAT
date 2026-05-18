//! DUST filter implementation for masking low-complexity regions in nucleotide sequences.
//!
//! This implements the symmetric DUST algorithm as described in NCBI BLAST.
//! The algorithm identifies low-complexity regions by calculating triplet-based
//! complexity scores within sliding windows.
//!
//! Reference: NCBI BLAST symdust.cpp/hpp

use std::collections::VecDeque;

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/util/random_gen.cpp:241-264
// ```c
// static const CRandom::TValue sm_State[CRandom::kStateSize] = {
//     0xd53f1852,  0xdfc78b83,  0x4f256096,  0xe643df7,
//     ...
//     0x6e37bd55
// };
// m_RJ = kStateOffset;
// m_RK = kStateSize - 1;
// ```
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/util/random_gen.hpp:217-240
// ```c
// r = m_State[m_RK] + m_State[m_RJ--];
// m_State[m_RK--] = r;
// ...
// return x_GetRand32Bits() >> 1;
// ```
struct NcbiLfgRandom {
    state: [u32; 33],
    rj: i32,
    rk: i32,
}

impl NcbiLfgRandom {
    fn new() -> Self {
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
        let level = if level >= 2 && level <= 64 { level } else { 20 };
        let window = if window >= 8 && window <= 64 {
            window
        } else {
            64
        };
        let linker = if linker >= 1 && linker <= 32 {
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

    /// Check if a position is within this masked interval
    #[inline]
    pub fn contains(&self, pos: usize) -> bool {
        pos >= self.start && pos < self.end
    }
}

/// Perfect interval - represents a region that exceeds the complexity threshold
#[derive(Debug, Clone)]
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
    #[allow(dead_code)]
    level: u32,
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
            level: params.level,
            window: params.window,
            linker: params.linker,
            low_k,
            thresholds,
        }
    }

    /// Create a DUST masker with default parameters
    pub fn with_defaults() -> Self {
        Self::new(20, 64, 1)
    }

    /// Convert a nucleotide base to 2-bit encoding (NCBI2NA format)
    /// A=0, C=1, G=2, T/U=3
    #[inline]
    fn encode_base(base: u8) -> Option<u8> {
        match base {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' | b'U' | b'u' => Some(3),
            _ => None, // Ambiguous bases
        }
    }

    #[inline]
    fn convert_iupac_to_ncbi2na(base: u8, rng: &mut NcbiLfgRandom) -> u8 {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/dustmask/symdust.hpp:75-84
        // ```c
        // switch( r )
        // {
        //     case 67: return 1;
        //     case 71: return 2;
        //     case 84: return 3;
        //     case 78: return (m_Random.GetRand() & 0x3);
        //     default: return 0;
        // }
        // ```
        match base {
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' | b'U' | b'u' => 3,
            b'N' | b'n' => (rng.next() & 0x3) as u8,
            _ => 0,
        }
    }

    #[inline]
    #[allow(dead_code)]
    fn encode_triplet(b1: u8, b2: u8, b3: u8) -> Option<u8> {
        let e1 = Self::encode_base(b1)?;
        let e2 = Self::encode_base(b2)?;
        let e3 = Self::encode_base(b3)?;
        Some((e1 << 4) | (e2 << 2) | e3)
    }

    /// Mask a sequence and return the list of masked intervals
    pub fn mask_sequence(&self, seq: &[u8]) -> Vec<MaskedInterval> {
        if seq.len() < 3 {
            return Vec::new();
        }
        self.mask_subsequence(seq, 0, seq.len())
    }

    /// Mask a subsequence and return the list of masked intervals
    pub fn mask_subsequence(&self, seq: &[u8], start: usize, stop: usize) -> Vec<MaskedInterval> {
        let mut result = Vec::new();

        if seq.is_empty() {
            return result;
        }

        let stop = stop.min(seq.len());
        let start = start.min(stop);

        // Need at least 3 bases for one triplet
        if stop <= start + 2 {
            return result;
        }

        let mut current_start = start;
        let mut rng = NcbiLfgRandom::new();

        while stop > current_start + 2 {
            // Initialize perfect intervals list for this window
            let mut perfect_list: VecDeque<PerfectInterval> = VecDeque::new();

            // Create triplet window tracker
            let mut window = TripletWindow::new(self.window, self.low_k, &self.thresholds);

            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/dustmask/symdust.cpp:235-249
            // ```c
            // seq_citer_type it(seq, start);
            // char c1 = *it, c2 = *++it;
            // triplet_type t = (converter_( c1 )<<2) + converter_( c2 );
            // it.SetPos(start + w.stop() + 2);
            // ...
            // t = ((t<<2)&TRIPLET_MASK) + (converter_( *it )&0x3);
            // ++it;
            // ```
            let mut current_triplet =
                (Self::convert_iupac_to_ncbi2na(seq[current_start], &mut rng) << 2)
                    + Self::convert_iupac_to_ncbi2na(seq[current_start + 1], &mut rng);
            let mut pos = current_start + 2;
            let mut done = false;

            while !done && pos < stop {
                // Save masked regions from previous window position
                self.save_masked_regions(
                    &mut result,
                    window.start(),
                    current_start,
                    &mut perfect_list,
                );

                // Shift window by adding new triplet
                let new_triplet = ((current_triplet << 2) & 0x3F)
                    + Self::convert_iupac_to_ncbi2na(seq[pos], &mut rng);
                current_triplet = new_triplet;
                pos += 1;

                if window.shift_window(new_triplet, &mut perfect_list) {
                    if window.needs_processing() {
                        window.find_perfect(&mut perfect_list);
                    }
                } else {
                    // Window contains only one triplet value - fast path
                    while pos < stop {
                        self.save_masked_regions(
                            &mut result,
                            window.start(),
                            current_start,
                            &mut perfect_list,
                        );

                        let new_triplet = ((current_triplet << 2) & 0x3F)
                            + Self::convert_iupac_to_ncbi2na(seq[pos], &mut rng);
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
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/dustmask/symdust.cpp:185-208
                // ```c
                // TMaskedInterval b = P.back().bounds_;
                // if( b.first < wstart ) {
                //     TMaskedInterval b1( b.first + start, b.second + start );
                //     ...
                //     if( s + linker_ >= b1.first ) {
                //         res.back().second = max( s, b1.second );
                //     }
                // ```
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

/// Triplet window tracker for DUST algorithm
/// Uses VecDeque to match NCBI BLAST's deque with push_front/pop_back semantics
// NCBI reference: ncbi-blast/c++/src/algo/dustmask/symdust.cpp:40-45
// ```c
// CSymDustMasker::triplets::triplets(
//     size_type window, Uint1 low_k,
//     perfect_list_type & perfect_list, thres_table_type & thresholds )
//     : ... P( perfect_list ), thresholds_( thresholds ) ...
// ```
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
        // NCBI: off = triplet_list_.size() - (L - start_) - 1
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

    fn find_perfect(&mut self, perfect_list: &mut VecDeque<PerfectInterval>) {
        let suffix_len = self.stop - self.l;

        if suffix_len >= self.triplet_list.len() {
            return;
        }

        let mut counts = self.c_v;
        let mut score = self.r_v;
        let mut max_perfect_score = 0u32;
        let mut max_len = 0usize;

        // NCBI: pos = L - 1, count starts at suffix_len and increments each iteration
        // it = triplet_list_.begin() + count (starts at suffix_len index)
        let mut pos = self.l.saturating_sub(1);
        let mut perfect_idx = 0usize;
        let mut count = suffix_len; // This is the candidate interval length variable

        // Iterate from suffix_len to end of triplet_list
        for idx in suffix_len..self.triplet_list.len() {
            let triplet = self.triplet_list[idx];
            let cnt = counts[triplet as usize];
            Self::add_triplet(&mut score, &mut counts, triplet);

            // Use count for threshold lookup (NCBI: thresholds_[count])
            if cnt > 0 && count < self.thresholds.len() && score * 10 > self.thresholds[count] {
                while perfect_idx < perfect_list.len() && pos <= perfect_list[perfect_idx].start {
                    let p = &perfect_list[perfect_idx];
                    if max_perfect_score == 0
                        || max_len * p.score as usize > max_perfect_score as usize * p.len
                    {
                        max_perfect_score = p.score;
                        max_len = p.len;
                    }
                    perfect_idx += 1;
                }

                if max_perfect_score == 0
                    || score as usize * max_len >= max_perfect_score as usize * count
                {
                    max_perfect_score = score;
                    max_len = count;
                    // NCBI reference: ncbi-blast/c++/src/algo/dustmask/symdust.cpp:156-163
                    // ```c
                    // if( max_perfect_score == 0 || score*max_len >= max_perfect_score*count ) {
                    //     max_perfect_score = score;
                    //     max_len = count;
                    //     perfect_iter = P.insert(
                    //             perfect_iter, perfect( pos, stop_ + 1,
                    //             max_perfect_score, count ) );
                    // }
                    // ```
                    let interval =
                        PerfectInterval::new(pos, self.stop + 1, max_perfect_score, count);
                    if perfect_idx < perfect_list.len() {
                        perfect_list.insert(perfect_idx, interval);
                    } else {
                        perfect_list.push_back(interval);
                    }
                }
            }

            // Increment count each iteration (NCBI: ++count in for-loop header)
            count += 1;
            if pos > 0 {
                pos -= 1;
            }
        }
    }
}

/// Check if a position is within any masked interval
pub fn is_position_masked(intervals: &[MaskedInterval], pos: usize) -> bool {
    intervals.iter().any(|interval| interval.contains(pos))
}

/// Check if a k-mer starting at position overlaps with any masked interval
pub fn is_kmer_masked(intervals: &[MaskedInterval], start: usize, kmer_len: usize) -> bool {
    let end = start + kmer_len;
    intervals.iter().any(|interval| {
        // Check if [start, end) overlaps with [interval.start, interval.end)
        start < interval.end && end > interval.start
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_base() {
        assert_eq!(DustMasker::encode_base(b'A'), Some(0));
        assert_eq!(DustMasker::encode_base(b'a'), Some(0));
        assert_eq!(DustMasker::encode_base(b'C'), Some(1));
        assert_eq!(DustMasker::encode_base(b'G'), Some(2));
        assert_eq!(DustMasker::encode_base(b'T'), Some(3));
        assert_eq!(DustMasker::encode_base(b'U'), Some(3));
        assert_eq!(DustMasker::encode_base(b'N'), None);
    }

    #[test]
    fn test_encode_triplet() {
        // AAA = 0b000000 = 0
        assert_eq!(DustMasker::encode_triplet(b'A', b'A', b'A'), Some(0));
        // TTT = 0b111111 = 63
        assert_eq!(DustMasker::encode_triplet(b'T', b'T', b'T'), Some(63));
        // ACG = 0b000110 = 6
        assert_eq!(DustMasker::encode_triplet(b'A', b'C', b'G'), Some(6));
    }

    #[test]
    fn test_simple_repeat() {
        let masker = DustMasker::with_defaults();

        // Simple repeat sequence should be masked
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let intervals = masker.mask_sequence(seq);

        // Should have at least one masked interval
        assert!(!intervals.is_empty(), "Poly-A sequence should be masked");
    }

    #[test]
    fn test_complex_sequence() {
        let masker = DustMasker::with_defaults();

        // Use a de Bruijn sequence B(4,3) over A/C/G/T: each 3-mer appears exactly once.
        // This is a more appropriate \"high complexity\" control than a short periodic repeat.
        fn debruijn_acgt_k3() -> Vec<u8> {
            fn db(t: usize, p: usize, k: usize, n: usize, a: &mut [usize], seq: &mut Vec<usize>) {
                if t > n {
                    if n % p == 0 {
                        for i in 1..=p {
                            seq.push(a[i]);
                        }
                    }
                } else {
                    a[t] = a[t - p];
                    db(t + 1, p, k, n, a, seq);
                    for j in (a[t - p] + 1)..k {
                        a[t] = j;
                        db(t + 1, t, k, n, a, seq);
                    }
                }
            }

            let alphabet: [u8; 4] = [b'A', b'C', b'G', b'T'];
            let k = alphabet.len();
            let n = 3usize;
            let mut a = vec![0usize; k * n + 1];
            let mut idx_seq: Vec<usize> = Vec::new();
            db(1, 1, k, n, &mut a, &mut idx_seq);

            let mut out: Vec<u8> = idx_seq.into_iter().map(|i| alphabet[i]).collect();
            // Linearize the cyclic de Bruijn sequence by appending the first n-1 symbols.
            let prefix: Vec<u8> = out[..(n - 1)].to_vec();
            out.extend_from_slice(&prefix);
            out
        }

        let seq = debruijn_acgt_k3();
        let intervals = masker.mask_sequence(&seq);

        // Complex sequence should have few or no masked regions
        let total_masked: usize = intervals.iter().map(|i| i.end - i.start).sum();
        assert!(
            total_masked < seq.len() / 2,
            "Complex sequence should not be heavily masked"
        );
    }

    #[test]
    fn test_short_sequence() {
        let masker = DustMasker::with_defaults();

        // Very short sequences should return empty
        let seq = b"AC";
        let intervals = masker.mask_sequence(seq);
        assert!(intervals.is_empty());
    }

    #[test]
    fn test_is_kmer_masked() {
        let intervals = vec![MaskedInterval::new(10, 20), MaskedInterval::new(30, 40)];

        // K-mer completely before masked region
        assert!(!is_kmer_masked(&intervals, 0, 5));

        // K-mer overlapping start of masked region
        assert!(is_kmer_masked(&intervals, 8, 5));

        // K-mer completely within masked region
        assert!(is_kmer_masked(&intervals, 12, 5));

        // K-mer overlapping end of masked region
        assert!(is_kmer_masked(&intervals, 18, 5));

        // K-mer between masked regions
        assert!(!is_kmer_masked(&intervals, 22, 5));
    }

    #[test]
    fn test_params_validation() {
        // Test parameter validation
        let params = DustParams::new(1, 5, 0);
        assert_eq!(params.level, 20); // Should default to 20 (out of range)
        assert_eq!(params.window, 64); // Should default to 64 (out of range)
        assert_eq!(params.linker, 1); // Should default to 1 (out of range)

        let params = DustParams::new(30, 32, 16);
        assert_eq!(params.level, 30);
        assert_eq!(params.window, 32);
        assert_eq!(params.linker, 16);
    }
}
