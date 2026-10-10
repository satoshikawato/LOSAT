//! DUST filter implementation for masking low-complexity regions in nucleotide sequences.
//!
//! This implements the symmetric DUST algorithm as described in NCBI BLAST.
//! The algorithm identifies low-complexity regions by calculating triplet-based
//! complexity scores within sliding windows.
//!
//! Reference: NCBI BLAST symdust.cpp/hpp

use std::collections::VecDeque;

// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:136-141
// ```c
// for( ; it != iend; ++it, ++count, --pos ) {
//     Uint1 cnt = counts[*it];
//     add_triplet_info( score, counts, *it );
//
//     if( cnt > 0 && score*10 > thresholds_[count] ) {
// ```
// Reads the LOSAT_X_DUSTFAST switch once. `TripletWindow::find_perfect` uses it to choose
// between the port of this loop and the one-pass merge that builds the same list.
fn x_dust_fast() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_DUSTFAST").is_some())
}

/// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:243-253
/// ```c
///         while( !done && it.GetPos() <= stop )
///         {
///             save_masked_regions( *res.get(), w.start(), start );
/// ...
///             if( w.shift_window( t ) ) {
///                 if( w.needs_processing() ) {
///                     w.find_perfect();
/// ```
/// Reads the LOSAT_X_DUSTRING / LOSAT_X_DUSTSHADOW switches once. The mode chooses which
/// window type runs this loop: the ported `TripletWindow` (0), the ring-buffer
/// `XTripletWindow` (1), or both with a comparison (2).
/// EXPERIMENT: 0 = reference window, 1 = LOSAT_X_DUSTRING (the same steps on
/// a fixed ring buffer), 2 = LOSAT_X_DUSTSHADOW (both, compared).
fn x_dust_ring() -> u8 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_DUSTSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_DUSTRING").is_some() {
            1
        } else {
            0
        }
    })
}

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
        // NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:229-233,280-281
        // ```c
        // while( stop > 2 + start )    // there must be at least one triplet
        // {
        //     // initializations
        //     P.clear();
        //     triplets w( window_, low_k_, P, thresholds_ );
        // ...
        //     if( w.start() > 0 ) start += w.start();
        //     else break;
        // ```
        // Dispatch point of LOSAT_X_DUSTRING / LOSAT_X_DUSTSHADOW. `mask_subsequence_reference`
        // ports this function. The ring-buffer version runs the same steps and only needs a
        // window of at most `X_RING` (64, the largest NCBI window); a larger window uses the
        // reference. In mode 2 both run and `assert!` compares the intervals.
        match x_dust_ring() {
            1 if self.window <= X_RING => self.x_mask_subsequence(seq, start, stop),
            2 if self.window <= X_RING => {
                let fast = self.x_mask_subsequence(seq, start, stop);
                let reference = self.mask_subsequence_reference(seq, start, stop);
                assert!(
                    fast == reference,
                    "LOSAT_X_DUSTSHADOW: ring-buffer DUST differs from the reference"
                );
                reference
            }
            _ => self.mask_subsequence_reference(seq, start, stop),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:229-282
    // ```c
    // while( stop > 2 + start )    // there must be at least one triplet
    // ...
    //         char c1 = *it, c2 = *++it;
    //         triplet_type t = (converter_( c1 )<<2) + converter_( c2 );
    // ...
    //             t = ((t<<2)&TRIPLET_MASK) + (converter_( *it )&0x3);
    // ...
    //     // append the rest of the perfect intervals to the result
    // ...
    //     if( w.start() > 0 ) start += w.start();
    // ```
    // This is `CSymDustMasker::operator()` with the window in `XTripletWindow` and the base
    // codes from `X_BASE_CODE`. The bases are converted in the same order, so a random
    // number for an `N` is drawn at the same point. It gives the same intervals.
    // EXPERIMENT (LOSAT_X_DUSTRING): `mask_subsequence_reference` step for
    // step, with the window in `XTripletWindow` and the base codes from a
    // table.
    fn x_mask_subsequence(&self, seq: &[u8], start: usize, stop: usize) -> Vec<MaskedInterval> {
        let mut result = Vec::new();

        if seq.is_empty() {
            return result;
        }

        let stop = stop.min(seq.len());
        let start = start.min(stop);

        if stop <= start + 2 {
            return result;
        }

        let mut current_start = start;
        let mut rng = NcbiLfgRandom::new();
        // NCBI reference: ncbi-blast/c++/include/algo/dustmask/symdust.hpp:75-84
        // ```c
        // case 67: return 1;
        // case 71: return 2;
        // case 84: return 3;
        // case 78: return (m_Random.GetRand() & 0x3);
        // default: return 0;
        // ```
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
                (code(seq[current_start], &mut rng) << 2) + code(seq[current_start + 1], &mut rng);
            let mut pos = current_start + 2;
            let mut done = false;

            while !done && pos < stop {
                if !perfect_list.is_empty() {
                    self.save_masked_regions(
                        &mut result,
                        window.start,
                        current_start,
                        &mut perfect_list,
                    );
                }

                let new_triplet = ((current_triplet << 2) & 0x3F) + code(seq[pos], &mut rng);
                current_triplet = new_triplet;
                pos += 1;

                if window.shift_window(new_triplet, &mut perfect_list) {
                    if window.needs_processing() {
                        window.find_perfect(&mut perfect_list);
                    }
                } else {
                    while pos < stop {
                        if !perfect_list.is_empty() {
                            self.save_masked_regions(
                                &mut result,
                                window.start,
                                current_start,
                                &mut perfect_list,
                            );
                        }

                        let new_triplet =
                            ((current_triplet << 2) & 0x3F) + code(seq[pos], &mut rng);
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

    fn mask_subsequence_reference(
        &self,
        seq: &[u8],
        start: usize,
        stop: usize,
    ) -> Vec<MaskedInterval> {
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

/// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:43,173
/// ```c
///       max_size_( window - 2 ), low_k_( low_k ),
/// ...
///       window_( (window >= 8 && window <= 64) ? window : DEFAULT_WINDOW ),
/// ```
/// NCBI accepts windows up to 64, so a window holds at most 62 triplets. 64 slots with
/// index arithmetic modulo 64 are enough.
/// Capacity of the ring buffer (the largest DUST window).
const X_RING: usize = 64;

/// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:75-84
/// ```c
/// Uint1 operator()( Uint1 r )
/// {
///     switch( r )
///     {
///         case 67: return 1;
///         case 71: return 2;
///         case 84: return 3;
///         case 78: return (m_Random.GetRand() & 0x3);
///         default: return 0;
/// ```
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
/// ```c
/// CSymDustMasker::triplets::triplets(
///     size_type window, Uint1 low_k,
///     perfect_list_type & perfect_list, thres_table_type & thresholds )
///     : start_( 0 ), stop_( 0 ), max_size_( window - 2 ), low_k_( low_k ),
///       L( 0 ), P( perfect_list ), thresholds_( thresholds ),
///       r_w( 0 ), r_v( 0 ), num_diff( 0 )
/// {
///     std::fill( c_w, c_w + 64, 0 );
/// ```
/// The same state as NCBI's `triplets` class. Only the container for `triplet_list_`
/// differs (a fixed ring instead of a `std::deque`). The perfect list `P` is passed to
/// each method instead of being held by reference.
/// EXPERIMENT (LOSAT_X_DUSTRING): `TripletWindow` with the triplet deque in a
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
}

impl<'a> XTripletWindow<'a> {
    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:43-48
    /// ```c
    ///     : start_( 0 ), stop_( 0 ), max_size_( window - 2 ), low_k_( low_k ),
    ///       L( 0 ), P( perfect_list ), thresholds_( thresholds ),
    ///       r_w( 0 ), r_v( 0 ), num_diff( 0 )
    /// {
    ///     std::fill( c_w, c_w + 64, 0 );
    /// ```
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
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:103
    /// ```c
    ///             rem_triplet_info( r_v, c_v, triplet_list_[off] );
    /// ```
    /// `triplet_list_[index]` of the deque: index 0 is the newest triplet.
    #[inline(always)]
    fn at(&self, index: usize) -> u8 {
        self.ring[(self.head + index) & (X_RING - 1)]
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:54-55
    /// ```c
    ///     triplet_type s = triplet_list_.back();
    ///     triplet_list_.pop_back();
    /// ```
    /// Returns the oldest triplet and removes it.
    #[inline(always)]
    fn pop_back(&mut self) -> u8 {
        self.len -= 1;
        self.ring[(self.head + self.len) & (X_RING - 1)]
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:60
    /// ```c
    ///     triplet_list_.push_front( t );
    /// ```
    /// Makes `triplet` the newest element (index 0).
    #[inline(always)]
    fn push_front(&mut self, triplet: u8) {
        self.head = (self.head + X_RING - 1) & (X_RING - 1);
        self.ring[self.head] = triplet;
        self.len += 1;
    }

    /// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:275-277
    /// ```c
    /// void add_triplet_info(
    ///         Uint4 & r, counts_type & c, triplet_type t )
    /// { r += c[t]; ++c[t]; }
    /// ```
    /// The same two steps in the same order.
    #[inline(always)]
    fn add_triplet(sum: &mut u32, counts: &mut [u8; 64], triplet: u8) {
        let idx = (triplet & 63) as usize;
        *sum += counts[idx] as u32;
        counts[idx] += 1;
    }

    /// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:287-289
    /// ```c
    /// void rem_triplet_info(
    ///         Uint4 & r, counts_type & c, triplet_type t )
    /// { --c[t]; r -= c[t]; }
    /// ```
    /// The same two steps in the same order.
    #[inline(always)]
    fn rem_triplet(sum: &mut u32, counts: &mut [u8; 64], triplet: u8) {
        let idx = (triplet & 63) as usize;
        counts[idx] -= 1;
        *sum -= counts[idx] as u32;
    }

    /// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:251-256
    /// ```c
    /// bool needs_processing() const
    /// {
    ///   Uint4 count = stop_ - L;
    ///   return count < triplet_list_.size() &&
    ///          10*r_w > thresholds_[count];
    /// }
    /// ```
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
    /// ```c
    ///     if( triplet_list_.size() >= max_size_ ) {
    ///         if( num_diff <= 1 ) {
    ///             return shift_high( t );
    ///         }
    /// ...
    ///         if( L == start_ ) {
    ///             ++L;
    ///             rem_triplet_info( r_v, c_v, s );
    ///         }
    /// ...
    ///     if( c_v[t] > low_k_ ) {
    ///         Uint4 off = triplet_list_.size() - (L - start_) - 1;
    /// ...
    ///     if( triplet_list_.size() >= max_size_ && num_diff <= 1 ) {
    /// ```
    /// The same steps in the same order on the same integers. The two `if`s that depend on
    /// the sequence (`c_w[s] == 0` and `L == start_`) are written as arithmetic on the 0/1
    /// value of the condition, which gives the same values.
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

            // if (L == start) { ++L; rem_triplet_info(r_v, c_v, old); }
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
    /// ```c
    /// bool CSymDustMasker::triplets::shift_high( triplet_type t )
    /// {
    ///     triplet_type s = triplet_list_.back();
    ///     triplet_list_.pop_back();
    ///     rem_triplet_info( r_w, c_w, s );
    ///     if( c_w[s] == 0 ) --num_diff;
    ///     ++start_;
    /// ...
    ///     if( num_diff <= 1 ) {
    ///         P.insert( P.begin(), perfect( start_, stop_ + 1, 0, 0 ) );
    /// ```
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
    /// ```c
    ///     Uint4 count = stop_ - L; // count is the suffix length
    /// ...
    ///     Uint4 score = r_v; // and of the partial sum
    ///     perfect_iter_type perfect_iter = P.begin();
    /// ...
    ///     for( ; it != iend; ++it, ++count, --pos ) {
    ///         Uint1 cnt = counts[*it];
    /// ...
    ///                 perfect_iter = P.insert(
    ///                         perfect( pos, stop_ + 1,
    ///                         max_perfect_score, count ) );
    /// ```
    /// The one-pass merge of `TripletWindow::x_find_perfect_merge` (see there), reading the
    /// triplets from the ring.
    /// `TripletWindow::x_find_perfect_merge` (the `LOSAT_X_DUSTFAST` form of
    /// `find_perfect`, which builds the same list) over the ring buffer.
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
        let mut merged: Vec<PerfectInterval> = Vec::new();
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
            if pos > 0 {
                pos -= 1;
            }
        }
        if inserted {
            merged.extend(perfect_list.iter().skip(read).copied());
            perfect_list.clear();
            perfect_list.extend(merged);
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
        // NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:118,136-141,161-163
        // ```c
        // inline void CSymDustMasker::triplets::find_perfect()
        // ...
        // for( ; it != iend; ++it, ++count, --pos ) {
        //     Uint1 cnt = counts[*it];
        // ...
        //         perfect_iter = P.insert(
        //                 perfect( pos, stop_ + 1,
        //                 max_perfect_score, count ) );
        // ```
        // Dispatch point of LOSAT_X_DUSTFAST. `find_perfect_reference` ports this function. The
        // merge builds the same list of perfect intervals in one pass.
        if x_dust_fast() {
            return self.x_find_perfect_merge(perfect_list);
        }
        self.find_perfect_reference(perfect_list)
    }

    // NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:123-166
    // ```c
    // Uint4 count = stop_ - L; // count is the suffix length
    // ...
    // perfect_iter_type perfect_iter = P.begin();
    // ...
    // for( ; it != iend; ++it, ++count, --pos ) {
    //     Uint1 cnt = counts[*it];
    //     add_triplet_info( score, counts, *it );
    //     if( cnt > 0 && score*10 > thresholds_[count] ) {
    // ...
    //         if(    max_perfect_score == 0
    //             || score*max_len >= max_perfect_score*count ) {
    // ```
    // Same walk over the triplets with the same tests and the same integer updates. Only
    // the list update differs: NCBI inserts into a `std::list`; here the old list and the
    // new elements are merged into a new vector in the order the inserts would give. The
    // argument is in the comment below. LOSAT_X_DUSTSHADOW and the random test compare it
    // with the port.
    // EXPERIMENT (LOSAT_X_DUSTFAST): the same list, built by one merge pass.
    //
    // NCBI reference: ncbi-blast/c++/src/algo/dustmask/symdust.cpp:137-166
    // ```c
    // for( impl_citer_type it( triplet_list_.begin() + count ),
    //      iend( triplet_list_.end() ); it != iend; ++it, ++count, --pos ) {
    //     Uint1 cnt( counts[*it] );
    //     add_triplet_info( score, counts, *it );
    //     if( cnt > 0 && score*10 > thresholds_[count] ) {
    //         while(    perfect_iter != P.end()
    //                && pos <= perfect_iter->bounds_.first ) {
    //             if(    max_perfect_score == 0
    //                 || max_len*perfect_iter->score_
    //                    > max_perfect_score*perfect_iter->len_ ) {
    //                 max_perfect_score = perfect_iter->score_;
    //                 max_len = perfect_iter->len_;
    //             }
    //             ++perfect_iter;
    //         }
    //         if( max_perfect_score == 0 || score*max_len >= max_perfect_score*count ) {
    //             max_perfect_score = score;
    //             max_len = count;
    //             perfect_iter = P.insert(
    //                     perfect_iter, perfect( pos, stop_ + 1,
    //                     max_perfect_score, count ) );
    //         }
    //     }
    // }
    // ```
    // NCBI's P is a std::list, so the insert is O(1); a VecDeque insert moves
    // every later element. Insert positions never move backwards during one
    // call, and stepping over an element this call inserted leaves
    // max_perfect_score/max_len unchanged (it compares the element with
    // itself), so the result is the old list with the new elements merged in
    // at the positions where the walk stood.
    fn x_find_perfect_merge(&mut self, perfect_list: &mut VecDeque<PerfectInterval>) {
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
            if pos > 0 {
                pos -= 1;
            }
        }
        if inserted {
            merged.extend(perfect_list.iter().skip(read).copied());
            perfect_list.clear();
            perfect_list.extend(merged);
        }
    }

    fn find_perfect_reference(&mut self, perfect_list: &mut VecDeque<PerfectInterval>) {
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

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:214-216,229-233
    /// ```c
    /// CSymDustMasker::operator()( const sequence_type & seq,
    ///                             size_type start, size_type stop )
    /// ...
    /// while( stop > 2 + start )    // there must be at least one triplet
    /// {
    ///     // initializations
    ///     P.clear();
    ///     triplets w( window_, low_k_, P, thresholds_ );
    /// ```
    /// The test compares the ring-buffer masker with the port of this function on random
    /// sequences, levels, windows and ranges.
    /// The ring-buffer window must give the reference intervals on random
    /// sequences of every kind the masker distinguishes: plain, biased,
    /// tandem repeats of every short period, homopolymer runs and Ns.
    #[test]
    fn x_ring_buffer_dust_matches_reference_on_random_sequences() {
        struct Rng(u64);
        impl Rng {
            fn next(&mut self) -> u64 {
                self.0 ^= self.0 << 13;
                self.0 ^= self.0 >> 7;
                self.0 ^= self.0 << 17;
                self.0
            }
            fn below(&mut self, n: u64) -> u64 {
                self.next() % n
            }
        }
        let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(3000);
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let bases = b"ACGT";
        let mut masked_cases = 0usize;
        for case in 0..cases {
            let len = 1 + rng.below(1500) as usize;
            let mut seq = Vec::with_capacity(len);
            while seq.len() < len {
                let piece = 1 + rng.below(200) as usize;
                match rng.below(7) {
                    0 | 1 => {
                        for _ in 0..piece {
                            seq.push(bases[rng.below(4) as usize]);
                        }
                    }
                    2 => {
                        // biased composition
                        let major = bases[rng.below(4) as usize];
                        for _ in 0..piece {
                            seq.push(if rng.below(5) == 0 {
                                bases[rng.below(4) as usize]
                            } else {
                                major
                            });
                        }
                    }
                    3 => {
                        // tandem repeat with a few substitutions
                        let period = 1 + rng.below(7) as usize;
                        let unit: Vec<u8> =
                            (0..period).map(|_| bases[rng.below(4) as usize]).collect();
                        for i in 0..piece {
                            seq.push(if rng.below(25) == 0 {
                                bases[rng.below(4) as usize]
                            } else {
                                unit[i % period]
                            });
                        }
                    }
                    4 => {
                        let base = bases[rng.below(4) as usize];
                        seq.extend(std::iter::repeat(base).take(piece));
                    }
                    5 => {
                        let n = 1 + rng.below(12) as usize;
                        seq.extend(std::iter::repeat(b'N').take(n));
                    }
                    _ => {
                        // other IUPAC letters and lower case
                        for _ in 0..piece.min(20) {
                            seq.push(b"acgtRYKMSWnBDHV"[rng.below(15) as usize]);
                        }
                    }
                }
            }
            seq.truncate(len);
            let (level, window, linker) = match case % 4 {
                0 => (20, 64, 1),
                1 => (
                    10 + rng.below(40) as u32,
                    8 + rng.below(57) as usize,
                    1 + rng.below(32) as usize,
                ),
                2 => (20, 64, 1),
                _ => (2 + rng.below(63) as u32, 64, 1),
            };
            let masker = DustMasker::new(level, window, linker);
            let start = if case % 5 == 0 {
                rng.below(len as u64) as usize
            } else {
                0
            };
            let stop = if case % 7 == 0 {
                rng.below(len as u64 + 1) as usize
            } else {
                len
            };
            let reference = masker.mask_subsequence_reference(&seq, start, stop);
            let fast = masker.x_mask_subsequence(&seq, start, stop);
            assert_eq!(
                fast, reference,
                "case {case}: level={level} window={window} linker={linker} start={start} stop={stop}"
            );
            if !reference.is_empty() {
                masked_cases += 1;
            }
        }
        // the generator must actually exercise the masking paths
        assert!(
            masked_cases * 2 > cases,
            "only {masked_cases} of {cases} cases masked anything"
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
