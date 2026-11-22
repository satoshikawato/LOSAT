//! SEG filter implementation for masking low-complexity regions in amino acid sequences.
//!
//! This implements the SEG algorithm as described in NCBI BLAST.
//! The algorithm identifies low-complexity regions by calculating the K2 complexity
//! score within sliding windows, using the probability-based method from
//! Wootton & Federhen.
//!
//! Reference:
//! - NCBI BLAST blast_seg.c/h
//! - Wootton & Federhen (1993) "Statistics of local complexity in amino acid sequences and sequence databases"
//! - Wootton & Federhen (1996) Methods Enzymol. 266:554-71

use crate::utils::dust::MaskedInterval;
use crate::utils::glibc_fma_libm::GlibcLibm;

use super::seg_lnfact::LNFAC;

/// Natural log of 20 (alphabet size for amino acids).
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2196
const K_LN20: f64 = 2.9957322735539909;

/// NCBI BLAST ln(2) constant.
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/ncbi_math.h:161
const NCBIMATH_LN2: f64 = std::f64::consts::LN_2; // the same f64 as NCBI's literal

/// Natural log values for 0.1, 0.2, 0.3, ... 1.0 used by NCBI when window total=10.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1305-1311
#[expect(
    clippy::approx_constant,
    reason = "NCBI's log_win10[] digits (blast_seg.c:1307-1311): -2.30258509 is not the f64 of -LN_10, and the std constant would change SEG's entropy"
)]
const LOG_WIN10: [f64; 11] = [
    0.0,
    -2.30258509,
    -1.60943791,
    -1.203982804,
    -0.91629073,
    -0.6931478,
    -0.510825623,
    -0.356674944,
    -0.22314355,
    -0.105360515,
    0.0,
];

/// NCBI s_lnfact: log(n!) using either tabulated data (lnfact[]) or Stirling's formula.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1854-1859
#[inline]
fn s_lnfact(n: usize) -> f64 {
    if n < LNFAC.len() {
        LNFAC[n]
    } else {
        let nf = n as f64;
        (nf + 0.5) * nf.glibc_ln() - nf + 0.9189385332
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1324-1350,2190-2217
// `lnalphasize` is the constant `K_LN20`, which s_GetProb (below) reads directly.
#[derive(Debug, Clone)]
struct SegAlpha {
    alphasize: usize,
    alphaindex: [i32; 128],
    alphaflag: [bool; 128],
}

impl SegAlpha {
    fn aa20alpha_std() -> Self {
        let mut alphaindex = [20i32; 128];
        let mut alphaflag = [true; 128];
        let mut next_index = 0i32;

        for letter in 0u8..128 {
            if letter == 1 || (3..=20).contains(&letter) || letter == 22 {
                alphaflag[letter as usize] = false;
                alphaindex[letter as usize] = next_index;
                next_index += 1;
            }
        }

        Self {
            alphasize: 20,
            alphaindex,
            alphaflag,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1335-1350,1541-1705
#[derive(Debug)]
struct SegWindow<'a> {
    storage: &'a [u8],
    alpha: &'a SegAlpha,
    start: usize,
    length: usize,
    bogus: i32,
    punctuation: bool,
    composition: Vec<i32>,
    state: Vec<i32>,
    entropy: f64,
}

impl<'a> SegWindow<'a> {
    fn open(storage: &'a [u8], start: usize, length: usize, alpha: &'a SegAlpha) -> Option<Self> {
        if start > storage.len() || length > storage.len().saturating_sub(start) {
            return None;
        }

        let mut win = Self {
            storage,
            alpha,
            start,
            length,
            bogus: 0,
            punctuation: false,
            composition: Vec::new(),
            state: Vec::new(),
            entropy: -2.0,
        };
        win.state_on();
        Some(win)
    }

    #[inline]
    fn slice(&self) -> &'a [u8] {
        &self.storage[self.start..self.start + self.length]
    }

    fn comp_on(&mut self) {
        self.composition = vec![0; self.alpha.alphasize];
        self.bogus = 0;

        let slice = &self.storage[self.start..self.start + self.length];
        for &letter in slice {
            let letter = letter as usize;
            if letter < self.alpha.alphaflag.len() && !self.alpha.alphaflag[letter] {
                let index = self.alpha.alphaindex[letter] as usize;
                self.composition[index] += 1;
            } else {
                self.bogus += 1;
            }
        }
    }

    fn state_on(&mut self) {
        if self.composition.is_empty() {
            self.comp_on();
        }

        self.state = vec![0; self.alpha.alphasize + 1];
        let mut nel = 0usize;
        for &count in &self.composition {
            if count != 0 {
                self.state[nel] = count;
                nel += 1;
            }
        }
        self.state[..nel].sort_unstable_by(|left, right| right.cmp(left));
    }

    fn entropy_on(&mut self) {
        if self.state.is_empty() {
            self.state_on();
        }
        self.entropy = entropy_from_state_vector(&self.state);
    }

    fn has_dash(&self) -> bool {
        self.slice().contains(&b'-')
    }

    fn shift_win1(&mut self) -> bool {
        if self.start + self.length >= self.storage.len() {
            return false;
        }

        let outgoing = self.storage[self.start] as usize;
        if outgoing < self.alpha.alphaflag.len() && !self.alpha.alphaflag[outgoing] {
            let index = self.alpha.alphaindex[outgoing] as usize;
            let class = self.composition[index];
            decrement_sv(&mut self.state, class);
            self.composition[index] -= 1;
        } else {
            self.bogus -= 1;
        }

        let incoming = self.storage[self.start + self.length] as usize;
        self.start += 1;

        if incoming < self.alpha.alphaflag.len() && !self.alpha.alphaflag[incoming] {
            let index = self.alpha.alphaindex[incoming] as usize;
            let class = self.composition[index];
            increment_sv(&mut self.state, class);
            self.composition[index] += 1;
        } else {
            self.bogus += 1;
        }

        if self.entropy > -2.0 {
            self.entropy = entropy_from_state_vector(&self.state);
        }

        true
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1631-1659
fn decrement_sv(state: &mut [i32], class: i32) {
    let mut index = 0usize;
    while index < state.len() && state[index] != 0 {
        let current = state[index];
        let next = state.get(index + 1).copied().unwrap_or(i32::MIN);
        if current == class && next < class {
            state[index] = current - 1;
            break;
        }
        index += 1;
    }
}

fn increment_sv(state: &mut [i32], class: i32) {
    let mut index = 0usize;
    while index < state.len() {
        if state[index] == class {
            state[index] += 1;
            break;
        }
        index += 1;
    }
}

/// Calculate ln(number of permutations) for a given state vector
/// This is equation 3 from Wootton & Federhen: ln(n!) - sum(ln(count_i!))
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c s_LnPerm()
fn ln_perm(sv: &[i32], window_length: i32) -> f64 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1867-1882
    let mut ans = s_lnfact(window_length as usize);
    for &count in sv {
        if count == 0 {
            break;
        }
        ans -= s_lnfact(count as usize);
    }
    ans
}

/// Calculate ln(number of compositions) for a given state vector
/// This is equation 1 from Wootton & Federhen
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c s_LnAss()
fn ln_ass(sv: &[i32], alphasize: i32) -> f64 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1892-1933
    let a = alphasize as usize;
    let mut ans = s_lnfact(a);

    // sv is expected to be length `alphasize` with trailing zeros; our `sv` is
    // typically the non-zero prefix. Treat missing entries as zeros.
    let sv0 = if !sv.is_empty() { sv[0] } else { 0 };
    if sv0 == 0 {
        return ans;
    }

    let mut total: i32 = alphasize;
    let mut class: i32 = 1;
    let mut svim1: i32 = sv0;

    // NCBI-style loop: increment `i` first, then stop when i == alphasize,
    // otherwise consume next sv entry and compare.
    let mut i: usize = 0;
    loop {
        i += 1;
        if i == a {
            ans -= s_lnfact(class as usize);
            break;
        }

        let svi: i32 = if i < sv.len() { sv[i] } else { 0 };
        if svi == svim1 {
            class += 1;
        } else {
            total -= class;
            ans -= s_lnfact(class as usize);
            if svi == 0 {
                ans -= s_lnfact(total as usize);
                break;
            }
            class = 1;
        }
        svim1 = svi;
    }

    ans
}

/// Calculate the K2 complexity score (probability-based)
/// This is the natural log of P_0 from equation 3 of Wootton & Federhen
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c s_GetProb()
fn get_prob(sv: &[i32], window_length: i32) -> f64 {
    let alphasize = 20; // Standard amino acid alphabet
    let totseq = (window_length as f64) * K_LN20;

    let ans1 = ln_ass(sv, alphasize);
    // NCBI guards `ans2` computation on `ans1 > -100000` and `sv[0]` not equal to `INT4_MIN`.
    let ans2 = if ans1 > -100000.0 {
        ln_perm(sv, window_length)
    } else {
        0.0
    };

    ans1 + ans2 - totseq
}

/// NCBI SEG entropy (Shannon entropy, bits) implementation.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1592-1624 (s_Entropy)
fn entropy_from_state_vector(state: &[i32]) -> f64 {
    // Memo (was LOSAT_X_SEGMEMO): the entropy is a pure function of the
    // state vector, and a 12-letter window has at most a few hundred distinct
    // state vectors, so remember the value computed for each one.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1598-1603,1616-1621
    // `entropy_from_state_vector_reference` is the port of this function. The memo is keyed by the sorted
    // counts of the state vector, which is all that `s_Entropy` reads, so a repeated vector gets the value
    // the first call returned. A vector that cannot be packed (a window of more than 15 letters) takes the
    // port directly.
    let mut key = 0u64;
    let mut classes = 0u32;
    let mut total = 0i32;
    let mut packable = true;
    for &count in state {
        if count == 0 {
            break;
        }
        total += count;
        // Windows of at most 15 letters: fewer than 700 possible keys.
        if count < 0 || total > 15 {
            packable = false;
            break;
        }
        key |= (count as u64) << (4 * classes);
        classes += 1;
    }
    if packable {
        if let Some(value) = x_memo_get(x_state_entropies(), key) {
            return value;
        }
    }
    entropy_from_state_vector_reference(state)
}

// No NCBI counterpart: the windows of at most 15 letters, as the count vectors NCBI's state
// vector holds for them: every non-increasing list of positive counts summing to at most 15 (684
// lists, the empty one among them), each passed to `visit`.
fn x_each_count_list(visit: &mut impl FnMut(&[i32])) {
    fn more(list: &mut Vec<i32>, left: i32, most: i32, visit: &mut impl FnMut(&[i32])) {
        visit(list);
        for count in (1..=most.min(left)).rev() {
            list.push(count);
            more(list, left - count, count, visit);
            list.pop();
        }
    }
    more(&mut Vec::new(), 15, 15, visit);
}

// No NCBI counterpart: an open-addressing table of `slots` (key, entropy bits) pairs, one for each
// count list of `x_each_count_list`, keyed by `key` and valued by `value`.
fn x_entropy_table(
    slots: usize,
    key: impl Fn(&[i32]) -> u64,
    value: impl Fn(&[i32], u64) -> f64,
) -> Box<[(u64, u64)]> {
    let mut table = vec![(u64::MAX, 0); slots].into_boxed_slice();
    let mask = slots - 1;
    x_each_count_list(&mut |counts| {
        let key = key(counts);
        let mut slot = (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 40) as usize & mask;
        while table[slot].0 != u64::MAX {
            slot = (slot + 1) & mask;
        }
        table[slot] = (key, value(counts, key).to_bits());
    });
    table
}

// No NCBI counterpart: the entropy stored for `key` in a table of `x_entropy_table`, if any.
#[inline]
fn x_memo_get(table: &[(u64, u64)], key: u64) -> Option<f64> {
    let mask = table.len() - 1;
    let mut slot = (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 40) as usize & mask;
    loop {
        let (stored, bits) = table[slot];
        if stored == key {
            return Some(f64::from_bits(bits));
        }
        if stored == u64::MAX {
            return None;
        }
        slot = (slot + 1) & mask;
    }
}

// No NCBI counterpart: memo table of entropy values per state vector (`s_Entropy` is a pure
// function of the vector); it does not change any value NCBI computes. Fewer than 700 state
// vectors sum to at most 15, so the table holds each of them, as NCBI keeps the vector (counts in
// non-increasing order), computed once for the process and read by every thread. A vector not in
// it is computed by the port.
fn x_state_entropies() -> &'static [(u64, u64)] {
    static TABLE: std::sync::OnceLock<Box<[(u64, u64)]>> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        x_entropy_table(
            4096,
            |counts| {
                counts.iter().enumerate().fold(0, |key, (class, &count)| {
                    key | ((count as u64) << (4 * class))
                })
            },
            |counts, _| entropy_from_state_vector_reference(counts),
        )
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1598-1606,1616-1621
// The original `entropy_from_state_vector`, renamed; it is the port of this function.
fn entropy_from_state_vector_reference(state: &[i32]) -> f64 {
    let mut total = 0i32;
    for &count in state {
        if count == 0 {
            break;
        }
        total += count;
    }
    if total == 0 {
        return 0.0;
    }

    let mut ent = 0.0f64;
    if total == 10 {
        for &count in state {
            if count == 0 {
                break;
            }
            ent += (count as f64) * LOG_WIN10[count as usize] / NCBIMATH_LN2;
        }
    } else {
        let total_f = total as f64;
        for &count in state {
            if count == 0 {
                break;
            }
            ent += (count as f64) * ((count as f64) / total_f).glibc_ln() / NCBIMATH_LN2;
        }
    }

    (ent / (total as f64)).abs()
}

/// SEG filter parameters
#[derive(Debug, Clone)]
pub struct SegParams {
    /// Window size (default: 12)
    pub window: usize,
    /// Low complexity threshold (default: 2.2)
    pub locut: f64,
    /// High complexity threshold (default: 2.5)
    pub hicut: f64,
    /// Maximum number of invalid amino acids allowed in window (default: 2)
    pub maxbogus: usize,
    /// Maximum trim size for boundary optimization (default: 50)
    pub maxtrim: usize,
}

impl Default for SegParams {
    fn default() -> Self {
        Self {
            window: 12,
            locut: 2.2,
            hicut: 2.5,
            maxbogus: 2,
            maxtrim: 50,
        }
    }
}

impl SegParams {
    pub fn new(window: usize, locut: f64, hicut: f64) -> Self {
        // NCBI parameter normalization (c++/src/algo/blast/core/blast_seg.c:2251-2258:s_SegParametersCheck):
        // (`< 0.0`, so a NaN stays: segmasker's `-locut -nan`.)
        let window = if window > 0 { window } else { 12 };
        let locut = if locut < 0.0 { 0.0 } else { locut };
        let mut hicut = if hicut < 0.0 { 0.0 } else { hicut };
        if locut > hicut {
            hicut = locut;
        }
        Self {
            window,
            locut,
            hicut,
            maxbogus: 2,
            maxtrim: 50,
        }
    }

    /// Create with all parameters including maxbogus and maxtrim
    pub fn with_all(
        window: usize,
        locut: f64,
        hicut: f64,
        maxbogus: usize,
        maxtrim: usize,
    ) -> Self {
        let mut params = Self::new(window, locut, hicut);
        params.maxbogus = maxbogus.min(window);
        params.maxtrim = maxtrim;
        params
    }
}

/// SEG masker implementation following NCBI BLAST's SEG algorithm
pub struct SegMasker {
    alpha: SegAlpha,
    window: usize,
    locut: f64,
    hicut: f64,
    maxbogus: usize,
    maxtrim: usize,
    downset: usize,
    upset: usize,
    /// LOSAT_X_SEGTRIMMEMO (graduated 2026-10-26): whether the fast path takes a segment's trim from the
    /// memo shared by all threads (`x_seg_trim_memo.rs`); set by the maskers that repeat segments.
    x_trim_memo: bool,
}

impl SegMasker {
    /// Create a new SEG masker with the given parameters
    pub fn new(window: usize, locut: f64, hicut: f64) -> Self {
        let params = SegParams::new(window, locut, hicut);

        // downset = (window+1)/2 - 1, upset = window - downset
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2050-2051
        let downset = params.window.div_ceil(2) - 1;
        let upset = params.window - downset;

        Self {
            alpha: SegAlpha::aa20alpha_std(),
            window: params.window,
            locut: params.locut,
            hicut: params.hicut,
            maxbogus: params.maxbogus,
            maxtrim: params.maxtrim,
            downset,
            upset,
            x_trim_memo: false,
        }
    }

    /// LOSAT_X_SEGTRIMMEMO (graduated 2026-10-26): the masker that takes the trim of a segment it (or
    /// another such masker, on any thread) trimmed before from the memo. For maskers
    /// that see the same segments again (TBLASTN's redo windows, BLASTX's query and its chunks);
    /// the result is the same either way.
    // No NCBI counterpart: chooses where a trim result comes from; it does not change any value
    // NCBI computes.
    pub(crate) fn x_with_trim_memo(mut self) -> Self {
        self.x_trim_memo = true;
        self
    }

    /// Create a SEG masker with default parameters (window=12, locut=2.2, hicut=2.5)
    #[cfg(test)]
    pub fn with_defaults() -> Self {
        Self::new(12, 2.2, 2.5)
    }

    /// Create with full parameters including maxbogus and maxtrim
    pub fn with_params(params: &SegParams) -> Self {
        let downset = params.window.div_ceil(2) - 1;
        let upset = params.window - downset;

        Self {
            alpha: SegAlpha::aa20alpha_std(),
            window: params.window,
            locut: params.locut,
            hicut: params.hicut,
            maxbogus: params.maxbogus,
            maxtrim: params.maxtrim,
            downset,
            upset,
            x_trim_memo: false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1541-1604,1745-1795
    #[cfg(test)]
    fn calculate_entropy(&self, window: &[u8]) -> (f64, usize) {
        let Some(mut win) = SegWindow::open(window, 0, window.len(), &self.alpha) else {
            return (0.0, 0);
        };
        win.entropy_on();
        (win.entropy, win.bogus as usize)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1745-1803
    fn calculate_entropy_array(&self, seq: &[u8]) -> Vec<f64> {
        let len = seq.len();
        let mut h = vec![-1.0; len];

        if len < self.window {
            return h;
        }

        let Some(mut win) = SegWindow::open(seq, 0, self.window, &self.alpha) else {
            return h;
        };
        win.entropy_on();

        let first = self.downset;
        let last = len - self.upset;

        for i in first..=last {
            if win.punctuation && win.has_dash() {
                h[i] = -1.0;
                win.shift_win1();
                continue;
            }
            if win.bogus > self.maxbogus as i32 {
                h[i] = -1.0;
                win.shift_win1();
                continue;
            }
            h[i] = win.entropy;
            win.shift_win1();
        }

        h
    }

    /// Find the left boundary of a low-complexity region
    /// Starting from position i, search left until entropy > hicut
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1813-1825 (s_FindLow)
    fn find_low(&self, i: usize, limit: usize, h: &[f64]) -> usize {
        let mut j: isize = i as isize;
        let limit: isize = limit as isize;
        while j >= limit {
            let v = h[j as usize];
            if v == -1.0 || v > self.hicut {
                break;
            }
            j -= 1;
        }
        (j + 1) as usize
    }

    /// Find the right boundary of a low-complexity region
    /// Starting from position i, search right until entropy > hicut
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1836-1848 (s_FindHigh)
    fn find_high(&self, i: usize, limit: usize, h: &[f64]) -> usize {
        let mut j: isize = i as isize;
        let limit: isize = limit as isize;
        while j <= limit {
            let v = h[j as usize];
            if v == -1.0 || v > self.hicut {
                break;
            }
            j += 1;
        }
        (j - 1) as usize
    }

    /// NCBI s_Trim: trim [leftend..=rightend] to minimize s_GetProb.
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1974-2018
    ///
    /// Note: `leftend` and `rightend` are **inclusive** indices in `seq`.
    fn trim_segment(&self, seq: &[u8], leftend: usize, rightend: usize) -> (usize, usize) {
        if leftend >= seq.len() || rightend >= seq.len() || leftend > rightend {
            return (leftend, rightend);
        }
        let seg_len = rightend - leftend + 1;
        if seg_len == 0 {
            return (leftend, rightend);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1984-2007
        let mut minlen: usize = 1;
        if seg_len.saturating_sub(self.maxtrim) > minlen {
            minlen = seg_len - self.maxtrim;
        }

        let segment = &seq[leftend..=rightend];
        let mut best_lend: usize = 0;
        let mut best_rend: usize = seg_len - 1;
        let mut minprob: f64 = 1.0;

        for cur_len in (minlen + 1..=seg_len).rev() {
            let Some(mut win) = SegWindow::open(segment, 0, cur_len, &self.alpha) else {
                continue;
            };

            let mut shift = true;
            let mut i = 0usize;
            while shift {
                let prob = get_prob(&win.state, cur_len as i32);
                if prob < minprob {
                    minprob = prob;
                    best_lend = i;
                    best_rend = i + cur_len - 1;
                }
                shift = win.shift_win1();
                i += 1;
            }
        }

        (leftend + best_lend, leftend + best_rend)
    }

    /// Mask a sequence and return the list of masked intervals, the overlapping segments merged:
    /// `SeqBufferSeg` as BLAST's filters call it (`sparamsp->overlaps` set to `TRUE`, blast_filter.c
    /// and the composition adjustment's `s_DoSegSequenceData`).
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2030-2116 (s_SegSeq)
    pub fn mask_sequence(&self, seq: &[u8]) -> Vec<MaskedInterval> {
        self.seq_buffer_seg::<true>(seq)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2224-2238,2311-2316,2320-2322
    /// `SeqBufferSeg` with `sparamsp->overlaps` = `OVERLAPS`: the intervals of the segments in
    /// NCBI's order (`s_SegsToBlastSeqLoc`), merged where they overlap when `OVERLAPS`. The
    /// masker of `losat segmasker` (`CSegMasker`) keeps `SegParametersNewAa`'s `FALSE`; BLAST's
    /// filters set `TRUE` (`mask_sequence`). Each value is its own copy of the code, so the
    /// filters' copy is the same machine code as before the masker of segmasker came.
    #[inline(always)]
    pub(crate) fn seq_buffer_seg<const OVERLAPS: bool>(&self, seq: &[u8]) -> Vec<MaskedInterval> {
        if seq.len() < self.window {
            return Vec::new();
        }
        // Fast path (was LOSAT_X_SEGFAST): see `x_mask_sequence`.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2050-2053,2062-2064
        // The fast path runs where `x_fast_applies` (window of at most 15 letters, 20-letter alphabet);
        // other settings run `mask_sequence_reference`, the port of `s_SegSeq`. The fast path keeps the
        // control flow of this loop and replaces only the entropy and trim computations.
        if self.x_fast_applies() {
            return self.x_mask_sequence::<OVERLAPS>(seq);
        }
        self.mask_sequence_reference::<OVERLAPS>(seq)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2030-2032,2053-2060
    // The original body of `mask_sequence`, moved here unchanged; it is the port of this function. It runs
    // where the fast path does not apply (`x_fast_applies`).
    fn mask_sequence_reference<const OVERLAPS: bool>(&self, seq: &[u8]) -> Vec<MaskedInterval> {
        let mut segs_inclusive: Vec<(usize, usize)> = Vec::new();

        fn seg_seq(masker: &SegMasker, seq: &[u8], offset: usize, segs: &mut Vec<(usize, usize)>) {
            if seq.len() < masker.window {
                return;
            }

            let h = masker.calculate_entropy_array(seq);
            let first = masker.downset;
            let last = seq.len().saturating_sub(masker.upset);
            let mut lowlim = first;

            let mut i = first;
            while i <= last {
                if h[i] != -1.0 && h[i] <= masker.locut {
                    let loi = masker.find_low(i, lowlim, &h);
                    let hii = masker.find_high(i, last, &h);

                    // NCBI:
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2070-2071
                    let mut leftend = loi - masker.downset;
                    let mut rightend = hii + masker.upset - 1; // inclusive

                    // Trim to minimize probability (NCBI s_Trim)
                    (leftend, rightend) = masker.trim_segment(seq, leftend, rightend);

                    // NCBI recursion for trigger window in left trim:
                    // when `i + upset - 1 < leftend`, recurse on `[lend..rend]`.
                    if i + masker.upset - 1 < leftend {
                        let lend = loi - masker.downset;
                        if lend < leftend {
                            let rend = leftend - 1;
                            if rend < seq.len() && lend <= rend {
                                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2086-2101
                                // The assignment of `*segs` to `leftsegs->next` overwrites the link of the head of
                                // the recursion's list, so when the recursion found more than
                                // one segment only its head (the one found last) is kept; the
                                // others are dropped (and leaked). This deterministic result
                                // is reproduced.
                                let mut left_segs: Vec<(usize, usize)> = Vec::new();
                                seg_seq(masker, &seq[lend..=rend], offset + lend, &mut left_segs);
                                if let Some(&head) = left_segs.first() {
                                    segs.insert(0, head);
                                }
                            }
                        }
                    }

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2104-2109
                    segs.insert(0, (leftend + offset, rightend + offset));

                    // NCBI:
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2110-2111
                    i = hii.min(rightend + masker.downset);
                    lowlim = i + 1;
                    i += 1;
                } else {
                    i += 1;
                }
            }
        }

        seg_seq(self, seq, 0, &mut segs_inclusive);

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1142-1157
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2321-2322
        if OVERLAPS {
            merge_segs_inclusive(seq.len(), &mut segs_inclusive);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2325
        segs_inclusive.reverse();

        segs_inclusive
            .into_iter()
            .map(|(b, e)| MaskedInterval::new(b, e.saturating_add(1)))
            .collect()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2125-2149
// ---------------------------------------------------------------------------
// Fast path (was LOSAT_X_SEGFAST): the same SEG, with the window kept as counts.
//
// NCBI keeps a window as a "state vector": the counts of its letters sorted in
// decreasing order (blast_seg.c s_StateOn / s_DecrementSV / s_IncrementSV), and
// derives the entropy (s_Entropy), ln(compositions) (s_LnAss) and
// ln(permutations) (s_LnPerm) from it.  All three are sums taken in the order of
// that vector, so they are functions of *how many letters occur k times* for each
// k.  The code below maintains exactly that (a histogram of counts) under a
// one-letter shift in O(1) and then performs NCBI's additions and subtractions in
// NCBI's order: the largest count first, each count once per letter.
//
// Everything else (find_low / find_high, the recursion, the merge) is shared with
// the reference implementation.
// ---------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1529-1533,1537-1537
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1642-1646,1658-1661
// The state vector is the sorted list of letter counts. A shift by one letter lowers the count of the
// leaving letter and raises the count of the entering one (`s_DecrementSV`, `s_IncrementSV`). The code
// below keeps a histogram of the counts instead, which holds the same information, and updates it per
// shift.

// No NCBI counterpart: size bound of the histogram arrays of `x_trim_segment`; longer segments use the reference trim; it does not change any value NCBI computes.
/// Longest segment handled by the count histogram of `x_trim_segment`.
pub(crate) const X_TRIM_MAX: usize = 127;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1593-1593,1599-1602,1621-1623
// No NCBI counterpart: memo table of `s_Entropy` values per count histogram. The entropy depends only on the
// counts, so reusing it changes no value.
// Entropy by count histogram (4 bits per count value, windows up to 15): the table holds every
// histogram of a window of at most 15 letters (fewer than 700), computed once for the process
// and read by every thread. A histogram not in it is computed by the port.
fn x_histogram_entropies() -> &'static [(u64, u64)] {
    static TABLE: std::sync::OnceLock<Box<[(u64, u64)]>> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        x_entropy_table(
            2048,
            |counts| {
                counts
                    .iter()
                    .fold(0, |key, &count| key + (1u64 << (4 * (count - 1))))
            },
            |_, key| x_entropy_of_histogram(key),
        )
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1598-1603,1616-1621
// The sorted state vector is rebuilt from the histogram (largest count first), and the port of `s_Entropy`
// is called on it.
/// NCBI's entropy of the window whose count histogram is `key` (nibble `k - 1`
/// holds the number of letters that occur `k` times).
fn x_entropy_of_histogram(key: u64) -> f64 {
    let mut state = [0i32; 16];
    let mut used = 0usize;
    for count in (1..=15usize).rev() {
        let letters = (key >> (4 * (count - 1))) & 0xF;
        for _ in 0..letters {
            state[used] = count as i32;
            used += 1;
        }
    }
    entropy_from_state_vector_reference(&state[..used])
}

// No NCBI counterpart: a bit set of the counts present in a window, the index into the histogram of `x_prob_of_histogram`; it does not change any value NCBI computes.
/// The counts that occur in a window: bit `k` is set when some letter occurs
/// `k` times (`k <= X_TRIM_MAX = 127`). Two 64-bit words rather than a `u128`,
/// whose shifts are library calls on wasm32.
#[derive(Clone, Copy, PartialEq, Eq)]
struct XCountMask([u64; 2]);

impl XCountMask {
    const EMPTY: Self = Self([0, 0]);

    #[inline(always)]
    fn set(&mut self, count: usize) {
        self.0[count >> 6] |= 1u64 << (count & 63);
    }

    #[inline(always)]
    fn clear(&mut self, count: usize) {
        self.0[count >> 6] &= !(1u64 << (count & 63));
    }

    #[inline(always)]
    fn is_empty(&self) -> bool {
        (self.0[0] | self.0[1]) == 0
    }

    /// The largest count in the mask (which must not be empty).
    #[inline(always)]
    fn largest(&self) -> usize {
        if self.0[1] != 0 {
            127 - self.0[1].leading_zeros() as usize
        } else {
            63 - self.0[0].leading_zeros() as usize
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1949-1955,1962-1962
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1901-1905,1909-1912,1919-1923
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1874-1879
// The terms of `s_LnAss` and `s_LnPerm` are taken in the order of the sorted state vector (largest count first),
// as here: one `s_lnfact(class)` per run of equal counts, and one `s_lnfact(count)` per letter.
/// `get_prob` from the histogram: `mask` has bit `k` set when some letter
/// occurs `k` times and `letters_with[k]` says how many.
#[inline]
fn x_prob_of_histogram(
    mask: XCountMask,
    letters_with: &[u8; X_TRIM_MAX + 1],
    window_length: usize,
) -> f64 {
    let totseq = (window_length as f64) * K_LN20;
    // s_LnAss: one term per run of equal counts, largest count first, then the
    // letters that do not occur (absent when all twenty do).
    let mut ans1 = s_lnfact(20);
    if !mask.is_empty() {
        let mut total = 20usize;
        let mut rest = mask;
        while !rest.is_empty() {
            let count = rest.largest();
            rest.clear(count);
            let class = letters_with[count] as usize;
            total -= class;
            ans1 -= s_lnfact(class);
        }
        if total > 0 {
            ans1 -= s_lnfact(total);
        }
    }
    // s_LnPerm: one term per letter, largest count first.
    let ans2 = if ans1 > -100000.0 {
        let mut ans = s_lnfact(window_length);
        let mut rest = mask;
        while !rest.is_empty() {
            let count = rest.largest();
            rest.clear(count);
            let term = s_lnfact(count);
            for _ in 0..letters_with[count] {
                ans -= term;
            }
        }
        ans
    } else {
        0.0
    };
    ans1 + ans2 - totseq
}

impl SegMasker {
    // No NCBI counterpart: the fast path needs windows of at most 15 letters and the 20-letter alphabet; other settings use the reference; it does not change any value NCBI computes.
    fn x_fast_applies(&self) -> bool {
        (1..=15).contains(&self.window) && self.alpha.alphasize == 20
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1499-1505
    // The table maps a residue to `alphaindex` when `alphaflag` is clear, and to 255 (counted as bogus) otherwise.
    /// 0..19 for a letter of the SEG alphabet, 255 for anything else.
    fn x_letter_index(&self) -> [u8; 256] {
        let mut table = [255u8; 256];
        for (letter, slot) in table
            .iter_mut()
            .enumerate()
            .take(self.alpha.alphaflag.len())
        {
            if !self.alpha.alphaflag[letter] {
                *slot = self.alpha.alphaindex[letter] as u8;
            }
        }
        table
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1775-1776,1778-1782,1789-1796
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1695-1704
    // The first window is counted once; each later window follows from the previous one by one removal and one
    // addition, as in `s_ShiftWin1`. The entropy of each window with at most `maxbogus` bogus letters comes from
    // the histogram.
    /// `calculate_entropy_array` with the window kept as a count histogram.
    fn x_entropy_array(&self, seq: &[u8], index: &[u8; 256]) -> Vec<f64> {
        let len = seq.len();
        let mut h = vec![-1.0; len];
        let window = self.window;
        if len < window {
            return h;
        }
        let mut counts = [0u8; 20];
        let mut bogus = 0i32;
        let mut key = 0u64;
        for &letter in &seq[..window] {
            let letter = index[letter as usize];
            if letter == 255 {
                bogus += 1;
            } else {
                let count = counts[letter as usize];
                if count > 0 {
                    key -= 1u64 << (4 * (count - 1));
                }
                key += 1u64 << (4 * count);
                counts[letter as usize] = count + 1;
            }
        }
        let maxbogus = self.maxbogus as i32;
        let first = self.downset;
        let windows = len - window + 1;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1621-1623
        // The entropy of a window depends only on its count histogram, so the value for a seen histogram is reused.
        let memo = x_histogram_entropies();
        for start in 0..windows {
            h[first + start] = if bogus > maxbogus {
                -1.0
            } else {
                x_memo_get(memo, key).unwrap_or_else(|| x_entropy_of_histogram(key))
            };
            if start + window < len {
                let outgoing = index[seq[start] as usize];
                if outgoing == 255 {
                    bogus -= 1;
                } else {
                    let count = counts[outgoing as usize];
                    key -= 1u64 << (4 * (count - 1));
                    if count > 1 {
                        key += 1u64 << (4 * (count - 2));
                    }
                    counts[outgoing as usize] = count - 1;
                }
                let incoming = index[seq[start + window] as usize];
                if incoming == 255 {
                    bogus += 1;
                } else {
                    let count = counts[incoming as usize];
                    if count > 0 {
                        key -= 1u64 << (4 * (count - 1));
                    }
                    key += 1u64 << (4 * count);
                    counts[incoming as usize] = count + 1;
                }
            }
        }
        h
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1988-1992,1998-2001,2007-2008,2013-2014
    // Windows of each length are visited left to right, with the smallest `prob` kept on a strict `<` as here.
    // Each window is updated from its neighbour by one removal and one addition, and the sums in `s_LnAss` and
    // `s_LnPerm` use the histogram in NCBI's order.
    // No NCBI counterpart for the attribute: on wasm32 (LOSAT Web) the trim stays a function of its
    // own, called once per segment, so that V8 runs the later segments in its optimized code (it
    // switches a function to that code only at a later call, never inside a running one); inlined
    // into `seg_seq`, the trims of a frame ran in V8's baseline code to the end of the frame. Other
    // targets keep the compiler's choice.
    /// `trim_segment` with every window kept as a count histogram.
    #[cfg_attr(target_arch = "wasm32", inline(never))]
    fn x_trim_segment(
        &self,
        seq: &[u8],
        leftend: usize,
        rightend: usize,
        index: &[u8; 256],
    ) -> (usize, usize) {
        if leftend >= seq.len() || rightend >= seq.len() || leftend > rightend {
            return (leftend, rightend);
        }
        let seg_len = rightend - leftend + 1;
        if seg_len > X_TRIM_MAX {
            return self.trim_segment(seq, leftend, rightend);
        }
        let mut minlen: usize = 1;
        if seg_len.saturating_sub(self.maxtrim) > minlen {
            minlen = seg_len - self.maxtrim;
        }
        let segment = &seq[leftend..=rightend];

        #[derive(Clone, Copy)]
        struct Window {
            counts: [u8; 20],
            letters_with: [u8; X_TRIM_MAX + 1],
            mask: XCountMask,
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1642-1646,1658-1661
        // `add` and `remove` change the count of one letter by one, as `s_IncrementSV` and `s_DecrementSV` do for
        // the state vector; here the histogram `letters_with` and its mask are updated.
        impl Window {
            #[inline(always)]
            fn add(&mut self, letter: u8) {
                if letter != 255 {
                    let count = self.counts[letter as usize] as usize;
                    if count > 0 {
                        self.letters_with[count] -= 1;
                        if self.letters_with[count] == 0 {
                            self.mask.clear(count);
                        }
                    }
                    self.letters_with[count + 1] += 1;
                    self.mask.set(count + 1);
                    self.counts[letter as usize] = (count + 1) as u8;
                }
            }
            #[inline(always)]
            fn remove(&mut self, letter: u8) {
                if letter != 255 {
                    let count = self.counts[letter as usize] as usize;
                    self.letters_with[count] -= 1;
                    if self.letters_with[count] == 0 {
                        self.mask.clear(count);
                    }
                    if count > 1 {
                        self.letters_with[count - 1] += 1;
                        self.mask.set(count - 1);
                    }
                    self.counts[letter as usize] = (count - 1) as u8;
                }
            }
        }

        // The window of the current length at the left end of the segment; it
        // loses its last letter each time the length drops by one.
        let mut leftmost = Window {
            counts: [0; 20],
            letters_with: [0; X_TRIM_MAX + 1],
            mask: XCountMask::EMPTY,
        };
        for &letter in segment {
            leftmost.add(index[letter as usize]);
        }

        let mut best_lend: usize = 0;
        let mut best_rend: usize = seg_len - 1;
        let mut minprob: f64 = 1.0;
        for cur_len in (minlen + 1..=seg_len).rev() {
            let mut win = leftmost;
            let mut i = 0usize;
            loop {
                let prob = x_prob_of_histogram(win.mask, &win.letters_with, cur_len);
                if prob < minprob {
                    minprob = prob;
                    best_lend = i;
                    best_rend = i + cur_len - 1;
                }
                if i + cur_len >= seg_len {
                    break;
                }
                win.remove(index[segment[i] as usize]);
                win.add(index[segment[i + cur_len] as usize]);
                i += 1;
            }
            leftmost.remove(index[segment[cur_len - 1] as usize]);
        }

        (leftend + best_lend, leftend + best_rend)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2062-2071,2073-2074
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2081-2084,2105-2111
    // The recursion, the left-trim handling and the order of the segments are the same as in
    // `mask_sequence_reference`; only the entropy array and the trim use the histogram.
    /// `mask_sequence` through `x_entropy_array` and `x_trim_segment`.
    fn x_mask_sequence<const OVERLAPS: bool>(&self, seq: &[u8]) -> Vec<MaskedInterval> {
        let index = self.x_letter_index();
        let mut segs_inclusive: Vec<(usize, usize)> = Vec::new();

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2030-2116 (s_SegSeq); the control flow of `seg_seq` in
        // `mask_sequence_reference`.
        fn seg_seq(
            masker: &SegMasker,
            index: &[u8; 256],
            seq: &[u8],
            offset: usize,
            segs: &mut Vec<(usize, usize)>,
        ) {
            if seq.len() < masker.window {
                return;
            }
            let h = masker.x_entropy_array(seq, index);
            let first = masker.downset;
            let last = seq.len().saturating_sub(masker.upset);
            let mut lowlim = first;

            let mut i = first;
            while i <= last {
                if h[i] != -1.0 && h[i] <= masker.locut {
                    let loi = masker.find_low(i, lowlim, &h);
                    let hii = masker.find_high(i, last, &h);

                    let mut leftend = loi - masker.downset;
                    let mut rightend = hii + masker.upset - 1; // inclusive

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2073-2074
                    // LOSAT_X_SEGTRIMMEMO (graduated 2026-10-26, always on): the trim's offsets
                    // within the segment are a function of its letters and `maxtrim`
                    // (x_seg_trim_memo.rs); for a masker that asks for it, a segment trimmed
                    // before (by any thread) takes the stored offsets.
                    (leftend, rightend) =
                        if masker.x_trim_memo && leftend <= rightend && rightend < seq.len() {
                            let base = leftend;
                            let (lend, rend) = crate::utils::x_seg_trim_memo::x_trimmed(
                                &seq[leftend..=rightend],
                                masker.maxtrim,
                                || {
                                    let (l, r) =
                                        masker.x_trim_segment(seq, leftend, rightend, index);
                                    (l - base, r - base)
                                },
                            );
                            (base + lend, base + rend)
                        } else {
                            masker.x_trim_segment(seq, leftend, rightend, index)
                        };

                    if i + masker.upset - 1 < leftend {
                        let lend = loi - masker.downset;
                        if lend < leftend {
                            let rend = leftend - 1;
                            if rend < seq.len() && lend <= rend {
                                let mut left_segs: Vec<(usize, usize)> = Vec::new();
                                seg_seq(
                                    masker,
                                    index,
                                    &seq[lend..=rend],
                                    offset + lend,
                                    &mut left_segs,
                                );
                                if let Some(&head) = left_segs.first() {
                                    segs.insert(0, head);
                                }
                            }
                        }
                    }

                    segs.insert(0, (leftend + offset, rightend + offset));

                    i = hii.min(rightend + masker.downset);
                    lowlim = i + 1;
                    i += 1;
                } else {
                    i += 1;
                }
            }
        }

        seg_seq(self, &index, seq, 0, &mut segs_inclusive);
        if OVERLAPS {
            merge_segs_inclusive(seq.len(), &mut segs_inclusive);
        }
        segs_inclusive.reverse();
        segs_inclusive
            .into_iter()
            .map(|(b, e)| MaskedInterval::new(b, e.saturating_add(1)))
            .collect()
    }
}

fn merge_segs_inclusive(seq_len: usize, segs_inclusive: &mut Vec<(usize, usize)>) {
    if segs_inclusive.is_empty() {
        return;
    }

    let max_end = seq_len.saturating_sub(1);
    segs_inclusive[0].1 = segs_inclusive[0].1.min(max_end);

    let mut index = 0usize;
    while index + 1 < segs_inclusive.len() {
        let (seg_begin, seg_end) = segs_inclusive[index];
        let (next_begin, next_end) = segs_inclusive[index + 1];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2130,2141-2145
        // Strict signed comparison merges overlaps; adjacent intervals stay
        // separate, including their distinct returned DNA-mask endpoints.
        if seg_begin <= next_end {
            segs_inclusive[index].0 = seg_begin.min(next_begin);
            segs_inclusive[index].1 = seg_end.max(next_end);
            segs_inclusive.remove(index + 1);
        } else {
            index += 1;
        }
    }
}

#[cfg(test)]
#[path = "seg_tests.rs"]
mod tests;
