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

use super::seg_lnfact::LNFAC;

/// Natural log of 20 (alphabet size for amino acids).
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:2196
///   const double kLn20 = 2.9957322735539909;
const K_LN20: f64 = 2.9957322735539909;

/// NCBI BLAST ln(2) constant.
/// Reference: ncbi-blast/c++/include/algo/blast/core/ncbi_math.h
/// #define NCBIMATH_LN2 0.69314718055994530941723212145818
const NCBIMATH_LN2: f64 = 0.69314718055994530941723212145818;

/// Natural log values for 0.1, 0.2, 0.3, ... 1.0 used by NCBI when window total=10.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1305-1311
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
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1854-1859
/// ```c
/// static double s_lnfact(Int4 n) {
///   if (n < sizeof(lnfact)/sizeof(*lnfact))
///      return lnfact[n];
///   else return ((n+0.5)*log(n) - n + 0.9189385332);
/// }
/// ```
#[inline]
fn s_lnfact(n: usize) -> f64 {
    if n < LNFAC.len() {
        LNFAC[n]
    } else {
        let nf = n as f64;
        (nf + 0.5) * nf.ln() - nf + 0.9189385332
    }
}

/// Map NCBISTDAA letter (0-27) to SEG's 20-letter alphabet index.
///
/// NCBI SEG uses ncbistdaa codes, but only a 20-letter subset is considered valid;
/// everything else (including X) is "bogus".
///
/// NCBI reference (verbatim, blast_seg.c:2197-2216):
/// ```c
/// palpha->alphasize = 20;
/// ...
/// for (c=0, i=0; c<kCharSet; c++)
/// {
///    if (c == 1 || (c >= 3 && c <= 20) || c == 22) {
///       alphaflag[c] = FALSE;
///       alphaindex[c] = i;
///       ++i;
///    } else {
///       alphaflag[c] = TRUE; alphaindex[c] = 20;
///    }
/// }
/// ```
///
/// Valid letters (NCBISTDAA codes): A(1), C..W(3..20), Y(22) => 20 letters.
#[inline]
fn seg_alpha_index_ncbistdaa(letter: u8) -> Option<usize> {
    match letter {
        1 => Some(0),                          // A
        3..=20 => Some((letter - 2) as usize), // C..W => 1..18
        22 => Some(19),                        // Y
        _ => None,                             // bogus: -, B, X, Z, U, *, O, J, etc.
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1324-1350,2190-2217
// ```c
// typedef struct Alpha {
//    Int4 alphabet;
//    Int4 alphasize;
//    double lnalphasize;
//    Int4* alphaindex;
//    unsigned char* alphaflag;
// } Alpha;
// ...
// palpha->alphasize = 20;
// palpha->lnalphasize = kLn20;
// for (c=0, i=0; c<kCharSet; c++) {
//    if (c == 1 || (c >= 3 && c <= 20) || c == 22) {
//       alphaflag[c] = FALSE;
//       alphaindex[c] = i;
//       ++i;
//    } else {
//       alphaflag[c] = TRUE;
//       alphaindex[c] = 20;
//    }
// }
// ```
#[derive(Debug, Clone)]
struct SegAlpha {
    alphasize: usize,
    lnalphasize: f64,
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
            lnalphasize: K_LN20,
            alphaindex,
            alphaflag,
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1335-1350,1541-1705
// ```c
// typedef struct SSequence {
//    struct SSequence* parent;
//    char* seq;
//    Alpha* palpha;
//    Int4 start;
//    Int4 length;
//    Int4 bogus;
//    Boolean punctuation;
//    Int4* composition;
//    Int4* state;
//    double entropy;
// } SSequence;
// ...
// win->start = start;
// win->length = length;
// win->seq = parent->seq + start;
// win->bogus = 0;
// win->punctuation = FALSE;
// win->entropy = -2.;
// s_StateOn(win);
// ```
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
        self.slice().iter().any(|&letter| letter == b'-')
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1631-1657
// ```c
// while ((svi = *sv++) != 0) {
//    if (svi == class && *sv < class) {
//       sv[-1] = svi - 1;
//       break;
//    }
// }
// ...
// for (;;) {
//    if (*sv++ == class) {
//       sv[-1]++;
//       break;
//    }
// }
// ```
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

/// State vector: sorted amino acid counts (descending order, non-zero only)
/// Reference: blast_seg.c uses this for entropy calculation
fn compute_state_vector(counts: &[u32; 20]) -> Vec<i32> {
    let mut sv: Vec<i32> = counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| c as i32)
        .collect();
    sv.sort_by(|a, b| b.cmp(a)); // Sort descending
    sv
}

/// Calculate ln(number of permutations) for a given state vector
/// This is equation 3 from Wootton & Federhen: ln(n!) - sum(ln(count_i!))
/// Reference: blast_seg.c s_LnPerm()
fn ln_perm(sv: &[i32], window_length: i32) -> f64 {
    // NCBI reference (verbatim, blast_seg.c:1867-1882):
    // ```c
    // ans = s_lnfact(window_length);
    // for (i=0; sv[i]!=0; i++) { ans -= s_lnfact(sv[i]); }
    // ```
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
/// Reference: blast_seg.c s_LnAss()
fn ln_ass(sv: &[i32], alphasize: i32) -> f64 {
    // NCBI reference (verbatim, blast_seg.c:1892-1933):
    // ```c
    // ans = lnfact[alphasize];
    // if (sv[0] == 0) return ans;
    // total = alphasize;
    // class = 1;
    // svi = *sv;
    // svim1 = sv[0];
    // for (i=0;; svim1 = svi) {
    //   if (++i==alphasize) { ans -= s_lnfact(class); break; }
    //   else if ((svi = *++sv) == svim1) { class++; continue; }
    //   else {
    //     total -= class;
    //     ans -= s_lnfact(class);
    //     if (svi == 0) { ans -= s_lnfact(total); break; }
    //     else { class = 1; continue; }
    //   }
    // }
    // return ans;
    // ```
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
/// Reference: blast_seg.c s_GetProb()
fn get_prob(sv: &[i32], window_length: i32) -> f64 {
    let alphasize = 20; // Standard amino acid alphabet
    let totseq = (window_length as f64) * K_LN20;

    let ans1 = ln_ass(sv, alphasize);
    // NCBI: guard ans2 computation with ans1 > -100000 and sv[0] != INT4_MIN
    let ans2 = if ans1 > -100000.0 {
        ln_perm(sv, window_length)
    } else {
        0.0
    };

    ans1 + ans2 - totseq
}

/// NCBI SEG entropy (Shannon entropy, bits) implementation.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1592-1624 (s_Entropy)
fn entropy_from_state_vector(state: &[i32]) -> f64 {
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
            ent += (count as f64) * ((count as f64) / total_f).ln() / NCBIMATH_LN2;
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
        // NCBI parameter normalization (blast_seg.c:s_SegParametersCheck):
        // ```c
        // if (sparamsp->window <= 0) sparamsp->window = 12;
        // if (sparamsp->locut < 0.0) sparamsp->locut = 0.0;
        // if (sparamsp->hicut < 0.0) sparamsp->hicut = 0.0;
        // if (sparamsp->locut > sparamsp->hicut)
        //     sparamsp->hicut = sparamsp->locut;
        // ```
        let window = if window > 0 { window } else { 12 };
        let mut locut = if locut >= 0.0 { locut } else { 0.0 };
        let mut hicut = if hicut >= 0.0 { hicut } else { 0.0 };
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
}

impl SegMasker {
    /// Create a new SEG masker with the given parameters
    pub fn new(window: usize, locut: f64, hicut: f64) -> Self {
        let params = SegParams::new(window, locut, hicut);

        // downset = (window+1)/2 - 1, upset = window - downset
        // Reference: blast_seg.c:2050-2051
        let downset = (params.window + 1) / 2 - 1;
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
        }
    }

    /// Create a SEG masker with default parameters (window=12, locut=2.2, hicut=2.5)
    pub fn with_defaults() -> Self {
        Self::new(12, 2.2, 2.5)
    }

    /// Create with full parameters including maxbogus and maxtrim
    pub fn with_params(params: &SegParams) -> Self {
        let downset = (params.window + 1) / 2 - 1;
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
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1541-1604,1745-1780
    // ```c
    // win = s_OpenWin(seq, 0, window);
    // s_EntropyOn(win);
    // ...
    // if (win->bogus > maxbogus) { H[i] = -1.; ... }
    // H[i] = win->entropy;
    // ```
    fn calculate_entropy(&self, window: &[u8]) -> (f64, usize) {
        let Some(mut win) = SegWindow::open(window, 0, window.len(), &self.alpha) else {
            return (0.0, 0);
        };
        win.entropy_on();
        (win.entropy, win.bogus as usize)
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1745-1803
    // ```c
    // win = s_OpenWin(seq, 0, window);
    // s_EntropyOn(win);
    // first = downset;
    // last = seq->length - upset;
    // for (i=first; i<=last; i++) {
    //    if (seq->punctuation && s_HasDash(win)) { H[i] = -1.; s_ShiftWin1(win); continue; }
    //    if (win->bogus > maxbogus) { H[i] = -1.; s_ShiftWin1(win); continue; }
    //    H[i] = win->entropy;
    //    s_ShiftWin1(win);
    // }
    // ```
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
    /// Reference: blast_seg.c:1813-1825 (s_FindLow)
    fn find_low(&self, i: usize, limit: usize, h: &[f64]) -> usize {
        // NCBI reference (verbatim):
        // ```c
        // for (j=i; j>=limit; j--) {
        //   if (H[j]==-1.0) break;
        //   if (H[j]>hicut) break;
        // }
        // return (j+1);
        // ```
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
    /// Reference: blast_seg.c:1836-1848 (s_FindHigh)
    fn find_high(&self, i: usize, limit: usize, h: &[f64]) -> usize {
        // NCBI reference (verbatim):
        // ```c
        // for (j=i; j<=limit; j++) {
        //   if (H[j]==-1.0) break;
        //   if (H[j]>hicut) break;
        // }
        // return (j-1);
        // ```
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

    /// Calculate probability for a subsequence
    /// Returns the K2 complexity probability (lower = more low-complexity)
    fn get_subseq_prob(&self, seq: &[u8]) -> f64 {
        if seq.is_empty() {
            return 1.0;
        }

        let mut counts = [0u32; 20];
        let mut valid_count = 0usize;

        for &aa in seq {
            if aa < 20 {
                counts[aa as usize] += 1;
                valid_count += 1;
            }
        }

        if valid_count == 0 {
            return 1.0;
        }

        let sv = compute_state_vector(&counts);
        get_prob(&sv, valid_count as i32)
    }

    /// NCBI s_Trim: trim [leftend..=rightend] to minimize s_GetProb.
    /// Reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1974-2018
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

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1987-2007
        // ```c
        // lend = 0;
        // rend = seq->length - 1;
        // minlen = 1;
        // maxtrim = sparamsp->maxtrim;
        // if ((seq->length-maxtrim)>minlen) minlen = seq->length-maxtrim;
        // minprob = 1.;
        // for (len=seq->length; len>minlen; len--) {
        //    Boolean shift = TRUE;
        //    Int4 i = 0;
        //    SSequence* win = s_OpenWin(seq, 0, len);
        //    while (shift) {
        //       prob = s_GetProb(win->state, len, win->palpha);
        //       if (prob<minprob) { minprob = prob; lend = i; rend = len + i - 1; }
        //       shift = s_ShiftWin1(win);
        //       i++;
        //    }
        // }
        // ```
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

    /// Mask a sequence and return the list of masked intervals
    /// Reference: blast_seg.c:2030-2116 (s_SegSeq)
    pub fn mask_sequence(&self, seq: &[u8]) -> Vec<MaskedInterval> {
        if seq.len() < self.window {
            return Vec::new();
        }

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
                    //   leftend = loi - downset;
                    //   rightend = hii + upset - 1;
                    let mut leftend = loi - masker.downset;
                    let mut rightend = hii + masker.upset - 1; // inclusive

                    // Trim to minimize probability (NCBI s_Trim)
                    (leftend, rightend) = masker.trim_segment(seq, leftend, rightend);

                    // NCBI recursion for trigger window in left trim:
                    //   if (i+upset-1 < leftend) { recurse on [lend..rend] }
                    if i + masker.upset - 1 < leftend {
                        let lend = loi - masker.downset;
                        if lend < leftend {
                            let rend = leftend - 1;
                            if rend < seq.len() && lend <= rend {
                                seg_seq(masker, &seq[lend..=rend], offset + lend, segs);
                            }
                        }
                    }

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:2104-2108
                    // ```c
                    // seg = (SSeg*) calloc(1, sizeof(SSeg));
                    // seg->begin = leftend + offset;
                    // seg->end = rightend + offset;
                    // seg->next = *segs;
                    // *segs = seg;
                    // ```
                    segs.insert(0, (leftend + offset, rightend + offset));

                    // NCBI:
                    //   i = MIN(hii, rightend+downset);
                    //   lowlim = i + 1;
                    i = hii.min(rightend + masker.downset);
                    lowlim = i + 1;
                    i += 1;
                } else {
                    i += 1;
                }
            }
        }

        seg_seq(self, seq, 0, &mut segs_inclusive);

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1142-1154
        // ```c
        // sparamsp = SegParametersNewAa();
        // sparamsp->overlaps = TRUE;
        // status = SeqBufferSeg(sequence, length, offset, sparamsp,
        //                       seqloc_retval);
        // ```
        //
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:2125-2149
        // ```c
        // if (sparamsp->overlaps)
        //    s_MergeSegs(seqwin, segs);
        // ```
        merge_segs_inclusive(seq.len(), &mut segs_inclusive);

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:2300-2319
        // ```c
        // s_SegsToBlastSeqLoc(segs, offset, seg_locs);
        // ```
        segs_inclusive.reverse();

        segs_inclusive
            .into_iter()
            .map(|(b, e)| MaskedInterval::new(b, e.saturating_add(1)))
            .collect()
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:2125-2149
// ```c
// s_MergeSegs(SSequence* seq, SSeg* segs)
// {
//    ...
//    while (nextseg!=NULL) {
//       if (seg->begin - nextseg->end - 1 < hilenmin) {
//          if (seg->end < nextseg->end) seg->end = nextseg->end;
//          if (seg->begin > nextseg->begin) seg->begin = nextseg->begin;
//          seg->next = nextseg->next;
//       } else {
//          seg = nextseg;
//       }
//       nextseg = seg->next;
//    }
//    ...
// }
// ```
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
        if seg_begin <= next_end.saturating_add(1) {
            segs_inclusive[index].0 = seg_begin.min(next_begin);
            segs_inclusive[index].1 = seg_end.max(next_end);
            segs_inclusive.remove(index + 1);
        } else {
            index += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::matrix::aa_char_to_ncbistdaa;

    fn fasta_sequence_by_id(contents: &str, wanted_id: &str) -> Vec<u8> {
        let mut current_id = None;
        let mut residues = Vec::new();
        for line in contents.lines() {
            if let Some(rest) = line.strip_prefix('>') {
                let found = rest.split_whitespace().next();
                if current_id == Some(wanted_id) {
                    break;
                }
                current_id = found;
                continue;
            }
            if current_id == Some(wanted_id) {
                residues.extend(
                    line.bytes()
                        .map(|aa| aa_char_to_ncbistdaa(aa.to_ascii_uppercase())),
                );
            }
        }
        residues
    }

    #[test]
    fn test_seg_params_default() {
        let params = SegParams::default();
        assert_eq!(params.window, 12);
        assert_eq!(params.locut, 2.2);
        assert_eq!(params.hicut, 2.5);
    }

    #[test]
    fn test_seg_params_validation() {
        // NCBI: negative values become 0.0, not defaults
        // if (sparamsp->locut < 0.0) sparamsp->locut = 0.0;
        // if (sparamsp->hicut < 0.0) sparamsp->hicut = 0.0;
        // if (sparamsp->locut > sparamsp->hicut) sparamsp->hicut = sparamsp->locut;
        let params = SegParams::new(0, -1.0, 1.0);
        assert_eq!(params.window, 12); // Should default to 12
        assert_eq!(params.locut, 0.0); // Negative becomes 0.0
        assert_eq!(params.hicut, 1.0); // 1.0 > 0.0, so stays at 1.0
    }

    #[test]
    fn test_calculate_entropy() {
        let masker = SegMasker::with_defaults();

        // Low complexity: all same amino acid (entropy = 0)
        // NCBISTDAA: A=1 (not 0, which is '-')
        let low_complex = vec![1u8; 12]; // All alanine (NCBISTDAA code 1)
        let (entropy_low, _) = masker.calculate_entropy(&low_complex);
        assert!(
            entropy_low < 1.0,
            "Low complexity should have low entropy: got {}",
            entropy_low
        );

        // High complexity: all different amino acids
        // Use valid NCBISTDAA codes: 1=A, 3-20=C..W, 22=Y
        let high_complex: Vec<u8> = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13].to_vec();
        let (entropy_high, _) = masker.calculate_entropy(&high_complex);
        assert!(
            entropy_high > 2.0,
            "High complexity should have high entropy: got {}",
            entropy_high
        );
    }

    #[test]
    fn test_mask_simple_repeat() {
        let masker = SegMasker::with_defaults();

        // Simple repeat sequence should be masked
        // NCBISTDAA: A=1 (not 0, which is '-')
        let seq: Vec<u8> = vec![1u8; 100]; // All alanine (NCBISTDAA code 1)
        let intervals = masker.mask_sequence(&seq);

        assert!(
            !intervals.is_empty(),
            "Poly-alanine sequence should be masked"
        );
        if !intervals.is_empty() {
            assert!(intervals[0].start < intervals[0].end);
        }
    }

    #[test]
    fn test_mask_short_sequence() {
        let masker = SegMasker::with_defaults();

        // Very short sequences should return empty (shorter than window)
        // NCBISTDAA: A=1
        let seq = vec![1u8; 5];
        let intervals = masker.mask_sequence(&seq);
        assert!(intervals.is_empty());
    }

    #[test]
    fn test_mask_complex_sequence() {
        let masker = SegMasker::with_defaults();

        // Complex sequence should not be heavily masked
        let seq: Vec<u8> = (0..100).map(|i| (i % 20) as u8).collect();
        let intervals = masker.mask_sequence(&seq);

        // Complex sequence should have few or no masked regions
        let total_masked: usize = intervals.iter().map(|i| i.end - i.start).sum();
        assert!(
            total_masked < seq.len() / 2,
            "Complex sequence should not be heavily masked"
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1000-1004
    // ```c
    // /** NCBIstdaa encoding for 'X' character */
    // #define BLASTP_MASK_RESIDUE 21
    // /** Default instructions and mask residue for SEG filtering */
    // #define BLASTP_MASK_INSTRUCTIONS "S 10 1.8 2.1"
    // ```
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_bdt63528() {
        let fasta = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/PajaWSV.faa"
        ));
        let seq = fasta_sequence_by_id(fasta, "BDT63528.1");
        let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![
                MaskedInterval::new(217, 223),
                MaskedInterval::new(568, 580),
                MaskedInterval::new(762, 774),
            ]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1000-1004
    // ```c
    // /** NCBIstdaa encoding for 'X' character */
    // #define BLASTP_MASK_RESIDUE 21
    // /** Default instructions and mask residue for SEG filtering */
    // #define BLASTP_MASK_INSTRUCTIONS "S 10 1.8 2.1"
    // ```
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_bdt63573() {
        let fasta = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/PajaWSV.faa"
        ));
        let seq = fasta_sequence_by_id(fasta, "BDT63573.1");
        let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![
                MaskedInterval::new(146, 176),
                MaskedInterval::new(192, 209),
                MaskedInterval::new(214, 227),
                MaskedInterval::new(408, 420),
                MaskedInterval::new(587, 601),
                MaskedInterval::new(605, 619),
            ]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1000-1004
    // ```c
    // /** NCBIstdaa encoding for 'X' character */
    // #define BLASTP_MASK_RESIDUE 21
    // /** Default instructions and mask residue for SEG filtering */
    // #define BLASTP_MASK_INSTRUCTIONS "S 10 1.8 2.1"
    // ```
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_bdt63533() {
        let fasta = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/PajaWSV.faa"
        ));
        let seq = fasta_sequence_by_id(fasta, "BDT63533.1");
        let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![
                MaskedInterval::new(293, 303),
                MaskedInterval::new(412, 427),
                MaskedInterval::new(873, 882),
                MaskedInterval::new(1077, 1087),
            ]
        );
    }

    // NCBI reference: /home/kawato/micromamba/bin/segmasker
    // Command:
    // `segmasker -in /tmp/BDT63510.faa -outfmt interval -window 10 -locut 1.8 -hicut 2.1`
    //
    // `segmasker` interval output is `0-based closed`; Rust stores merged
    // `MaskedInterval` values as `0-based half-open`.
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1415-1438
    // ```c
    // #define BLASTP_MASK_INSTRUCTIONS "S 10 1.8 2.1"
    // static int
    // s_DoSegSequenceData(BlastCompo_SequenceData * seqData,
    //                     EBlastProgramType program,
    //                     Boolean * pSequenceIsBiased)
    // ```
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_bdt63510() {
        let fasta = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/PajaWSV.faa"
        ));
        let seq = fasta_sequence_by_id(fasta, "BDT63510.1");
        let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![
                MaskedInterval::new(207, 228),
                MaskedInterval::new(257, 266),
                MaskedInterval::new(591, 610),
                MaskedInterval::new(789, 797),
            ]
        );
    }

    // NCBI reference: /home/kawato/micromamba/bin/segmasker
    // Command:
    // `segmasker -in - -outfmt interval < tests/fasta/WSSV.faa`
    //
    // `segmasker` interval output is `0-based closed`; Rust stores merged
    // `MaskedInterval` values as `0-based half-open`.
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:337-370
    // ```c
    // BlastSetUp_Filter(..., query_blk, ...);
    // query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
    // query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
    // ```
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_yp_009220537_1() {
        let fasta = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fasta/WSSV.faa"));
        let seq = fasta_sequence_by_id(fasta, "YP_009220537.1");
        let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![MaskedInterval::new(319, 335), MaskedInterval::new(920, 931),]
        );
    }

    // NCBI reference: /home/kawato/micromamba/bin/segmasker
    // Command:
    // `segmasker -in - -outfmt interval < tests/fasta/WSSV.faa`
    //
    // `segmasker` interval output is `0-based closed`; Rust stores merged
    // `MaskedInterval` values as `0-based half-open`.
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:337-370
    // ```c
    // BlastSetUp_Filter(..., query_blk, ...);
    // query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
    // query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
    // ```
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_yp_009220568_1() {
        let fasta = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fasta/WSSV.faa"));
        let seq = fasta_sequence_by_id(fasta, "YP_009220568.1");
        let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![
                MaskedInterval::new(12, 28),
                MaskedInterval::new(94, 129),
                MaskedInterval::new(280, 293),
                MaskedInterval::new(452, 536),
                MaskedInterval::new(649, 655),
                MaskedInterval::new(709, 727),
                MaskedInterval::new(1041, 1054),
                MaskedInterval::new(1139, 1149),
            ]
        );
    }

    // NCBI reference: /home/kawato/micromamba/bin/segmasker
    // Command:
    // `segmasker -in - -outfmt interval < tests/fasta/WSSV.faa`
    //
    // `segmasker` interval output is `0-based closed`; Rust stores merged
    // `MaskedInterval` values as `0-based half-open`.
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:337-370
    // ```c
    // BlastSetUp_Filter(..., query_blk, ...);
    // query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
    // query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
    // ```
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_yp_009220574_1() {
        let fasta = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fasta/WSSV.faa"));
        let seq = fasta_sequence_by_id(fasta, "YP_009220574.1");
        let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![
                MaskedInterval::new(32, 39),
                MaskedInterval::new(364, 382),
                MaskedInterval::new(774, 791),
                MaskedInterval::new(827, 839),
                MaskedInterval::new(915, 923),
                MaskedInterval::new(1103, 1112),
                MaskedInterval::new(1227, 1246),
                MaskedInterval::new(1267, 1279),
            ]
        );
    }

    // NCBI reference: /home/kawato/micromamba/bin/segmasker
    // Command:
    // `segmasker -in /tmp/bdv02435_query.faa -outfmt interval`
    //
    // `segmasker` interval output is `0-based closed`; Rust stores merged
    // `MaskedInterval` values as `0-based half-open`.
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_bdv02435_1_query_defaults() {
        let fasta = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/AP027131.faa"
        ));
        let seq = fasta_sequence_by_id(fasta, "BDV02435.1");
        let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![
                MaskedInterval::new(3, 26),
                MaskedInterval::new(263, 278),
                MaskedInterval::new(348, 362),
                MaskedInterval::new(427, 443),
                MaskedInterval::new(452, 466),
                MaskedInterval::new(521, 542),
                MaskedInterval::new(623, 640),
                MaskedInterval::new(697, 708),
                MaskedInterval::new(732, 760),
                MaskedInterval::new(1021, 1036),
                MaskedInterval::new(1214, 1227),
                MaskedInterval::new(1273, 1295),
                MaskedInterval::new(1418, 1431),
                MaskedInterval::new(1477, 1499),
                MaskedInterval::new(2489, 2509),
                MaskedInterval::new(2593, 2613),
                MaskedInterval::new(2697, 2717),
                MaskedInterval::new(3814, 3841),
                MaskedInterval::new(4087, 4107),
                MaskedInterval::new(4124, 4139),
                MaskedInterval::new(4142, 4157),
                MaskedInterval::new(4204, 4224),
                MaskedInterval::new(4241, 4256),
                MaskedInterval::new(4347, 4365),
                MaskedInterval::new(4516, 4533),
                MaskedInterval::new(4568, 4578),
            ]
        );
    }

    // NCBI reference: /home/kawato/micromamba/bin/segmasker
    // Command:
    // `segmasker -in /tmp/BDT63134.faa -outfmt interval`
    //
    // `segmasker` interval output is `0-based closed`; Rust stores merged
    // `MaskedInterval` values as `0-based half-open`.
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1415-1438
    // ```c
    // #define BLASTP_MASK_INSTRUCTIONS "S 10 1.8 2.1"
    // static int
    // s_DoSegSequenceData(BlastCompo_SequenceData * seqData,
    //                     EBlastProgramType program,
    //                     Boolean * pSequenceIsBiased)
    // ```
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_bdt63134() {
        let fasta = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/SicyWSV.faa"
        ));
        let seq = fasta_sequence_by_id(fasta, "BDT63134.1");
        let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![
                MaskedInterval::new(251, 264),
                MaskedInterval::new(311, 321),
                MaskedInterval::new(348, 363),
                MaskedInterval::new(382, 402),
            ]
        );
    }

    // NCBI reference: /home/kawato/micromamba/bin/segmasker
    // Command:
    // `segmasker -in /tmp/ap_subjects.faa -outfmt interval -window 10 -locut 1.8 -hicut 2.1`
    //
    // `segmasker` interval output is `0-based closed`; Rust stores merged
    // `MaskedInterval` values as `0-based half-open`.
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_wp_025208654_1() {
        let fasta = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/NZ_CP006932.faa"
        ));
        let seq = fasta_sequence_by_id(fasta, "WP_025208654.1");
        let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(intervals, vec![MaskedInterval::new(137, 156),]);
    }

    // NCBI reference: /home/kawato/micromamba/bin/segmasker
    // Command:
    // `segmasker -in /tmp/ap_subjects.faa -outfmt interval -window 10 -locut 1.8 -hicut 2.1`
    //
    // `segmasker` interval output is `0-based closed`; Rust stores merged
    // `MaskedInterval` values as `0-based half-open`.
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_wp_025208655_1() {
        let fasta = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/NZ_CP006932.faa"
        ));
        let seq = fasta_sequence_by_id(fasta, "WP_025208655.1");
        let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![
                MaskedInterval::new(148, 167),
                MaskedInterval::new(1270, 1280),
            ]
        );
    }

    // NCBI reference: /home/kawato/micromamba/bin/segmasker
    // Command:
    // `segmasker -in /tmp/ap_targets.faa -outfmt interval -window 10 -locut 1.8 -hicut 2.1`
    //
    // `segmasker` interval output is `0-based closed`; Rust stores merged
    // `MaskedInterval` values as `0-based half-open`.
    #[test]
    fn test_mask_sequence_matches_ncbi_segmasker_bdv02435_1() {
        let fasta = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/AP027131.faa"
        ));
        let seq = fasta_sequence_by_id(fasta, "BDV02435.1");
        let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);

        assert_eq!(
            intervals,
            vec![
                MaskedInterval::new(8, 19),
                MaskedInterval::new(174, 183),
                MaskedInterval::new(629, 640),
                MaskedInterval::new(697, 708),
                MaskedInterval::new(908, 918),
                MaskedInterval::new(1273, 1287),
                MaskedInterval::new(1477, 1491),
                MaskedInterval::new(1946, 1956),
                MaskedInterval::new(4091, 4105),
                MaskedInterval::new(4208, 4222),
                MaskedInterval::new(4241, 4256),
                MaskedInterval::new(4347, 4365),
                MaskedInterval::new(4568, 4578),
            ]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:2125-2149
    // ```c
    // if (seg->begin - nextseg->end - 1 < hilenmin) {
    //     ...
    // }
    // ```
    #[test]
    fn test_merge_segs_inclusive_matches_ncbi_overlap_merge() {
        let mut segs = vec![(20, 30), (10, 25), (0, 5)];
        merge_segs_inclusive(64, &mut segs);
        assert_eq!(segs, vec![(10, 30), (0, 5)]);
    }
}
