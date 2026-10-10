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
    // EXPERIMENT (LOSAT_X_SEGMEMO): the entropy is a pure function of the
    // state vector, and a 12-letter window has at most a few hundred distinct
    // state vectors, so remember the value computed for each one.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1598-1603,1616-1621
    // ```c
    //    total = 0;
    //    for (i=0; sv[i]!=0; i++)
    //      {
    //       total += sv[i];
    //      }
    //    if (total==0) return(0.);
    // ...
    //     for (i=0; sv[i]!=0; i++)
    //         {
    //             ent += ((double)sv[i])*log(((double)sv[i])/(double)total)/NCBIMATH_LN2;
    //         }
    //    }
    //    ent = fabs(ent/(double)total);
    // ```
    // Dispatch point of LOSAT_X_SEGMEMO: the reference path is `entropy_from_state_vector_reference`, the port
    // of this function. The memo is keyed by the sorted counts of the state vector, which is all that `s_Entropy`
    // reads, so a repeated vector gets the value the first call returned.
    if x_seg_memo() {
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
            return X_ENTROPY_MEMO.with(|memo| {
                let mut memo = memo.borrow_mut();
                let mask = memo.len() - 1;
                let mut slot = (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 40) as usize & mask;
                loop {
                    let (stored_key, bits) = memo[slot];
                    if stored_key == key {
                        return f64::from_bits(bits);
                    }
                    if stored_key == u64::MAX {
                        let value = entropy_from_state_vector_reference(state);
                        memo[slot] = (key, value.to_bits());
                        return value;
                    }
                    slot = (slot + 1) & mask;
                }
            });
        }
    }
    entropy_from_state_vector_reference(state)
}

// No NCBI counterpart: reads LOSAT_X_SEGMEMO once; it does not change any value NCBI computes.
fn x_seg_memo() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_SEGMEMO").is_some())
}

thread_local! {
    // No NCBI counterpart: memo table of entropy values per state vector (`s_Entropy` is a pure function of the vector); it does not change any value NCBI computes.
    // Open-addressing table. Only state vectors summing to at most 15 are
    // stored, and there are fewer than 700 of those, so it never fills.
    static X_ENTROPY_MEMO: std::cell::RefCell<Vec<(u64, u64)>> =
        std::cell::RefCell::new(vec![(u64::MAX, 0); 4096]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1598-1606,1616-1621
// ```c
//    total = 0;
//    for (i=0; sv[i]!=0; i++)
//      {
//       total += sv[i];
//      }
//    if (total==0) return(0.);
//    ent = 0.0;
//    if (total == 10)
// ...
//     for (i=0; sv[i]!=0; i++)
//         {
//             ent += ((double)sv[i])*log(((double)sv[i])/(double)total)/NCBIMATH_LN2;
//         }
//    }
//    ent = fabs(ent/(double)total);
// ```
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
    /// Whether a recursion on the left part of a trimmed segment keeps only the head of
    /// its result, as NCBI's `s_SegSeq` does (true, the default; see `mask_sequence`).
    ncbi_left_segments: bool,
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
            ncbi_left_segments: true,
        }
    }

    /// The masker that keeps every segment of the left recursion (LOSAT's SEG before S08).
    /// BLASTX keeps it until it is integrated (plan DW-10, SX); every other program uses
    /// NCBI's behavior.
    pub fn keeping_all_left_segments(mut self) -> Self {
        self.ncbi_left_segments = false;
        self
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
            ncbi_left_segments: true,
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
        // EXPERIMENT (LOSAT_X_SEGFAST / LOSAT_X_SEGSHADOW): see `x_mask_sequence`.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2050-2053,2062-2064
        // ```c
        //    downset = (window+1)/2 - 1;
        //    upset = window - downset;
        //    H = s_SeqEntropy(seq, window, sparamsp->maxbogus);
        // ...
        //    for (i=first; i<=last; i++)
        //    {
        //       if (H[i] <= locut && H[i] != -1.0)
        // ```
        // Dispatch point of LOSAT_X_SEGFAST / LOSAT_X_SEGSHADOW: the reference path is `mask_sequence_reference`, the
        // port of `s_SegSeq`. The fast path keeps the control flow of this loop and replaces only the entropy
        // and trim computations.
        match x_seg_fast() {
            1 if self.x_fast_applies() => return self.x_mask_sequence(seq),
            2 if self.x_fast_applies() => {
                let fast = self.x_mask_sequence(seq);
                let reference = self.mask_sequence_reference(seq);
                assert!(
                    fast == reference,
                    "LOSAT_X_SEGSHADOW: fast SEG differs from the reference on a {}-residue sequence",
                    seq.len()
                );
                return fast;
            }
            _ => {}
        }
        self.mask_sequence_reference(seq)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2030-2032,2053-2060
    // ```c
    // s_SegSeq(SSequence* seq, SegParameters* sparamsp, SSeg **segs,
    //                    Int4 offset)
    // {
    // ...
    //    H = s_SeqEntropy(seq, window, sparamsp->maxbogus);
    //    if (H == NULL)
    //       return status;
    //    first = downset;
    //    last = seq->length - upset;
    //    lowlim = first;
    // ```
    // The original body of `mask_sequence`, moved here unchanged; it is the port of this function.
    fn mask_sequence_reference(&self, seq: &[u8]) -> Vec<MaskedInterval> {
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
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:2086-2101
                                // ```c
                                //             SSeg *leftsegs = (SSeg*) NULL;
                                // ...
                                //             status = s_SegSeq(leftseq, sparamsp, &leftsegs, offset+lend);
                                // ...
                                //             if (leftsegs!=NULL)
                                //             {
                                //                leftsegs->next = *segs;
                                //                *segs = leftsegs;
                                //             }
                                // ```
                                // `leftsegs->next = *segs` overwrites the link of the head of
                                // the recursion's list, so when the recursion found more than
                                // one segment only its head (the one found last) is kept; the
                                // others are dropped (and leaked). This deterministic result
                                // is reproduced.
                                if masker.ncbi_left_segments {
                                    let mut left_segs: Vec<(usize, usize)> = Vec::new();
                                    seg_seq(
                                        masker,
                                        &seq[lend..=rend],
                                        offset + lend,
                                        &mut left_segs,
                                    );
                                    if let Some(&head) = left_segs.first() {
                                        segs.insert(0, head);
                                    }
                                } else {
                                    seg_seq(masker, &seq[lend..=rend], offset + lend, segs);
                                }
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
// ---------------------------------------------------------------------------
// EXPERIMENT (LOSAT_X_SEGFAST): the same SEG, with the window kept as counts.
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
// ```c
//     for (letter = nel = 0; letter < alphasize; ++letter) {
//         if ((c = win->composition[letter]) == 0)
//             continue;
//         win->state[nel++] = c;
//     }
// ...
//     qsort(win->state, nel, sizeof(win->state[0]), s_StateCmp);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1642-1646,1658-1661
// ```c
//     while ((svi = *sv++) != 0) {
//         if (svi == class && *sv < class) {
//             sv[-1] = svi - 1;
//             break;
//         }
// ...
//     for (;;) {
//         if (*sv++ == class) {
//             sv[-1]++;
//             break;
// ```
// The state vector is the sorted list of letter counts. A shift by one letter lowers the count of the
// leaving letter and raises the count of the entering one (`s_DecrementSV`, `s_IncrementSV`). The code below keeps a histogram of the counts instead, which holds the
// same information, and updates it per shift. `x_seg_fast` itself reads LOSAT_X_SEGFAST and LOSAT_X_SEGSHADOW once.
fn x_seg_fast() -> u8 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_SEGSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_SEGFAST").is_some() {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: size bound of the histogram arrays of `x_trim_segment`; longer segments use the reference trim; it does not change any value NCBI computes.
/// Longest segment handled by the count histogram of `x_trim_segment`.
const X_TRIM_MAX: usize = 127;

thread_local! {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1593-1593,1599-1602,1621-1623
    // ```c
    // s_Entropy(Int4* sv)
    // ...
    //    for (i=0; sv[i]!=0; i++)
    //      {
    //       total += sv[i];
    //      }
    // ...
    //    ent = fabs(ent/(double)total);
    //    return(ent);
    // ```
    // No NCBI counterpart: memo table of `s_Entropy` values per count histogram. The entropy depends only on the
    // counts, so reusing it changes no value.
    // Entropy by count histogram (4 bits per count value, windows up to 15).
    static X_ENTROPY_BY_HISTOGRAM: std::cell::RefCell<Vec<(u64, u64)>> =
        std::cell::RefCell::new(vec![(u64::MAX, 0); 2048]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1598-1603,1616-1621
// ```c
//    total = 0;
//    for (i=0; sv[i]!=0; i++)
//      {
//       total += sv[i];
//      }
//    if (total==0) return(0.);
// ...
//     for (i=0; sv[i]!=0; i++)
//         {
//             ent += ((double)sv[i])*log(((double)sv[i])/(double)total)/NCBIMATH_LN2;
//         }
//    }
//    ent = fabs(ent/(double)total);
// ```
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
// ```c
//    totseq = ((double) total) * palpha->lnalphasize;
//    ans1 = s_LnAss(sv, palpha->alphasize);
//    if (ans1 > -100000.0 && sv[0] != INT4_MIN)
//    {
//     ans2 = s_LnPerm(sv, total);
//    }
// ...
//    ans = ans1 + ans2 - totseq;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1901-1905,1909-1912,1919-1923
// ```c
//     ans = lnfact[alphasize];
//     if (sv[0] == 0)
//         return ans;
//     total = alphasize;
// ...
//     for (i=0;; svim1 = svi) {
//             if (++i==alphasize) {
//                 ans -= s_lnfact(class);
//             break;
// ...
//             total -= class;
//             ans -= s_lnfact(class);
//             if (svi == 0) {
//                 ans -= s_lnfact(total);
//                 break;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1874-1879
// ```c
//    ans = s_lnfact(window_length);
//    for (i=0; sv[i]!=0; i++)
//      {
//       ans -= s_lnfact(sv[i]);
//      }
// ```
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
    // ```c
    //     while (seq < seqmax) {
    //         letter = *seq++;
    //         if (!alphaflag[letter])
    //             comp[alphaindex[letter]]++;
    //                 else
    //                         win->bogus++;
    //     }
    // ```
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
    // ```c
    //    win = s_OpenWin(seq, 0, window);
    //    s_EntropyOn(win);
    // ...
    //    first = downset;
    //    last = seq->length - upset;
    //    for (i=first; i<=last; i++)
    //      {
    // ...
    //       if (win->bogus > maxbogus)
    //         {
    //          H[i] = -1.;
    //          s_ShiftWin1(win);
    //          continue;
    //         }
    //       H[i] = win->entropy;
    //       s_ShiftWin1(win);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1695-1704
    // ```c
    //     if (!alphaflag[j])
    //         s_DecrementSV(win->state, comp[alphaindex[j]]--);
    //     else win->bogus--;
    //     j = (Uint1) win->seq[length];   /* prevent sign-extension */
    //     ++win->seq;
    //     if (!alphaflag[j])
    //         s_IncrementSV(win->state, comp[alphaindex[j]]++);
    //     else win->bogus++;
    // ```
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
        // ```c
        //    ent = fabs(ent/(double)total);
        //    return(ent);
        // ```
        // The entropy of a window depends only on its count histogram, so the value for a seen histogram is reused.
        X_ENTROPY_BY_HISTOGRAM.with(|memo| {
            let mut memo = memo.borrow_mut();
            let slots = memo.len() - 1;
            let mut filled = 0usize;
            for start in 0..windows {
                h[first + start] = if bogus > maxbogus {
                    -1.0
                } else {
                    let mut slot = (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 40) as usize & slots;
                    loop {
                        let (stored, bits) = memo[slot];
                        if stored == key {
                            break f64::from_bits(bits);
                        }
                        if stored == u64::MAX {
                            let value = x_entropy_of_histogram(key);
                            // Fewer than 700 histograms exist for windows up to
                            // 15 letters; the table is never close to full.
                            memo[slot] = (key, value.to_bits());
                            filled += 1;
                            debug_assert!(filled < slots);
                            break value;
                        }
                        slot = (slot + 1) & slots;
                    }
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
        });
        h
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1988-1992,1998-2001,2007-2008,2013-2014
    // ```c
    //    if ((seq->length-maxtrim)>minlen)
    //         minlen = seq->length-maxtrim;
    //    minprob = 1.;
    //    for (len=seq->length; len>minlen; len--)
    // ...
    //       while (shift)
    //       {
    //          prob = s_GetProb(win->state, len, win->palpha);
    //          if (prob<minprob)
    // ...
    //          shift = s_ShiftWin1(win);
    //          i++;
    // ...
    //    *leftend = *leftend + lend;
    //    *rightend = *rightend - (seq->length - rend - 1);
    // ```
    // Windows of each length are visited left to right, with the smallest `prob` kept on a strict `<` as here.
    // Each window is updated from its neighbour by one removal and one addition, and the sums in `s_LnAss` and
    // `s_LnPerm` use the histogram in NCBI's order.
    /// `trim_segment` with every window kept as a count histogram.
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
        // ```c
        //     while ((svi = *sv++) != 0) {
        //         if (svi == class && *sv < class) {
        //             sv[-1] = svi - 1;
        //             break;
        //         }
        // ...
        //     for (;;) {
        //         if (*sv++ == class) {
        //             sv[-1]++;
        //             break;
        // ```
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
    // ```c
    //    for (i=first; i<=last; i++)
    //    {
    //       if (H[i] <= locut && H[i] != -1.0)
    //         {
    //          Int4 loi = s_FindLow(i, lowlim, hicut, H);
    //          Int4 hii = s_FindHigh(i, last, hicut, H);
    //          SSequence* temp_seq = NULL;
    //          leftend = loi - downset;
    //          rightend = hii + upset - 1;
    // ...
    //          temp_seq = s_OpenWin(seq, leftend, rightend-leftend+1);
    //          status = s_Trim(temp_seq, &leftend, &rightend, sparamsp);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2081-2084,2105-2111
    // ```c
    //          if (i+upset-1<leftend)   /* check for trigger window in left trim */
    //          {
    //             Int4 lend = loi - downset;
    //             Int4 rend = leftend - 1;
    // ...
    //          seg = (SSeg*) calloc(1, sizeof(SSeg));
    //          seg->begin = leftend + offset;
    //          seg->end = rightend + offset;
    //          seg->next = *segs;
    //          *segs = seg;
    //          i = MIN(hii, rightend+downset);
    //          lowlim = i + 1;
    // ```
    // The recursion, the left-trim handling and the order of the segments are the same as in
    // `mask_sequence_reference`; only the entropy array and the trim use the histogram.
    /// `mask_sequence` through `x_entropy_array` and `x_trim_segment`.
    fn x_mask_sequence(&self, seq: &[u8]) -> Vec<MaskedInterval> {
        let index = self.x_letter_index();
        let mut segs_inclusive: Vec<(usize, usize)> = Vec::new();

        // NCBI blast_seg.c:2030-2116 (s_SegSeq); the control flow of `seg_seq` in
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

                    (leftend, rightend) = masker.x_trim_segment(seq, leftend, rightend, index);

                    if i + masker.upset - 1 < leftend {
                        let lend = loi - masker.downset;
                        if lend < leftend {
                            let rend = leftend - 1;
                            if rend < seq.len() && lend <= rend {
                                if masker.ncbi_left_segments {
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
                                } else {
                                    seg_seq(masker, index, &seq[lend..=rend], offset + lend, segs);
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
        merge_segs_inclusive(seq.len(), &mut segs_inclusive);
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
        // ```c
        //    hilenmin = 0;               /* hilenmin - temporary default */
        // ```
        // ```c
        //       if (seg->begin - nextseg->end - 1 < hilenmin) {
        //          if (seg->end < nextseg->end) seg->end = nextseg->end;
        //          if (seg->begin > nextseg->begin) seg->begin = nextseg->begin;
        //          seg->next = nextseg->next;
        //          sfree(nextseg);
        // ```
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

    // EXPERIMENT (LOSAT_X_SEGFAST): the histogram path against the reference on
    // random protein-like sequences with planted low-complexity stretches,
    // ambiguity codes and stop codons.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2030-2032,2053-2060
    // ```c
    // s_SegSeq(SSequence* seq, SegParameters* sparamsp, SSeg **segs,
    //                    Int4 offset)
    // {
    // ...
    //    H = s_SeqEntropy(seq, window, sparamsp->maxbogus);
    //    if (H == NULL)
    //       return status;
    //    first = downset;
    //    last = seq->length - upset;
    //    lowlim = first;
    // ```
    // The test masks random sequences with both implementations and compares the intervals.
    #[test]
    fn x_fast_seg_matches_reference_on_random_sequences() {
        let mut state = 0x243F_6A88_85A3_08D3u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let valid: Vec<u8> = (0u8..28)
            .filter(|&c| c == 1 || (3..=20).contains(&c) || c == 22)
            .collect();
        let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(4000);
        let mut masked = 0usize;
        for case in 0..cases {
            let params = match case % 5 {
                0 => SegParams::new(12, 2.2, 2.5),
                1 => SegParams::new(10, 1.8, 2.1),
                2 => SegParams::with_all(12, 2.2, 2.5, 2, 50),
                3 => SegParams::with_all(15, 2.6, 2.9, 15, 100),
                _ => SegParams::with_all(8, 1.5, 1.9, 1, 10),
            };
            let masker = if case % 7 == 3 {
                SegMasker::with_params(&params).keeping_all_left_segments()
            } else {
                SegMasker::with_params(&params)
            };
            assert!(masker.x_fast_applies());
            let len = 1 + (next() % if case % 11 == 0 { 3000 } else { 400 }) as usize;
            let mut seq = Vec::with_capacity(len);
            while seq.len() < len {
                match next() % 10 {
                    // a run over a tiny alphabet (possibly a homopolymer)
                    0 | 1 => {
                        let letters = 1 + (next() % 4) as usize;
                        let pool: Vec<u8> = (0..letters)
                            .map(|_| valid[(next() % valid.len() as u64) as usize])
                            .collect();
                        let run = 1 + (next() % if case % 13 == 0 { 300 } else { 60 }) as usize;
                        for _ in 0..run {
                            seq.push(pool[(next() % pool.len() as u64) as usize]);
                        }
                    }
                    // letters outside the SEG alphabet (X, B, Z, *, gap, ...)
                    2 => {
                        let run = 1 + (next() % 6) as usize;
                        for _ in 0..run {
                            seq.push([0u8, 2, 21, 23, 24, 25, 26, 27, 200][(next() % 9) as usize]);
                        }
                    }
                    _ => {
                        let run = 1 + (next() % 40) as usize;
                        for _ in 0..run {
                            seq.push(valid[(next() % valid.len() as u64) as usize]);
                        }
                    }
                }
            }
            seq.truncate(len);
            let reference = if seq.len() < masker.window {
                Vec::new()
            } else {
                masker.mask_sequence_reference(&seq)
            };
            let fast = if seq.len() < masker.window {
                Vec::new()
            } else {
                masker.x_mask_sequence(&seq)
            };
            assert_eq!(fast, reference, "case {case}, length {len}");
            masked += reference.len();
        }
        assert!(
            masked > cases / 4,
            "only {masked} intervals in {cases} cases"
        );
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1140-1157
    // ```c
    //        return status;
    //
    // 	if (filter_options->segOptions)
    // 	{
    //         SSegOptions* seg_options = filter_options->segOptions;
    //         SegParameters* sparamsp=NULL;
    //
    //         sparamsp = SegParametersNewAa();
    //         sparamsp->overlaps = TRUE;
    //         if (seg_options->window > 0)
    //             sparamsp->window = seg_options->window;
    //         if (seg_options->locut > 0.0)
    //             sparamsp->locut = seg_options->locut;
    //         if (seg_options->hicut > 0.0)
    //             sparamsp->hicut = seg_options->hicut;
    //
    // 		status = SeqBufferSeg(sequence, length, offset, sparamsp,
    //                               seqloc_retval);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2125-2149
    // ```c
    // s_MergeSegs(SSequence* seq, SSeg* segs)
    // {
    //    SSeg* seg,* nextseg;
    //    Int4 hilenmin;              /* hilenmin yet unset */
    //
    //    hilenmin = 0;               /* hilenmin - temporary default */
    //
    //    if (segs==NULL) return;
    //
    //    if (seq->length -1 - segs->end < hilenmin)
    //        segs->end = seq->length -1;
    //
    //    seg = segs;
    //    nextseg = seg->next;
    //
    //    while (nextseg!=NULL) {
    //       if (seg->begin - nextseg->end - 1 < hilenmin) {
    //          if (seg->end < nextseg->end) seg->end = nextseg->end;
    //          if (seg->begin > nextseg->begin) seg->begin = nextseg->begin;
    //          seg->next = nextseg->next;
    //          sfree(nextseg);
    //       } else {
    //          seg = nextseg;
    //       }
    //       nextseg = seg->next;
    // ```
    // Independent SeqBufferSeg API output with this actual caller's overlaps=TRUE.
    // The former assertion combined adjacent locations, unlike this core caller.
    #[test]
    fn test_mask_sequence_matches_ncbi_core_bdv02435_1_query_defaults() {
        let fasta = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/AP027131.faa"
        ));
        let seq = fasta_sequence_by_id(fasta, "BDV02435.1");
        let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));
        let mut intervals = masker.mask_sequence(&seq);
        intervals.sort_by_key(|interval| interval.start);
        let expected: Vec<_> =
            include_str!("../../tests/unit/blastx_stage_e_seg_sequence_expected.tsv")
                .lines()
                .map(|line| {
                    let (begin, end) = line.split_once('\t').unwrap();
                    MaskedInterval::new(begin.parse().unwrap(), end.parse().unwrap())
                })
                .collect();
        assert_eq!(expected.len(), 27);
        assert_eq!(intervals, expected);
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2140-2151
    // ```c
    //    while (nextseg!=NULL) {
    //       if (seg->begin - nextseg->end - 1 < hilenmin) {
    //          if (seg->end < nextseg->end) seg->end = nextseg->end;
    //          if (seg->begin > nextseg->begin) seg->begin = nextseg->begin;
    //          seg->next = nextseg->next;
    //          sfree(nextseg);
    //       } else {
    //          seg = nextseg;
    //       }
    //       nextseg = seg->next;
    //    }
    //
    // ```
    // Expected survivors come from the verbatim pinned C function, recorded in
    // Session E/oracles/seg_overlap; the Rust helper never generates expected data.
    #[test]
    fn test_merge_segs_independent_ncbi_boundary_oracle() {
        let parse = |text: &str| -> Vec<(usize, usize)> {
            if text == "-" {
                return Vec::new();
            }
            text.split(',')
                .map(|range| {
                    let (begin, end) = range.split_once(':').unwrap();
                    (begin.parse().unwrap(), end.parse().unwrap())
                })
                .collect()
        };
        for row in include_str!("../../tests/unit/blastx_stage_e_seg_expected.tsv").lines() {
            let fields: Vec<_> = row.split('\t').collect();
            let mut intervals = parse(fields[2]);
            merge_segs_inclusive(fields[1].parse().unwrap(), &mut intervals);
            assert_eq!(intervals, parse(fields[3]), "{}", fields[0]);
        }
    }
}
