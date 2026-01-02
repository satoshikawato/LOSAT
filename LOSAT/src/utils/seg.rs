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
        1 => Some(0),              // A
        3..=20 => Some((letter - 2) as usize), // C..W => 1..18
        22 => Some(19),            // Y
        _ => None,                 // bogus: -, B, X, Z, U, *, O, J, etc.
    }
}

/// State vector: sorted amino acid counts (descending order, non-zero only)
/// Reference: blast_seg.c uses this for entropy calculation
fn compute_state_vector(counts: &[u32; 20]) -> Vec<i32> {
    let mut sv: Vec<i32> = counts.iter()
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
    let ans2 = if ans1 > -100000.0 { ln_perm(sv, window_length) } else { 0.0 };
    
    ans1 + ans2 - totseq
}

/// NCBI SEG entropy (Shannon entropy, bits) implementation.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:1592-1624 (s_Entropy)
fn entropy_ncbi_from_counts(counts: &[u32; 20]) -> f64 {
    // Build state vector (sorted counts, terminating at 0)
    let mut sv: Vec<i32> = counts.iter().copied().filter(|&c| c > 0).map(|c| c as i32).collect();
    sv.sort_by(|a, b| b.cmp(a));

    let total: i32 = sv.iter().sum();
    if total == 0 {
        return 0.0;
    }

    let mut ent = 0.0f64;
    if total == 10 {
        for &c in &sv {
            // c in 1..=10
            ent += (c as f64) * LOG_WIN10[c as usize] / NCBIMATH_LN2;
        }
    } else {
        let total_f = total as f64;
        for &c in &sv {
            ent += (c as f64) * ((c as f64) / total_f).ln() / NCBIMATH_LN2;
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
    pub fn with_all(window: usize, locut: f64, hicut: f64, maxbogus: usize, maxtrim: usize) -> Self {
        let mut params = Self::new(window, locut, hicut);
        params.maxbogus = maxbogus.min(window);
        params.maxtrim = maxtrim;
        params
    }
}

/// SEG masker implementation following NCBI BLAST's SEG algorithm
pub struct SegMasker {
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
            window: params.window,
            locut: params.locut,
            hicut: params.hicut,
            maxbogus: params.maxbogus,
            maxtrim: params.maxtrim,
            downset,
            upset,
        }
    }

    /// Calculate Shannon entropy for a window of amino acids (in bits)
    /// Reference: Original SEG algorithm (Wootton & Federhen)
    /// 
    /// Entropy range: 0 (all same amino acid) to ~4.32 (log2(20), uniform distribution)
    /// Thresholds: locut=2.2, hicut=2.5 are designed for this scale
    fn calculate_entropy(&self, window: &[u8]) -> (f64, usize) {
        if window.is_empty() {
            return (0.0, 0);
        }

        // NCBI SEG uses ncbistdaa codes but a 20-letter subset. Everything else is bogus.
        // Reference: blast_seg.c Alpha init + s_CompOn bogus counting.
        let mut counts = [0u32; 20];
        let mut bogus_count = 0usize;

        for &aa in window {
            if let Some(idx) = seg_alpha_index_ncbistdaa(aa) {
                counts[idx] += 1;
            } else {
                bogus_count += 1;
            }
        }

        // Calculate entropy in bits (NCBI s_Entropy)
        let entropy = entropy_ncbi_from_counts(&counts);
        
        (entropy, bogus_count)
    }

    /// Calculate entropy array for the entire sequence
    /// Returns array of entropy values, with -1.0 for invalid positions
    /// Reference: blast_seg.c s_SeqEntropy()
    fn calculate_entropy_array(&self, seq: &[u8]) -> Vec<f64> {
        let len = seq.len();
        let mut h = vec![-1.0; len];

        if len < self.window {
            return h;
        }

        // Calculate entropy for each valid window position
        // Valid positions: downset to (len - upset)
        // Reference: blast_seg.c:2058-2059, 1781
        let first = self.downset;
        let last = len.saturating_sub(self.upset);

        for i in first..=last {
            let window_start = i - self.downset;
            let window_end = window_start + self.window;
            if window_end <= len {
                let (entropy, bogus) = self.calculate_entropy(&seq[window_start..window_end]);
                
                // Skip windows with too many invalid amino acids (NCBI line 1789-1793)
                if bogus > self.maxbogus {
                    h[i] = -1.0;
                } else {
                    h[i] = entropy;
                }
            }
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

        // NCBI:
        //   minlen = 1;
        //   if ((seq->length-maxtrim)>minlen) minlen = seq->length-maxtrim;
        let mut minlen: usize = 1;
        if seg_len.saturating_sub(self.maxtrim) > minlen {
            minlen = seg_len - self.maxtrim;
        }

        let mut best_lend: usize = 0;
        let mut best_rend: usize = seg_len - 1;
        let mut minprob: f64 = 1.0;

        // For each candidate length, slide a window across the segment and find min prob.
        // This is a direct port of NCBI's nested loop (len desc, shift by 1).
        for cur_len in (minlen + 1..=seg_len).rev() {
            // Initialize composition for window [0..cur_len)
            let mut counts = [0u32; 20];
            for &aa in &seq[leftend..leftend + cur_len] {
                if let Some(idx) = seg_alpha_index_ncbistdaa(aa) {
                    counts[idx] += 1;
                }
            }

            // Helper to compute prob from counts (K2 probability)
            let mut prob_from_counts = |counts: &[u32; 20], total: i32| -> f64 {
                if total <= 0 {
                    return 1.0;
                }
                let sv = compute_state_vector(counts);
                get_prob(&sv, total)
            };

            let mut i = 0usize;
            loop {
                let total: i32 = counts.iter().map(|&c| c as i32).sum();
                let prob = prob_from_counts(&counts, cur_len as i32);
                if prob < minprob {
                    minprob = prob;
                    best_lend = i;
                    best_rend = i + cur_len - 1;
                }

                if i + cur_len >= seg_len {
                    break;
                }

                // Slide by 1: remove outgoing, add incoming
                let out_aa = seq[leftend + i];
                if let Some(idx) = seg_alpha_index_ncbistdaa(out_aa) {
                    counts[idx] = counts[idx].saturating_sub(1);
                }
                let in_aa = seq[leftend + i + cur_len];
                if let Some(idx) = seg_alpha_index_ncbistdaa(in_aa) {
                    counts[idx] += 1;
                }
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

        // NCBI s_SegSeq returns inclusive [begin,end] segments; we collect those and then
        // convert to [start,end) MaskedInterval at the end.
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

                    // Record segment in original coordinates
                    segs.push((leftend + offset, rightend + offset));

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

        // NCBI s_MergeSegs (hilenmin=0) semantics: merge overlapping segments only.
        if !segs_inclusive.is_empty() {
            // NCBI list order is descending begin due to prepend; mimic before merge.
            segs_inclusive.sort_by(|a, b| b.0.cmp(&a.0));
            let mut merged: Vec<(usize, usize)> = Vec::new();
            let mut cur = segs_inclusive[0];
            for &next in &segs_inclusive[1..] {
                // overlap if cur.begin <= next.end (inclusive coords)
                if cur.0 <= next.1 {
                    cur.0 = cur.0.min(next.0);
                    cur.1 = cur.1.max(next.1);
                } else {
                    merged.push(cur);
                    cur = next;
                }
            }
            merged.push(cur);
            segs_inclusive = merged;
        }

        // Convert to MaskedInterval [start, end) with end exclusive
        segs_inclusive
            .into_iter()
            .map(|(b, e)| MaskedInterval::new(b, e.saturating_add(1)))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seg_params_default() {
        let params = SegParams::default();
        assert_eq!(params.window, 12);
        assert_eq!(params.locut, 2.2);
        assert_eq!(params.hicut, 2.5);
    }

    #[test]
    fn test_seg_params_validation() {
        let params = SegParams::new(0, -1.0, 1.0);
        assert_eq!(params.window, 12); // Should default
        assert_eq!(params.locut, 2.2); // Should default
        assert_eq!(params.hicut, 2.5); // Should default to max(locut, 2.5)
    }

    #[test]
    fn test_calculate_entropy() {
        let masker = SegMasker::with_defaults();
        
        // Low complexity: all same amino acid (entropy = 0)
        let low_complex = vec![0u8; 12]; // All alanine (0)
        let (entropy_low, _) = masker.calculate_entropy(&low_complex);
        assert!(entropy_low < 1.0, "Low complexity should have low entropy: got {}", entropy_low);
        
        // High complexity: all different amino acids
        let high_complex: Vec<u8> = (0..12).map(|i| (i % 20) as u8).collect();
        let (entropy_high, _) = masker.calculate_entropy(&high_complex);
        assert!(entropy_high > 2.0, "High complexity should have high entropy: got {}", entropy_high);
    }

    #[test]
    fn test_mask_simple_repeat() {
        let masker = SegMasker::with_defaults();
        
        // Simple repeat sequence should be masked
        let seq: Vec<u8> = vec![0u8; 100]; // All alanine
        let intervals = masker.mask_sequence(&seq);
        
        assert!(!intervals.is_empty(), "Poly-alanine sequence should be masked");
        if !intervals.is_empty() {
            assert!(intervals[0].start < intervals[0].end);
        }
    }

    #[test]
    fn test_mask_short_sequence() {
        let masker = SegMasker::with_defaults();
        
        // Very short sequences should return empty
        let seq = vec![0u8; 5];
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
        assert!(total_masked < seq.len() / 2, "Complex sequence should not be heavily masked");
    }
}
