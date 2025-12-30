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

/// Natural log of 20 (alphabet size for amino acids)
const LN_20: f64 = 2.995732273553991;

use std::sync::OnceLock;

/// Precomputed ln(n!) for n = 0..256 (NCBI style lookup table)
/// Initialized on first access for fast subsequent lookups
static LN_FACT: OnceLock<[f64; 256]> = OnceLock::new();

fn init_ln_fact() -> [f64; 256] {
    let mut arr = [0.0; 256];
    let mut sum = 0.0;
    for i in 0..256 {
        if i > 1 {
            sum += (i as f64).ln();
        }
        arr[i] = sum;
    }
    arr
}

/// Calculate ln(n!) using precomputed table or Stirling's approximation
/// Reference: blast_seg.c s_lnfact()
#[inline]
fn ln_factorial(n: usize) -> f64 {
    let table = LN_FACT.get_or_init(init_ln_fact);
    if n < 256 {
        table[n]
    } else {
        // Stirling's approximation for large n
        let nf = n as f64;
        (nf + 0.5) * nf.ln() - nf + 0.9189385332
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
    let mut ans = ln_factorial(window_length as usize);
    for &count in sv {
        if count > 0 {
            ans -= ln_factorial(count as usize);
        }
    }
    ans
}

/// Calculate ln(number of compositions) for a given state vector
/// This is equation 1 from Wootton & Federhen
/// Reference: blast_seg.c s_LnAss()
fn ln_ass(sv: &[i32], alphasize: i32) -> f64 {
    if sv.is_empty() {
        return ln_factorial(alphasize as usize);
    }
    
    let mut ans = ln_factorial(alphasize as usize);
    let mut total = alphasize;
    let mut class_count = 1;
    let mut prev_val = sv[0];
    
    for i in 1..sv.len() {
        if sv[i] == prev_val {
            class_count += 1;
        } else {
            total -= class_count;
            ans -= ln_factorial(class_count as usize);
            class_count = 1;
            prev_val = sv[i];
        }
    }
    
    // Handle the last class
    total -= class_count;
    ans -= ln_factorial(class_count as usize);
    
    // Account for unused letters
    if total > 0 {
        ans -= ln_factorial(total as usize);
    }
    
    ans
}

/// Calculate the K2 complexity score (probability-based)
/// This is the natural log of P_0 from equation 3 of Wootton & Federhen
/// Reference: blast_seg.c s_GetProb()
fn get_prob(sv: &[i32], window_length: i32) -> f64 {
    let alphasize = 20; // Standard amino acid alphabet
    let totseq = (window_length as f64) * LN_20;
    
    let ans1 = ln_ass(sv, alphasize);
    let ans2 = ln_perm(sv, window_length);
    
    ans1 + ans2 - totseq
}

/// Calculate Shannon entropy for a window (in bits)
/// Reference: Original SEG algorithm (Wootton & Federhen)
/// 
/// The SEG thresholds (locut=2.2, hicut=2.5) are designed for Shannon entropy
/// measured in bits (log base 2). Range is 0 (homogeneous) to ~4.32 (log2(20)).
/// 
/// Low complexity: entropy < locut (e.g., < 2.2)
/// High complexity: entropy > hicut (e.g., > 2.5)
fn calculate_shannon_entropy(counts: &[u32; 20], window_length: usize) -> f64 {
    if window_length == 0 {
        return 0.0;
    }
    
    let total = window_length as f64;
    let mut entropy = 0.0;
    
    for &count in counts {
        if count > 0 {
            let p = (count as f64) / total;
            entropy -= p * p.log2();
        }
    }
    
    entropy
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
        // Validate and clamp parameters to NCBI BLAST ranges
        let window = if window > 0 { window } else { 12 };
        let locut = if locut >= 0.0 { locut } else { 2.2 };
        let hicut = if hicut >= 0.0 && hicut >= locut { hicut } else { locut.max(2.5) };
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

        // Count amino acid frequencies (20 standard amino acids: 0-19 in NCBI encoding)
        // Also count "bogus" (invalid) amino acids: stop codon (24), X (23), etc.
        let mut counts = [0u32; 20];
        let mut bogus_count = 0usize;
        let mut valid_count = 0usize;

        for &aa in window {
            if aa < 20 {
                counts[aa as usize] += 1;
                valid_count += 1;
            } else {
                bogus_count += 1;
            }
        }

        if valid_count == 0 {
            return (-1.0, bogus_count); // Invalid window (all stop codons or invalid)
        }

        // Calculate Shannon entropy in bits
        let entropy = calculate_shannon_entropy(&counts, valid_count);
        
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
        let mut j = i;
        while j >= limit && j < h.len() {
            if h[j] == -1.0 {
                break;
            }
            if h[j] > self.hicut {
                break;
            }
            if j == 0 {
                break;
            }
            j -= 1;
        }
        j + 1
    }

    /// Find the right boundary of a low-complexity region
    /// Starting from position i, search right until entropy > hicut
    /// Reference: blast_seg.c:1836-1848 (s_FindHigh)
    fn find_high(&self, i: usize, limit: usize, h: &[f64]) -> usize {
        let mut j = i;
        while j <= limit && j < h.len() {
            if h[j] == -1.0 {
                break;
            }
            if h[j] > self.hicut {
                break;
            }
            j += 1;
        }
        if j > 0 {
            j - 1
        } else {
            0
        }
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
    
    /// Trim a segment to minimize the probability
    /// This finds the optimal boundaries within the given segment
    /// Reference: blast_seg.c:1974-2018 (s_Trim)
    /// 
    /// NCBI uses s_OpenWin/s_ShiftWin1/s_CloseWin for efficient incremental calculation.
    /// This implementation uses a simplified but still efficient approach:
    /// - Only search when segment is larger than window size
    /// - Limit iterations by stepping through positions
    fn trim_segment(&self, seq: &[u8], leftend: usize, rightend: usize) -> (usize, usize) {
        let seg_len = rightend.saturating_sub(leftend);
        if seg_len == 0 || seg_len < self.window {
            return (leftend, rightend);
        }
        
        // If segment is small enough, use original boundaries
        // Reference: NCBI's maxtrim limits how much trimming can occur
        let max_trim_total = self.maxtrim.min(seg_len / 2);
        if max_trim_total == 0 {
            return (leftend, rightend);
        }
        
        // Calculate probability for original segment
        let orig_prob = self.get_subseq_prob(&seq[leftend..rightend]);
        
        let mut best_left = leftend;
        let mut best_right = rightend;
        let mut min_prob = orig_prob;
        
        // Try trimming from left and right by small amounts
        // This is a simplified version of NCBI's recursive trim
        let step = 1.max(max_trim_total / 10); // Limit iterations
        
        for trim_left in (0..=max_trim_total).step_by(step) {
            let new_left = leftend + trim_left;
            let remaining = seg_len - trim_left;
            let max_right_trim = max_trim_total.saturating_sub(trim_left).min(remaining.saturating_sub(self.window));
            
            for trim_right in (0..=max_right_trim).step_by(step) {
                let new_right = rightend - trim_right;
                
                if new_left >= new_right || new_right > seq.len() || new_right - new_left < self.window {
                    continue;
                }
                
                let prob = self.get_subseq_prob(&seq[new_left..new_right]);
                
                if prob < min_prob {
                    min_prob = prob;
                    best_left = new_left;
                    best_right = new_right;
                }
            }
        }
        
        (best_left, best_right)
    }

    /// Mask a sequence and return the list of masked intervals
    /// Reference: blast_seg.c:2030-2116 (s_SegSeq)
    pub fn mask_sequence(&self, seq: &[u8]) -> Vec<MaskedInterval> {
        if seq.len() < self.window {
            return Vec::new();
        }

        let mut result: Vec<MaskedInterval> = Vec::new();
        let h = self.calculate_entropy_array(seq);

        let first = self.downset;
        let last = seq.len().saturating_sub(self.upset);
        let mut lowlim = first;

        // Scan for low-complexity regions
        // Reference: blast_seg.c:2062-2112
        let mut i = first;
        while i <= last {
            // Check if entropy is below locut threshold
            if h[i] != -1.0 && h[i] <= self.locut {
                // Find boundaries of low-complexity region
                let loi = self.find_low(i, lowlim, &h);
                let hii = self.find_high(i, last, &h);

                // Convert window-relative positions to sequence positions
                // Reference: blast_seg.c:2070-2071
                let mut leftend = loi.saturating_sub(self.downset);
                let mut rightend = (hii + self.upset - 1).min(seq.len());

                // Apply trimming to optimize boundaries (NCBI s_Trim)
                // Reference: blast_seg.c:2073-2079
                if leftend < rightend && rightend <= seq.len() {
                    let (trimmed_left, trimmed_right) = self.trim_segment(seq, leftend, rightend);
                    leftend = trimmed_left;
                    rightend = trimmed_right;
                }

                if leftend < rightend {
                    // Try to merge with previous interval if overlapping or adjacent
                    if let Some(last_interval) = result.last_mut() {
                        if leftend <= last_interval.end {
                            last_interval.end = last_interval.end.max(rightend);
                        } else {
                            result.push(MaskedInterval::new(leftend, rightend));
                        }
                    } else {
                        result.push(MaskedInterval::new(leftend, rightend));
                    }
                }

                // Skip to after this region
                // Reference: blast_seg.c:2110
                i = hii.min(rightend.saturating_sub(self.downset)) + 1;
                lowlim = i;
            } else {
                i += 1;
            }
        }

        result
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

