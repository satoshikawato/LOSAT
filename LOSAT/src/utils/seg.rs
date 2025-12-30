//! SEG filter implementation for masking low-complexity regions in amino acid sequences.
//!
//! This implements the SEG algorithm as described in NCBI BLAST.
//! The algorithm identifies low-complexity regions by calculating entropy-based
//! complexity scores within sliding windows.
//!
//! Reference: 
//! - NCBI BLAST blast_seg.c/h
//! - Wootton & Federhen (1993) "Statistics of local complexity in amino acid sequences and sequence databases"
//! - Wootton & Federhen (1996) Methods Enzymol. 266:554-71

use crate::utils::dust::MaskedInterval;

/// SEG filter parameters
#[derive(Debug, Clone)]
pub struct SegParams {
    /// Window size (default: 12)
    pub window: usize,
    /// Low complexity threshold (default: 2.2)
    pub locut: f64,
    /// High complexity threshold (default: 2.5)
    pub hicut: f64,
}

impl Default for SegParams {
    fn default() -> Self {
        Self {
            window: 12,
            locut: 2.2,
            hicut: 2.5,
        }
    }
}

impl SegParams {
    pub fn new(window: usize, locut: f64, hicut: f64) -> Self {
        // Validate and clamp parameters to NCBI BLAST ranges
        let window = if window > 0 { window } else { 12 };
        let locut = if locut >= 0.0 { locut } else { 2.2 };
        let hicut = if hicut >= 0.0 && hicut >= locut { hicut } else { locut.max(2.5) };
        Self { window, locut, hicut }
    }
}

/// SEG masker implementation following NCBI BLAST's SEG algorithm
pub struct SegMasker {
    window: usize,
    locut: f64,
    hicut: f64,
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
            downset,
            upset,
        }
    }

    /// Create a SEG masker with default parameters (window=12, locut=2.2, hicut=2.5)
    pub fn with_defaults() -> Self {
        Self::new(12, 2.2, 2.5)
    }

    /// Calculate entropy for a window of amino acids
    /// Uses Shannon entropy: H = -sum(p_i * log2(p_i))
    /// where p_i is the frequency of amino acid i in the window
    fn calculate_entropy(&self, window: &[u8]) -> f64 {
        if window.is_empty() {
            return 0.0;
        }

        // Count amino acid frequencies (20 standard amino acids: 0-19 in NCBI encoding)
        // Stop codon (24) and invalid (>= 25) are ignored
        let mut counts = [0u32; 20];
        let mut total = 0u32;

        for &aa in window {
            if aa < 20 {
                counts[aa as usize] += 1;
                total += 1;
            }
        }

        if total == 0 {
            return -1.0; // Invalid window (all stop codons or invalid)
        }

        // Calculate Shannon entropy
        let mut entropy = 0.0;
        for &count in &counts {
            if count > 0 {
                let p = count as f64 / total as f64;
                entropy -= p * p.log2();
            }
        }

        entropy
    }

    /// Calculate entropy array for the entire sequence
    /// Returns array of entropy values, with -1.0 for invalid positions
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
                let entropy = self.calculate_entropy(&seq[window_start..window_end]);
                h[i] = entropy;
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
                let leftend = loi.saturating_sub(self.downset);
                let rightend = (hii + self.upset - 1).min(seq.len());

                if leftend < rightend {
                    // Try to merge with previous interval if close
                    if let Some(last_interval) = result.last_mut() {
                        // Merge if intervals are close (within window size)
                        if leftend <= last_interval.end + self.window {
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
        let entropy_low = masker.calculate_entropy(&low_complex);
        assert!(entropy_low < 1.0, "Low complexity should have low entropy");
        
        // High complexity: all different amino acids
        let high_complex: Vec<u8> = (0..12).map(|i| (i % 20) as u8).collect();
        let entropy_high = masker.calculate_entropy(&high_complex);
        assert!(entropy_high > 2.0, "High complexity should have high entropy");
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

