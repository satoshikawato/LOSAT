//! Lookup table construction for TBLASTX
//!
//! This module handles building the k-mer lookup table for amino acid sequences.

use crate::utils::dust::MaskedInterval;
use super::constants::MAX_HITS_PER_KMER;
use super::translation::QueryFrame;

/// Direct lookup table type: Vec<Vec<(query_idx, frame_idx, aa_pos)>>
pub type DirectLookup = Vec<Vec<(u32, u8, u32)>>;

/// Encode a 3-amino-acid k-mer to an index (0-15624 = 25^3 - 1)
/// Uses NCBI BLAST amino acid order: ARNDCQEGHILKMFPSTWYVBJZX* (25 amino acids)
/// Returns None if the k-mer contains a stop codon or invalid amino acid
/// 
/// Reference: NCBI BLAST uses 25 amino acids for k-mer encoding
pub fn encode_aa_kmer(seq: &[u8], pos: usize) -> Option<usize> {
    if pos + 3 > seq.len() {
        return None;
    }
    let c1 = unsafe { *seq.get_unchecked(pos) } as usize;
    let c2 = unsafe { *seq.get_unchecked(pos + 1) } as usize;
    let c3 = unsafe { *seq.get_unchecked(pos + 2) } as usize;

    // Stop codon (24 = '*')を含むk-merはシードにしない
    // NCBI BLAST order: indices 0-23 are valid, 24 is stop codon
    if c1 >= 25 || c2 >= 25 || c3 >= 25 {
        return None;
    }
    // 25^3 = 15,625 possible k-mers
    Some(c1 * 625 + c2 * 25 + c3)
}

/// Check if a DNA position range overlaps with any masked interval
fn is_dna_pos_masked(intervals: &[MaskedInterval], start: usize, end: usize) -> bool {
    for interval in intervals {
        // Check if [start, end) overlaps with [interval.start, interval.end)
        if start < interval.end && end > interval.start {
            return true;
        }
    }
    false
}

/// Convert amino acid position to DNA position range for a given frame
fn aa_pos_to_dna_range(aa_pos: usize, kmer_len: usize, frame: i8, orig_len: usize) -> (usize, usize) {
    let aa_end = aa_pos + kmer_len;
    if frame > 0 {
        // Forward frames: 1, 2, 3
        let offset = (frame - 1) as usize;
        let dna_start = aa_pos * 3 + offset;
        let dna_end = aa_end * 3 + offset;
        (dna_start, dna_end.min(orig_len))
    } else {
        // Reverse frames: -1, -2, -3
        let offset = (-frame - 1) as usize;
        let dna_end = orig_len - (aa_pos * 3 + offset);
        let dna_start = orig_len - (aa_end * 3 + offset);
        (dna_start.max(0), dna_end)
    }
}

/// Build a direct lookup table for amino acid k-mers
///
/// # Arguments
/// * `queries` - Vector of query frames for each query sequence
/// * `query_masks` - Vector of masked intervals for each query sequence
///
/// # Returns
/// Direct lookup table mapping k-mer codes to (query_idx, frame_idx, aa_pos) tuples
pub fn build_direct_lookup(
    queries: &[Vec<QueryFrame>],
    query_masks: &[Vec<MaskedInterval>],
) -> DirectLookup {
    let table_size = 15625; // 25^3 (NCBI BLAST uses 25 amino acids: ARNDCQEGHILKMFPSTWYVBJZX*)
    let mut lookup = vec![Vec::new(); table_size];

    for (q_idx, frames) in queries.iter().enumerate() {
        let masks = &query_masks[q_idx];
        for (f_idx, frame) in frames.iter().enumerate() {
            let seq = &frame.aa_seq;
            if seq.len() < 3 {
                continue;
            }
            for i in 0..=(seq.len() - 3) {
                // Check if this amino acid position corresponds to a masked DNA region
                if !masks.is_empty() {
                    let (dna_start, dna_end) = aa_pos_to_dna_range(i, 3, frame.frame, frame.orig_len);
                    if is_dna_pos_masked(masks, dna_start, dna_end) {
                        continue;
                    }
                }
                if let Some(code) = encode_aa_kmer(seq, i) {
                    lookup[code].push((q_idx as u32, f_idx as u8, i as u32));
                }
            }
        }
    }

    // Mask buckets that have too many hits (prevents excessive extension)
    for bucket in &mut lookup {
        if bucket.len() > MAX_HITS_PER_KMER {
            bucket.clear();
        }
    }
    lookup
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_aa_kmer() {
        // Amino acid indices: A=0, B=1, C=2 (not ASCII characters)
        let seq = [0u8, 1u8, 2u8];
        let code = encode_aa_kmer(&seq, 0);
        assert_eq!(code, Some(0 * 676 + 1 * 26 + 2)); // 28
        
        // Test with stop codon
        let seq_with_stop = [0u8, 1u8, 25u8]; // Contains stop codon
        let code = encode_aa_kmer(&seq_with_stop, 0);
        assert_eq!(code, None);
    }
}

