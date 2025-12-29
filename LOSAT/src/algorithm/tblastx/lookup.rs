//! Lookup table construction for TBLASTX
//!
//! This module handles building the k-mer lookup table for amino acid sequences.
//! Uses exact-match lookup for speed (f593910 approach).

use crate::utils::dust::MaskedInterval;
use super::diagnostics::diagnostics_enabled;
use super::translation::QueryFrame;

/// Direct lookup table type: Vec<Vec<(query_idx, frame_idx, aa_pos)>>
pub type DirectLookup = Vec<Vec<(u32, u8, u32)>>;

// Use MAX_HITS_PER_KMER from constants (200, matching f593910)
use super::constants::MAX_HITS_PER_KMER;

/// Encode a 3-amino-acid k-mer to an index (0-17575 = 26^3 - 1)
/// Returns None if the k-mer contains a stop codon or invalid amino acid
#[inline]
pub fn encode_aa_kmer(seq: &[u8], pos: usize) -> Option<usize> {
    if pos + 3 > seq.len() {
        return None;
    }
    let c1 = unsafe { *seq.get_unchecked(pos) } as usize;
    let c2 = unsafe { *seq.get_unchecked(pos + 1) } as usize;
    let c3 = unsafe { *seq.get_unchecked(pos + 2) } as usize;

    // Stop codon (25) makes k-mer invalid
    if c1 >= 25 || c2 >= 25 || c3 >= 25 {
        return None;
    }
    Some(c1 * 676 + c2 * 26 + c3)
}

/// Check if a DNA position range overlaps with any masked interval
fn is_dna_pos_masked(intervals: &[MaskedInterval], start: usize, end: usize) -> bool {
    for interval in intervals {
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
        let offset = (frame - 1) as usize;
        let dna_start = aa_pos * 3 + offset;
        let dna_end = aa_end * 3 + offset;
        (dna_start, dna_end.min(orig_len))
    } else {
        let offset = (-frame - 1) as usize;
        let dna_end = orig_len - (aa_pos * 3 + offset);
        let dna_start = orig_len - (aa_end * 3 + offset);
        (dna_start.max(0), dna_end)
    }
}

/// Build a direct lookup table for amino acid k-mers (exact match only)
///
/// Fast approach: each query k-mer is stored at its exact code position.
/// No neighborhood expansion - seed extension handles similarity.
pub fn build_direct_lookup(
    queries: &[Vec<QueryFrame>],
    query_masks: &[Vec<MaskedInterval>],
) -> DirectLookup {
    let table_size = 17576; // 26^3
    let mut lookup: DirectLookup = vec![Vec::new(); table_size];

    let diag_enabled = diagnostics_enabled();
    let mut diag_positions_used = 0usize;

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
                    if diag_enabled {
                        diag_positions_used += 1;
                    }
                }
            }
        }
    }

    // Mask buckets that have too many hits (low-complexity filter)
    let mut buckets_cleared = 0usize;
    for bucket in &mut lookup {
        if bucket.len() > MAX_HITS_PER_KMER {
            bucket.clear();
            buckets_cleared += 1;
        }
    }

    if diag_enabled {
        let buckets_nonempty = lookup.iter().filter(|b| !b.is_empty()).count();
        let total_entries: usize = lookup.iter().map(|b| b.len()).sum();
        eprintln!("\n=== TBLASTX Lookup Table ===");
        eprintln!("Query positions indexed: {}", diag_positions_used);
        eprintln!("Buckets cleared (>{} hits): {}", MAX_HITS_PER_KMER, buckets_cleared);
        eprintln!("Non-empty buckets: {} / {}", buckets_nonempty, table_size);
        eprintln!("Total entries: {}", total_entries);
        eprintln!("============================\n");
    }

    lookup
}

/// Build lookup table (compatibility wrapper)
pub fn build_direct_lookup_with_threshold(
    queries: &[Vec<QueryFrame>],
    query_masks: &[Vec<MaskedInterval>],
    _threshold: i32,
) -> DirectLookup {
    // Threshold is ignored - using exact match for speed
    build_direct_lookup(queries, query_masks)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_aa_kmer() {
        let seq = [0u8, 1u8, 2u8];
        let code = encode_aa_kmer(&seq, 0);
        assert_eq!(code, Some(0 * 676 + 1 * 26 + 2));
        
        let seq_with_stop = [0u8, 1u8, 25u8];
        let code = encode_aa_kmer(&seq_with_stop, 0);
        assert_eq!(code, None);
    }
}
