//! Lookup table construction for TBLASTX
//!
//! This module handles building the k-mer lookup table for amino acid sequences.
//! Uses NCBI BLAST compatible neighborhood word expansion.
//! Amino acids are encoded in NCBI matrix order: ARNDCQEGHILKMFPSTWYVBJZX* (0-24)
//!
//! IMPORTANT: QueryFrame.aa_seq contains sentinel bytes at positions 0 and len-1.
//! When scanning for k-mers, we skip sentinel positions and report logical
//! amino acid positions (0-indexed from the first actual amino acid).

use crate::utils::dust::MaskedInterval;
use crate::utils::matrix::MATRIX;
use super::constants::STOP_CODON;
use super::diagnostics::diagnostics_enabled;
use super::translation::QueryFrame;

/// Direct lookup table type: Vec<Vec<(query_idx, frame_idx, aa_pos)>>
/// aa_pos is the LOGICAL position (0-indexed, not counting sentinels)
pub type DirectLookup = Vec<Vec<(u32, u8, u32)>>;

// Use MAX_HITS_PER_KMER from constants
use super::constants::MAX_HITS_PER_KMER;

/// Table size for 3-mer amino acid k-mers excluding stop codon
/// NCBI uses 24 amino acids (0-23) for k-mers, excluding stop codon (24)
/// 24^3 = 13824 possible k-mers
const TABLE_SIZE: usize = 24 * 24 * 24;

/// Encode a 3-amino-acid k-mer to an index (0-13823 = 24^3 - 1)
/// Returns None if the k-mer contains a stop codon (24), sentinel (255), or invalid amino acid
/// Amino acids must be in NCBI matrix order (0-24)
#[inline]
pub fn encode_aa_kmer(seq: &[u8], pos: usize) -> Option<usize> {
    if pos + 3 > seq.len() {
        return None;
    }
    let c1 = unsafe { *seq.get_unchecked(pos) };
    let c2 = unsafe { *seq.get_unchecked(pos + 1) };
    let c3 = unsafe { *seq.get_unchecked(pos + 2) };

    // Stop codon (24) or sentinel (255) makes k-mer invalid
    if c1 >= STOP_CODON || c2 >= STOP_CODON || c3 >= STOP_CODON {
        return None;
    }
    Some((c1 as usize) * 576 + (c2 as usize) * 24 + (c3 as usize))
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

/// Convert LOGICAL amino acid position to DNA position range for a given frame.
/// aa_pos is 0-indexed (not counting sentinels).
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
///
/// NOTE: aa_seq contains sentinels at positions 0 and len-1.
/// We scan from position 1 to len-4 (inclusive) and report logical positions (raw_pos - 1).
pub fn build_direct_lookup(
    queries: &[Vec<QueryFrame>],
    query_masks: &[Vec<MaskedInterval>],
) -> DirectLookup {
    let mut lookup: DirectLookup = vec![Vec::new(); TABLE_SIZE];

    let diag_enabled = diagnostics_enabled();
    let mut diag_positions_used = 0usize;

    for (q_idx, frames) in queries.iter().enumerate() {
        let masks = &query_masks[q_idx];
        for (f_idx, frame) in frames.iter().enumerate() {
            let seq = &frame.aa_seq;
            // aa_seq layout: [SENTINEL, aa0, aa1, ..., aaN-1, SENTINEL]
            // Valid k-mer positions in raw array: 1 to len-4 (so k-mer doesn't touch end sentinel)
            // For a sequence of N actual amino acids, aa_seq.len() = N + 2
            // Valid raw positions: 1, 2, ..., N-2 (i.e., 1 to len-4)
            if seq.len() < 5 {
                // Need at least: SENTINEL + 3 AAs + SENTINEL = 5
                continue;
            }
            // Scan from raw position 1 to len-4 (inclusive)
            for raw_pos in 1..=(seq.len() - 4) {
                // Logical AA position (0-indexed, not counting sentinels)
                let logical_pos = raw_pos - 1;
                
                // Check if this amino acid position corresponds to a masked DNA region
                if !masks.is_empty() {
                    let (dna_start, dna_end) = aa_pos_to_dna_range(logical_pos, 3, frame.frame, frame.orig_len);
                    if is_dna_pos_masked(masks, dna_start, dna_end) {
                        continue;
                    }
                }
                if let Some(code) = encode_aa_kmer(seq, raw_pos) {
                    // Store logical position (not raw position)
                    lookup[code].push((q_idx as u32, f_idx as u8, logical_pos as u32));
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
        eprintln!("Non-empty buckets: {} / {}", buckets_nonempty, TABLE_SIZE);
        eprintln!("Total entries: {}", total_entries);
        eprintln!("============================\n");
    }

    lookup
}

/// Calculate score between two 3-mers using BLOSUM62 matrix
/// Amino acids must be in NCBI matrix order (0-23, excluding stop codon)
#[inline]
fn kmer_score(kmer1: &[u8; 3], kmer2: &[u8; 3]) -> i32 {
    let mut score = 0i32;
    for i in 0..3 {
        let idx = (kmer1[i] as usize) * 25 + (kmer2[i] as usize);
        score += MATRIX[idx] as i32;
    }
    score
}

/// Decode a k-mer index back to amino acid codes (NCBI matrix order 0-23)
#[inline]
fn decode_kmer(idx: usize) -> [u8; 3] {
    let c1 = (idx / 576) as u8;
    let c2 = ((idx % 576) / 24) as u8;
    let c3 = (idx % 24) as u8;
    [c1, c2, c3]
}

/// Build lookup table with NCBI BLAST compatible neighborhood word expansion
///
/// For each query k-mer, we add entries for ALL subject k-mers that would
/// score >= threshold when aligned to that query k-mer.
/// This is the reverse of what NeighborhoodIndex does - we need to find
/// which subject k-mers hit each query k-mer.
///
/// NOTE: aa_seq contains sentinels at positions 0 and len-1.
/// We scan from position 1 to len-4 (inclusive) and report logical positions (raw_pos - 1).
pub fn build_direct_lookup_with_threshold(
    queries: &[Vec<QueryFrame>],
    query_masks: &[Vec<MaskedInterval>],
    threshold: i32,
) -> DirectLookup {
    let mut lookup: DirectLookup = vec![Vec::new(); TABLE_SIZE];

    let diag_enabled = diagnostics_enabled();
    let mut diag_positions_used = 0usize;
    let mut diag_neighbor_entries = 0usize;

    for (q_idx, frames) in queries.iter().enumerate() {
        let masks = &query_masks[q_idx];
        for (f_idx, frame) in frames.iter().enumerate() {
            let seq = &frame.aa_seq;
            // aa_seq layout: [SENTINEL, aa0, aa1, ..., aaN-1, SENTINEL]
            if seq.len() < 5 {
                continue;
            }
            // Scan from raw position 1 to len-4 (inclusive)
            for raw_pos in 1..=(seq.len() - 4) {
                let logical_pos = raw_pos - 1;
                
                // Check if this amino acid position corresponds to a masked DNA region
                if !masks.is_empty() {
                    let (dna_start, dna_end) = aa_pos_to_dna_range(logical_pos, 3, frame.frame, frame.orig_len);
                    if is_dna_pos_masked(masks, dna_start, dna_end) {
                        continue;
                    }
                }
                
                // Get the query k-mer (from raw position in aa_seq)
                if let Some(_code) = encode_aa_kmer(seq, raw_pos) {
                    let query_kmer = [
                        unsafe { *seq.get_unchecked(raw_pos) },
                        unsafe { *seq.get_unchecked(raw_pos + 1) },
                        unsafe { *seq.get_unchecked(raw_pos + 2) },
                    ];
                    
                    if diag_enabled {
                        diag_positions_used += 1;
                    }
                    
                    // For each possible subject k-mer (0 to 24^3-1), check if it scores >= threshold
                    // against this query k-mer, and if so, add this query position
                    // to that subject k-mer's bucket
                    for subject_code in 0..TABLE_SIZE {
                        let subject_kmer = decode_kmer(subject_code);
                        
                        let score = kmer_score(&query_kmer, &subject_kmer);
                        if score >= threshold {
                            // Store logical position
                            lookup[subject_code].push((q_idx as u32, f_idx as u8, logical_pos as u32));
                            if diag_enabled {
                                diag_neighbor_entries += 1;
                            }
                        }
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
        eprintln!("\n=== TBLASTX Lookup Table (Neighborhood) ===");
        eprintln!("Query positions indexed: {}", diag_positions_used);
        eprintln!("Neighbor entries added: {}", diag_neighbor_entries);
        eprintln!("Avg neighbors per position: {:.1}", diag_neighbor_entries as f64 / diag_positions_used.max(1) as f64);
        eprintln!("Buckets cleared (>{} hits): {}", MAX_HITS_PER_KMER, buckets_cleared);
        eprintln!("Non-empty buckets: {} / {}", buckets_nonempty, TABLE_SIZE);
        eprintln!("Total entries: {}", total_entries);
        eprintln!("===========================================\n");
    }

    lookup
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_aa_kmer() {
        // Test with valid k-mer (A=0, R=1, N=2 in NCBI order)
        let seq = [0u8, 1u8, 2u8];
        let code = encode_aa_kmer(&seq, 0);
        assert_eq!(code, Some(0 * 576 + 1 * 24 + 2));
        
        // Test with stop codon (24) - should return None
        let seq_with_stop = [0u8, 1u8, STOP_CODON];
        let code = encode_aa_kmer(&seq_with_stop, 0);
        assert_eq!(code, None);
        
        // Test with X (23) - should be valid (X is a valid amino acid in NCBI)
        let seq_with_x = [0u8, 1u8, 23u8];
        let code = encode_aa_kmer(&seq_with_x, 0);
        assert_eq!(code, Some(0 * 576 + 1 * 24 + 23));
        
        // Test with sentinel (255) - should return None
        let seq_with_sentinel = [0u8, 255u8, 2u8];
        let code = encode_aa_kmer(&seq_with_sentinel, 0);
        assert_eq!(code, None);
    }
    
    #[test]
    fn test_decode_kmer() {
        // Test round-trip encoding/decoding
        for c1 in 0..24u8 {
            for c2 in 0..24u8 {
                for c3 in 0..24u8 {
                    let seq = [c1, c2, c3];
                    if let Some(code) = encode_aa_kmer(&seq, 0) {
                        let decoded = decode_kmer(code);
                        assert_eq!(decoded, seq, "Failed for {:?}", seq);
                    }
                }
            }
        }
    }
}
