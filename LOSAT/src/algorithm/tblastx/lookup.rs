//! Lookup table construction for TBLASTX
//!
//! This module handles building the k-mer lookup table for amino acid sequences.
//! Implements NCBI BLAST's "neighborhood word" approach for finding similar seeds.

use crate::utils::dust::MaskedInterval;
use crate::utils::matrix::MATRIX;
use super::constants::MAX_HITS_PER_KMER;
use super::translation::QueryFrame;

/// Direct lookup table type: Vec<Vec<(query_idx, frame_idx, aa_pos)>>
pub type DirectLookup = Vec<Vec<(u32, u8, u32)>>;

/// Encode a 3-amino-acid k-mer to an index (0-17575 = 26^3 - 1)
/// Returns None if the k-mer contains a stop codon or invalid amino acid
pub fn encode_aa_kmer(seq: &[u8], pos: usize) -> Option<usize> {
    if pos + 3 > seq.len() {
        return None;
    }
    let c1 = unsafe { *seq.get_unchecked(pos) } as usize;
    let c2 = unsafe { *seq.get_unchecked(pos + 1) } as usize;
    let c3 = unsafe { *seq.get_unchecked(pos + 2) } as usize;

    // 終止コドン(25)を含むk-merはシードにしない
    if c1 >= 25 || c2 >= 25 || c3 >= 25 {
        return None;
    }
    Some(c1 * 676 + c2 * 26 + c3)
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

/// Generate all neighboring k-mer codes that score >= threshold against a query k-mer
/// 
/// NCBI BLAST's neighborhood word approach: for each query word, find ALL words
/// in the alphabet that score >= threshold when aligned against the query word.
/// This allows finding similar (not just exact) matches.
///
/// # Arguments
/// * `query_kmer` - The query k-mer (3 amino acid indices)
/// * `threshold` - Minimum score threshold (NCBI BLAST BLOSUM62 default = 11)
///
/// # Returns
/// Vector of (neighbor_code, score) pairs for all neighbors scoring >= threshold
fn generate_neighborhood_words(query_kmer: [u8; 3], threshold: i32) -> Vec<(usize, i32)> {
    let mut neighbors = Vec::new();
    
    // Iterate through all possible 3-mers (26^3 = 17576 combinations)
    // For efficiency, we prune early based on partial scores
    for aa1 in 0u8..25 {
        // Matrix indexing: MATRIX[a * 27 + b]
        let score1 = MATRIX[(query_kmer[0] as usize) * 27 + (aa1 as usize)] as i32;
        // Early termination: if first position alone can't reach threshold
        // with best possible scores for positions 2 and 3 (max ~11 each)
        if score1 + 22 < threshold {
            continue;
        }
        
        for aa2 in 0u8..25 {
            let score2 = MATRIX[(query_kmer[1] as usize) * 27 + (aa2 as usize)] as i32;
            let partial_score = score1 + score2;
            // Early termination: if first two positions can't reach threshold
            // with best possible score for position 3
            if partial_score + 11 < threshold {
                continue;
            }
            
            for aa3 in 0u8..25 {
                let score3 = MATRIX[(query_kmer[2] as usize) * 27 + (aa3 as usize)] as i32;
                let total_score = partial_score + score3;
                
                if total_score >= threshold {
                    let code = (aa1 as usize) * 676 + (aa2 as usize) * 26 + (aa3 as usize);
                    neighbors.push((code, total_score));
                }
            }
        }
    }
    
    neighbors
}

/// Build a direct lookup table for amino acid k-mers with neighborhood words
///
/// Implements NCBI BLAST's neighborhood word approach: for each query k-mer,
/// add ALL neighboring words that score >= threshold to the lookup table.
/// This significantly increases sensitivity compared to exact-match only.
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
    // NCBI BLAST uses threshold 11 for BLOSUM62 (reference: blast_options.c kB62_threshold)
    build_direct_lookup_with_threshold(queries, query_masks, 11)
}

/// Build a direct lookup table with configurable threshold
pub fn build_direct_lookup_with_threshold(
    queries: &[Vec<QueryFrame>],
    query_masks: &[Vec<MaskedInterval>],
    threshold: i32,
) -> DirectLookup {
    let table_size = 17576; // 26^3
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
                
                // Get the query k-mer
                let c1 = unsafe { *seq.get_unchecked(i) };
                let c2 = unsafe { *seq.get_unchecked(i + 1) };
                let c3 = unsafe { *seq.get_unchecked(i + 2) };
                
                // Skip if contains stop codon
                if c1 >= 25 || c2 >= 25 || c3 >= 25 {
                    continue;
                }
                
                let query_kmer = [c1, c2, c3];
                
                // Generate all neighboring words that score >= threshold
                // and add them to the lookup table
                let neighbors = generate_neighborhood_words(query_kmer, threshold);
                for (code, _score) in neighbors {
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

