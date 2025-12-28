//! Lookup table construction for TBLASTX
//!
//! This module handles building the k-mer lookup table for amino acid sequences.
//! Implements NCBI BLAST's "neighborhood word" approach for finding similar seeds.

use crate::utils::dust::MaskedInterval;
use crate::utils::matrix::MATRIX;
use super::constants::MAX_HITS_PER_KMER;
use super::diagnostics::diagnostics_enabled;
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
    // NCBI BLAST+ tblastx default: Neighboring words threshold = 13
    // (visible in pairwise output header: "Neighboring words threshold: 13")
    build_direct_lookup_with_threshold(queries, query_masks, 13)
}

/// Build a direct lookup table with configurable threshold
pub fn build_direct_lookup_with_threshold(
    queries: &[Vec<QueryFrame>],
    query_masks: &[Vec<MaskedInterval>],
    threshold: i32,
) -> DirectLookup {
    let table_size = 17576; // 26^3
    let mut lookup = vec![Vec::new(); table_size];

    let diag_enabled = diagnostics_enabled();
    let mut diag_query_positions_total: usize = 0;
    let mut diag_query_positions_masked: usize = 0;
    let mut diag_query_positions_stop: usize = 0;
    let mut diag_query_positions_used: usize = 0;
    let mut diag_unique_query_kmers: usize = 0;
    let mut diag_neighbors_generated: usize = 0;
    let mut diag_neighbors_max_per_pos: usize = 0;
    let mut query_kmers_seen: Vec<bool> = if diag_enabled {
        vec![false; table_size]
    } else {
        Vec::new()
    };

    for (q_idx, frames) in queries.iter().enumerate() {
        let masks = &query_masks[q_idx];
        for (f_idx, frame) in frames.iter().enumerate() {
            let seq = &frame.aa_seq;
            if seq.len() < 3 {
                continue;
            }
            for i in 0..=(seq.len() - 3) {
                if diag_enabled {
                    diag_query_positions_total += 1;
                }
                // Check if this amino acid position corresponds to a masked DNA region
                if !masks.is_empty() {
                    let (dna_start, dna_end) = aa_pos_to_dna_range(i, 3, frame.frame, frame.orig_len);
                    if is_dna_pos_masked(masks, dna_start, dna_end) {
                        if diag_enabled {
                            diag_query_positions_masked += 1;
                        }
                        continue;
                    }
                }
                
                // Get the query k-mer
                let c1 = unsafe { *seq.get_unchecked(i) };
                let c2 = unsafe { *seq.get_unchecked(i + 1) };
                let c3 = unsafe { *seq.get_unchecked(i + 2) };
                
                // Skip if contains stop codon
                if c1 >= 25 || c2 >= 25 || c3 >= 25 {
                    if diag_enabled {
                        diag_query_positions_stop += 1;
                    }
                    continue;
                }
                
                let query_kmer = [c1, c2, c3];

                if diag_enabled {
                    diag_query_positions_used += 1;
                    let q_code = (c1 as usize) * 676 + (c2 as usize) * 26 + (c3 as usize);
                    if !query_kmers_seen[q_code] {
                        query_kmers_seen[q_code] = true;
                        diag_unique_query_kmers += 1;
                    }
                }
                
                // Generate all neighboring words that score >= threshold
                // and add them to the lookup table
                let neighbors = generate_neighborhood_words(query_kmer, threshold);
                if diag_enabled {
                    diag_neighbors_generated += neighbors.len();
                    diag_neighbors_max_per_pos = diag_neighbors_max_per_pos.max(neighbors.len());
                }
                for (code, _score) in neighbors {
                    lookup[code].push((q_idx as u32, f_idx as u8, i as u32));
                }
            }
        }
    }

    // Compute and print lookup statistics before masking (diagnostics only)
    let (mut buckets_nonempty_pre, mut total_entries_pre, mut max_bucket_pre) = (0usize, 0usize, 0usize);
    let mut top_buckets_pre: Vec<(usize, usize)> = Vec::new(); // (len, code)
    if diag_enabled {
        top_buckets_pre.reserve(table_size);
        for (code, bucket) in lookup.iter().enumerate() {
            let len = bucket.len();
            if len > 0 {
                buckets_nonempty_pre += 1;
                total_entries_pre += len;
                if len > max_bucket_pre {
                    max_bucket_pre = len;
                }
            }
            top_buckets_pre.push((len, code));
        }
        top_buckets_pre.sort_by(|a, b| b.0.cmp(&a.0));
    }

    // Mask buckets that have too many hits (prevents excessive extension)
    let mut buckets_cleared: usize = 0;
    let mut entries_cleared: usize = 0;
    let mut max_cleared_bucket: usize = 0;
    for bucket in &mut lookup {
        let len = bucket.len();
        if len > MAX_HITS_PER_KMER {
            if diag_enabled {
                buckets_cleared += 1;
                entries_cleared += len;
                max_cleared_bucket = max_cleared_bucket.max(len);
            }
            bucket.clear();
        }
    }

    // Compute and print lookup statistics after masking (diagnostics only)
    if diag_enabled {
        let (mut buckets_nonempty_post, mut total_entries_post, mut max_bucket_post) = (0usize, 0usize, 0usize);
        for bucket in &lookup {
            let len = bucket.len();
            if len > 0 {
                buckets_nonempty_post += 1;
                total_entries_post += len;
                if len > max_bucket_post {
                    max_bucket_post = len;
                }
            }
        }

        eprintln!("\n=== TBLASTX Lookup Table Diagnostics ===");
        eprintln!("Lookup construction:");
        eprintln!("  Neighborhood threshold:     {}", threshold);
        eprintln!("  Table size (26^3):          {}", table_size);
        eprintln!("  Query positions scanned:    {}", diag_query_positions_total);
        eprintln!("  Query positions masked:     {}", diag_query_positions_masked);
        eprintln!("  Query positions stop-codon: {}", diag_query_positions_stop);
        eprintln!("  Query positions used:       {}", diag_query_positions_used);
        eprintln!("  Unique query k-mers:        {}", diag_unique_query_kmers);
        if diag_query_positions_used > 0 {
            eprintln!(
                "  Neighbors generated:        {} (avg {:.1}, max {})",
                diag_neighbors_generated,
                diag_neighbors_generated as f64 / diag_query_positions_used as f64,
                diag_neighbors_max_per_pos
            );
        } else {
            eprintln!("  Neighbors generated:        {}", diag_neighbors_generated);
        }

        eprintln!("Bucket statistics (pre-mask):");
        eprintln!(
            "  Non-empty buckets:          {} / {}",
            buckets_nonempty_pre, table_size
        );
        eprintln!("  Total entries:              {}", total_entries_pre);
        eprintln!("  Max bucket size:            {}", max_bucket_pre);

        eprintln!("Frequency masking (MAX_HITS_PER_KMER = {}):", MAX_HITS_PER_KMER);
        eprintln!("  Buckets cleared:            {}", buckets_cleared);
        eprintln!("  Entries removed:            {}", entries_cleared);
        eprintln!("  Max cleared bucket size:    {}", max_cleared_bucket);

        eprintln!("Bucket statistics (post-mask):");
        eprintln!(
            "  Non-empty buckets:          {} / {}",
            buckets_nonempty_post, table_size
        );
        eprintln!("  Total entries:              {}", total_entries_post);
        eprintln!("  Max bucket size:            {}", max_bucket_post);

        // Print a small snapshot of the largest buckets to spot low-complexity dominance
        eprintln!("Top buckets by size (pre-mask, top 10):");
        for (len, code) in top_buckets_pre.iter().filter(|(l, _)| *l > 0).take(10) {
            let aa1 = code / 676;
            let aa2 = (code / 26) % 26;
            let aa3 = code % 26;
            eprintln!(
                "  code {:5} [{:02},{:02},{:02}]  len {}",
                code, aa1, aa2, aa3, len
            );
        }
        eprintln!("========================================\n");
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

