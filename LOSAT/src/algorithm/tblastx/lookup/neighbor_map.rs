//! Neighbor word map for efficient k-mer matching
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c

use crate::utils::matrix::blosum62_score;
use crate::stats::{KarlinParams, lookup_protein_params_ungapped};
use crate::config::ScoringMatrix;
use crate::stats::karlin_calc::{
    compute_aa_composition, compute_std_aa_composition, compute_score_freq_profile,
    compute_karlin_params_ungapped, apply_check_ideal,
};
use super::super::diagnostics::diagnostics_enabled;
use super::super::translation::QueryFrame;
use super::{
    QueryContext, LOOKUP_WORD_LENGTH, LOOKUP_ALPHABET_SIZE,
    ilog2, compute_backbone_size, compute_mask, encode_kmer_3, pv_test, pv_set,
};
use super::backbone::build_direct_lookup;

// Presence Vector bucket bits
const PV_BUCKET_BITS: usize = 64;

/// Pre-computed neighbor map: for each k-mer, list of k-mers that are neighbors (score >= threshold)
/// This is the key optimization: compute once, use during every scan.
pub struct NeighborMap {
    /// For each k-mer index, the list of k-mer indices that are neighbors
    /// neighbor_map[subject_kmer] = [query_kmer_1, query_kmer_2, ...] such that
    /// score(query_kmer, subject_kmer) >= threshold
    pub map: Vec<Vec<u16>>,
    pub threshold: i32,
    pub backbone_size: usize,
}

impl NeighborMap {
    /// Pre-compute all neighbor relationships for a given threshold.
    /// This is O(alphabet^6) but only done once.
    pub fn new(threshold: i32) -> Self {
        let diag_enabled = diagnostics_enabled();
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        let charsize = ilog2(alphabet_size) + 1;
        let backbone_size = compute_backbone_size(LOOKUP_WORD_LENGTH, alphabet_size, charsize);

        if diag_enabled {
            eprintln!("Pre-computing neighbor map (threshold={})...", threshold);
        }
        let start = std::time::Instant::now();

        // For each possible k-mer (subject side), find all k-mers (query side)
        // that would be neighbors
        let mut map: Vec<Vec<u16>> = vec![Vec::new(); backbone_size];

        // Row max for pruning
        let row_max: Vec<i32> = (0..alphabet_size)
            .map(|i| {
                (0..alphabet_size)
                    .map(|j| blosum62_score(i as u8, j as u8))
                    .max()
                    .unwrap_or(-4)
            })
            .collect();

        let mask = compute_mask(LOOKUP_WORD_LENGTH, charsize);
        let residue_mask = (1usize << charsize) - 1;

        // For each subject k-mer S
        for s_idx in 0..backbone_size {
            let s0 = (s_idx >> (2 * charsize)) & residue_mask;
            let s1 = (s_idx >> charsize) & residue_mask;
            let s2 = s_idx & residue_mask;

            if s0 >= alphabet_size || s1 >= alphabet_size || s2 >= alphabet_size {
                continue;
            }

            // Find all query k-mers Q such that score(Q, S) >= threshold
            let rm12 = row_max[s1] + row_max[s2];
            let rm2 = row_max[s2];

            for q0 in 0..alphabet_size {
                let sc0 = blosum62_score(q0 as u8, s0 as u8);
                if sc0 + rm12 < threshold {
                    continue;
                }
                for q1 in 0..alphabet_size {
                    let sc1 = sc0 + blosum62_score(q1 as u8, s1 as u8);
                    if sc1 + rm2 < threshold {
                        continue;
                    }
                    for q2 in 0..alphabet_size {
                        if sc1 + blosum62_score(q2 as u8, s2 as u8) >= threshold {
                            let q_idx = encode_kmer_3(q0, q1, q2, charsize) & mask;
                            map[s_idx].push(q_idx as u16);
                        }
                    }
                }
            }
        }

        let total_neighbors: usize = map.iter().map(|v| v.len()).sum();
        let nonempty = map.iter().filter(|v| !v.is_empty()).count();
        let elapsed = start.elapsed();

        if diag_enabled {
            eprintln!("Neighbor map computed in {:?}", elapsed);
            eprintln!("  Non-empty entries: {}", nonempty);
            eprintln!("  Total neighbor pairs: {}", total_neighbors);
            eprintln!(
                "  Avg neighbors per k-mer: {:.1}",
                total_neighbors as f64 / nonempty.max(1) as f64
            );
        }

        Self {
            map,
            threshold,
            backbone_size,
        }
    }

    /// Get neighbors for a subject k-mer
    #[inline]
    pub fn get_neighbors(&self, subject_kmer: usize) -> &[u16] {
        &self.map[subject_kmer]
    }

    /// Build expanded lookup table using neighbor map.
    /// This is more efficient than generating neighbors per query position.
    pub fn build_expanded_lookup(
        &self,
        query_lookup: &[Vec<(u32, u8, u32)>],
    ) -> Vec<Vec<(u32, u8, u32)>> {
        let diag_enabled = diagnostics_enabled();
        let start = std::time::Instant::now();

        // For each subject k-mer, collect all query positions that would match
        let mut expanded: Vec<Vec<(u32, u8, u32)>> = vec![Vec::new(); self.backbone_size];

        // Build reverse neighbor map: query_kmer -> [subject_kmer, ...]
        // This is essentially: for which subject k-mers would this query k-mer be a neighbor?
        // In our current map: subject_kmer -> [query_kmer, ...]
        // Reversed: query_kmer -> [subject_kmer, ...]
        let mut reverse_map: Vec<Vec<u16>> = vec![Vec::new(); self.backbone_size];
        for (s_kmer, neighbors) in self.map.iter().enumerate() {
            for &q_kmer in neighbors {
                reverse_map[q_kmer as usize].push(s_kmer as u16);
            }
        }

        // Now, for each query k-mer with positions, add those positions to all
        // subject k-mers that are its neighbors
        for (q_kmer, positions) in query_lookup.iter().enumerate() {
            if positions.is_empty() {
                continue;
            }
            // Get all subject k-mers that would match this query k-mer
            let subject_kmers = &reverse_map[q_kmer];
            for &s_kmer in subject_kmers {
                expanded[s_kmer as usize].extend_from_slice(positions);
            }
        }

        let total_entries: usize = expanded.iter().map(|v| v.len()).sum();
        let nonempty = expanded.iter().filter(|v| !v.is_empty()).count();
        let elapsed = start.elapsed();

        if diag_enabled {
            eprintln!("Expanded lookup built in {:?}", elapsed);
            eprintln!("  Total entries: {}", total_entries);
            eprintln!("  Non-empty buckets: {}", nonempty);
        }

        expanded
    }

    /// Pre-compute all neighbor relationships, but only keep those where
    /// the query k-mer actually exists in the query lookup table.
    /// This dramatically reduces the number of neighbors to check during scan.
    pub fn new_filtered(threshold: i32, query_pv: &[u64]) -> Self {
        let diag_enabled = diagnostics_enabled();
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        let charsize = ilog2(alphabet_size) + 1;
        let backbone_size = compute_backbone_size(LOOKUP_WORD_LENGTH, alphabet_size, charsize);

        if diag_enabled {
            eprintln!("Pre-computing filtered neighbor map (threshold={})...", threshold);
        }
        let start = std::time::Instant::now();

        // For each possible k-mer (subject side), find all k-mers (query side)
        // that would be neighbors AND exist in the query
        let mut map: Vec<Vec<u16>> = vec![Vec::new(); backbone_size];

        // Row max for pruning
        let row_max: Vec<i32> = (0..alphabet_size)
            .map(|i| {
                (0..alphabet_size)
                    .map(|j| blosum62_score(i as u8, j as u8))
                    .max()
                    .unwrap_or(-4)
            })
            .collect();

        let mask = compute_mask(LOOKUP_WORD_LENGTH, charsize);
        let residue_mask = (1usize << charsize) - 1;

        // For each subject k-mer S
        for s_idx in 0..backbone_size {
            let s0 = (s_idx >> (2 * charsize)) & residue_mask;
            let s1 = (s_idx >> charsize) & residue_mask;
            let s2 = s_idx & residue_mask;

            if s0 >= alphabet_size || s1 >= alphabet_size || s2 >= alphabet_size {
                continue;
            }

            // Find all query k-mers Q such that score(Q, S) >= threshold
            // AND Q exists in the query
            let rm12 = row_max[s1] + row_max[s2];
            let rm2 = row_max[s2];

            for q0 in 0..alphabet_size {
                let sc0 = blosum62_score(q0 as u8, s0 as u8);
                if sc0 + rm12 < threshold {
                    continue;
                }
                for q1 in 0..alphabet_size {
                    let sc1 = sc0 + blosum62_score(q1 as u8, s1 as u8);
                    if sc1 + rm2 < threshold {
                        continue;
                    }
                    for q2 in 0..alphabet_size {
                        if sc1 + blosum62_score(q2 as u8, s2 as u8) >= threshold {
                            let q_idx = encode_kmer_3(q0, q1, q2, charsize) & mask;
                            // Only add if query k-mer exists
                            if pv_test(query_pv, q_idx) {
                                map[s_idx].push(q_idx as u16);
                            }
                        }
                    }
                }
            }
        }

        let total_neighbors: usize = map.iter().map(|v| v.len()).sum();
        let nonempty = map.iter().filter(|v| !v.is_empty()).count();
        let elapsed = start.elapsed();

        if diag_enabled {
            eprintln!("Filtered neighbor map computed in {:?}", elapsed);
            eprintln!("  Non-empty entries: {}", nonempty);
            eprintln!("  Total neighbor pairs: {}", total_neighbors);
            eprintln!(
                "  Avg neighbors per k-mer: {:.1}",
                total_neighbors as f64 / nonempty.max(1) as f64
            );
        }

        Self {
            map,
            threshold,
            backbone_size,
        }
    }
}

/// Lookup table with pre-computed neighbor map for fast scanning
pub struct NeighborLookup {
    /// Exact query k-mer positions: query_lookup[kmer] = [(q_idx, f_idx, aa_pos), ...]
    pub query_lookup: Vec<Vec<(u32, u8, u32)>>,
    /// Presence vector for query k-mers
    pub query_pv: Vec<u64>,
    /// Pre-computed neighbor map
    pub neighbor_map: NeighborMap,
    /// Frame bases for context lookup
    pub frame_bases: Vec<i32>,
    pub contexts: Vec<QueryContext>,
}

impl NeighborLookup {
    pub fn build(
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        _karlin_params: &KarlinParams, // Unused - computed per context
    ) -> Self {
        let diag_enabled = diagnostics_enabled();
        // Build exact query lookup
        let query_lookup = build_direct_lookup(queries);

        // Build presence vector
        let pv_size = (query_lookup.len() + PV_BUCKET_BITS - 1) / PV_BUCKET_BITS;
        let mut query_pv = vec![0u64; pv_size];
        for (idx, entries) in query_lookup.iter().enumerate() {
            if !entries.is_empty() {
                pv_set(&mut query_pv, idx);
            }
        }

        // Compute ideal Karlin parameters (kbp_ideal) - used for check_ideal logic
        // Reference: NCBI blast_stat.c:2754 Blast_ScoreBlkKbpIdealCalc
        let ideal_params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);

        // Compute standard amino acid composition (for database/subject)
        // Reference: NCBI blast_stat.c:2759 Blast_ResFreqStdComp
        let std_comp = compute_std_aa_composition();

        // Build contexts and frame bases
        let mut frame_bases = Vec::new();
        let mut contexts = Vec::new();
        let mut base: i32 = 0;

        for (q_idx, frames) in queries.iter().enumerate() {
            for (f_idx, frame) in frames.iter().enumerate() {
                frame_bases.push(base);

                // Compute context-specific Karlin parameters
                // Reference: NCBI blast_stat.c:2778-2797
                // 1. Compute amino acid composition for this context
                let ctx_comp = compute_aa_composition(&frame.aa_seq, frame.aa_len);

                // 2. Compute score frequency profile
                // Use standard composition for subject (database composition)
                let score_min = -4; // BLOSUM62 minimum
                let score_max = 11; // BLOSUM62 maximum
                let sfp = compute_score_freq_profile(&ctx_comp, &std_comp, score_min, score_max);

                // 3. Compute Karlin parameters
                let computed_params = compute_karlin_params_ungapped(&sfp)
                    .unwrap_or_else(|_| {
                        // Fallback to ideal if computation fails
                        ideal_params
                    });

                // 4. Apply check_ideal logic (tblastx uses check_ideal = TRUE)
                // Reference: NCBI blast_stat.c:2796-2797
                let final_params = apply_check_ideal(computed_params, ideal_params);

                contexts.push(QueryContext {
                    q_idx: q_idx as u32,
                    f_idx: f_idx as u8,
                    frame: frame.frame,
                    aa_seq: frame.aa_seq.clone(),
                    aa_seq_nomask: frame.aa_seq_nomask.clone(),
                    aa_len: frame.aa_len,
                    orig_len: frame.orig_len,
                    frame_base: base,
                    karlin_params: final_params,
                });
                // See build_ncbi_lookup() above: share the trailing NULLB sentinel between frames.
                base += frame.aa_seq.len() as i32 - 1;
            }
        }

        // Pre-compute neighbor map and filter by query presence
        let neighbor_map = NeighborMap::new_filtered(threshold, &query_pv);

        let total_query_entries: usize = query_lookup.iter().map(|v| v.len()).sum();
        if diag_enabled {
            eprintln!("Query lookup: {} entries (exact matches only)", total_query_entries);
        }

        Self {
            query_lookup,
            query_pv,
            neighbor_map,
            frame_bases,
            contexts,
        }
    }

    /// Check if a query k-mer exists
    #[inline]
    pub fn has_query_kmer(&self, kmer: usize) -> bool {
        pv_test(&self.query_pv, kmer)
    }

    /// Get query positions for a k-mer
    #[inline]
    pub fn get_query_hits(&self, kmer: usize) -> &[(u32, u8, u32)] {
        &self.query_lookup[kmer]
    }
}
