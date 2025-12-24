use crate::utils::matrix::MATRIX;

/// Size of the amino acid k-mer table (26^3 = 17,576 for 3-mers)
pub const AA_TABLE_SIZE: usize = 26 * 26 * 26;

/// Maximum hits per k-mer before filtering
pub const MAX_HITS_PER_KMER: usize = 200;

/// Stop codon index
pub const STOP_CODON: u8 = 25;

/// Type alias for amino acid k-mer lookup table
/// Direct-address table: index = encoded k-mer, value = list of (query_idx, frame, position)
pub type AaKmerLookup = Vec<Vec<(u32, usize, usize)>>;

/// Encode a 3-amino acid k-mer to an index
///
/// Returns None if the k-mer contains stop codons or invalid characters
#[inline]
pub fn encode_aa_kmer(seq: &[u8], pos: usize) -> Option<usize> {
    if pos + 3 > seq.len() {
        return None;
    }

    let c1 = seq[pos];
    let c2 = seq[pos + 1];
    let c3 = seq[pos + 2];

    // Check for stop codons
    if c1 == STOP_CODON || c2 == STOP_CODON || c3 == STOP_CODON {
        return None;
    }

    // Check bounds (amino acids are 0-25)
    if c1 >= 26 || c2 >= 26 || c3 >= 26 {
        return None;
    }

    Some((c1 as usize) * 676 + (c2 as usize) * 26 + (c3 as usize))
}

/// Decode an index back to a 3-amino acid k-mer
pub fn decode_aa_kmer(index: usize) -> [u8; 3] {
    let c1 = (index / 676) as u8;
    let c2 = ((index / 26) % 26) as u8;
    let c3 = (index % 26) as u8;
    [c1, c2, c3]
}

/// Build a direct-address lookup table for amino acid k-mers
pub fn build_aa_lookup(
    frames: &[(Vec<u8>, usize)], // (aa_sequence, frame_number)
    query_idx: u32,
) -> AaKmerLookup {
    let mut lookup: AaKmerLookup = vec![Vec::new(); AA_TABLE_SIZE];

    for (frame_idx, (aa_seq, _frame_num)) in frames.iter().enumerate() {
        if aa_seq.len() < 3 {
            continue;
        }

        for pos in 0..=(aa_seq.len() - 3) {
            if let Some(kmer_idx) = encode_aa_kmer(aa_seq, pos) {
                let hits = &mut lookup[kmer_idx];
                if hits.len() < MAX_HITS_PER_KMER {
                    hits.push((query_idx, frame_idx, pos));
                }
            }
        }
    }

    lookup
}

/// Configuration for amino acid word finding
#[derive(Debug, Clone, Copy)]
pub struct AaWordFinderConfig {
    /// Word size (k-mer length, typically 3)
    pub word_size: usize,
    /// Threshold for neighborhood word finding
    pub threshold: i32,
    /// Whether to use neighborhood word finding (NCBI style)
    pub use_neighborhood: bool,
    /// Maximum hits per k-mer
    pub max_hits_per_kmer: usize,
}

impl Default for AaWordFinderConfig {
    fn default() -> Self {
        Self {
            word_size: 3,
            threshold: 13,
            use_neighborhood: false,
            max_hits_per_kmer: MAX_HITS_PER_KMER,
        }
    }
}

impl AaWordFinderConfig {
    /// Configuration for NCBI BLAST compatible seeding
    pub fn ncbi_compat(threshold: i32) -> Self {
        Self {
            word_size: 3,
            threshold,
            use_neighborhood: true,
            max_hits_per_kmer: MAX_HITS_PER_KMER,
        }
    }
}

/// Neighborhood word index for finding similar k-mers
///
/// Pre-computes all k-mers that score above threshold for each possible k-mer
pub struct NeighborhoodIndex {
    /// For each k-mer index, list of neighbor k-mer indices that score >= threshold
    neighbors: Vec<Vec<u16>>,
    /// Threshold used for neighbor computation
    threshold: i32,
}

impl NeighborhoodIndex {
    /// Create a new neighborhood index with the given threshold
    ///
    /// This pre-computes neighbors for all 17,576 possible 3-mers
    pub fn new(threshold: i32) -> Self {
        let mut neighbors = vec![Vec::new(); AA_TABLE_SIZE];

        // For each possible k-mer, find all neighbors that score >= threshold
        for kmer_idx in 0..AA_TABLE_SIZE {
            let kmer = decode_aa_kmer(kmer_idx);

            // Skip k-mers with stop codons
            if kmer[0] == STOP_CODON || kmer[1] == STOP_CODON || kmer[2] == STOP_CODON {
                continue;
            }

            // Check all possible neighbor k-mers
            for neighbor_idx in 0..AA_TABLE_SIZE {
                let neighbor = decode_aa_kmer(neighbor_idx);

                // Skip neighbors with stop codons
                if neighbor[0] == STOP_CODON
                    || neighbor[1] == STOP_CODON
                    || neighbor[2] == STOP_CODON
                {
                    continue;
                }

                // Calculate score using substitution matrix
                let score = get_kmer_score(&kmer, &neighbor);

                if score >= threshold {
                    neighbors[kmer_idx].push(neighbor_idx as u16);
                }
            }
        }

        Self {
            neighbors,
            threshold,
        }
    }

    /// Get all neighbor k-mers for a given k-mer index
    #[inline]
    pub fn get_neighbors(&self, kmer_idx: usize) -> &[u16] {
        if kmer_idx < AA_TABLE_SIZE {
            &self.neighbors[kmer_idx]
        } else {
            &[]
        }
    }

    /// Get the threshold used for this index
    pub fn threshold(&self) -> i32 {
        self.threshold
    }

    /// Get statistics about the neighborhood index
    pub fn stats(&self) -> NeighborhoodStats {
        let mut total_neighbors = 0;
        let mut max_neighbors = 0;
        let mut non_empty = 0;

        for neighbors in &self.neighbors {
            let count = neighbors.len();
            total_neighbors += count;
            if count > max_neighbors {
                max_neighbors = count;
            }
            if count > 0 {
                non_empty += 1;
            }
        }

        NeighborhoodStats {
            total_kmers: AA_TABLE_SIZE,
            non_empty_kmers: non_empty,
            total_neighbors,
            max_neighbors,
            avg_neighbors: if non_empty > 0 {
                total_neighbors as f64 / non_empty as f64
            } else {
                0.0
            },
        }
    }
}

/// Statistics about a neighborhood index
#[derive(Debug, Clone)]
pub struct NeighborhoodStats {
    pub total_kmers: usize,
    pub non_empty_kmers: usize,
    pub total_neighbors: usize,
    pub max_neighbors: usize,
    pub avg_neighbors: f64,
}

/// Calculate the score of aligning two 3-mers using the substitution matrix
#[inline]
fn get_kmer_score(kmer1: &[u8; 3], kmer2: &[u8; 3]) -> i32 {
    let mut score = 0i32;
    for i in 0..3 {
        let a = kmer1[i] as usize;
        let b = kmer2[i] as usize;
        if a < 27 && b < 27 {
            score += MATRIX[a * 27 + b] as i32;
        }
    }
    score
}

/// Get score for two amino acids from the substitution matrix
#[inline]
pub fn get_aa_score(a: u8, b: u8) -> i8 {
    let a_idx = a as usize;
    let b_idx = b as usize;
    if a_idx < 27 && b_idx < 27 {
        MATRIX[a_idx * 27 + b_idx]
    } else {
        -4 // Default penalty for unknown
    }
}

/// Amino acid word finder for seed detection
pub struct AaWordFinder {
    /// Direct-address lookup table
    lookup: AaKmerLookup,
    /// Optional neighborhood index for NCBI-style seeding
    neighborhood: Option<NeighborhoodIndex>,
    /// Configuration
    config: AaWordFinderConfig,
}

impl AaWordFinder {
    /// Create a new word finder with exact matching only
    pub fn new(config: AaWordFinderConfig) -> Self {
        Self {
            lookup: vec![Vec::new(); AA_TABLE_SIZE],
            neighborhood: if config.use_neighborhood {
                Some(NeighborhoodIndex::new(config.threshold))
            } else {
                None
            },
            config,
        }
    }

    /// Add query frames to the index
    pub fn add_query(&mut self, query_idx: u32, frames: &[(Vec<u8>, usize)]) {
        for (frame_idx, (aa_seq, _frame_num)) in frames.iter().enumerate() {
            if aa_seq.len() < 3 {
                continue;
            }

            for pos in 0..=(aa_seq.len() - 3) {
                if let Some(kmer_idx) = encode_aa_kmer(aa_seq, pos) {
                    let hits = &mut self.lookup[kmer_idx];
                    if hits.len() < self.config.max_hits_per_kmer {
                        hits.push((query_idx, frame_idx, pos));
                    }
                }
            }
        }
    }

    /// Find all seed hits for a subject frame
    ///
    /// If neighborhood finding is enabled, also returns hits for similar k-mers
    pub fn find_seeds(&self, subject_aa: &[u8], subject_frame: usize) -> Vec<AaSeedHit> {
        let mut hits = Vec::new();

        if subject_aa.len() < 3 {
            return hits;
        }

        for s_pos in 0..=(subject_aa.len() - 3) {
            if let Some(kmer_idx) = encode_aa_kmer(subject_aa, s_pos) {
                // Exact matches
                for &(query_idx, q_frame, q_pos) in &self.lookup[kmer_idx] {
                    hits.push(AaSeedHit {
                        query_idx: query_idx as usize,
                        query_frame: q_frame,
                        query_pos: q_pos,
                        subject_frame,
                        subject_pos: s_pos,
                        is_exact: true,
                    });
                }

                // Neighborhood matches (if enabled)
                if let Some(ref neighborhood) = self.neighborhood {
                    for &neighbor_idx in neighborhood.get_neighbors(kmer_idx) {
                        let neighbor_idx = neighbor_idx as usize;
                        if neighbor_idx == kmer_idx {
                            continue; // Skip exact match (already added)
                        }

                        for &(query_idx, q_frame, q_pos) in &self.lookup[neighbor_idx] {
                            hits.push(AaSeedHit {
                                query_idx: query_idx as usize,
                                query_frame: q_frame,
                                query_pos: q_pos,
                                subject_frame,
                                subject_pos: s_pos,
                                is_exact: false,
                            });
                        }
                    }
                }
            }
        }

        hits
    }

    /// Get the configuration
    pub fn config(&self) -> &AaWordFinderConfig {
        &self.config
    }

    /// Check if neighborhood finding is enabled
    pub fn uses_neighborhood(&self) -> bool {
        self.neighborhood.is_some()
    }
}

/// A seed hit for amino acid sequences
#[derive(Debug, Clone, Copy)]
pub struct AaSeedHit {
    /// Index of the query sequence
    pub query_idx: usize,
    /// Frame of the query (0-5 for 6 frames)
    pub query_frame: usize,
    /// Position in the query amino acid sequence (0-based)
    pub query_pos: usize,
    /// Frame of the subject (0-5 for 6 frames)
    pub subject_frame: usize,
    /// Position in the subject amino acid sequence (0-based)
    pub subject_pos: usize,
    /// Whether this is an exact match or neighborhood match
    pub is_exact: bool,
}

impl AaSeedHit {
    /// Get the diagonal of this hit (query_pos - subject_pos)
    pub fn diagonal(&self) -> isize {
        self.query_pos as isize - self.subject_pos as isize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_aa_kmer() {
        // Test encoding and decoding
        let kmer = [0u8, 1, 2]; // A, C, D in amino acid indices
        let idx = (0 * 676) + (1 * 26) + 2;

        assert_eq!(encode_aa_kmer(&kmer, 0), Some(idx));
        assert_eq!(decode_aa_kmer(idx), kmer);
    }

    #[test]
    fn test_encode_with_stop_codon() {
        let kmer = [0u8, STOP_CODON, 2];
        assert!(encode_aa_kmer(&kmer, 0).is_none());
    }

    #[test]
    fn test_kmer_score() {
        // Same k-mer should have positive score
        let kmer = [0u8, 1, 2];
        let score = get_kmer_score(&kmer, &kmer);
        assert!(score > 0);

        // Different k-mers should have lower score
        let kmer2 = [3u8, 4, 5];
        let score2 = get_kmer_score(&kmer, &kmer2);
        assert!(score2 < score);
    }

    #[test]
    fn test_neighborhood_index() {
        // Create with a moderate threshold
        let index = NeighborhoodIndex::new(11);

        // Each k-mer should at least have itself as a neighbor
        let kmer = [0u8, 1, 2];
        let kmer_idx = encode_aa_kmer(&kmer, 0).unwrap();
        let neighbors = index.get_neighbors(kmer_idx);

        // Should have at least one neighbor (itself)
        assert!(!neighbors.is_empty());

        // Self should be a neighbor
        assert!(neighbors.contains(&(kmer_idx as u16)));
    }

    #[test]
    fn test_aa_word_finder_exact() {
        let config = AaWordFinderConfig {
            use_neighborhood: false,
            ..Default::default()
        };

        let mut finder = AaWordFinder::new(config);

        // Add a query with a simple sequence
        let query_frames = vec![
            (vec![0u8, 1, 2, 3, 4], 1), // Frame 1
        ];
        finder.add_query(0, &query_frames);

        // Search with matching subject
        let subject = vec![0u8, 1, 2, 3, 4];
        let hits = finder.find_seeds(&subject, 0);

        // Should find exact matches
        assert!(!hits.is_empty());
        assert!(hits.iter().all(|h| h.is_exact));
    }

    #[test]
    fn test_aa_word_finder_neighborhood() {
        let config = AaWordFinderConfig {
            use_neighborhood: true,
            threshold: 11,
            ..Default::default()
        };

        let mut finder = AaWordFinder::new(config);

        // Add a query
        let query_frames = vec![(vec![0u8, 1, 2, 3, 4], 1)];
        finder.add_query(0, &query_frames);

        // Search with similar subject
        let subject = vec![0u8, 1, 2, 3, 4];
        let hits = finder.find_seeds(&subject, 0);

        // Should find matches (exact and possibly neighborhood)
        assert!(!hits.is_empty());
    }
}
