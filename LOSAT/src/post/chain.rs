use crate::common::Hit;
use rustc_hash::FxHashMap;
use std::sync::Arc;

/// Configuration for HSP chaining
#[derive(Debug, Clone, Copy)]
pub struct ChainConfig {
    /// Whether to enable HSP chaining
    pub enable_chaining: bool,
    /// Maximum gap allowed between HSPs for chaining
    pub max_gap: usize,
    /// Maximum diagonal drift allowed for chaining
    pub max_diag_drift: isize,
    /// Overlap threshold for filtering redundant HSPs
    pub overlap_threshold: f64,
}

impl Default for ChainConfig {
    fn default() -> Self {
        Self {
            enable_chaining: true,
            max_gap: 1000,
            max_diag_drift: 100,
            overlap_threshold: 0.5,
        }
    }
}

impl ChainConfig {
    /// Configuration for NCBI BLAST compatible behavior (no chaining)
    pub fn ncbi_compat() -> Self {
        Self {
            enable_chaining: false,
            max_gap: 0,
            max_diag_drift: 0,
            overlap_threshold: 0.0,
        }
    }
}

/// Chain and filter HSPs for a query-subject pair
///
/// This implements the BLEMIR chaining algorithm that merges nearby HSPs
/// into longer alignments. When `enable_chaining` is false, HSPs are
/// returned without merging (NCBI-compatible behavior).
pub fn chain_and_filter_hsps(hits: Vec<Hit>, config: &ChainConfig) -> Vec<Hit> {
    if hits.is_empty() {
        return hits;
    }

    if !config.enable_chaining {
        // NCBI-compatible: just filter overlapping hits without chaining
        return filter_overlapping_hsps(hits, config.overlap_threshold);
    }

    // Group hits by query-subject pair
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
    //    Int4 hspcnt; /**< Number of HSPs saved */
    //    ...
    // } BlastHSPList;
    // ```
    let mut grouped: FxHashMap<(Arc<str>, Arc<str>), Vec<Hit>> = FxHashMap::default();
    for hit in hits {
        let key = (hit.query_id.clone(), hit.subject_id.clone());
        grouped.entry(key).or_default().push(hit);
    }

    let mut result = Vec::new();

    for ((_query_id, _subject_id), mut pair_hits) in grouped {
        // Sort by query start position
        pair_hits.sort_by_key(|h| h.q_start);

        // Chain nearby HSPs
        let chained = chain_hsps_for_pair(pair_hits, config);

        // Filter overlapping chains
        let filtered = filter_overlapping_hsps(chained, config.overlap_threshold);

        result.extend(filtered);
    }

    // Sort final results by bit score descending
    result.sort_by(|a, b| {
        b.bit_score
            .partial_cmp(&a.bit_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    result
}

/// Chain HSPs for a single query-subject pair
fn chain_hsps_for_pair(hits: Vec<Hit>, config: &ChainConfig) -> Vec<Hit> {
    if hits.is_empty() {
        return hits;
    }

    let mut chains: Vec<HspChain> = Vec::new();

    for hit in hits {
        let hit_diag = hit.q_start as isize - hit.s_start as isize;

        // Try to extend an existing chain
        let mut merged = false;
        for chain in &mut chains {
            // Check if this hit can extend the chain
            let diag_diff = (hit_diag - chain.diagonal).abs();
            let q_gap = hit.q_start.saturating_sub(chain.q_end);
            let s_gap = hit.s_start.saturating_sub(chain.s_end);
            let gap = q_gap.max(s_gap);

            if diag_diff <= config.max_diag_drift as isize && gap <= config.max_gap {
                // Merge into this chain
                chain.extend(&hit);
                merged = true;
                break;
            }
        }

        if !merged {
            // Start a new chain
            chains.push(HspChain::from_hit(&hit));
        }
    }

    // Convert chains back to hits
    chains.into_iter().map(|c| c.to_hit()).collect()
}

/// Filter overlapping HSPs, keeping higher-scoring ones
fn filter_overlapping_hsps(hits: Vec<Hit>, overlap_threshold: f64) -> Vec<Hit> {
    if hits.is_empty() || overlap_threshold <= 0.0 {
        return hits;
    }

    // Sort by bit score descending
    let mut sorted = hits;
    sorted.sort_by(|a, b| {
        b.bit_score
            .partial_cmp(&a.bit_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut kept: Vec<Hit> = Vec::new();

    for hit in sorted {
        let dominated = kept.iter().any(|kept_hit| {
            // Must be same query-subject pair
            if hit.query_id != kept_hit.query_id || hit.subject_id != kept_hit.subject_id {
                return false;
            }

            // Calculate overlap fractions
            let q_overlap =
                calculate_overlap(hit.q_start, hit.q_end, kept_hit.q_start, kept_hit.q_end);
            let s_overlap =
                calculate_overlap(hit.s_start, hit.s_end, kept_hit.s_start, kept_hit.s_end);

            let hit_q_len = (hit.q_end - hit.q_start + 1) as f64;
            let hit_s_len = (hit.s_end - hit.s_start + 1) as f64;

            let q_frac = q_overlap as f64 / hit_q_len;
            let s_frac = s_overlap as f64 / hit_s_len;

            q_frac > overlap_threshold && s_frac > overlap_threshold
        });

        if !dominated {
            kept.push(hit);
        }
    }

    kept
}

/// Calculate overlap between two intervals
fn calculate_overlap(start1: usize, end1: usize, start2: usize, end2: usize) -> usize {
    let overlap_start = start1.max(start2);
    let overlap_end = end1.min(end2);

    if overlap_start <= overlap_end {
        overlap_end - overlap_start + 1
    } else {
        0
    }
}

/// Internal representation of an HSP chain during merging
struct HspChain {
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
    //    Int4 hspcnt; /**< Number of HSPs saved */
    //    ...
    // } BlastHSPList;
    // ```
    query_id: Arc<str>,
    subject_id: Arc<str>,
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    diagonal: isize,
    total_score: f64,
    total_matches: usize,
    total_mismatches: usize,
    total_gap_opens: usize,
    total_length: usize,
    min_evalue: f64,
}

impl HspChain {
    fn from_hit(hit: &Hit) -> Self {
        let diagonal = hit.q_start as isize - hit.s_start as isize;
        Self {
            query_id: hit.query_id.clone(),
            subject_id: hit.subject_id.clone(),
            q_start: hit.q_start,
            q_end: hit.q_end,
            s_start: hit.s_start,
            s_end: hit.s_end,
            diagonal,
            total_score: hit.bit_score,
            total_matches: (hit.identity * hit.length as f64 / 100.0).round() as usize,
            total_mismatches: hit.mismatch,
            total_gap_opens: hit.gapopen,
            total_length: hit.length,
            min_evalue: hit.e_value,
        }
    }

    fn extend(&mut self, hit: &Hit) {
        // Update boundaries
        self.q_end = self.q_end.max(hit.q_end);
        self.s_end = self.s_end.max(hit.s_end);

        // Accumulate statistics
        self.total_score += hit.bit_score;
        self.total_matches += (hit.identity * hit.length as f64 / 100.0).round() as usize;
        self.total_mismatches += hit.mismatch;
        self.total_gap_opens += hit.gapopen;
        self.total_length += hit.length;
        self.min_evalue = self.min_evalue.min(hit.e_value);

        // Update diagonal (use weighted average)
        let hit_diag = hit.q_start as isize - hit.s_start as isize;
        self.diagonal = (self.diagonal + hit_diag) / 2;
    }

    fn to_hit(self) -> Hit {
        let length = (self.q_end - self.q_start + 1).max(self.s_end - self.s_start + 1);
        let identity = if length > 0 {
            100.0 * self.total_matches as f64 / length as f64
        } else {
            0.0
        };

        Hit {
            query_id: self.query_id,
            subject_id: self.subject_id,
            identity,
            length,
            mismatch: self.total_mismatches,
            gapopen: self.total_gap_opens,
            q_start: self.q_start,
            q_end: self.q_end,
            s_start: self.s_start,
            s_end: self.s_end,
            e_value: self.min_evalue,
            bit_score: self.total_score,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
            // ```c
            // if (hsp->query.frame != hsp->subject.frame) {
            //    *q_end = query_length - hsp->query.offset;
            //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
            // }
            // ```
            query_frame: 1,
            query_length: 0,
            q_idx: 0,
            s_idx: 0,
            raw_score: (self.total_score * 2.0) as i32,
            gap_info: None,
        }
    }
}

/// Cluster HSPs by diagonal for efficient processing
pub fn cluster_by_diagonal(hits: &[Hit], max_diag_drift: isize) -> Vec<Vec<usize>> {
    if hits.is_empty() {
        return Vec::new();
    }

    // Calculate diagonals
    let diagonals: Vec<isize> = hits
        .iter()
        .map(|h| h.q_start as isize - h.s_start as isize)
        .collect();

    // Sort indices by diagonal
    let mut indices: Vec<usize> = (0..hits.len()).collect();
    indices.sort_by_key(|&i| diagonals[i]);

    // Cluster by diagonal proximity
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut current_cluster: Vec<usize> = vec![indices[0]];
    let mut current_diag = diagonals[indices[0]];

    for &idx in &indices[1..] {
        let diag = diagonals[idx];
        if (diag - current_diag).abs() <= max_diag_drift {
            current_cluster.push(idx);
        } else {
            if !current_cluster.is_empty() {
                clusters.push(current_cluster);
            }
            current_cluster = vec![idx];
        }
        current_diag = diag;
    }

    if !current_cluster.is_empty() {
        clusters.push(current_cluster);
    }

    clusters
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hit(q_start: usize, q_end: usize, s_start: usize, s_end: usize, bit_score: f64) -> Hit {
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
        //                       Set to 0 if not applicable */
        //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
        //    Int4 hspcnt; /**< Number of HSPs saved */
        //    ...
        // } BlastHSPList;
        // ```
        Hit {
            query_id: "q1".into(),
            subject_id: "s1".into(),
            identity: 90.0,
            length: q_end - q_start + 1,
            mismatch: 0,
            gapopen: 0,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value: 1e-10,
            bit_score,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
            // ```c
            // if (hsp->query.frame != hsp->subject.frame) {
            //    *q_end = query_length - hsp->query.offset;
            //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
            // }
            // ```
            query_frame: 1,
            query_length: 0,
            q_idx: 0,
            s_idx: 0,
            raw_score: (bit_score * 2.0) as i32,
            gap_info: None,
        }
    }

    #[test]
    fn test_chain_nearby_hsps() {
        let hits = vec![
            make_hit(1, 100, 1, 100, 50.0),
            make_hit(110, 200, 110, 200, 50.0), // Close, same diagonal
        ];

        let config = ChainConfig::default();
        let chained = chain_and_filter_hsps(hits, &config);

        // Should be merged into one chain
        assert_eq!(chained.len(), 1);
        assert_eq!(chained[0].q_start, 1);
        assert_eq!(chained[0].q_end, 200);
    }

    #[test]
    fn test_no_chain_distant_hsps() {
        let hits = vec![
            make_hit(1, 100, 1, 100, 50.0),
            make_hit(2000, 2100, 2000, 2100, 50.0), // Far apart
        ];

        let config = ChainConfig::default();
        let chained = chain_and_filter_hsps(hits, &config);

        // Should remain separate
        assert_eq!(chained.len(), 2);
    }

    #[test]
    fn test_ncbi_compat_no_chaining() {
        let hits = vec![
            make_hit(1, 100, 1, 100, 50.0),
            make_hit(110, 200, 110, 200, 50.0),
        ];

        let config = ChainConfig::ncbi_compat();
        let result = chain_and_filter_hsps(hits, &config);

        // Should not be merged (NCBI mode)
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_filter_overlapping() {
        let hits = vec![
            make_hit(1, 100, 1, 100, 100.0),
            make_hit(10, 90, 10, 90, 50.0), // Significantly overlaps
        ];

        let filtered = filter_overlapping_hsps(hits, 0.5);

        // Lower-scoring overlapping hit should be removed
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].bit_score, 100.0);
    }

    #[test]
    fn test_cluster_by_diagonal() {
        let hits = vec![
            make_hit(10, 20, 10, 20, 10.0),   // diag = 0
            make_hit(30, 40, 30, 40, 10.0),   // diag = 0
            make_hit(100, 110, 50, 60, 10.0), // diag = 50
        ];

        let clusters = cluster_by_diagonal(&hits, 10);

        // Should have 2 clusters (diag 0 and diag 50)
        assert_eq!(clusters.len(), 2);
    }
}
