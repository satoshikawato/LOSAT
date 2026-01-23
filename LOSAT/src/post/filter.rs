use crate::common::Hit;

/// Filter hits by E-value threshold
pub fn filter_by_evalue(hits: Vec<Hit>, max_evalue: f64) -> Vec<Hit> {
    hits.into_iter()
        .filter(|h| h.e_value <= max_evalue)
        .collect()
}

/// Filter hits by minimum bit score
pub fn filter_by_bit_score(hits: Vec<Hit>, min_bit_score: f64) -> Vec<Hit> {
    hits.into_iter()
        .filter(|h| h.bit_score >= min_bit_score)
        .collect()
}

/// Filter hits by minimum identity percentage
pub fn filter_by_identity(hits: Vec<Hit>, min_identity: f64) -> Vec<Hit> {
    hits.into_iter()
        .filter(|h| h.identity >= min_identity)
        .collect()
}

/// Filter hits by minimum alignment length
pub fn filter_by_length(hits: Vec<Hit>, min_length: usize) -> Vec<Hit> {
    hits.into_iter()
        .filter(|h| h.length >= min_length)
        .collect()
}

/// Remove hits that overlap significantly with higher-scoring hits
///
/// Two hits are considered overlapping if they share more than `overlap_threshold`
/// fraction of their aligned positions on both query and subject
pub fn filter_overlapping(hits: Vec<Hit>, overlap_threshold: f64) -> Vec<Hit> {
    if hits.is_empty() || overlap_threshold <= 0.0 {
        return hits;
    }

    // Sort by bit score descending
    let mut sorted_hits = hits;
    sorted_hits.sort_by(|a, b| {
        b.bit_score
            .partial_cmp(&a.bit_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut kept_hits: Vec<Hit> = Vec::new();

    for hit in sorted_hits {
        let dominated = kept_hits.iter().any(|kept| {
            // Check if same query-subject pair
            // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
            // ```c
            // typedef struct BlastHSPList {
            //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
            //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
            //                       Set to 0 if not applicable */
            // } BlastHSPList;
            // ```
            if hit.q_idx != kept.q_idx || hit.s_idx != kept.s_idx {
                return false;
            }

            // Calculate overlap on query
            let q_overlap = calculate_overlap(hit.q_start, hit.q_end, kept.q_start, kept.q_end);

            // Calculate overlap on subject
            let s_overlap = calculate_overlap(hit.s_start, hit.s_end, kept.s_start, kept.s_end);

            // Check if overlap exceeds threshold on both
            let hit_q_len = (hit.q_end - hit.q_start + 1) as f64;
            let hit_s_len = (hit.s_end - hit.s_start + 1) as f64;

            let q_overlap_frac = q_overlap as f64 / hit_q_len;
            let s_overlap_frac = s_overlap as f64 / hit_s_len;

            q_overlap_frac > overlap_threshold && s_overlap_frac > overlap_threshold
        });

        if !dominated {
            kept_hits.push(hit);
        }
    }

    kept_hits
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

/// Keep only the top N hits per query-subject pair
pub fn keep_top_n_per_pair(hits: Vec<Hit>, n: usize) -> Vec<Hit> {
    use std::collections::HashMap;

    if n == 0 {
        return Vec::new();
    }

    // Sort by bit score descending
    let mut sorted_hits = hits;
    sorted_hits.sort_by(|a, b| {
        b.bit_score
            .partial_cmp(&a.bit_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Count hits per query-subject pair
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
    let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
    let mut kept_hits = Vec::new();

    for hit in sorted_hits {
        let key = (hit.q_idx, hit.s_idx);
        let count = counts.entry(key).or_insert(0);

        if *count < n {
            *count += 1;
            kept_hits.push(hit);
        }
    }

    kept_hits
}

/// Keep only the top N hits overall
pub fn keep_top_n(hits: Vec<Hit>, n: usize) -> Vec<Hit> {
    if n == 0 {
        return Vec::new();
    }

    let mut sorted_hits = hits;
    sorted_hits.sort_by(|a, b| {
        b.bit_score
            .partial_cmp(&a.bit_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    sorted_hits.truncate(n);
    sorted_hits
}

/// Configuration for hit filtering
#[derive(Debug, Clone)]
pub struct FilterConfig {
    /// Maximum E-value (None = no filter)
    pub max_evalue: Option<f64>,
    /// Minimum bit score (None = no filter)
    pub min_bit_score: Option<f64>,
    /// Minimum identity percentage (None = no filter)
    pub min_identity: Option<f64>,
    /// Minimum alignment length (None = no filter)
    pub min_length: Option<usize>,
    /// Overlap threshold for filtering redundant hits (0.0 = disabled)
    pub overlap_threshold: f64,
    /// Maximum hits per query-subject pair (None = no limit)
    pub max_hits_per_pair: Option<usize>,
    /// Maximum total hits (None = no limit)
    pub max_hits: Option<usize>,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            max_evalue: Some(10.0),
            min_bit_score: None,
            min_identity: None,
            min_length: None,
            overlap_threshold: 0.5,
            max_hits_per_pair: None,
            max_hits: Some(500),
        }
    }
}

impl FilterConfig {
    /// Configuration for NCBI BLAST compatible filtering (no overlap filtering)
    pub fn ncbi_compat() -> Self {
        Self {
            max_evalue: Some(10.0),
            min_bit_score: None,
            min_identity: None,
            min_length: None,
            overlap_threshold: 0.0, // Disabled
            max_hits_per_pair: None,
            max_hits: Some(500),
        }
    }
}

/// Apply all configured filters to hits
pub fn apply_filters(hits: Vec<Hit>, config: &FilterConfig) -> Vec<Hit> {
    let mut filtered = hits;

    if let Some(max_ev) = config.max_evalue {
        filtered = filter_by_evalue(filtered, max_ev);
    }

    if let Some(min_bs) = config.min_bit_score {
        filtered = filter_by_bit_score(filtered, min_bs);
    }

    if let Some(min_id) = config.min_identity {
        filtered = filter_by_identity(filtered, min_id);
    }

    if let Some(min_len) = config.min_length {
        filtered = filter_by_length(filtered, min_len);
    }

    if config.overlap_threshold > 0.0 {
        filtered = filter_overlapping(filtered, config.overlap_threshold);
    }

    if let Some(n) = config.max_hits_per_pair {
        filtered = keep_top_n_per_pair(filtered, n);
    }

    if let Some(n) = config.max_hits {
        filtered = keep_top_n(filtered, n);
    }

    filtered
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
    fn test_filter_by_evalue() {
        let hits = vec![
            Hit {
                e_value: 1e-10,
                ..make_hit(1, 100, 1, 100, 100.0)
            },
            Hit {
                e_value: 1e-5,
                ..make_hit(1, 100, 1, 100, 50.0)
            },
            Hit {
                e_value: 1.0,
                ..make_hit(1, 100, 1, 100, 10.0)
            },
        ];

        let filtered = filter_by_evalue(hits, 1e-6);
        assert_eq!(filtered.len(), 1);
    }

    #[test]
    fn test_filter_overlapping() {
        let hits = vec![
            make_hit(1, 100, 1, 100, 100.0),    // Best hit
            make_hit(10, 90, 10, 90, 80.0),     // Overlaps significantly
            make_hit(200, 300, 200, 300, 90.0), // No overlap
        ];

        let filtered = filter_overlapping(hits, 0.5);
        assert_eq!(filtered.len(), 2); // First and third should be kept
    }

    #[test]
    fn test_calculate_overlap() {
        assert_eq!(calculate_overlap(1, 100, 50, 150), 51);
        assert_eq!(calculate_overlap(1, 100, 101, 200), 0);
        assert_eq!(calculate_overlap(1, 100, 1, 100), 100);
    }

    #[test]
    fn test_keep_top_n() {
        let hits = vec![
            make_hit(1, 100, 1, 100, 100.0),
            make_hit(1, 100, 1, 100, 80.0),
            make_hit(1, 100, 1, 100, 90.0),
        ];

        let top = keep_top_n(hits, 2);
        assert_eq!(top.len(), 2);
        assert_eq!(top[0].bit_score, 100.0);
        assert_eq!(top[1].bit_score, 90.0);
    }
}
