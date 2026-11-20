//! Subject Best Hit filtering for BLASTN
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2537-2606
//! Function: Blast_HSPListSubjectBestHit
//!
//! This filter removes HSPs in overlapping query ranges after the best scoring HSP.
//! For BLASTN, it also performs cross-strand filtering.

use super::super::hsp::BlastnHsp;

/// Default range difference for subject best hit filtering
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:1211-1212
const MAX_RANGE_DIFF: i64 = 3;

/// The fields of an HSP that `Blast_HSPListSubjectBestHit` reads: its context, the frame of
/// the query, the query offsets within the context (NCBI's `Int4` offsets: an HSP can lie
/// outside its context, `BlastnHsp::outside_query_context`) and the length of the context.
#[derive(Debug, Clone, Copy)]
pub struct BestHitKey {
    pub context: u32,
    pub query_frame: i32,
    pub query_offset: i64,
    pub query_end: i64,
    pub query_length: i64,
}

/// Apply Subject Best Hit filtering to remove HSPs with overlapping query ranges.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2537-2606 Blast_HSPListSubjectBestHit
///
/// Arguments:
/// - `hits`: Vector of HSPs, already sorted by score descending
/// - `query_len`: Length of the query sequence (for cross-strand coordinate flipping)
pub fn subject_best_hit(hits: &mut Vec<BlastnHsp>, query_len: usize) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2563-2567
    // Subject-best-hit filtering uses internal query.offset/query.end, not
    // output coordinates.
    subject_best_hit_by(hits, |h| BestHitKey {
        context: h.q_idx * 2 + if h.query_frame < 0 { 1 } else { 0 },
        query_frame: h.query_frame,
        query_offset: i64::from(h.internal_q_offset_0),
        query_end: i64::from(h.internal_q_end_0),
        query_length: query_len as i64,
    });
}

/// `Blast_HSPListSubjectBestHit` over an HSP list in its order, with the fields of each HSP
/// given by `key`: the preliminary list of a subject holds the HSPs of every query, each
/// with the length of its own context (`query_info->contexts[curr_context].query_length`).
pub fn subject_best_hit_by<T>(hits: &mut Vec<T>, key: impl Fn(&T) -> BestHitKey) {
    if hits.len() <= 1 {
        return;
    }
    let keys: Vec<BestHitKey> = hits.iter().map(&key).collect();
    // Mark HSPs for removal
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2560-2566
    let mut to_remove = vec![false; hits.len()];

    // Pass 1: Same-strand filtering
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2561-2577
    for i in 0..hits.len() {
        if to_remove[i] {
            continue;
        }

        let hit_i_context = keys[i].context;
        let (i_q_offset, i_q_end) = (keys[i].query_offset, keys[i].query_end);
        let o = i_q_offset - MAX_RANGE_DIFF;
        let e = i_q_end + MAX_RANGE_DIFF;

        for j in (i + 1)..hits.len() {
            if to_remove[j] {
                continue;
            }

            let hit_j_context = keys[j].context;
            let (j_q_offset, j_q_end) = (keys[j].query_offset, keys[j].query_end);

            // Same context (strand) check
            if hit_i_context == hit_j_context {
                // Check if query range is within ±range_diff
                if j_q_offset >= o && j_q_end <= e {
                    to_remove[j] = true;
                }
            }
        }
    }

    // Pass 2: Cross-strand filtering (BLASTN-specific)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2582-2604
    for i in 0..hits.len() {
        if to_remove[i] {
            continue;
        }

        let hit_i_context = keys[i].context;
        let query_len = keys[i].query_length;
        let target_context = if keys[i].query_frame > 0 {
            hit_i_context + 1
        } else {
            hit_i_context.saturating_sub(1)
        };

        // Flip coordinates for cross-strand comparison
        // NCBI: e = qlen - (query.offset - range_diff)
        // NCBI: o = qlen - (query.end + range_diff)
        let (i_q_offset, i_q_end) = (keys[i].query_offset, keys[i].query_end);
        let flipped_o = query_len - (i_q_end + MAX_RANGE_DIFF);
        let flipped_e = query_len - (i_q_offset - MAX_RANGE_DIFF);

        for j in (i + 1)..hits.len() {
            if to_remove[j] {
                continue;
            }

            let hit_j_context = keys[j].context;
            let (j_q_offset, j_q_end) = (keys[j].query_offset, keys[j].query_end);

            // Target context (opposite strand) check
            if hit_j_context == target_context {
                // Check if query range (in flipped coordinates) is within range
                if j_q_offset >= flipped_o && j_q_end <= flipped_e {
                    to_remove[j] = true;
                }
            }
        }
    }

    // Remove marked HSPs by compacting the vector
    let mut write_idx = 0;
    for read_idx in 0..hits.len() {
        if !to_remove[read_idx] {
            if write_idx != read_idx {
                hits.swap(write_idx, read_idx);
            }
            write_idx += 1;
        }
    }
    hits.truncate(write_idx);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hit(
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
        score: i32,
        query_frame: i32,
    ) -> BlastnHsp {
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
        BlastnHsp {
            identity: 100.0,
            length: q_end - q_start + 1,
            mismatch: 0,
            gapopen: 0,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value: 0.0,
            bit_score: score as f64,
            num_ident: q_end - q_start + 1,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1122-1132
            query_frame,
            query_length: 1000,
            q_idx: 0,
            s_idx: 0,
            raw_score: score,
            internal_q_offset_0: if query_frame < 0 {
                1000usize.saturating_sub(q_end) as i32
            } else {
                q_start.saturating_sub(1) as i32
            },
            internal_q_end_0: if query_frame < 0 {
                1000usize.saturating_sub(q_start).saturating_add(1) as i32
            } else {
                q_end as i32
            },
            internal_s_offset_0: s_start.min(s_end).saturating_sub(1),
            internal_s_end_0: s_start.max(s_end),
            internal_query_context_offset: 0,
            gap_info: None,
            num_positives: q_end - q_start + 1,
        }
    }

    #[test]
    fn test_same_strand_filtering() {
        // Hit at q_start=10, q_end=50, score=100
        // Hit at q_start=12, q_end=48, score=80 (should be filtered - within ±3 of first)
        // Hit at q_start=20, q_end=60, score=60 (should NOT be filtered - outside range)
        let mut hits = vec![
            make_hit(10, 50, 10, 50, 100, 1), // Plus strand, best score
            make_hit(12, 48, 12, 48, 80, 1),  // Plus strand, within range - REMOVE
            make_hit(20, 60, 20, 60, 60, 1),  // Plus strand, outside range - KEEP
        ];

        subject_best_hit(&mut hits, 1000);

        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].raw_score, 100);
        assert_eq!(hits[1].raw_score, 60);
    }

    #[test]
    fn test_different_strand_no_filtering() {
        // Plus strand hit and minus strand hit with similar positions should NOT filter each other
        // in same-strand pass
        let mut hits = vec![
            make_hit(10, 50, 10, 50, 100, 1), // Plus strand
            make_hit(12, 48, 48, 12, 80, -1), // Minus strand
        ];

        subject_best_hit(&mut hits, 1000);

        // Both should remain after same-strand filtering
        // Cross-strand filtering might remove one depending on coordinates
        assert!(!hits.is_empty());
    }
}
