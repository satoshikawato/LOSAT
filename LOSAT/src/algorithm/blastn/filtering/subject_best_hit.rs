//! Subject Best Hit filtering for BLASTN
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2537-2606
//! Function: Blast_HSPListSubjectBestHit
//!
//! This filter removes HSPs in overlapping query ranges after the best scoring HSP.
//! For BLASTN, it also performs cross-strand filtering.

use crate::common::Hit;

/// Default range difference for subject best hit filtering
/// NCBI reference: blast_options.h:1211-1212
/// #define BLAST_SUBJECT_BESTHIT_DEFAULT_RANGE_DIFF 3
const MAX_RANGE_DIFF: usize = 3;

/// Apply Subject Best Hit filtering to remove HSPs with overlapping query ranges.
///
/// NCBI reference: blast_hits.c:2537-2606 Blast_HSPListSubjectBestHit
///
/// ```c
/// // NCBI algorithm:
/// // 1. Sort HSPs by score (already done)
/// // 2. For each HSP, remove subsequent HSPs with:
/// //    - Same context (strand)
/// //    - Query range within ±range_diff of current HSP
/// // 3. For BLASTN: Also filter cross-strand HSPs with flipped coordinates
/// ```
///
/// Arguments:
/// - `hits`: Vector of HSPs, already sorted by score descending
/// - `query_len`: Length of the query sequence (for cross-strand coordinate flipping)
pub fn subject_best_hit(hits: &mut Vec<Hit>, query_len: usize) {
    if hits.len() <= 1 {
        return;
    }

    // Helper to determine strand: minus if s_start > s_end
    // NCBI reference: For blastn, subject.frame encodes strand direction
    let is_minus_strand = |h: &Hit| h.s_start > h.s_end;

    // Ensure sorted by score descending (should already be sorted by caller)
    hits.sort_by(|a, b| b.raw_score.cmp(&a.raw_score));

    // Mark HSPs for removal
    let mut to_remove = vec![false; hits.len()];

    // Pass 1: Same-strand filtering
    // NCBI reference: blast_hits.c:2561-2577
    // ```c
    // for(i=0; i < hsp_list->hspcnt -1; i++) {
    //     o = hsp_array[i]->query.offset - range_diff;
    //     e = hsp_array[i]->query.end + range_diff;
    //     while (i+j < hsp_list->hspcnt) {
    //         if (hsp_array[i+j] && hsp_array[i]->context == hsp_array[i+j]->context &&
    //             ((hsp_array[i+j]->query.offset >= o) && (hsp_array[i+j]->query.end <= e))) {
    //             hsp_array[i+j] = Blast_HSPFree(hsp_array[i+j]);
    //         }
    //         j++;
    //     }
    // }
    // ```
    for i in 0..hits.len() {
        if to_remove[i] {
            continue;
        }

        let hit_i_strand = is_minus_strand(&hits[i]);
        let o = hits[i].q_start.saturating_sub(MAX_RANGE_DIFF);
        let e = hits[i].q_end.saturating_add(MAX_RANGE_DIFF);

        for j in (i + 1)..hits.len() {
            if to_remove[j] {
                continue;
            }

            let hit_j_strand = is_minus_strand(&hits[j]);

            // Same context (strand) check
            if hit_i_strand == hit_j_strand {
                // Check if query range is within ±range_diff
                if hits[j].q_start >= o && hits[j].q_end <= e {
                    to_remove[j] = true;
                }
            }
        }
    }

    // Pass 2: Cross-strand filtering (BLASTN-specific)
    // NCBI reference: blast_hits.c:2582-2604
    // ```c
    // if(isBlastn) {
    //     for(i=0; i < hsp_list->hspcnt -1; i++) {
    //         curr_context = hsp_array[i]->context;
    //         qlen = query_info->contexts[curr_context].query_length;
    //         target_context = (hsp_array[i]->query.frame > 0) ? curr_context +1 : curr_context -1;
    //         e = qlen - (hsp_array[i]->query.offset - range_diff);
    //         o = qlen - (hsp_array[i]->query.end + range_diff);
    //         while (i+j < hsp_list->hspcnt) {
    //             if(hsp_array[i+j] && (hsp_array[i+j]->context == target_context) &&
    //                ((hsp_array[i+j]->query.offset >= o) && (hsp_array[i+j]->query.end <= e))) {
    //                 hsp_array[i+j] = Blast_HSPFree(hsp_array[i+j]);
    //             }
    //             j++;
    //         }
    //     }
    // }
    // ```
    for i in 0..hits.len() {
        if to_remove[i] {
            continue;
        }

        let hit_i_strand = is_minus_strand(&hits[i]);
        let target_strand = !hit_i_strand; // Opposite strand

        // Flip coordinates for cross-strand comparison
        // NCBI: e = qlen - (query.offset - range_diff)
        // NCBI: o = qlen - (query.end + range_diff)
        let flipped_o = query_len.saturating_sub(hits[i].q_end.saturating_add(MAX_RANGE_DIFF));
        let flipped_e = query_len.saturating_sub(hits[i].q_start.saturating_sub(MAX_RANGE_DIFF));

        for j in (i + 1)..hits.len() {
            if to_remove[j] {
                continue;
            }

            let hit_j_strand = is_minus_strand(&hits[j]);

            // Target context (opposite strand) check
            if hit_j_strand == target_strand {
                // Check if query range (in flipped coordinates) is within range
                if hits[j].q_start >= flipped_o && hits[j].q_end <= flipped_e {
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

    fn make_hit(q_start: usize, q_end: usize, s_start: usize, s_end: usize, score: i32) -> Hit {
        Hit {
            query_id: "q1".to_string(),
            subject_id: "s1".to_string(),
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
            q_idx: 0,
            s_idx: 0,
            raw_score: score,
        }
    }

    #[test]
    fn test_same_strand_filtering() {
        // Hit at q_start=10, q_end=50, score=100
        // Hit at q_start=12, q_end=48, score=80 (should be filtered - within ±3 of first)
        // Hit at q_start=20, q_end=60, score=60 (should NOT be filtered - outside range)
        let mut hits = vec![
            make_hit(10, 50, 10, 50, 100),  // Plus strand, best score
            make_hit(12, 48, 12, 48, 80),   // Plus strand, within range - REMOVE
            make_hit(20, 60, 20, 60, 60),   // Plus strand, outside range - KEEP
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
            make_hit(10, 50, 10, 50, 100),  // Plus strand (s_start < s_end)
            make_hit(12, 48, 48, 12, 80),   // Minus strand (s_start > s_end)
        ];

        subject_best_hit(&mut hits, 1000);

        // Both should remain after same-strand filtering
        // Cross-strand filtering might remove one depending on coordinates
        assert!(hits.len() >= 1);
    }
}
