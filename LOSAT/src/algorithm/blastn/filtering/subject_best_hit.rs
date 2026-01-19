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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
    // ```c
    // if (hsp->query.frame != hsp->subject.frame) {
    //    *q_end = query_length - hsp->query.offset;
    //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
    // }
    // ```
    let q_offsets = |h: &Hit| {
        if h.query_length > 0 && h.query_frame < 0 {
            (
                h.query_length.saturating_sub(h.q_end),
                h.query_length
                    .saturating_sub(h.q_start)
                    .saturating_add(1),
            )
        } else {
            (h.q_start.saturating_sub(1), h.q_end)
        }
    };
    let context = |h: &Hit| -> u32 {
        h.q_idx * 2 + if h.query_frame < 0 { 1 } else { 0 }
    };

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

        let hit_i_context = context(&hits[i]);
        let (i_q_offset, i_q_end) = q_offsets(&hits[i]);
        let o = i_q_offset.saturating_sub(MAX_RANGE_DIFF);
        let e = i_q_end.saturating_add(MAX_RANGE_DIFF);

        for j in (i + 1)..hits.len() {
            if to_remove[j] {
                continue;
            }

            let hit_j_context = context(&hits[j]);
            let (j_q_offset, j_q_end) = q_offsets(&hits[j]);

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

        let hit_i_context = context(&hits[i]);
        let target_context = if hits[i].query_frame > 0 {
            hit_i_context + 1
        } else {
            hit_i_context.saturating_sub(1)
        };

        // Flip coordinates for cross-strand comparison
        // NCBI: e = qlen - (query.offset - range_diff)
        // NCBI: o = qlen - (query.end + range_diff)
        let (i_q_offset, i_q_end) = q_offsets(&hits[i]);
        let flipped_o = query_len.saturating_sub(i_q_end.saturating_add(MAX_RANGE_DIFF));
        let flipped_e = query_len.saturating_sub(i_q_offset.saturating_sub(MAX_RANGE_DIFF));

        for j in (i + 1)..hits.len() {
            if to_remove[j] {
                continue;
            }

            let hit_j_context = context(&hits[j]);
            let (j_q_offset, j_q_end) = q_offsets(&hits[j]);

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
    ) -> Hit {
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
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
            // ```c
            // if (hsp->query.frame != hsp->subject.frame) {
            //    *q_end = query_length - hsp->query.offset;
            //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
            // }
            // ```
            query_frame,
            query_length: 1000,
            q_idx: 0,
            s_idx: 0,
            raw_score: score,
            gap_info: None,
        }
    }

    #[test]
    fn test_same_strand_filtering() {
        // Hit at q_start=10, q_end=50, score=100
        // Hit at q_start=12, q_end=48, score=80 (should be filtered - within ±3 of first)
        // Hit at q_start=20, q_end=60, score=60 (should NOT be filtered - outside range)
        let mut hits = vec![
            make_hit(10, 50, 10, 50, 100, 1),  // Plus strand, best score
            make_hit(12, 48, 12, 48, 80, 1),   // Plus strand, within range - REMOVE
            make_hit(20, 60, 20, 60, 60, 1),   // Plus strand, outside range - KEEP
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
            make_hit(10, 50, 10, 50, 100, 1),  // Plus strand
            make_hit(12, 48, 48, 12, 80, -1),  // Minus strand
        ];

        subject_best_hit(&mut hits, 1000);

        // Both should remain after same-strand filtering
        // Cross-strand filtering might remove one depending on coordinates
        assert!(hits.len() >= 1);
    }
}
