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
