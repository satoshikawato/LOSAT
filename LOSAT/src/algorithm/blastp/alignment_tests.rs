use crate::utils::matrix::ncbistdaa;

use super::*;

#[test]
fn test_build_alignment_view_ungapped() {
    let hit = Hit {
        identity: 75.0,
        length: 4,
        mismatch: 1,
        gapopen: 0,
        q_start: 1,
        q_end: 4,
        s_start: 1,
        s_end: 4,
        e_value: 1e-10,
        bit_score: 50.0,
        num_ident: 3,
        query_frame: 0,
        query_length: 4,
        q_idx: 0,
        s_idx: 0,
        raw_score: 20,
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-143
        // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: 4,
    };

    let query = [1u8, 2, 3, 4];
    let subject = [1u8, 2, 5, 4];
    let view = build_alignment_view_with_matrix(&query, &subject, &hit, ScoringMatrix::Blosum62)
        .expect("alignment");
    assert_eq!(view.alignment_len, 4);
    assert_eq!(view.num_ident, 3);
    assert_eq!(view.gap_letters, 0);
}

#[test]
fn test_build_alignment_view_gapped() {
    let hit = Hit {
        identity: 75.0,
        length: 4,
        mismatch: 0,
        gapopen: 1,
        q_start: 1,
        q_end: 3,
        s_start: 1,
        s_end: 4,
        e_value: 1e-10,
        bit_score: 50.0,
        num_ident: 3,
        query_frame: 0,
        query_length: 3,
        q_idx: 0,
        s_idx: 0,
        raw_score: 20,
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-143
        // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: Some(vec![
            GapEditOp::Sub(2),
            GapEditOp::Del(1),
            GapEditOp::Sub(1),
        ]),
        num_positives: 3,
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:115-121
    let query = [ncbistdaa::A, ncbistdaa::R, ncbistdaa::N];
    let subject = [ncbistdaa::A, ncbistdaa::R, ncbistdaa::D, ncbistdaa::C];
    let view = build_alignment_view_with_matrix(&query, &subject, &hit, ScoringMatrix::Blosum62)
        .expect("alignment");
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/gapinfo.h:45-51
    assert_eq!(view.query_aln, "AR-N");
    assert_eq!(view.subject_aln, "ARDC");
    assert_eq!(view.gap_letters, 1);
    assert_eq!(view.gap_opens, 1);
}

// An ungapped HSP of length 0 (kept by `Blast_HSPTest`) has no column, also at the first letter
// of the subject, whose range then ends before it.
#[test]
fn test_build_alignment_view_of_length_zero() {
    let hit = Hit {
        identity: 0.0,
        length: 0,
        mismatch: 0,
        gapopen: 0,
        q_start: 3,
        q_end: 2,
        s_start: 1,
        s_end: 0,
        e_value: 1.0,
        bit_score: 10.0,
        num_ident: 0,
        query_frame: 0,
        query_length: 4,
        q_idx: 0,
        s_idx: 0,
        raw_score: 20,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: Some(vec![GapEditOp::Sub(0)]),
        num_positives: 0,
    };
    let view =
        build_alignment_view_with_matrix(&[1, 2, 3, 4], &[1, 2], &hit, ScoringMatrix::Blosum62)
            .expect("alignment");
    assert_eq!(view.query_aln, "");
    assert_eq!(view.subject_aln, "");
    assert_eq!(
        (view.num_ident, view.num_positives, view.alignment_len),
        (0, 0, 0)
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:282-300
// A Seq-align that runs past the end of its query (a PHI-BLAST HSP of a long subject's merged
// parts) shows the gap character for the letters the query lacks, in a substitution's column and
// in an insertion's.
#[test]
fn test_build_alignment_view_past_the_query_end_shows_gaps() {
    let hit = Hit {
        identity: 40.0,
        length: 5,
        mismatch: 2,
        gapopen: 1,
        q_start: 1,
        q_end: 3,
        s_start: 1,
        s_end: 4,
        e_value: 1e-10,
        bit_score: 50.0,
        num_ident: 2,
        query_frame: 0,
        query_length: 3,
        q_idx: 0,
        s_idx: 0,
        raw_score: 20,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: Some(vec![GapEditOp::Sub(4), GapEditOp::Ins(1)]),
        num_positives: 2,
    };
    let query = [ncbistdaa::A, ncbistdaa::R, ncbistdaa::N];
    let subject = [ncbistdaa::A, ncbistdaa::R, ncbistdaa::D, ncbistdaa::C];
    let view = build_alignment_view_with_matrix(&query, &subject, &hit, ScoringMatrix::Blosum62)
        .expect("alignment");
    assert_eq!(view.query_aln, "ARN--");
    assert_eq!(view.subject_aln, "ARDC-");
    assert_eq!(view.alignment_len, 5);
    assert_eq!(view.gap_letters, 1);
}
