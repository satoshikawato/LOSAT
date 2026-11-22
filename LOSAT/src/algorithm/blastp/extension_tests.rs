use super::*;
use crate::utils::matrix::ncbistdaa;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1078-1127
#[test]
fn test_extend_one_hit_returns_ncbi_blast_ungapped_data_payload() {
    let query = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
    ];
    let subject = query.clone();

    let result = extend_one_hit(ScoringMatrix::Blosum62, &query, &subject, 1, 1, 7, 3)
        .expect("one-hit ungapped extension");

    assert_eq!(result.ungapped_data.q_start, 0);
    assert_eq!(result.ungapped_data.s_start, 0);
    assert_eq!(result.ungapped_data.length, 5);
    assert_eq!(
        result.ungapped_data.score,
        5 * protein_score(ScoringMatrix::Blosum62, ncbistdaa::A, ncbistdaa::A)
    );
    assert_eq!(result.s_last_off, 5);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1124-1158
#[test]
fn test_extend_two_hit_returns_ncbi_payload_when_left_reaches_first_hit() {
    let query = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
    ];
    let subject = query.clone();

    let result = extend_two_hit(ScoringMatrix::Blosum62, &query, &subject, 0, 3, 3, 7, 3)
        .expect("two-hit ungapped extension");

    assert!(result.right_extend);
    assert_eq!(result.ungapped_data.q_start, 0);
    assert_eq!(result.ungapped_data.s_start, 0);
    assert_eq!(result.ungapped_data.length, 6);
    assert_eq!(
        result.ungapped_data.score,
        6 * protein_score(ScoringMatrix::Blosum62, ncbistdaa::A, ncbistdaa::A)
    );
    assert_eq!(result.s_last_off, 6);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1129-1158
#[test]
fn test_extend_two_hit_returns_ncbi_payload_when_left_stops_before_first_hit() {
    let query = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::D,
        ncbistdaa::D,
        ncbistdaa::D,
        ncbistdaa::C,
        ncbistdaa::C,
        ncbistdaa::C,
    ];
    let subject = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::W,
        ncbistdaa::W,
        ncbistdaa::W,
        ncbistdaa::C,
        ncbistdaa::C,
        ncbistdaa::C,
    ];

    let result = extend_two_hit(ScoringMatrix::Blosum62, &query, &subject, 0, 6, 6, 2, 3)
        .expect("two-hit ungapped extension");

    assert!(!result.right_extend);
    assert_eq!(result.ungapped_data.q_start, 6);
    assert_eq!(result.ungapped_data.s_start, 6);
    assert_eq!(result.ungapped_data.length, 3);
    assert_eq!(
        result.ungapped_data.score,
        3 * protein_score(ScoringMatrix::Blosum62, ncbistdaa::C, ncbistdaa::C)
    );
    assert_eq!(result.s_last_off, 9);
}
