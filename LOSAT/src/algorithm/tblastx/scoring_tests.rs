use super::*;
use crate::utils::matrix::BLASTAA_SIZE;

// No NCBI counterpart: the table of BLOSUM62 holds the scores the search read before
// `-matrix`, for every two bytes.
#[test]
fn the_table_of_blosum62_holds_blosum62_score() {
    let scorer = MatrixScore::new(ScoringMatrix::Blosum62);
    for a in 0..=u8::MAX {
        for b in 0..=u8::MAX {
            assert_eq!(scorer.score(a, b), Blosum62Score.score(a, b), "{a} {b}");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/util/tables/sm_pam30.c (the PAM30 scores of A/A and W/W) and
// blast_stat.c:1559-1592 (the gap letter scores BLAST_SCORE_MIN).
#[test]
fn the_table_of_another_matrix_holds_its_scores() {
    use crate::utils::matrix::ncbistdaa::{A, GAP, W};
    let scorer = MatrixScore::new(ScoringMatrix::Pam30);
    assert_eq!(scorer.score(A, A), 6);
    assert_eq!(scorer.score(W, W), 13);
    assert_eq!(scorer.score(A, W), -13);
    assert_eq!(scorer.score(GAP, A), crate::utils::matrix::BLAST_SCORE_MIN);
    for a in 0..BLASTAA_SIZE as u8 {
        for b in 0..BLASTAA_SIZE as u8 {
            assert_eq!(
                scorer.score(a, b),
                protein_score(ScoringMatrix::Pam30, a, b)
            );
        }
    }
}
