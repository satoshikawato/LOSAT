use super::*;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::core::blast_util::BLASTAA_SIZE;

// A PSSM that scores letter 1 (A) +5 and every other letter -4 at every position: the
// extensions then behave as on a sequence of A against BLOSUM62-like scores.
fn rows(len: usize) -> QueryPssm {
    let mut row = [-4i32; BLASTAA_SIZE];
    row[1] = 5;
    QueryPssm {
        rows: vec![row; len],
    }
}

// Hand-checked against `s_BlastAaExtendTwoHit` (aa_ungapped.c:1104-1157): the subject has
// A at 2..=9 and W elsewhere; hits at s 3 and s 6 on the diagonal q = s.
#[test]
fn two_hit_extends_both_ways_to_the_run_of_matches() {
    let pssm = rows(12);
    let mut subject = [20u8; 12];
    subject[2..10].fill(1);
    // s_left_off = last_hit + word_size = 3 + 3; the second hit at s 6, q 6.
    let hit = pssm_extend_two_hit(&pssm, &subject, 6, 6, 6, 16, 3);
    assert!(hit.right_extend);
    // The run of eight A from 2 to 9 scores 40.
    assert_eq!(
        hit.ungapped_data,
        BlastpUngappedData {
            q_start: 2,
            s_start: 2,
            length: 8,
            score: 40,
        }
    );
    // The right walk from s 9 reaches 40 at s 9 and drops by 4 at s 10 and 11 (8 < 16): it
    // ends at the subject's end, s_last_off 9 + 3 = 12.
    assert_eq!(hit.s_last_off, 12);
}

// The left extension that does not reach the first hit leaves the right side alone.
#[test]
fn two_hit_without_reaching_the_first_hit_does_not_extend_right() {
    let pssm = rows(20);
    let mut subject = [20u8; 20];
    subject[10..13].fill(1);
    // First hit far to the left (s_left_off 3): the walk left from s 9 drops after the
    // word.
    let hit = pssm_extend_two_hit(&pssm, &subject, 3, 10, 10, 7, 3);
    assert!(!hit.right_extend);
    assert_eq!(hit.ungapped_data.q_start, 10);
    assert_eq!(hit.ungapped_data.length, 3);
    assert_eq!(hit.ungapped_data.score, 15);
    assert_eq!(hit.s_last_off, 13);
}

#[test]
fn one_hit_extends_from_the_best_part_of_the_word() {
    let pssm = rows(10);
    let mut subject = [20u8; 10];
    subject[3..7].fill(1);
    // Word at q 2 / s 2: W A A -> sum -4 (reset), 5, 10.
    let hit = pssm_extend_one_hit(&pssm, &subject, 2, 2, 16, 3);
    assert_eq!(
        hit.ungapped_data,
        BlastpUngappedData {
            q_start: 3,
            s_start: 3,
            length: 4,
            score: 20,
        }
    );
}
