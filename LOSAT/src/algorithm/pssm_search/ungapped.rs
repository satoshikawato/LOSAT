//! The ungapped extensions of a PSSM search: NCBI's `s_BlastPSSMExtendLeft`,
//! `s_BlastPSSMExtendRight` and the `use_pssm` branches of `s_BlastAaExtendOneHit` and
//! `s_BlastAaExtendTwoHit`, which score `matrix[query position][subject residue]`.
//!
//! The subject is `subject->sequence` (its letters without the sentinels, `subject.len()` is
//! `subject->length`); the query length is the PSSM's (`query->length` of the one query).

use super::rows::PssmRows;
use crate::algorithm::blastp::extension::{
    BlastpOneHitUngappedResult, BlastpTwoHitUngappedResult, BlastpUngappedData,
};

/// `s_BlastPSSMExtendRight`: (the best score, its length, `s_last_off`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:938-975
/// ```c
///     n = MIN(subject->length - s_off, query_size - q_off);
///     s = subject->sequence + s_off;
///
///     for (i = 0; i < n; i++) {
///         score += matrix[q_off + i][s[i]];
///
///         if (score > maxscore) {
///             maxscore = score;
///             best_i = i;
///         }
///
///         /* The comparison below is really >= and is different than the old
///            code (e.g., blast.c:BlastWordExtend_prelim). In the old code the
///            loop continued as long as sum > X (X being negative).  The loop
///            control here is different and we *break out* when the if statement
///            below is true. */
///         if (score <= 0 || (maxscore - score) >= dropoff)
///             break;
///     }
///
///     *length = best_i + 1;
///     *s_last_off = s_off + i;
///     return maxscore;
/// ```
#[inline]
fn pssm_extend_right<R: PssmRows>(
    rows: &R,
    subject: &[u8],
    s_off: i32,
    q_off: i32,
    dropoff: i32,
    mut maxscore: i32,
) -> (i32, i32, i32) {
    let n = (subject.len() as i32 - s_off).min(rows.len() as i32 - q_off);
    let mut score = maxscore;
    let mut best_i = -1i32;
    let mut i = 0i32;
    while i < n {
        score += rows.row((q_off + i) as usize)[subject[(s_off + i) as usize] as usize];
        if score > maxscore {
            maxscore = score;
            best_i = i;
        }
        if score <= 0 || (maxscore - score) >= dropoff {
            break;
        }
        i += 1;
    }
    (maxscore, best_i + 1, s_off + i)
}

/// `s_BlastPSSMExtendLeft`: (the best score, its length).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:986-1017
/// ```c
///     n = MIN(s_off, q_off);
///     best_i = n + 1;
///     s = subject->sequence + s_off - n;
///
///     for (i = n; i >= 0; i--) {
///         score += matrix[q_off - n + i][s[i]];
///
///         if (score > maxscore) {
///             maxscore = score;
///             best_i = i;
///         }
///         /* The comparison below is really >= and is different than the old
///            code (e.g., blast.c:BlastWordExtend_prelim). In the old code the
///            loop continued as long as sum > X (X being negative).  The loop
///            control here is different and we *break out* when the if statement
///            below is true. */
///         if ((maxscore - score) >= dropoff)
///             break;
///     }
///
///     *length = n - best_i + 1;
///     return maxscore;
/// ```
#[inline]
fn pssm_extend_left<R: PssmRows>(
    rows: &R,
    subject: &[u8],
    s_off: i32,
    q_off: i32,
    dropoff: i32,
    mut maxscore: i32,
) -> (i32, i32) {
    let n = s_off.min(q_off);
    let mut best_i = n + 1;
    let mut score = maxscore;
    let mut i = n;
    while i >= 0 {
        score += rows.row((q_off - n + i) as usize)[subject[(s_off - n + i) as usize] as usize];
        if score > maxscore {
            maxscore = score;
            best_i = i;
        }
        if (maxscore - score) >= dropoff {
            break;
        }
        i -= 1;
    }
    (maxscore, n - best_i + 1)
}

/// `s_BlastAaExtendOneHit` with `use_pssm`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1036-1083
/// ```c
///     for (i = 0; i < word_size; i++) {
///         if (use_pssm)
///             sum += matrix[q_off + i][s[s_off + i]];
///         else
///             sum += matrix[q[q_off + i]][s[s_off + i]];
///
///         if (sum > score) {
///             score = sum;
///             q_best_left_off = q_left_off;
///             q_right_off = q_off + i;
///         } else if (sum <= 0) {
///             sum = 0;
///             q_left_off = q_off + i + 1;
///         }
///     }
///
///     init_hit_width = q_right_off - q_left_off + 1;
///
///     q_left_off = q_best_left_off;
///
///     s_left_off = q_left_off + (s_off - q_off);
///     s_right_off = q_right_off + (s_off - q_off);
///
///     if (use_pssm) {
///         left_score = s_BlastPSSMExtendLeft(matrix, subject,
///                                            s_left_off - 1, q_left_off - 1,
///                                            dropoff, &left_disp, score);
///
///         total_score = s_BlastPSSMExtendRight(matrix, subject, query->length,
///                                              s_right_off + 1, q_right_off + 1,
///                                              dropoff, &right_disp, left_score,
///                                              s_last_off);
///     } else {
///     ...
///     *hsp_q = q_left_off - left_disp;
///     *hsp_s = s_left_off - left_disp;
///     *hsp_len = left_disp + right_disp + init_hit_width;
///
///     return total_score;
/// ```
pub(crate) fn pssm_extend_one_hit<R: PssmRows>(
    rows: &R,
    subject: &[u8],
    s_off: i32,
    q_off: i32,
    dropoff: i32,
    word_size: i32,
) -> BlastpOneHitUngappedResult {
    let mut score = 0i32;
    let mut sum = 0i32;
    let mut q_left_off = q_off;
    let mut q_right_off = q_off + word_size;
    let mut q_best_left_off = q_off;
    for i in 0..word_size {
        sum += rows.row((q_off + i) as usize)[subject[(s_off + i) as usize] as usize];
        if sum > score {
            score = sum;
            q_best_left_off = q_left_off;
            q_right_off = q_off + i;
        } else if sum <= 0 {
            sum = 0;
            q_left_off = q_off + i + 1;
        }
    }
    let init_hit_width = q_right_off - q_left_off + 1;
    let q_left_off = q_best_left_off;
    let s_left_off = q_left_off + (s_off - q_off);
    let s_right_off = q_right_off + (s_off - q_off);
    let (left_score, left_disp) = pssm_extend_left(
        rows,
        subject,
        s_left_off - 1,
        q_left_off - 1,
        dropoff,
        score,
    );
    let (total_score, right_disp, s_last_off) = pssm_extend_right(
        rows,
        subject,
        s_right_off + 1,
        q_right_off + 1,
        dropoff,
        left_score,
    );
    BlastpOneHitUngappedResult {
        ungapped_data: BlastpUngappedData {
            q_start: q_left_off - left_disp,
            s_start: s_left_off - left_disp,
            length: left_disp + right_disp + init_hit_width,
            score: total_score,
        },
        s_last_off,
    }
}

/// `s_BlastAaExtendTwoHit` with `use_pssm`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1104-1157
/// ```c
///     for (i = 0; i < word_size; i++) {
///         if (use_pssm)
///             score += matrix[q_right_off + i][s[s_right_off + i]];
///         else
///             score += matrix[q[q_right_off + i]][s[s_right_off + i]];
///
///         if (score > left_score) {
///             left_score = score;
///             right_d = i + 1;    /* Position is one beyond the end of the
///                                    word. */
///         }
///     }
///     q_right_off += right_d;
///     s_right_off += right_d;
///
///     right_d = 0;
///     *right_extend = FALSE;
///     *s_last_off = s_right_off;
///
///     /* first, try to extend left, from the second hit to the first hit. */
///     if (use_pssm)
///         left_score = s_BlastPSSMExtendLeft(matrix, subject,
///                                            s_right_off - 1, q_right_off - 1,
///                                            dropoff, &left_d, 0);
///     ...
///     /* Extend to the right only if left extension reached the first hit. */
///     if (left_d >= (s_right_off - s_left_off)) {
///
///         *right_extend = TRUE;
///         if (use_pssm)
///             right_score = s_BlastPSSMExtendRight(matrix, subject, query->length,
///                                                  s_right_off, q_right_off,
///                                                  dropoff, &right_d, left_score,
///                                                  s_last_off);
///     ...
///     *hsp_q = q_right_off - left_d;
///     *hsp_s = s_right_off - left_d;
///     *hsp_len = left_d + right_d;
///     return MAX(left_score, right_score);
/// ```
pub(crate) fn pssm_extend_two_hit<R: PssmRows>(
    rows: &R,
    subject: &[u8],
    s_left_off: i32,
    s_right_off: i32,
    q_right_off: i32,
    dropoff: i32,
    word_size: i32,
) -> BlastpTwoHitUngappedResult {
    let mut left_score = 0i32;
    let mut right_d = 0i32;
    let mut score = 0i32;
    for i in 0..word_size {
        score += rows.row((q_right_off + i) as usize)[subject[(s_right_off + i) as usize] as usize];
        if score > left_score {
            left_score = score;
            right_d = i + 1;
        }
    }
    let q_right_off = q_right_off + right_d;
    let s_right_off = s_right_off + right_d;
    let mut right_extend = false;
    let mut s_last_off = s_right_off;
    let (left_score, left_d) =
        pssm_extend_left(rows, subject, s_right_off - 1, q_right_off - 1, dropoff, 0);
    let mut right_score = 0i32;
    let mut right_d = 0i32;
    if left_d >= (s_right_off - s_left_off) {
        right_extend = true;
        let (score, length, last) =
            pssm_extend_right(rows, subject, s_right_off, q_right_off, dropoff, left_score);
        right_score = score;
        right_d = length;
        s_last_off = last;
    }
    BlastpTwoHitUngappedResult {
        ungapped_data: BlastpUngappedData {
            q_start: q_right_off - left_d,
            s_start: s_right_off - left_d,
            length: left_d + right_d,
            score: left_score.max(right_score),
        },
        right_extend,
        s_last_off,
    }
}

#[cfg(test)]
mod tests {
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
}
