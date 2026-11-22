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
    // NCBI's bound past the rows of this PSSM (RPS-BLAST's concatenated length above its rows,
    // `PssmRows::past_end`): the loop above it scores the next row for the next letter as well.
    if i == n {
        if let Some(past) = rows.past_end(q_off + i, subject.get((s_off + i) as usize).copied()) {
            score += past;
            if score > maxscore {
                maxscore = score;
                best_i = i;
            }
            if !(score <= 0 || (maxscore - score) >= dropoff) {
                i += 1;
                rows.past_end_continues(q_off + i, s_off + i < subject.len() as i32);
            }
        }
    }
    (maxscore, best_i + 1, s_off + i)
}

/// `s_BlastPSSMExtendLeft`: (the best score, its length).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:986-1017
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1036-1085
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
#[path = "ungapped_tests.rs"]
mod tests;
