//! BLASTP ungapped extension helpers.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c

use crate::config::ScoringMatrix;
use crate::utils::matrix::protein_score;

#[inline(always)]
fn residue_score(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    q_off: i32,
    s_off: i32,
) -> i32 {
    protein_score(matrix, query[q_off as usize], subject[s_off as usize])
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:804-835
// ```c
// static Int4 s_BlastAaExtendRight(Int4 ** matrix,
//                         const BLAST_SequenceBlk * subject,
//                         const BLAST_SequenceBlk * query,
//                         Int4 s_off, Int4 q_off, Int4 dropoff,
//                         Int4 * length, Int4 maxscore, Int4 * s_last_off)
// {
//     Int4 i, n, best_i = -1;
//     Int4 score = maxscore;
//     ...
//     for (i = 0; i < n; i++) {
//         score += matrix[q[i]][s[i]];
//         if (score > maxscore) { ... }
//         if (score <= 0 || (maxscore - score) >= dropoff)
//             break;
//     }
//     *length = best_i + 1;
//     *s_last_off = s_off + i;
//     return maxscore;
// }
// ```
fn extend_right(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    q_off: i32,
    s_off: i32,
    dropoff: i32,
    maxscore: i32,
) -> (i32, i32, i32) {
    let n = (subject.len() as i32 - s_off).min(query.len() as i32 - q_off);
    let mut best_i = -1i32;
    let mut score = maxscore;
    let mut best_score = maxscore;
    let mut i = 0i32;

    while i < n {
        score += residue_score(matrix, query, subject, q_off + i, s_off + i);

        if score > best_score {
            best_score = score;
            best_i = i;
        }

        if score <= 0 || (best_score - score) >= dropoff {
            break;
        }
        i += 1;
    }

    (best_score, best_i + 1, s_off + i)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:853-883
// ```c
// static Int4 s_BlastAaExtendLeft(Int4 ** matrix,
//                        const BLAST_SequenceBlk * subject,
//                        const BLAST_SequenceBlk * query,
//                        Int4 s_off, Int4 q_off, Int4 dropoff,
//                        Int4 * length, Int4 maxscore)
// {
//     Int4 i, n, best_i;
//     Int4 score = maxscore;
//     ...
//     n = MIN(s_off, q_off);
//     best_i = n + 1;
//     ...
//     for (i = n; i >= 0; i--) {
//         score += matrix[q[i]][s[i]];
//         if (score > maxscore) { ... }
//         if ((maxscore - score) >= dropoff)
//             break;
//     }
//     *length = n - best_i + 1;
//     return maxscore;
// }
// ```
fn extend_left(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    q_off: i32,
    s_off: i32,
    dropoff: i32,
    maxscore: i32,
) -> (i32, i32) {
    let n = s_off.min(q_off);
    let mut best_i = n + 1;
    let mut score = maxscore;
    let mut best_score = maxscore;

    if n >= 0 {
        let q_base = q_off - n;
        let s_base = s_off - n;
        let mut i = n;

        loop {
            score += residue_score(matrix, query, subject, q_base + i, s_base + i);

            if score > best_score {
                best_score = score;
                best_i = i;
            }

            if (best_score - score) >= dropoff {
                break;
            }
            if i == 0 {
                break;
            }
            i -= 1;
        }
    }

    (best_score, n - best_i + 1)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1078-1127
// ```c
// static Int4 s_BlastAaExtendOneHit(...)
// {
//     ...
//     for (i = 0; i < word_size; i++) {
//         sum += matrix[q[q_off + i]][s[s_off + i]];
//         if (sum > score) {
//             score = sum;
//             q_best_left_off = q_left_off;
//             q_right_off = q_off + i;
//         } else if (sum <= 0) {
//             sum = 0;
//             q_left_off = q_off + i + 1;
//         }
//     }
//     ...
//     left_score = s_BlastAaExtendLeft(..., s_left_off - 1, q_left_off - 1, ...);
//     total_score = s_BlastAaExtendRight(..., s_right_off + 1, q_right_off + 1, ...);
//     *hsp_q = q_left_off - left_disp;
//     *hsp_s = s_left_off - left_disp;
//     *hsp_len = left_disp + right_disp + init_hit_width;
//     return total_score;
// }
// ```
pub fn extend_one_hit(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<(usize, usize, usize, usize, i32, usize)> {
    if q_off + word_size > query.len() || s_off + word_size > subject.len() {
        return None;
    }

    let q_off = q_off as i32;
    let s_off = s_off as i32;
    let word_size = word_size as i32;

    let mut score = 0i32;
    let mut sum = 0i32;
    let mut q_left_off = q_off;
    let mut q_right_off = q_off + word_size;
    let mut q_best_left_off = q_off;

    for i in 0..word_size {
        sum += residue_score(matrix, query, subject, q_off + i, s_off + i);

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
    q_left_off = q_best_left_off;

    let s_left_off = q_left_off + (s_off - q_off);
    let s_right_off = q_right_off + (s_off - q_off);

    let (left_score, left_disp) = extend_left(
        matrix,
        query,
        subject,
        q_left_off - 1,
        s_left_off - 1,
        dropoff,
        score,
    );
    let (total_score, right_disp, s_last_off) = extend_right(
        matrix,
        query,
        subject,
        q_right_off + 1,
        s_right_off + 1,
        dropoff,
        left_score,
    );

    let q_start = (q_left_off - left_disp) as usize;
    let q_end = (q_start as i32 + left_disp + right_disp + init_hit_width) as usize;
    let s_start = (s_left_off - left_disp) as usize;
    let s_end = (s_start as i32 + left_disp + right_disp + init_hit_width) as usize;

    Some((
        q_start,
        q_end,
        s_start,
        s_end,
        total_score,
        s_last_off as usize,
    ))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1129-1177
// ```c
// static Int4 s_BlastAaExtendTwoHit(...)
// {
//     ...
//     for (i = 0; i < word_size; i++) {
//         score += matrix[q[q_right_off + i]][s[s_right_off + i]];
//         if (score > left_score) {
//             left_score = score;
//             right_d = i + 1;
//         }
//     }
//     q_right_off += right_d;
//     s_right_off += right_d;
//     ...
//     left_score = s_BlastAaExtendLeft(..., s_right_off - 1, q_right_off - 1, ...);
//     if (left_d >= (s_right_off - s_left_off)) {
//         *right_extend = TRUE;
//         right_score = s_BlastAaExtendRight(..., s_right_off, q_right_off, ...);
//     }
//     *hsp_q = q_right_off - left_d;
//     *hsp_s = s_right_off - left_d;
//     *hsp_len = left_d + right_d;
//     return MAX(left_score, right_score);
// }
// ```
pub fn extend_two_hit(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<(usize, usize, usize, usize, i32, bool, usize)> {
    if q_right_off + word_size > query.len() || s_right_off + word_size > subject.len() {
        return None;
    }

    let mut q_right_off = q_right_off as i32;
    let s_left_off = s_left_off as i32;
    let mut s_right_off = s_right_off as i32;
    let word_size = word_size as i32;

    let mut left_d = 0i32;
    let mut right_d = 0i32;
    let mut left_score = 0i32;
    let mut right_score = 0i32;
    let mut score = 0i32;

    for i in 0..word_size {
        score += residue_score(matrix, query, subject, q_right_off + i, s_right_off + i);
        if score > left_score {
            left_score = score;
            right_d = i + 1;
        }
    }
    q_right_off += right_d;
    s_right_off += right_d;
    right_d = 0;

    let mut right_extend = false;
    let mut s_last_off = s_right_off;

    (left_score, left_d) = extend_left(
        matrix,
        query,
        subject,
        q_right_off - 1,
        s_right_off - 1,
        dropoff,
        0,
    );

    if left_d >= (s_right_off - s_left_off) {
        right_extend = true;
        (right_score, right_d, s_last_off) = extend_right(
            matrix,
            query,
            subject,
            q_right_off,
            s_right_off,
            dropoff,
            left_score,
        );
    }

    let q_start = (q_right_off - left_d) as usize;
    let q_end = (q_right_off + right_d) as usize;
    let s_start = (s_right_off - left_d) as usize;
    let s_end = (s_right_off + right_d) as usize;
    let best_score = left_score.max(right_score);

    Some((
        q_start,
        q_end,
        s_start,
        s_end,
        best_score,
        right_extend,
        s_last_off as usize,
    ))
}
