//! The one-hit ungapped extension (`-window_size 0`).
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c

use super::query_batch::InBatch;
use super::two_hit::{extend_left_ungapped, extend_right_ungapped};
use crate::algorithm::tblastx::scoring::TbxScore;

/// What `s_BlastAaExtendOneHit` returns and writes: `*hsp_q`, `*hsp_s`, `*hsp_len` (which can be 0 or
/// negative), the score, and `*s_last_off`, in the coordinates of the slices it was given.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct OneHit {
    pub(crate) q_start: i32,
    pub(crate) s_start: i32,
    pub(crate) length: i32,
    pub(crate) score: i32,
    pub(crate) s_last_off: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1037-1083
// ```c
//     for (i = 0; i < word_size; i++) {
//         if (use_pssm)
//             sum += matrix[q_off + i][s[s_off + i]];
//         else
//             sum += matrix[q[q_off + i]][s[s_off + i]];
//
//         if (sum > score) {
//             score = sum;
//             q_best_left_off = q_left_off;
//             q_right_off = q_off + i;
//         } else if (sum <= 0) {
//             sum = 0;
//             q_left_off = q_off + i + 1;
//         }
//     }
//
//     init_hit_width = q_right_off - q_left_off + 1;
//
//     q_left_off = q_best_left_off;
//
//     s_left_off = q_left_off + (s_off - q_off);
//     s_right_off = q_right_off + (s_off - q_off);
// ...
//         left_score = s_BlastAaExtendLeft(matrix, subject, query,
//                                          s_left_off - 1, q_left_off - 1,
//                                          dropoff, &left_disp, score);
//
//         total_score = s_BlastAaExtendRight(matrix, subject, query,
//                                            s_right_off + 1, q_right_off + 1,
//                                            dropoff, &right_disp, left_score,
//                                            s_last_off);
//
//     *hsp_q = q_left_off - left_disp;
//     *hsp_s = s_left_off - left_disp;
//     *hsp_len = left_disp + right_disp + init_hit_width;
//
//     return total_score;
// ```
/// `s_BlastAaExtendOneHit` with the scores of `sc` (a matrix, not a PSSM) on a context's letters
/// `q_seq` and a chunk's `s_seq`, from the word at `q_off` and `s_off`; `None` if the word does not
/// fit in both. `init_hit_width` is taken from the last reset of the running sum, so the length is 0
/// or negative when the sum is reset after its best part. The left and right walks are the two-hit
/// extension's (`s_BlastAaExtendLeft` and `s_BlastAaExtendRight`): they stop at the subject
/// slice's ends, and at the context's ends they go on into the query batch `qb` as NCBI's walks do
/// (query_batch.rs).
pub(crate) fn extend_hit_one_hit<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    x_drop: i32,
) -> Option<OneHit> {
    // Not the first pass: the walks give their result (`OneHitAlone::Again` is the first pass's).
    match one_hit_walks(sc, q_seq, qb, s_seq, q_off, s_off, x_drop, false)? {
        OneHitAlone::Hit(hit) => Some(hit),
        OneHitAlone::Again => None,
    }
}

/// What the one-hit extension on a context alone gives (`extend_hit_one_hit_alone`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum OneHitAlone {
    /// The extension's result, which NCBI's walks on the query batch also give.
    Hit(OneHit),
    /// A walk ended at the context's end where NCBI's would go on past it into the query batch
    /// (the NULLB not breaking it, `InBatch::left_at_bound`, `InBatch::right_at_bound`): the search
    /// runs the extension again on the batch (`extend_hit_one_hit`).
    Again,
}

/// `extend_hit_one_hit` on the context alone (`InBatch::ALONE`: the walks stop at its ends), up to
/// a walk that ends at the context's end where NCBI's would go on into the query batch
/// (`OneHitAlone::Again`); elsewhere the two give the same result. The search's loop so keeps no
/// values of the batch through the walks.
pub(crate) fn extend_hit_one_hit_alone<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    x_drop: i32,
) -> Option<OneHitAlone> {
    one_hit_walks(sc, q_seq, InBatch::ALONE, s_seq, q_off, s_off, x_drop, true)
}

/// The steps of `extend_hit_one_hit` with the walks of `qb`. `first`: the search's first pass on
/// the context alone, which ends at a walk that asks for the batch (`OneHitAlone::Again`; no other
/// pass gives it).
#[inline(always)]
fn one_hit_walks<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    x_drop: i32,
    first: bool,
) -> Option<OneHitAlone> {
    let word = sc.word();
    if q_off + word > q_seq.len() || s_off + word > s_seq.len() {
        return None;
    }
    let (mut score, mut sum) = (0i32, 0i32);
    let (mut q_left_off, mut q_right_off, mut q_best_left_off) = (q_off, q_off + word, q_off);
    for i in 0..word {
        sum += sc.score(q_seq[q_off + i], s_seq[s_off + i]);
        if sum > score {
            score = sum;
            q_best_left_off = q_left_off;
            q_right_off = q_off + i;
        } else if sum <= 0 {
            sum = 0;
            q_left_off = q_off + i + 1;
        }
    }
    let init_hit_width = q_right_off as i32 - q_left_off as i32 + 1;
    let q_left_off = q_best_left_off;
    let s_left_off = q_left_off + s_off - q_off;
    let s_right_off = q_right_off + s_off - q_off;
    // `extend_left_ungapped` walks from the letter before the offsets it is given.
    let (left_score, left_disp, left_again) =
        extend_left_ungapped(sc, q_seq, qb, s_seq, q_left_off, s_left_off, score, x_drop);
    if first && left_again {
        return Some(OneHitAlone::Again);
    }
    // `last` is the walk's `i` at its end: `*s_last_off = s_off + i`.
    let (total_score, right_disp, last, right_again) = extend_right_ungapped(
        sc,
        q_seq,
        qb,
        s_seq,
        q_right_off + 1,
        s_right_off + 1,
        left_score,
        x_drop,
    );
    if first && right_again {
        return Some(OneHitAlone::Again);
    }
    Some(OneHitAlone::Hit(OneHit {
        q_start: (q_left_off - left_disp) as i32,
        s_start: (s_left_off - left_disp) as i32,
        length: left_disp as i32 + right_disp as i32 + init_hit_width,
        score: total_score,
        s_last_off: (s_right_off + 1 + last) as i32,
    }))
}
