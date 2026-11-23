//! The one-hit ungapped extension (`-window_size 0`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c

use super::query_batch::{Again, InBatch, Pass, WalkEnd};
use super::two_hit::{extend_left_ungapped, extend_right_ungapped};
use crate::algorithm::common::translated_chunks::ChunkLetters;
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1037-1085
/// `s_BlastAaExtendOneHit` with the scores of `sc` (a matrix, not a PSSM) on a context's letters
/// `q_seq` and a chunk's `s_seq`, from the word at `q_off` and `s_off`; `None` if the word does not
/// fit in both. `init_hit_width` is taken from the last reset of the running sum, so the length is 0
/// or negative when the sum is reset after its best part. The left and right walks are the two-hit
/// extension's (`s_BlastAaExtendLeft` and `s_BlastAaExtendRight`): they stop at the subject
/// slice's ends, and at the context's ends they go on into the query batch `qb` as NCBI's walks do
/// (query_batch.rs). The two passes as the one-hit scan runs them (`extend_hit_one_hit_alone`,
/// then `extend_hit_one_hit_resume` where the first stops; `word_finder::scan_chunk_one_hit`); the
/// tests' entry.
#[cfg(test)]
pub(crate) fn extend_hit_one_hit<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    x_drop: i32,
) -> Option<OneHit> {
    match extend_hit_one_hit_alone(sc, q_seq, s_seq, q_off, s_off, x_drop)? {
        OneHitAlone::Hit(hit) => Some(hit),
        OneHitAlone::Again(again) => {
            extend_hit_one_hit_resume(sc, q_seq, qb, s_seq, q_off, s_off, x_drop, again)
        }
    }
}

/// What the one-hit extension on a context alone gives (`extend_hit_one_hit_alone`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum OneHitAlone {
    /// The extension's result, which NCBI's walks on the query batch also give.
    Hit(OneHit),
    /// A walk ended at the context's end where NCBI's could go on past it into the query batch
    /// (the NULLB not breaking it, `InBatch::left_at_bound`, `InBatch::right_at_bound`): the search
    /// goes on from there on the batch (`extend_hit_one_hit_resume`).
    Again(Again),
}

/// The first pass of `extend_hit_one_hit`, on the context alone (`InBatch::ALONE`: the walks stop
/// at its ends), up to a walk that ends at the context's end where NCBI's could go on into the
/// query batch (`OneHitAlone::Again`); elsewhere its result is NCBI's. The search's loop so keeps
/// no values of the batch through the walks.
#[inline(always)]
pub(crate) fn extend_hit_one_hit_alone<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    x_drop: i32,
) -> Option<OneHitAlone> {
    Some(
        match one_hit_walks(
            sc,
            q_seq,
            InBatch::ALONE,
            s_seq,
            q_off,
            s_off,
            x_drop,
            Pass::First,
        )? {
            Ok(hit) => OneHitAlone::Hit(hit),
            Err(again) => OneHitAlone::Again(again),
        },
    )
}

/// The second pass of `extend_hit_one_hit`: from where the first stopped (`again`), on the batch
/// `qb` (the word's offsets and scores are taken again; its walks are not).
pub(crate) fn extend_hit_one_hit_resume<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    x_drop: i32,
    again: Again,
) -> Option<OneHit> {
    match one_hit_walks(
        sc,
        q_seq,
        qb,
        s_seq,
        q_off,
        s_off,
        x_drop,
        Pass::Resume(again),
    )? {
        Ok(hit) => Some(hit),
        Err(_) => unreachable!("only the first pass stops at a walk"),
    }
}

/// The steps of `extend_hit_one_hit` with the walks of `qb` in `pass` (`Pass::First` or
/// `Pass::Resume`): `None` where the word does not fit, else the result, or where the first pass
/// stopped.
#[inline(always)]
fn one_hit_walks<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    x_drop: i32,
    pass: Pass,
) -> Option<Result<OneHit, Again>> {
    let word = sc.word();
    // The word's letters (`None` where it does not fit in both).
    let (Some(q_word), Some(s_word)) = (
        q_seq.get(q_off..q_off.checked_add(word)?),
        s_seq.get(s_off..s_off.checked_add(word)?),
    ) else {
        return None;
    };
    let (mut score, mut sum) = (0i32, 0i32);
    let (mut q_left_off, mut q_right_off, mut q_best_left_off) = (q_off, q_off + word, q_off);
    for (i, (&q, &s)) in q_word.iter().zip(s_word).enumerate() {
        sum += sc.score(q, s);
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
    let left = match pass {
        // NCBI's walk goes on past the context's start (which leaves its best and length).
        Pass::Resume(Again::Left(left)) => {
            qb.left_at_bound(
                sc,
                s_seq,
                q_left_off,
                s_left_off,
                left.score,
                left.max_score,
                x_drop,
            );
            left
        }
        Pass::Resume(Again::Right(left, _)) => left,
        Pass::First => {
            let (max_score, disp, again) =
                extend_left_ungapped(sc, q_seq, qb, s_seq, q_left_off, s_left_off, score, x_drop);
            let left = WalkEnd {
                max_score,
                disp,
                steps: 0,
                score: again.unwrap_or(max_score),
            };
            if let (Pass::First, Some(_)) = (pass, again) {
                return Some(Err(Again::Left(left)));
            }
            left
        }
    };
    let (left_score, left_disp) = (left.max_score, left.disp);
    // `last` is the walk's `i` at its end: `*s_last_off = s_off + i`.
    let (q_start, s_start) = (q_right_off + 1, s_right_off + 1);
    let (total_score, right_disp, last) = match pass {
        // NCBI's walk goes on past the context's end (which leaves its best and length).
        Pass::Resume(Again::Right(_, right)) => {
            let (j, _) = qb.right_at_bound(
                sc,
                s_seq,
                s_start + right.steps,
                q_start + right.steps,
                q_start > q_seq.len(),
                right.steps,
                right.score,
                right.max_score,
                x_drop,
            );
            (right.max_score, right.disp, j)
        }
        _ => {
            let (max_score, disp, j, again) =
                extend_right_ungapped(sc, q_seq, qb, s_seq, q_start, s_start, left_score, x_drop);
            if let (Pass::First, Some(score)) = (pass, again) {
                return Some(Err(Again::Right(
                    left,
                    WalkEnd {
                        max_score,
                        disp,
                        steps: j,
                        score,
                    },
                )));
            }
            (max_score, disp, j)
        }
    };
    Some(Ok(OneHit {
        q_start: (q_left_off - left_disp) as i32,
        s_start: (s_left_off - left_disp) as i32,
        length: left_disp as i32 + right_disp as i32 + init_hit_width,
        score: total_score,
        s_last_off: (s_right_off + 1 + last) as i32,
    }))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1020-1083, :897 and :841
/// `s_BlastAaExtendOneHit` with the scores of `sc` of a word (at `s_off`, below 0) that starts
/// before its chunk, on the letters around the chunk (`ChunkLetters`, whose type doc gives NCBI's
/// reads there): the word on the frame's letters, the left walk down to the chunk's first letter
/// (none where the word's best part starts at or before it), the right walk from the frame's
/// letters after the word up to the chunk's last letter. The walks go on into the query batch
/// `qb` past the context's ends as NCBI's do; the subject offsets are the chunk's (`Int4`).
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn extend_hit_one_hit_outside<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    chunk: ChunkLetters<'_>,
    q_off: usize,
    s_off: i32,
    x_drop: i32,
) -> OneHit {
    let word = sc.word();
    let (mut score, mut sum) = (0i32, 0i32);
    let (mut q_left_off, mut q_right_off, mut q_best_left_off) = (q_off, q_off + word, q_off);
    for i in 0..word {
        sum += sc.score(q_seq[q_off + i], chunk.at(s_off + i as i32));
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
    let s_left_off = q_left_off as i32 + s_off - q_off as i32;
    let s_right_off = q_right_off as i32 + s_off - q_off as i32;
    // `extend_left_ungapped` walks from the letter before the offsets it is given.
    let (left_score, left_disp) = if s_left_off > 0 {
        let (max_score, disp, _) = extend_left_ungapped(
            sc,
            q_seq,
            qb,
            chunk.from_chunk(),
            q_left_off,
            s_left_off as usize,
            score,
            x_drop,
        );
        (max_score, disp)
    } else {
        (score, 0)
    };
    let (total_score, right_disp, last, _) = extend_right_ungapped(
        sc,
        q_seq,
        qb,
        chunk.to_chunk_end(),
        q_right_off + 1,
        (chunk.lead as i32 + s_right_off + 1) as usize,
        left_score,
        x_drop,
    );
    OneHit {
        q_start: (q_left_off - left_disp) as i32,
        s_start: s_left_off - left_disp as i32,
        length: left_disp as i32 + right_disp as i32 + init_hit_width,
        score: total_score,
        s_last_off: s_right_off + 1 + last as i32,
    }
}
