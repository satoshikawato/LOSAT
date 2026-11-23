//! BLASTP ungapped extension helpers.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c

use super::x_ungapped_fast;
pub(crate) use crate::algorithm::common::translated_chunks::ChunkLetters;
use crate::config::ScoringMatrix;
use crate::utils::matrix::{blosum62_score_ncbistdaa_direct, protein_score};

#[inline(always)]
fn residue_score<const BLOSUM62: bool>(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    q_off: i32,
    s_off: i32,
) -> i32 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883
    let query_residue = query[q_off as usize];
    let subject_residue = subject[s_off as usize];
    if BLOSUM62 {
        blosum62_score_ncbistdaa_direct(query_residue, subject_residue)
    } else {
        protein_score(matrix, query_residue, subject_residue)
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:142-154
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct BlastpUngappedData {
    pub(crate) q_start: i32,
    pub(crate) s_start: i32,
    pub(crate) length: i32,
    pub(crate) score: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1078-1127
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct BlastpOneHitUngappedResult {
    pub(crate) ungapped_data: BlastpUngappedData,
    pub(crate) s_last_off: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1124-1158
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct BlastpTwoHitUngappedResult {
    pub(crate) ungapped_data: BlastpUngappedData,
    pub(crate) right_extend: bool,
    pub(crate) s_last_off: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-865
fn extend_right<const BLOSUM62: bool>(
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
        score += residue_score::<BLOSUM62>(matrix, query, subject, q_off + i, s_off + i);

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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:853-919
fn extend_left<const BLOSUM62: bool>(
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
            score += residue_score::<BLOSUM62>(matrix, query, subject, q_base + i, s_base + i);

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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1038-1127
#[allow(dead_code)]
pub fn extend_one_hit(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpOneHitUngappedResult> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1039-1044
    if matrix == ScoringMatrix::Blosum62 {
        extend_one_hit_impl::<true>(matrix, query, subject, q_off, s_off, dropoff, word_size)
    } else {
        extend_one_hit_impl::<false>(matrix, query, subject, q_off, s_off, dropoff, word_size)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1039-1044
/// [`extend_one_hit`] of a matrix other than BLOSUM62 (BLASTP's other matrices, item SC2b): the
/// reference walk with `protein_score`, without the dispatch, so that BLOSUM62's instance keeps the
/// one caller (`extend_one_hit_blosum62`'s reference) that inlines it.
pub(crate) fn extend_one_hit_of(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpOneHitUngappedResult> {
    debug_assert_ne!(matrix, ScoringMatrix::Blosum62);
    extend_one_hit_impl::<false>(matrix, query, subject, q_off, s_off, dropoff, word_size)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1039-1044
#[inline]
pub(crate) fn extend_one_hit_blosum62(
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpOneHitUngappedResult> {
    let reference = || {
        extend_one_hit_impl::<true>(
            ScoringMatrix::Blosum62,
            query,
            subject,
            q_off,
            s_off,
            dropoff,
            word_size,
        )
    };
    // x_ungapped_fast.rs: the same result, the X-drop walks with a byte score table and no bounds
    // checks; the reference runs where the fast walks do not cover the call.
    x_ungapped_fast::one_hit(query, subject, q_off, s_off, dropoff, word_size)
        .unwrap_or_else(reference)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:785-790,1040-1044
/// `s_BlastAaExtendOneHit` with the score matrix `matrix` (not a PSSM), chosen at each call:
/// BLOSUM62's walks (`extend_one_hit_blosum62`) or the reference walk of another matrix
/// (`extend_one_hit_of`). TBLASTN's one-hit word finder (`-window_size 0`) calls it with its matrix
/// (TBLASTX, below BLASTP in the layering, has its own on its two-hit walks).
#[inline]
pub(crate) fn extend_one_hit_for(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpOneHitUngappedResult> {
    if matrix == ScoringMatrix::Blosum62 {
        extend_one_hit_blosum62(query, subject, q_off, s_off, dropoff, word_size)
    } else {
        extend_one_hit_of(matrix, query, subject, q_off, s_off, dropoff, word_size)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1020-1083
fn extend_one_hit_impl<const BLOSUM62: bool>(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpOneHitUngappedResult> {
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
        sum += residue_score::<BLOSUM62>(matrix, query, subject, q_off + i, s_off + i);

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

    let (left_score, left_disp) = extend_left::<BLOSUM62>(
        matrix,
        query,
        subject,
        q_left_off - 1,
        s_left_off - 1,
        dropoff,
        score,
    );
    let (total_score, right_disp, s_last_off) = extend_right::<BLOSUM62>(
        matrix,
        query,
        subject,
        q_right_off + 1,
        s_right_off + 1,
        dropoff,
        left_score,
    );

    Some(BlastpOneHitUngappedResult {
        ungapped_data: BlastpUngappedData {
            q_start: q_left_off - left_disp,
            s_start: s_left_off - left_disp,
            length: left_disp + right_disp + init_hit_width,
            score: total_score,
        },
        s_last_off,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1108-1158
#[allow(dead_code)]
pub fn extend_two_hit(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpTwoHitUngappedResult> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1109-1116
    if matrix == ScoringMatrix::Blosum62 {
        extend_two_hit_impl::<true>(
            matrix,
            query,
            subject,
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            word_size,
        )
    } else {
        extend_two_hit_impl::<false>(
            matrix,
            query,
            subject,
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            word_size,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1109-1116
#[inline]
pub(crate) fn extend_two_hit_blosum62(
    query: &[u8],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpTwoHitUngappedResult> {
    let reference = || {
        extend_two_hit_impl::<true>(
            ScoringMatrix::Blosum62,
            query,
            subject,
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            word_size,
        )
    };
    // x_ungapped_fast.rs: the same result, the X-drop walks with a byte score table and no bounds
    // checks; the reference runs where the fast walks do not cover the call.
    x_ungapped_fast::two_hit(
        query,
        subject,
        s_left_off,
        s_right_off,
        q_right_off,
        dropoff,
        word_size,
    )
    .unwrap_or_else(reference)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1109-1116
/// [`extend_two_hit`] of a matrix other than BLOSUM62 (item SC2b), with the byte-table walks of
/// the matrix (`x_ungapped_fast::two_hit_of`; `of.1` is the matrix's `matrix_table`): the same
/// result. The reference runs where the matrix has no table or the walks do not cover the call.
#[inline]
pub(crate) fn extend_two_hit_table(
    of: (ScoringMatrix, Option<&'static x_ungapped_fast::MatrixTable>),
    query: &[u8],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpTwoHitUngappedResult> {
    let (matrix, table) = of;
    table
        .and_then(|table| {
            x_ungapped_fast::two_hit_of(
                table,
                query,
                subject,
                s_left_off,
                s_right_off,
                q_right_off,
                dropoff,
                word_size,
            )
        })
        .unwrap_or_else(|| {
            extend_two_hit(
                matrix,
                query,
                subject,
                s_left_off,
                s_right_off,
                q_right_off,
                dropoff,
                word_size,
            )
        })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155
fn extend_two_hit_impl<const BLOSUM62: bool>(
    matrix: ScoringMatrix,
    query: &[u8],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpTwoHitUngappedResult> {
    if q_right_off + word_size > query.len() || s_right_off + word_size > subject.len() {
        return None;
    }

    let mut q_right_off = q_right_off as i32;
    let s_left_off = s_left_off as i32;
    let mut s_right_off = s_right_off as i32;
    let word_size = word_size as i32;

    let left_d: i32;
    let mut right_d = 0i32;
    let mut left_score = 0i32;
    let mut right_score = 0i32;
    let mut score = 0i32;

    for i in 0..word_size {
        score +=
            residue_score::<BLOSUM62>(matrix, query, subject, q_right_off + i, s_right_off + i);
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

    (left_score, left_d) = extend_left::<BLOSUM62>(
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
        (right_score, right_d, s_last_off) = extend_right::<BLOSUM62>(
            matrix,
            query,
            subject,
            q_right_off,
            s_right_off,
            dropoff,
            left_score,
        );
    }

    Some(BlastpTwoHitUngappedResult {
        ungapped_data: BlastpUngappedData {
            q_start: q_right_off - left_d,
            s_start: s_right_off - left_d,
            length: left_d + right_d,
            score: left_score.max(right_score),
        },
        right_extend,
        s_last_off,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 (the two-hit
// extension), :897 (the left walk counts the letters down to the subject's first one), :841 (the
// right walk up to the subject's length) and :1020-1083 (the one-hit extension), on the letters
// around a chunk as NCBI's extensions read them there (`ChunkLetters`).
/// `s_BlastAaExtendTwoHit` with the matrix `matrix` of a hit whose first hit (`s_left_off`) lies
/// before its chunk, on the letters around the chunk (`ChunkLetters`).
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn extend_two_hit_outside(
    matrix: ScoringMatrix,
    query: &[u8],
    chunk: ChunkLetters<'_>,
    s_left_off: i32,
    s_right_off: i32,
    q_right_off: i32,
    dropoff: i32,
    word_size: i32,
) -> BlastpTwoHitUngappedResult {
    if matrix == ScoringMatrix::Blosum62 {
        extend_two_hit_outside_impl::<true>(
            matrix,
            query,
            chunk,
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            word_size,
        )
    } else {
        extend_two_hit_outside_impl::<false>(
            matrix,
            query,
            chunk,
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            word_size,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn extend_two_hit_outside_impl<const BLOSUM62: bool>(
    matrix: ScoringMatrix,
    query: &[u8],
    chunk: ChunkLetters<'_>,
    s_left_off: i32,
    s_right_off: i32,
    q_right_off: i32,
    dropoff: i32,
    word_size: i32,
) -> BlastpTwoHitUngappedResult {
    let mut right_d = 0i32;
    let mut left_score = 0i32;
    let mut score = 0i32;
    for i in 0..word_size {
        score += residue_score::<BLOSUM62>(
            matrix,
            query,
            chunk.letters,
            q_right_off + i,
            chunk.lead as i32 + s_right_off + i,
        );
        if score > left_score {
            left_score = score;
            right_d = i + 1;
        }
    }
    let q_right_off = q_right_off + right_d;
    let s_right_off = s_right_off + right_d;
    // From the letter before the best end down to the chunk's first letter: none where that end is
    // at or before the chunk's start (the walk's count is then below 0).
    let (left_score, left_d) = extend_left::<BLOSUM62>(
        matrix,
        query,
        chunk.chunk_onward(),
        q_right_off - 1,
        s_right_off - 1,
        dropoff,
        0,
    );
    let mut right_score = 0i32;
    let mut right_d = 0i32;
    let mut right_extend = false;
    let mut s_last_off = s_right_off;
    if left_d >= (s_right_off - s_left_off) {
        right_extend = true;
        let lead = chunk.lead as i32;
        let (best, length, last) = extend_right::<BLOSUM62>(
            matrix,
            query,
            chunk.to_chunk_end(),
            q_right_off,
            lead + s_right_off,
            dropoff,
            left_score,
        );
        (right_score, right_d, s_last_off) = (best, length, last - lead);
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

/// `s_BlastAaExtendOneHit` with the matrix `matrix` of a word (at `s_off`) that starts before its
/// chunk, on the letters around the chunk (`ChunkLetters`).
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn extend_one_hit_outside(
    matrix: ScoringMatrix,
    query: &[u8],
    chunk: ChunkLetters<'_>,
    q_off: i32,
    s_off: i32,
    dropoff: i32,
    word_size: i32,
) -> BlastpOneHitUngappedResult {
    if matrix == ScoringMatrix::Blosum62 {
        extend_one_hit_outside_impl::<true>(matrix, query, chunk, q_off, s_off, dropoff, word_size)
    } else {
        extend_one_hit_outside_impl::<false>(matrix, query, chunk, q_off, s_off, dropoff, word_size)
    }
}

fn extend_one_hit_outside_impl<const BLOSUM62: bool>(
    matrix: ScoringMatrix,
    query: &[u8],
    chunk: ChunkLetters<'_>,
    q_off: i32,
    s_off: i32,
    dropoff: i32,
    word_size: i32,
) -> BlastpOneHitUngappedResult {
    let lead = chunk.lead as i32;
    let mut score = 0i32;
    let mut sum = 0i32;
    let mut q_left_off = q_off;
    let mut q_right_off = q_off + word_size;
    let mut q_best_left_off = q_off;
    for i in 0..word_size {
        sum += residue_score::<BLOSUM62>(matrix, query, chunk.letters, q_off + i, lead + s_off + i);
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
    let (left_score, left_disp) = extend_left::<BLOSUM62>(
        matrix,
        query,
        chunk.chunk_onward(),
        q_left_off - 1,
        s_left_off - 1,
        dropoff,
        score,
    );
    let (total_score, right_disp, s_last_off) = extend_right::<BLOSUM62>(
        matrix,
        query,
        chunk.to_chunk_end(),
        q_right_off + 1,
        lead + s_right_off + 1,
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
        s_last_off: s_last_off - lead,
    }
}

#[cfg(test)]
#[path = "extension_tests.rs"]
mod tests;
