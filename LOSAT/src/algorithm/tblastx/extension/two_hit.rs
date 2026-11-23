//! Two-hit ungapped extension algorithm
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c

use crate::algorithm::tblastx::scoring::{Blosum62Score, TbxScore};
// Where the context's letters lie in the query batch: NCBI's walks go on past the context's ends
// (query_batch.rs).
use super::query_batch::InBatch;
// The letters around a chunk that NCBI's extensions read for a hit before it or past its end.
use crate::algorithm::common::translated_chunks::ChunkLetters;
// The table walks of the two-hit extension (x_two_hit_fast.rs).
use super::x_two_hit_fast;

// Import the internal extension functions from ungapped module
// These are used for left/right extension
// The wasm simd128 build extends with `extend_left_ungapped_simd128` instead.
#[cfg(not(all(
    target_arch = "wasm32",
    target_feature = "simd128",
    not(feature = "tblastx-wasm-scalar")
)))]
#[inline(always)]
fn extend_left_ungapped_scalar<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize, Option<i32>) {
    // The left walk of the two-hit and one-hit extensions, in NCBI's form: one count `k` down over
    // the `n` letters before the offsets (the walk's step `i` reads the letters at
    // `k = n - 1 - i`), and the best's `k`; the length is `n - best_k`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:897-920
    // NCBI's offsets are those of the letters before the word (`s_off - 1`, `q_off - 1`), LOSAT's
    // one past them: LOSAT's `n` is NCBI's `n + 1`, its `k` NCBI's `i` (from the same first
    // letter), its `best_k` NCBI's `best_i` and `n - best_k` NCBI's length.
    let n = q_off.min(s_off);
    debug_assert!(q_off <= q_seq.len() && s_off <= s_seq.len());
    // SAFETY: both callers' offsets are at most the word's end in each sequence, where the word
    // fits (`two_hit_walks`, `one_hit_walks`); the walk reads the `n` letters before them.
    let (q, s) = unsafe {
        (
            q_seq.get_unchecked(q_off - n..q_off),
            s_seq.get_unchecked(s_off - n..s_off),
        )
    };

    let mut score = initial_score;
    let mut max_score = initial_score;
    let mut best_k = n;
    let mut k = n;

    while k > 0 {
        k -= 1;
        // SAFETY: `k < n`, the length of both slices.
        let (q_char, s_char) = unsafe { (*q.get_unchecked(k), *s.get_unchecked(k)) };

        score += sc.score(q_char, s_char);

        if score > max_score {
            max_score = score;
            best_k = k;
        }

        if (max_score - score) >= x_drop {
            return (max_score, n - best_k, None);
        }
    }
    let left_disp = n - best_k;
    // The walk read the context's first letter, or the subject's, without a break. NCBI's
    // `n`, the smaller of `s_off` and `q_off`, counts the batch's letters before the walk (aa_ungapped.c:897): where
    // the subject has letters left, the walk goes on before the context (query_batch.rs).
    let again = qb.left_at_bound(sc, s_seq, q_off, s_off, score, max_score, x_drop);

    (max_score, left_disp, again.then_some(score))
}

// NCBI BLAST reference (c++/src/algo/blast/core/aa_ungapped.c:886-921):
#[cfg(all(
    target_arch = "wasm32",
    target_feature = "simd128",
    not(feature = "tblastx-wasm-scalar")
))]
#[target_feature(enable = "simd128")]
unsafe fn extend_left_ungapped_simd128<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize, Option<i32>) {
    use std::arch::wasm32::*;

    let max_left = q_off.min(s_off);

    let mut score = initial_score;
    let mut max_score = initial_score;
    let mut left_disp = 0usize;
    let mut i = 0usize;

    #[expect(
        unused_assignments,
        reason = "the lane macro steps pos after every lane, the last one too; the next block starts again from its own offset"
    )]
    while i + 16 <= max_left {
        let q_ptr = q_seq.as_ptr().add(q_off - i - 16);
        let s_ptr = s_seq.as_ptr().add(s_off - i - 16);
        // Safety: q_ptr/s_ptr are in-bounds by loop condition.
        let qv = v128_load(q_ptr as *const v128);
        let sv = v128_load(s_ptr as *const v128);
        let mut pos = i;
        macro_rules! step_lane_rev {
            ($lane:expr, $qv:ident, $sv:ident, $pos:ident) => {{
                let q = i8x16_extract_lane::<$lane>($qv) as u8;
                let s = i8x16_extract_lane::<$lane>($sv) as u8;
                score += sc.score(q, s);
                let ext_len = $pos + 1;
                if score > max_score {
                    max_score = score;
                    left_disp = ext_len;
                }
                if (max_score - score) >= x_drop {
                    return (max_score, left_disp, None);
                }
                $pos += 1;
            }};
        }
        step_lane_rev!(15, qv, sv, pos);
        step_lane_rev!(14, qv, sv, pos);
        step_lane_rev!(13, qv, sv, pos);
        step_lane_rev!(12, qv, sv, pos);
        step_lane_rev!(11, qv, sv, pos);
        step_lane_rev!(10, qv, sv, pos);
        step_lane_rev!(9, qv, sv, pos);
        step_lane_rev!(8, qv, sv, pos);
        step_lane_rev!(7, qv, sv, pos);
        step_lane_rev!(6, qv, sv, pos);
        step_lane_rev!(5, qv, sv, pos);
        step_lane_rev!(4, qv, sv, pos);
        step_lane_rev!(3, qv, sv, pos);
        step_lane_rev!(2, qv, sv, pos);
        step_lane_rev!(1, qv, sv, pos);
        step_lane_rev!(0, qv, sv, pos);
        i += 16;
    }

    while i < max_left {
        let q_char = *q_seq.get_unchecked(q_off - 1 - i);
        let s_char = *s_seq.get_unchecked(s_off - 1 - i);

        score += sc.score(q_char, s_char);
        let ext_len = i + 1;

        if score > max_score {
            max_score = score;
            left_disp = ext_len;
        }

        if (max_score - score) >= x_drop {
            return (max_score, left_disp, None);
        }
        i += 1;
    }
    // As in `extend_left_ungapped_scalar` (aa_ungapped.c:897): both loops return at a break.
    let again = qb.left_at_bound(sc, s_seq, q_off, s_off, score, max_score, x_drop);

    (max_score, left_disp, again.then_some(score))
}

/// NCBI's left walk (`s_BlastAaExtendLeft`) from the letters before `q_off` and `s_off`: its best
/// sum and length, and where the search is to go on on the query batch, the walk's sum at its end
/// (on `InBatch::ALONE`, which stops at the context's start, where NCBI's walk could go on past
/// it; a walk of `QueryBatch::at` goes on into the batch as NCBI's, `InBatch::left_at_bound`).
#[inline(always)]
pub(super) fn extend_left_ungapped<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize, Option<i32>) {
    #[cfg(all(
        target_arch = "wasm32",
        target_feature = "simd128",
        not(feature = "tblastx-wasm-scalar")
    ))]
    {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:886-921
        // The `tblastx-wasm-scalar` feature is diagnostic-only: it keeps this
        // exact scalar NCBI control flow available for native-vs-Wasm parity
        // isolation while leaving production SIMD builds unchanged.
        unsafe {
            extend_left_ungapped_simd128(sc, q_seq, qb, s_seq, q_off, s_off, initial_score, x_drop)
        }
    }
    #[cfg(not(all(
        target_arch = "wasm32",
        target_feature = "simd128",
        not(feature = "tblastx-wasm-scalar")
    )))]
    {
        extend_left_ungapped_scalar(sc, q_seq, qb, s_seq, q_off, s_off, initial_score, x_drop)
    }
}

// The wasm simd128 build extends with `extend_right_ungapped_simd128` instead.
#[cfg(not(all(
    target_arch = "wasm32",
    target_feature = "simd128",
    not(feature = "tblastx-wasm-scalar")
)))]
#[inline(always)]
fn extend_right_ungapped_scalar<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize, usize, Option<i32>) {
    let q_limit = q_seq.len();

    // NCBI's `n`, the smaller of `subject->length` less `s_off` and `query->length` less `q_off` (c++/src/algo/blast/core/aa_ungapped.c:841), on the
    // context's letters: one test a step. None where the walk starts past them
    // (`extend_hit_one_hit`).
    let q_rest = q_seq.get(q_start..).unwrap_or_default();
    let s_rest = s_seq.get(s_start..).unwrap_or_default();
    let n = q_rest.len().min(s_rest.len());

    let mut score = initial_score;
    let mut max_score = initial_score;
    let mut right_disp = 0usize;
    let mut j = 0usize;

    while j < n {
        let q_char = q_rest[j];
        let s_char = s_rest[j];

        score += sc.score(q_char, s_char);

        if score > max_score {
            max_score = score;
            right_disp = j + 1;
        }

        if score <= 0 || (max_score - score) >= x_drop {
            return (max_score, right_disp, j, None);
        }
        j += 1;
    }
    // The walk read the context's last letter, or the subject's, without a break; where the
    // subject has letters left, it read the context's. NCBI's `n`, the smaller of `subject->length` less `s_off` and
    // `query->length` less `q_off`, counts the batch's letters after the walk (aa_ungapped.c:841): the
    // walk goes on after the context (query_batch.rs), and its count is NCBI's `i` there
    // (`*s_last_off = s_off + i`). The one-hit extension can start the walk one past the
    // context's end (`extend_hit_one_hit`): NCBI's walk starts on the next context's first letter
    // there.
    // A walk that ends here started past the context's letters only where it took no step, and
    // then `q_start > q_limit` (`InBatch::right_at_bound`).
    let (j, again) = qb.right_at_bound(
        sc,
        s_seq,
        s_start + j,
        q_start + j,
        q_start > q_limit,
        j,
        score,
        max_score,
        x_drop,
    );

    (max_score, right_disp, j, again.then_some(score))
}

// NCBI BLAST reference (c++/src/algo/blast/core/aa_ungapped.c:846-866):
#[cfg(all(
    target_arch = "wasm32",
    target_feature = "simd128",
    not(feature = "tblastx-wasm-scalar")
))]
#[target_feature(enable = "simd128")]
unsafe fn extend_right_ungapped_simd128<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize, usize, Option<i32>) {
    use std::arch::wasm32::*;

    let q_limit = q_seq.len();
    let s_limit = s_seq.len();
    // The one-hit extension can start the walk one past either slice's end (`extend_hit_one_hit`):
    // no step there, as in the scalar walk's loop test.
    let n = q_limit
        .saturating_sub(q_start)
        .min(s_limit.saturating_sub(s_start));

    let mut score = initial_score;
    let mut max_score = initial_score;
    let mut right_disp = 0usize;
    let mut j = 0usize;

    #[expect(
        unused_assignments,
        reason = "the lane macro steps pos after every lane, the last one too; the next block starts again from its own offset"
    )]
    while j + 16 <= n {
        let base = j;
        let q_ptr = q_seq.as_ptr().add(q_start + base);
        let s_ptr = s_seq.as_ptr().add(s_start + base);
        // Safety: q_ptr/s_ptr are in-bounds by loop condition.
        let qv = v128_load(q_ptr as *const v128);
        let sv = v128_load(s_ptr as *const v128);
        let mut pos = base;
        macro_rules! step_lane {
            ($lane:expr, $qv:ident, $sv:ident, $pos:ident) => {{
                let q = i8x16_extract_lane::<$lane>($qv) as u8;
                let s = i8x16_extract_lane::<$lane>($sv) as u8;
                score += sc.score(q, s);
                if score > max_score {
                    max_score = score;
                    right_disp = $pos + 1;
                }
                if score <= 0 || (max_score - score) >= x_drop {
                    return (max_score, right_disp, $pos, None);
                }
                $pos += 1;
            }};
        }
        step_lane!(0, qv, sv, pos);
        step_lane!(1, qv, sv, pos);
        step_lane!(2, qv, sv, pos);
        step_lane!(3, qv, sv, pos);
        step_lane!(4, qv, sv, pos);
        step_lane!(5, qv, sv, pos);
        step_lane!(6, qv, sv, pos);
        step_lane!(7, qv, sv, pos);
        step_lane!(8, qv, sv, pos);
        step_lane!(9, qv, sv, pos);
        step_lane!(10, qv, sv, pos);
        step_lane!(11, qv, sv, pos);
        step_lane!(12, qv, sv, pos);
        step_lane!(13, qv, sv, pos);
        step_lane!(14, qv, sv, pos);
        step_lane!(15, qv, sv, pos);
        j += 16;
    }

    while j < n {
        let q_char = *q_seq.get_unchecked(q_start + j);
        let s_char = *s_seq.get_unchecked(s_start + j);

        score += sc.score(q_char, s_char);

        if score > max_score {
            max_score = score;
            right_disp = j + 1;
        }

        if score <= 0 || (max_score - score) >= x_drop {
            return (max_score, right_disp, j, None);
        }
        j += 1;
    }
    // As in `extend_right_ungapped_scalar` (aa_ungapped.c:841): both loops return at a break, so
    // `j == n`, and where the subject has letters left, `n` is the context's count.
    // A walk that ends here started past the context's letters only where it took no step, and
    // then `q_start > q_limit` (`InBatch::right_at_bound`).
    let (j, again) = qb.right_at_bound(
        sc,
        s_seq,
        s_start + j,
        q_start + j,
        q_start > q_limit,
        j,
        score,
        max_score,
        x_drop,
    );

    (max_score, right_disp, j, again.then_some(score))
}

/// NCBI's right walk (`s_BlastAaExtendRight`) from `q_start` and `s_start`: its best sum and
/// length, its count (`*s_last_off = s_off + i`), and where the search is to go on on the query
/// batch, the walk's sum at its end (on `InBatch::ALONE`, which stops at the context's end, where
/// NCBI's walk could go on past it; a walk of `QueryBatch::at` goes on into the batch as NCBI's,
/// `InBatch::right_at_bound`).
#[inline(always)]
pub(super) fn extend_right_ungapped<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize, usize, Option<i32>) {
    #[cfg(all(
        target_arch = "wasm32",
        target_feature = "simd128",
        not(feature = "tblastx-wasm-scalar")
    ))]
    {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:846-866
        // The `tblastx-wasm-scalar` feature is diagnostic-only: it keeps this
        // exact scalar NCBI control flow available for native-vs-Wasm parity
        // isolation while leaving production SIMD builds unchanged.
        unsafe {
            extend_right_ungapped_simd128(
                sc,
                q_seq,
                qb,
                s_seq,
                q_start,
                s_start,
                initial_score,
                x_drop,
            )
        }
    }
    #[cfg(not(all(
        target_arch = "wasm32",
        target_feature = "simd128",
        not(feature = "tblastx-wasm-scalar")
    )))]
    {
        extend_right_ungapped_scalar(
            sc,
            q_seq,
            qb,
            s_seq,
            q_start,
            s_start,
            initial_score,
            x_drop,
        )
    }
}

/// The two-hit extension for the searches, without a test of the bounds: the word at `q_right_off` and
/// `s_right_off` must fit in both sequences. The searches pass lookup hits on the letters; no
/// lookup word holds a frame's NULLB (its matrix row is BLAST_SCORE_MIN).
// The values come from the table walks of x_two_hit_fast.rs, which hand back to
// `extend_hit_two_hit_reference` what they do not cover; the debug output (`debug_ext`) stays
// with the reference.
#[inline]
pub(crate) fn extend_hit_two_hit_fitting(
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
    debug_ext: bool,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    if !debug_ext {
        if let Some(r) =
            x_two_hit_fast::two_hit(q_seq, s_seq, s_left_off, s_right_off, q_right_off, x_drop)
        {
            return r;
        }
    }
    extend_hit_two_hit_reference(
        Blosum62Score,
        q_seq,
        qb,
        s_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
        debug_ext,
    )
}

/// `extend_hit_two_hit_fitting` with the scores of `sc` (the word must fit, as there): BLOSUM62's
/// instance is `extend_hit_two_hit_fitting`,
/// another matrix's takes the table walks of TBXUNGAP with its byte table
/// (`TbxScore::table_two_hit`), and the reference walks with its scores for the calls those do
/// not cover.
#[inline(always)]
pub(crate) fn extend_hit_two_hit_with<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
    debug_ext: bool,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    if S::BLOSUM62 {
        return extend_hit_two_hit_fitting(
            q_seq,
            qb,
            s_seq,
            s_left_off,
            s_right_off,
            q_right_off,
            x_drop,
            debug_ext,
        );
    }
    if !debug_ext {
        if let Some(r) =
            sc.table_two_hit(q_seq, s_seq, s_left_off, s_right_off, q_right_off, x_drop)
        {
            return r;
        }
    }
    extend_hit_two_hit_reference(
        sc,
        q_seq,
        qb,
        s_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
        debug_ext,
    )
}

/// EXPERIMENT (LOSAT_X_TBXUNGAPSENT): `extend_hit_two_hit_fitting(&q_full[1..q_full.len() - 1],
/// &s_full[1..s_full.len() - 1], s_left_off, s_right_off, q_right_off, x_drop, false)` with the
/// walks that the sentinels of the two buffers stop (x_two_hit_fast.rs).
///
/// SAFETY: `q_full` is the buffer of a context of a batch `x_two_hit_fast::sent_contexts`
/// accepted, and `s_full` the buffer of a subject frame `x_two_hit_fast::sent_subject` accepted.
/// The word at `q_right_off` and `s_right_off` lies in the letters (the buffers without their
/// sentinels), as for `extend_hit_two_hit_fitting`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit`.
#[inline(always)]
pub(crate) unsafe fn extend_hit_two_hit_sent(
    q_full: &[u8],
    qb: InBatch<'_>,
    s_full: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    // SAFETY: the caller's.
    match unsafe {
        x_two_hit_fast::two_hit_sent(q_full, s_full, s_left_off, s_right_off, q_right_off, x_drop)
    } {
        Some(r) => r,
        None => extend_hit_two_hit_sent_refused(
            q_full,
            qb,
            s_full,
            s_left_off,
            s_right_off,
            q_right_off,
            x_drop,
        ),
    }
}

/// The calls `extend_hit_two_hit_sent` leaves to `extend_hit_two_hit_fitting` (an X-drop outside
/// 1..=127; `two_hit_sent` also refuses a word past the letters, which a lookup hit's is not), out
/// of the callers' code.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit`.
#[cold]
#[inline(never)]
fn extend_hit_two_hit_sent_refused(
    q_full: &[u8],
    qb: InBatch<'_>,
    s_full: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    extend_hit_two_hit_fitting(
        &q_full[1..q_full.len() - 1],
        qb,
        &s_full[1..s_full.len() - 1],
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
        false,
    )
}

/// The two-hit extension as written before the table walks: the calls they do not cover (an X-drop
/// outside 1..=127, a word that does not fit) and the debug output, with the scores of `sc`.
///
/// Its walks go on past the context's ends into the query batch `qb` as NCBI's do (query_batch.rs).
/// They first run on the context alone, and only where one of them ends at the context's end
/// where NCBI's could go on past it (the NULLB not breaking it, `InBatch::left_at_bound`,
/// `InBatch::right_at_bound`) does the extension run again on the batch, both walks
/// (`extend_hit_two_hit_in_batch`): elsewhere the two give the same result, and the walks' loops
/// keep no values of the batch. (The one-hit extension instead goes on from where its first pass
/// stopped, `Pass::Resume`; here that form cost the first pass's left loop a register, X2c
/// `-xdrop_ungap 100` +0.4 % instructions, and the two walks again stay below NCBI's instructions
/// where they run: X2c s6i `f48_t13` 0.49x.)
pub(super) fn extend_hit_two_hit_reference<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
    debug_ext: bool,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    if !debug_ext {
        let (r, again) = two_hit_walks(
            sc,
            q_seq,
            InBatch::ALONE,
            s_seq,
            s_left_off,
            s_right_off,
            q_right_off,
            x_drop,
            false,
            true,
        );
        if !again {
            return r;
        }
    }
    extend_hit_two_hit_in_batch(
        sc,
        q_seq,
        qb,
        s_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
        debug_ext,
    )
}

/// `extend_hit_two_hit_reference` on the batch: its walks past the context's ends, and the debug
/// output.
#[cold]
#[inline(never)]
fn extend_hit_two_hit_in_batch<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
    debug_ext: bool,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    two_hit_walks(
        sc,
        q_seq,
        qb,
        s_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
        debug_ext,
        false,
    )
    .0
}

/// The steps of `extend_hit_two_hit_reference` with the walks of `qb`, and whether the extension
/// is to run again on the query batch (`InBatch::ALONE` only, where a walk ends at the context's
/// end and NCBI's could go on). `first`: the search's first pass on the context alone, which ends
/// at that walk (its result is then not used: no other value is kept past the walk).
#[inline(always)]
pub(super) fn two_hit_walks<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    // Where `q_seq` lies in the query batch: the walks go on past its ends as NCBI's do on the
    // batch (aa_ungapped.c:841, 897; query_batch.rs).
    qb: InBatch<'_>,
    s_seq: &[u8],
    s_left_off: usize,  // Position of first hit (L) in subject
    s_right_off: usize, // Position of second hit (R) in subject
    q_right_off: usize, // Position of second hit (R) in query
    x_drop: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1096
    // Existing LOSAT diagnostics only; supplied by the current search owner.
    debug_ext: bool,
    first: bool,
) -> ((usize, usize, usize, usize, i32, bool, usize), bool) {
    // DEBUG: Print sequence segment for tracing
    static DEBUG_PRINTED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1104-1108 (the word length of the
    // lookup table, `-word_size`: 3 for BLOSUM62's instance, `TbxScore::word`)
    let k_size = sc.word();

    // First, find the best scoring position within the word at R
    // NCBI BLAST assumes caller provides valid coordinates within bounds.
    // The lookup table and offset generation guarantee valid positions.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c does NOT have bounds checks in word scanning loop.
    let mut score = 0i32;
    let mut seed_score = 0i32;
    let mut right_d = 0usize;

    for i in 0..k_size {
        // SAFETY: Lookup table guarantees q_right_off and s_right_off have at least k_size valid positions
        let q_char = unsafe { *q_seq.get_unchecked(q_right_off + i) };
        let s_char = unsafe { *s_seq.get_unchecked(s_right_off + i) };
        score += sc.score(q_char, s_char);

        if score > seed_score {
            seed_score = score;
            right_d = i + 1; // Position is one beyond the end of the word
        }
    }

    let q_right_off = q_right_off + right_d;
    let s_right_off = s_right_off + right_d;

    // Extend to the left from R, trying to reach L
    // Note: NCBI BLAST does NOT explicitly check for stop codons in s_BlastAaExtendLeft.
    // Stop codons have very negative scores in BLOSUM62 (-4 to -5), so X-drop handles them.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:886-921
    let (left_score, left_disp, left_again) =
        extend_left_ungapped(sc, q_seq, qb, s_seq, q_right_off, s_right_off, 0, x_drop);
    if first && left_again.is_some() {
        return ((0, 0, 0, 0, 0, false, 0), true);
    }

    // Check if left extension reached the first hit
    // The distance from R to L is (s_right_off - s_left_off)
    // We need left_disp >= this distance for the extension to reach L
    let distance_to_first_hit = s_right_off.saturating_sub(s_left_off);
    let reached_first_hit = left_disp >= distance_to_first_hit;

    // Only extend to the right if left extension reached the first hit (NCBI BLAST behavior)
    let mut right_disp = 0usize;
    let mut right_score = 0i32;
    let mut right_extended = false;
    let mut s_last_off = s_right_off; // Default to current position if no right extension

    if reached_first_hit {
        // NCBI BLAST: right_extend is set to TRUE when we enter this block,
        // regardless of how many characters are actually extended.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1140
        right_extended = true;

        // NCBI BLAST extends right starting from left_score as initial score
        let (new_max_score, new_right_disp, j, right_again) = extend_right_ungapped(
            sc,
            q_seq,
            qb,
            s_seq,
            q_right_off,
            s_right_off,
            left_score,
            x_drop,
        );
        if first && right_again.is_some() {
            return ((0, 0, 0, 0, 0, false, 0), true);
        }
        right_score = new_max_score;
        right_disp = new_right_disp;
        // s_last_off is the rightmost subject position scanned (NCBI BLAST style)
        s_last_off = s_right_off + j;
    }

    // NCBI BLAST returns MAX(left_score, right_score)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1157
    let max_score_total = left_score.max(right_score);

    let q_start = q_right_off - left_disp;
    let q_end = q_right_off + right_disp;
    let s_start = s_right_off - left_disp;
    let s_end = s_right_off + right_disp;

    // DEBUG: Print alignment details for first few high-scoring extensions
    if debug_ext
        && max_score_total >= 100
        && !DEBUG_PRINTED.swap(true, std::sync::atomic::Ordering::Relaxed)
    {
        eprintln!(
            "[DEBUG_EXT] left_score={} right_score={} max_score_total={}",
            left_score, right_score, max_score_total
        );
        eprintln!(
            "[DEBUG_EXT] q_start={} q_end={} s_start={} s_end={} score={}",
            q_start, q_end, s_start, s_end, max_score_total
        );
        eprintln!(
            "[DEBUG_EXT] left_disp={} right_disp={} x_drop={}",
            left_disp, right_disp, x_drop
        );
        // Print position-by-position scores
        let align_len = q_end - q_start;
        eprintln!("[DEBUG_EXT] alignment_len={}", align_len);
        let mut total = 0i32;
        for i in 0..align_len.min(100) {
            let q_char = q_seq[q_start + i];
            let s_char = s_seq[s_start + i];
            let cell = sc.score(q_char, s_char);
            total += cell;
            if i < 10 || i >= align_len - 5 {
                eprintln!(
                    "[DEBUG_EXT] pos={:2} q={:2} s={:2} score={:3} cumul={:4}",
                    i, q_char, s_char, cell, total
                );
            }
        }
        eprintln!(
            "[DEBUG_EXT] manual_total={} reported_score={}",
            total, max_score_total
        );
    }

    (
        (
            q_start,
            q_end,
            s_start,
            s_end,
            max_score_total,
            right_extended,
            s_last_off,
        ),
        false,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158, :897 and :841
/// NCBI's two-hit extension (`s_BlastAaExtendTwoHit`) of a hit whose first hit (`s_left_off`)
/// lies before its chunk, with the scores of `sc`, on the letters around the chunk (`ChunkLetters`, whose type doc gives NCBI's
/// reads there): the word on the frame's letters, the left walk down to the chunk's first letter
/// (none where the word's best end is at or before it), the right walk up to the chunk's last
/// letter. The walks go on into the query batch `qb` past the context's ends as NCBI's do. The
/// subject offsets are the chunk's (`Int4`, below 0 before it): `(q_start, q_end, s_start, s_end,
/// score, right_extend, s_last_off)`.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn extend_hit_two_hit_outside<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    qb: InBatch<'_>,
    chunk: ChunkLetters<'_>,
    s_left_off: i32,
    s_right_off: i32,
    q_right_off: usize,
    x_drop: i32,
) -> (usize, usize, i32, i32, i32, bool, i32) {
    let mut score = 0i32;
    let mut seed_score = 0i32;
    let mut right_d = 0usize;
    for i in 0..sc.word() {
        score += sc.score(q_seq[q_right_off + i], chunk.at(s_right_off + i as i32));
        if score > seed_score {
            seed_score = score;
            right_d = i + 1;
        }
    }
    let q_right_off = q_right_off + right_d;
    let s_right_off = s_right_off + right_d as i32;
    let (left_score, left_disp) = if s_right_off > 0 {
        let (max_score, disp, _) = extend_left_ungapped(
            sc,
            q_seq,
            qb,
            chunk.chunk_onward(),
            q_right_off,
            s_right_off as usize,
            0,
            x_drop,
        );
        (max_score, disp as i32)
    } else {
        (0, 0)
    };
    let (mut right_score, mut right_disp, mut right_extend, mut s_last_off) =
        (0i32, 0i32, false, s_right_off);
    if left_disp >= s_right_off - s_left_off {
        right_extend = true;
        let (max_score, disp, j, _) = extend_right_ungapped(
            sc,
            q_seq,
            qb,
            chunk.to_chunk_end(),
            q_right_off,
            (chunk.lead as i32 + s_right_off) as usize,
            left_score,
            x_drop,
        );
        (right_score, right_disp, s_last_off) = (max_score, disp as i32, s_right_off + j as i32);
    }
    (
        q_right_off - left_disp as usize,
        q_right_off + right_disp as usize,
        s_right_off - left_disp,
        s_right_off + right_disp,
        left_score.max(right_score),
        right_extend,
        s_last_off,
    )
}

// No NCBI counterpart (LOSAT_X_STATS): the cells one call of `extend_hit_two_hit_reference` scores,
// with the search's scores and word length, for the statistics of the unit path
// (tblastx/x_unit_probe.rs). It repeats the call's scalar loops (the same sums and the same stops,
// aa_ungapped.c:846-866 and 886-921) and counts their steps; the search never calls it. The test
// below checks that it returns what the call returns. It counts the cells of the context alone
// (`InBatch::ALONE`): a walk that goes on into the query batch past the context's ends (a raw
// X-drop above 32,768) takes more cells and may end elsewhere than it counts.
/// What one call of `extend_hit_two_hit_reference` scored, and its result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct XTwoHitCells {
    /// Cells of the left extension (from the word's best end towards the first hit).
    pub(crate) left: usize,
    /// Cells of the right extension; 0 when it did not run.
    pub(crate) right: usize,
    /// First subject position scored (the word's cells included).
    pub(crate) s_lo: usize,
    /// One past the last subject position scored.
    pub(crate) s_hi: usize,
    /// The seven values `extend_hit_two_hit_reference` returns for the same arguments.
    pub(crate) result: (usize, usize, usize, usize, i32, bool, usize),
}

/// The cells of `extend_hit_two_hit_reference(sc, q_seq, InBatch::ALONE, s_seq, s_left_off,
/// s_right_off, q_right_off, x_drop, false)`: the scores of `sc` and its word length
/// (`TbxScore::word`).
pub(crate) fn x_two_hit_cells<S: TbxScore>(
    sc: S,
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> XTwoHitCells {
    let k_size = sc.word();
    let (mut score, mut seed_score, mut right_d) = (0i32, 0i32, 0usize);
    for i in 0..k_size {
        score += sc.score(q_seq[q_right_off + i], s_seq[s_right_off + i]);
        if score > seed_score {
            seed_score = score;
            right_d = i + 1;
        }
    }
    let (q_r, s_r) = (q_right_off + right_d, s_right_off + right_d);

    let max_left = q_r.min(s_r);
    let (mut score, mut left_score, mut left_disp, mut left) = (0i32, 0i32, 0usize, 0usize);
    while left < max_left {
        score += sc.score(q_seq[q_r - 1 - left], s_seq[s_r - 1 - left]);
        left += 1;
        if score > left_score {
            left_score = score;
            left_disp = left;
        }
        if (left_score - score) >= x_drop {
            break;
        }
    }

    let reached_first_hit = left_disp >= s_r.saturating_sub(s_left_off);
    let (mut right, mut right_disp, mut right_score, mut s_last_off) = (0usize, 0usize, 0i32, s_r);
    if reached_first_hit {
        let (mut score, mut max_score, mut j) = (left_score, left_score, 0usize);
        while q_r + j < q_seq.len() && s_r + j < s_seq.len() {
            score += sc.score(q_seq[q_r + j], s_seq[s_r + j]);
            right += 1;
            if score > max_score {
                max_score = score;
                right_disp = j + 1;
            }
            if score <= 0 || (max_score - score) >= x_drop {
                break;
            }
            j += 1;
        }
        right_score = max_score;
        s_last_off = s_r + j;
    }

    XTwoHitCells {
        left,
        right,
        s_lo: s_right_off.min(s_r - left),
        s_hi: (s_right_off + k_size).max(s_r + right),
        result: (
            q_r - left_disp,
            q_r + right_disp,
            s_r - left_disp,
            s_r + right_disp,
            left_score.max(right_score),
            reached_first_hit,
            s_last_off,
        ),
    }
}

#[cfg(test)]
#[path = "two_hit_x_cells_tests.rs"]
mod x_cells_tests;
