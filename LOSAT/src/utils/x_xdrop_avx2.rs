//! EXPERIMENT (LOSAT_X_DPAVX2): 256-bit AVX2 kernels for the X-drop dynamic programming of
//! `xdrop_simd.rs`: 16 lanes of 16 bits and 32 lanes of 8 bits, with (`TB`) and without the
//! traceback. x86_64 only, chosen at run time when the CPU has AVX2. They run the recurrence of
//! NCBI `ALIGN_EX` / `Blast_SemiGappedAlign` row by row, with the value representation of the
//! 128-bit kernels, and return exactly the `XdropResult` and edit operations of those kernels.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:563-635
//! ```c
//!         for (b_index = first_b_index; b_index < b_size; b_index++) {
//! ...
//!             if (matrix_index == FENCE_SENTRY) {
//! ...
//!             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
//! ...
//!             if (score < score_gap_col) {
//!                 script = SCRIPT_GAP_IN_B;
//!                 score = score_gap_col;
//!             }
//!             if (score < score_gap_row) {
//!                 script = SCRIPT_GAP_IN_A;
//!                 score = score_gap_row;
//!             }
//!
//!             if (best_score - score > x_dropoff) {
//!
//!                 if (first_b_index == b_index)
//!                     first_b_index++;
//!                 else
//!                     score_array[b_index].best = MININT;
//!             }
//!             else {
//!                 last_b_index = b_index;
//!                 if (score > best_score) {
//!                     best_score = score;
//!                     *a_offset = a_index;
//! ...
//!                 if (score_gap_row < (score - gap_open_extend))
//!                     score_gap_row = score - gap_open_extend;
//!                 else
//!                     script += script_row;
//!
//!                 score_array[b_index].best = score;
//! ...
//!             edit_script_row[b_index] = script;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:652-673
//! ```c
//!         if (last_b_index < b_size - 1) {
//!             b_size = last_b_index + 1;
//!         }
//!         else {
//!             while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {
//!
//!                 score_array[b_size].best = score_gap_row;
//!                 score_array[b_size].best_gap = score_gap_row - gap_open_extend;
//!                 score_gap_row -= gap_extend;
//!                 edit_script_row[b_size] = SCRIPT_GAP_IN_A;
//!                 b_size++;
//! ...
//!         if (b_size <= N) {
//!             score_array[b_size].best = MININT;
//!             score_array[b_size].best_gap = MININT;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:862-872,903-908
//! ```c
//!         for (b_index = first_b_index; b_index < b_size; b_index++) {
//!
//!             b_ptr += b_increment;
//!             score_gap_col = score_array[b_index].best_gap;
//!             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
//!
//!             if (score < score_gap_col)
//!                 score = score_gap_col;
//!
//!             if (score < score_gap_row)
//!                 score = score_gap_row;
//! ...
//!                 score_gap_row -= gap_extend;
//!                 score_gap_col -= gap_extend;
//!                 score_array[b_index].best_gap = MAX(score - gap_open_extend,
//!                                                     score_gap_col);
//!                 score_gap_row = MAX(score - gap_open_extend, score_gap_row);
//!                 score_array[b_index].best = score;
//! ```
//!
//! Why the result is exact. The argument of `xdrop_simd.rs` holds unchanged: with
//! `thr = best_score - x_dropoff` at the moment a cell is evaluated (never decreasing), a cell
//! that survives has `score >= thr` and its score is decided by candidates `>= thr`; whatever
//! the vector code computes differently from the C loop (gap states behind a dropped cell,
//! saturated values, minus infinity plus a score) is `< thr` in both computations, so it cannot
//! win the `max` of a surviving cell, change which cells are dropped, or change a script bit the
//! traceback walk reads (the walk only visits surviving cells). Each lane here holds exactly the
//! value the 128-bit kernel holds for the same cell:
//!
//! - Same stored values (`score - base + bias`, 0 = minus infinity), the same `bias`, rebase
//!   limit and acceptance tests, and lookup tables with the bytes of `build_row_table`.
//! - Diagonal, column gap, X-drop test and script bits: the same lane operations, 16 or 32
//!   lanes at a time instead of 8 or 16.
//! - Row gap (`score_gap_row`): both kernels compute, with saturating arithmetic,
//!   `F[i] = max(max_{j<i} (v[j] - (i-1-j)*gap_extend), carry - i*gap_extend)` with
//!   `v = max(next_score, score_gap_col) - gap_open_extend` and `carry` the row gap entering the
//!   block. The 128-bit kernel uses shifted maxima with decay; this file uses the identity
//!   `max_j sat(v[j] - (i-1-j)g) = sat(max_j (v[j] + (j+1)g) - i*g)` (a plain prefix maximum),
//!   exact because `v + 15g` (16-bit) or `v + 31g` (8-bit) fits the lane, which is checked
//!   (otherwise None). The carry to the next block is `max(carry - L*g, max(F_in[L-1] - g,
//!   v[L-1]))`, the same value as the 128-bit kernel's `max(H - gap_open_extend, F - g)`.
//! - Running best: the per-lane exclusive prefix maximum of the row, as in the 128-bit kernel.
//! - Band steps: the blocks of a row start at `first_b_index + 16k`, as in the 128-bit kernels,
//!   and past the band (`b_size`) a row stops at the first such point where the row gap entering
//!   the next block is below the threshold, or at the end of `B` (or the fence). A 32-lane step
//!   is taken while more than 16 lanes of the band are left, so a row may run up to 16 lanes past
//!   the point where the 128-bit kernel stops. Those lanes lie past the band with a row gap below
//!   the threshold and nothing above it from the previous row, so they are dropped (H = 0, E
//!   cleared right of the last survivor) and change nothing. The first and last survivors, hence
//!   the band of the next row, the cell count, the best cell and the traceback rows the walk
//!   reads are those of the 128-bit kernel. The traceback walk is the 128-bit kernels' `finish`.
//!
//! Less work per row and per call than the 128-bit kernels: wider steps, first and last
//! survivor tracked while the row runs (no survivor bitmap), the first lane's missing diagonal
//! as a lane mask instead of a store under the next load, the column residues copied 32 at a
//! time, row tables of `Rows28` / `Flat` matrices built with vector instructions, and the row
//! gap carried from block to block with two operations.

#![allow(clippy::too_many_arguments)]

use super::kernel::finish;
use super::{
    build_row_table, fits_8bit, ColSeq, RowSeq, Scores, XdropResult, XdropScratch, FENCE_SENTRY,
    MAX_M, PAD, SCRIPT_EXTEND_GAP_A, SCRIPT_EXTEND_GAP_B, SCRIPT_GAP_IN_A, SCRIPT_SUB, TAB, VMAX,
};
use core::arch::x86_64::*;

type W = __m256i;
type X = __m128i;

// No NCBI counterpart: shuffle and lane masks of the 256-bit lanes; they do not change any value
// NCBI computes.
/// pshufb mask: the last 16-bit lane of each 128-bit half to the whole half.
static KB7X2: [u8; 32] = [
    14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15,
    14, 15, 14, 15, 14, 15, 14, 15,
];
/// `LANE_MASK16[16 - n..]` has the first `n` 16-bit lanes set.
static LANE_MASK16: [i16; 32] = [
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
];
/// `LANE_MASK8[32 - n..]` has the first `n` bytes set.
static LANE_MASK8: [u8; 64] = [
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];
/// `LANE_MASK8X[16 - n..]` has the first `n` bytes set.
static LANE_MASK8X: [u8; 32] = [
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];
/// pshufb mask: low bytes then high bytes of eight 16-bit lanes, per 128-bit half.
static KSPLIT: [u8; 32] = [
    0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15, 0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9,
    11, 13, 15,
];

/// No NCBI counterpart: LOSAT_X_DPAVX2 (any value) and a CPU that runs these kernels.
pub(super) fn switch_on() -> bool {
    std::env::var_os("LOSAT_X_DPAVX2").is_some() && available()
}

/// No NCBI counterpart: AVX2 and the features of the 128-bit AVX build (LOSAT_X_DPNOAVX, which
/// turns that build off, turns these kernels off too).
pub(super) fn available() -> bool {
    super::cpu_level() == 2 && is_x86_feature_detected!("avx2")
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:766-770
/// ```c
///     if (!score_only) {
///         return ALIGN_EX(A, B, M, N, a_offset, b_offset, edit_block, gap_align,
/// ```
/// Dispatch point of the 256-bit kernels: `tb == false` is the score-only recurrence
/// (`Blast_SemiGappedAlign`, `s_BlastAlignPackedNucl`), `tb == true` is `ALIGN_EX`. The lane width
/// is the one the 128-bit dispatch would choose (`fits_8bit`). None: not handled here; the caller
/// runs the 128-bit dispatch.
///
/// # Safety
/// `available()` must be true.
pub(super) unsafe fn align(
    q: &RowSeq<'_>,
    s: &ColSeq<'_>,
    scores: &Scores<'_>,
    len1: usize,
    len2: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    tb: bool,
    check_fence: bool,
    sc: &mut XdropScratch,
) -> Option<XdropResult> {
    let m8 = fits_8bit(scores, gap_open, gap_extend, x_drop);
    match (tb, m8) {
        (true, Some(m)) => align8::<true>(
            q,
            s,
            scores,
            m,
            len1,
            len2,
            gap_open,
            gap_extend,
            x_drop,
            check_fence,
            sc,
        ),
        (false, Some(m)) => align8::<false>(
            q,
            s,
            scores,
            m,
            len1,
            len2,
            gap_open,
            gap_extend,
            x_drop,
            check_fence,
            sc,
        ),
        (true, None) => align16::<true>(
            q,
            s,
            scores,
            len1,
            len2,
            gap_open,
            gap_extend,
            x_drop,
            check_fence,
            sc,
        ),
        (false, None) => align16::<false>(
            q,
            s,
            scores,
            len1,
            len2,
            gap_open,
            gap_extend,
            x_drop,
            check_fence,
            sc,
        ),
    }
}

// ---------------------------------------------------------------------------
// Lane helpers (inlined into the `target_feature` entry points)
// ---------------------------------------------------------------------------

#[inline(always)]
unsafe fn ld<T>(p: *const T) -> W {
    _mm256_loadu_si256(p as *const W)
}
#[inline(always)]
unsafe fn st<T>(p: *mut T, v: W) {
    _mm256_storeu_si256(p as *mut W, v)
}
#[inline(always)]
unsafe fn ldx<T>(p: *const T) -> X {
    _mm_loadu_si128(p as *const X)
}
#[inline(always)]
unsafe fn stx<T>(p: *mut T, v: X) {
    _mm_storeu_si128(p as *mut X, v)
}
/// 16 bytes to both 128-bit halves.
#[inline(always)]
unsafe fn bcast128(p: *const u8) -> W {
    _mm256_broadcastsi128_si256(ldx(p))
}
/// 16-bit lanes moved up by one across the whole vector; lane 0 becomes 0.
#[inline(always)]
unsafe fn up1_16(v: W) -> W {
    _mm256_alignr_epi8::<14>(v, _mm256_permute2x128_si256::<0x08>(v, v))
}
/// Inclusive prefix maximum (unsigned 16-bit) over the 16 lanes.
#[inline(always)]
unsafe fn pmax16(w: W, kb7: W) -> W {
    let p = _mm256_max_epu16(w, _mm256_bslli_epi128::<2>(w));
    let p = _mm256_max_epu16(p, _mm256_bslli_epi128::<4>(p));
    let p = _mm256_max_epu16(p, _mm256_bslli_epi128::<8>(p));
    let c = _mm256_shuffle_epi8(p, kb7);
    _mm256_max_epu16(p, _mm256_permute2x128_si256::<0x08>(c, c))
}
/// 16-bit lane 15 to every lane.
#[inline(always)]
unsafe fn last16(v: W, kb7: W) -> W {
    _mm256_permute4x64_epi64::<0xFF>(_mm256_shuffle_epi8(v, kb7))
}
/// 8-bit lanes moved up by one across the whole vector; lane 0 becomes 0.
#[inline(always)]
unsafe fn up1_8(v: W) -> W {
    _mm256_alignr_epi8::<15>(v, _mm256_permute2x128_si256::<0x08>(v, v))
}
/// Inclusive prefix maximum (unsigned 8-bit) over the 32 lanes.
#[inline(always)]
unsafe fn pmax8(w: W, kb15: W) -> W {
    let p = _mm256_max_epu8(w, _mm256_bslli_epi128::<1>(w));
    let p = _mm256_max_epu8(p, _mm256_bslli_epi128::<2>(p));
    let p = _mm256_max_epu8(p, _mm256_bslli_epi128::<4>(p));
    let p = _mm256_max_epu8(p, _mm256_bslli_epi128::<8>(p));
    let c = _mm256_shuffle_epi8(p, kb15);
    _mm256_max_epu8(p, _mm256_permute2x128_si256::<0x08>(c, c))
}
/// 8-bit lane 31 to every lane.
#[inline(always)]
unsafe fn last8(v: W, kb15: W) -> W {
    _mm256_permute4x64_epi64::<0xFF>(_mm256_shuffle_epi8(v, kb15))
}
/// Inclusive prefix maximum (unsigned 8-bit) over 16 lanes.
#[inline(always)]
unsafe fn pmax8x(w: X) -> X {
    let p = _mm_max_epu8(w, _mm_bslli_si128::<1>(w));
    let p = _mm_max_epu8(p, _mm_bslli_si128::<2>(p));
    let p = _mm_max_epu8(p, _mm_bslli_si128::<4>(p));
    _mm_max_epu8(p, _mm_bslli_si128::<8>(p))
}
/// 32-entry byte lookup, per 128-bit half (an index >= 32 gives one of the table bytes or 0).
#[inline(always)]
unsafe fn lut32(t0: W, t1: W, i0: W, i1: W) -> W {
    _mm256_or_si256(_mm256_shuffle_epi8(t0, i0), _mm256_shuffle_epi8(t1, i1))
}
#[inline(always)]
unsafe fn lut32x(t0: X, t1: X, idx: X) -> X {
    let i0 = _mm_add_epi8(idx, _mm_set1_epi8(0x70));
    let i1 = _mm_sub_epi8(idx, _mm_set1_epi8(0x10));
    _mm_or_si128(_mm_shuffle_epi8(t0, i0), _mm_shuffle_epi8(t1, i1))
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543,578
/// ```c
///                 matrix_row = matrix[ A[ M - a_index ] ];
/// ...
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
/// ```
/// The same lookup planes as `build_row_table` (same bytes, same rejection of scores above the
/// 16-bit path's range), built with vector instructions for `Rows28` and `Flat` matrices; other
/// sources use `build_row_table`.
#[inline(always)]
unsafe fn build_row(scores: &Scores<'_>, q: usize, out: &mut [u8]) -> bool {
    let row: &[i32] = match *scores {
        Scores::Rows28(m) => {
            if q >= 28 {
                return false;
            }
            &m[q]
        }
        Scores::Flat { data, n } => {
            if q >= n || n > 32 {
                return false;
            }
            &data[q * n..q * n + n]
        }
        _ => return build_row_table(scores, q, out),
    };
    assert!(out.len() >= TAB);
    let n = row.len();
    let lanes = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    let mut v = [_mm256_setzero_si256(); 4];
    for (k, slot) in v.iter_mut().enumerate() {
        if 8 * k < n {
            // lanes 8k + i with 8k + i < n; the others read nothing and are 0
            let m = _mm256_cmpgt_epi32(_mm256_set1_epi32((n - 8 * k) as i32), lanes);
            *slot = _mm256_maskload_epi32(row.as_ptr().add(8 * k), m);
        }
    }
    let mx = _mm256_set1_epi32(MAX_M);
    let over = _mm256_or_si256(
        _mm256_or_si256(_mm256_cmpgt_epi32(v[0], mx), _mm256_cmpgt_epi32(v[1], mx)),
        _mm256_or_si256(_mm256_cmpgt_epi32(v[2], mx), _mm256_cmpgt_epi32(v[3], mx)),
    );
    if _mm256_movemask_epi8(over) != 0 {
        return false;
    }
    // i32 -> i16 saturates below at -32768 (the clamp of `build_row_table`)
    let w0 = _mm256_permute4x64_epi64::<0xD8>(_mm256_packs_epi32(v[0], v[1]));
    let w1 = _mm256_permute4x64_epi64::<0xD8>(_mm256_packs_epi32(v[2], v[3]));
    // i16 -> i8 saturates to -128..=127 (the 8-bit plane)
    let b8 = _mm256_permute4x64_epi64::<0xD8>(_mm256_packs_epi16(w0, w1));
    let ks = ld(KSPLIT.as_ptr());
    let x0 = _mm256_permute4x64_epi64::<0xD8>(_mm256_shuffle_epi8(w0, ks));
    let x1 = _mm256_permute4x64_epi64::<0xD8>(_mm256_shuffle_epi8(w1, ks));
    let p = out.as_mut_ptr();
    st(p, _mm256_permute2x128_si256::<0x20>(x0, x1));
    st(p.add(32), _mm256_permute2x128_si256::<0x31>(x0, x1));
    st(p.add(64), b8);
    true
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:554-557,566-576
/// ```c
///         if(reverse_sequence)
///             b_ptr = &B[N - first_b_index];
///         else
///             b_ptr = &B[first_b_index];
/// ...
///             b_ptr += b_increment;
/// ...
///             matrix_index = *b_ptr;
///
///             if (matrix_index == FENCE_SENTRY) {
/// ```
/// The column residues `s.get(k)` (`*b_ptr`) for `k` in `sc.sb_filled..to` into the lane buffer,
/// with the fence and bad-residue bookkeeping of the 128-bit kernels (first column of each).
/// 32 columns at a time while they are a plain forward or backward run of `s.data` with no
/// residue `>= ncols`; the rest, and any block with such a residue, by the scalar loop of the
/// 128-bit kernels.
#[inline(always)]
unsafe fn fill_cols(
    s: &ColSeq<'_>,
    sc: &mut XdropScratch,
    to: usize,
    ncols: usize,
    check_fence: bool,
) {
    let mut k = sc.sb_filled;
    if ncols > 0 && (s.step == 1 || s.step == -1) {
        let lim = _mm256_set1_epi8((ncols - 1) as i8);
        let rev = _mm256_setr_epi8(
            15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7,
            6, 5, 4, 3, 2, 1, 0,
        );
        let dl = s.data.len() as isize;
        while k + 32 <= to && k + 32 <= s.zero_from {
            let v = if s.step == 1 {
                // columns k..k+32 are s.data[i..i + 32]
                let i = s.base + k as isize;
                if i < 0 || i + 32 > dl {
                    break;
                }
                ld(s.data.as_ptr().offset(i))
            } else {
                // columns k..k+32 are s.data[i], s.data[i - 1], ..., s.data[i - 31]
                let i = s.base - k as isize;
                if i - 31 < 0 || i >= dl {
                    break;
                }
                _mm256_permute4x64_epi64::<0x4E>(_mm256_shuffle_epi8(
                    ld(s.data.as_ptr().offset(i - 31)),
                    rev,
                ))
            };
            if _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epu8(v, lim), lim)) != -1 {
                break;
            }
            st(sc.sb.as_mut_ptr().add(PAD + k + 1), v);
            k += 32;
        }
    }
    for k in k..to {
        let ch = s.get(k);
        if ch as usize >= ncols {
            if check_fence && ch == FENCE_SENTRY {
                if sc.fence_k == usize::MAX {
                    sc.fence_k = k;
                }
            } else if sc.bad_k == usize::MAX {
                sc.bad_k = k;
            }
        }
        *sc.sb.get_unchecked_mut(PAD + k + 1) = ch;
    }
    sc.sb_filled = to;
}

// ---------------------------------------------------------------------------
// 16-bit lanes, 16 per vector
// ---------------------------------------------------------------------------

// No NCBI counterpart: the gap costs and constants of one call on 16-bit lanes.
struct C16 {
    goe: W,
    ge1: W,
    ge16: W,
    /// `i * gap_extend` in lane `i`.
    ramp: W,
    xv: W,
    kb7: W,
    c3: W,
    c10: W,
    c40: W,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:601,610-612
// ```c
//             if (best_score - score > x_dropoff) {
// ...
//                 if (score > best_score) {
//                     best_score = score;
//                     *a_offset = a_index;
// ```
// The row gap entering the next block (`score_gap_row`), the running `best_score`, the threshold
// `best_score - x_dropoff` and the cell where the best was last raised, carried across the blocks
// of a row (all lanes hold the same value).
struct Row {
    carry: W,
    bestv: W,
    thrv: W,
    rec_b: usize,
    slow: bool,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:578-601,616-631
/// ```c
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
/// ...
///             if (score < score_gap_col) {
///                 script = SCRIPT_GAP_IN_B;
///                 score = score_gap_col;
///             }
///             if (score < score_gap_row) {
///                 script = SCRIPT_GAP_IN_A;
///                 score = score_gap_row;
///             }
///
///             if (best_score - score > x_dropoff) {
/// ...
///                 score_gap_row -= gap_extend;
///                 score_gap_col -= gap_extend;
///                 if (score_gap_col < (score - gap_open_extend)) {
///                     score_array[b_index].best_gap = score - gap_open_extend;
/// ...
///                 score_array[b_index].best = score;
/// ```
/// One block of 16 cells of one row: `D` = `next_score` of the diagonal, `E` = `score_gap_col`,
/// `F` = `score_gap_row`, `H` = `score`, the same lane operations as the 128-bit `block`, with the
/// row gap as described at the head of this file. Returns (drop mask, script).
#[inline(always)]
unsafe fn block16<const TB: bool>(
    c: &C16,
    rs: &mut Row,
    hr: *const i16,
    hw: *mut i16,
    e: *mut i16,
    b: usize,
    s: W,
    lane_mask: Option<W>,
    hpm: W,
) -> (W, W) {
    // `hpm`: lanes with no diagonal (the first lane of the band), as a 0 in the H row
    let hp = _mm256_andnot_si256(hpm, ld(hr.add(b).sub(1)));
    let ep = ld(e.add(b));
    let d = _mm256_adds_epi16(hp, s);
    let h0 = _mm256_max_epi16(d, ep);
    let v = _mm256_subs_epu16(h0, c.goe);
    // row gap opened inside the block: sat(max_j (v[j] + (j+1)g) - i*g) over j < i
    let fi = _mm256_subs_epu16(pmax16(_mm256_adds_epu16(up1_16(v), c.ramp), c.kb7), c.ramp);
    let f = _mm256_max_epi16(fi, _mm256_subs_epu16(rs.carry, c.ramp));
    let mut hh = _mm256_max_epi16(h0, f);
    if let Some(m) = lane_mask {
        hh = _mm256_and_si256(hh, m);
    }
    let hg = _mm256_subs_epu16(hh, c.goe);
    let fs = _mm256_subs_epu16(f, c.ge1);
    let out = last16(_mm256_max_epi16(_mm256_subs_epu16(fi, c.ge1), v), c.kb7);
    rs.carry = _mm256_max_epi16(_mm256_subs_epu16(rs.carry, c.ge16), out);
    let es = _mm256_subs_epu16(ep, c.ge1);
    let en = _mm256_max_epi16(hg, es);
    let gt = _mm256_cmpgt_epi16(hh, rs.bestv);
    let thr;
    if _mm256_movemask_epi8(gt) != 0 {
        // A new best inside this block: exact running best per lane.
        let pfx = _mm256_max_epi16(pmax16(up1_16(hh), c.kb7), rs.bestv);
        thr = _mm256_subs_epu16(pfx, c.xv);
        let rm = _mm256_movemask_epi8(_mm256_cmpgt_epi16(hh, pfx)) as u32;
        rs.rec_b = b + ((31 - rm.leading_zeros()) >> 1) as usize;
        rs.bestv = last16(_mm256_max_epi16(pfx, hh), c.kb7);
        rs.thrv = _mm256_subs_epu16(rs.bestv, c.xv);
        rs.slow = true;
    } else {
        thr = rs.thrv;
    }
    let drop = _mm256_cmpgt_epi16(thr, hh);
    st(hw.add(b), _mm256_andnot_si256(drop, hh));
    st(e.add(b), en);
    let scr = if TB {
        let me = _mm256_cmpgt_epi16(ep, d);
        let mf = _mm256_cmpgt_epi16(f, h0);
        let op = _mm256_andnot_si256(mf, _mm256_add_epi16(c.c3, _mm256_and_si256(me, c.c3)));
        let fb = _mm256_andnot_si256(_mm256_cmpgt_epi16(hg, es), c.c40);
        let fa = _mm256_andnot_si256(_mm256_cmpgt_epi16(hg, fs), c.c10);
        _mm256_or_si256(op, _mm256_or_si256(fb, fa))
    } else {
        _mm256_setzero_si256()
    };
    (drop, scr)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:476-500,563,638,652-673
// ```c
//     score_array[0].best = 0;
//     score_array[0].best_gap = -gap_open_extend;
// ...
//     for (a_index = 1; a_index <= M; a_index++) {
// ...
//         for (b_index = first_b_index; b_index < b_size; b_index++) {
// ...
//         if (first_b_index == b_size || (fence_hit && *fence_hit))
//             break;
// ...
//         if (last_b_index < b_size - 1) {
//             b_size = last_b_index + 1;
//         }
// ...
//         if (b_size <= N) {
//             score_array[b_size].best = MININT;
// ```
// The whole DP on 16-bit lanes, the steps of the 128-bit `align_body` in the same order (row 0,
// then per row: column residues, fence, row table, rebase, blocks, band update). None for a
// problem the lanes do not hold (the caller then runs the 128-bit dispatch).
#[target_feature(enable = "avx2,bmi1,lzcnt,popcnt")]
unsafe fn align16<const TB: bool>(
    q: &RowSeq<'_>,
    s: &ColSeq<'_>,
    scores: &Scores<'_>,
    len1: usize,
    len2: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    check_fence: bool,
    sc: &mut XdropScratch,
) -> Option<XdropResult> {
    sc.ops.clear();
    if len1 == 0 || len2 == 0 {
        return Some(XdropResult {
            a_offset: 0,
            b_offset: 0,
            score: 0,
            fence_hit: false,
            cells: 0,
        });
    }
    let goe = gap_open.checked_add(gap_extend)?;
    let mut x = x_drop;
    if x < goe {
        x = goe;
    }
    if gap_extend < 1 || gap_open < 0 {
        return None;
    }
    let bias = x as i64 + goe as i64 + MAX_M as i64 + 1;
    let limit_rebase = VMAX - bias - MAX_M as i64;
    if limit_rebase < MAX_M as i64 || gap_extend as i64 * 8 > VMAX {
        return None;
    }
    // No NCBI counterpart: the prefix-maximum form of the row gap needs v + 15 * gap_extend to
    // fit 16 bits (v <= 32767).
    if gap_extend > 2048 {
        return None;
    }
    let bias = bias as i32;
    let limit_rebase = limit_rebase as i32;
    let num_extra = (x / gap_extend + 3) as usize;
    sc.reset();

    // Row 0 goes to buffer A (as in `align_body`).
    sc.ensure(num_extra.min(len2 + 2) + 16);
    let mut b_size = 1usize;
    {
        let h = sc.ha.as_mut_ptr().add(PAD);
        let e = sc.e.as_mut_ptr().add(PAD);
        *h = bias as i16;
        *e = (bias - goe) as i16;
        let mut score = -goe;
        for i in 1..=len2 {
            if score < -x {
                break;
            }
            *h.add(i) = (bias + score) as i16;
            *e.add(i) = (bias + score - goe) as i16;
            score -= gap_extend;
            b_size = i + 1;
        }
        // In row 1 the scalar loop never forms the diagonal into column b_size.
        *h.add(b_size - 1) = 0;
    }
    sc.dirty = b_size;
    if TB {
        sc.trace.resize(b_size, SCRIPT_GAP_IN_A);
        sc.rows.push((0, 0, b_size as u32));
    }

    let mut ramp = [0i16; 16];
    for (i, r) in ramp.iter_mut().enumerate() {
        *r = (i as i32 * gap_extend) as i16;
    }
    let c = C16 {
        goe: _mm256_set1_epi16(goe as i16),
        ge1: _mm256_set1_epi16(gap_extend as i16),
        ge16: _mm256_set1_epi16((16 * gap_extend) as i16),
        ramp: ld(ramp.as_ptr()),
        xv: _mm256_set1_epi16(x as i16),
        kb7: ld(KB7X2.as_ptr()),
        c3: _mm256_set1_epi16(SCRIPT_SUB as i16),
        c10: _mm256_set1_epi16(SCRIPT_EXTEND_GAP_A as i16),
        c40: _mm256_set1_epi16(SCRIPT_EXTEND_GAP_B as i16),
    };
    let zero = _mm256_setzero_si256();
    let lane0 = _mm256_setr_epi16(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    let k70 = _mm256_set1_epi8(0x70);
    let k10 = _mm256_set1_epi8(0x10);
    let static_tabs: Option<*const u8> = match scores {
        Scores::Tables { t, .. } => Some(t.data.as_ptr()),
        _ => None,
    };
    let ncols = scores.ncols();
    let mut cells: u64 = 0;

    let mut base: i32 = 0;
    let mut best_rel: i32 = 0;
    let mut a_off = 0usize;
    let mut b_off = 0usize;
    let mut first = 0usize;
    let mut fence_hit = false;
    // exclusive end of the lanes last written in each H buffer
    let mut wend = [b_size, 0usize];
    let mut cur = 0usize; // buffer holding the previous row

    for a in 1..=len1 {
        sc.ensure(b_size + num_extra + 64);
        if sc.sb_filled < b_size {
            fill_cols(s, sc, (b_size + 64).min(len2 + 1), ncols, check_fence);
        }
        if sc.bad_k < b_size {
            return None;
        }
        let mut limit = len2 + 1;
        let mut bend = b_size;
        let mut fence_row = false;
        if sc.fence_k < b_size {
            limit = sc.fence_k;
            bend = sc.fence_k;
            fence_row = true;
            if bend <= first {
                fence_hit = true;
                break;
            }
        }
        let qa = q.get(a)? as usize;
        if qa >= ncols {
            return None;
        }
        cells += (b_size - first) as u64;
        let t: *const u8 = if let Some(tabs) = static_tabs {
            tabs.add(qa * TAB)
        } else {
            if sc.tab_ok & (1u32 << qa) == 0 {
                if !build_row(scores, qa, &mut sc.tabs[qa * TAB..(qa + 1) * TAB]) {
                    return None;
                }
                sc.tab_ok |= 1u32 << qa;
            }
            sc.tabs.as_ptr().add(qa * TAB)
        };
        let (hr, hw) = if cur == 0 {
            (sc.ha.as_mut_ptr().add(PAD), sc.hb.as_mut_ptr().add(PAD))
        } else {
            (sc.hb.as_mut_ptr().add(PAD), sc.ha.as_mut_ptr().add(PAD))
        };
        let e = sc.e.as_mut_ptr().add(PAD);
        let sb = sc.sb.as_ptr().add(PAD);
        let t0lo = bcast128(t);
        let t1lo = bcast128(t.add(16));
        let t0hi = bcast128(t.add(32));
        let t1hi = bcast128(t.add(48));

        if best_rel > limit_rebase {
            let dv = _mm256_set1_epi16(best_rel as i16);
            let mut b = first;
            while b < b_size {
                st(hr.add(b), _mm256_subs_epu16(ld(hr.add(b)), dv));
                st(e.add(b), _mm256_subs_epu16(ld(e.add(b)), dv));
                b += 16;
            }
            base += best_rel;
            best_rel = 0;
        }
        // no diagonal into the first lane of the band: lane 0 of the row's first block (`hpm`);
        // a store of 0 into the H row would sit right under the first block's load
        let mut hpm = lane0;

        let row_off = sc.trace.len();
        if TB {
            sc.trace.reserve(bend - first + num_extra + 96);
        }
        let tr = sc.trace.as_mut_ptr().add(row_off);

        let bestv0 = _mm256_set1_epi16((bias + best_rel) as i16);
        let mut rs = Row {
            carry: zero,
            bestv: bestv0,
            thrv: _mm256_subs_epu16(bestv0, c.xv),
            rec_b: 0,
            slow: false,
        };
        let mut lo_surv = usize::MAX;
        let mut hi_surv = 0usize;
        let mut b0 = first;
        // Two blocks (32 lanes) while more than 16 lanes of the band are left, then single blocks
        // for the rest of the band and the row-gap extension. The stop test runs on the grid
        // first + 16k once the band is passed; lanes processed past the point where the 128-bit
        // kernel stops are computed exactly and dropped (row gap below the threshold), so the
        // result is the same.
        while b0 + 16 < bend {
            // scores of lanes b0..b0+32 (lane b reads column b-1)
            let idx = _mm256_permute4x64_epi64::<0xD8>(ld(sb.add(b0)));
            let i0 = _mm256_add_epi8(idx, k70);
            let i1 = _mm256_sub_epi8(idx, k10);
            let lo = lut32(t0lo, t1lo, i0, i1);
            let hi = lut32(t0hi, t1hi, i0, i1);
            // the first block lies inside the band (b0 + 16 < bend <= limit); the second may
            // reach past the band and past `limit`
            let mask1 = if b0 + 32 <= limit {
                None
            } else {
                Some(ld(LANE_MASK16.as_ptr().add(32 - (limit - b0))))
            };
            let (d0, r0) = block16::<TB>(
                &c,
                &mut rs,
                hr,
                hw,
                e,
                b0,
                _mm256_unpacklo_epi8(lo, hi),
                None,
                hpm,
            );
            hpm = zero;
            let (d1, r1) = block16::<TB>(
                &c,
                &mut rs,
                hr,
                hw,
                e,
                b0 + 16,
                _mm256_unpackhi_epi8(lo, hi),
                mask1,
                zero,
            );
            if TB {
                st(
                    tr.add(b0 - first),
                    _mm256_permute4x64_epi64::<0xD8>(_mm256_packus_epi16(r0, r1)),
                );
            }
            // 2 bits per lane
            let live = !((_mm256_movemask_epi8(d0) as u32 as u64)
                | ((_mm256_movemask_epi8(d1) as u32 as u64) << 32));
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + (live.trailing_zeros() as usize >> 1);
                }
                hi_surv = b0 + ((63 - live.leading_zeros() as usize) >> 1);
            }
            b0 += 32;
        }
        loop {
            if b0 >= bend {
                if b0 >= limit {
                    break;
                }
                if (_mm256_cvtsi256_si32(rs.carry) as i16) < (_mm256_cvtsi256_si32(rs.thrv) as i16)
                {
                    break;
                }
            }
            let idx = _mm256_permute4x64_epi64::<0xD8>(ld(sb.add(b0)));
            let i0 = _mm256_add_epi8(idx, k70);
            let i1 = _mm256_sub_epi8(idx, k10);
            let lo = lut32(t0lo, t1lo, i0, i1);
            let hi = lut32(t0hi, t1hi, i0, i1);
            let mask = if b0 + 16 <= limit {
                None
            } else {
                Some(ld(LANE_MASK16.as_ptr().add(16 - (limit - b0))))
            };
            let (d0, r0) = block16::<TB>(
                &c,
                &mut rs,
                hr,
                hw,
                e,
                b0,
                _mm256_unpacklo_epi8(lo, hi),
                mask,
                hpm,
            );
            hpm = zero;
            if TB {
                stx(
                    tr.add(b0 - first),
                    _mm256_castsi256_si128(_mm256_permute4x64_epi64::<0xD8>(_mm256_packus_epi16(
                        r0, r0,
                    ))),
                );
            }
            let live = !(_mm256_movemask_epi8(d0) as u32);
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + (live.trailing_zeros() as usize >> 1);
                }
                hi_surv = b0 + ((31 - live.leading_zeros() as usize) >> 1);
            }
            b0 += 16;
        }
        let proc_end = b0;
        if proc_end > sc.dirty {
            sc.dirty = proc_end;
        }
        // lanes this buffer held beyond what the row wrote are stale
        let wi = 1 - cur;
        let mut b = proc_end;
        while b < wend[wi] {
            st(hw.add(b), zero);
            b += 16;
        }
        wend[wi] = proc_end;
        cur = wi;
        if rs.slow {
            b_off = rs.rec_b;
            a_off = a;
            best_rel = (_mm256_cvtsi256_si32(rs.bestv) as i16 as i32) - bias;
        }
        if TB {
            sc.trace.set_len(row_off + proc_end - first);
            sc.rows
                .push((row_off as u32, first as u32, (proc_end - first) as u32));
        }
        if fence_row {
            fence_hit = true;
            break;
        }
        if lo_surv == usize::MAX {
            break;
        }
        // columns right of the last survivor leave the band: E = -inf
        let last = hi_surv;
        let mut b = last + 1;
        while b < proc_end {
            st(e.add(b), zero);
            b += 16;
        }
        first = lo_surv;
        b_size = last + 1;
        if b_size <= len2 {
            b_size += 1;
        }
    }

    let best = base + best_rel;
    finish::<TB>(sc, a_off, b_off, best, fence_hit, cells)
}

// ---------------------------------------------------------------------------
// 8-bit lanes, 32 per vector (16 per block in the row tails)
// ---------------------------------------------------------------------------

// No NCBI counterpart: the gap costs and constants of one call on 8-bit lanes.
struct C8 {
    goe: W,
    ge1: W,
    ge16: W,
    ge32: W,
    /// `i * gap_extend` in lane `i`.
    ramp: W,
    xv: W,
    kb15: W,
    c3: W,
    c10: W,
    c40: W,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:578-601,616-631
/// ```c
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
/// ...
///             if (score < score_gap_col) {
///                 script = SCRIPT_GAP_IN_B;
///                 score = score_gap_col;
///             }
///             if (score < score_gap_row) {
///                 script = SCRIPT_GAP_IN_A;
///                 score = score_gap_row;
///             }
///
///             if (best_score - score > x_dropoff) {
/// ...
///                 score_array[b_index].best = score;
/// ```
/// One block of 32 cells on 8-bit lanes: the lane operations of the 128-bit `block8`, row gap as
/// at the head of this file. Returns (drop mask, script).
#[inline(always)]
unsafe fn block8<const TB: bool>(
    c: &C8,
    rs: &mut Row,
    hr: *const u8,
    hw: *mut u8,
    e: *mut u8,
    b: usize,
    s: W,
    lane_mask: Option<W>,
    hpm: W,
) -> (W, W) {
    let hp = _mm256_andnot_si256(hpm, ld(hr.add(b).sub(1)));
    let ep = ld(e.add(b));
    let d = _mm256_adds_epi8(hp, s);
    let h0 = _mm256_max_epi8(d, ep);
    let v = _mm256_subs_epu8(h0, c.goe);
    let fi = _mm256_subs_epu8(pmax8(_mm256_adds_epu8(up1_8(v), c.ramp), c.kb15), c.ramp);
    let f = _mm256_max_epu8(fi, _mm256_subs_epu8(rs.carry, c.ramp));
    let mut hh = _mm256_max_epi8(h0, f);
    if let Some(m) = lane_mask {
        hh = _mm256_and_si256(hh, m);
    }
    let hg = _mm256_subs_epu8(hh, c.goe);
    let fs = _mm256_subs_epu8(f, c.ge1);
    let out = last8(_mm256_max_epu8(_mm256_subs_epu8(fi, c.ge1), v), c.kb15);
    rs.carry = _mm256_max_epu8(_mm256_subs_epu8(rs.carry, c.ge32), out);
    let es = _mm256_subs_epu8(ep, c.ge1);
    let en = _mm256_max_epi8(hg, es);
    let gt = _mm256_cmpgt_epi8(hh, rs.bestv);
    let thr;
    if _mm256_movemask_epi8(gt) != 0 {
        let pfx = _mm256_max_epu8(pmax8(up1_8(hh), c.kb15), rs.bestv);
        thr = _mm256_subs_epu8(pfx, c.xv);
        let rm = _mm256_movemask_epi8(_mm256_cmpgt_epi8(hh, pfx)) as u32;
        rs.rec_b = b + (31 - rm.leading_zeros()) as usize;
        rs.bestv = last8(_mm256_max_epu8(pfx, hh), c.kb15);
        rs.thrv = _mm256_subs_epu8(rs.bestv, c.xv);
        rs.slow = true;
    } else {
        thr = rs.thrv;
    }
    let drop = _mm256_cmpgt_epi8(thr, hh);
    st(hw.add(b), _mm256_andnot_si256(drop, hh));
    st(e.add(b), en);
    let scr = if TB {
        let me = _mm256_cmpgt_epi8(ep, d);
        let mf = _mm256_cmpgt_epi8(f, h0);
        let op = _mm256_andnot_si256(mf, _mm256_add_epi8(c.c3, _mm256_and_si256(me, c.c3)));
        let fb = _mm256_andnot_si256(_mm256_cmpgt_epi8(hg, es), c.c40);
        let fa = _mm256_andnot_si256(_mm256_cmpgt_epi8(hg, fs), c.c10);
        _mm256_or_si256(op, _mm256_or_si256(fb, fa))
    } else {
        _mm256_setzero_si256()
    };
    (drop, scr)
}

// No NCBI counterpart: the row state of `Row` in 128-bit registers, for the 16-lane tail blocks.
struct RowX {
    carry: X,
    bestv: X,
    thrv: X,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:578-601,616-631
/// ```c
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
/// ...
///             if (best_score - score > x_dropoff) {
/// ...
///                 score_array[b_index].best = score;
/// ```
/// `block8` on 16 lanes (the tail of a row and the row-gap extension), with the same values.
#[inline(always)]
unsafe fn block8x<const TB: bool>(
    c: &C8,
    rs: &mut RowX,
    rec_b: &mut usize,
    slow: &mut bool,
    hr: *const u8,
    hw: *mut u8,
    e: *mut u8,
    b: usize,
    s: X,
    lane_mask: Option<X>,
    hpm: X,
) -> (X, X) {
    let goe = _mm256_castsi256_si128(c.goe);
    let ge1 = _mm256_castsi256_si128(c.ge1);
    let ramp = _mm256_castsi256_si128(c.ramp);
    let xv = _mm256_castsi256_si128(c.xv);
    let k15 = _mm256_castsi256_si128(c.kb15);
    let hp = _mm_andnot_si128(hpm, ldx(hr.add(b).sub(1)));
    let ep = ldx(e.add(b));
    let d = _mm_adds_epi8(hp, s);
    let h0 = _mm_max_epi8(d, ep);
    let v = _mm_subs_epu8(h0, goe);
    let fi = _mm_subs_epu8(pmax8x(_mm_adds_epu8(_mm_bslli_si128::<1>(v), ramp)), ramp);
    let f = _mm_max_epu8(fi, _mm_subs_epu8(rs.carry, ramp));
    let mut hh = _mm_max_epi8(h0, f);
    if let Some(m) = lane_mask {
        hh = _mm_and_si128(hh, m);
    }
    let hg = _mm_subs_epu8(hh, goe);
    let fs = _mm_subs_epu8(f, ge1);
    let out = _mm_shuffle_epi8(_mm_max_epu8(_mm_subs_epu8(fi, ge1), v), k15);
    rs.carry = _mm_max_epu8(_mm_subs_epu8(rs.carry, _mm256_castsi256_si128(c.ge16)), out);
    let es = _mm_subs_epu8(ep, ge1);
    let en = _mm_max_epi8(hg, es);
    let gt = _mm_cmpgt_epi8(hh, rs.bestv);
    let thr;
    if _mm_movemask_epi8(gt) != 0 {
        let pfx = _mm_max_epu8(pmax8x(_mm_bslli_si128::<1>(hh)), rs.bestv);
        thr = _mm_subs_epu8(pfx, xv);
        let rm = _mm_movemask_epi8(_mm_cmpgt_epi8(hh, pfx)) as u32;
        *rec_b = b + (31 - rm.leading_zeros()) as usize;
        rs.bestv = _mm_shuffle_epi8(_mm_max_epu8(pfx, hh), k15);
        rs.thrv = _mm_subs_epu8(rs.bestv, xv);
        *slow = true;
    } else {
        thr = rs.thrv;
    }
    let drop = _mm_cmpgt_epi8(thr, hh);
    stx(hw.add(b), _mm_andnot_si128(drop, hh));
    stx(e.add(b), en);
    let scr = if TB {
        let c3 = _mm256_castsi256_si128(c.c3);
        let me = _mm_cmpgt_epi8(ep, d);
        let mf = _mm_cmpgt_epi8(f, h0);
        let op = _mm_andnot_si128(mf, _mm_add_epi8(c3, _mm_and_si128(me, c3)));
        let fb = _mm_andnot_si128(_mm_cmpgt_epi8(hg, es), _mm256_castsi256_si128(c.c40));
        let fa = _mm_andnot_si128(_mm_cmpgt_epi8(hg, fs), _mm256_castsi256_si128(c.c10));
        _mm_or_si128(op, _mm_or_si128(fb, fa))
    } else {
        _mm_setzero_si128()
    };
    (drop, scr)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:476-500,563,638,652-673
// ```c
//     score_array[0].best = 0;
//     score_array[0].best_gap = -gap_open_extend;
// ...
//     for (a_index = 1; a_index <= M; a_index++) {
// ...
//         for (b_index = first_b_index; b_index < b_size; b_index++) {
// ...
//         if (first_b_index == b_size || (fence_hit && *fence_hit))
//             break;
// ...
//         if (last_b_index < b_size - 1) {
//             b_size = last_b_index + 1;
//         }
// ```
// The whole DP on 8-bit lanes, the steps of the 128-bit `align_body8` in the same order. `max_m`
// is the largest score of the matrix (`fits_8bit`). None for a problem the lanes do not hold.
#[target_feature(enable = "avx2,bmi1,lzcnt,popcnt")]
unsafe fn align8<const TB: bool>(
    q: &RowSeq<'_>,
    s: &ColSeq<'_>,
    scores: &Scores<'_>,
    max_m: i32,
    len1: usize,
    len2: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    check_fence: bool,
    sc: &mut XdropScratch,
) -> Option<XdropResult> {
    sc.ops.clear();
    if len1 == 0 || len2 == 0 {
        return Some(XdropResult {
            a_offset: 0,
            b_offset: 0,
            score: 0,
            fence_hit: false,
            cells: 0,
        });
    }
    let goe = gap_open.checked_add(gap_extend)?;
    let mut x = x_drop;
    if x < goe {
        x = goe;
    }
    if gap_extend < 1 || gap_open < 0 || max_m < 0 {
        return None;
    }
    let bias = x as i64 + goe as i64 + max_m as i64 + 1;
    let limit_rebase = 127 - bias - max_m as i64;
    if limit_rebase < 0 {
        return None;
    }
    // No NCBI counterpart: the prefix-maximum form of the row gap needs v + 31 * gap_extend to
    // fit 8 bits (v <= 127).
    if gap_extend > 4 {
        return None;
    }
    let bias = bias as i32;
    let limit_rebase = limit_rebase as i32;
    let num_extra = (x / gap_extend + 3) as usize;
    let ncols = scores.ncols();
    sc.reset();

    // Row 0 goes to buffer A (as in `align_body8`).
    sc.ensure(num_extra.min(len2 + 2) + 16);
    let mut b_size = 1usize;
    {
        let h = sc.h8a.as_mut_ptr().add(PAD);
        let e = sc.e8.as_mut_ptr().add(PAD);
        *h = bias as u8;
        *e = (bias - goe) as u8;
        let mut score = -goe;
        for i in 1..=len2 {
            if score < -x {
                break;
            }
            *h.add(i) = (bias + score) as u8;
            *e.add(i) = (bias + score - goe) as u8;
            score -= gap_extend;
            b_size = i + 1;
        }
        *h.add(b_size - 1) = 0;
    }
    sc.dirty8 = b_size;
    if TB {
        sc.trace.resize(b_size, SCRIPT_GAP_IN_A);
        sc.rows.push((0, 0, b_size as u32));
    }

    let sat = |v: i32| v.min(255) as u8 as i8;
    let mut ramp = [0u8; 32];
    for (i, r) in ramp.iter_mut().enumerate() {
        *r = (i as i32 * gap_extend) as u8;
    }
    let c = C8 {
        goe: _mm256_set1_epi8(sat(goe)),
        ge1: _mm256_set1_epi8(sat(gap_extend)),
        ge16: _mm256_set1_epi8(sat(16 * gap_extend)),
        ge32: _mm256_set1_epi8(sat(32 * gap_extend)),
        ramp: ld(ramp.as_ptr()),
        xv: _mm256_set1_epi8(sat(x)),
        kb15: _mm256_set1_epi8(15),
        c3: _mm256_set1_epi8(SCRIPT_SUB as i8),
        c10: _mm256_set1_epi8(SCRIPT_EXTEND_GAP_A as i8),
        c40: _mm256_set1_epi8(SCRIPT_EXTEND_GAP_B as i8),
    };
    let zero = _mm256_setzero_si256();
    let lane0 = _mm256_setr_epi8(
        -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0,
    );
    let k70 = _mm256_set1_epi8(0x70);
    let k10 = _mm256_set1_epi8(0x10);
    let static_tabs: Option<*const u8> = match scores {
        Scores::Tables { t, .. } => Some(t.data.as_ptr()),
        _ => None,
    };
    let small_alphabet = ncols <= 16;
    let mut cells: u64 = 0;

    let mut base: i32 = 0;
    let mut best_rel: i32 = 0;
    let mut a_off = 0usize;
    let mut b_off = 0usize;
    let mut first = 0usize;
    let mut fence_hit = false;
    let mut wend = [b_size, 0usize];
    let mut cur = 0usize;

    for a in 1..=len1 {
        sc.ensure(b_size + num_extra + 64);
        if sc.sb_filled < b_size {
            fill_cols(s, sc, (b_size + 64).min(len2 + 1), ncols, check_fence);
        }
        if sc.bad_k < b_size {
            return None;
        }
        let mut limit = len2 + 1;
        let mut bend = b_size;
        let mut fence_row = false;
        if sc.fence_k < b_size {
            limit = sc.fence_k;
            bend = sc.fence_k;
            fence_row = true;
            if bend <= first {
                fence_hit = true;
                break;
            }
        }
        let qa = q.get(a)? as usize;
        if qa >= ncols {
            return None;
        }
        cells += (b_size - first) as u64;
        let t: *const u8 = if let Some(tabs) = static_tabs {
            tabs.add(qa * TAB)
        } else {
            if sc.tab_ok & (1u32 << qa) == 0 {
                if !build_row(scores, qa, &mut sc.tabs[qa * TAB..(qa + 1) * TAB]) {
                    return None;
                }
                sc.tab_ok |= 1u32 << qa;
            }
            sc.tabs.as_ptr().add(qa * TAB)
        };
        let (hr, hw) = if cur == 0 {
            (sc.h8a.as_mut_ptr().add(PAD), sc.h8b.as_mut_ptr().add(PAD))
        } else {
            (sc.h8b.as_mut_ptr().add(PAD), sc.h8a.as_mut_ptr().add(PAD))
        };
        let e = sc.e8.as_mut_ptr().add(PAD);
        let sb = sc.sb.as_ptr().add(PAD);
        let t0 = bcast128(t.add(64));
        let t1 = bcast128(t.add(80));

        if best_rel > limit_rebase {
            let dv = _mm256_set1_epi8(best_rel as u8 as i8);
            let mut b = first;
            while b < b_size {
                st(hr.add(b), _mm256_subs_epu8(ld(hr.add(b)), dv));
                st(e.add(b), _mm256_subs_epu8(ld(e.add(b)), dv));
                b += 32;
            }
            base += best_rel;
            best_rel = 0;
        }
        let mut hpm = lane0;

        let row_off = sc.trace.len();
        if TB {
            sc.trace.reserve(bend - first + num_extra + 96);
        }
        let tr = sc.trace.as_mut_ptr().add(row_off);

        let bestv0 = _mm256_set1_epi8((bias + best_rel) as u8 as i8);
        let mut rs = Row {
            carry: zero,
            bestv: bestv0,
            thrv: _mm256_subs_epu8(bestv0, c.xv),
            rec_b: 0,
            slow: false,
        };
        let mut lo_surv = usize::MAX;
        let mut hi_surv = 0usize;
        let mut b0 = first;
        // 32-lane blocks while more than 16 lanes of the band are left, then 16-lane blocks for the
        // rest of the band and the row-gap extension (stop test as in `align16`).
        while b0 + 16 < bend {
            let idx = ld(sb.add(b0));
            let sv = if small_alphabet {
                _mm256_shuffle_epi8(t0, idx)
            } else {
                lut32(t0, t1, _mm256_add_epi8(idx, k70), _mm256_sub_epi8(idx, k10))
            };
            let mask = if b0 + 32 <= limit {
                None
            } else {
                Some(ld(LANE_MASK8.as_ptr().add(32 - (limit - b0))))
            };
            let (d0, r0) = block8::<TB>(&c, &mut rs, hr, hw, e, b0, sv, mask, hpm);
            hpm = zero;
            if TB {
                st(tr.add(b0 - first), r0);
            }
            let live = !(_mm256_movemask_epi8(d0) as u32);
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + live.trailing_zeros() as usize;
                }
                hi_surv = b0 + (31 - live.leading_zeros() as usize);
            }
            b0 += 32;
        }
        let mut rx = RowX {
            carry: _mm256_castsi256_si128(rs.carry),
            bestv: _mm256_castsi256_si128(rs.bestv),
            thrv: _mm256_castsi256_si128(rs.thrv),
        };
        let (mut rec_b, mut slow) = (rs.rec_b, rs.slow);
        let t0x = _mm256_castsi256_si128(t0);
        let t1x = _mm256_castsi256_si128(t1);
        loop {
            if b0 >= bend {
                if b0 >= limit {
                    break;
                }
                if (_mm_cvtsi128_si32(rx.carry) & 0xFF) < (_mm_cvtsi128_si32(rx.thrv) & 0xFF) {
                    break;
                }
            }
            let idx = ldx(sb.add(b0));
            let sv = if small_alphabet {
                _mm_shuffle_epi8(t0x, idx)
            } else {
                lut32x(t0x, t1x, idx)
            };
            let mask = if b0 + 16 <= limit {
                None
            } else {
                Some(ldx(LANE_MASK8X.as_ptr().add(16 - (limit - b0))))
            };
            let (d0, r0) = block8x::<TB>(
                &c,
                &mut rx,
                &mut rec_b,
                &mut slow,
                hr,
                hw,
                e,
                b0,
                sv,
                mask,
                _mm256_castsi256_si128(hpm),
            );
            hpm = zero;
            if TB {
                stx(tr.add(b0 - first), r0);
            }
            let live = (!(_mm_movemask_epi8(d0) as u32)) & 0xFFFF;
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + live.trailing_zeros() as usize;
                }
                hi_surv = b0 + (31 - live.leading_zeros() as usize);
            }
            b0 += 16;
        }
        let proc_end = b0;
        if proc_end > sc.dirty8 {
            sc.dirty8 = proc_end;
        }
        let wi = 1 - cur;
        let mut b = proc_end;
        while b < wend[wi] {
            st(hw.add(b), zero);
            b += 32;
        }
        wend[wi] = proc_end;
        cur = wi;
        if slow {
            b_off = rec_b;
            a_off = a;
            best_rel = (_mm_cvtsi128_si32(rx.bestv) & 0xFF) - bias;
        }
        if TB {
            sc.trace.set_len(row_off + proc_end - first);
            sc.rows
                .push((row_off as u32, first as u32, (proc_end - first) as u32));
        }
        if fence_row {
            fence_hit = true;
            break;
        }
        if lo_surv == usize::MAX {
            break;
        }
        let last = hi_surv;
        let mut b = last + 1;
        while b < proc_end {
            st(e.add(b), zero);
            b += 32;
        }
        first = lo_surv;
        b_size = last + 1;
        if b_size <= len2 {
            b_size += 1;
        }
    }

    let best = base + best_rel;
    finish::<TB>(sc, a_off, b_off, best, fence_hit, cells)
}

#[cfg(test)]
mod tests {
    //! The 256-bit kernels against the scalar transliteration of the NCBI loops and against the
    //! 128-bit kernels: random problems (forward and reversed rows and columns, both lane widths,
    //! traceback and score-only, fences, scale 1 and 32, bands wider than 256 cells), and the
    //! replay of captured problems (`LOSAT_X_DPREPLAY`).
    use super::super::kernel;
    use super::super::tests::{random_problem, reference, Problem, Rng};
    use super::super::x_capture::load_dir;
    use super::*;

    #[derive(Clone, Copy, PartialEq, Eq, Debug)]
    enum Kern {
        V128,
        V256,
    }

    /// One call of the chosen kernel; `m8` selects the 8-bit lanes (as `fits_8bit` would).
    ///
    /// # Safety
    /// `available()` must be true.
    unsafe fn run(
        k: Kern,
        m8: Option<i32>,
        tb: bool,
        q: &RowSeq<'_>,
        s: &ColSeq<'_>,
        scores: &Scores<'_>,
        len1: usize,
        len2: usize,
        go: i32,
        ge: i32,
        x: i32,
        cf: bool,
        sc: &mut XdropScratch,
    ) -> Option<XdropResult> {
        match (k, tb, m8) {
            (Kern::V128, true, Some(m)) => {
                kernel::align8_avx::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (Kern::V128, false, Some(m)) => {
                kernel::align8_avx::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (Kern::V128, true, None) => {
                kernel::align_avx::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
            (Kern::V128, false, None) => {
                kernel::align_avx::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
            (Kern::V256, true, Some(m)) => {
                align8::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (Kern::V256, false, Some(m)) => {
                align8::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (Kern::V256, true, None) => {
                align16::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
            (Kern::V256, false, None) => {
                align16::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
        }
    }

    /// The vector row tables are byte for byte those of `build_row_table`, and are rejected in
    /// the same cases.
    #[test]
    fn row_tables_match_build_row_table() {
        if !available() {
            eprintln!("no AVX2: skipped");
            return;
        }
        let special = [
            i32::MIN,
            i32::MIN / 2,
            -40000,
            -32769,
            -32768,
            -32767,
            -129,
            -128,
            -127,
            -1,
            0,
            1,
            127,
            128,
            1999,
            2000,
            2001,
            32767,
            i32::MAX,
        ];
        let mut rng = Rng(0x5151_7A7A_0101_2323);
        let mut accepted = 0usize;
        for _ in 0..3000 {
            let rare = rng.below(8) == 0;
            let mut pick = |rng: &mut Rng| {
                if rare && rng.below(50) == 0 {
                    special[rng.below(special.len() as u64) as usize]
                } else if rng.below(10) == 0 {
                    special[rng.below(14) as usize]
                } else {
                    rng.below(300) as i32 - 150
                }
            };
            let mut m = [[0i32; 28]; 28];
            for row in m.iter_mut() {
                for v in row.iter_mut() {
                    *v = pick(&mut rng);
                }
            }
            let q = rng.below(30) as usize;
            let (mut a, mut b) = ([0u8; TAB], [0u8; TAB]);
            let rows = Scores::Rows28(&m);
            let ra = build_row_table(&rows, q, &mut a);
            // SAFETY: AVX2 detected above.
            let rb = unsafe { build_row(&rows, q, &mut b) };
            assert_eq!(ra, rb, "Rows28 q={q}");
            if ra {
                assert_eq!(a, b, "Rows28 q={q}");
                accepted += 1;
            }
            let n = [1usize, 5, 8, 9, 16, 20, 28, 31, 32, 33][rng.below(10) as usize];
            let data: Vec<i32> = (0..n * n).map(|_| pick(&mut rng)).collect();
            let q = rng.below(n as u64 + 2) as usize;
            let flat = Scores::Flat { data: &data, n };
            let (mut a, mut b) = ([0u8; TAB], [0u8; TAB]);
            let ra = build_row_table(&flat, q, &mut a);
            // SAFETY: AVX2 detected above.
            let rb = unsafe { build_row(&flat, q, &mut b) };
            assert_eq!(ra, rb, "Flat n={n} q={q}");
            if ra {
                assert_eq!(a, b, "Flat n={n} q={q}");
                accepted += 1;
            }
        }
        assert!(accepted > 1000, "only {accepted} tables accepted");
    }

    /// Repetitive sequences: every diagonal scores alike, so each row keeps the widest band
    /// x_drop allows.
    fn repeat_problem(rng: &mut Rng, n: usize, scale: i32, len: usize) -> Problem {
        let mut p = random_problem(rng, n, scale, 1);
        let letters = (n - 1) as u64;
        let period = 1 + rng.below(3) as usize;
        let motif: Vec<u8> = (0..period).map(|_| 1 + rng.below(letters) as u8).collect();
        p.q = (0..len).map(|i| motif[i % period]).collect();
        p.s = (0..len + rng.below(50) as usize)
            .map(|i| motif[i % period])
            .collect();
        for _ in 0..rng.below(len as u64 / 20 + 1) {
            let k = rng.below(p.s.len() as u64) as usize;
            p.s[k] = 1 + rng.below(letters) as u8;
        }
        p
    }

    fn junk(rng: &mut Rng, k: usize) -> Vec<u8> {
        (0..k)
            .map(|_| match rng.below(4) {
                0 => FENCE_SENTRY,
                1 => 255,
                _ => rng.below(32) as u8,
            })
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:541-557
    // ```c
    //         if (!(gap_align->positionBased)) {
    //             if(reverse_sequence)
    //                 matrix_row = matrix[ A[ M - a_index ] ];
    //             else
    //                 matrix_row = matrix[ A[ a_index ] ];
    // ...
    //         if(reverse_sequence)
    //             b_ptr = &B[N - first_b_index];
    //         else
    //             b_ptr = &B[first_b_index];
    // ```
    // Random problems run through the scalar loops (`reference`), the 128-bit kernels and the
    // 256-bit kernels, with the rows and columns given forward and reversed (as the call sites
    // pass them for the left extension), on 16-bit lanes and, where the problem fits, 8-bit lanes,
    // with and without traceback. Offsets, score, fence flag, edit script and cell count must be
    // equal. `LOSAT_FUZZ_CASES` sets the number of problems.
    #[test]
    fn avx2_kernels_match_scalar_loops_on_random_problems() {
        if !available() {
            eprintln!("no AVX2: skipped");
            return;
        }
        let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(600);
        let mut rng = Rng(0xD1B5_4A32_D192_ED03);
        let (mut sc128, mut sc256) = (XdropScratch::new(), XdropScratch::new());
        // [width 8 / 16][tb][reverse]
        let mut handled = [[[0usize; 2]; 2]; 2];
        let (mut fences, mut ops_total, mut declined, mut wide) = (0usize, 0usize, 0usize, 0usize);
        for case in 0..cases {
            let n = if rng.below(2) == 0 { 16 } else { 28 };
            let scale = if rng.below(3) == 0 { 32 } else { 1 };
            let long = rng.below(6) == 0;
            let len = 1 + rng.below(if long { 2500 } else { 300 }) as usize;
            let repeats = rng.below(4) == 0;
            let p = if repeats {
                repeat_problem(&mut rng, n, scale, len)
            } else {
                random_problem(&mut rng, n, scale, len)
            };
            let (go, ge, x) = match rng.below(8) {
                0 => (11 * scale, scale, 38 * scale),
                1 => (11 * scale, scale, 64 * scale + rng.below(5) as i32),
                2 => (5 * scale, 2 * scale, 33 * scale),
                3 => (5 * scale, 2 * scale, 110 * scale),
                4 => (
                    rng.below(12) as i32 * scale,
                    (1 + rng.below(3) as i32) * scale,
                    (3 + rng.below(150) as i32) * scale,
                ),
                5 => (9 * scale, scale, 3 * scale),
                // wide bands (more than 256 cells on 16-bit lanes)
                6 => (
                    rng.below(12) as i32 * scale,
                    scale,
                    (150 + rng.below(450) as i32) * scale,
                ),
                // small costs for the 8-bit lanes, gap_extend up to 5 (5 is left to 128 bits)
                _ => (
                    rng.below(8) as i32,
                    1 + rng.below(5) as i32,
                    10 + rng.below(80) as i32,
                ),
            };
            let len1 = p.q.len();
            if p.s.is_empty() {
                continue;
            }
            for tb in [false, true] {
                let mut s = p.s.clone();
                let mut len2 = s.len();
                let mut sentinel = false;
                if tb {
                    match rng.below(6) {
                        0 => {
                            s.push(FENCE_SENTRY);
                            sentinel = true;
                        }
                        1 if len2 > 3 => {
                            let k = 1 + rng.below(len2 as u64 - 1) as usize;
                            s[k] = FENCE_SENTRY;
                        }
                        2 if len2 > 3 => {
                            len2 = 1 + rng.below(len2 as u64 - 1) as usize;
                            sentinel = rng.below(2) == 0;
                        }
                        _ => {}
                    }
                }
                let prob = Problem {
                    q: p.q.clone(),
                    s: s.clone(),
                    n,
                    m: p.m.clone(),
                };
                let expected = reference(&prob, len1, len2, go, ge, x, tb, sentinel);
                let scores = Scores::Flat { data: &p.m, n };
                // forward forms, at an offset inside a larger buffer
                let pad_q = rng.below(5) as usize;
                let pad_s = rng.below(5) as usize;
                let mut fq = junk(&mut rng, pad_q);
                fq.extend_from_slice(&p.q);
                let mut fs = junk(&mut rng, pad_s);
                fs.extend_from_slice(&s);
                // reversed forms (`A[M - a_index]`, `B[N - b]`): the columns stop at len2
                let mut rq = junk(&mut rng, pad_q);
                rq.extend(p.q.iter().rev());
                let mut rs = junk(&mut rng, pad_s);
                rs.extend(s[..len2].iter().rev());
                let tail = junk(&mut rng, 3);
                rs.extend_from_slice(&tail);
                for reverse in [false, true] {
                    if reverse && sentinel {
                        continue;
                    }
                    let (rows, cols) = if reverse {
                        (
                            RowSeq::Bytes {
                                data: &rq,
                                base: (pad_q + len1) as isize,
                                step: -1,
                            },
                            ColSeq {
                                data: &rs,
                                base: (pad_s + len2) as isize - 1,
                                step: -1,
                                zero_from: len2,
                            },
                        )
                    } else {
                        (
                            RowSeq::Bytes {
                                data: &fq,
                                base: pad_q as isize - 1,
                                step: 1,
                            },
                            ColSeq {
                                data: &fs,
                                base: pad_s as isize,
                                step: 1,
                                zero_from: if sentinel { usize::MAX } else { len2 },
                            },
                        )
                    };
                    let m8 = fits_8bit(&scores, go, ge, x);
                    for wide8 in [false, true] {
                        let m = if wide8 {
                            match m8 {
                                Some(m) => Some(m),
                                None => continue,
                            }
                        } else {
                            None
                        };
                        let tag = format!(
                            "case {case} tb={tb} rev={reverse} w8={wide8} n={n} scale={scale} \
                             go={go} ge={ge} x={x} len1={len1} len2={len2} rep={repeats}"
                        );
                        // SAFETY: AVX2 detected above.
                        let r128 = unsafe {
                            run(
                                Kern::V128,
                                m,
                                tb,
                                &rows,
                                &cols,
                                &scores,
                                len1,
                                len2,
                                go,
                                ge,
                                x,
                                tb,
                                &mut sc128,
                            )
                        };
                        // SAFETY: AVX2 detected above.
                        let r256 = unsafe {
                            run(
                                Kern::V256,
                                m,
                                tb,
                                &rows,
                                &cols,
                                &scores,
                                len1,
                                len2,
                                go,
                                ge,
                                x,
                                tb,
                                &mut sc256,
                            )
                        };
                        match (r128, r256) {
                            (Some(a), Some(b)) => {
                                assert_eq!(a, b, "{tag}: 128 vs 256");
                                assert_eq!(sc128.ops, sc256.ops, "{tag}: 128 vs 256 ops");
                            }
                            (None, None) => {}
                            (Some(_), None) => declined += 1,
                            (None, Some(_)) => {
                                panic!("{tag}: 256-bit kernel handled a call the 128-bit declines")
                            }
                        }
                        let Some(r) = r256 else {
                            continue;
                        };
                        assert_eq!(
                            (r.a_offset, r.b_offset, r.score, r.fence_hit),
                            (expected.0, expected.1, expected.2, expected.3),
                            "{tag}: 256 vs scalar"
                        );
                        assert_eq!(sc256.ops, expected.4, "{tag}: 256 vs scalar ops");
                        if !r.fence_hit {
                            assert_eq!(r.cells, expected.5, "{tag}: cells");
                        }
                        handled[usize::from(!wide8)][usize::from(tb)][usize::from(reverse)] += 1;
                        fences += usize::from(r.fence_hit);
                        ops_total += sc256.ops.len();
                        if !wide8 && r.cells > 256 * len1.min(len2) as u64 && len1.min(len2) > 500 {
                            wide += 1;
                        }
                    }
                }
            }
        }
        eprintln!(
            "avx2 fuzz: cases={cases} handled[w8/w16][so/tb][fwd/rev]={handled:?} fences={fences} \
             wide={wide} declined={declined}"
        );
        if cases >= 200 {
            for (w, by_tb) in handled.iter().enumerate() {
                for (t, by_dir) in by_tb.iter().enumerate() {
                    for (d, &count) in by_dir.iter().enumerate() {
                        assert!(
                            count > 0,
                            "nothing handled for width {w} tb {t} reverse {d}"
                        );
                    }
                }
            }
            assert!(fences > 0 && ops_total > 0 && wide > 0);
        }
    }

    /// Replay of captured problems (`LOSAT_X_DPCAPTURE` output): every problem through the scalar
    /// loops, the 128-bit and the 256-bit kernels and the `xdrop_align` dispatch, all compared
    /// with the recorded result; then the 128-bit and the 256-bit kernels timed over all
    /// problems (`LOSAT_X_DPREPLAY_REPS` repetitions, default 7, median), per lane width and
    /// traceback / score-only. `LOSAT_X_DPREPLAY_OUT=<file>` appends the table to a file.
    #[test]
    #[ignore = "needs LOSAT_X_DPREPLAY=<capture directory>"]
    fn replay_captured_problems() {
        use std::time::Instant;
        let Some(dir) = std::env::var_os("LOSAT_X_DPREPLAY") else {
            eprintln!("LOSAT_X_DPREPLAY is not set: nothing to replay");
            return;
        };
        assert!(available(), "the replay needs AVX2");
        let dir = std::path::PathBuf::from(dir);
        let probs = load_dir(&dir).expect("read the capture files");
        let reps: usize = std::env::var("LOSAT_X_DPREPLAY_REPS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(7)
            .max(1);
        let (mut sa, mut sb, mut sd) = (
            XdropScratch::new(),
            XdropScratch::new(),
            XdropScratch::new(),
        );
        // categories: [tb][8-bit] -> (problem indices, cells)
        let mut cats: [[(Vec<usize>, u64); 2]; 2] = Default::default();
        // [tb][8-bit]: rows and lanes the 256-bit kernel processes (a traceback run of the problem
        // records them; the band evolves the same way without traceback)
        let mut shape = [[(0u64, 0u64); 2]; 2];
        let (mut unhandled, mut declined) = (0usize, 0usize);
        let mut m8s = Vec::with_capacity(probs.len());
        for (i, p) in probs.iter().enumerate() {
            let (q, s, scores) = (p.row_seq(), p.col_seq(), p.scores());
            let m8 = fits_8bit(&scores, p.gap_open, p.gap_extend, p.x_drop);
            m8s.push(m8);
            let Some(want) = p.result else {
                unhandled += 1;
                continue;
            };
            assert_eq!(p.check_fence, p.tb, "problem {i}: check_fence != tb");
            let args = (p.len1, p.len2, p.gap_open, p.gap_extend, p.x_drop);
            // SAFETY: AVX2 checked above.
            let r128 = unsafe {
                run(
                    Kern::V128,
                    m8,
                    p.tb,
                    &q,
                    &s,
                    &scores,
                    args.0,
                    args.1,
                    args.2,
                    args.3,
                    args.4,
                    p.check_fence,
                    &mut sa,
                )
            };
            assert_eq!(r128, Some(want), "problem {i}: 128-bit vs recorded");
            assert_eq!(sa.ops, p.ops, "problem {i}: 128-bit ops vs recorded");
            // SAFETY: AVX2 checked above.
            let r256 = unsafe {
                run(
                    Kern::V256,
                    m8,
                    p.tb,
                    &q,
                    &s,
                    &scores,
                    args.0,
                    args.1,
                    args.2,
                    args.3,
                    args.4,
                    p.check_fence,
                    &mut sb,
                )
            };
            // SAFETY: AVX2 checked above.
            if unsafe {
                run(
                    Kern::V256,
                    m8,
                    true,
                    &q,
                    &s,
                    &scores,
                    args.0,
                    args.1,
                    args.2,
                    args.3,
                    args.4,
                    p.check_fence,
                    &mut sd,
                )
            }
            .is_some()
            {
                let sh = &mut shape[usize::from(p.tb)][usize::from(m8.is_some())];
                sh.0 += sd.rows.len().saturating_sub(1) as u64;
                sh.1 += (sd.trace.len() - sd.rows.first().map_or(0, |r| r.2 as usize)) as u64;
            }
            match r256 {
                Some(r) => {
                    assert_eq!(r, want, "problem {i}: 256-bit vs recorded");
                    assert_eq!(sb.ops, p.ops, "problem {i}: 256-bit ops vs recorded");
                }
                None => declined += 1,
            }
            assert_eq!(p.replay(&mut sd), Some(want), "problem {i}: dispatch");
            assert_eq!(sd.ops, p.ops, "problem {i}: dispatch ops");
            let prob = Problem {
                q: p.rows.clone(),
                s: p.cols.clone(),
                n: p.n,
                m: p.matrix.clone(),
            };
            let exp = reference(&prob, args.0, args.1, args.2, args.3, args.4, p.tb, true);
            assert_eq!(
                (want.a_offset, want.b_offset, want.score, want.fence_hit),
                (exp.0, exp.1, exp.2, exp.3),
                "problem {i}: recorded vs scalar"
            );
            assert_eq!(p.ops, exp.4, "problem {i}: recorded ops vs scalar");
            if !want.fence_hit {
                assert_eq!(want.cells, exp.5, "problem {i}: cells vs scalar");
            }
            let cat = &mut cats[usize::from(p.tb)][usize::from(m8.is_some())];
            cat.0.push(i);
            cat.1 += want.cells;
        }
        let checked: usize = cats.iter().flatten().map(|c| c.0.len()).sum();
        let mut report = format!(
            "replay {}: problems={} checked={checked} (scalar, 128-bit, 256-bit, dispatch: 0 \
             mismatches) unhandled={unhandled} declined_by_256={declined} reps={reps}\n",
            dir.display(),
            probs.len()
        );
        // One call as `xdrop_align` makes it: a call the 256-bit kernel declines goes to the
        // 128-bit kernel.
        let call = |k: Kern, i: usize, sc: &mut XdropScratch| {
            let p = &probs[i];
            let (q, s, scores) = (p.row_seq(), p.col_seq(), p.scores());
            let go = |k: Kern, sc: &mut XdropScratch| {
                // SAFETY: AVX2 checked above.
                unsafe {
                    run(
                        k,
                        m8s[i],
                        p.tb,
                        &q,
                        &s,
                        &scores,
                        p.len1,
                        p.len2,
                        p.gap_open,
                        p.gap_extend,
                        p.x_drop,
                        p.check_fence,
                        sc,
                    )
                }
            };
            let r = go(k, sc);
            if r.is_none() && k == Kern::V256 {
                go(Kern::V128, sc)
            } else {
                r
            }
        };
        // Reads a problem's residues and matrix, so that a "warm" call finds them in cache as the
        // search's calls do (the matrix and sequences are in use there).
        let touch = |i: usize| {
            let p = &probs[i];
            let mut acc = 0i64;
            for &b in p.rows.iter().chain(p.cols.iter()) {
                acc += b as i64;
            }
            for &v in &p.matrix {
                acc += v as i64;
            }
            if let Scores::Rows28(m) = p.scores() {
                for &v in m.iter().flatten() {
                    acc += v as i64;
                }
            }
            std::hint::black_box(acc);
        };
        // cost of one pair of `Instant::now()` reads, subtracted from the warm per-call times
        let timer_ns = {
            let mut v: Vec<f64> = (0..15)
                .map(|_| {
                    let t0 = Instant::now();
                    for _ in 0..10_000 {
                        std::hint::black_box(Instant::now().elapsed());
                    }
                    t0.elapsed().as_nanos() as f64 / 10_000.0
                })
                .collect();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            v[v.len() / 2]
        };
        report += &format!(
            "  cold = problems back to back in capture order; warm = each call timed alone after its \
             residues and matrix were read (timer cost {timer_ns:.1} ns/call subtracted)\n"
        );
        for (tb, by_w) in cats.iter().enumerate() {
            for (w8, (list, cells)) in by_w.iter().enumerate() {
                if list.is_empty() {
                    continue;
                }
                let (rows, lanes) = shape[tb][w8];
                report += &format!(
                    "  {}-{}: rows per call {:.1}, cells per row {:.1}, lanes processed per cell \
                     {:.3} (256-bit)\n",
                    if tb == 1 { "traceback" } else { "score-only" },
                    if w8 == 1 { 8 } else { 16 },
                    rows as f64 / list.len() as f64,
                    *cells as f64 / rows.max(1) as f64,
                    lanes as f64 / *cells as f64
                );
                for warm in [false, true] {
                    let mut times = [Vec::new(), Vec::new()];
                    // warm: fastest time of each call over the repetitions
                    let mut best = [vec![f64::MAX; list.len()], vec![f64::MAX; list.len()]];
                    for rep in 0..reps {
                        let order = if rep % 2 == 0 {
                            [Kern::V128, Kern::V256]
                        } else {
                            [Kern::V256, Kern::V128]
                        };
                        for k in order {
                            let sc = if k == Kern::V128 { &mut sa } else { &mut sb };
                            let total = if warm {
                                let mut ns = 0f64;
                                let bk = &mut best[usize::from(k == Kern::V256)];
                                for (j, &i) in list.iter().enumerate() {
                                    touch(i);
                                    let t0 = Instant::now();
                                    std::hint::black_box(call(k, i, sc));
                                    let t = t0.elapsed().as_nanos() as f64 - timer_ns;
                                    ns += t;
                                    bk[j] = bk[j].min(t);
                                }
                                ns
                            } else {
                                let t0 = Instant::now();
                                for &i in list {
                                    std::hint::black_box(call(k, i, sc));
                                }
                                t0.elapsed().as_nanos() as f64
                            };
                            times[usize::from(k == Kern::V256)].push(total);
                        }
                    }
                    let stat = |v: &mut Vec<f64>| {
                        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                        (v[v.len() / 2], v[0], v[v.len() - 1])
                    };
                    let (m128, lo128, hi128) = stat(&mut times[0]);
                    let (m256, lo256, hi256) = stat(&mut times[1]);
                    let (nc, nl) = (*cells as f64, list.len() as f64);
                    if warm {
                        let (b128, b256): (f64, f64) = (best[0].iter().sum(), best[1].iter().sum());
                        report += &format!(
                            "  {}-{} warm-min (fastest of {reps} per call, summed): 128-bit \
                             ns/cell={:.3} ns/call={:.0} | 256-bit ns/cell={:.3} ns/call={:.0} | \
                             128/256={:.2}\n",
                            if tb == 1 { "traceback" } else { "score-only" },
                            if w8 == 1 { 8 } else { 16 },
                            b128 / nc,
                            b128 / nl,
                            b256 / nc,
                            b256 / nl,
                            b128 / b256
                        );
                    }
                    report += &format!(
                        "  {}-{} {}: problems={} cells={} | 128-bit ns/cell={:.3} ns/call={:.0} \
                         (min {:.3} max {:.3}) | 256-bit ns/cell={:.3} ns/call={:.0} (min {:.3} \
                         max {:.3}) | 128/256={:.2}\n",
                        if tb == 1 { "traceback" } else { "score-only" },
                        if w8 == 1 { 8 } else { 16 },
                        if warm { "warm" } else { "cold" },
                        list.len(),
                        cells,
                        m128 / nc,
                        m128 / nl,
                        lo128 / nc,
                        hi128 / nc,
                        m256 / nc,
                        m256 / nl,
                        lo256 / nc,
                        hi256 / nc,
                        m128 / m256
                    );
                }
            }
        }
        eprint!("{report}");
        if let Some(out) = std::env::var_os("LOSAT_X_DPREPLAY_OUT") {
            use std::io::Write;
            let mut f = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(out)
                .expect("open LOSAT_X_DPREPLAY_OUT");
            f.write_all(report.as_bytes()).expect("write report");
        }
    }
}
