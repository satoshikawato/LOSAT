//! EXPERIMENT (LOSAT_X_DPW128): 128-bit kernels for the X-drop dynamic programming of
//! `xdrop_simd.rs` with the per-row savings of the 256-bit kernels (`x_xdrop_avx2.rs`): 8 lanes
//! of 16 bits and 16 lanes of 8 bits, with (`TB`) and without the traceback. Written once over
//! the 128-bit `arch` layer of `xdrop_simd.rs`, so the same code is built for x86_64 (SSE4.1 and
//! AVX-128 entry points, chosen at run time like the 128-bit kernels) and for wasm32 simd128,
//! where no 256-bit kernel exists. They run the recurrence of NCBI `ALIGN_EX` /
//! `Blast_SemiGappedAlign` row by row, with the value representation of the 128-bit kernels, and
//! return exactly the `XdropResult` and edit operations of those kernels.
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
//! `thr = best_score - x_dropoff` at the moment a cell is evaluated (never decreasing), a cell that
//! survives has `score >= thr` and its score is decided by candidates `>= thr`; whatever the vector
//! code computes differently from the C loop is `< thr` in both computations and cannot win the
//! `max` of a surviving cell, change which cells are dropped, or change a script bit the traceback
//! walk reads. Each lane here holds exactly the value the 128-bit kernel holds for the same cell:
//! the same stored values (`score - base + bias`, 0 = minus infinity), the same `bias`, rebase
//! limit and acceptance tests (plus the precondition of the row gap below, otherwise None), the
//! same blocks (8 lanes of 16 bits or 16 lanes of 8 bits starting at `first_b_index + L*k`, two
//! blocks per step on 16-bit lanes), the same end of a row (the end test on the grid
//! `first_b_index + 16k`), the same diagonal, column gap, X-drop test, running best and script
//! bits (the same lane operations), and lookup tables with the bytes of `build_row_table`. What
//! changes, item by item:
//!
//! 1. First and last survivor. The 128-bit kernel writes a survivor bitmap word per block and
//!    scans it after the row; here the first and the last lane that passed the X-drop test are kept
//!    in two scalars while the row runs (`lo_surv`, `hi_surv`). Both are the first and last set bit
//!    of the same drop masks, so `first_b_index`, `last_b_index` and the new `b_size` are the same.
//! 2. Row gap (`score_gap_row`). Both kernels compute, with saturating arithmetic,
//!    `F[i] = max(max_{j<i} (v[j] - (i-1-j)*gap_extend), carry - i*gap_extend)` with
//!    `v = max(next_score, score_gap_col) - gap_open_extend` and `carry` the row gap entering the
//!    block. The 128-bit kernel uses shifted maxima with decay and carries `max(H - goe, F - g)` of
//!    the last lane. Here `w[j] = v[j] + (j+1)*g` (no saturation: `v + L*g` fits the unsigned lane,
//!    `v <= max - gap_open_extend`, checked, otherwise None), `pw` its inclusive prefix maximum,
//!    `F[i] = sat(max(pw[i-1], carry) - i*g)` (`pw[-1]` = 0) and the next carry
//!    `sat(max(carry, pw[L-1]) - L*g)`. Since `sat(max(a, b) - c) = max(sat(a - c), sat(b - c))`,
//!    `F` is the same value, and the carry is `max(v[L-1], F[L-1] - g)`, which equals
//!    `max(H[L-1] - gap_open_extend, F[L-1] - g)` because `gap_open >= 0`. The carry depends on the
//!    previous one through two operations (a maximum and a subtraction) instead of the whole block.
//!    A block that is cut by a lane mask (the end of `B` or the fence) is the last of its row in
//!    both kernels, so the carry it leaves is never used.
//! 3. The first lane's missing diagonal. The 128-bit kernel stores 0 at `H[first_b_index - 1]` of
//!    the previous row and then loads the block over it (a store under an overlapping load every
//!    row); here the first block of the row clears lane 0 of the loaded `H` with a lane mask
//!    (`hpm`). That location is read by no other lane of this or a later row (bands only move
//!    right), so the store and the mask give every lane the same input.
//! 4. Column residues (`*b_ptr`) are copied 16 at a time while they are a plain forward or backward
//!    run of `B` with no residue `>= ncols` (`fill_cols`); a block holding such a residue, and the
//!    rest, take the 128-bit kernels' scalar loop, which records the first fence and the first bad
//!    column. The lane buffer holds the same bytes.
//! 5. Row tables of `Rows28` / `Flat` matrices are built with vector instructions (`build_row`):
//!    the same rejection of scores above the 16-bit range, `i32 -> i16` saturation (the clamp at
//!    -32768 of `build_row_table`) and `i16 -> i8` saturation (its clamp to -128..=127), so the
//!    same 96 bytes.
//! 6. Lane shifts. The byte shifts across the vector (`up1b` .. `up8b`; the same values as the
//!    `bsl*` of `xdrop_simd.rs`) are written for wasm32 as a window of the concatenation `v ++ 0`,
//!    which V8 compiles to one `palignr`; the form of `bsl*` is compiled as a general two-input
//!    shuffle (two `pshufb` and two masks built in registers). The exclusive prefix maximum of the
//!    running best is formed as the inclusive one moved up by a lane (the same values), so that no
//!    shift is applied to a shifted vector, a pair the compiler merges into one general shuffle.
//!
//! Unchanged: the 8 / 16-bit choice (`fits_8bit`), the per-row `sc.ensure`, the rebase test, the
//! traceback walk (`finish` of the 128-bit kernels) and the cell count.

#![allow(clippy::too_many_arguments)]

use super::arch::*;
use super::kernel::finish;
use super::{
    build_row_table, fits_8bit, ColSeq, RowSeq, Scores, XdropResult, XdropScratch, FENCE_SENTRY,
    MAX_M, PAD, SCRIPT_EXTEND_GAP_A, SCRIPT_EXTEND_GAP_B, SCRIPT_GAP_IN_A, SCRIPT_SUB, TAB, VMAX,
};

// No NCBI counterpart: lane operations the `arch` layer of `xdrop_simd.rs` does not have (unsigned
// maximum and saturating add, any-bit test, byte reversal, 32-bit compare and narrowing, 64-bit
// halves); they do not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
mod ext {
    use core::arch::x86_64::*;
    type V = __m128i;
    #[inline(always)]
    pub unsafe fn max_u16(a: V, b: V) -> V {
        _mm_max_epu16(a, b)
    }
    #[inline(always)]
    pub unsafe fn adds_u16(a: V, b: V) -> V {
        _mm_adds_epu16(a, b)
    }
    #[inline(always)]
    pub unsafe fn max_u8(a: V, b: V) -> V {
        _mm_max_epu8(a, b)
    }
    #[inline(always)]
    pub unsafe fn adds_u8(a: V, b: V) -> V {
        _mm_adds_epu8(a, b)
    }
    /// Any bit set.
    #[inline(always)]
    pub unsafe fn any(v: V) -> bool {
        _mm_testz_si128(v, v) == 0
    }
    /// The 16 bytes in reverse order.
    #[inline(always)]
    pub unsafe fn rev8(v: V) -> V {
        _mm_shuffle_epi8(
            v,
            _mm_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0),
        )
    }
    #[inline(always)]
    pub unsafe fn splat32(x: i32) -> V {
        _mm_set1_epi32(x)
    }
    #[inline(always)]
    pub unsafe fn gt32(a: V, b: V) -> V {
        _mm_cmpgt_epi32(a, b)
    }
    /// i32 to i16 with signed saturation: the lanes of `a`, then those of `b`.
    #[inline(always)]
    pub unsafe fn narrow32(a: V, b: V) -> V {
        _mm_packs_epi32(a, b)
    }
    /// i16 to i8 with signed saturation: the lanes of `a`, then those of `b`.
    #[inline(always)]
    pub unsafe fn narrow16(a: V, b: V) -> V {
        _mm_packs_epi16(a, b)
    }
    /// The low bytes of the eight 16-bit lanes, then their high bytes.
    #[inline(always)]
    pub unsafe fn split16(v: V) -> V {
        _mm_shuffle_epi8(
            v,
            _mm_setr_epi8(0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15),
        )
    }
    /// The low 64 bits of `a`, then the low 64 bits of `b`.
    #[inline(always)]
    pub unsafe fn lo64(a: V, b: V) -> V {
        _mm_unpacklo_epi64(a, b)
    }
    /// The high 64 bits of `a`, then the high 64 bits of `b`.
    #[inline(always)]
    pub unsafe fn hi64(a: V, b: V) -> V {
        _mm_unpackhi_epi64(a, b)
    }
    /// Bytes moved up by 1 / 2 / 4 / 8 lanes across the vector, zeros shifted in (`pslldq`).
    #[inline(always)]
    pub unsafe fn up1b(v: V) -> V {
        _mm_slli_si128::<1>(v)
    }
    #[inline(always)]
    pub unsafe fn up2b(v: V) -> V {
        _mm_slli_si128::<2>(v)
    }
    #[inline(always)]
    pub unsafe fn up4b(v: V) -> V {
        _mm_slli_si128::<4>(v)
    }
    #[inline(always)]
    pub unsafe fn up8b(v: V) -> V {
        _mm_slli_si128::<8>(v)
    }
}

// No NCBI counterpart: the same lane operations for wasm32 simd128.
#[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
mod ext {
    use core::arch::wasm32::*;
    type V = v128;
    #[inline(always)]
    pub unsafe fn max_u16(a: V, b: V) -> V {
        u16x8_max(a, b)
    }
    #[inline(always)]
    pub unsafe fn adds_u16(a: V, b: V) -> V {
        u16x8_add_sat(a, b)
    }
    #[inline(always)]
    pub unsafe fn max_u8(a: V, b: V) -> V {
        u8x16_max(a, b)
    }
    #[inline(always)]
    pub unsafe fn adds_u8(a: V, b: V) -> V {
        u8x16_add_sat(a, b)
    }
    /// Any bit set.
    #[inline(always)]
    pub unsafe fn any(v: V) -> bool {
        v128_any_true(v)
    }
    /// The 16 bytes in reverse order.
    #[inline(always)]
    pub unsafe fn rev8(v: V) -> V {
        i8x16_shuffle::<15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0>(v, v)
    }
    #[inline(always)]
    pub unsafe fn splat32(x: i32) -> V {
        i32x4_splat(x)
    }
    #[inline(always)]
    pub unsafe fn gt32(a: V, b: V) -> V {
        i32x4_gt(a, b)
    }
    /// i32 to i16 with signed saturation: the lanes of `a`, then those of `b`.
    #[inline(always)]
    pub unsafe fn narrow32(a: V, b: V) -> V {
        i16x8_narrow_i32x4(a, b)
    }
    /// i16 to i8 with signed saturation: the lanes of `a`, then those of `b`.
    #[inline(always)]
    pub unsafe fn narrow16(a: V, b: V) -> V {
        i8x16_narrow_i16x8(a, b)
    }
    /// The low bytes of the eight 16-bit lanes, then their high bytes.
    #[inline(always)]
    pub unsafe fn split16(v: V) -> V {
        i8x16_shuffle::<0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15>(v, v)
    }
    /// The low 64 bits of `a`, then the low 64 bits of `b`.
    #[inline(always)]
    pub unsafe fn lo64(a: V, b: V) -> V {
        i64x2_shuffle::<0, 2>(a, b)
    }
    /// The high 64 bits of `a`, then the high 64 bits of `b`.
    #[inline(always)]
    pub unsafe fn hi64(a: V, b: V) -> V {
        i64x2_shuffle::<1, 3>(a, b)
    }
    // Bytes moved up by 1 / 2 / 4 / 8 lanes across the vector, zeros shifted in. Written as a
    // window of the concatenation `v ++ 0` (the zero bytes taken from its end, indices 31, 30, ..),
    // which V8 compiles to one `palignr`; the same shift with the zero bytes taken from the start
    // of the zero vector (16, 17, .., as `bsl*` of `xdrop_simd.rs`) is compiled as a general
    // two-input shuffle (two `pshufb` and two built masks).
    #[inline(always)]
    pub unsafe fn up1b(v: V) -> V {
        i8x16_shuffle::<31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14>(v, u8x16_splat(0))
    }
    #[inline(always)]
    pub unsafe fn up2b(v: V) -> V {
        i8x16_shuffle::<30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13>(v, u8x16_splat(0))
    }
    #[inline(always)]
    pub unsafe fn up4b(v: V) -> V {
        i8x16_shuffle::<28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11>(v, u8x16_splat(0))
    }
    #[inline(always)]
    pub unsafe fn up8b(v: V) -> V {
        i8x16_shuffle::<24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7>(v, u8x16_splat(0))
    }
}

use ext::*;

// No NCBI counterpart: shuffle and lane masks; they do not change any value NCBI computes.
/// Byte shuffle: 16-bit lane 7 to every lane (x86 `bcast_last16`; wasm ignores it).
static KB7: [u8; 16] = [
    14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15,
];
/// Byte shuffle: 8-bit lane 15 to every lane.
static KB15: [u8; 16] = [15; 16];
/// Lane 0 of 16-bit lanes (the first lane of the band: no diagonal).
static HPM16: [i16; 8] = [-1, 0, 0, 0, 0, 0, 0, 0];
/// Lane 0 of 8-bit lanes.
static HPM8: [u8; 16] = [255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
/// `LANE_MASK16[8 - n..]` has the first `n` 16-bit lanes set.
static LANE_MASK16: [i16; 16] = [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0];
/// `LANE_MASK8[16 - n..]` has the first `n` bytes set.
static LANE_MASK8: [u8; 32] = [
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// No NCBI counterpart: LOSAT_X_DPW128 (any value) and a CPU that runs these kernels (SSE4.1 on
/// x86_64, as the 128-bit kernels; always on wasm32 simd128).
pub(super) fn switch_on() -> bool {
    std::env::var_os("LOSAT_X_DPW128").is_some() && available()
}

/// No NCBI counterpart: the 128-bit kernels' CPU check (LOSAT_X_DPNOAVX selects the SSE4.1 build).
pub(super) fn available() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        super::cpu_level() >= 1
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        true
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:766-770
/// ```c
///     if (!score_only) {
///         return ALIGN_EX(A, B, M, N, a_offset, b_offset, edit_block, gap_align,
/// ```
/// Dispatch point of these kernels: `tb == false` is the score-only recurrence
/// (`Blast_SemiGappedAlign`, `s_BlastAlignPackedNucl`), `tb == true` is `ALIGN_EX`. The lane width
/// is the one the 128-bit dispatch chooses (`fits_8bit`). None: not handled here; the caller runs
/// the 128-bit dispatch.
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
        (true, Some(m)) => run8::<true>(
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
        (false, Some(m)) => run8::<false>(
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
        (true, None) => run16::<true>(
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
        (false, None) => run16::<false>(
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

/// No NCBI counterpart: the 16-bit kernel in the build the CPU runs (the AVX-128 or the SSE4.1
/// build on x86_64, as the 128-bit dispatch picks; the simd128 build on wasm32).
///
/// # Safety
/// `available()` must be true.
#[inline]
pub(super) unsafe fn run16<const TB: bool>(
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
    #[cfg(target_arch = "x86_64")]
    {
        if super::cpu_level() == 2 {
            w16_avx::<TB>(
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
            )
        } else {
            w16_sse::<TB>(
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
            )
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        body16::<TB>(
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
        )
    }
}

/// No NCBI counterpart: the 8-bit kernel in the build the CPU runs (see `run16`). `max_m` is the
/// largest matrix score (`fits_8bit`).
///
/// # Safety
/// `available()` must be true.
#[inline]
pub(super) unsafe fn run8<const TB: bool>(
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
    #[cfg(target_arch = "x86_64")]
    {
        if super::cpu_level() == 2 {
            w8_avx::<TB>(
                q,
                s,
                scores,
                max_m,
                len1,
                len2,
                gap_open,
                gap_extend,
                x_drop,
                check_fence,
                sc,
            )
        } else {
            w8_sse::<TB>(
                q,
                s,
                scores,
                max_m,
                len1,
                len2,
                gap_open,
                gap_extend,
                x_drop,
                check_fence,
                sc,
            )
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        body8::<TB>(
            q,
            s,
            scores,
            max_m,
            len1,
            len2,
            gap_open,
            gap_extend,
            x_drop,
            check_fence,
            sc,
        )
    }
}

// No NCBI counterpart: entry points compiled with CPU features enabled (AVX or SSE4.1) around
// `body16` and `body8`, as the 128-bit kernels' `x86_entry`. The same code runs; only the
// instructions the compiler may use change.
macro_rules! x86_entry {
    ($n16:ident, $n8:ident, $feat:literal) => {
        #[cfg(target_arch = "x86_64")]
        #[target_feature(enable = $feat)]
        unsafe fn $n16<const TB: bool>(
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
            body16::<TB>(
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
            )
        }

        #[cfg(target_arch = "x86_64")]
        #[target_feature(enable = $feat)]
        unsafe fn $n8<const TB: bool>(
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
            body8::<TB>(
                q,
                s,
                scores,
                max_m,
                len1,
                len2,
                gap_open,
                gap_extend,
                x_drop,
                check_fence,
                sc,
            )
        }
    };
}
x86_entry!(w16_avx, w8_avx, "avx,ssse3,sse4.1,bmi1,lzcnt,popcnt");
x86_entry!(w16_sse, w8_sse, "ssse3,sse4.1");

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543,578
/// ```c
///                 matrix_row = matrix[ A[ M - a_index ] ];
/// ...
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
/// ```
/// The lookup planes of `build_row_table` (the same bytes, the same rejection of scores above the
/// 16-bit path's range), built with vector instructions for `Rows28` and `Flat` matrices; other
/// sources use `build_row_table` (item 5 of the file header).
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
    // the row, then 0 up to 32 entries (as `build_row_table`)
    let mut tmp = [0i32; 32];
    tmp[..row.len()].copy_from_slice(row);
    let p = tmp.as_ptr();
    let mx = splat32(MAX_M);
    let mut v = [zero(); 8];
    let mut over = zero();
    for (k, slot) in v.iter_mut().enumerate() {
        *slot = load(p.add(4 * k));
        over = or(over, gt32(*slot, mx));
    }
    if any(over) {
        return false;
    }
    // i32 -> i16 saturates below at -32768 (the clamp of `build_row_table`; nothing is above 2000)
    let w0 = narrow32(v[0], v[1]);
    let w1 = narrow32(v[2], v[3]);
    let w2 = narrow32(v[4], v[5]);
    let w3 = narrow32(v[6], v[7]);
    let (x0, x1, x2, x3) = (split16(w0), split16(w1), split16(w2), split16(w3));
    let o = out.as_mut_ptr();
    // low bytes of entries 0..32, high bytes of entries 0..32
    store(o, lo64(x0, x1));
    store(o.add(16), lo64(x2, x3));
    store(o.add(32), hi64(x0, x1));
    store(o.add(48), hi64(x2, x3));
    // i16 -> i8 saturates to -128..=127 (the 8-bit plane)
    store(o.add(64), narrow16(w0, w1));
    store(o.add(80), narrow16(w2, w3));
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
/// with the fence and bad-residue bookkeeping of the 128-bit kernels (first column of each). 16
/// columns at a time while they are a plain forward or backward run of `s.data` with no residue
/// `>= ncols`; the rest, and any block with such a residue, by the scalar loop of the 128-bit
/// kernels (item 4 of the file header).
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
        let lim = splat8((ncols - 1) as i8);
        let dl = s.data.len() as isize;
        while k + 16 <= to && k + 16 <= s.zero_from {
            let v = if s.step == 1 {
                // columns k..k+16 are s.data[i..i + 16]
                let i = s.base + k as isize;
                if i < 0 || i + 16 > dl {
                    break;
                }
                load(s.data.as_ptr().offset(i))
            } else {
                // columns k..k+16 are s.data[i], s.data[i - 1], ..., s.data[i - 15]
                let i = s.base - k as isize;
                if i - 15 < 0 || i >= dl {
                    break;
                }
                rev8(load(s.data.as_ptr().offset(i - 15)))
            };
            // some residue > ncols - 1 (unsigned)
            if any(subsu8(v, lim)) {
                break;
            }
            store(sc.sb.as_mut_ptr().add(PAD + k + 1), v);
            k += 16;
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:601,610-612,616,627
// ```c
//             if (best_score - score > x_dropoff) {
// ...
//                 if (score > best_score) {
//                     best_score = score;
//                     *a_offset = a_index;
// ...
//                 score_gap_row -= gap_extend;
// ...
//                     score_gap_row = score - gap_open_extend;
// ```
// The row gap entering the next block (`score_gap_row`, broadcast), the running `best_score`, the
// threshold `best_score - x_dropoff` and the cell where the best was last raised, carried across
// the blocks of a row.
struct Row {
    carry: V,
    bestv: V,
    thrv: V,
    rec_b: usize,
    slow: bool,
}

// ---------------------------------------------------------------------------
// 16-bit lanes, 8 per block
// ---------------------------------------------------------------------------

// No NCBI counterpart: the gap costs and constants of one call on 16-bit lanes.
struct C16 {
    goe: V,
    ge1: V,
    ge8: V,
    /// `i * gap_extend` in lane `i`.
    ramp: V,
    /// `(i + 1) * gap_extend` in lane `i`.
    ramp1: V,
    xv: V,
    kb7: V,
    c3: V,
    c10: V,
    c40: V,
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
/// One block of 8 cells of one row: `D` = `next_score` of the diagonal, `E` = `score_gap_col`,
/// `F` = `score_gap_row`, `H` = `score`, the lane operations of the 128-bit `block`, with the row
/// gap as in item 2 of the file header and the first lane's diagonal cleared by `hpm` (item 3).
/// Returns (drop mask, script).
#[inline(always)]
unsafe fn block16<const TB: bool>(
    c: &C16,
    rs: &mut Row,
    hr: *const i16,
    hw: *mut i16,
    e: *mut i16,
    b: usize,
    s: V,
    lane_mask: Option<V>,
    hpm: V,
) -> (V, V) {
    let hp = andnot(hpm, load(hr.add(b).sub(1)));
    let ep = load(e.add(b));
    let d = adds16(hp, s);
    let h0 = max16(d, ep);
    let v = subsu16(h0, c.goe);
    // row gap: pw = inclusive prefix maximum of v[j] + (j+1)g (unsigned, exact)
    let w = adds_u16(v, c.ramp1);
    let p = max_u16(w, up2b(w));
    let p = max_u16(p, up4b(p));
    let pw = max_u16(p, up8b(p));
    let f = subsu16(max_u16(up2b(pw), rs.carry), c.ramp);
    rs.carry = subsu16(max_u16(rs.carry, bcast_last16(pw, c.kb7)), c.ge8);
    let mut hh = max16(h0, f);
    if let Some(m) = lane_mask {
        hh = and(hh, m);
    }
    let hg = subsu16(hh, c.goe);
    let es = subsu16(ep, c.ge1);
    let en = max16(hg, es);
    let gt = gt16(hh, rs.bestv);
    let thr;
    if mask8(gt) != 0 {
        // A new best inside this block: exact running best per lane (the exclusive prefix
        // maximum of `H`, formed as the inclusive one moved up by a lane, so that no shift is
        // applied to a shifted vector: such a pair is merged into one general shuffle).
        let y1 = max16(hh, up2b(hh));
        let y2 = max16(y1, up4b(y1));
        let y3 = max16(y2, up8b(y2));
        let pfx = max16(up2b(y3), rs.bestv);
        thr = subsu16(pfx, c.xv);
        let rm = mask8(gt16(hh, pfx));
        rs.rec_b = b + ((31 - rm.leading_zeros()) >> 1) as usize;
        rs.bestv = bcast_last16(max16(pfx, hh), c.kb7);
        rs.thrv = subsu16(rs.bestv, c.xv);
        rs.slow = true;
    } else {
        thr = rs.thrv;
    }
    let drop = gt16(thr, hh);
    store(hw.add(b), andnot(drop, hh));
    store(e.add(b), en);
    let scr = if TB {
        let fs = subsu16(f, c.ge1);
        let me = gt16(ep, d);
        let mf = gt16(f, h0);
        let op = andnot(mf, add16(c.c3, and(me, c.c3)));
        let fb = andnot(gt16(hg, es), c.c40);
        let fa = andnot(gt16(hg, fs), c.c10);
        or(op, or(fb, fa))
    } else {
        zero()
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
#[inline(always)]
unsafe fn body16<const TB: bool>(
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
    // No NCBI counterpart: the prefix-maximum form of the row gap needs v + 8 * gap_extend to fit
    // the unsigned 16-bit lane, with v <= 32767 - gap_open_extend (item 2). The test above already
    // implies it; it is kept here as the stated precondition.
    if VMAX - goe as i64 + 8 * gap_extend as i64 > 65535 {
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

    let mut ramp = [0i16; 8];
    let mut ramp1 = [0i16; 8];
    for (i, (r, r1)) in ramp.iter_mut().zip(ramp1.iter_mut()).enumerate() {
        *r = (i as i32 * gap_extend) as i16;
        *r1 = ((i as i32 + 1) * gap_extend) as i16;
    }
    let c = C16 {
        goe: splat16(goe as i16),
        ge1: splat16(gap_extend as i16),
        ge8: splat16((8 * gap_extend) as i16),
        ramp: load(ramp.as_ptr()),
        ramp1: load(ramp1.as_ptr()),
        xv: splat16(x as i16),
        kb7: load(core::hint::black_box(&KB7).as_ptr()),
        c3: splat16(SCRIPT_SUB as i16),
        c10: splat16(SCRIPT_EXTEND_GAP_A as i16),
        c40: splat16(SCRIPT_EXTEND_GAP_B as i16),
    };
    let zero_v = zero();
    let hpm0 = load(HPM16.as_ptr());
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
        let t0lo = load(t);
        let t1lo = load(t.add(16));
        let t0hi = load(t.add(32));
        let t1hi = load(t.add(48));

        if best_rel > limit_rebase {
            let dv = splat16(best_rel as i16);
            let mut b = first;
            while b < b_size {
                store(hr.add(b), subsu16(load(hr.add(b)), dv));
                store(e.add(b), subsu16(load(e.add(b)), dv));
                b += 8;
            }
            base += best_rel;
            best_rel = 0;
        }
        // no diagonal into the first lane of the band: lane 0 of the row's first block (item 3)
        let mut hpm = hpm0;

        let row_off = sc.trace.len();
        if TB {
            sc.trace.reserve(bend - first + num_extra + 96);
        }
        let tr = sc.trace.as_mut_ptr().add(row_off);

        let bestv0 = splat16((bias + best_rel) as i16);
        let mut rs = Row {
            carry: zero_v,
            bestv: bestv0,
            thrv: subsu16(bestv0, c.xv),
            rec_b: 0,
            slow: false,
        };
        let mut lo_surv = usize::MAX;
        let mut hi_surv = 0usize;
        let mut b0 = first;
        // Two blocks (16 lanes) per step from `first_b_index`, the end test on the grid
        // first + 16k: the steps of the 128-bit kernel.
        loop {
            // scores for the 16 lanes b0..b0+16 (lane b reads column b-1)
            let idx = load(sb.add(b0));
            let lo = lut32(t0lo, t1lo, idx);
            let hi = lut32(t0hi, t1hi, idx);
            let (m0, m1) = if b0 + 16 <= limit {
                (None, None)
            } else {
                let n = limit - b0;
                (
                    Some(load(LANE_MASK16.as_ptr().add(8 - n.min(8)))),
                    Some(load(LANE_MASK16.as_ptr().add(8 - n.saturating_sub(8)))),
                )
            };
            let (d0, r0) = block16::<TB>(&c, &mut rs, hr, hw, e, b0, unpacklo8(lo, hi), m0, hpm);
            hpm = zero_v;
            let (d1, r1) = block16::<TB>(
                &c,
                &mut rs,
                hr,
                hw,
                e,
                b0 + 8,
                unpackhi8(lo, hi),
                m1,
                zero_v,
            );
            if TB {
                store(tr.add(b0 - first), packus16(r0, r1));
            }
            // 2 bits per lane (item 1)
            let live = !(mask8(d0) | (mask8(d1) << 16));
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + (live.trailing_zeros() as usize >> 1);
                }
                hi_surv = b0 + ((31 - live.leading_zeros() as usize) >> 1);
            }
            b0 += 16;
            if b0 >= bend {
                if b0 >= limit {
                    break;
                }
                if lane0(rs.carry) < lane0(rs.thrv) {
                    break;
                }
            }
        }
        let proc_end = b0;
        if proc_end > sc.dirty {
            sc.dirty = proc_end;
        }
        // lanes this buffer held beyond what the row wrote are stale
        let wi = 1 - cur;
        let mut b = proc_end;
        while b < wend[wi] {
            store(hw.add(b), zero_v);
            b += 8;
        }
        wend[wi] = proc_end;
        cur = wi;
        if rs.slow {
            b_off = rs.rec_b;
            a_off = a;
            best_rel = lane0(rs.bestv) - bias;
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
            store(e.add(b), zero_v);
            b += 8;
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
// 8-bit lanes, 16 per block
// ---------------------------------------------------------------------------

// No NCBI counterpart: the gap costs and constants of one call on 8-bit lanes.
struct C8 {
    goe: V,
    ge1: V,
    ge16: V,
    /// `i * gap_extend` in lane `i`.
    ramp: V,
    /// `(i + 1) * gap_extend` in lane `i`.
    ramp1: V,
    xv: V,
    kb15: V,
    c3: V,
    c10: V,
    c40: V,
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
/// One block of 16 cells on 8-bit lanes: the lane operations of the 128-bit `block8`, row gap as
/// in item 2 of the file header, first lane's diagonal cleared by `hpm`. Returns (drop mask,
/// script).
#[inline(always)]
unsafe fn block8<const TB: bool>(
    c: &C8,
    rs: &mut Row,
    hr: *const u8,
    hw: *mut u8,
    e: *mut u8,
    b: usize,
    s: V,
    lane_mask: Option<V>,
    hpm: V,
) -> (V, V) {
    let hp = andnot(hpm, load(hr.add(b).sub(1)));
    let ep = load(e.add(b));
    let d = adds8(hp, s);
    let h0 = max8(d, ep);
    let v = subsu8(h0, c.goe);
    // row gap: pw = inclusive prefix maximum of v[j] + (j+1)g (unsigned, exact)
    let w = adds_u8(v, c.ramp1);
    let p = max_u8(w, up1b(w));
    let p = max_u8(p, up2b(p));
    let p = max_u8(p, up4b(p));
    let pw = max_u8(p, up8b(p));
    let f = subsu8(max_u8(up1b(pw), rs.carry), c.ramp);
    rs.carry = subsu8(max_u8(rs.carry, bcast_last8(pw, c.kb15)), c.ge16);
    let mut hh = max8(h0, f);
    if let Some(m) = lane_mask {
        hh = and(hh, m);
    }
    let hg = subsu8(hh, c.goe);
    let es = subsu8(ep, c.ge1);
    let en = max8(hg, es);
    let gt = gt8(hh, rs.bestv);
    let thr;
    if mask8(gt) != 0 {
        // A new best inside this block: exact running best per lane (exclusive prefix maximum
        // as in `block16`).
        let y1 = max8(hh, up1b(hh));
        let y2 = max8(y1, up2b(y1));
        let y3 = max8(y2, up4b(y2));
        let y4 = max8(y3, up8b(y3));
        let pfx = max8(up1b(y4), rs.bestv);
        thr = subsu8(pfx, c.xv);
        let rm = mask8(gt8(hh, pfx));
        rs.rec_b = b + (31 - rm.leading_zeros()) as usize;
        rs.bestv = bcast_last8(max8(pfx, hh), c.kb15);
        rs.thrv = subsu8(rs.bestv, c.xv);
        rs.slow = true;
    } else {
        thr = rs.thrv;
    }
    let drop = gt8(thr, hh);
    store(hw.add(b), andnot(drop, hh));
    store(e.add(b), en);
    let scr = if TB {
        let fs = subsu8(f, c.ge1);
        let me = gt8(ep, d);
        let mf = gt8(f, h0);
        let op = andnot(mf, add8(c.c3, and(me, c.c3)));
        let fb = andnot(gt8(hg, es), c.c40);
        let fa = andnot(gt8(hg, fs), c.c10);
        or(op, or(fb, fa))
    } else {
        zero()
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
// The whole DP on 8-bit lanes, the steps of the 128-bit `align_body8` in the same order (blocks of
// 16 lanes on the same grid). `max_m` is the largest score of the matrix (`fits_8bit`). None for a
// problem the lanes do not hold.
#[inline(always)]
unsafe fn body8<const TB: bool>(
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
    // No NCBI counterpart: the prefix-maximum form of the row gap needs v + 16 * gap_extend to fit
    // the unsigned 8-bit lane, with v <= 127 - gap_open_extend (item 2).
    if 127 - goe as i64 + 16 * gap_extend as i64 > 255 {
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

    // 16 * gap_extend <= 128 + gap_open_extend <= 254 (checked above and by the bias test)
    let mut ramp = [0u8; 16];
    let mut ramp1 = [0u8; 16];
    for (i, (r, r1)) in ramp.iter_mut().zip(ramp1.iter_mut()).enumerate() {
        *r = (i as i32 * gap_extend) as u8;
        *r1 = ((i as i32 + 1) * gap_extend) as u8;
    }
    let sat = |v: i32| v.min(255) as u8 as i8;
    let c = C8 {
        goe: splat8(sat(goe)),
        ge1: splat8(sat(gap_extend)),
        ge16: splat8(sat(16 * gap_extend)),
        ramp: load(ramp.as_ptr()),
        ramp1: load(ramp1.as_ptr()),
        xv: splat8(sat(x)),
        kb15: load(core::hint::black_box(&KB15).as_ptr()),
        c3: splat8(SCRIPT_SUB as i8),
        c10: splat8(SCRIPT_EXTEND_GAP_A as i8),
        c40: splat8(SCRIPT_EXTEND_GAP_B as i8),
    };
    let zero_v = zero();
    let hpm0 = load(HPM8.as_ptr());
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
        let t0 = load(t.add(64));
        let t1 = load(t.add(80));

        if best_rel > limit_rebase {
            let dv = splat8(best_rel as u8 as i8);
            let mut b = first;
            while b < b_size {
                store(hr.add(b), subsu8(load(hr.add(b)), dv));
                store(e.add(b), subsu8(load(e.add(b)), dv));
                b += 16;
            }
            base += best_rel;
            best_rel = 0;
        }
        let mut hpm = hpm0;

        let row_off = sc.trace.len();
        if TB {
            sc.trace.reserve(bend - first + num_extra + 96);
        }
        let tr = sc.trace.as_mut_ptr().add(row_off);

        let bestv0 = splat8((bias + best_rel) as u8 as i8);
        let mut rs = Row {
            carry: zero_v,
            bestv: bestv0,
            thrv: subsu8(bestv0, c.xv),
            rec_b: 0,
            slow: false,
        };
        let mut lo_surv = usize::MAX;
        let mut hi_surv = 0usize;
        let mut b0 = first;
        loop {
            // scores for the 16 lanes b0..b0+16 (lane b reads column b-1)
            let idx = load(sb.add(b0));
            let sv = if small_alphabet {
                pick(t0, idx)
            } else {
                lut32(t0, t1, idx)
            };
            let mask = if b0 + 16 <= limit {
                None
            } else {
                Some(load(LANE_MASK8.as_ptr().add(16 - (limit - b0))))
            };
            let (d0, r0) = block8::<TB>(&c, &mut rs, hr, hw, e, b0, sv, mask, hpm);
            hpm = zero_v;
            if TB {
                store(tr.add(b0 - first), r0);
            }
            // 1 bit per lane (item 1)
            let live = (!mask8(d0)) & 0xFFFF;
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + live.trailing_zeros() as usize;
                }
                hi_surv = b0 + (31 - live.leading_zeros() as usize);
            }
            b0 += 16;
            if b0 >= bend {
                if b0 >= limit {
                    break;
                }
                if lane0_u8(rs.carry) < lane0_u8(rs.thrv) {
                    break;
                }
            }
        }
        let proc_end = b0;
        if proc_end > sc.dirty8 {
            sc.dirty8 = proc_end;
        }
        let wi = 1 - cur;
        let mut b = proc_end;
        while b < wend[wi] {
            store(hw.add(b), zero_v);
            b += 16;
        }
        wend[wi] = proc_end;
        cur = wi;
        if rs.slow {
            b_off = rs.rec_b;
            a_off = a;
            best_rel = lane0_u8(rs.bestv) - bias;
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
            store(e.add(b), zero_v);
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

#[cfg(test)]
mod tests {
    //! These kernels against the scalar transliteration of the NCBI loops and against the 128-bit
    //! kernels of `xdrop_simd.rs` (random problems: forward and reversed rows and columns, both
    //! lane widths, traceback and score-only, fences, scale 1 and 32, bands wider than 256 cells,
    //! gap costs at the row-gap precondition), the vector row tables and the column copy. They run
    //! natively on the x86_64 build of the same code. The replay of captured problems
    //! (`LOSAT_X_DPREPLAY`) is in `x_xdrop_avx2.rs`.
    use super::super::kernel;
    use super::super::tests::{random_problem, reference, Problem, Rng};
    use super::super::x_capture::load_dir;
    use super::*;

    /// The 128-bit kernel of `xdrop_simd.rs` in the build its dispatch picks; `m8` selects the
    /// 8-bit lanes (as `fits_8bit` would).
    ///
    /// # Safety
    /// `available()` must be true.
    unsafe fn xs128(
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
        #[cfg(target_arch = "x86_64")]
        {
            let avx = super::super::cpu_level() == 2;
            match (avx, tb, m8) {
                (true, true, Some(m)) => {
                    kernel::align8_avx::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
                }
                (true, false, Some(m)) => {
                    kernel::align8_avx::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
                }
                (false, true, Some(m)) => {
                    kernel::align8_sse::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
                }
                (false, false, Some(m)) => {
                    kernel::align8_sse::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
                }
                (true, true, None) => {
                    kernel::align_avx::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
                }
                (true, false, None) => {
                    kernel::align_avx::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
                }
                (false, true, None) => {
                    kernel::align_sse::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
                }
                (false, false, None) => {
                    kernel::align_sse::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
                }
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            match (tb, m8) {
                (true, Some(m)) => {
                    kernel::align_body8::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
                }
                (false, Some(m)) => {
                    kernel::align_body8::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
                }
                (true, None) => {
                    kernel::align_body::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
                }
                (false, None) => {
                    kernel::align_body::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
                }
            }
        }
    }

    /// The kernels of this file; `m8` selects the 8-bit lanes.
    ///
    /// # Safety
    /// `available()` must be true.
    unsafe fn w128(
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
        match (tb, m8) {
            (true, Some(m)) => run8::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc),
            (false, Some(m)) => run8::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc),
            (true, None) => run16::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc),
            (false, None) => run16::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc),
        }
    }

    /// The vector row tables are byte for byte those of `build_row_table`, and are rejected in
    /// the same cases.
    #[test]
    fn w128_row_tables_match_build_row_table() {
        if !available() {
            eprintln!("no SSE4.1: skipped");
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
        let mut rng = Rng(0x5151_7A7A_0101_2324);
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
            // SAFETY: SSE4.1 checked above (or simd128).
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
            // SAFETY: as above.
            let rb = unsafe { build_row(&flat, q, &mut b) };
            assert_eq!(ra, rb, "Flat n={n} q={q}");
            if ra {
                assert_eq!(a, b, "Flat n={n} q={q}");
                accepted += 1;
            }
        }
        assert!(accepted > 1000, "only {accepted} tables accepted");
    }

    /// The column copy fills the lane buffer, the first fence column and the first bad column as
    /// the scalar loop of the 128-bit kernels does, forward and backward, near both ends of the
    /// data and with `zero_from` inside the range.
    #[test]
    fn w128_fill_cols_matches_scalar_loop() {
        if !available() {
            eprintln!("no SSE4.1: skipped");
            return;
        }
        let mut rng = Rng(0x0BAD_C0DE_1234_5678);
        let mut vector_blocks = 0usize;
        for case in 0..4000 {
            let len = rng.below(200) as usize;
            let ncols = [4usize, 16, 25, 28, 32][rng.below(5) as usize];
            let data: Vec<u8> = (0..len)
                .map(|_| match rng.below(60) {
                    0 => FENCE_SENTRY,
                    1 => ncols as u8 + rng.below(20) as u8,
                    2 => 255,
                    _ => rng.below(ncols as u64) as u8,
                })
                .collect();
            let step: isize = if rng.below(2) == 0 { 1 } else { -1 };
            let base = rng.below(len as u64 + 10) as isize - 5;
            let zero_from = if rng.below(3) == 0 {
                rng.below(len as u64 + 20) as usize
            } else {
                usize::MAX
            };
            let s = ColSeq {
                data: &data,
                base,
                step,
                zero_from,
            };
            let check_fence = rng.below(2) == 0;
            let start = rng.below(40) as usize;
            let to = start + rng.below(150) as usize;
            let mut a = XdropScratch::new();
            let mut b = XdropScratch::new();
            a.ensure(to + 64);
            b.ensure(to + 64);
            for sc in [&mut a, &mut b] {
                sc.sb_filled = start;
                sc.fence_k = usize::MAX;
                sc.bad_k = usize::MAX;
            }
            // the scalar loop of `align_body` / `align_body8`
            for k in a.sb_filled..to {
                let ch = s.get(k);
                if ch as usize >= ncols {
                    if check_fence && ch == FENCE_SENTRY {
                        if a.fence_k == usize::MAX {
                            a.fence_k = k;
                        }
                    } else if a.bad_k == usize::MAX {
                        a.bad_k = k;
                    }
                }
                a.sb[PAD + k + 1] = ch;
            }
            a.sb_filled = to;
            // SAFETY: SSE4.1 checked above (or simd128); the buffer holds `to + 64` lanes.
            unsafe { fill_cols(&s, &mut b, to, ncols, check_fence) };
            assert_eq!(
                (a.sb_filled, a.fence_k, a.bad_k),
                (b.sb_filled, b.fence_k, b.bad_k),
                "case {case}"
            );
            assert_eq!(a.sb, b.sb, "case {case}");
            if (start..to).filter(|k| s.get(*k) < ncols as u8).count() >= 16 {
                vector_blocks += 1;
            }
        }
        assert!(vector_blocks > 500, "{vector_blocks}");
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
    // Random problems run through the scalar loops (`reference`), the 128-bit kernels of
    // `xdrop_simd.rs` and the kernels of this file, with the rows and columns given forward and
    // reversed (as the call sites pass them for the left extension), on 16-bit lanes and, where
    // the problem fits, 8-bit lanes, with and without traceback. Offsets, score, fence flag, edit
    // script and cell count must be equal; a call these kernels decline must be one that fails the
    // row-gap precondition. `LOSAT_FUZZ_CASES` sets the number of problems.
    #[test]
    fn w128_kernels_match_scalar_loops_on_random_problems() {
        if !available() {
            eprintln!("no SSE4.1: skipped");
            return;
        }
        let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(600);
        let mut rng = Rng(0x7F4A_7C15_D1B5_4A33);
        let (mut sc128, mut scw) = (XdropScratch::new(), XdropScratch::new());
        // [width 8 / 16][tb][reverse]
        let mut handled = [[[0usize; 2]; 2]; 2];
        let (mut fences, mut ops_total, mut wide) = (0usize, 0usize, 0usize);
        // 8-bit calls with 16 * gap_extend - gap_open_extend in 113..=128 (handled) and above 128
        // (declined, run by the 128-bit kernel)
        let (mut edge_in, mut edge_out, mut big_ge) = (0usize, 0usize, 0usize);
        for case in 0..cases {
            let n = if rng.below(2) == 0 { 16 } else { 28 };
            let mut scale = if rng.below(3) == 0 { 32 } else { 1 };
            let long = rng.below(6) == 0;
            let len = 1 + rng.below(if long { 2500 } else { 300 }) as usize;
            let kind = rng.below(10);
            if kind == 8 {
                scale = 1;
            }
            let repeats = rng.below(4) == 0;
            let p = if repeats {
                repeat_problem(&mut rng, n, scale, len)
            } else {
                random_problem(&mut rng, n, scale, len)
            };
            let (go, ge, x) = match kind {
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
                // small costs for the 8-bit lanes
                7 => (
                    rng.below(8) as i32,
                    1 + rng.below(5) as i32,
                    10 + rng.below(80) as i32,
                ),
                // 8-bit lanes at the row-gap precondition (16 * ge - goe around 128)
                8 => {
                    let ge = 8 + rng.below(4) as i32;
                    let go = rng.below(45) as i32;
                    (go, ge, go + ge + rng.below(30) as i32)
                }
                // 16-bit lanes with gap_extend up to the 128-bit kernel's limit (8 * ge <= 32767)
                _ => {
                    let ge = 1000 + rng.below(3096) as i32;
                    let go = rng.below(2000) as i32;
                    (go, ge, go + ge + rng.below(4000) as i32)
                }
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
                        // SAFETY: SSE4.1 checked above (or simd128).
                        let (r128, rw) = unsafe {
                            (
                                xs128(
                                    m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb,
                                    &mut sc128,
                                ),
                                w128(
                                    m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb,
                                    &mut scw,
                                ),
                            )
                        };
                        let edge = 16 * ge - (go + ge);
                        match (r128, rw) {
                            (Some(a), Some(b)) => {
                                assert_eq!(a, b, "{tag}: 128 vs w128");
                                assert_eq!(sc128.ops, scw.ops, "{tag}: 128 vs w128 ops");
                            }
                            (None, None) => {}
                            (Some(_), None) => {
                                assert!(wide8 && edge > 128, "{tag}: w128 declined");
                                edge_out += 1;
                            }
                            (None, Some(_)) => {
                                panic!("{tag}: w128 handled a call the 128-bit kernel declines")
                            }
                        }
                        let Some(r) = rw else {
                            continue;
                        };
                        assert_eq!(
                            (r.a_offset, r.b_offset, r.score, r.fence_hit),
                            (expected.0, expected.1, expected.2, expected.3),
                            "{tag}: w128 vs scalar"
                        );
                        assert_eq!(scw.ops, expected.4, "{tag}: w128 vs scalar ops");
                        if !r.fence_hit {
                            assert_eq!(r.cells, expected.5, "{tag}: cells");
                        }
                        handled[usize::from(!wide8)][usize::from(tb)][usize::from(reverse)] += 1;
                        fences += usize::from(r.fence_hit);
                        ops_total += scw.ops.len();
                        if wide8 && (113..=128).contains(&edge) {
                            edge_in += 1;
                        }
                        if !wide8 && ge >= 1000 {
                            big_ge += 1;
                        }
                        if !wide8 && r.cells > 256 * len1.min(len2) as u64 && len1.min(len2) > 500 {
                            wide += 1;
                        }
                    }
                }
            }
        }
        eprintln!(
            "w128 fuzz: cases={cases} handled[w8/w16][so/tb][fwd/rev]={handled:?} fences={fences} \
             wide={wide} edge_in={edge_in} edge_out={edge_out} big_ge={big_ge}"
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
            assert!(edge_in > 0 && edge_out > 0 && big_ge > 0);
        }
    }

    /// Replay of captured problems (`LOSAT_X_DPCAPTURE` output) through the 128-bit kernel of
    /// `xdrop_simd.rs` and the kernels of this file, both compared with the recorded result, then
    /// both timed over all problems per lane width and traceback / score-only: "cold" = the
    /// problems back to back in capture order (median of `LOSAT_X_DPREPLAY_REPS` repetitions,
    /// default 7), "warm-min" = each call timed alone after its residues and matrix were read,
    /// fastest of the repetitions, summed. It runs on every target with the 128-bit kernels, the
    /// wasm32 simd128 build included (the replay in `x_xdrop_avx2.rs` needs AVX2).
    /// `LOSAT_X_DPREPLAY_OUT=<file>` appends the table to a file.
    #[test]
    #[ignore = "needs LOSAT_X_DPREPLAY=<capture directory>"]
    fn w128_replay_captured_problems() {
        use std::time::Instant;
        let Some(dir) = std::env::var_os("LOSAT_X_DPREPLAY") else {
            eprintln!("LOSAT_X_DPREPLAY is not set: nothing to replay");
            return;
        };
        assert!(available(), "the replay needs the 128-bit kernels");
        let dir = std::path::PathBuf::from(dir);
        let probs = load_dir(&dir).expect("read the capture files");
        let reps: usize = std::env::var("LOSAT_X_DPREPLAY_REPS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(7)
            .max(1);
        let (mut sa, mut sw) = (XdropScratch::new(), XdropScratch::new());
        // categories: [tb][8-bit] -> (problem indices, cells)
        let mut cats: [[(Vec<usize>, u64); 2]; 2] = Default::default();
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
            let args = (p.len1, p.len2, p.gap_open, p.gap_extend, p.x_drop);
            // SAFETY: `available()` checked above.
            let (r128, rw) = unsafe {
                (
                    xs128(
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
                    ),
                    w128(
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
                        &mut sw,
                    ),
                )
            };
            assert_eq!(r128, Some(want), "problem {i}: 128-bit vs recorded");
            assert_eq!(sa.ops, p.ops, "problem {i}: 128-bit ops vs recorded");
            match rw {
                Some(r) => {
                    assert_eq!(r, want, "problem {i}: w128 vs recorded");
                    assert_eq!(sw.ops, p.ops, "problem {i}: w128 ops vs recorded");
                }
                None => declined += 1,
            }
            let cat = &mut cats[usize::from(p.tb)][usize::from(m8.is_some())];
            cat.0.push(i);
            cat.1 += want.cells;
        }
        let checked: usize = cats.iter().flatten().map(|c| c.0.len()).sum();
        let mut report = format!(
            "w128 replay {} ({}): problems={} checked={checked} (128-bit, w128: 0 mismatches) \
             unhandled={unhandled} declined_by_w128={declined} reps={reps}\n",
            dir.display(),
            std::env::consts::ARCH,
            probs.len()
        );
        // One call as `xdrop_align` makes it: a call w128 declines goes to the 128-bit kernel.
        let call = |w: bool, i: usize, sc: &mut XdropScratch| {
            let p = &probs[i];
            let (q, s, scores) = (p.row_seq(), p.col_seq(), p.scores());
            let a = (
                p.len1,
                p.len2,
                p.gap_open,
                p.gap_extend,
                p.x_drop,
                p.check_fence,
            );
            // SAFETY: `available()` checked above.
            unsafe {
                let r = if w {
                    w128(
                        m8s[i], p.tb, &q, &s, &scores, a.0, a.1, a.2, a.3, a.4, a.5, sc,
                    )
                } else {
                    None
                };
                if r.is_none() {
                    xs128(
                        m8s[i], p.tb, &q, &s, &scores, a.0, a.1, a.2, a.3, a.4, a.5, sc,
                    )
                } else {
                    r
                }
            }
        };
        let touch = |i: usize| {
            let p = &probs[i];
            let mut acc = 0i64;
            for &b in p.rows.iter().chain(p.cols.iter()) {
                acc += b as i64;
            }
            for &v in &p.matrix {
                acc += v as i64;
            }
            std::hint::black_box(acc);
        };
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
        report += &format!("  timer cost {timer_ns:.1} ns/call (subtracted from warm calls)\n");
        for (tb, by_w) in cats.iter().enumerate() {
            for (w8, (list, cells)) in by_w.iter().enumerate() {
                if list.is_empty() {
                    continue;
                }
                let cat = format!(
                    "{}-{}",
                    if tb == 1 { "traceback" } else { "score-only" },
                    if w8 == 1 { 8 } else { 16 }
                );
                let mut cold = [Vec::new(), Vec::new()];
                let mut best = [vec![f64::MAX; list.len()], vec![f64::MAX; list.len()]];
                for rep in 0..reps {
                    for slot in 0..2 {
                        let w = (slot + rep) % 2 == 1;
                        let sc = if w { &mut sw } else { &mut sa };
                        let t0 = Instant::now();
                        for &i in list {
                            std::hint::black_box(call(w, i, sc));
                        }
                        cold[usize::from(w)].push(t0.elapsed().as_nanos() as f64);
                        let bk = &mut best[usize::from(w)];
                        for (j, &i) in list.iter().enumerate() {
                            touch(i);
                            let t0 = Instant::now();
                            std::hint::black_box(call(w, i, sc));
                            bk[j] = bk[j].min(t0.elapsed().as_nanos() as f64 - timer_ns);
                        }
                    }
                }
                let (nc, nl) = (*cells as f64, list.len() as f64);
                let med = |v: &mut Vec<f64>| {
                    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    v[v.len() / 2]
                };
                let (c128, cw) = (med(&mut cold[0]), med(&mut cold[1]));
                let (b128, bw): (f64, f64) = (best[0].iter().sum(), best[1].iter().sum());
                report += &format!(
                    "  {cat}: problems={} cells={cells} | cold 128-bit ns/cell={:.3} ns/call={:.0} \
                     w128 ns/cell={:.3} ns/call={:.0} 128/w128={:.2} | warm-min 128-bit \
                     ns/cell={:.3} ns/call={:.0} w128 ns/cell={:.3} ns/call={:.0} 128/w128={:.2}\n",
                    list.len(),
                    c128 / nc,
                    c128 / nl,
                    cw / nc,
                    cw / nl,
                    c128 / cw,
                    b128 / nc,
                    b128 / nl,
                    bw / nc,
                    bw / nl,
                    b128 / bw
                );
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
