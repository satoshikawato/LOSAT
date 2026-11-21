//! EXPERIMENT: 128-bit kernels for the X-drop dynamic programming of
//! `xdrop_simd.rs` with the per-row savings of the 256-bit kernels (`x_xdrop_avx2.rs`): 8 lanes
//! of 16 bits and 16 lanes of 8 bits, with (`TB`) and without the traceback. Written once over
//! the 128-bit `arch` layer of `xdrop_simd.rs`, so the same code is built for x86_64 (SSE4.1 and
//! AVX-128 entry points, chosen at run time like the 128-bit kernels), for wasm32 simd128 and for
//! aarch64 NEON, where no 256-bit kernel exists. They run the recurrence of NCBI `ALIGN_EX` /
//! `Blast_SemiGappedAlign` row by row, with the value representation of the 128-bit kernels, and
//! return exactly the `XdropResult` and edit operations of those kernels.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:563-635
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:652-673
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:862-872,903-908
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
//! Unchanged: the 8 / 16-bit choice (`fits_8bit`), the row-start `sc.fits` / `x_room`, the rebase
//! test, the traceback walk (`finish` of the 128-bit kernels) and the cell count.

#![allow(clippy::too_many_arguments)]

use super::arch::*;
use super::kernel::finish;
use super::{
    build_row_table, fits_8bit, script_codes, ColSeq, RowSeq, Scores, XdropResult, XdropScratch,
    FENCE_SENTRY, MAX_M, PAD, TAB, VMAX,
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

// No NCBI counterpart: the same lane operations for aarch64 NEON, on the register type of the
// `arch` layer (`uint8x16_t`).
#[cfg(target_arch = "aarch64")]
mod ext {
    use core::arch::aarch64::*;
    type V = uint8x16_t;
    #[inline(always)]
    pub unsafe fn max_u16(a: V, b: V) -> V {
        vreinterpretq_u8_u16(vmaxq_u16(vreinterpretq_u16_u8(a), vreinterpretq_u16_u8(b)))
    }
    #[inline(always)]
    pub unsafe fn adds_u16(a: V, b: V) -> V {
        vreinterpretq_u8_u16(vqaddq_u16(vreinterpretq_u16_u8(a), vreinterpretq_u16_u8(b)))
    }
    #[inline(always)]
    pub unsafe fn max_u8(a: V, b: V) -> V {
        vmaxq_u8(a, b)
    }
    #[inline(always)]
    pub unsafe fn adds_u8(a: V, b: V) -> V {
        vqaddq_u8(a, b)
    }
    /// Any bit set.
    #[inline(always)]
    pub unsafe fn any(v: V) -> bool {
        vmaxvq_u32(vreinterpretq_u32_u8(v)) != 0
    }
    /// The 16 bytes in reverse order: the bytes of each 64-bit half reversed, then the halves
    /// swapped.
    #[inline(always)]
    pub unsafe fn rev8(v: V) -> V {
        let r = vrev64q_u8(v);
        vextq_u8::<8>(r, r)
    }
    #[inline(always)]
    pub unsafe fn splat32(x: i32) -> V {
        vreinterpretq_u8_s32(vdupq_n_s32(x))
    }
    #[inline(always)]
    pub unsafe fn gt32(a: V, b: V) -> V {
        vreinterpretq_u8_u32(vcgtq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b)))
    }
    /// i32 to i16 with signed saturation: the lanes of `a`, then those of `b`.
    #[inline(always)]
    pub unsafe fn narrow32(a: V, b: V) -> V {
        let lo = vqmovn_s32(vreinterpretq_s32_u8(a));
        vreinterpretq_u8_s16(vqmovn_high_s32(lo, vreinterpretq_s32_u8(b)))
    }
    /// i16 to i8 with signed saturation: the lanes of `a`, then those of `b`.
    #[inline(always)]
    pub unsafe fn narrow16(a: V, b: V) -> V {
        let lo = vqmovn_s16(vreinterpretq_s16_u8(a));
        vreinterpretq_u8_s8(vqmovn_high_s16(lo, vreinterpretq_s16_u8(b)))
    }
    /// The low bytes of the eight 16-bit lanes, then their high bytes: a truncating narrow, then
    /// a narrow of the lanes shifted right by 8 into the high half.
    #[inline(always)]
    pub unsafe fn split16(v: V) -> V {
        let w = vreinterpretq_u16_u8(v);
        vshrn_high_n_u16::<8>(vmovn_u16(w), w)
    }
    /// The low 64 bits of `a`, then the low 64 bits of `b`.
    #[inline(always)]
    pub unsafe fn lo64(a: V, b: V) -> V {
        vreinterpretq_u8_u64(vzip1q_u64(vreinterpretq_u64_u8(a), vreinterpretq_u64_u8(b)))
    }
    /// The high 64 bits of `a`, then the high 64 bits of `b`.
    #[inline(always)]
    pub unsafe fn hi64(a: V, b: V) -> V {
        vreinterpretq_u8_u64(vzip2q_u64(vreinterpretq_u64_u8(a), vreinterpretq_u64_u8(b)))
    }
    /// Bytes moved up by 1 / 2 / 4 / 8 lanes across the vector, zeros shifted in: the window of
    /// `0 ++ v` that starts 16 - n bytes in.
    #[inline(always)]
    pub unsafe fn up1b(v: V) -> V {
        vextq_u8::<15>(vdupq_n_u8(0), v)
    }
    #[inline(always)]
    pub unsafe fn up2b(v: V) -> V {
        vextq_u8::<14>(vdupq_n_u8(0), v)
    }
    #[inline(always)]
    pub unsafe fn up4b(v: V) -> V {
        vextq_u8::<12>(vdupq_n_u8(0), v)
    }
    #[inline(always)]
    pub unsafe fn up8b(v: V) -> V {
        vextq_u8::<8>(vdupq_n_u8(0), v)
    }
}

use ext::*;

// No NCBI counterpart: shuffle and lane masks; they do not change any value NCBI computes.
/// Byte shuffle: 16-bit lane 7 to every lane (x86 `bcast_last16`; wasm and aarch64 ignore it).
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

/// No NCBI counterpart: the 128-bit kernels' CPU check (LOSAT_X_DPNOAVX selects the SSE4.1 build):
/// SSE4.1 on x86_64, always true on wasm32 simd128 and aarch64.
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
/// build on x86_64, as the 128-bit dispatch picks; the simd128 build on wasm32, NEON on aarch64).
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
/// The lookup planes of `build_row_table` (the same bytes, the same rejection of scores above the
/// 16-bit path's range), built with vector instructions for `Rows28` and `Flat` matrices; other
/// sources use `build_row_table` (item 5 of the file header).
#[inline(always)]
pub(super) unsafe fn build_row(scores: &Scores<'_>, q: usize, out: &mut [u8]) -> bool {
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
/// The column residues `s.get(k)` (`*b_ptr`) for `k` in `sc.sb_filled..to` into the lane buffer,
/// with the fence and bad-residue bookkeeping of the 128-bit kernels (first column of each). 16
/// columns at a time while they are a plain forward or backward run of `s.data` with no residue
/// `>= ncols`; the rest, and any block with such a residue, by the scalar loop of the 128-bit
/// kernels (item 4 of the file header).
#[inline(always)]
pub(super) unsafe fn fill_cols(
    s: &ColSeq<'_>,
    sc: &mut XdropScratch,
    to: usize,
    ncols: usize,
    check_fence: bool,
) {
    // `sb` lane k + 1 holds column k + sh (rows over the band, `XdropScratch::x_room`)
    let sh = sc.x_sh;
    let mut k = sc.sb_filled;
    if ncols > 0 && (s.step == 1 || s.step == -1) {
        let lim = splat8((ncols - 1) as i8);
        let dl = s.data.len() as isize;
        while k + 16 <= to && k + sh + 16 <= s.zero_from {
            let v = if s.step == 1 {
                // columns k..k+16 are s.data[i..i + 16]
                let i = s.base + (k + sh) as isize;
                if i < 0 || i + 16 > dl {
                    break;
                }
                load(s.data.as_ptr().offset(i))
            } else {
                // columns k..k+16 are s.data[i], s.data[i - 1], ..., s.data[i - 15]
                let i = s.base - (k + sh) as isize;
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
        let ch = s.get(k + sh);
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
    if !TB || sc.x_nib {
        body16_rows::<TB, TB>(
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
        body16_bytes(
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

// No NCBI counterpart: the traceback of `body16` with one byte per cell (`LOSAT_X_DPNIBBLE=0`). Off
// x86_64 it is a function of its own, so that the default entry keeps one traceback body.
#[cfg_attr(target_arch = "x86_64", inline(always))]
#[cfg_attr(not(target_arch = "x86_64"), inline(never))]
unsafe fn body16_bytes(
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
    body16_rows::<true, false>(
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

// No NCBI counterpart: the kernel of `body16` for one trace layout (`NIB`: four-bit codes,
// `LOSAT_X_DPNIBBLE`).
#[inline(always)]
unsafe fn body16_rows<const TB: bool, const NIB: bool>(
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
        sc.trace.x_row_zero(NIB, b_size);
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
        c3: splat16(const { script_codes(NIB).0 } as i16),
        c10: splat16(const { script_codes(NIB).1 } as i16),
        c40: splat16(const { script_codes(NIB).2 } as i16),
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

    // the rows over the band (`XdropScratch::x_room`): columns below are relative to `sc.x_sh`
    let mut len2 = len2;
    for a in 1..=len1 {
        if !sc.fits(b_size + num_extra + 64) {
            let d = sc.x_room(first, b_size + num_extra + 64, true);
            first -= d;
            b_size -= d;
            len2 -= d;
            wend[0] = wend[0].saturating_sub(d);
            wend[1] = wend[1].saturating_sub(d);
        }
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

        let tr = if TB {
            // SAFETY: room for every byte of the row loop (see `XTrace::x_row`)
            sc.trace.x_row(bend - first + num_extra + 96)
        } else {
            core::ptr::null_mut()
        };

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
            if TB && NIB {
                store64(tr.add((b0 - first) >> 1), nibs(packus16(r0, r1)));
            } else if TB {
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
        let rec = if TB {
            sc.trace.x_row_done(NIB, first + sc.x_sh, proc_end - first)
        } else {
            core::ptr::null()
        };
        if rs.slow {
            b_off = rs.rec_b + sc.x_sh;
            a_off = a;
            best_rel = lane0(rs.bestv) - bias;
            if TB {
                // No NCBI counterpart: the walk starts at this row's record (`XTrace::x_best`)
                sc.trace.x_best(rec);
            }
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
    finish::<TB, NIB>(sc, a_off, b_off, best, fence_hit, cells)
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
    if !TB || sc.x_nib {
        body8_rows::<TB, TB>(
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
        body8_bytes(
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

// No NCBI counterpart: the traceback of `body8` with one byte per cell (`LOSAT_X_DPNIBBLE=0`). Off
// x86_64 it is a function of its own, so that the default entry keeps one traceback body.
#[cfg_attr(target_arch = "x86_64", inline(always))]
#[cfg_attr(not(target_arch = "x86_64"), inline(never))]
unsafe fn body8_bytes(
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
    body8_rows::<true, false>(
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

// No NCBI counterpart: the kernel of `body8` for one trace layout (`NIB`: four-bit codes,
// `LOSAT_X_DPNIBBLE`), with its row loop compiled once for a matrix of at most 16 columns (one
// table lookup per block) and once for a larger one (two), so that no block computes the lookup
// it does not use (wasm32 had both, merged by a select).
#[inline(always)]
unsafe fn body8_rows<const TB: bool, const NIB: bool>(
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
    if scores.ncols() <= 16 {
        body8_rows_t::<TB, NIB, true>(
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
        body8_rows_t::<TB, NIB, false>(
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

// No NCBI counterpart: `body8_rows` for a matrix of at most 16 columns (`SMALL`) or more.
#[inline(always)]
unsafe fn body8_rows_t<const TB: bool, const NIB: bool, const SMALL: bool>(
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
        sc.trace.x_row_zero(NIB, b_size);
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
        // (opaque, as `kb15`: wasm32's engines build a constant vector again at each use, here
        // in every block of the row loop; a vector made once from a value is kept or reloaded)
        c3: splat8(core::hint::black_box(const { script_codes(NIB).0 }) as i8),
        c10: splat8(core::hint::black_box(const { script_codes(NIB).1 }) as i8),
        c40: splat8(core::hint::black_box(const { script_codes(NIB).2 }) as i8),
    };
    let zero_v = zero();
    let hpm0 = load(HPM8.as_ptr());
    let static_tabs: Option<*const u8> = match scores {
        Scores::Tables { t, .. } => Some(t.data.as_ptr()),
        _ => None,
    };
    debug_assert_eq!(SMALL, ncols <= 16);
    let mut cells: u64 = 0;

    let mut base: i32 = 0;
    let mut best_rel: i32 = 0;
    let mut a_off = 0usize;
    let mut b_off = 0usize;
    let mut first = 0usize;
    let mut fence_hit = false;
    let mut wend = [b_size, 0usize];
    let mut cur = 0usize;

    // the rows over the band (`XdropScratch::x_room`): columns below are relative to `sc.x_sh`
    let mut len2 = len2;
    for a in 1..=len1 {
        if !sc.fits(b_size + num_extra + 64) {
            let d = sc.x_room(first, b_size + num_extra + 64, false);
            first -= d;
            b_size -= d;
            len2 -= d;
            wend[0] = wend[0].saturating_sub(d);
            wend[1] = wend[1].saturating_sub(d);
        }
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
        let tr = if TB {
            // SAFETY: room for every byte of the row loop (see `XTrace::x_row`)
            sc.trace.x_row(bend - first + num_extra + 96)
        } else {
            core::ptr::null_mut()
        };

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
        // One block of 16 lanes at `b0` (lane b reads column b-1): `$mask` the lanes of the row's
        // last block up to `limit`, `$hpm` the left edge of the first block. The first block, the
        // blocks inside the band and the last block are written out apart (each with its own
        // constant mask and edge), so the blocks between carry neither from one to the next.
        macro_rules! step {
            ($mask:expr, $hpm:expr) => {{
                let idx = load(sb.add(b0));
                let sv = if SMALL {
                    pick(t0, idx)
                } else {
                    lut32(t0, t1, idx)
                };
                let (d0, r0) = block8::<TB>(&c, &mut rs, hr, hw, e, b0, sv, $mask, $hpm);
                if TB && NIB {
                    store64(tr.add((b0 - first) >> 1), nibs(r0));
                } else if TB {
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
            }};
        }
        let last_mask = |b0: usize| load(LANE_MASK8.as_ptr().add(16 - (limit - b0)));
        if b0 + 16 <= limit {
            step!(None, hpm0);
        } else {
            step!(Some(last_mask(b0)), hpm0);
        }
        loop {
            if b0 >= bend {
                if b0 >= limit {
                    break;
                }
                if lane0_u8(rs.carry) < lane0_u8(rs.thrv) {
                    break;
                }
            }
            if b0 + 16 <= limit {
                step!(None, zero_v);
            } else {
                step!(Some(last_mask(b0)), zero_v);
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
        let rec = if TB {
            sc.trace.x_row_done(NIB, first + sc.x_sh, proc_end - first)
        } else {
            core::ptr::null()
        };
        if rs.slow {
            b_off = rs.rec_b + sc.x_sh;
            a_off = a;
            best_rel = lane0_u8(rs.bestv) - bias;
            if TB {
                // No NCBI counterpart: the walk starts at this row's record (`XTrace::x_best`)
                sc.trace.x_best(rec);
            }
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
    finish::<TB, NIB>(sc, a_off, b_off, best, fence_hit, cells)
}

#[cfg(test)]
#[path = "x_xdrop_w128_tests.rs"]
mod tests;

// No NCBI counterpart: the aarch64 `ext` layer against a scalar model of each operation (the x86
// and wasm definition), on random vectors.
#[cfg(all(test, target_arch = "aarch64"))]
#[path = "x_xdrop_w128_neon_ext_tests.rs"]
mod neon_ext_tests;
