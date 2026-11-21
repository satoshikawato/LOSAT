//! EXPERIMENT: 256-bit AVX2 kernels for the X-drop dynamic programming of
//! `xdrop_simd.rs`: 16 lanes of 16 bits and 32 lanes of 8 bits, with (`TB`) and without the
//! traceback. x86_64 only, chosen at run time when the CPU has AVX2. They run the recurrence of
//! NCBI `ALIGN_EX` / `Blast_SemiGappedAlign` row by row, with the value representation of the
//! 128-bit kernels, and return exactly the `XdropResult` and edit operations of those kernels.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:563-635
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:652-673
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:862-872,903-908
//!
//! Why the result is exact. The argument of `xdrop_simd.rs` holds unchanged: with
//! `thr = best_score - x_dropoff` at the moment a cell is evaluated (never decreasing), a cell
//! that survives has `score >= thr` and its score is decided by candidates `>= thr`; whatever
//! the vector code computes differently from the C loop (gap states behind a dropped cell,
//! saturated values, minus infinity plus a score) is `< thr` in both computations, so it cannot
//! win the `max` of a surviving cell, change which cells are dropped, or change a script bit the
//! traceback walk reads (the walk only visits surviving cells). Each lane here holds exactly the
//! value the 128-bit kernel holds for the same cell, except the values below the threshold that
//! the last item changes:
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
//!   the threshold and nothing above it from the previous row, so they are dropped (H = 0; E 0
//!   from the clearing at the row end, or from the block itself in the 8-bit `EM` instances) and
//!   change nothing. The first and last survivors, hence the band of the next row, the cell
//!   count, the best cell and the traceback rows the walk reads are those of the 128-bit kernel.
//!   The traceback walk is the 128-bit kernels' `finish`.
//! - Row end of the 8-bit lanes (`LOSAT_X_DPEMASK`, the `EM` instances of `align8`; on by
//!   default): each block stores the E of a dropped cell as 0 (minus infinity) instead of
//!   `max(H - gap_open_extend, E - gap_extend)`, and the row no longer clears E right of its last
//!   survivor (it also keeps the high-water mark `dirty8` as a max, the same value). A dropped
//!   cell has `H < thr` and `E <= H`, so the E it passes down is below the threshold in both
//!   forms. What differs downstream is gap states below the threshold, and the gap-extension
//!   script bits that only such gap states decide; the walk reads such a bit only where it enters
//!   a cell through a gap whose state is at least the threshold, so no surviving cell's score,
//!   drop decision or script bit the walk reads changes (the argument above). Every lane right of
//!   the last survivor was dropped in the row that last wrote it (a lane the next row reads but
//!   this row did not write lies right of the band of the row after the one that wrote it), so it
//!   holds 0 without the clearing. Without the clearing the next row's E loads no longer overlap
//!   a store whose address the end of the row decides, which made them wait until that store was
//!   written (line R, R1). The 16-bit lanes keep the clearing: there the store-forward blocks
//!   were 0.1-0.2 per row and the change gained nothing.
//!
//! Less work per row and per call than the 128-bit kernels: wider steps, first and last
//! survivor tracked while the row runs (no survivor bitmap), the first lane's missing diagonal
//! as a lane mask instead of a store under the next load, the column residues copied 32 at a
//! time, row tables of `Rows28` / `Flat` matrices built with vector instructions, and the row
//! gap carried from block to block with two operations.

#![allow(clippy::too_many_arguments)]

use super::kernel::finish;
use super::{
    build_row_table, fits_8bit, script_codes, ColSeq, RowSeq, Scores, XdropResult, XdropScratch,
    FENCE_SENTRY, MAX_M, PAD, TAB, VMAX,
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

/// No NCBI counterpart: LOSAT_X_DPEMASK on (`utils/xswitch.rs`): the `EM` instances of `align8`
/// (E of dropped cells stored as minus infinity, no E clearing at the row end).
pub(super) fn emask_on() -> bool {
    crate::utils::xswitch::on("LOSAT_X_DPEMASK")
}

/// EXPERIMENT (v4-BPF s19 probe, not for dev): `LOSAT_BPF_S19` (decimal bits, read once), the arms
/// of `align16`'s traceback instance (`V`): bit 1 = the diagonal H from the previous block's load
/// (no load at `b - 1`); bit 2 = the scores of the row's band in a pass of their own before the
/// blocks (`srow`); bit 4 = four identity operations on `hh` (a latency slope, not a candidate).
pub(crate) fn bpf_s19() -> u32 {
    static V: std::sync::OnceLock<u32> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("LOSAT_BPF_S19")
            .ok()
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(0)
    })
}

/// No NCBI counterpart: AVX2 and the features of the 128-bit AVX build (LOSAT_X_DPNOAVX, which
/// turns that build off, turns these kernels off too).
pub(super) fn available() -> bool {
    super::cpu_level() == 2 && is_x86_feature_detected!("avx2")
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:766-770
/// Dispatch point of the 256-bit kernels: `tb == false` is the score-only recurrence
/// (`Blast_SemiGappedAlign`, `s_BlastAlignPackedNucl`), `tb == true` is `ALIGN_EX`. The lane width
/// is the one the 128-bit dispatch would choose (`fits_8bit`). `em`: the `EM` instances of
/// `align8` (LOSAT_X_DPEMASK). None: not handled here; the caller runs the 128-bit dispatch.
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
    em: bool,
    sc: &mut XdropScratch,
) -> Option<XdropResult> {
    if em {
        dispatch::<true>(
            q,
            s,
            scores,
            len1,
            len2,
            gap_open,
            gap_extend,
            x_drop,
            tb,
            check_fence,
            sc,
        )
    } else {
        dispatch::<false>(
            q,
            s,
            scores,
            len1,
            len2,
            gap_open,
            gap_extend,
            x_drop,
            tb,
            check_fence,
            sc,
        )
    }
}

/// No NCBI counterpart: `align` for one choice of `EM`.
///
/// # Safety
/// `available()` must be true.
#[inline(always)]
unsafe fn dispatch<const EM: bool>(
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
        (true, Some(m)) if sc.x_nib => align8::<true, EM, true>(
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
        (true, Some(m)) => align8::<true, EM, false>(
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
        (false, Some(m)) => align8::<false, EM, false>(
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
        (true, None) if sc.x_nib => align16::<true, true, 0>(
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
        (true, None) if bpf_s19() & 7 != 0 => {
            let f = match bpf_s19() & 7 {
                1 => align16::<true, false, 1>,
                2 => align16::<true, false, 2>,
                3 => align16::<true, false, 3>,
                4 => align16::<true, false, 4>,
                5 => align16::<true, false, 5>,
                6 => align16::<true, false, 6>,
                _ => align16::<true, false, 7>,
            };
            f(
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
        (true, None) => align16::<true, false, 0>(
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
        (false, None) => align16::<false, false, 0>(
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
/// The low 8 bytes of `v`.
#[inline(always)]
unsafe fn stl<T>(p: *mut T, v: X) {
    _mm_storel_epi64(p as *mut X, v)
}
/// 32 four-bit codes (bytes <= 15) to 16 bytes, byte k = code 2k | code 2k+1 << 4
/// (`LOSAT_X_DPNIBBLE`).
#[inline(always)]
unsafe fn nibs32(v: W) -> X {
    let m = _mm256_maddubs_epi16(v, _mm256_set1_epi16(0x1001));
    _mm_packus_epi16(_mm256_castsi256_si128(m), _mm256_extracti128_si256::<1>(m))
}
/// 16 four-bit codes to 8 bytes in the low half (`nibs32`).
#[inline(always)]
unsafe fn nibs16(v: X) -> X {
    let m = _mm_maddubs_epi16(v, _mm_set1_epi16(0x1001));
    _mm_packus_epi16(m, m)
}
/// EXPERIMENT (s19): `v` through four dependent identity operations (+4 cycles of latency).
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn lat4(mut v: W) -> W {
    core::arch::asm!(
        "vpor {0}, {0}, {0}",
        "vpor {0}, {0}, {0}",
        "vpor {0}, {0}, {0}",
        "vpor {0}, {0}, {0}",
        inout(ymm_reg) v,
        options(pure, nomem, nostack)
    );
    v
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
    // `sb` lane k + 1 holds column k + sh (rows over the band, `XdropScratch::x_room`)
    let sh = sc.x_sh;
    let mut k = sc.sb_filled;
    if ncols > 0 && (s.step == 1 || s.step == -1) {
        let lim = _mm256_set1_epi8((ncols - 1) as i8);
        let rev = _mm256_setr_epi8(
            15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7,
            6, 5, 4, 3, 2, 1, 0,
        );
        let dl = s.data.len() as isize;
        while k + 32 <= to && k + sh + 32 <= s.zero_from {
            let v = if s.step == 1 {
                // columns k..k+32 are s.data[i..i + 32]
                let i = s.base + (k + sh) as isize;
                if i < 0 || i + 32 > dl {
                    break;
                }
                ld(s.data.as_ptr().offset(i))
            } else {
                // columns k..k+32 are s.data[i], s.data[i - 1], ..., s.data[i - 31]
                let i = s.base - (k + sh) as isize;
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
/// One block of 16 cells of one row: `D` = `next_score` of the diagonal, `E` = `score_gap_col`,
/// `F` = `score_gap_row`, `H` = `score`, the same lane operations as the 128-bit `block`, with the
/// row gap as described at the head of this file. Returns (drop mask, script).
#[inline(always)]
unsafe fn block16<const TB: bool, const V: u32>(
    c: &C16,
    rs: &mut Row,
    hr: *const i16,
    hw: *mut i16,
    e: *mut i16,
    b: usize,
    s: W,
    lane_mask: Option<W>,
    hpm: W,
    prev: &mut W,
) -> (W, W) {
    let hp = if V & 1 != 0 {
        // EXPERIMENT (s19): lanes b..b+16 of the H row, moved up one lane with lane b-1 from the
        // previous block's load (`prev`, 0 before the row's first block: no diagonal there)
        let cur = ld(hr.add(b));
        let t = _mm256_permute2x128_si256::<0x03>(cur, *prev);
        *prev = cur;
        _mm256_alignr_epi8::<14>(cur, t)
    } else {
        // `hpm`: lanes with no diagonal (the first lane of the band), as a 0 in the H row
        _mm256_andnot_si256(hpm, ld(hr.add(b).sub(1)))
    };
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
    if V & 4 != 0 {
        // EXPERIMENT (s19): +4 cycles on the block's chain after `hh` (latency slope)
        hh = lat4(hh);
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
// The whole DP on 16-bit lanes, the steps of the 128-bit `align_body` in the same order (row 0,
// then per row: column residues, fence, row table, rebase, blocks, band update). None for a
// problem the lanes do not hold (the caller then runs the 128-bit dispatch).
#[target_feature(enable = "avx2,bmi1,lzcnt,popcnt")]
unsafe fn align16<const TB: bool, const NIB: bool, const V: u32>(
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
    // No NCBI counterpart (memory placement only): see `x_frame_align`.
    crate::utils::x_frame::x_frame_align();
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
        sc.trace.x_row_zero(NIB, b_size);
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
        c3: _mm256_set1_epi16(const { script_codes(NIB).0 } as i16),
        c10: _mm256_set1_epi16(const { script_codes(NIB).1 } as i16),
        c40: _mm256_set1_epi16(const { script_codes(NIB).2 } as i16),
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
    // EXPERIMENT (s19, V & 2): the scores of the band's lanes of one row, lane `first + k` at `k`
    let mut srow: Vec<i16> = Vec::new();

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
        // EXPERIMENT (s19, V & 2): the scores of every block the band loops take (the two blocks
        // of each 32-lane step, and the first single block below `bend`) before the blocks; the
        // blocks past `bend` (the row-gap extension) look theirs up as before
        let sp: *mut i16 = if V & 2 != 0 {
            let need = bend - first + 64;
            if srow.len() < need {
                srow.resize(need.next_power_of_two(), 0);
            }
            let sp = srow.as_mut_ptr();
            let mut k = first;
            while k < bend + 16 {
                let idx = _mm256_permute4x64_epi64::<0xD8>(ld(sb.add(k)));
                let i0 = _mm256_add_epi8(idx, k70);
                let i1 = _mm256_sub_epi8(idx, k10);
                let lo = lut32(t0lo, t1lo, i0, i1);
                let hi = lut32(t0hi, t1hi, i0, i1);
                st(sp.add(k - first), _mm256_unpacklo_epi8(lo, hi));
                st(sp.add(k - first + 16), _mm256_unpackhi_epi8(lo, hi));
                k += 32;
            }
            sp
        } else {
            core::ptr::null_mut()
        };
        // no diagonal into the first lane of the band: lane 0 of the row's first block (`hpm`);
        // a store of 0 into the H row would sit right under the first block's load
        let mut hpm = lane0;
        let mut hprev = zero;

        let tr = if TB {
            // SAFETY: room for every byte of the row loop (see `XTrace::x_row`)
            sc.trace.x_row(bend - first + num_extra + 96)
        } else {
            core::ptr::null_mut()
        };

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
            let (s0, s1) = if V & 2 != 0 {
                (ld(sp.add(b0 - first)), ld(sp.add(b0 - first + 16)))
            } else {
                let idx = _mm256_permute4x64_epi64::<0xD8>(ld(sb.add(b0)));
                let i0 = _mm256_add_epi8(idx, k70);
                let i1 = _mm256_sub_epi8(idx, k10);
                let lo = lut32(t0lo, t1lo, i0, i1);
                let hi = lut32(t0hi, t1hi, i0, i1);
                (_mm256_unpacklo_epi8(lo, hi), _mm256_unpackhi_epi8(lo, hi))
            };
            // the first block lies inside the band (b0 + 16 < bend <= limit); the second may
            // reach past the band and past `limit`
            let mask1 = if b0 + 32 <= limit {
                None
            } else {
                Some(ld(LANE_MASK16.as_ptr().add(32 - (limit - b0))))
            };
            let (d0, r0) = block16::<TB, V>(
                &c,
                &mut rs,
                hr,
                hw,
                e,
                b0,
                s0,
                None,
                hpm,
                &mut hprev,
            );
            hpm = zero;
            let (d1, r1) = block16::<TB, V>(
                &c,
                &mut rs,
                hr,
                hw,
                e,
                b0 + 16,
                s1,
                mask1,
                zero,
                &mut hprev,
            );
            if TB {
                let v = _mm256_permute4x64_epi64::<0xD8>(_mm256_packus_epi16(r0, r1));
                if NIB {
                    stx(tr.add((b0 - first) >> 1), nibs32(v));
                } else {
                    st(tr.add(b0 - first), v);
                }
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
            let s0 = if V & 2 != 0 && b0 < bend {
                ld(sp.add(b0 - first))
            } else {
                let idx = _mm256_permute4x64_epi64::<0xD8>(ld(sb.add(b0)));
                let i0 = _mm256_add_epi8(idx, k70);
                let i1 = _mm256_sub_epi8(idx, k10);
                let lo = lut32(t0lo, t1lo, i0, i1);
                let hi = lut32(t0hi, t1hi, i0, i1);
                _mm256_unpacklo_epi8(lo, hi)
            };
            let mask = if b0 + 16 <= limit {
                None
            } else {
                Some(ld(LANE_MASK16.as_ptr().add(16 - (limit - b0))))
            };
            let (d0, r0) = block16::<TB, V>(
                &c,
                &mut rs,
                hr,
                hw,
                e,
                b0,
                s0,
                mask,
                hpm,
                &mut hprev,
            );
            hpm = zero;
            if TB {
                let v = _mm256_castsi256_si128(_mm256_permute4x64_epi64::<0xD8>(
                    _mm256_packus_epi16(r0, r0),
                ));
                if NIB {
                    stl(tr.add((b0 - first) >> 1), nibs16(v));
                } else {
                    stx(tr.add(b0 - first), v);
                }
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
        let rec = if TB {
            sc.trace.x_row_done(NIB, first + sc.x_sh, proc_end - first)
        } else {
            core::ptr::null()
        };
        if rs.slow {
            b_off = rs.rec_b + sc.x_sh;
            a_off = a;
            best_rel = (_mm256_cvtsi256_si32(rs.bestv) as i16 as i32) - bias;
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
    finish::<TB, NIB>(sc, a_off, b_off, best, fence_hit, cells)
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
/// One block of 32 cells on 8-bit lanes: the lane operations of the 128-bit `block8`, row gap as
/// at the head of this file. `EM`: E of a dropped cell stored as 0. Returns (drop mask, script).
#[inline(always)]
unsafe fn block8<const TB: bool, const EM: bool>(
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
    // LOSAT_X_DPEMASK: the E of a dropped cell as minus infinity
    let en = if EM {
        _mm256_andnot_si256(drop, en)
    } else {
        en
    };
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
/// `block8` on 16 lanes (the tail of a row and the row-gap extension), with the same values.
#[inline(always)]
unsafe fn block8x<const TB: bool, const EM: bool>(
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
    // LOSAT_X_DPEMASK: the E of a dropped cell as minus infinity
    let en = if EM { _mm_andnot_si128(drop, en) } else { en };
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
// The whole DP on 8-bit lanes, the steps of the 128-bit `align_body8` in the same order. `max_m`
// is the largest score of the matrix (`fits_8bit`). None for a problem the lanes do not hold.
#[target_feature(enable = "avx2,bmi1,lzcnt,popcnt")]
unsafe fn align8<const TB: bool, const EM: bool, const NIB: bool>(
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
    // No NCBI counterpart (memory placement only): see `x_frame_align`.
    crate::utils::x_frame::x_frame_align();
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
        sc.trace.x_row_zero(NIB, b_size);
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
        c3: _mm256_set1_epi8(const { script_codes(NIB).0 } as i8),
        c10: _mm256_set1_epi8(const { script_codes(NIB).1 } as i8),
        c40: _mm256_set1_epi8(const { script_codes(NIB).2 } as i8),
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

        let tr = if TB {
            // SAFETY: room for every byte of the row loop (see `XTrace::x_row`)
            sc.trace.x_row(bend - first + num_extra + 96)
        } else {
            core::ptr::null_mut()
        };

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
            let (d0, r0) = block8::<TB, EM>(&c, &mut rs, hr, hw, e, b0, sv, mask, hpm);
            hpm = zero;
            if TB && NIB {
                stx(tr.add((b0 - first) >> 1), nibs32(r0));
            } else if TB {
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
            let (d0, r0) = block8x::<TB, EM>(
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
            if TB && NIB {
                stl(tr.add((b0 - first) >> 1), nibs16(r0));
            } else if TB {
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
        if EM {
            // the high-water mark without a branch (its mispredictions were half of the row's)
            sc.dirty8 = sc.dirty8.max(proc_end);
        } else if proc_end > sc.dirty8 {
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
        let rec = if TB {
            sc.trace.x_row_done(NIB, first + sc.x_sh, proc_end - first)
        } else {
            core::ptr::null()
        };
        if slow {
            b_off = rec_b + sc.x_sh;
            a_off = a;
            best_rel = (_mm_cvtsi128_si32(rx.bestv) & 0xFF) - bias;
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
        if !EM {
            // columns right of the last survivor leave the band: E = -inf
            let mut b = last + 1;
            while b < proc_end {
                st(e.add(b), zero);
                b += 32;
            }
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
#[path = "x_xdrop_avx2_tests.rs"]
mod tests;
