//! EXPERIMENT (LOSAT_X_DPNEON2): the 256-bit kernels of `x_xdrop_avx2.rs` on aarch64 NEON, each
//! 256-bit vector held in two 128-bit registers: 16 lanes of 16 bits and 32 lanes of 8 bits, with
//! (`TB`) and without the traceback. aarch64 only (NEON is a baseline feature of the targets).
//! They run the recurrence of NCBI `ALIGN_EX` / `Blast_SemiGappedAlign` row by row, with the value
//! representation of the 128-bit kernels, and return exactly the `XdropResult` and edit operations
//! of those kernels.
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
//!
//! Why the result is exact. The argument at the head of `x_xdrop_avx2.rs` holds unchanged: the
//! same stored values (`score - base + bias`, 0 = minus infinity), the same `bias`, rebase limit,
//! acceptance tests and lookup tables; the same lanes per block (16 x i16, 32 x i8) and the same
//! band steps (16-bit: two 16-lane blocks per step while more than 16 lanes of the band are left,
//! then 16-lane blocks with the stop test on the grid `first_b_index + 16k`; 8-bit: 32-lane blocks,
//! then 16-lane blocks), so the same first and last survivors, cell count, best cell and traceback
//! rows; the traceback walk is the 128-bit kernels' `finish`. What differs from that file:
//!
//! 1. Registers. AVX2's low 128-bit half is `.0` of a `W`, its high half `.1`. Lane operations act
//!    on both registers; the operations that cross the halves (one lane up, prefix maximum, last
//!    lane to all lanes) are written for the pair, and the AVX2 shuffles that only undo the
//!    in-half behaviour of `unpack` / `packus` (`permute4x64::<0xD8>`) are not needed. Values that
//!    are the same in every lane (row-gap carry, running best, threshold, gap constants) are held
//!    in one register and used for both halves.
//! 2. Row gap (`score_gap_row`) in the form of item 2 of `x_xdrop_w128.rs` with blocks of `L` = 16
//!    lanes (16-bit) and 32 / 16 lanes (8-bit): `w[j] = v[j] + (j+1)*g` (no saturation: `v + L*g`
//!    fits the unsigned lane, `v <= max - gap_open_extend`, checked, otherwise None), `pw` its
//!    inclusive prefix maximum, `F[i] = sat(max(pw[i-1], carry) - i*g)` (`pw[-1]` = 0) and the next
//!    carry `sat(max(carry, pw[L-1]) - L*g)`. That item shows these are the values of the 128-bit
//!    kernel (and of the AVX2 file, whose form differs only in how the same maxima are grouped).
//! 3. Lane masks. The any-lane test is `umaxv` over the two registers; lane positions (the cell of
//!    a new best, the first and last survivor) come from a mask with 4 bits per lane (`shrn #4`),
//!    the 16-bit lanes first narrowed to one byte each (`uzp1`).
//! 4. Column residues and row tables: the 128-bit NEON `fill_cols` and `build_row` of
//!    `x_xdrop_w128.rs` (16 columns at a time; the same bytes as the AVX2 file's 32 at a time).
//!
//! Index bytes of 32 and up (the fence sentinel, letters the kernel declines, columns past the
//! band) give 0 through `tbl`; they reach only lanes the kernels mask or do not read, as in the
//! 128-bit NEON kernels.

#![allow(clippy::too_many_arguments)]

use super::arch::*;
use super::kernel::finish;
use super::x_w128::{build_row, fill_cols};
use super::{
    fits_8bit, ColSeq, RowSeq, Scores, XdropResult, XdropScratch, MAX_M, PAD, SCRIPT_EXTEND_GAP_A,
    SCRIPT_EXTEND_GAP_B, SCRIPT_GAP_IN_A, SCRIPT_SUB, TAB, VMAX,
};
use core::arch::aarch64::*;

// No NCBI counterpart: lane masks of the 256-bit lanes; they do not change any value NCBI
// computes.
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
/// Lane 0 of 16-bit lanes (the first lane of the band: no diagonal).
static HPM16: [i16; 8] = [-1, 0, 0, 0, 0, 0, 0, 0];
/// Lane 0 of 8-bit lanes.
static HPM8: [u8; 16] = [255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

/// No NCBI counterpart: LOSAT_X_DPNEON2 set (any value; not in the strict set, off unless set).
/// `neon_on` in `xdrop_simd.rs` lets it turn the NEON kernels on as LOSAT_X_DPNEON does.
pub(super) fn switch_on() -> bool {
    std::env::var_os("LOSAT_X_DPNEON2").is_some()
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:766-770
/// ```c
///     if (!score_only) {
///         return ALIGN_EX(A, B, M, N, a_offset, b_offset, edit_block, gap_align,
/// ```
/// Dispatch point of these kernels: `tb == false` is the score-only recurrence
/// (`Blast_SemiGappedAlign`, `s_BlastAlignPackedNucl`), `tb == true` is `ALIGN_EX`. The lane width
/// is the one the 128-bit dispatch would choose (`fits_8bit`). None: not handled here; the caller
/// runs the next dispatch.
///
/// # Safety
/// NEON (a baseline feature of the aarch64 targets).
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
// Lane helpers
// ---------------------------------------------------------------------------

/// No NCBI counterpart: one 256-bit vector of the AVX2 file in two registers (`.0` = lanes of the
/// low half, `.1` = lanes of the high half).
#[derive(Clone, Copy)]
struct W(V, V);

#[inline(always)]
unsafe fn ldw<T>(p: *const T) -> W {
    let p = p as *const u8;
    W(vld1q_u8(p), vld1q_u8(p.add(16)))
}
#[inline(always)]
unsafe fn stw<T>(p: *mut T, v: W) {
    let p = p as *mut u8;
    vst1q_u8(p, v.0);
    vst1q_u8(p.add(16), v.1);
}
/// A value that is the same in every lane, for both registers.
#[inline(always)]
fn bc(v: V) -> W {
    W(v, v)
}

#[inline(always)]
unsafe fn max_u16(a: V, b: V) -> V {
    vreinterpretq_u8_u16(vmaxq_u16(vreinterpretq_u16_u8(a), vreinterpretq_u16_u8(b)))
}
#[inline(always)]
unsafe fn adds_u16(a: V, b: V) -> V {
    vreinterpretq_u8_u16(vqaddq_u16(vreinterpretq_u16_u8(a), vreinterpretq_u16_u8(b)))
}
#[inline(always)]
unsafe fn max_u8(a: V, b: V) -> V {
    vmaxq_u8(a, b)
}
#[inline(always)]
unsafe fn adds_u8(a: V, b: V) -> V {
    vqaddq_u8(a, b)
}

// The lane operations of the `arch` layer (and the four above) on both registers.
macro_rules! pairwise {
    ($($w:ident = $f:ident;)*) => {
        $(
            #[inline(always)]
            unsafe fn $w(a: W, b: W) -> W {
                W($f(a.0, b.0), $f(a.1, b.1))
            }
        )*
    };
}
pairwise! {
    adds16w = adds16;
    add16w = add16;
    subsu16w = subsu16;
    max16w = max16;
    gt16w = gt16;
    max_u16w = max_u16;
    adds_u16w = adds_u16;
    adds8w = adds8;
    add8w = add8;
    subsu8w = subsu8;
    max8w = max8;
    gt8w = gt8;
    max_u8w = max_u8;
    adds_u8w = adds_u8;
    andw = and;
    orw = or;
    andnotw = andnot;
}

/// Bytes moved up by `N` lanes within one register, zeros shifted in.
#[inline(always)]
unsafe fn upb<const N: i32>(v: V) -> V {
    match N {
        1 => vextq_u8::<15>(vdupq_n_u8(0), v),
        2 => vextq_u8::<14>(vdupq_n_u8(0), v),
        4 => vextq_u8::<12>(vdupq_n_u8(0), v),
        _ => vextq_u8::<8>(vdupq_n_u8(0), v),
    }
}
/// 16-bit lanes moved up by one across the pair; lane 0 becomes 0 (`up1_16` of the AVX2 file).
#[inline(always)]
unsafe fn up1_16(v: W) -> W {
    W(vextq_u8::<14>(vdupq_n_u8(0), v.0), vextq_u8::<14>(v.0, v.1))
}
/// 8-bit lanes moved up by one across the pair; lane 0 becomes 0 (`up1_8` of the AVX2 file).
#[inline(always)]
unsafe fn up1_8(v: W) -> W {
    W(vextq_u8::<15>(vdupq_n_u8(0), v.0), vextq_u8::<15>(v.0, v.1))
}
/// 16-bit lane 7 of one register to every lane.
#[inline(always)]
unsafe fn dup7_16(v: V) -> V {
    vreinterpretq_u8_u16(vdupq_laneq_u16::<7>(vreinterpretq_u16_u8(v)))
}
/// 16-bit lane 15 of the pair to every lane of one register (`last16` of the AVX2 file).
#[inline(always)]
unsafe fn last16(v: W) -> V {
    dup7_16(v.1)
}
/// 8-bit lane 31 of the pair to every lane of one register (`last8` of the AVX2 file).
#[inline(always)]
unsafe fn last8(v: W) -> V {
    vdupq_laneq_u8::<15>(v.1)
}
/// Inclusive prefix maximum over the 16 lanes of 16 bits, unsigned (`pmax16` of the AVX2 file).
#[inline(always)]
unsafe fn pmax_u16w(w: W) -> W {
    let a = max_u16(w.0, upb::<2>(w.0));
    let a = max_u16(a, upb::<4>(a));
    let a = max_u16(a, upb::<8>(a));
    let b = max_u16(w.1, upb::<2>(w.1));
    let b = max_u16(b, upb::<4>(b));
    let b = max_u16(b, upb::<8>(b));
    W(a, max_u16(b, dup7_16(a)))
}
/// The same, signed (for `H`, whose lanes are 0..=32767, where both orders agree).
#[inline(always)]
unsafe fn pmax_s16w(w: W) -> W {
    let a = max16(w.0, upb::<2>(w.0));
    let a = max16(a, upb::<4>(a));
    let a = max16(a, upb::<8>(a));
    let b = max16(w.1, upb::<2>(w.1));
    let b = max16(b, upb::<4>(b));
    let b = max16(b, upb::<8>(b));
    W(a, max16(b, dup7_16(a)))
}
/// Inclusive prefix maximum over the 16 lanes of one register, 8 bits, unsigned.
#[inline(always)]
unsafe fn pmax_u8x(w: V) -> V {
    let p = max_u8(w, upb::<1>(w));
    let p = max_u8(p, upb::<2>(p));
    let p = max_u8(p, upb::<4>(p));
    max_u8(p, upb::<8>(p))
}
/// The same, signed (for `H`, whose lanes are 0..=127).
#[inline(always)]
unsafe fn pmax_s8x(w: V) -> V {
    let p = max8(w, upb::<1>(w));
    let p = max8(p, upb::<2>(p));
    let p = max8(p, upb::<4>(p));
    max8(p, upb::<8>(p))
}
/// Inclusive prefix maximum over the 32 lanes of 8 bits, unsigned (`pmax8` of the AVX2 file).
#[inline(always)]
unsafe fn pmax_u8w(w: W) -> W {
    let a = pmax_u8x(w.0);
    let b = pmax_u8x(w.1);
    W(a, max_u8(b, vdupq_laneq_u8::<15>(a)))
}
/// The same, signed.
#[inline(always)]
unsafe fn pmax_s8w(w: W) -> W {
    let a = pmax_s8x(w.0);
    let b = pmax_s8x(w.1);
    W(a, max8(b, vdupq_laneq_u8::<15>(a)))
}
/// Any bit set in one register.
#[inline(always)]
unsafe fn any(v: V) -> bool {
    vmaxvq_u32(vreinterpretq_u32_u8(v)) != 0
}
/// Any bit set in the pair (`movemask != 0` of the AVX2 file).
#[inline(always)]
unsafe fn anyw(v: W) -> bool {
    any(vorrq_u8(v.0, v.1))
}
/// A byte mask (each byte 0 or 255) as 4 bits per byte: byte `i` at bits `4i..4i + 4`.
#[inline(always)]
unsafe fn nib8(m: V) -> u64 {
    vget_lane_u64::<0>(vreinterpret_u64_u8(vshrn_n_u16::<4>(vreinterpretq_u16_u8(
        m,
    ))))
}
/// A 16-bit lane mask of the pair (each lane 0 or 0xFFFF) as 4 bits per lane, lanes 0..16.
#[inline(always)]
unsafe fn nib16w(m: W) -> u64 {
    nib8(vuzp1q_u8(m.0, m.1))
}
/// An 8-bit lane mask of the pair as 4 bits per lane, lanes 0..32.
#[inline(always)]
unsafe fn nib8w(m: W) -> u128 {
    nib8(m.0) as u128 | (nib8(m.1) as u128) << 64
}
/// The script bytes of 16 lanes of 16 bits: the low byte of each lane (the values are 0..=0x56,
/// so this is the unsigned saturation `packus_epi16` of the AVX2 file).
#[inline(always)]
unsafe fn pack16(r: W) -> V {
    vuzp1q_u8(r.0, r.1)
}
/// Scores of 16 lanes of 16 bits from the low and high byte planes of a row table.
#[inline(always)]
unsafe fn scores16(tl: uint8x16x2_t, th: uint8x16x2_t, idx: V) -> W {
    let lo = vqtbl2q_u8(tl, idx);
    let hi = vqtbl2q_u8(th, idx);
    W(vzip1q_u8(lo, hi), vzip2q_u8(lo, hi))
}

// ---------------------------------------------------------------------------
// 16-bit lanes, 16 per block
// ---------------------------------------------------------------------------

// No NCBI counterpart: the gap costs and constants of one call on 16-bit lanes.
struct C16 {
    goe: V,
    ge1: V,
    ge16: V,
    /// `i * gap_extend` in lane `i`.
    ramp: W,
    /// `(i + 1) * gap_extend` in lane `i`.
    ramp1: W,
    xv: V,
    c3: V,
    c10: V,
    c40: V,
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
// `best_score - x_dropoff` (each the same in every lane of one register) and the cell where the
// best was last raised, carried across the blocks of a row.
struct Row {
    carry: V,
    bestv: V,
    thrv: V,
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
/// `F` = `score_gap_row`, `H` = `score`, the lane operations of `block16` of the AVX2 file with the
/// row gap of item 2 at the head of this file; `hpm` clears the first lane's diagonal (lanes of
/// `.0`). Returns (drop mask, script).
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
    hpm: V,
) -> (W, W) {
    let hl = ldw(hr.add(b).sub(1));
    let hp = W(andnot(hpm, hl.0), hl.1);
    let ep = ldw(e.add(b));
    let d = adds16w(hp, s);
    let h0 = max16w(d, ep);
    let v = subsu16w(h0, bc(c.goe));
    // row gap: pw = inclusive prefix maximum of v[j] + (j+1)g (unsigned, exact)
    let pw = pmax_u16w(adds_u16w(v, c.ramp1));
    let f = subsu16w(max_u16w(up1_16(pw), bc(rs.carry)), c.ramp);
    rs.carry = subsu16(max_u16(rs.carry, last16(pw)), c.ge16);
    let mut hh = max16w(h0, f);
    if let Some(m) = lane_mask {
        hh = andw(hh, m);
    }
    let hg = subsu16w(hh, bc(c.goe));
    let es = subsu16w(ep, bc(c.ge1));
    let en = max16w(hg, es);
    let gt = gt16w(hh, bc(rs.bestv));
    let thr;
    if anyw(gt) {
        // A new best inside this block: exact running best per lane (the exclusive prefix
        // maximum of `H`, the inclusive one moved up by a lane).
        let pfx = max16w(up1_16(pmax_s16w(hh)), bc(rs.bestv));
        thr = subsu16w(pfx, bc(c.xv));
        let rm = nib16w(gt16w(hh, pfx));
        rs.rec_b = b + ((63 - rm.leading_zeros()) >> 2) as usize;
        rs.bestv = last16(max16w(pfx, hh));
        rs.thrv = subsu16(rs.bestv, c.xv);
        rs.slow = true;
    } else {
        thr = bc(rs.thrv);
    }
    let drop = gt16w(thr, hh);
    stw(hw.add(b), andnotw(drop, hh));
    stw(e.add(b), en);
    let scr = if TB {
        let fs = subsu16w(f, bc(c.ge1));
        let me = gt16w(ep, d);
        let mf = gt16w(f, h0);
        let op = andnotw(mf, add16w(bc(c.c3), andw(me, bc(c.c3))));
        let fb = andnotw(gt16w(hg, es), bc(c.c40));
        let fa = andnotw(gt16w(hg, fs), bc(c.c10));
        orw(op, orw(fb, fa))
    } else {
        bc(zero())
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
// The whole DP on 16-bit lanes, the steps of `align16` of the AVX2 file in the same order (row 0,
// then per row: column residues, fence, row table, rebase, blocks, band update). None for a
// problem the lanes do not hold (the caller then runs the next dispatch).
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
    // No NCBI counterpart: the prefix-maximum form of the row gap needs v + 16 * gap_extend to
    // fit the unsigned 16-bit lane, with v <= 32767 - gap_open_extend (item 2).
    if VMAX - goe as i64 + 16 * gap_extend as i64 > 65535 {
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
    let mut ramp1 = [0i16; 16];
    for (i, (r, r1)) in ramp.iter_mut().zip(ramp1.iter_mut()).enumerate() {
        *r = (i as i32 * gap_extend) as i16;
        *r1 = ((i as i32 + 1) * gap_extend) as i16;
    }
    let c = C16 {
        goe: splat16(goe as i16),
        ge1: splat16(gap_extend as i16),
        ge16: splat16((16 * gap_extend) as i16),
        ramp: ldw(ramp.as_ptr()),
        ramp1: ldw(ramp1.as_ptr()),
        xv: splat16(x as i16),
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
        // low bytes of entries 0..32, then high bytes (`build_row_table`)
        let tl = uint8x16x2_t(load(t), load(t.add(16)));
        let th = uint8x16x2_t(load(t.add(32)), load(t.add(48)));

        if best_rel > limit_rebase {
            let dv = bc(splat16(best_rel as i16));
            let mut b = first;
            while b < b_size {
                stw(hr.add(b), subsu16w(ldw(hr.add(b)), dv));
                stw(e.add(b), subsu16w(ldw(e.add(b)), dv));
                b += 16;
            }
            base += best_rel;
            best_rel = 0;
        }
        // no diagonal into the first lane of the band: lane 0 of the row's first block (`hpm`)
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
        // Two blocks (32 lanes) while more than 16 lanes of the band are left, then single blocks
        // for the rest of the band and the row-gap extension, as in the AVX2 file.
        while b0 + 16 < bend {
            // scores of lanes b0..b0+32 (lane b reads column b-1)
            let s0 = scores16(tl, th, load(sb.add(b0)));
            let s1 = scores16(tl, th, load(sb.add(b0 + 16)));
            // the first block lies inside the band (b0 + 16 < bend <= limit); the second may
            // reach past the band and past `limit`
            let mask1 = if b0 + 32 <= limit {
                None
            } else {
                Some(ldw(LANE_MASK16.as_ptr().add(32 - (limit - b0))))
            };
            let (d0, r0) = block16::<TB>(&c, &mut rs, hr, hw, e, b0, s0, None, hpm);
            hpm = zero_v;
            let (d1, r1) = block16::<TB>(&c, &mut rs, hr, hw, e, b0 + 16, s1, mask1, zero_v);
            if TB {
                store(tr.add(b0 - first), pack16(r0));
                store(tr.add(b0 - first + 16), pack16(r1));
            }
            // 4 bits per lane
            let live = !(nib16w(d0) as u128 | (nib16w(d1) as u128) << 64);
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + (live.trailing_zeros() as usize >> 2);
                }
                hi_surv = b0 + ((127 - live.leading_zeros() as usize) >> 2);
            }
            b0 += 32;
        }
        loop {
            if b0 >= bend {
                if b0 >= limit {
                    break;
                }
                if lane0(rs.carry) < lane0(rs.thrv) {
                    break;
                }
            }
            let s0 = scores16(tl, th, load(sb.add(b0)));
            let mask = if b0 + 16 <= limit {
                None
            } else {
                Some(ldw(LANE_MASK16.as_ptr().add(16 - (limit - b0))))
            };
            let (d0, r0) = block16::<TB>(&c, &mut rs, hr, hw, e, b0, s0, mask, hpm);
            hpm = zero_v;
            if TB {
                store(tr.add(b0 - first), pack16(r0));
            }
            let live = !nib16w(d0);
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + (live.trailing_zeros() as usize >> 2);
                }
                hi_surv = b0 + ((63 - live.leading_zeros() as usize) >> 2);
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
            stw(hw.add(b), bc(zero_v));
            b += 16;
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
            stw(e.add(b), bc(zero_v));
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
// 8-bit lanes, 32 per block (16 per block in the row tails)
// ---------------------------------------------------------------------------

// No NCBI counterpart: the gap costs and constants of one call on 8-bit lanes.
struct C8 {
    goe: V,
    ge1: V,
    ge16: V,
    ge32: V,
    /// `i * gap_extend` in lane `i`.
    ramp: W,
    /// `(i + 1) * gap_extend` in lane `i`.
    ramp1: W,
    xv: V,
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
/// One block of 32 cells on 8-bit lanes: the lane operations of `block8` of the AVX2 file, row
/// gap as in item 2 at the head of this file. Returns (drop mask, script).
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
    hpm: V,
) -> (W, W) {
    let hl = ldw(hr.add(b).sub(1));
    let hp = W(andnot(hpm, hl.0), hl.1);
    let ep = ldw(e.add(b));
    let d = adds8w(hp, s);
    let h0 = max8w(d, ep);
    let v = subsu8w(h0, bc(c.goe));
    let pw = pmax_u8w(adds_u8w(v, c.ramp1));
    let f = subsu8w(max_u8w(up1_8(pw), bc(rs.carry)), c.ramp);
    rs.carry = subsu8(max_u8(rs.carry, last8(pw)), c.ge32);
    let mut hh = max8w(h0, f);
    if let Some(m) = lane_mask {
        hh = andw(hh, m);
    }
    let hg = subsu8w(hh, bc(c.goe));
    let es = subsu8w(ep, bc(c.ge1));
    let en = max8w(hg, es);
    let gt = gt8w(hh, bc(rs.bestv));
    let thr;
    if anyw(gt) {
        let pfx = max8w(up1_8(pmax_s8w(hh)), bc(rs.bestv));
        thr = subsu8w(pfx, bc(c.xv));
        let rm = nib8w(gt8w(hh, pfx));
        rs.rec_b = b + ((127 - rm.leading_zeros()) >> 2) as usize;
        rs.bestv = last8(max8w(pfx, hh));
        rs.thrv = subsu8(rs.bestv, c.xv);
        rs.slow = true;
    } else {
        thr = bc(rs.thrv);
    }
    let drop = gt8w(thr, hh);
    stw(hw.add(b), andnotw(drop, hh));
    stw(e.add(b), en);
    let scr = if TB {
        let fs = subsu8w(f, bc(c.ge1));
        let me = gt8w(ep, d);
        let mf = gt8w(f, h0);
        let op = andnotw(mf, add8w(bc(c.c3), andw(me, bc(c.c3))));
        let fb = andnotw(gt8w(hg, es), bc(c.c40));
        let fa = andnotw(gt8w(hg, fs), bc(c.c10));
        orw(op, orw(fb, fa))
    } else {
        bc(zero())
    };
    (drop, scr)
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:578-601,616-631
/// ```c
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
/// ...
///             if (best_score - score > x_dropoff) {
/// ...
///                 score_array[b_index].best = score;
/// ```
/// `block8` on 16 lanes (one register: the tail of a row and the row-gap extension), with the
/// same values (`block8x` of the AVX2 file).
#[inline(always)]
unsafe fn block8x<const TB: bool>(
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
    let pw = pmax_u8x(adds_u8(v, c.ramp1.0));
    let f = subsu8(max_u8(upb::<1>(pw), rs.carry), c.ramp.0);
    rs.carry = subsu8(max_u8(rs.carry, vdupq_laneq_u8::<15>(pw)), c.ge16);
    let mut hh = max8(h0, f);
    if let Some(m) = lane_mask {
        hh = and(hh, m);
    }
    let hg = subsu8(hh, c.goe);
    let es = subsu8(ep, c.ge1);
    let en = max8(hg, es);
    let gt = gt8(hh, rs.bestv);
    let thr;
    if any(gt) {
        let pfx = max8(upb::<1>(pmax_s8x(hh)), rs.bestv);
        thr = subsu8(pfx, c.xv);
        let rm = nib8(gt8(hh, pfx));
        rs.rec_b = b + ((63 - rm.leading_zeros()) >> 2) as usize;
        rs.bestv = vdupq_laneq_u8::<15>(max8(pfx, hh));
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
// The whole DP on 8-bit lanes, the steps of `align8` of the AVX2 file in the same order. `max_m`
// is the largest score of the matrix (`fits_8bit`). None for a problem the lanes do not hold.
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
    // No NCBI counterpart: the prefix-maximum form of the row gap needs v + 32 * gap_extend to
    // fit the unsigned 8-bit lane, with v <= 127 - gap_open_extend (item 2).
    if 127 - goe as i64 + 32 * gap_extend as i64 > 255 {
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
    let mut ramp1 = [0u8; 32];
    for (i, (r, r1)) in ramp.iter_mut().zip(ramp1.iter_mut()).enumerate() {
        *r = (i as i32 * gap_extend) as u8;
        *r1 = ((i as i32 + 1) * gap_extend) as u8;
    }
    let c = C8 {
        goe: splat8(sat(goe)),
        ge1: splat8(sat(gap_extend)),
        ge16: splat8(sat(16 * gap_extend)),
        ge32: splat8(sat(32 * gap_extend)),
        ramp: ldw(ramp.as_ptr()),
        ramp1: ldw(ramp1.as_ptr()),
        xv: splat8(sat(x)),
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
        let lookup = |idx: V| {
            if small_alphabet {
                pick(t0, idx)
            } else {
                lut32(t0, t1, idx)
            }
        };

        if best_rel > limit_rebase {
            let dv = bc(splat8(best_rel as u8 as i8));
            let mut b = first;
            while b < b_size {
                stw(hr.add(b), subsu8w(ldw(hr.add(b)), dv));
                stw(e.add(b), subsu8w(ldw(e.add(b)), dv));
                b += 32;
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
        // 32-lane blocks while more than 16 lanes of the band are left, then 16-lane blocks for the
        // rest of the band and the row-gap extension (stop test as in `align16`).
        while b0 + 16 < bend {
            let sv = W(lookup(load(sb.add(b0))), lookup(load(sb.add(b0 + 16))));
            let mask = if b0 + 32 <= limit {
                None
            } else {
                Some(ldw(LANE_MASK8.as_ptr().add(32 - (limit - b0))))
            };
            let (d0, r0) = block8::<TB>(&c, &mut rs, hr, hw, e, b0, sv, mask, hpm);
            hpm = zero_v;
            if TB {
                stw(tr.add(b0 - first), r0);
            }
            let live = !nib8w(d0);
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + (live.trailing_zeros() as usize >> 2);
                }
                hi_surv = b0 + ((127 - live.leading_zeros() as usize) >> 2);
            }
            b0 += 32;
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
            let sv = lookup(load(sb.add(b0)));
            let mask = if b0 + 16 <= limit {
                None
            } else {
                Some(load(LANE_MASK8X.as_ptr().add(16 - (limit - b0))))
            };
            let (d0, r0) = block8x::<TB>(&c, &mut rs, hr, hw, e, b0, sv, mask, hpm);
            hpm = zero_v;
            if TB {
                store(tr.add(b0 - first), r0);
            }
            let live = !nib8(d0);
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + (live.trailing_zeros() as usize >> 2);
                }
                hi_surv = b0 + ((63 - live.leading_zeros() as usize) >> 2);
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
            stw(hw.add(b), bc(zero_v));
            b += 32;
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
            stw(e.add(b), bc(zero_v));
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
    //! The two-register kernels against the scalar transliteration of the NCBI loops and against
    //! the 128-bit NEON kernels of `xdrop_simd.rs` (random problems: forward and reversed rows and
    //! columns, both lane widths, traceback and score-only, fences, scale 1 and 32, bands wider than
    //! 256 cells, gap costs at both row-gap preconditions), and the pair helpers against scalar
    //! models of their lanes.
    use super::super::kernel;
    use super::super::tests::{random_problem, reference, Problem, Rng};
    use super::super::FENCE_SENTRY;
    use super::*;

    unsafe fn v(b: [u8; 16]) -> V {
        vld1q_u8(b.as_ptr())
    }
    unsafe fn b(x: V) -> [u8; 16] {
        let mut o = [0u8; 16];
        vst1q_u8(o.as_mut_ptr(), x);
        o
    }
    unsafe fn wv(x: [u8; 32]) -> W {
        ldw(x.as_ptr())
    }
    unsafe fn wb(x: W) -> [u8; 32] {
        let mut o = [0u8; 32];
        stw(o.as_mut_ptr(), x);
        o
    }
    fn u16s(x: [u8; 32]) -> [u16; 16] {
        core::array::from_fn(|i| u16::from_le_bytes([x[2 * i], x[2 * i + 1]]))
    }
    fn bytes(rng: &mut Rng) -> [u8; 32] {
        // mostly small values, some near the ends of the lanes
        core::array::from_fn(|_| match rng.below(4) {
            0 => 255 - rng.below(3) as u8,
            1 => rng.below(3) as u8,
            _ => rng.below(256) as u8,
        })
    }

    /// Each pair helper against the scalar definition of its lanes, on random inputs.
    #[test]
    fn neon2_pair_helpers_match_scalar_models() {
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        for _ in 0..20_000 {
            let x = bytes(&mut rng);
            let x16 = u16s(x);
            // SAFETY: NEON is a baseline feature of the aarch64 targets.
            unsafe {
                // one lane up
                let up = u16s(wb(up1_16(wv(x))));
                for i in 0..16 {
                    assert_eq!(up[i], if i == 0 { 0 } else { x16[i - 1] });
                }
                let up = wb(up1_8(wv(x)));
                for i in 0..32 {
                    assert_eq!(up[i], if i == 0 { 0 } else { x[i - 1] });
                }
                // prefix maxima
                // the signed forms are applied to `H`, which is never negative: inputs 0..=i16::MAX
                // (zeros are shifted in, so negative lanes would give a mix of 0 and the inputs)
                let xp: [u8; 32] =
                    core::array::from_fn(|i| if i % 2 == 1 { x[i] & 0x7F } else { x[i] });
                let xp16 = u16s(xp);
                let pu = u16s(wb(pmax_u16w(wv(x))));
                let ps = u16s(wb(pmax_s16w(wv(xp))));
                let (mut mu, mut ms) = (0u16, 0u16);
                for i in 0..16 {
                    mu = mu.max(x16[i]);
                    ms = ms.max(xp16[i]);
                    assert_eq!(pu[i], mu);
                    assert_eq!(ps[i], ms);
                }
                let x7: [u8; 32] = core::array::from_fn(|i| x[i] & 0x7F);
                let pu = wb(pmax_u8w(wv(x)));
                let ps = wb(pmax_s8w(wv(x7)));
                let (mut mu, mut ms) = (0u8, 0u8);
                for i in 0..32 {
                    mu = mu.max(x[i]);
                    ms = ms.max(x7[i]);
                    assert_eq!(pu[i], mu);
                    assert_eq!(ps[i], ms);
                }
                // last lane to all lanes
                let l = b(last16(wv(x)));
                for i in 0..8 {
                    assert_eq!(u16::from_le_bytes([l[2 * i], l[2 * i + 1]]), x16[15]);
                }
                assert_eq!(b(last8(wv(x))), [x[31]; 16]);
                // any bit
                let sparse: [u8; 32] =
                    core::array::from_fn(|i| if rng.below(40) == 0 { x[i] } else { 0 });
                assert_eq!(anyw(wv(sparse)), sparse.iter().any(|&c| c != 0));
                // lane masks to 4 bits per lane
                let m8: [u8; 32] = core::array::from_fn(|i| if x[i] & 1 != 0 { 255 } else { 0 });
                let want8: u128 = (0..32).fold(0, |a, i| {
                    a | if m8[i] != 0 { 0xFu128 << (4 * i) } else { 0 }
                });
                assert_eq!(nib8w(wv(m8)), want8);
                let lo: [u8; 16] = core::array::from_fn(|i| m8[i]);
                assert_eq!(nib8(v(lo)), want8 as u64);
                let m16: [u8; 32] = core::array::from_fn(|i| m8[i & !1]);
                let want16: u64 = (0..16).fold(0, |a, i| {
                    a | if m16[2 * i] != 0 {
                        0xFu64 << (4 * i)
                    } else {
                        0
                    }
                });
                assert_eq!(nib16w(wv(m16)), want16);
                // script bytes
                let r: [u8; 32] =
                    core::array::from_fn(|i| if i % 2 == 0 { (x[i] % 0x57) as u8 } else { 0 });
                let p = b(pack16(wv(r)));
                for i in 0..16 {
                    assert_eq!(p[i], r[2 * i]);
                }
                // 16-bit scores from the two byte planes (an index of 32 and up gives 0)
                let tl = bytes(&mut rng);
                let th = bytes(&mut rng);
                let idx: [u8; 16] =
                    core::array::from_fn(|i| if rng.below(8) == 0 { x[i] } else { x[i] % 32 });
                let s = u16s(wb(scores16(
                    uint8x16x2_t(ldw(tl.as_ptr()).0, ldw(tl.as_ptr()).1),
                    uint8x16x2_t(ldw(th.as_ptr()).0, ldw(th.as_ptr()).1),
                    v(idx),
                )));
                for i in 0..16 {
                    let k = idx[i] as usize;
                    let want = if k < 32 {
                        u16::from_le_bytes([tl[k], th[k]])
                    } else {
                        0
                    };
                    assert_eq!(s[i], want, "scores16 lane {i} index {k}");
                }
            }
        }
    }

    /// The 128-bit NEON kernel of `xdrop_simd.rs` or the kernels of this file; `m8` selects the
    /// 8-bit lanes (as `fits_8bit` would).
    unsafe fn run(
        two: bool,
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
        match (two, tb, m8) {
            (false, true, Some(m)) => {
                kernel::align_body8::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (false, false, Some(m)) => {
                kernel::align_body8::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (false, true, None) => {
                kernel::align_body::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
            (false, false, None) => {
                kernel::align_body::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
            (true, true, Some(m)) => align8::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc),
            (true, false, Some(m)) => {
                align8::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (true, true, None) => align16::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc),
            (true, false, None) => align16::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc),
        }
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
    // Random problems run through the scalar loops (`reference`), the 128-bit NEON kernels and the
    // kernels of this file, with the rows and columns given forward and reversed (as the call sites
    // pass them for the left extension), on 16-bit lanes and, where the problem fits, 8-bit lanes,
    // with and without traceback. Offsets, score, fence flag, edit script and cell count must be
    // equal; a call these kernels decline must be one that fails a row-gap precondition.
    // `LOSAT_FUZZ_CASES` sets the number of problems.
    #[test]
    fn neon2_kernels_match_scalar_loops_on_random_problems() {
        let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(600);
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let (mut sc128, mut sc2) = (XdropScratch::new(), XdropScratch::new());
        // [width 8 / 16][tb][reverse]
        let mut handled = [[[0usize; 2]; 2]; 2];
        let (mut fences, mut ops_total, mut wide) = (0usize, 0usize, 0usize);
        // calls at the row-gap preconditions: 8-bit with 32 * gap_extend - gap_open_extend in
        // 113..=128 (handled), above 128 or 16-bit above 32768 (declined, run by the 128-bit
        // kernel); 16-bit with gap_extend >= 1000 handled
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
                // 8-bit lanes at the row-gap precondition (32 * ge - goe around 128)
                8 => {
                    let ge = 4 + rng.below(2) as i32;
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
                        // SAFETY: NEON is a baseline feature of the aarch64 targets.
                        let (r128, r2) = unsafe {
                            (
                                run(
                                    false, m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb,
                                    &mut sc128,
                                ),
                                run(
                                    true, m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb,
                                    &mut sc2,
                                ),
                            )
                        };
                        let edge = if wide8 {
                            32 * ge - (go + ge)
                        } else {
                            16 * ge - (go + ge)
                        };
                        match (r128, r2) {
                            (Some(a), Some(b)) => {
                                assert_eq!(a, b, "{tag}: 128 vs neon2");
                                assert_eq!(sc128.ops, sc2.ops, "{tag}: 128 vs neon2 ops");
                            }
                            (None, None) => {}
                            (Some(_), None) => {
                                assert!(
                                    edge > if wide8 { 128 } else { 32768 },
                                    "{tag}: neon2 declined"
                                );
                                edge_out += 1;
                            }
                            (None, Some(_)) => {
                                panic!("{tag}: neon2 handled a call the 128-bit kernel declines")
                            }
                        }
                        let Some(r) = r2 else {
                            continue;
                        };
                        assert_eq!(
                            (r.a_offset, r.b_offset, r.score, r.fence_hit),
                            (expected.0, expected.1, expected.2, expected.3),
                            "{tag}: neon2 vs scalar"
                        );
                        assert_eq!(sc2.ops, expected.4, "{tag}: neon2 vs scalar ops");
                        if !r.fence_hit {
                            assert_eq!(r.cells, expected.5, "{tag}: cells");
                        }
                        handled[usize::from(!wide8)][usize::from(tb)][usize::from(reverse)] += 1;
                        fences += usize::from(r.fence_hit);
                        ops_total += sc2.ops.len();
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
            "neon2 fuzz: cases={cases} handled[w8/w16][so/tb][fwd/rev]={handled:?} fences={fences} \
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
}
