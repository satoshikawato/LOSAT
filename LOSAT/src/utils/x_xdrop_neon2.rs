//! LOSAT_X_DPNEON2 (a switch of the strict set: on unless `=0`, which leaves every call to the
//! DPW128 / 128-bit NEON kernels): the 256-bit kernels of `x_xdrop_avx2.rs` on aarch64 NEON, each
//! 256-bit vector held in two 128-bit registers: 16 lanes of 16 bits and 32 lanes of 8 bits, with
//! (`TB`) and without the traceback. aarch64 only (NEON is a baseline feature of the targets).
//! They run the recurrence of NCBI `ALIGN_EX` / `Blast_SemiGappedAlign` row by row, with the value
//! representation of the 128-bit kernels, and return exactly the `XdropResult` and edit operations
//! of those kernels.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:563-635
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:652-673
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
//! 5. Rows over the band (No NCBI counterpart): as in the AVX2 file, a kernel whose rows do not
//!    hold the band's end calls `XdropScratch::x_room`, which moves its rows down (by a
//!    multiple of 64 lanes, so every lane keeps its place within the 64-byte-aligned rows; the
//!    blocks start at `first` in every row either way, and at least 64 lanes stay below `first`)
//!    and lowers `first`, `b_size`, `len2` and the buffers' write ends by
//!    the same amount; the best cell and the traceback rows add `x_sh` back, so they keep absolute
//!    columns, and the result is the same as with rows that hold every column.
//!
//! Index bytes of 32 and up (the fence sentinel, letters the kernel declines, columns past the
//! band) give 0 through `tbl`; they reach only lanes the kernels mask or do not read, as in the
//! 128-bit NEON kernels.

#![allow(clippy::too_many_arguments)]

use super::arch::*;
use super::kernel::finish;
use super::x_w128::{build_row, fill_cols};
use super::{
    fits_8bit, script_codes, ColSeq, RowSeq, Scores, XdropResult, XdropScratch, MAX_M, PAD, TAB,
    VMAX,
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

/// No NCBI counterpart: LOSAT_X_DPNEON2, a switch of the strict set: on unless `LOSAT_X_DPNEON2=0`,
/// which leaves the calls to the DPW128 / 128-bit NEON kernels. The two-register kernels run first
/// on aarch64 (`x_extra_switches` in `xdrop_simd.rs`, which reads this once).
pub(super) fn switch_on() -> bool {
    crate::utils::xswitch::on("LOSAT_X_DPNEON2")
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:766-770
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
        (true, Some(m)) if sc.x_nib => align8::<true, true>(
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
        (true, Some(m)) => align8::<true, false>(
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
        (false, Some(m)) => align8::<false, false>(
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
        (true, None) if sc.x_nib => align16::<true, true>(
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
        (true, None) => align16::<true, false>(
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
        (false, None) => align16::<false, false>(
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
/// Any bit set in the pair (the nonzero test of `movemask` in the AVX2 file).
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
/// The script bytes of 16 lanes of 16 bits: the low byte of each lane (the values are 0..=0x56, or
/// four-bit codes, so this is the unsigned saturation `packus_epi16` of the AVX2 file).
#[inline(always)]
unsafe fn pack16(r: W) -> V {
    vuzp1q_u8(r.0, r.1)
}
/// 32 four-bit codes (bytes <= 15) to 16 bytes, byte k = code 2k | code 2k+1 << 4
/// (`LOSAT_X_DPNIBBLE`; `nibs` of the `arch` layer for 16).
#[inline(always)]
unsafe fn nibw(v: W) -> V {
    let a = vreinterpretq_u16_u8(v.0);
    let b = vreinterpretq_u16_u8(v.1);
    vuzp1q_u8(
        vreinterpretq_u8_u16(vsraq_n_u16::<4>(a, a)),
        vreinterpretq_u8_u16(vsraq_n_u16::<4>(b, b)),
    )
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
// The whole DP on 16-bit lanes, the steps of `align16` of the AVX2 file in the same order (row 0,
// then per row: column residues, fence, row table, rebase, blocks, band update). None for a
// problem the lanes do not hold (the caller then runs the next dispatch).
unsafe fn align16<const TB: bool, const NIB: bool>(
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
        sc.trace.x_row_zero(NIB, b_size);
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
            if TB && NIB {
                store(tr.add((b0 - first) >> 1), nibw(W(pack16(r0), pack16(r1))));
            } else if TB {
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
            if TB && NIB {
                store64(tr.add((b0 - first) >> 1), nibs(pack16(r0)));
            } else if TB {
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
    finish::<TB, NIB>(sc, a_off, b_off, best, fence_hit, cells)
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
// The whole DP on 8-bit lanes, the steps of `align8` of the AVX2 file in the same order. `max_m`
// is the largest score of the matrix (`fits_8bit`). None for a problem the lanes do not hold.
unsafe fn align8<const TB: bool, const NIB: bool>(
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
        sc.trace.x_row_zero(NIB, b_size);
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
        c3: splat8(const { script_codes(NIB).0 } as i8),
        c10: splat8(const { script_codes(NIB).1 } as i8),
        c40: splat8(const { script_codes(NIB).2 } as i8),
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
            if TB && NIB {
                store(tr.add((b0 - first) >> 1), nibw(r0));
            } else if TB {
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
            if TB && NIB {
                store64(tr.add((b0 - first) >> 1), nibs(r0));
            } else if TB {
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
    finish::<TB, NIB>(sc, a_off, b_off, best, fence_hit, cells)
}

#[cfg(test)]
#[path = "x_xdrop_neon2_tests.rs"]
mod tests;
