//! The two-hit loop of tblastn's unit path (`XTwoHit::pairs`)
//! with its two common outcomes, too far and overlap, taken without a branch.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-551
//! On a genome searched with proteins most hits end in one of the two `continue`s of a flag-0 cell:
//! too far (83-87 % of the hits of the four ledger cases, v4-A10 probe 1) or overlap (12-15 %).
//! Which one is hard to predict: taking both without a branch removed 0.12-0.17 mispredicted
//! branches per hit from the whole run (0.25-0.35 before). A flag-1 cell (0.5-0.7 %) or a hit past
//! both tests, a candidate for the extension (1.2-1.6 %), is rare.
//!
//! `x_hits_bf` reads the cell of each hit and computes `d = (s + diag_offset) - last_hit`. That is
//! the reference's `diff = s - (last_hit - diag_offset)` in wrapping `i32` operations, so in a release
//! build (no overflow checks) the same value whatever the operands; and the operands stay far from
//! overflow (`Diagonals::finish_frame` starts again from `window` once the offset reaches
//! INT4_MAX / 4). A candidate (`wordsize <= d < window`) is
//! `(d - wordsize) as u32 < window - wordsize`, one unsigned compare: for `d >= window` the
//! difference is at least `window - wordsize`; for `d < wordsize` it is negative, or, when it wraps
//! (`d < i32::MIN + wordsize`), at least `2^31 - wordsize`, which is not below
//! `window - wordsize` because `window <= i32::MAX`. So the test is exact whenever
//! `0 <= wordsize < window` (`fits`). A flag-1 cell or a candidate stops the loop at that hit; the
//! caller (`XTwoHit::x_pairs_bf`) runs the body of `pairs` on that one hit and calls the loop again
//! after it. Every other hit is a flag-0 cell whose next `last_hit` is selected without a branch:
//! `s + diag_offset` when too far (`d >= window`), the value it read when it overlaps. The cell's
//! flag stays 0. So every cell holds after each hit what the reference leaves there, and the
//! reference's body sees the same table for every hit it runs and saves the same init HSPs in the
//! same order.
//!
//! The cell index is the hit's diagonal masked by the table's length minus one. The table's length
//! is a power of two and its `mask` is that length minus one (`Diagonals::new`), so this is the
//! reference's index, and it is never past the end of the table.
//!
//! The loop runs for every offset array of the unit path (the former switch `LOSAT_X_TBNHITBF`
//! graduated on 2026-10-26). TBLASTX's unit scan has a loop of this form under its own switch,
//! LOSAT_X_TBXHITBF.
//!
//! EXPERIMENT (LOSAT_X_TBNHITSTEP, in the strict set): `x_hits_steps`, the same outcomes from one
//! compare, two hits a step. With `dw = diff - wordsize` and `k = window - wordsize`, a candidate is
//! `0 <= dw < k` (`dw` below `k` as an unsigned number) and too far is `diff >= window`, which is
//! `dw >= k` whenever `diff - wordsize` is exact; so one compare of `dw` with `k` gives both, the
//! unsigned outcome for the stop and the signed one for the cell's next value. The loop above
//! takes two compares (`diff - wordsize` with `k`, `diff` with `window`) and a copy of `diff`.
//! `diff - wordsize` is exact:
//! - on a packed table (`u32` cells, `XDiagCell::PACKED`) in 32 bits. A flag-0 cell holds
//!   last_hit in [0, 2^31) (31 bits), and the table keeps `s + diag_offset` in [window, 2^31)
//!   (`search_init::x_packs`: offsets below INT4_MAX / 4 + window + the longest unit, a window of
//!   at most 2^28, below 8.3 x 10^8 in all). So `diff` lies in (window - 2^31, 2^31) and
//!   `diff - wordsize` above `window - wordsize - 2^31 >= 1 - 2^31`: neither wraps. A flag-1 cell
//!   stops the loop whatever `dw` is.
//! - on a `(last_hit, flag)` table (any window, offsets that may wrap) in 64 bits: `diff` is the
//!   reference's wrapped `i32`, and `diff as i64 - wordsize` never wraps.
//!
//! Two hits a step: the step's control (one add, one decrement and its branch) is shared by the
//! two hits, which are taken one after the other in order, as above. The hits after the last whole
//! step are taken one by one.

use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:51-55
/// A cell of the two-hit diagonal table (`search_init::XDiagCells`): `u32`, NCBI's `DiagStruct` in
/// its four bytes (`LOSAT_X_TBNDIAGPACK`: last_hit in bits 0-30, flag in bit 31; the table's
/// values fit 31 bits without a sign, `search_init::x_packs`), or `(last_hit, flag)`.
pub(super) trait XDiagCell: Copy {
    /// Whether the table is packed (`u32` cells): a flag-0 cell holds last_hit in [0, 2^31), and
    /// the table keeps `s + diag_offset` in [window, 2^31) (`search_init::x_packs`), so
    /// `x_hits_steps` takes `diff - wordsize` in 32 bits (module doc).
    const PACKED: bool;
    /// `(last_hit, flag)`.
    fn get(self) -> (i32, bool);
    /// Stores both `last_hit` and `flag`.
    fn set(&mut self, last_hit: i32, flag: bool);
    /// Stores `last_hit` on a cell whose flag is 0.
    fn set_last_hit(&mut self, last_hit: i32);
    /// `(value, flag)` for `x_hits_bf`: `value` is `last_hit` when `flag` is false.
    fn bf_read(self) -> (i32, bool);
}

impl XDiagCell for (i32, bool) {
    const PACKED: bool = false;

    #[inline(always)]
    fn get(self) -> (i32, bool) {
        self
    }

    #[inline(always)]
    fn set(&mut self, last_hit: i32, flag: bool) {
        *self = (last_hit, flag);
    }

    #[inline(always)]
    fn set_last_hit(&mut self, last_hit: i32) {
        self.0 = last_hit;
    }

    #[inline(always)]
    fn bf_read(self) -> (i32, bool) {
        self
    }
}

impl XDiagCell for u32 {
    const PACKED: bool = true;

    #[inline(always)]
    fn get(self) -> (i32, bool) {
        ((self & 0x7fff_ffff) as i32, self >> 31 != 0)
    }

    #[inline(always)]
    fn set(&mut self, last_hit: i32, flag: bool) {
        debug_assert!(last_hit >= 0, "a packed cell holds 31 bits without a sign");
        *self = last_hit as u32 | u32::from(flag) << 31;
    }

    #[inline(always)]
    fn set_last_hit(&mut self, last_hit: i32) {
        debug_assert!(last_hit >= 0, "a packed cell holds 31 bits without a sign");
        *self = last_hit as u32;
    }

    // A set flag is the sign bit; a clear one leaves last_hit, which is never negative.
    #[inline(always)]
    fn bf_read(self) -> (i32, bool) {
        (self as i32, (self as i32) < 0)
    }
}

/// Whether `x_hits_bf` gives the reference's outcomes for this
/// window and word size (module doc): `0 <= wordsize < window`.
pub(super) fn fits(window: i32, word_size: usize) -> bool {
    i32::try_from(word_size).is_ok_and(|word_size| word_size < window)
}

// No NCBI counterpart: reads LOSAT_X_TBNHITSTEP once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNHITSTEP): whether the two-hit loop of the units is `x_hits_steps` (module
/// doc). A switch of the strict set, on unless it is `0`.
pub(super) fn x_hit_steps() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNHITSTEP"))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-552 (module doc).
/// `x_hits_bf`, or `x_hits_steps` when `STEPS` (LOSAT_X_TBNHITSTEP): the same hits taken, the same
/// cells left.
#[inline(always)]
pub(super) fn x_hits<C: XDiagCell, const STEPS: bool>(
    entries: &mut [C],
    offset: i32,
    window: i32,
    word_size: i32,
    pairs: &[BlastOffsetPair],
) -> usize {
    if STEPS {
        x_hits_steps(entries, offset, window, word_size, pairs)
    } else {
        x_hits_bf(entries, offset, window, word_size, pairs)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-552 (module doc).
/// One hit of `x_hits_steps` on the cell `entry`: false (nothing written) for a flag-1 cell or a
/// candidate, else the cell's next `last_hit` written and true. `k` is `window - wordsize`.
///
/// # Safety
/// `entry` is a cell of the table.
#[inline(always)]
unsafe fn x_hit_step<C: XDiagCell>(
    entry: *mut C,
    s_off: u32,
    offset: i32,
    word_size: i32,
    k: u32,
) -> bool {
    // SAFETY: (caller) a cell of the table.
    let (last_hit, flag) = unsafe { entry.read() }.bf_read();
    // Unbias the stored `last_hit` by `diag_offset` and compute distance from `subject_offset`.
    let s_diag = (s_off as i32).wrapping_add(offset);
    let diff = s_diag.wrapping_sub(last_hit);
    // diff - wordsize, exact in either width (module doc): below k unsigned = a candidate, at
    // least k signed = too far.
    let (candidate, too_far) = if C::PACKED {
        let dw = diff.wrapping_sub(word_size);
        ((dw as u32) < k, dw >= k as i32)
    } else {
        let dw = i64::from(diff) - i64::from(word_size);
        ((dw as u64) < u64::from(k), dw >= i64::from(k))
    };
    // Stop on a flagged diagonal or when a hit passes both distance tests.
    if flag | candidate {
        return false;
    }
    // Update `last_hit` when distance reaches or exceeds `window`; skip overlapping hits below `wordsize`.
    let next = if too_far { s_diag } else { last_hit };
    // SAFETY: as above.
    unsafe { (*entry).set_last_hit(next) };
    true
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-552 (module doc).
/// EXPERIMENT (LOSAT_X_TBNHITSTEP): `x_hits_bf` with each hit's outcomes from one compare
/// (`x_hit_step`), two hits a step (module doc). The same hits taken, the same cells left; the
/// same conditions (`fits(window, word_size)`, `entries.len()` a power of two).
#[inline(always)]
pub(super) fn x_hits_steps<C: XDiagCell>(
    entries: &mut [C],
    offset: i32,
    window: i32,
    word_size: i32,
    pairs: &[BlastOffsetPair],
) -> usize {
    debug_assert!(entries.is_empty() || entries.len().is_power_of_two());
    debug_assert!(0 <= word_size && word_size < window);
    let Some(last) = entries.len().checked_sub(1) else {
        return 0;
    };
    // Diagonal table index: query offset minus subject offset, masked by length minus 1.
    let mask = last as u32;
    // window - wordsize: the length of the candidate range of diff - wordsize.
    let k = window.wrapping_sub(word_size) as u32;
    let base = entries.as_mut_ptr();
    let start = pairs.as_ptr();
    let mut p = start;
    // SAFETY: every pair read is one of `pairs` (the steps cover `2 * steps <= pairs.len()` of
    // them, the rest loop stops at the end); every cell index is at most `mask` =
    // `entries.len() - 1`.
    unsafe {
        let end = start.add(pairs.len());
        let mut steps = pairs.len() / 2;
        while steps != 0 {
            for j in 0..2 {
                let pair = p.add(j).read();
                let entry = base.add((pair.q_off.wrapping_sub(pair.s_off) & mask) as usize);
                if !x_hit_step(entry, pair.s_off, offset, word_size, k) {
                    return p.add(j).offset_from(start) as usize;
                }
            }
            p = p.add(2);
            steps -= 1;
        }
        while p < end {
            let pair = p.read();
            let entry = base.add((pair.q_off.wrapping_sub(pair.s_off) & mask) as usize);
            if !x_hit_step(entry, pair.s_off, offset, word_size, k) {
                break;
            }
            p = p.add(1);
        }
        p.offset_from(start) as usize
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-552 (module doc).
/// The hits of `pairs` from the first, each on its cell of
/// `entries` (`(last_hit, flag)`, `XDiagCell`), up to the first that is a flag-1 cell or a
/// candidate for the extension; returns the number of hits taken (`pairs.len()` when there is no
/// such hit). The hits taken are too far (`last_hit = s + offset`) or overlap (no change). The
/// outcomes are exact when `fits(window, word_size)` (module doc); `entries.len()` is a power of
/// two.
#[inline(always)]
pub(super) fn x_hits_bf<C: XDiagCell>(
    entries: &mut [C],
    offset: i32,
    window: i32,
    word_size: i32,
    pairs: &[BlastOffsetPair],
) -> usize {
    debug_assert!(entries.is_empty() || entries.len().is_power_of_two());
    debug_assert!(0 <= word_size && word_size < window);
    let Some(last) = entries.len().checked_sub(1) else {
        return 0;
    };
    // Diagonal table index: query offset minus subject offset, masked by length minus 1.
    let mask = last as u32;
    // window - wordsize: the length of the candidate range of diff - wordsize.
    let k = window.wrapping_sub(word_size) as u32;
    let base = entries.as_mut_ptr();
    let start = pairs.as_ptr();
    let mut p = start;
    // SAFETY: one past the last element of `pairs`.
    let end = unsafe { start.add(pairs.len()) };
    while p < end {
        // SAFETY: start <= p < end.
        let pair = unsafe { p.read() };
        let index = (pair.q_off.wrapping_sub(pair.s_off) & mask) as usize;
        // SAFETY: index <= mask = entries.len() - 1.
        let entry = unsafe { base.add(index) };
        // SAFETY: as above.
        let (last_hit, flag) = unsafe { entry.read() }.bf_read();
        // Unbias the stored `last_hit` by `diag_offset` and compute distance from `subject_offset`.
        let s_diag = (pair.s_off as i32).wrapping_add(offset);
        let diff = s_diag.wrapping_sub(last_hit);
        let candidate = (diff.wrapping_sub(word_size) as u32) < k;
        // Stop on a flagged diagonal or when a hit passes both distance tests.
        if flag | candidate {
            break;
        }
        // Update `last_hit` when distance reaches or exceeds `window`; skip overlapping hits below `wordsize`.
        let next = if diff >= window { s_diag } else { last_hit };
        // SAFETY: as above.
        unsafe { (*entry).set_last_hit(next) };
        // SAFETY: p < end, so p + 1 <= end.
        p = unsafe { p.add(1) };
    }
    // SAFETY: start <= p <= end, both in `pairs`.
    unsafe { p.offset_from(start) as usize }
}

#[cfg(test)]
#[path = "x_tbn_hitbf_tests.rs"]
mod tests;
