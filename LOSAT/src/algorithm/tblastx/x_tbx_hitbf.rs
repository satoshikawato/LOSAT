//! EXPERIMENT (LOSAT_X_TBXHITBF): the hit loop of a unit of the unit path
//! (x_pair_par.rs `scan_unit`, x_tbx_extbatch.rs `scan_unit_batched`, and the scan in slices of
//! LOSAT_X_TBXQSLICE, x_tbx_qslice.rs `scan_unit_sliced`) with the too-far and overlap outcomes of
//! the two-hit test taken without a branch. A child module of x_pair_par.rs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:518-552
//! In the reference loop most hits end in one of these two `continue`s (on AP027131 x AP027133,
//! 80 % too far and 16 % overlap), and which one is taken depends on the data: the too-far test
//! was mispredicted 0.14 times per hit there (v4-R6 REPORT, workstation). `x_hits_bf` computes
//! both tests from one subtraction and writes the cell's next value, chosen without a branch:
//! `s + window` (too far) or the value it read (overlap: the same value back). A flag-1 cell or a
//! hit that passes both tests (a candidate) goes to the per-hit body of the scan (`x_hit`, or
//! `x_hit_batched`), a branch that is rarely taken. The loop walks the batch by pointer, and the
//! prefetch of LOSAT_X_TBXPREFETCH (`PF`) reads the entry `X_PREFETCH_AHEAD` ahead without a bound
//! test (the offset-pair vector has that many entries more; `scan_unit_hits`).
//!
//! Why the result is that of the reference loop. Let L be a cell's `last_hit` (the 31-bit field,
//! sign-extended) and s the hit's subject offset. A unit's table starts at raw 0 with
//! `diag_offset == window` (x_pair_par.rs module doc, item 1), so in the flag-0 branch
//! `diff = s - L + window`, and in exact integers
//! - too far: `diff >= window` <=> `s - L >= 0`;
//! - overlap: `diff < wordsize` <=> `s - L < wordsize - window`;
//! - otherwise `0 <= diff - wordsize < window - wordsize`, the candidate.
//!
//! The rule `fits` (checked per unit): `1 <= wordsize < window` and chunk length + window < 2^29.
//! Every s is below the chunk length. Every value a unit stores with flag 0 is 0 (the fresh cell)
//! or `s + window` (too far, the flag reset, the context check, an extension without a right
//! extension), so a flag-0 cell holds L in [0, 2^29) as its raw word (bit 31 clear), with nothing
//! cut by the 31-bit field. The values stored with flag 1 (a right extension; LOSAT_X_TBXEXTBATCH's
//! pending marks) are never decoded here: the sign of the raw word sends them to the per-hit body
//! first. So for a flag-0 cell `d = s - raw` (wrapping i32) is exactly s - L, |d| < 2^29, and with
//! `k = window - wordsize` (in [1, 2^29)) `t = d + k` is exactly `diff - wordsize`, in
//! (-2^29, 2^30):
//! - too far <=> `d >= 0`;
//! - overlap <=> `t < 0`;
//! - candidate <=> `0 <= t < k` <=> `t as u32 < k`: a negative t is 2^31 or more as u32.
//!
//! The reference's too-far write stores `s + diag_offset` truncated to 31 bits with the flag kept
//! (0); `s + window` is below 2^29, so the raw word `s + window` is that value. The overlap path
//! writes back the word it read. A flag-1 cell or a candidate goes to the per-hit body unchanged,
//! which then reads and writes the cell as in the reference. The loop reads and writes only the cell
//! of the hit at hand, after every call of the per-hit body for the hits before it, and keeps no
//! cell value across such a call (the batched body writes other cells too: it resolves held-back
//! extensions, x_tbx_extbatch.rs). So every cell goes through the reference's states in scan order,
//! every call of the per-hit body sees the table the loop before the switch gives it, and the init
//! HSPs are saved by the same calls in the same order. The scan in slices (x_tbx_qslice.rs) runs the
//! slices of a unit one after the other on the unit's one table, set up as `scan_unit` sets it up
//! (the group's bounds: the same `lo` and width) and not cleared between slices; its hits are those
//! of the unit's chunk, so the above holds for them in the order the slices give them, against the
//! scan in slices with the reference loop. A unit outside `fits` runs the reference loop, which
//! keeps the 31-bit field's truncation (memory `losat-no-ncbi-bug-reproduction`, 2026-10-12); the
//! batch guard `last_hit_fits` already keeps the units of the unit path inside it.
//!
//! The rule: the switch is in the strict set (on unless `LOSAT_X_TBXHITBF` is `0`);
//! the probes of LOSAT_TIMING / LOSAT_X_STATS run the
//! reference loop (they count every outcome).
//!
//! On wasm32 (LOSAT Web) `x_hits_bf` runs as `x_hits_bf_split`: a function of its own,
//! `x_hits_bf_run`, takes the hits up to the first one for the per-hit body (flag 1 or a candidate)
//! and returns it with its cell not written; the caller runs the per-hit body on it and goes on
//! after it. The same reads and writes of the cells and the same calls of the per-hit body, in the
//! same order (the caller reads that hit and computes its cell a second time), so the same states
//! and calls as above. V8's optimized code (TurboFan) kept the table, the bounds and the loop
//! pointer of the one loop in stack slots (the inlined per-hit body needs every register) and read
//! or wrote six of them for every hit; no value of `x_hits_bf_run`'s loop lives across a call.
//! Native builds keep the one loop.

use super::x_tbx_qslice::XSliced;
use super::*;

/// `fits`: chunk length + window stays below this (2^29; module doc).
const X_HITBF_LIMIT: i64 = 1 << 29;

// No NCBI counterpart: which hit loop the units of the unit path use (scheduling only). Read once
// per process, as the other modes of x_pair_par.rs.
/// LOSAT_X_TBXHITBF: a strict-set switch (`utils/xswitch.rs`).
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXHITBF"))
}

// No NCBI counterpart: the range in which `x_hits_bf` takes the reference's outcomes (module doc).
/// Whether a unit with a subject chunk of `chunk_length` letters may run `x_hits_bf`:
/// `1 <= wordsize < window` and `chunk_length + window < 2^29`.
pub(super) fn fits(chunk_length: usize, window: i32, wordsize: i32) -> bool {
    wordsize >= 1
        && window > wordsize
        && chunk_length < X_HITBF_LIMIT as usize
        && (chunk_length as i64) + (window as i64) < X_HITBF_LIMIT
}

// No NCBI counterpart: the switch `on` for one unit (module doc).
/// Whether a unit runs the branch-free loop: the switch is on and the unit is inside `fits`
/// (otherwise the reference loop).
fn unit_bf(on: bool, chunk_length: usize, window: i32, wordsize: i32) -> bool {
    if !on {
        return false;
    }
    let inside = fits(chunk_length, window, wordsize);
    crate::utils::xreach::x_reach!("LOSAT_X_TBXHITBF", inside);
    inside
}

// No NCBI counterpart: records that a batch left with the reference loop (the unit path does not
// run for it) has the switch on and does not use it (`LOSAT_X_REACH`).
/// LOSAT_X_TBXHITBF: a batch the unit path does not take runs the reference loop.
pub(super) fn x_reach_not_taken() {
    if on() {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXHITBF", false);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-552 (module doc)
// The hit loop of `scan_unit` over one batch: the cell of each hit as there (the raw diagonal minus
// the unit's lowest one), the too-far and overlap outcomes of the flag-0 branch from one
// subtraction with the cell's next value selected (module doc), and `rare` for a flag-1 cell or a
// candidate. One constant (`k`) besides the window, so that both can stay in registers.
/// The two-hit test of the `n` hits at `pairs` on the unit table at `diag_ptr` (cell of a hit
/// (q, s): `cells.cell(q, s)`, q - s - lo or, with LOSAT_X_TBXUNITMASK, that masked), with
/// `rare(cell, q, s)` for a hit whose cell has flag 1 or that passes the test; `PF` prefetches the
/// cell of the hit `X_PREFETCH_AHEAD` entries ahead (x86-64).
///
/// SAFETY: `pairs` points to `n` entries, and with `PF` to `X_PREFETCH_AHEAD` initialized entries
/// after them; every hit's cell is below `width`, the length of the table `diag_ptr` points into;
/// `fits(chunk length, window, wordsize)` holds, every subject offset is below the chunk length and
/// every flag-0 cell holds a value the unit stored (module doc); `rare` writes only cells of the
/// table, values the unit's per-hit body stores (the loop keeps no cell value across a call).
#[cfg(not(target_arch = "wasm32"))]
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) unsafe fn x_hits_bf<C: XCells, const PF: bool>(
    diag_ptr: *mut DiagStruct,
    pairs: *const OffsetPair,
    n: usize,
    cells: C,
    width: usize,
    window: i32,
    wordsize: i32,
    mut rare: impl FnMut(usize, u32, u32),
) {
    let window_u = window as u32;
    // window - wordsize, at least 1 (`fits`): the length of the candidate range of diff - wordsize.
    let k = window.wrapping_sub(wordsize) as u32;
    let mut p = pairs;
    // SAFETY: the caller's bound (n entries).
    let end = unsafe { pairs.add(n) };
    while p < end {
        // EXPERIMENT (LOSAT_X_TBXPREFETCH): the cell of the hit X_PREFETCH_AHEAD entries ahead,
        // computed as below; past the batch the entry is a stale one and the address is only
        // prefetched (no architectural effect, `wrapping_add`).
        #[cfg(target_arch = "x86_64")]
        if PF {
            // SAFETY: the caller's bound (X_PREFETCH_AHEAD initialized entries after the n).
            let ahead = unsafe { p.add(X_PREFETCH_AHEAD).read() };
            let cell = cells.cell(ahead.q_off, ahead.s_off);
            // SAFETY: a prefetch never faults.
            unsafe {
                core::arch::x86_64::_mm_prefetch(
                    diag_ptr.wrapping_add(cell) as *const i8,
                    core::arch::x86_64::_MM_HINT_T0,
                )
            };
        }
        // SAFETY: p < end.
        let OffsetPair {
            q_off: query_offset,
            s_off: subject_offset,
        } = unsafe { p.read() };
        // SAFETY: p < end, so p + 1 <= end.
        p = unsafe { p.add(1) };
        // Map the hit diagonal (`query_offset` minus `subject_offset`) into the unit table.
        // Unit cell: the raw diagonal minus the unit's lowest one (`XUnmasked`: no mask;
        // LOSAT_X_TBXUNITMASK's `XMasked`: masked to the table).
        let cell = cells.cell(query_offset, subject_offset);
        debug_assert!(cell < width, "unit diagonal outside its table");
        // SAFETY: cell < width (the caller's bound; a masked cell is at most the mask, width - 1).
        let entry = unsafe { diag_ptr.add(cell) };
        // SAFETY: as above.
        let raw = unsafe { (*entry).raw_bits() };
        // s - L and diff - wordsize for a flag-0 cell (module doc); any values for a flag-1 cell,
        // which goes to `rare` below.
        let d = subject_offset.wrapping_sub(raw) as i32;
        let candidate = (d as u32).wrapping_add(k) < k;
        // Dispatch flagged cells or hits passing both distance tests to `rare`.
        if ((raw as i32) < 0) | candidate {
            rare(cell, query_offset, subject_offset);
            continue;
        }
        // Update `last_hit` when outside the window, or keep the existing entry on an overlapping hit.
        let next = if d >= 0 {
            subject_offset.wrapping_add(window_u)
        } else {
            raw
        };
        // SAFETY: as above.
        unsafe { (*entry).set_raw_bits(next) };
    }
}

// No NCBI counterpart (scheduling only): `x_hits_bf` on wasm32 (LOSAT Web), the same statements in
// two parts (`x_hits_bf_split`).
/// `x_hits_bf` on wasm32: `x_hits_bf_split` (`PF` prefetches on x86-64 only).
///
/// SAFETY: as `x_hits_bf` above.
#[cfg(target_arch = "wasm32")]
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) unsafe fn x_hits_bf<C: XCells, const PF: bool>(
    diag_ptr: *mut DiagStruct,
    pairs: *const OffsetPair,
    n: usize,
    cells: C,
    width: usize,
    window: i32,
    wordsize: i32,
    rare: impl FnMut(usize, u32, u32),
) {
    // SAFETY: the caller's bounds, passed on.
    unsafe { x_hits_bf_split(diag_ptr, pairs, n, cells, width, window, wordsize, rare) }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:510-552 (as `x_hits_bf`)
// No NCBI counterpart for the function boundary (scheduling of the statements of `x_hits_bf`
// only): on wasm32 (LOSAT Web) the loop of `x_hits_bf` runs as `x_hits_bf_run`, a function of its
// own that takes the hits up to the first one for `rare`, and this loop, which calls `rare` for
// that hit and goes on after it. In `x_hits_bf`'s one loop, V8's optimized code (TurboFan) kept
// the table, the bounds and the loop pointer in stack slots and read or wrote six of them for every
// hit, as the inlined per-hit body (`x_hit`, its extension and its calls) needs every register. No
// value of `x_hits_bf_run`'s loop lives across a call: V8 inlines it into this loop and unrolls it,
// without stack slots. Elsewhere `x_hits_bf` keeps its one loop.
/// `x_hits_bf` (same arguments, same safety conditions, without `PF`): the same reads and writes of
/// the cells and the same calls of `rare`, in the same order.
#[cfg(any(target_arch = "wasm32", test))]
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn x_hits_bf_split<C: XCells>(
    diag_ptr: *mut DiagStruct,
    pairs: *const OffsetPair,
    n: usize,
    cells: C,
    width: usize,
    window: i32,
    wordsize: i32,
    mut rare: impl FnMut(usize, u32, u32),
) {
    let window_u = window as u32;
    // window - wordsize, at least 1 (`fits`), as in `x_hits_bf`.
    let k = window.wrapping_sub(wordsize) as u32;
    // SAFETY: the caller's bound (n entries).
    let end = unsafe { pairs.add(n) };
    let mut p = pairs;
    loop {
        // SAFETY: pairs <= p <= end (the bounds `x_hits_bf_run` keeps).
        p = unsafe { x_hits_bf_run(diag_ptr, p, end, cells, width, window_u, k) };
        if p == end {
            break;
        }
        // SAFETY: p < end: a hit `x_hits_bf_run` left for `rare` (its cell neither read into a
        // value kept nor written).
        let OffsetPair {
            q_off: query_offset,
            s_off: subject_offset,
        } = unsafe { p.read() };
        // SAFETY: p < end, so p + 1 <= end.
        p = unsafe { p.add(1) };
        // Map the hit diagonal (`query_offset` minus `subject_offset`) into the unit table.
        let cell = cells.cell(query_offset, subject_offset);
        debug_assert!(cell < width, "unit diagonal outside its table");
        // Dispatch flagged cells or hits passing both distance tests to `rare`.
        rare(cell, query_offset, subject_offset);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:518-552 (module doc)
// The flag-0 branch's too-far and overlap outcomes of `x_hits_bf`, statement for statement, up to
// the first hit that goes to `rare` (`x_hits_bf_split`).
/// From `p` up to `end`: each hit whose cell has flag 0 and that does not pass the test gets the
/// cell's next value (`s + window`, or the word it read), as in `x_hits_bf`. Returns the first hit
/// whose cell has flag 1 or that passes the test, its cell not written, or `end`.
///
/// SAFETY: `x_hits_bf`'s conditions for the hits from `p` to `end` (`p <= end`, `window_u` the
/// window, `k` = window - wordsize).
#[cfg(any(target_arch = "wasm32", test))]
#[inline(never)]
unsafe fn x_hits_bf_run<C: XCells>(
    diag_ptr: *mut DiagStruct,
    mut p: *const OffsetPair,
    end: *const OffsetPair,
    cells: C,
    width: usize,
    window_u: u32,
    k: u32,
) -> *const OffsetPair {
    while p < end {
        // SAFETY: p < end.
        let OffsetPair {
            q_off: query_offset,
            s_off: subject_offset,
        } = unsafe { p.read() };
        // Map the hit diagonal (`query_offset` minus `subject_offset`) into the unit table.
        let cell = cells.cell(query_offset, subject_offset);
        debug_assert!(cell < width, "unit diagonal outside its table");
        // SAFETY: cell < width (the caller's bound; a masked cell is at most the mask, width - 1).
        let entry = unsafe { diag_ptr.add(cell) };
        // SAFETY: as above.
        let raw = unsafe { (*entry).raw_bits() };
        // s - L and diff - wordsize for a flag-0 cell (module doc).
        let d = subject_offset.wrapping_sub(raw) as i32;
        let candidate = (d as u32).wrapping_add(k) < k;
        // Return to the caller (`rare`) on flagged cells or hits passing both distance tests.
        if ((raw as i32) < 0) | candidate {
            return p;
        }
        // Update `last_hit` when outside the window, or keep the existing entry on an overlapping hit.
        let next = if d >= 0 {
            subject_offset.wrapping_add(window_u)
        } else {
            raw
        };
        // SAFETY: as above.
        unsafe { (*entry).set_raw_bits(next) };
        // SAFETY: p < end, so p + 1 <= end.
        p = unsafe { p.add(1) };
    }
    end
}

impl<M: TbxScore> XPairPar<M> {
    // EXPERIMENT (LOSAT_X_TBXHITBF): dispatch of the units without probes (`run_unit`,
    // `run_unit_compact`) by the switch and `fits`. EXPERIMENT (LOSAT_X_TBXUNGAPSENT): the
    // branch-free loop's units by `x_unit_sent`, between the instances of `scan_unit_bf` (`SENT`:
    // the search's word length, the walks `TbxScore::sent_two_hit_of` of 2 or 3 letters, or 0),
    // so that the instance without the switch keeps the code it had before it.
    /// `scan_unit` without probes, or its branch-free hit loop.
    #[inline(always)]
    pub(super) fn scan_unit_or_bf<S: XInitSink, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        if unit_bf(on(), unit.chunk.length, inp.window, inp.wordsize) {
            if self.x_unit_sent(unit, inp) {
                // The walks of the search's word length (`XPairPar::sent`: 2 or 3 letters); 2
                // only in another instance than BLOSUM62's (whose words have 3).
                if !M::BLOSUM62 && self.sc.word() == 2 {
                    self.scan_unit_bf::<S, PF, 2>(unit, inp, scratch, init_hsps);
                } else {
                    assert!(M::BLOSUM62 || self.sc.word() == 3);
                    self.scan_unit_bf::<S, PF, 3>(unit, inp, scratch, init_hsps);
                }
            } else {
                self.scan_unit_bf::<S, PF, 0>(unit, inp, scratch, init_hsps);
            }
        } else {
            self.scan_unit::<S, XNoProbe, PF>(unit, inp, scratch, init_hsps, &mut XNoProbe);
        }
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF): its own function, so that the loop gets its own registers; the
    // units outside the switch run the instance of `scan_unit` without it.
    /// `scan_unit` with the branch-free hit loop.
    #[inline(never)]
    fn scan_unit_bf<S: XInitSink, const PF: bool, const SENT: usize>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        self.scan_unit_hits::<S, XNoProbe, PF, true, SENT>(
            unit,
            inp,
            scratch,
            init_hsps,
            &mut XNoProbe,
        );
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF): dispatch of the batched scan (LOSAT_X_TBXEXTBATCH) by the
    // switch and `fits`; EXPERIMENT (LOSAT_X_TBXUNGAPSENT): and by `x_unit_sent`, as
    // `scan_unit_or_bf`.
    /// `scan_unit_batched` with full batches, or its branch-free hit loop.
    #[inline(always)]
    pub(super) fn scan_unit_batched_or_bf<S: XInitSink, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        if unit_bf(on(), unit.chunk.length, inp.window, inp.wordsize) {
            if self.x_unit_sent(unit, inp) {
                // The walks of the search's word length (`XPairPar::sent`: 2 or 3 letters); 2
                // only in another instance than BLOSUM62's (whose words have 3).
                if !M::BLOSUM62 && self.sc.word() == 2 {
                    self.scan_unit_batched_bf::<S, PF, 2>(unit, inp, scratch, init_hsps);
                } else {
                    assert!(M::BLOSUM62 || self.sc.word() == 3);
                    self.scan_unit_batched_bf::<S, PF, 3>(unit, inp, scratch, init_hsps);
                }
            } else {
                self.scan_unit_batched_bf::<S, PF, 0>(unit, inp, scratch, init_hsps);
            }
        } else {
            self.scan_unit_batched::<S, PF>(unit, inp, scratch, init_hsps, x_tbx_extbatch::BATCH);
        }
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF): its own function, as `scan_unit_bf`.
    /// `scan_unit_batched` with full batches and the branch-free hit loop.
    #[inline(never)]
    fn scan_unit_batched_bf<S: XInitSink, const PF: bool, const SENT: usize>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        self.scan_unit_batched_hits::<S, PF, true, SENT>(
            unit,
            inp,
            scratch,
            init_hsps,
            x_tbx_extbatch::BATCH,
        );
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF): dispatch of the scan in slices (LOSAT_X_TBXQSLICE) by the
    // switch and `fits`; EXPERIMENT (LOSAT_X_TBXUNGAPSENT): and by `x_unit_sent`, as
    // `scan_unit_or_bf`.
    /// `scan_unit_sliced` (batched extensions with `B`, the prefetch with `PF`), or its
    /// branch-free hit loop.
    #[inline(always)]
    pub(super) fn scan_unit_sliced_or_bf<S: XInitSink, const B: bool, const PF: bool>(
        &self,
        group: &XSliced,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        if unit_bf(on(), unit.chunk.length, inp.window, inp.wordsize) {
            if self.x_unit_sent(unit, inp) {
                // The walks of the search's word length (`XPairPar::sent`: 2 or 3 letters); 2
                // only in another instance than BLOSUM62's (whose words have 3).
                if !M::BLOSUM62 && self.sc.word() == 2 {
                    self.scan_unit_sliced_bf::<S, B, PF, 2>(group, unit, inp, scratch, init_hsps);
                } else {
                    assert!(M::BLOSUM62 || self.sc.word() == 3);
                    self.scan_unit_sliced_bf::<S, B, PF, 3>(group, unit, inp, scratch, init_hsps);
                }
            } else {
                self.scan_unit_sliced_bf::<S, B, PF, 0>(group, unit, inp, scratch, init_hsps);
            }
        } else {
            self.scan_unit_sliced::<S, B, PF>(group, unit, inp, scratch, init_hsps);
        }
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF): its own function, as `scan_unit_bf`.
    /// `scan_unit_sliced` with the branch-free hit loop.
    #[inline(never)]
    fn scan_unit_sliced_bf<S: XInitSink, const B: bool, const PF: bool, const SENT: usize>(
        &self,
        group: &XSliced,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        self.scan_unit_sliced_hits::<S, B, PF, true, SENT>(group, unit, inp, scratch, init_hsps);
    }
}

#[cfg(test)]
#[path = "x_tbx_hitbf_tests.rs"]
mod tests;
