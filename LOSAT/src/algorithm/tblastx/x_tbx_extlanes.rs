//! EXPERIMENT (LOSAT_X_TBXEXTLANES / LOSAT_X_TBXEXTLANESSHADOW): the two-hit stage of a unit of the
//! unit path (x_pair_par.rs) with its extensions held back and resolved in batches. A child module
//! of x_pair_par.rs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:575-605
//! ```c
//! score = s_BlastAaExtendTwoHit(matrix, subject, query,
//!                               last_hit + wordsize,
//!                               subject_offset, query_offset,
//!                               cutoffs->x_dropoff,
//!                               &hsp_q, &hsp_s,
//!                               &hsp_len, use_pssm,
//!                               wordsize, &right_extend,
//!                               &s_last_off);
//! ...
//! if (score >= cutoffs->cutoff_score)
//!     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
//!                      query_offset, subject_offset, hsp_len,
//!                      score);
//! ...
//! if (right_extend) {
//!     diag_array[diag_coord].flag = 1;
//!     diag_array[diag_coord].last_hit =
//!         s_last_off - (wordsize - 1) + diag_offset;
//! }
//! else {
//!     diag_array[diag_coord].last_hit =
//!         subject_offset + diag_offset;
//! }
//! ```
//! `x_hit` runs an extension as soon as a hit passes the two-hit test, so every extension starts
//! with loads at a query position the scan did not touch before (the query offset of a hit is a
//! random position), and its walks end on a mispredicted exit. Here a hit that passes the test
//! becomes a candidate: its inputs are recorded in a slot of a batch, the cache lines its
//! extension reads first are prefetched, and its cell gets a pending mark. The batch is resolved
//! when it holds `batch()` candidates and at the end of the unit; until then the scan goes on.
//!
//! The pending mark is `flag = 1` with the slot number (`0..LANES_MAX`) in `last_hit`. A real cell
//! with `flag = 1` holds `s_last_off - (wordsize - 1) + diag_offset`, and the unit's diag_offset
//! is the window, so with `s_last_off >= 0` its `last_hit` is at least `window - wordsize + 1`;
//! the module runs only where `window >= wordsize + LANES_MAX` (`active`), so a real flag-1 cell
//! never reads as a mark. The unit path's guard `last_hit_fits` keeps every stored value below
//! 2^29, so no value wraps into the slot range either. The mark is tested only in the `flag = 1`
//! branch of the per-hit body.
//!
//! Why the result is that of `scan_unit`.
//! 1. `x_hit` reads and writes only the cell of its hit; the extension reads only the two
//!    sequences and the matrix, and the read-only values of the unit (`XHitCtx`) are fixed. So an
//!    extension may run at any time before the next hit on its cell is tested, with the same
//!    result, provided the cell then holds what `x_hit` would have written.
//! 2. A hit on a pending cell resolves that one candidate first (the extension and the cell
//!    write of `x_hit`), then takes the per-hit body on the cell as written. A cell holds at most
//!    one pending candidate (a second hit on it resolves the first before it can record another),
//!    so every cell goes through the states of `scan_unit`, in its order.
//! 3. The init HSPs are saved in slot order when the batch is resolved, and the slots are filled
//!    in scan order, never reused within a batch; a batch is resolved before the next one starts.
//!    So the unit's list is that of `scan_unit`, element for element.
//!
//! Step 1 (this version): the batch resolves its candidates one after the other with
//! `extend_hit_two_hit`, the function `x_hit` calls (LOSAT_X_TBXUNGAP's table walks where that
//! switch is on); only the prefetch at record time and the place of the call change. The probes of
//! LOSAT_TIMING / LOSAT_X_STATS run `scan_unit` (the module has none).
//!
//! `LOSAT_X_TBXEXTLANESSHADOW=1` runs `scan_unit` and this path on every unit of the module (each
//! on a fresh table) and stops at the first difference of the saved records or of the final table.

use super::*;

// No NCBI counterpart: counters of LOSAT_X_TBXEXTLANES (units, candidates, early resolutions,
// batches) and of its shadow; they do not change any value NCBI computes.
static EL_UNITS: AtomicU64 = AtomicU64::new(0);
static EL_CANDS: AtomicU64 = AtomicU64::new(0);
static EL_EARLY: AtomicU64 = AtomicU64::new(0);
static EL_BATCHES: AtomicU64 = AtomicU64::new(0);
static EL_SHADOW_UNITS: AtomicU64 = AtomicU64::new(0);
static EL_SHADOW_RECORDS: AtomicU64 = AtomicU64::new(0);

/// The most candidates a batch holds; a pending mark keeps the slot number, below this, in the
/// cell's `last_hit` (module doc).
pub(super) const LANES_MAX: usize = 16;
/// Candidates per batch when LOSAT_X_TBXEXTLANES_BATCH is not set.
const BATCH_DEFAULT: usize = 8;
/// The packed cell's flag bit (`DiagStruct`: last_hit in bits 0-30, flag in bit 31).
const FLAG_BIT: u32 = 0x8000_0000;
/// The packed cell's `last_hit` bits.
const LAST_HIT_BITS: u32 = 0x7fff_ffff;

// No NCBI counterpart: which scan the units of the unit path use (scheduling only).
/// 0 = `scan_unit`, 1 = LOSAT_X_TBXEXTLANES, 2 = LOSAT_X_TBXEXTLANESSHADOW (both, compared).
pub(super) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_TBXEXTLANESSHADOW").is_some() {
            2
        } else if std::env::var("LOSAT_X_TBXEXTLANES").is_ok_and(|v| v != "0") {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: a tuning value of the batches (prototype; it changes speed only).
/// Candidates per batch: LOSAT_X_TBXEXTLANES_BATCH (1 to `LANES_MAX`), else `BATCH_DEFAULT`.
fn batch() -> usize {
    static V: OnceLock<usize> = OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("LOSAT_X_TBXEXTLANES_BATCH")
            .ok()
            .and_then(|v| v.trim().parse().ok())
            .filter(|n| (1..=LANES_MAX).contains(n))
            .unwrap_or(BATCH_DEFAULT)
    })
}

// No NCBI counterpart: the rule of the module for one subject (module doc, the pending mark).
/// Whether the units of a subject with this two-hit `window` and `wordsize` take the module.
pub(super) fn active(window: i32, wordsize: i32) -> bool {
    mode() != 0 && window as i64 >= wordsize as i64 + LANES_MAX as i64
}

/// The summary line (with LOSAT_X_STATS, LOSAT_TIMING or the shadow).
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub(super) fn print_summary() {
    if mode() == 2 {
        eprintln!(
            "[X_SHADOW] TBXEXTLANES units_compared={} records_compared={} (records and unit tables identical)",
            EL_SHADOW_UNITS.load(AtomicOrdering::Relaxed),
            EL_SHADOW_RECORDS.load(AtomicOrdering::Relaxed)
        );
    }
    if mode() != 0 && (std::env::var_os("LOSAT_X_STATS").is_some() || timing()) {
        eprintln!(
            "[X_STATS] TBXEXTLANES units={} candidates={} early={} batches={} batch={}",
            EL_UNITS.load(AtomicOrdering::Relaxed),
            EL_CANDS.load(AtomicOrdering::Relaxed),
            EL_EARLY.load(AtomicOrdering::Relaxed),
            EL_BATCHES.load(AtomicOrdering::Relaxed),
            batch()
        );
    }
}

/// One held-back extension: the inputs of `extend_hit_two_hit` and of the init HSP, then its
/// result once resolved.
#[derive(Clone, Copy, Default)]
struct XCand {
    /// The query offset of the second hit (the concatenated query).
    query_offset: u32,
    /// `subject_offset` of the second hit (the chunk).
    s_right: u32,
    /// `last_hit + wordsize`.
    s_left: u32,
    /// The query offset of the second hit within its context.
    q_raw: u32,
    /// The hit's unit cell.
    cell: u32,
    /// The context of the query offset.
    ctx_idx: u32,
    /// The context's X-drop.
    x_drop: i32,
    /// Whether the extension below has run.
    done: bool,
    /// `hsp_q` (within the context), `hsp_s`, `hsp_len` and the score of the extension.
    hsp_q: i32,
    hsp_s: i32,
    hsp_len: i32,
    score: i32,
}

/// The candidates of the batch that is filling, and the unit's counters.
struct XLaneBatch {
    c: [XCand; LANES_MAX],
    /// Slots filled.
    n: usize,
    /// Slots of a full batch (`batch()`).
    len: usize,
    cands: u64,
    early: u64,
    batches: u64,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:575-605 (module doc)
// The extension of one candidate and the cell write of `x_hit` (the statements after its call of
// `extend_hit_two_hit`), the flag written with the value, since the cell holds the pending mark.
/// Runs the extension of `c` and writes its cell as `x_hit` does.
///
/// SAFETY: `c.cell` is below the length of the table `diag_ptr` points into.
#[inline(always)]
unsafe fn resolve(c: &mut XCand, diag_ptr: *mut DiagStruct, h: &XHitCtx<'_>) {
    // SAFETY: the context index was returned by `get_context_idx` (below num_contexts).
    let ctx = unsafe { h.contexts.get_unchecked(c.ctx_idx as usize) };
    let query_full = &ctx.aa_seq;
    let query = &query_full[1..query_full.len() - 1];
    let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) =
        extend_hit_two_hit(
            query,
            h.subject,
            c.s_left as usize,
            c.s_right as usize,
            c.q_raw as usize,
            c.x_drop,
            false,
        );
    // [C] if (right_extend) { flag = 1; last_hit = s_last_off - (wordsize - 1) + diag_offset; }
    // [C] else { last_hit = subject_offset + diag_offset; }  (flag stays 0)
    let raw = if right_extend {
        FLAG_BIT | ((s_last_off_u as i32 - (h.wordsize - 1) + h.diag_offset) as u32 & LAST_HIT_BITS)
    } else {
        c.s_right.wrapping_add(h.diag_offset as u32) & LAST_HIT_BITS
    };
    // SAFETY: the caller's bound.
    unsafe { (*diag_ptr.add(c.cell as usize)).set_raw_bits(raw) };
    c.hsp_q = hsp_q_u as i32;
    c.hsp_s = hsp_s_u as i32;
    c.hsp_len = (hsp_qe_u - hsp_q_u) as i32;
    c.score = score;
    c.done = true;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:575-605 (module doc)
// A hit on a pending cell: the candidate of the cell's slot first (module doc, item 2).
/// Resolves the candidate in `slot` before a later hit on its cell is tested.
///
/// SAFETY: as `resolve`.
#[cold]
#[inline(never)]
unsafe fn resolve_early(
    b: &mut XLaneBatch,
    slot: usize,
    diag_ptr: *mut DiagStruct,
    h: &XHitCtx<'_>,
) {
    debug_assert!(
        slot < b.n && !b.c[slot].done,
        "a pending mark without its candidate"
    );
    b.early += 1;
    // SAFETY: the caller's bound; the slot holds the candidate that wrote the mark.
    unsafe { resolve(b.c.get_unchecked_mut(slot), diag_ptr, h) };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:599-602
// ```c
// if (score >= cutoffs->cutoff_score)
//     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
//                      query_offset, subject_offset, hsp_len,
//                      score);
// ```
// The batch, slot by slot: the candidates not resolved yet, then the init HSPs in slot order
// (module doc, item 3), by the push statement of `x_hit`.
/// Resolves every candidate of the batch, saves their init HSPs in slot order and empties it.
/// Its own function, so that the scan loop keeps its values in registers.
///
/// SAFETY: every candidate's cell is below the length of the table `diag_ptr` points into.
#[inline(never)]
unsafe fn flush<S: XInitSink>(
    b: &mut XLaneBatch,
    diag_ptr: *mut DiagStruct,
    h: &XHitCtx<'_>,
    init_hsps: &mut S,
) {
    b.batches += 1;
    for c in &mut b.c[..b.n] {
        if !c.done {
            // SAFETY: the caller's bound.
            unsafe { resolve(c, diag_ptr, h) };
        }
        // SAFETY: ctx_idx < contexts.len() == cutoff_scores.len().
        let cutoff = unsafe { *h.cutoff_scores.get_unchecked(c.ctx_idx as usize) };
        if c.score >= cutoff {
            // SAFETY: as in `resolve`.
            let ctx = unsafe { h.contexts.get_unchecked(c.ctx_idx as usize) };
            init_hsps.push(InitHSP {
                q_start_absolute: ctx.frame_base + c.hsp_q,
                q_end_absolute: ctx.frame_base + (c.hsp_q + c.hsp_len),
                s_start: c.hsp_s,
                s_end: c.hsp_s + c.hsp_len,
                q_seed_absolute: c.query_offset as i32,
                s_seed: c.s_right as i32,
                score: c.score,
                ctx_idx: c.ctx_idx as usize,
                s_f_idx: h.s_f_idx,
                q_idx: ctx.q_idx,
                s_idx: h.s_idx as u32,
                q_frame: ctx.frame,
                s_frame: h.s_frame,
                q_orig_len: ctx.orig_len,
                s_orig_len: h.s_len,
            });
        }
    }
    b.n = 0;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:518-606 (at `x_hit`)
// The per-hit body of `x_hit`, statement for statement, with two changes (module doc): a hit on a
// pending cell resolves that candidate first, and a hit that passes the two-hit test is recorded
// in the batch, its first query and subject lines prefetched and its cell marked, instead of
// extended at once.
/// The two-hit test of one hit on its cell, recording the extension it starts.
///
/// SAFETY: `cell` and the cell of every candidate of `b` are below the length of the table
/// `diag_ptr` points into.
#[inline(always)]
unsafe fn x_hit_lanes<S: XInitSink>(
    diag_ptr: *mut DiagStruct,
    cell: usize,
    query_offset: u32,
    subject_offset: u32,
    h: &XHitCtx<'_>,
    b: &mut XLaneBatch,
    init_hsps: &mut S,
) {
    let diag_offset = h.diag_offset;
    let window = h.window;
    let wordsize = h.wordsize;
    // SAFETY: the caller's bound. A raw pointer, not a reference: `resolve_early` writes the cell.
    let entry = unsafe { diag_ptr.add(cell) };
    loop {
        // SAFETY: as above, for every read and write of `entry` below.
        let raw = unsafe { (*entry).raw_bits() };
        // [C] if (diag_array[diag_coord].flag)
        if raw & FLAG_BIT != 0 {
            let last_hit = raw & LAST_HIT_BITS;
            if last_hit < LANES_MAX as u32 {
                // The pending mark of slot `last_hit` (module doc).
                // SAFETY: the caller's bounds.
                unsafe { resolve_early(b, last_hit as usize, diag_ptr, h) };
                continue;
            }
            // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
            let subject_plus_offset = subject_offset.wrapping_add(diag_offset as u32);
            if subject_plus_offset < unsafe { (*entry).last_hit() } as u32 {
                return;
            }
            // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
            // [C] diag_array[diag_coord].flag = 0;
            unsafe {
                (*entry).set_last_hit(subject_plus_offset as i32);
                (*entry).set_flag(0);
            }
            return;
        }
        // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
        let last_hit = unsafe { (*entry).last_hit() } - diag_offset;
        // [C] diff = subject_offset - last_hit;
        let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;
        // [C] if (diff >= window)
        if diff >= window {
            unsafe {
                (*entry).set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32)
            };
            return;
        }
        // [C] if (diff < wordsize)
        if diff < wordsize {
            return;
        }
        // [C] curr_context = BSearchContextInfo(query_offset, query_info);
        let ctx_idx = h.lookup.get_context_idx(query_offset as i32);
        // SAFETY: get_context_idx returns an index below num_contexts == contexts.len().
        let ctx = unsafe { h.contexts.get_unchecked(ctx_idx) };
        // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
        let q_minus_diff = query_offset.wrapping_sub(diff as u32);
        if q_minus_diff < ctx.frame_base as u32 {
            unsafe {
                (*entry).set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32)
            };
            return;
        }
        let q_raw = query_offset.wrapping_sub(ctx.frame_base as u32);
        // [C] cutoffs = word_params->cutoffs + curr_context;
        // SAFETY: ctx_idx < contexts.len() == x_dropoff_per_context.len().
        let x_drop = unsafe { *h.x_dropoff_per_context.get_unchecked(ctx_idx) };
        // The extension's first reads: the word at the second hit and the residues before it in
        // the query (the left walk), the word in the subject. A prefetch never faults and has no
        // architectural effect (`wrapping_add` / `wrapping_sub`: no bound needed).
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use core::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            let q_at = ctx.aa_seq.as_ptr().wrapping_add(1 + q_raw as usize);
            _mm_prefetch(q_at as *const i8, _MM_HINT_T0);
            _mm_prefetch(q_at.wrapping_sub(32) as *const i8, _MM_HINT_T0);
            _mm_prefetch(
                h.subject.as_ptr().wrapping_add(subject_offset as usize) as *const i8,
                _MM_HINT_T0,
            );
        }
        let slot = b.n;
        debug_assert!(slot < b.len && b.len <= LANES_MAX);
        // SAFETY: slot < b.len <= LANES_MAX (a full batch is flushed below).
        unsafe {
            *b.c.get_unchecked_mut(slot) = XCand {
                query_offset,
                s_right: subject_offset,
                s_left: (last_hit + wordsize) as u32,
                q_raw,
                cell: cell as u32,
                ctx_idx: ctx_idx as u32,
                x_drop,
                ..XCand::default()
            };
            (*entry).set_raw_bits(FLAG_BIT | slot as u32);
        }
        b.n = slot + 1;
        b.cands += 1;
        if b.n == b.len {
            // SAFETY: the caller's bounds.
            unsafe { flush(b, diag_ptr, h, init_hsps) };
        }
        return;
    }
}

/// A probe that records nothing, for the shadow's run of `scan_unit`: its own type, so that the
/// shadow adds no caller to the instance of `scan_unit` the units run (`XNoProbe`).
#[derive(Default)]
struct XShadowRef;

impl XProbe for XShadowRef {
    type Mark = ();
    fn unit_begin(&mut self, _width: usize) {}
    fn batch(&mut self, _hits: usize) {}
    fn at(&mut self, _cell: usize, _s: u32) {}
    fn hit(&mut self, _outcome: XHit, _window: i32) {}
    fn batch_end(&mut self) {}
    fn ext_begin(&mut self) {}
    fn ext_end(
        &mut self,
        _mark: (),
        _q: &[u8],
        _s: &[u8],
        _s_left: usize,
        _s_right: usize,
        _q_right: usize,
        _x_drop: i32,
    ) {
    }
    fn unit_end(&mut self) {}
}

/// The lists the shadow compares (`Vec<InitHSP>`, `Vec<XUnitHsp>`).
trait XShadowList: XInitSink + Default + PartialEq + std::fmt::Debug {
    fn len_of(&self) -> usize;
    fn first_difference(&self, other: &Self) -> usize;
}

impl<T: PartialEq + std::fmt::Debug> XShadowList for Vec<T>
where
    Vec<T>: XInitSink,
{
    fn len_of(&self) -> usize {
        self.len()
    }
    fn first_difference(&self, other: &Self) -> usize {
        self.iter()
            .zip(other)
            .position(|(a, b)| a != b)
            .unwrap_or(self.len().min(other.len()))
    }
}

impl XPairPar {
    // EXPERIMENT (LOSAT_X_TBXEXTLANES): `run_unit` with the batched extensions (module doc);
    // LOSAT_X_TBXEXTLANESSHADOW runs both scans on every unit and compares.
    /// The saved init HSPs of one unit, in its scan order.
    #[inline(never)]
    pub(super) fn run_unit_lanes<const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<InitHSP> {
        let mut init_hsps: Vec<InitHSP> = Vec::new();
        self.scan_unit_lanes_or_shadow::<_, PF>(unit, inp, scratch, &mut init_hsps);
        init_hsps
    }

    // EXPERIMENT (LOSAT_X_TBXEXTLANES): `run_unit_compact` with the batched extensions;
    // LOSAT_X_SEEDUNMAP acts as in `run_unit_compact` (`x_unit_records`, `x_unit_copy`).
    /// The saved init HSPs of one unit as compact records, in its scan order.
    #[inline(never)]
    pub(super) fn run_unit_compact_lanes<const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<XUnitHsp> {
        let mut records = std::mem::take(&mut scratch.records);
        records.clear();
        let mapped = self.x_unit_records(unit, scratch, &mut records);
        self.scan_unit_lanes_or_shadow::<_, PF>(unit, inp, scratch, &mut records);
        let exact = x_unit_copy(&records, mapped);
        scratch.records = records;
        exact
    }

    // EXPERIMENT (LOSAT_X_TBXEXTLANESSHADOW): `scan_unit` (on a fresh table, into its own list)
    // and then this module's scan, compared.
    /// The module's scan of one unit, compared with `scan_unit` in shadow mode.
    #[inline(always)]
    fn scan_unit_lanes_or_shadow<S: XShadowList, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        if mode() != 2 {
            self.scan_unit_lanes::<S, PF>(unit, inp, scratch, init_hsps, batch());
            return;
        }
        let mut reference = S::default();
        self.scan_unit::<S, XShadowRef, PF>(unit, inp, scratch, &mut reference, &mut XShadowRef);
        let reference_diag = scratch.diag.clone();
        self.scan_unit_lanes::<S, PF>(unit, inp, scratch, init_hsps, batch());
        if *init_hsps != reference || scratch.diag != reference_diag {
            let first = init_hsps.first_difference(&reference);
            let cell = scratch
                .diag
                .iter()
                .zip(&reference_diag)
                .position(|(a, b)| a != b);
            panic!(
                "LOSAT_X_TBXEXTLANESSHADOW: subject {} frame {} chunk {} group {}: records {} vs {} (first difference at {}), first table difference at {:?}",
                inp.s_idx,
                unit.f,
                unit.k,
                unit.g,
                init_hsps.len_of(),
                reference.len_of(),
                first,
                cell
            );
        }
        EL_SHADOW_UNITS.fetch_add(1, AtomicOrdering::Relaxed);
        EL_SHADOW_RECORDS.fetch_add(init_hsps.len_of() as u64, AtomicOrdering::Relaxed);
    }

    // The scan loop of `scan_unit` (RI:3082-3168 on one unit) with `x_hit_lanes` for `x_hit` and
    // the batch resolved at the end of the unit (module doc). LOSAT_X_TBXPREFETCH's hint (`PF`) as
    // in `scan_unit`.
    /// Scans one unit and saves its init HSPs, in scan order, into `init_hsps`, with batches of
    /// `len` candidates (1 to `LANES_MAX`). Returns the candidates and the early resolutions.
    #[inline(always)]
    fn scan_unit_lanes<S: XInitSink, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
        len: usize,
    ) -> (u64, u64) {
        debug_assert!((1..=LANES_MAX).contains(&len));
        let group = &self.groups[unit.g];
        let sub = &group.lookup;
        let wordsize = inp.wordsize;
        let window = inp.window;
        let s_frame = &inp.s_frames_preliminary[unit.f];
        let subject_full = &s_frame.aa_seq;
        let subject_all = &subject_full[1..subject_full.len() - 1];
        let chunk = unit.chunk;
        let chunk_end = chunk.offset.saturating_add(chunk.length);
        let subject = &subject_all[chunk.offset..chunk_end];
        if subject.len() < wordsize as usize {
            return (0, 0);
        }

        let offset_array_size: i32 = OFFSET_ARRAY_SIZE + sub.longest_chain.max(0);
        let offset_pairs_len =
            crate::algorithm::tblastx::x_scan_bf::x_offset_pairs_len(offset_array_size);
        if scratch.offset_pairs.len() < offset_pairs_len {
            scratch
                .offset_pairs
                .resize(offset_pairs_len, OffsetPair::default());
        }
        let width = x_unit_width(group, &chunk);
        let lo = group.q_min.wrapping_sub(chunk.length as u32 - 1);
        if scratch.diag.capacity() < width
            && x_unit_mapped(scratch.unmap, x_unit_span(group), chunk.length)
        {
            scratch.diag = Vec::with_capacity(x_mapped_len::<DiagStruct>(width));
            SU_DIAG.fetch_add(1, AtomicOrdering::Relaxed);
        }
        scratch.diag.clear();
        scratch.diag.resize(width, DiagStruct::default());

        let h = XHitCtx {
            diag_offset: window,
            window,
            wordsize,
            lookup: sub,
            contexts: inp.contexts,
            cutoff_scores: inp.cutoff_scores,
            x_dropoff_per_context: inp.x_dropoff_per_context,
            subject,
            s_f_idx: unit.f,
            s_frame: s_frame.frame,
            s_idx: inp.s_idx,
            s_len: inp.s_len,
        };
        let mut b = XLaneBatch {
            c: [XCand::default(); LANES_MAX],
            n: 0,
            len,
            cands: 0,
            early: 0,
            batches: 0,
        };

        let base_seq_ranges: [(i32, i32); 1] = [(0, chunk.length as i32)];
        let scan_interiors = tblastx_scan_interiors(chunk.length, None);
        for (interior_start, interior_end) in scan_interiors {
            let seq_ranges = clip_tblastx_seq_ranges_for_scan_interior(
                &base_seq_ranges,
                interior_start,
                interior_end,
                wordsize as usize,
                subject.len(),
            );
            if seq_ranges.is_empty() {
                continue;
            }
            let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];
            while scan_range[1] <= scan_range[2] {
                let prev_scan_left = scan_range[1];
                let hits = crate::algorithm::tblastx::x_scan_bf::x_scan_subject(
                    group.scan_bf,
                    sub,
                    subject,
                    &seq_ranges,
                    &mut scratch.offset_pairs[..],
                    offset_array_size,
                    &mut scan_range,
                );
                if hits == 0 && scan_range[1] == prev_scan_left {
                    break;
                }
                let diag_ptr = scratch.diag.as_mut_ptr();
                let offset_pairs_ptr = scratch.offset_pairs.as_ptr();
                for i in 0..hits as usize {
                    #[cfg(target_arch = "x86_64")]
                    if PF && i + X_PREFETCH_AHEAD < hits as usize {
                        // SAFETY: as in `scan_unit`.
                        unsafe {
                            let ahead = &*offset_pairs_ptr.add(i + X_PREFETCH_AHEAD);
                            let cell =
                                ahead.q_off.wrapping_sub(ahead.s_off).wrapping_sub(lo) as usize;
                            core::arch::x86_64::_mm_prefetch(
                                diag_ptr.wrapping_add(cell) as *const i8,
                                core::arch::x86_64::_MM_HINT_T0,
                            );
                        }
                    }
                    // SAFETY: i < hits <= offset_array_size <= offset_pairs.len().
                    let pair = unsafe { &*offset_pairs_ptr.add(i) };
                    let query_offset = pair.q_off;
                    let subject_offset = pair.s_off;
                    let diag_coord =
                        query_offset.wrapping_sub(subject_offset).wrapping_sub(lo) as usize;
                    debug_assert!(diag_coord < width, "unit diagonal outside its table");
                    // SAFETY: diag_coord < width = scratch.diag.len() (as in `scan_unit`), and
                    // every candidate's cell is such a cell of this unit.
                    unsafe {
                        x_hit_lanes(
                            diag_ptr,
                            diag_coord,
                            query_offset,
                            subject_offset,
                            &h,
                            &mut b,
                            init_hsps,
                        )
                    };
                }
            }
        }
        if b.n > 0 {
            // SAFETY: the candidates' cells are cells of this unit's table.
            unsafe { flush(&mut b, scratch.diag.as_mut_ptr(), &h, init_hsps) };
        }
        EL_UNITS.fetch_add(1, AtomicOrdering::Relaxed);
        EL_CANDS.fetch_add(b.cands, AtomicOrdering::Relaxed);
        EL_EARLY.fetch_add(b.early, AtomicOrdering::Relaxed);
        EL_BATCHES.fetch_add(b.batches, AtomicOrdering::Relaxed);
        (b.cands, b.early)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A frame of amino acids (NCBISTDAA codes) with its sentinels.
    fn frame(frame: i8, residues: &[u8]) -> QueryFrame {
        let mut aa_seq = Vec::with_capacity(residues.len() + 2);
        aa_seq.push(0);
        aa_seq.extend_from_slice(residues);
        aa_seq.push(0);
        QueryFrame {
            frame,
            aa_seq,
            aa_seq_nomask: None,
            aa_len: residues.len(),
            orig_len: residues.len() * 3 + 2,
            seg_masks: Vec::new(),
        }
    }

    /// xorshift64*, for reproducible random residues.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// Residues from `alphabet`, with stretches copied from `source` (shared diagonals: long
    /// extensions, right extensions and later hits on their cells).
    fn residues(rng: &mut Rng, len: usize, alphabet: &[u8], source: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(len);
        while out.len() < len {
            if !source.is_empty() && rng.below(3) == 0 {
                let at = rng.below(source.len() as u64) as usize;
                let take = (rng.below(60) as usize + 5).min(source.len() - at);
                out.extend_from_slice(&source[at..at + take]);
            } else {
                out.push(alphabet[rng.below(alphabet.len() as u64) as usize]);
            }
        }
        out.truncate(len);
        out
    }

    // No NCBI counterpart: compares two LOSAT scans of the same units (module doc, items 1-3).
    #[test]
    fn batched_extensions_equal_scan_unit() {
        let alphabet: Vec<u8> = "ACDEKLMNPQRSTVWY"
            .bytes()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect();
        let narrow: Vec<u8> = "AKLE"
            .bytes()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect();
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let mut rng = Rng(0x5eed_a5a5_0123_4567);
        let mut total_cands = 0u64;
        let mut total_early = 0u64;
        for case in 0..12 {
            let abc = if case % 3 == 2 { &narrow } else { &alphabet };
            let q_len = 200 + rng.below(400) as usize;
            let frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
                .iter()
                .map(|&f| frame(f, &residues(&mut rng, q_len, abc, &[])))
                .collect();
            let source: Vec<u8> = frames
                .iter()
                .flat_map(|f| f.aa_seq[1..f.aa_len + 1].to_vec())
                .collect();
            let s_len = 300 + rng.below(500) as usize;
            let s_frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
                .iter()
                .map(|&f| frame(f, &residues(&mut rng, s_len, abc, &source)))
                .collect();
            let (lookup, contexts) = build_ncbi_lookup(&[frames.clone()], 11, &params, true);
            let cutoff_scores: Vec<i32> = (0..contexts.len())
                .map(|i| 10 + (i as i32 % 4) * 3)
                .collect();
            // X-drops inside and outside LOSAT_X_TBXUNGAP's 1..=127 (its fallback).
            let x_dropoff_per_context: Vec<i32> = (0..contexts.len())
                .map(|i| [7, 16, 30, 200][i % 4])
                .collect();
            let pp = XPairPar {
                groups: build_groups(&lookup, &contexts),
                shadow: false,
                stub: false,
            };
            let window = [40, 19, 64][case % 3];
            let inp = XSubject {
                s_frames_preliminary: &s_frames,
                s_frames: &s_frames,
                contexts: &contexts,
                cutoff_scores: &cutoff_scores,
                x_dropoff_per_context: &x_dropoff_per_context,
                window,
                wordsize: 3,
                s_idx: 0,
                s_len: 999,
            };
            let mut scratch = XScratch::new(false);
            for (f, s_frame) in s_frames.iter().enumerate() {
                let chunk = SubjectChunk {
                    offset: 0,
                    length: s_frame.aa_len,
                    overlap: 0,
                };
                for (g, group) in pp.groups.iter().enumerate() {
                    if group.lookup.longest_chain <= 0 {
                        continue;
                    }
                    let unit = XUnit { f, k: 0, chunk, g };
                    let mut reference: Vec<InitHSP> = Vec::new();
                    pp.scan_unit::<_, _, false>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut reference,
                        &mut XNoProbe,
                    );
                    let reference_diag = scratch.diag.clone();
                    for len in [1, 2, 3, 8, 16] {
                        let mut lanes: Vec<InitHSP> = Vec::new();
                        let (cands, early) = if len == 8 {
                            pp.scan_unit_lanes::<_, true>(
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut lanes,
                                len,
                            )
                        } else {
                            pp.scan_unit_lanes::<_, false>(
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut lanes,
                                len,
                            )
                        };
                        assert_eq!(
                            lanes, reference,
                            "case {case} frame {f} group {g} batch {len}: records"
                        );
                        assert!(
                            scratch.diag == reference_diag,
                            "case {case} frame {f} group {g} batch {len}: table"
                        );
                        if len == 16 {
                            total_cands += cands;
                            total_early += early;
                        }
                    }
                    let mut compact: Vec<XUnitHsp> = Vec::new();
                    let mut compact_ref: Vec<XUnitHsp> = Vec::new();
                    pp.scan_unit::<_, _, false>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut compact_ref,
                        &mut XNoProbe,
                    );
                    pp.scan_unit_lanes::<_, false>(&unit, &inp, &mut scratch, &mut compact, 5);
                    assert_eq!(
                        compact, compact_ref,
                        "case {case} frame {f} group {g}: compact records"
                    );
                }
            }
        }
        // The cases reach both resolution paths: batches and hits on pending cells.
        assert!(total_cands > 1000, "candidates {total_cands}");
        assert!(
            total_early > 10,
            "early resolutions {total_early} of {total_cands}"
        );
    }
}
