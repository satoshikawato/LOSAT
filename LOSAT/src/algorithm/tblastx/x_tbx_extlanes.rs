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
static EL_SHADOW_CANDS: AtomicU64 = AtomicU64::new(0);
static EL_LANE_CANDS: AtomicU64 = AtomicU64::new(0);
static EL_LANE_BLOCKS: AtomicU64 = AtomicU64::new(0);
static EL_LANE_STEPS: AtomicU64 = AtomicU64::new(0);

/// The most candidates a batch holds; a pending mark keeps the slot number, below this, in the
/// cell's `last_hit` (module doc).
pub(super) const LANES_MAX: usize = 32;
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

/// The table of LOSAT_X_TBXUNGAP's walks when the batches run the lane kernel (`lanes`): x86-64
/// with AVX2, the table built, and LOSAT_X_TBXEXTLANES_SCALAR not set (a prototype choice that
/// keeps step 1, the scalar resolution, for comparison). `None`: every candidate resolves with
/// `extend_hit_two_hit`.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
fn lane_table() -> Option<&'static [i8; 1 << 16]> {
    static T: OnceLock<Option<&'static [i8; 1 << 16]>> = OnceLock::new();
    *T.get_or_init(|| {
        #[cfg(target_arch = "x86_64")]
        if std::env::var_os("LOSAT_X_TBXEXTLANES_SCALAR").is_none()
            && std::arch::is_x86_feature_detected!("avx2")
        {
            return crate::algorithm::tblastx::extension::x_tbxungap_table();
        }
        None
    })
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
    if mode() == 2 || (mode() != 0 && (std::env::var_os("LOSAT_X_STATS").is_some() || timing())) {
        eprintln!(
            "[X_STATS] TBXEXTLANES units={} candidates={} early={} batches={} batch={} kernel={} lane_candidates={} lane_blocks={} lane_steps={} shadow_candidates={}",
            EL_UNITS.load(AtomicOrdering::Relaxed),
            EL_CANDS.load(AtomicOrdering::Relaxed),
            EL_EARLY.load(AtomicOrdering::Relaxed),
            EL_BATCHES.load(AtomicOrdering::Relaxed),
            batch(),
            if lane_table().is_some() { "avx2" } else { "scalar" },
            EL_LANE_CANDS.load(AtomicOrdering::Relaxed),
            EL_LANE_BLOCKS.load(AtomicOrdering::Relaxed),
            EL_LANE_STEPS.load(AtomicOrdering::Relaxed),
            EL_SHADOW_CANDS.load(AtomicOrdering::Relaxed)
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
    /// The lane kernel's candidates, blocks and lane steps (`lanes`).
    lane_cands: u64,
    lane_blocks: u64,
    lane_steps: u64,
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
    // SAFETY: the caller's bound.
    unsafe {
        apply(
            c,
            diag_ptr,
            h,
            (
                hsp_q_u,
                hsp_qe_u,
                hsp_s_u,
                score,
                right_extend,
                s_last_off_u,
            ),
        )
    };
}

/// The result of `c`'s extension, (q_start, q_end, s_start, score, right_extend, s_last_off) of
/// `extend_hit_two_hit`, into `c` and its cell, as `x_hit` writes them.
///
/// SAFETY: as `resolve`.
#[inline(always)]
unsafe fn apply(
    c: &mut XCand,
    diag_ptr: *mut DiagStruct,
    h: &XHitCtx<'_>,
    (hsp_q_u, hsp_qe_u, hsp_s_u, score, right_extend, s_last_off_u): (
        usize,
        usize,
        usize,
        i32,
        bool,
        usize,
    ),
) {
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
    #[cfg(target_arch = "x86_64")]
    if let Some(tab) = lane_table() {
        // SAFETY: AVX2 is present (`lane_table`); the caller's bound.
        unsafe { lanes::resolve_batch(b, diag_ptr, h, tab) };
    }
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
            lane_cands: 0,
            lane_blocks: 0,
            lane_steps: 0,
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
        EL_LANE_CANDS.fetch_add(b.lane_cands, AtomicOrdering::Relaxed);
        EL_LANE_BLOCKS.fetch_add(b.lane_blocks, AtomicOrdering::Relaxed);
        EL_LANE_STEPS.fetch_add(b.lane_steps, AtomicOrdering::Relaxed);
        (b.cands, b.early)
    }
}

// EXPERIMENT (LOSAT_X_TBXEXTLANES, step 2): the lane kernel of a batch (x86-64, AVX2).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:846-866, 886-921 (the walks;
// module comment of tblastx/extension/x_two_hit_fast.rs)
// ```c
// for (i = 0; i < n; i++) {
//     score += matrix[q[i]][s[i]];
//     if (score > maxscore) { maxscore = score; best_i = i; }
//     if (score <= 0 || (maxscore - score) >= dropoff) break;
// }
// ```
// Eight walks run side by side, one per 32-bit lane: the left walk of a candidate, then its
// right walk in the same lane. A block writes the scores of the next K steps of every lane into
// a step-major byte buffer (one table load per lane and step, `T[q << 8 | s]` of
// LOSAT_X_TBXUNGAP), then runs the recurrence of the walks K steps on all lanes: per lane the
// reference's statements in its order (add, improve, stop test), on i32 as the reference. A lane
// that stopped keeps its values (its `active` lane is cleared). At the end of a block the lanes
// whose walk stopped hand its result on (a left walk to the right walk, a right walk to the
// candidate's result), and the free lanes take the next candidates of the batch.
//
// The control avoids branches on one lane's data (the walks of `two_hit` pay a misprediction at
// each exit): the lanes that stopped are taken from a bit mask, every left walk is followed by a
// right walk (of no steps when the left walk did not reach the first hit), and a walk with no
// steps left, like a free lane, reads from a block of code-0 residues (`ZEROS`), whose scores are
// the sentinel.
//
// Why the result is that of `two_hit` (x_two_hit_fast.rs), which is that of the reference.
// * The word loop is `two_hit`'s, in scalar, before the walks.
// * A walk reads the positions of `two_hit`'s walk, in its order: the left walk `min(q_r, s_r)`
//   steps back, the right walk `min(q_len - q_r, s_len - s_r)` steps forward (none if the left walk
//   did not reach the first hit); every step after the last one reads the sentinel -128 (the
//   buffer's fill, or `ZEROS`: code 0 against anything is `BLAST_SCORE_MIN`).
// * The sentinel: x_drop is 1..=127 (else the candidate is left to `extend_hit_two_hit`), and the
//   score before a step is never above the best, so a step of -128 stops the walk without raising
//   the best (two_hit's argument). Stopping one step after the last real step changes nothing on
//   the left (only the best and its step are returned) and gives the right walk's stop position
//   `n` (taken - 1 = n), the value of `two_hit`'s walk that runs out. A right walk of no steps
//   stops at its first step with its start (the left walk's best) as the best, length 0 and stop
//   position 0, which are the values `two_hit` uses when it does not walk right: score =
//   max(left, 0) = left, q_end = q_r, s_last_off = s_r.
// * Counters: `t` counts the steps of the walk (1-based); the best's step is `t` at its last strict
//   improvement (0 if none): the left walk's length and the right walk's `best_len`; the right
//   walk's stop position is the number of steps before the stopping one (`alive`).
#[cfg(target_arch = "x86_64")]
mod lanes {
    use super::*;
    use crate::algorithm::tblastx::extension::X_TBXUNGAP_MAX_X_DROP;
    use core::arch::x86_64::*;

    /// Lanes of a 256-bit vector of i32.
    const NL: usize = 8;
    /// Steps of a block.
    const K: usize = 16;
    /// `BLAST_SCORE_MIN` in the table (x_two_hit_fast.rs), written where a lane has no step.
    const SENTINEL: i8 = -128;
    /// Residues of the word loop (`k_size` of two_hit.rs).
    const WORD: usize = 3;

    /// Code-0 residues (the sentinel against anything): where a walk with no steps left reads,
    /// `K` bytes on either side of `ZEROS.0[K]`.
    #[repr(C, align(32))]
    struct Zeros([u8; 2 * K]);
    static ZEROS: Zeros = Zeros([0; 2 * K]);

    /// `pshufb` masks: the 16 bytes in order (a right walk) and reversed (a left walk, whose step
    /// `i` reads the byte `K - 1 - i` of the 16 before its pointer).
    #[repr(C, align(16))]
    struct Masks([[u8; K]; 2]);
    static MASKS: Masks = Masks([
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
    ]);

    type Table = [i8; 1 << 16];

    /// A candidate the lanes take, after the word loop.
    #[derive(Clone, Copy)]
    struct Prep {
        slot: usize,
        /// The end of the best word at R (`two_hit`'s `q_r`, `s_r`).
        q_r: usize,
        s_r: usize,
        /// `q_seq + q_r`, `s_seq + s_r`.
        q_at: *const u8,
        s_at: *const u8,
        /// The query of its context, for the shadow.
        q_base: *const u8,
        q_len: usize,
        /// The left walk's steps, the right walk's (if it runs), `s_r - s_left` (saturating).
        n_left: usize,
        n_right: usize,
        reach: usize,
        xd1: i32,
    }

    const PREP_NONE: Prep = Prep {
        slot: 0,
        q_r: 0,
        s_r: 0,
        q_at: core::ptr::null(),
        s_at: core::ptr::null(),
        q_base: core::ptr::null(),
        q_len: 0,
        n_left: 0,
        n_right: 0,
        reach: 0,
        xd1: 0,
    };

    /// The scalar state of the eight lanes; lane `l` runs the prepared candidate `cur[l]`.
    struct Lanes {
        cur: [usize; NL],
        /// Where the next step reads (`q[-1]`, `s[-1]` on the left, `q[0]`, `s[0]` on the right).
        q: [*const u8; NL],
        s: [*const u8; NL],
        /// Steps left.
        r: [usize; NL],
        /// The left walk's best and length, and whether it reached the first hit.
        left_score: [i32; NL],
        left_disp: [usize; NL],
        reached: [bool; NL],
    }

    /// The i32 values of the eight lanes, as stored between blocks.
    #[repr(C, align(32))]
    #[derive(Clone, Copy)]
    struct V([i32; NL]);

    /// The scores of a block, step-major: step `k` of lane `l` at `k * NL + l`.
    #[repr(C, align(32))]
    struct Buf([i8; K * NL]);

    /// The vector state of the eight walks.
    struct Walks {
        score: V,
        best: V,
        /// `t` at the best's last strict improvement.
        pos: V,
        /// Steps taken in the walk.
        t: V,
        /// x_drop - 1.
        xd1: V,
        /// 1 for a right walk (its `score <= 0` stop), `i32::MIN` for a left walk (none).
        floor: V,
        /// The steps of the last block after which the lane was still going.
        alive: V,
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:846-866 (module comment
    // above): `K` steps of the eight walks.
    /// Runs one block on the lanes of `busy`; returns those that stopped in it.
    #[target_feature(enable = "avx2")]
    unsafe fn run_block(buf: &Buf, w: &mut Walks, busy: u32) -> u32 {
        let ld = |v: &V| unsafe { _mm256_load_si256(v.0.as_ptr() as *const __m256i) };
        let mut sc = ld(&w.score);
        let mut bs = ld(&w.best);
        let mut ps = ld(&w.pos);
        let mut tt = ld(&w.t);
        let xd = ld(&w.xd1);
        let fl = ld(&w.floor);
        let bits = _mm256_setr_epi32(1, 2, 4, 8, 16, 32, 64, 128);
        let act = _mm256_cmpeq_epi32(_mm256_and_si256(_mm256_set1_epi32(busy as i32), bits), bits);
        let mut ac = act;
        let ones = _mm256_set1_epi32(1);
        let mut alive = _mm256_setzero_si256();
        for k in 0..K {
            // SAFETY: k * NL + 8 <= K * NL bytes of `buf`.
            let v = _mm256_cvtepi8_epi32(unsafe {
                _mm_loadl_epi64(buf.0.as_ptr().add(k * NL) as *const __m128i)
            });
            tt = _mm256_add_epi32(tt, ones);
            // [C] score += matrix[q[i]][s[i]];
            sc = _mm256_add_epi32(sc, v);
            // [C] if (score > maxscore) { maxscore = score; best_i = i; }
            let gt = _mm256_and_si256(_mm256_cmpgt_epi32(sc, bs), ac);
            bs = _mm256_blendv_epi8(bs, sc, gt);
            ps = _mm256_blendv_epi8(ps, tt, gt);
            // [C] if (score <= 0 || (maxscore - score) >= dropoff) break;
            let stop = _mm256_or_si256(
                _mm256_cmpgt_epi32(_mm256_sub_epi32(bs, sc), xd),
                _mm256_cmpgt_epi32(fl, sc),
            );
            ac = _mm256_andnot_si256(stop, ac);
            alive = _mm256_sub_epi32(alive, ac);
        }
        let st = |v: &mut V, x: __m256i| unsafe {
            _mm256_store_si256(v.0.as_mut_ptr() as *mut __m256i, x)
        };
        st(&mut w.score, sc);
        st(&mut w.best, bs);
        st(&mut w.pos, ps);
        st(&mut w.t, tt);
        st(&mut w.alive, alive);
        // Busy at the start and not at the end.
        _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_andnot_si256(ac, act))) as u32
    }

    /// Writes the scores of the next `K` steps of every lane. A lane with `K` steps or more left,
    /// or none (it reads `ZEROS`), takes its indices `q << 8 | s` from two 16-byte loads (reversed
    /// on the left) interleaved; a walk's last steps near the start or end of a sequence are read
    /// one by one, after the sentinel fill.
    ///
    /// SAFETY: every lane can read its `r` steps from its pointers, and a lane with no steps left
    /// points at `ZEROS.0[K]`.
    #[target_feature(enable = "avx2")]
    unsafe fn fill(buf: &mut Buf, ln: &Lanes, left: u32, tab: &Table) {
        buf.0 = [SENTINEL; K * NL];
        let mut idx = [0u16; K];
        for l in 0..NL {
            let back = (left >> l) & 1;
            let (q, s, r) = (ln.q[l], ln.s[l], ln.r[l]);
            if r >= K || r == 0 {
                let off = K * back as usize;
                // SAFETY: K positions before the pointers (left) or from them (right) are the
                // walk's next steps, or `ZEROS`; the mask table has two entries.
                unsafe {
                    let mk = _mm_load_si128(MASKS.0.as_ptr().add(back as usize) as *const __m128i);
                    let qv = _mm_shuffle_epi8(
                        _mm_loadu_si128(q.wrapping_sub(off) as *const __m128i),
                        mk,
                    );
                    let sv = _mm_shuffle_epi8(
                        _mm_loadu_si128(s.wrapping_sub(off) as *const __m128i),
                        mk,
                    );
                    _mm_storeu_si128(idx.as_mut_ptr() as *mut __m128i, _mm_unpacklo_epi8(sv, qv));
                    _mm_storeu_si128(
                        idx.as_mut_ptr().add(8) as *mut __m128i,
                        _mm_unpackhi_epi8(sv, qv),
                    );
                    for (i, &x) in idx.iter().enumerate() {
                        *buf.0.get_unchecked_mut(i * NL + l) = *tab.get_unchecked(x as usize);
                    }
                }
            } else if back != 0 {
                for i in 0..r {
                    // SAFETY: the walk's bound; q << 8 | s < 2^16.
                    unsafe {
                        let x = (*q.sub(1 + i) as usize) << 8 | *s.sub(1 + i) as usize;
                        *buf.0.get_unchecked_mut(i * NL + l) = *tab.get_unchecked(x);
                    }
                }
            } else {
                for i in 0..r {
                    // SAFETY: as above.
                    unsafe {
                        let x = (*q.add(i) as usize) << 8 | *s.add(i) as usize;
                        *buf.0.get_unchecked_mut(i * NL + l) = *tab.get_unchecked(x);
                    }
                }
            }
        }
    }

    /// Starts a walk on lane `l`: `n` steps from `q`, `s` (no steps: `ZEROS`), from `start`.
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    fn set_walk(
        ln: &mut Lanes,
        w: &mut Walks,
        l: usize,
        q: *const u8,
        s: *const u8,
        n: usize,
        start: i32,
        xd1: i32,
        floor: i32,
    ) {
        let zero = ZEROS.0[K..].as_ptr();
        ln.q[l] = if n == 0 { zero } else { q };
        ln.s[l] = if n == 0 { zero } else { s };
        ln.r[l] = n;
        w.score.0[l] = start;
        w.best.0[l] = start;
        w.pos.0[l] = 0;
        w.t.0[l] = 0;
        w.xd1.0[l] = xd1;
        w.floor.0[l] = floor;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158
    // `s_BlastAaExtendTwoHit` (two_hit.rs and x_two_hit_fast.rs follow it): the word loop, the
    // left walk, the test that it reached the first hit, the right walk, the result.
    /// Resolves the batch's candidates that the lanes can take; the others stay for
    /// `extend_hit_two_hit` (`flush`).
    ///
    /// SAFETY: AVX2 is available; every candidate's cell is below the length of the table
    /// `diag_ptr` points into.
    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn resolve_batch(
        b: &mut XLaneBatch,
        diag_ptr: *mut DiagStruct,
        h: &XHitCtx<'_>,
        tab: &Table,
    ) {
        let shadow = mode() == 2;
        let s_seq = h.subject;
        // The candidates the lanes take, in slot order, after the word loop.
        let mut prep = [PREP_NONE; LANES_MAX];
        let mut k = 0usize;
        for slot in 0..b.n {
            let c = b.c[slot];
            // SAFETY: the context index was returned by `get_context_idx` (below num_contexts).
            let ctx = unsafe { h.contexts.get_unchecked(c.ctx_idx as usize) };
            let query_full = &ctx.aa_seq;
            let q = &query_full[1..query_full.len() - 1];
            let (q_raw, s_right) = (c.q_raw as usize, c.s_right as usize);
            if c.done
                || !(1..=X_TBXUNGAP_MAX_X_DROP).contains(&c.x_drop)
                || q_raw >= q.len()
                || q.len() - q_raw < WORD
                || s_right >= s_seq.len()
                || s_seq.len() - s_right < WORD
            {
                continue;
            }
            // The best end of the word at R (`two_hit`).
            let mut score = 0i32;
            let mut seed_score = 0i32;
            let mut right_d = 0usize;
            for i in 0..WORD {
                // SAFETY: the word fits in both sequences (tested above).
                let (qc, sc) = unsafe {
                    (
                        *q.get_unchecked(q_raw + i),
                        *s_seq.get_unchecked(s_right + i),
                    )
                };
                score += tab[(qc as usize) << 8 | sc as usize] as i32;
                if score > seed_score {
                    seed_score = score;
                    right_d = i + 1;
                }
            }
            let (q_r, s_r) = (q_raw + right_d, s_right + right_d);
            prep[k] = Prep {
                slot,
                q_r,
                s_r,
                // SAFETY: q_r <= q.len() and s_r <= s.len() (the word fits).
                q_at: unsafe { q.as_ptr().add(q_r) },
                s_at: unsafe { s_seq.as_ptr().add(s_r) },
                q_base: q.as_ptr(),
                q_len: q.len(),
                n_left: q_r.min(s_r),
                n_right: (q.len() - q_r).min(s_seq.len() - s_r),
                reach: s_r.saturating_sub(c.s_left as usize),
                xd1: c.x_drop - 1,
            };
            k += 1;
        }
        b.lane_cands += k as u64;

        let zero = ZEROS.0[K..].as_ptr();
        let mut ln = Lanes {
            cur: [0; NL],
            q: [zero; NL],
            s: [zero; NL],
            r: [0; NL],
            left_score: [0; NL],
            left_disp: [0; NL],
            reached: [false; NL],
        };
        let mut w = Walks {
            score: V([0; NL]),
            best: V([0; NL]),
            pos: V([0; NL]),
            t: V([0; NL]),
            xd1: V([0; NL]),
            floor: V([0; NL]),
            alive: V([0; NL]),
        };
        let mut buf = Buf([SENTINEL; K * NL]);
        // Lanes with a candidate, and those of them on its left walk.
        let mut busy = 0u32;
        let mut left = 0u32;
        let mut next = 0usize;
        loop {
            let mut free = !busy & 0xff;
            while free != 0 && next < k {
                let l = free.trailing_zeros() as usize;
                free &= free - 1;
                let p = &prep[next];
                ln.cur[l] = next;
                set_walk(
                    &mut ln,
                    &mut w,
                    l,
                    p.q_at,
                    p.s_at,
                    p.n_left,
                    0,
                    p.xd1,
                    i32::MIN,
                );
                busy |= 1 << l;
                left |= 1 << l;
                next += 1;
            }
            if busy == 0 {
                break;
            }
            // SAFETY: the lanes' pointers (`set_walk`, the advance below).
            unsafe { fill(&mut buf, &ln, left, tab) };
            // SAFETY: AVX2 (the caller's bound).
            let stopped = unsafe { run_block(&buf, &mut w, busy) };
            b.lane_blocks += 1;
            for l in 0..NL {
                let r = ln.r[l];
                let m = r.min(K);
                b.lane_steps += m as u64;
                let d = if (left >> l) & 1 != 0 {
                    m.wrapping_neg()
                } else {
                    m
                };
                ln.r[l] = r - m;
                ln.q[l] = if r == m {
                    zero
                } else {
                    ln.q[l].wrapping_add(d)
                };
                ln.s[l] = if r == m {
                    zero
                } else {
                    ln.s[l].wrapping_add(d)
                };
            }
            // Right walks that stopped: the candidate's result (module comment above).
            let mut done = stopped & !left;
            while done != 0 {
                let l = done.trailing_zeros() as usize;
                done &= done - 1;
                let p = &prep[ln.cur[l]];
                let right_best = w.best.0[l];
                let right_disp = w.pos.0[l] as usize;
                let end = w.t.0[l] as usize - K + w.alive.0[l] as usize;
                let disp = ln.left_disp[l];
                let result = (
                    p.q_r - disp,
                    p.q_r + right_disp,
                    p.s_r - disp,
                    ln.left_score[l].max(right_best),
                    ln.reached[l],
                    p.s_r + end,
                );
                // SAFETY: the slot is below b.n.
                let c = unsafe { b.c.get_unchecked_mut(p.slot) };
                if shadow {
                    // SAFETY: q_base points at the candidate's query of q_len residues.
                    let q = unsafe { core::slice::from_raw_parts(p.q_base, p.q_len) };
                    let expect = extend_hit_two_hit(
                        q,
                        s_seq,
                        c.s_left as usize,
                        c.s_right as usize,
                        c.q_raw as usize,
                        c.x_drop,
                        false,
                    );
                    let got = (
                        result.0,
                        result.1,
                        result.2,
                        result.2 + (result.1 - result.0),
                        result.3,
                        result.4,
                        result.5,
                    );
                    assert_eq!(
                        got, expect,
                        "LOSAT_X_TBXEXTLANESSHADOW: lane kernel s_left={} s_right={} q_raw={} x_drop={} q_len={} s_len={}",
                        c.s_left, c.s_right, c.q_raw, c.x_drop, p.q_len, s_seq.len()
                    );
                    EL_SHADOW_CANDS.fetch_add(1, AtomicOrdering::Relaxed);
                }
                // SAFETY: the caller's bound.
                unsafe { apply(c, diag_ptr, h, result) };
                ln.q[l] = zero;
                ln.s[l] = zero;
                ln.r[l] = 0;
                busy &= !(1 << l);
            }
            // Left walks that stopped: the right walk, of no steps if it did not reach.
            // [C] if (left_d >= (s_right_off - s_left_off)) { *right_extend = TRUE; ... }
            let mut done = stopped & left;
            while done != 0 {
                let l = done.trailing_zeros() as usize;
                done &= done - 1;
                let p = &prep[ln.cur[l]];
                let best = w.best.0[l];
                let disp = w.pos.0[l] as usize;
                let reached = disp >= p.reach;
                ln.left_score[l] = best;
                ln.left_disp[l] = disp;
                ln.reached[l] = reached;
                let n = if reached { p.n_right } else { 0 };
                let xd1 = w.xd1.0[l];
                set_walk(&mut ln, &mut w, l, p.q_at, p.s_at, n, best, xd1, 1);
                left &= !(1 << l);
            }
        }
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
        for case in 0..16 {
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
            // Every fourth case scans the query against itself: long walks, walks that reach the
            // ends of both sequences.
            let s_frames: Vec<QueryFrame> = if case % 4 == 3 {
                frames.clone()
            } else {
                [1i8, 2, 3, -1, -2, -3]
                    .iter()
                    .map(|&f| frame(f, &residues(&mut rng, s_len, abc, &source)))
                    .collect()
            };
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
            // Windows the module takes (`active`): from its smallest, wordsize + LANES_MAX.
            let window = [40, 3 + LANES_MAX as i32, 64][case % 3];
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

    // No NCBI counterpart: the batch resolution (the lane kernel where AVX2 is present, and
    // `extend_hit_two_hit` for the rest) against `extend_hit_two_hit` on every candidate, with the
    // cell writes and the init HSPs of `x_hit`.
    #[test]
    fn batch_resolution_equals_extend_hit_two_hit() {
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let mut rng = Rng(0x0dd_ba11_cafe_f00d);
        let mut lane_cands = 0u64;
        let mut cands = 0u64;
        for round in 0..400 {
            let q_len = 3 + rng.below(if round % 5 == 0 { 40 } else { 300 }) as usize;
            // Codes 1..=24, a code 0 (BLAST_SCORE_MIN) now and then.
            let pick = |rng: &mut Rng| {
                if rng.below(97) == 0 {
                    0u8
                } else {
                    1 + rng.below(24) as u8
                }
            };
            let frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
                .iter()
                .map(|&f| frame(f, &(0..q_len).map(|_| pick(&mut rng)).collect::<Vec<u8>>()))
                .collect();
            let (lookup, contexts) = build_ncbi_lookup(&[frames.clone()], 11, &params, true);
            // The subject: pieces of the query frames with substitutions, and random residues.
            let mut subject: Vec<u8> = Vec::new();
            let s_len = 3 + rng.below(400) as usize;
            while subject.len() < s_len {
                if rng.below(2) == 0 {
                    let f = &frames[rng.below(6) as usize].aa_seq[1..q_len + 1];
                    let at = rng.below(q_len as u64) as usize;
                    let take = (1 + rng.below(80) as usize).min(q_len - at);
                    let rate = [3u64, 8, 30, 1_000_000][rng.below(4) as usize];
                    subject.extend(f[at..at + take].iter().map(|&c| {
                        if rng.below(rate) == 0 {
                            pick(&mut rng)
                        } else {
                            c
                        }
                    }));
                } else {
                    subject.push(pick(&mut rng));
                }
            }
            subject.truncate(s_len);
            let cutoff_scores: Vec<i32> = (0..contexts.len()).map(|i| 8 + i as i32).collect();
            let x_dropoff_per_context: Vec<i32> = (0..contexts.len()).map(|_| 0).collect();
            let window = 40;
            let h = XHitCtx {
                diag_offset: window,
                window,
                wordsize: 3,
                lookup: &lookup,
                contexts: &contexts,
                cutoff_scores: &cutoff_scores,
                x_dropoff_per_context: &x_dropoff_per_context,
                subject: &subject,
                s_f_idx: 2,
                s_frame: -1,
                s_idx: 7,
                s_len: 999,
            };
            let n = 1 + rng.below(LANES_MAX as u64) as usize;
            let mut b = XLaneBatch {
                c: [XCand::default(); LANES_MAX],
                n,
                len: LANES_MAX,
                cands: 0,
                early: 0,
                batches: 0,
                lane_cands: 0,
                lane_blocks: 0,
                lane_steps: 0,
            };
            let mut diag = vec![DiagStruct::default(); n];
            let mut expect_cells = Vec::new();
            let mut expect_hsps: Vec<InitHSP> = Vec::new();
            for i in 0..n {
                let ctx_idx = rng.below(contexts.len() as u64) as usize;
                let ctx = &contexts[ctx_idx];
                let q = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
                // A word that fits (a lookup hit), as `two_hit.rs` reads it unchecked.
                let q_raw = rng.below((q.len() - 2) as u64) as usize;
                let s_right = rng.below((subject.len() - 2) as u64) as usize;
                let diff = 3 + rng.below(37) as usize;
                let s_left = (s_right + 3).saturating_sub(diff);
                let x_drop = [1, 2, 7, 16, 30, 127, 128, 300][rng.below(8) as usize];
                b.c[i] = XCand {
                    query_offset: ctx.frame_base as u32 + q_raw as u32,
                    s_right: s_right as u32,
                    s_left: s_left as u32,
                    q_raw: q_raw as u32,
                    cell: i as u32,
                    ctx_idx: ctx_idx as u32,
                    x_drop,
                    ..XCand::default()
                };
                diag[i].set_raw_bits(FLAG_BIT | i as u32);
                let (q0, qe, s0, _se, score, right, s_last) =
                    extend_hit_two_hit(q, &subject, s_left, s_right, q_raw, x_drop, false);
                let mut cell = DiagStruct::default();
                if right {
                    cell.set_flag(1);
                    cell.set_last_hit(s_last as i32 - 2 + window);
                } else {
                    cell.set_last_hit(s_right as i32 + window);
                }
                expect_cells.push(cell);
                if score >= cutoff_scores[ctx_idx] {
                    let (q0, qe, s0) = (q0 as i32, qe as i32, s0 as i32);
                    expect_hsps.push(InitHSP {
                        q_start_absolute: ctx.frame_base + q0,
                        q_end_absolute: ctx.frame_base + qe,
                        s_start: s0,
                        s_end: s0 + (qe - q0),
                        q_seed_absolute: ctx.frame_base + q_raw as i32,
                        s_seed: s_right as i32,
                        score,
                        ctx_idx,
                        s_f_idx: 2,
                        q_idx: ctx.q_idx,
                        s_idx: 7,
                        q_frame: ctx.frame,
                        s_frame: -1,
                        q_orig_len: ctx.orig_len,
                        s_orig_len: 999,
                    });
                }
            }
            // A few candidates resolved already, as a hit on their cell does.
            for i in 0..n {
                if rng.below(10) == 0 {
                    // SAFETY: the cell is below n = diag.len().
                    unsafe { resolve(&mut b.c[i], diag.as_mut_ptr(), &h) };
                }
            }
            let mut got: Vec<InitHSP> = Vec::new();
            // SAFETY: every candidate's cell is below n = diag.len().
            unsafe { flush(&mut b, diag.as_mut_ptr(), &h, &mut got) };
            assert_eq!(diag, expect_cells, "round {round}: cells");
            assert_eq!(got, expect_hsps, "round {round}: init HSPs");
            assert_eq!(b.n, 0);
            lane_cands += b.lane_cands;
            cands += n as u64;
        }
        if lane_table().is_some() {
            assert!(
                lane_cands * 2 > cands,
                "the lane kernel took {lane_cands} of {cands} candidates"
            );
        }
    }
}
