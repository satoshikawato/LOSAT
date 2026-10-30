//! EXPERIMENT (LOSAT_X_TBXEXTBATCH): the two-hit stage of a unit of the
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
//! `x_hit` runs an extension as soon as a hit passes the two-hit test, and the extension starts
//! with loads at the hit's query position, a random position of the query group. When the group's
//! frames do not fit the core's caches, those loads wait on memory. Here a hit that passes the
//! test becomes a candidate: its inputs are recorded in a slot of a batch, the cache lines its
//! extension reads first are prefetched, and its cell gets a pending mark. The batch is resolved,
//! by `extend_hit_two_hit` (the function `x_hit` calls), when it holds `BATCH` candidates and at
//! the end of the unit; until then the scan goes on, and the prefetched lines arrive meanwhile.
//!
//! The pending mark is `flag = 1` with the slot number (`0..BATCH`) in `last_hit`. A real cell
//! with `flag = 1` holds `s_last_off - (wordsize - 1) + diag_offset`, and the unit's diag_offset
//! is the window, so with `s_last_off >= 0` its `last_hit` is at least `window - wordsize + 1`;
//! the module runs only where `window >= wordsize + BATCH` (`active`), so a real flag-1 cell
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
//! The rule (`active`): the switch is a value switch, the smallest largest-unit query-group span,
//! in Ki letters, of a subject whose units take this path (built-in `X_EXTBATCH_MIN_KI`); the
//! other subjects run `scan_unit`, whose machine code does not change. x86-64 only (the prefetch;
//! elsewhere the path would only add work). The probes of LOSAT_TIMING / LOSAT_X_STATS run
//! `scan_unit` (this path has none).

use super::*;

/// The candidates of a batch; a pending mark keeps the slot number, below this, in the cell's
/// `last_hit` (module doc). 16: b8 and b16 were alike on the large query, b16 a little ahead on a
/// small one (v4-A DESIGN_A5 §4.1).
pub(super) const BATCH: usize = 16;
/// The packed cell's flag bit (`DiagStruct`: last_hit in bits 0-30, flag in bit 31).
const FLAG_BIT: u32 = 0x8000_0000;
/// The packed cell's `last_hit` bits.
const LAST_HIT_BITS: u32 = 0x7fff_ffff;
/// The switch's built-in: subjects whose largest unit query-group span has at least 256 Ki letters
/// (262,144). Measured (workstation, 8 threads, Sakai prefixes against MG1655's first 1 Mb): frames
/// of 133 k letters x1.005, 267 k x0.985, 533 k x0.977, 1.07 M x0.956, 1.83 M x0.918; the small
/// genomes below (67-96 k letters) +0.4 to +1.5 % with the path on every unit.
const X_EXTBATCH_MIN_KI: usize = 256;
/// The largest number the switch takes (4 Gi letters); a value outside 2..=X_EXTBATCH_MAX_KI other
/// than `0`, or a non-numeric one, gives the built-in (`xswitch::value`).
const X_EXTBATCH_MAX_KI: usize = 1 << 22;

// No NCBI counterpart: which scan the units of the unit path use (scheduling only). Read once per
// process, as the other modes of x_pair_par.rs.
/// LOSAT_X_TBXEXTBATCH: a strict-set value switch (`utils/xswitch.rs`). `None` when it is off
/// (`0`) and on targets other than x86-64; otherwise the
/// smallest largest-unit query-group span, in Ki letters, of a subject whose units take the path.
fn min_span_ki() -> Option<usize> {
    static MIN: OnceLock<Option<usize>> = OnceLock::new();
    *MIN.get_or_init(|| {
        if cfg!(target_arch = "x86_64") {
            crate::utils::xswitch::value(
                "LOSAT_X_TBXEXTBATCH",
                2..=X_EXTBATCH_MAX_KI,
                X_EXTBATCH_MIN_KI,
            )
        } else {
            None
        }
    })
}

// No NCBI counterpart: the rule of the module for one subject (module doc).
/// Whether the units of a subject take the path: the two-hit `window` leaves room for the pending
/// marks, and the largest of the units' query-group `spans` reaches the switch's value.
pub(super) fn active(window: i32, wordsize: i32, spans: impl Iterator<Item = usize>) -> bool {
    let active =
        window as i64 >= wordsize as i64 + BATCH as i64 && x_extbatch_units(min_span_ki(), spans);
    if min_span_ki().is_some() {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXEXTBATCH", active);
    }
    active
}

// No NCBI counterpart: the size rule (module doc), apart from the switch so that a test can give it
// values.
/// Whether a subject whose units have the query-group `spans` takes the path at the switch's value
/// `min_ki`: one of them has at least `min_ki` Ki letters.
fn x_extbatch_units(min_ki: Option<usize>, mut spans: impl Iterator<Item = usize>) -> bool {
    min_ki.is_some_and(|ki| spans.any(|span| span / 1024 >= ki))
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
pub(super) struct XPendBatch {
    c: [XCand; BATCH],
    /// Slots filled.
    n: usize,
    /// Slots of a full batch (`BATCH`; tests take others).
    len: usize,
    cands: u64,
    early: u64,
}

// No NCBI counterpart: the batch of held-back extensions (module doc).
impl XPendBatch {
    /// An empty batch of `len` slots (1 to `BATCH`).
    pub(super) fn new(len: usize) -> Self {
        debug_assert!((1..=BATCH).contains(&len));
        Self {
            c: [XCand::default(); BATCH],
            n: 0,
            len,
            cands: 0,
            early: 0,
        }
    }

    /// Candidates recorded and not yet resolved by `flush`.
    pub(super) fn pending(&self) -> usize {
        self.n
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:575-605 (module doc)
// The extension of one candidate and the cell write of `x_hit` (the statements after its call of
// `extend_hit_two_hit`), the flag written with the value, since the cell holds the pending mark.
/// Runs the extension of `c` and writes its cell as `x_hit` does (`SENT`: LOSAT_X_TBXUNGAPSENT, as
/// in `x_hit`).
///
/// SAFETY: `c.cell` is below the length of the table `diag_ptr` points into; with `SENT`, the unit
/// is one `x_unit_sent` accepts.
#[inline(always)]
unsafe fn resolve<M: TbxScore, const SENT: bool>(
    c: &mut XCand,
    diag_ptr: *mut DiagStruct,
    h: &XHitCtx<'_, M>,
) {
    // SAFETY: the context index was returned by `get_context_idx` (below num_contexts).
    let ctx = unsafe { h.contexts.get_unchecked(c.ctx_idx as usize) };
    let query_full = &ctx.aa_seq;
    // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): as in `x_hit`; another matrix's extension as in `x_hit`.
    let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) = if !M::BLOSUM62
    {
        extend_hit_two_hit_with(
            h.sc,
            &query_full[1..query_full.len() - 1],
            h.subject,
            c.s_left as usize,
            c.s_right as usize,
            c.q_raw as usize,
            c.x_drop,
            false,
        )
    } else if !SENT {
        extend_hit_two_hit(
            &query_full[1..query_full.len() - 1],
            h.subject,
            c.s_left as usize,
            c.s_right as usize,
            c.q_raw as usize,
            c.x_drop,
            false,
        )
    } else {
        // SAFETY: as in `x_hit` (the caller's).
        unsafe {
            extend_hit_two_hit_sent(
                query_full,
                h.subject_sent,
                c.s_left as usize,
                c.s_right as usize,
                c.q_raw as usize,
                c.x_drop,
            )
        }
    };
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
unsafe fn apply<M: TbxScore>(
    c: &mut XCand,
    diag_ptr: *mut DiagStruct,
    h: &XHitCtx<'_, M>,
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
unsafe fn resolve_early<M: TbxScore, const SENT: bool>(
    b: &mut XPendBatch,
    slot: usize,
    diag_ptr: *mut DiagStruct,
    h: &XHitCtx<'_, M>,
) {
    debug_assert!(
        slot < b.n && !b.c[slot].done,
        "a pending mark without its candidate"
    );
    b.early += 1;
    // SAFETY: the caller's bound; the slot holds the candidate that wrote the mark.
    unsafe { resolve::<_, SENT>(b.c.get_unchecked_mut(slot), diag_ptr, h) };
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
pub(super) unsafe fn flush<S: XInitSink, M: TbxScore, const SENT: bool>(
    b: &mut XPendBatch,
    diag_ptr: *mut DiagStruct,
    h: &XHitCtx<'_, M>,
    init_hsps: &mut S,
) {
    for c in &mut b.c[..b.n] {
        if !c.done {
            // SAFETY: the caller's bound.
            unsafe { resolve::<_, SENT>(c, diag_ptr, h) };
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
pub(super) unsafe fn x_hit_batched<S: XInitSink, M: TbxScore, const SENT: bool>(
    diag_ptr: *mut DiagStruct,
    cell: usize,
    query_offset: u32,
    subject_offset: u32,
    h: &XHitCtx<'_, M>,
    b: &mut XPendBatch,
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
            if last_hit < BATCH as u32 {
                // The pending mark of slot `last_hit` (module doc).
                // SAFETY: the caller's bounds.
                unsafe { resolve_early::<_, SENT>(b, last_hit as usize, diag_ptr, h) };
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
        debug_assert!(slot < b.len && b.len <= BATCH);
        // SAFETY: slot < b.len <= BATCH (a full batch is flushed below).
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
            unsafe { flush::<S, _, SENT>(b, diag_ptr, h, init_hsps) };
        }
        return;
    }
}

impl<M: TbxScore> XPairPar<M> {
    // EXPERIMENT (LOSAT_X_TBXEXTBATCH): `run_unit` with the batched extensions (module doc).
    /// The saved init HSPs of one unit, in its scan order.
    #[inline(never)]
    pub(super) fn run_unit_batched<const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<InitHSP> {
        let mut init_hsps: Vec<InitHSP> = Vec::new();
        // EXPERIMENT (LOSAT_X_TBXHITBF): the batched scan or its branch-free hit loop
        // (x_tbx_hitbf.rs).
        self.scan_unit_batched_or_bf::<_, PF>(unit, inp, scratch, &mut init_hsps);
        init_hsps
    }

    // EXPERIMENT (LOSAT_X_TBXEXTBATCH): `run_unit_compact` with the batched extensions;
    // LOSAT_X_SEEDUNMAP acts as in `run_unit_compact` (`x_unit_records`, `x_unit_copy`).
    /// The saved init HSPs of one unit as compact records, in its scan order.
    #[inline(never)]
    pub(super) fn run_unit_compact_batched<const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<XUnitHsp> {
        let mut records = std::mem::take(&mut scratch.records);
        records.clear();
        let mapped = self.x_unit_records(unit, scratch, &mut records);
        self.scan_unit_batched_or_bf::<_, PF>(unit, inp, scratch, &mut records);
        let exact = x_unit_copy(&records, mapped);
        scratch.records = records;
        exact
    }

    // The scan loop of `scan_unit` (RI:3082-3168 on one unit) with `x_hit_batched` for `x_hit` and
    // the batch resolved at the end of the unit (module doc). LOSAT_X_TBXPREFETCH's hint (`PF`) as
    // in `scan_unit`.
    /// Scans one unit and saves its init HSPs, in scan order, into `init_hsps`, with batches of
    /// `len` candidates (1 to `BATCH`). Returns the candidates and the early resolutions.
    #[inline(always)]
    pub(super) fn scan_unit_batched<S: XInitSink, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
        len: usize,
    ) -> (u64, u64) {
        self.scan_unit_batched_hits::<S, PF, false, false>(unit, inp, scratch, init_hsps, len)
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF, x_tbx_hitbf.rs): with `BF` the hits of a batch go through
    // `x_hits_bf`, which takes the too-far and overlap outcomes without a branch and hands the
    // others (a flag, which includes a pending mark, or a candidate) to `x_hit_batched`; without
    // it this is the loop before the switch. `BF` runs only where `x_tbx_hitbf::fits` holds.
    /// `scan_unit_batched`, with `BF` the branch-free hit loop of LOSAT_X_TBXHITBF and with `SENT`
    /// the extensions of `extend_hit_two_hit_sent` (LOSAT_X_TBXUNGAPSENT).
    #[inline(always)]
    pub(super) fn scan_unit_batched_hits<
        S: XInitSink,
        const PF: bool,
        const BF: bool,
        const SENT: bool,
    >(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
        len: usize,
    ) -> (u64, u64) {
        // EXPERIMENT (LOSAT_X_TBXUNITMASK): the masked table where the rule holds (x_pair_par.rs
        // module doc), as in `scan_unit_hits`.
        match self.x_unit_mask(unit, inp.window) {
            None => self.scan_unit_batched_cells::<S, XUnmasked, PF, BF, SENT>(
                unit, inp, scratch, init_hsps, len, 0,
            ),
            Some(mask) => self.scan_unit_batched_masked::<S, PF, BF, SENT>(
                unit, inp, scratch, init_hsps, len, mask,
            ),
        }
    }

    // EXPERIMENT (LOSAT_X_TBXUNITMASK): the masked form in a function of its own, as
    // `scan_unit_masked`.
    /// `scan_unit_batched_cells` with `XMasked`.
    #[inline(never)]
    fn scan_unit_batched_masked<S: XInitSink, const PF: bool, const BF: bool, const SENT: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
        len: usize,
        mask: u32,
    ) -> (u64, u64) {
        self.scan_unit_batched_cells::<S, XMasked, PF, BF, SENT>(
            unit, inp, scratch, init_hsps, len, mask,
        )
    }

    // The statements of `scan_unit_batched_hits` for one form of the unit's cells (`XCells`), as
    // `scan_unit_cells` is for `scan_unit_hits`.
    /// `scan_unit_batched_hits` with the cells `C` (`mask`: the masked table's cells minus one).
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    fn scan_unit_batched_cells<
        S: XInitSink,
        C: XCells,
        const PF: bool,
        const BF: bool,
        const SENT: bool,
    >(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
        len: usize,
        mask: u32,
    ) -> (u64, u64) {
        debug_assert!((1..=BATCH).contains(&len));
        // LOSAT_X_TBXUNGAPSENT: as in `scan_unit_cells`.
        assert!(!SENT || self.x_unit_sent(unit, inp));
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
        // EXPERIMENT (LOSAT_X_TBXHITBF): with `BF` and `PF`, X_PREFETCH_AHEAD entries more, as in
        // `scan_unit_hits`.
        let offset_pairs_len =
            crate::algorithm::tblastx::x_scan_bf::x_offset_pairs_len(offset_array_size)
                + if BF && PF { X_PREFETCH_AHEAD } else { 0 };
        if scratch.offset_pairs.len() < offset_pairs_len {
            scratch
                .offset_pairs
                .resize(offset_pairs_len, OffsetPair::default());
        }
        // LOSAT_X_TBXUNITMASK: a masked table has mask + 1 cells (`scan_unit_cells`).
        let lo = group.q_min.wrapping_sub(chunk.length as u32 - 1);
        let cells = C::of(lo, mask);
        let width = if C::MASKED {
            mask as usize + 1
        } else {
            x_unit_width(group, &chunk)
        };
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
            subject_sent: if SENT { subject_full } else { &[] },
            s_f_idx: unit.f,
            s_frame: s_frame.frame,
            s_idx: inp.s_idx,
            s_len: inp.s_len,
            sc: self.sc,
        };
        let mut b = XPendBatch {
            c: [XCand::default(); BATCH],
            n: 0,
            len,
            cands: 0,
            early: 0,
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
                let hits = group.scan(
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
                // EXPERIMENT (LOSAT_X_TBXHITBF): the loop below, with the too-far and overlap
                // outcomes taken without a branch (x_tbx_hitbf.rs).
                if BF {
                    // SAFETY: as in `scan_unit_hits`; every candidate's cell is a cell of this
                    // unit's table.
                    unsafe {
                        x_tbx_hitbf::x_hits_bf::<C, PF>(
                            diag_ptr,
                            offset_pairs_ptr,
                            hits as usize,
                            cells,
                            width,
                            window,
                            wordsize,
                            |cell, query_offset, subject_offset| {
                                x_hit_batched::<S, _, SENT>(
                                    diag_ptr,
                                    cell,
                                    query_offset,
                                    subject_offset,
                                    &h,
                                    &mut b,
                                    init_hsps,
                                );
                            },
                        )
                    };
                    continue;
                }
                for i in 0..hits as usize {
                    #[cfg(target_arch = "x86_64")]
                    if PF && i + X_PREFETCH_AHEAD < hits as usize {
                        // SAFETY: as in `scan_unit`.
                        unsafe {
                            let ahead = &*offset_pairs_ptr.add(i + X_PREFETCH_AHEAD);
                            let cell = cells.cell(ahead.q_off, ahead.s_off);
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
                    let diag_coord = cells.cell(query_offset, subject_offset);
                    debug_assert!(diag_coord < width, "unit diagonal outside its table");
                    // SAFETY: diag_coord < width = scratch.diag.len() (as in `scan_unit`), and
                    // every candidate's cell is such a cell of this unit.
                    unsafe {
                        x_hit_batched::<S, _, SENT>(
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
            unsafe { flush::<S, _, SENT>(&mut b, scratch.diag.as_mut_ptr(), &h, init_hsps) };
        }
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
            // X-drops inside and outside the 1..=127 of the two-hit table walks (their fallback).
            let x_dropoff_per_context: Vec<i32> = (0..contexts.len())
                .map(|i| [7, 16, 30, 200][i % 4])
                .collect();
            let pp = XPairPar {
                groups: build_groups(&lookup, &contexts, 8),
                stub: false,
                sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
                sc: Blosum62Score,
                sliced: None,
                mask_rule: None,
                extents: Vec::new(),
            };
            // Windows the module takes (`active`): from its smallest, wordsize + BATCH.
            let window = [40, 3 + BATCH as i32, 64][case % 3];
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
                // Every other case scans a chunk that starts inside the frame (offset > 0).
                let offset = if case % 2 == 1 { s_frame.aa_len / 3 } else { 0 };
                let chunk = SubjectChunk {
                    offset,
                    length: s_frame.aa_len - offset,
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
                        let mut batched: Vec<InitHSP> = Vec::new();
                        let (cands, early) = if len == 8 {
                            pp.scan_unit_batched::<_, true>(
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut batched,
                                len,
                            )
                        } else {
                            pp.scan_unit_batched::<_, false>(
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut batched,
                                len,
                            )
                        };
                        assert_eq!(
                            batched, reference,
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
                    pp.scan_unit_batched::<_, false>(&unit, &inp, &mut scratch, &mut compact, 5);
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

    // No NCBI counterpart: the size rule of the switch (module doc).
    #[test]
    fn the_rule_takes_the_largest_span() {
        let spans = [100 * 1024, 300 * 1024 - 1, 50];
        assert!(!x_extbatch_units(None, spans.into_iter()));
        assert!(x_extbatch_units(Some(256), spans.into_iter()));
        assert!(x_extbatch_units(Some(299), spans.into_iter()));
        assert!(!x_extbatch_units(Some(300), spans.into_iter()));
        assert!(x_extbatch_units(Some(2), [2048].into_iter()));
        assert!(!x_extbatch_units(Some(2), [2047].into_iter()));
        assert!(!x_extbatch_units(Some(2), std::iter::empty()));
    }

    // No NCBI counterpart: the batch resolution against `extend_hit_two_hit` on every candidate,
    // with the cell writes and the init HSPs of `x_hit`, some candidates resolved early.
    #[test]
    fn batch_resolution_equals_extend_hit_two_hit() {
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let mut rng = Rng(0x0dd_ba11_cafe_f00d);
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
                subject_sent: &[],
                s_f_idx: 2,
                s_frame: -1,
                s_idx: 7,
                s_len: 999,
                sc: Blosum62Score,
            };
            let n = 1 + rng.below(BATCH as u64) as usize;
            let mut b = XPendBatch {
                c: [XCand::default(); BATCH],
                n,
                len: BATCH,
                cands: 0,
                early: 0,
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
                    unsafe { resolve::<_, false>(&mut b.c[i], diag.as_mut_ptr(), &h) };
                }
            }
            let mut got: Vec<InitHSP> = Vec::new();
            // SAFETY: every candidate's cell is below n = diag.len().
            unsafe { flush::<_, _, false>(&mut b, diag.as_mut_ptr(), &h, &mut got) };
            assert_eq!(diag, expect_cells, "round {round}: cells");
            assert_eq!(got, expect_hsps, "round {round}: init HSPs");
            assert_eq!(b.n, 0);
            cands += n as u64;
        }
        assert!(cands > 1000, "candidates {cands}");
    }
}
