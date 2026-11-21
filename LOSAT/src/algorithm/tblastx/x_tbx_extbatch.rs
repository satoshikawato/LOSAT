//! EXPERIMENT (LOSAT_X_TBXEXTBATCH): the two-hit stage of a unit of the
//! unit path (x_pair_par.rs) with its extensions held back and resolved in batches. A child module
//! of x_pair_par.rs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:575-605
//! `x_hit` runs an extension as soon as a hit passes the two-hit test, and the extension starts
//! with loads at the hit's query position, a random position of the query group. When the group's
//! frames do not fit the core's caches, those loads wait on memory. Here a hit that passes the
//! test becomes a candidate: its inputs are recorded in a slot of a batch, the cache lines its
//! extension reads first are prefetched, and its cell gets a pending mark. The batch is resolved,
//! by `extend_hit_two_hit_fitting` (the function `x_hit` calls), when it holds `BATCH` candidates and at
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

/// One held-back extension: the inputs of `extend_hit_two_hit_fitting` and of the init HSP, then its
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
// `extend_hit_two_hit_fitting`), the flag written with the value, since the cell holds the pending mark.
/// Runs the extension of `c` and writes its cell as `x_hit` does (`SENT`: LOSAT_X_TBXUNGAPSENT, as
/// in `x_hit`).
///
/// SAFETY: `c.cell` is below the length of the table `diag_ptr` points into; with `SENT`, the unit
/// is one `x_unit_sent` accepts.
#[inline(always)]
unsafe fn resolve<M: TbxScore, const SENT: usize>(
    c: &mut XCand,
    diag_ptr: *mut DiagStruct,
    h: &XHitCtx<'_, M>,
) {
    // SAFETY: the context index was returned by `get_context_idx` (below num_contexts).
    let ctx = unsafe { h.contexts.get_unchecked(c.ctx_idx as usize) };
    let query_full = &ctx.aa_seq;
    // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): as in `x_hit`; another matrix's extension as in `x_hit`.
    // The walks go on into the query batch past the context's ends, as in `x_hit`.
    let in_batch = h.query_batch.at(c.ctx_idx as usize, h.s_idx);
    let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) =
        if !M::BLOSUM62 && SENT == 0 {
            extend_hit_two_hit_with(
                h.sc,
                &query_full[1..query_full.len() - 1],
                in_batch,
                h.subject,
                c.s_left as usize,
                c.s_right as usize,
                c.q_raw as usize,
                c.x_drop,
                false,
            )
        } else if SENT == 0 {
            extend_hit_two_hit_fitting(
                &query_full[1..query_full.len() - 1],
                in_batch,
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
                h.sc.sent_two_hit_of::<SENT>(
                    query_full,
                    in_batch,
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
/// `extend_hit_two_hit_fitting`, into `c` and its cell, as `x_hit` writes them.
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
unsafe fn resolve_early<M: TbxScore, const SENT: usize>(
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-591
// The batch, slot by slot: the candidates not resolved yet, then the init HSPs in slot order
// (module doc, item 3), by the push statement of `x_hit`.
/// Resolves every candidate of the batch, saves their init HSPs in slot order and empties it.
/// Its own function, so that the scan loop keeps its values in registers.
///
/// SAFETY: every candidate's cell is below the length of the table `diag_ptr` points into.
#[inline(never)]
pub(super) unsafe fn flush<S: XInitSink, M: TbxScore, const SENT: usize>(
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
pub(super) unsafe fn x_hit_batched<S: XInitSink, M: TbxScore, const SENT: usize>(
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
            // SAFETY: `entry` (cell `cell`) is in the table by this function's contract.
            if subject_plus_offset < unsafe { (*entry).last_hit() } as u32 {
                return;
            }
            // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
            // [C] diag_array[diag_coord].flag = 0;
            // SAFETY: `entry` (cell `cell`) is in the table by this function's contract.
            unsafe {
                (*entry).set_last_hit(subject_plus_offset as i32);
                (*entry).set_flag(0);
            }
            return;
        }
        // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
        // SAFETY: `entry` (cell `cell`) is in the table by this function's contract.
        let last_hit = unsafe { (*entry).last_hit() } - diag_offset;
        // [C] diff = subject_offset - last_hit;
        let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;
        // [C] if (diff >= window)
        if diff >= window {
            // SAFETY: `entry` (cell `cell`) is in the table by this function's contract.
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
        let ctx_idx = h.context_of(query_offset);
        // SAFETY: context_of returns an index below num_contexts == contexts.len().
        let ctx = unsafe { h.contexts.get_unchecked(ctx_idx) };
        // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
        let q_minus_diff = query_offset.wrapping_sub(diff as u32);
        if q_minus_diff < ctx.frame_base as u32 {
            // SAFETY: `entry` (cell `cell`) is in the table by this function's contract.
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
        // SAFETY: `_mm_prefetch` needs only SSE, which every x86_64 target has, and never faults.
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

    // The scan loop of `scan_unit` (the reference's, `scan_chunk_plain` with
    // `word_finder::scan_next_hits`, on one unit) with `x_hit_batched` for `x_hit` and
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
        self.scan_unit_batched_hits::<S, PF, false, 0>(unit, inp, scratch, init_hsps, len)
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF, x_tbx_hitbf.rs): with `BF` the hits of a batch go through
    // `x_hits_bf`, which takes the too-far and overlap outcomes without a branch and hands the
    // others (a flag, which includes a pending mark, or a candidate) to `x_hit_batched`; without
    // it this is the loop before the switch. `BF` runs only where `x_tbx_hitbf::fits` holds.
    /// `scan_unit_batched`, with `BF` the branch-free hit loop of LOSAT_X_TBXHITBF and with `SENT`
    /// (a word length, 0: none) the extensions of `TbxScore::sent_two_hit_of`
    /// (LOSAT_X_TBXUNGAPSENT).
    #[inline(always)]
    pub(super) fn scan_unit_batched_hits<
        S: XInitSink,
        const PF: bool,
        const BF: bool,
        const SENT: usize,
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
    fn scan_unit_batched_masked<S: XInitSink, const PF: bool, const BF: bool, const SENT: usize>(
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
        const SENT: usize,
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
        let sent = SENT != 0;
        assert!(!sent || self.x_unit_sent(unit, inp));
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
            ctx_of: &self.ctx_of,
            contexts: inp.contexts,
            cutoff_scores: inp.cutoff_scores,
            x_dropoff_per_context: inp.x_dropoff_per_context,
            subject,
            subject_sent: if sent { subject_full } else { &[] },
            s_f_idx: unit.f,
            s_frame: s_frame.frame,
            s_idx: inp.s_idx,
            s_len: inp.s_len,
            sc: self.sc,
            query_batch: inp.query_batch(),
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
#[path = "x_tbx_extbatch_tests.rs"]
mod tests;
