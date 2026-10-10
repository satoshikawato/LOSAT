//! EXPERIMENT (LOSAT_X_PAIRPAR / LOSAT_X_PAIRPARSHADOW): the seed stage of one TBLASTX subject
//! (scan, two-hit test, ungapped extension) run as independent units, one per
//! (subject frame, subject chunk, query context), on the search pool.
//!
//! This file is a child module of `blast_engine/run_impl.rs` (declared there with `#[path]`, as
//! `lookup/x_lut_direct.rs` is for `backbone.rs`), so it uses the reference's chunk helpers
//! (`SubjectSplitState`, `merge_tblastx_subject_chunk_hits`, ...) without changing them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-516,614
//! ```c
//! while (scan_range[1] <= scan_range[2]) {
//!     /* scan the subject sequence for hits */
//!     hits = scansub(lookup_wrap, subject,
//!                               offset_pairs, array_size, scan_range);
//! ...
//!     for (i = 0; i < hits; ++i) {
//!         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//!         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//! ...
//!         diag_coord = (query_offset - subject_offset) & diag_mask;
//! ...
//! /* increment the offset in the diagonal array */
//! Blast_ExtendWordExit(ewp, subject->length);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805,835-844
//! ```c
//! for (context=first_context; context<=last_context; context++) {
//! ...
//!     status = s_BlastSearchEngineOneContext(program_number, query, query_info,
//! ...
//!     if (Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax)) {
//! ```
//! NCBI runs the two-hit word finder once per subject frame (and per chunk of a frame longer than
//! MAX_DBSEQ_LEN), with one scan stream that emits the hits of every query context, on one
//! diagonal table shared by the frames. This path runs the same per-hit statements
//! (`x_hit`, a copy of the reference loop body RI:3204-3546 without its trace, debug and
//! diagnostics branches) on units: for each (frame, chunk), one scan per query context with a
//! lookup table that holds only that context's entries (`XGroup`), on a fresh diagonal table
//! without masking. The unit lists are concatenated in context order and then go through the
//! reference's own steps in the reference's order: the stable init-HSP score sort
//! (`Blast_InitHitListSortByScore`), `BLAST_GetUngappedHSPList`, the chunk offset adjustment and
//! merge, and the per-frame append and sort (RI:4556-4581). Only the placement of work changes.
//!
//! # Why the units give the reference's `combined_ungapped_hits`
//!
//! Notation: QL = `query_length` (concatenated contexts with their shared sentinels), 2^k the
//! reference table (`diag_array_size`, next power of two >= QL + window), slack = 2^k - QL
//! (>= window), D the reference `diag_offset` of a chunk, len_c a context length, s_len the
//! chunk length. A hit is (q, s); its raw diagonal is d = q - s; the reference cell is
//! d & (2^k - 1), the unit cell is d itself (offset by the unit's lowest diagonal).
//!
//! 1. Cell-local state. The two-hit test of a hit reads and writes only its own cell
//!    (`last_hit` relative to D, and `flag`); cutoffs, x-drop and D are fixed during a chunk;
//!    the extension reads only the two sequences (aa_ungapped.c:516-606; RI:3685-3687). So a hit's
//!    decision depends only on the earlier hits of its cell. All comparisons are of `last_hit - D`
//!    with `s`, so a unit starting at D = window and the reference at any D agree whenever their
//!    cells hold the same `last_hit - D` and flag.
//! 2. Same hits per context, same order. A sub-lookup chain is the reference chain filtered to
//!    one context with the order kept (`build_groups`; NCBI blast_aalookup.c:337-360 layout), so a
//!    unit's scan emits exactly the reference's hit stream restricted to that context, in the same
//!    order (s ascending, chain order). Where scan calls end does not matter (per-hit state only).
//! 3. Fresh start. Every reference cell, at the first hit that a unit cell would see first,
//!    takes the same branch as a fresh cell (`DiagStruct::default()`, D = window: flag 0,
//!    `last_hit - D = -window`, so "diff >= window", `last_hit = s + D`) or a branch that ends in
//!    the same state (flag reset, `last_hit = s + D`, flag 0) without an extension:
//!    - stale values of earlier frames or chunks: `Blast_ExtendWordExit` moves D past every old
//!      value (`D += len + window`, blast_extend.c:162-173), so `last_hit - D <= -window - 1`
//!      (flag 0: too far) or `last_hit < s + D` (flag 1: reset). After the reset at
//!      D >= INT4_MAX/4 every cell is `-window` (`s_BlastDiagClear`, blast_extend.c:87-104) with
//!      D = window: `last_hit - D = -2 window`, too far. The value is read through the
//!      sign-extending `DiagStruct::last_hit` (round-1 bug: `-window` read as `2^31 - window`).
//!    - cross-sentinel diagonals: one raw diagonal can run from context c into c+1 (q grows with
//!      s, so all c hits come first). Words do not span the sentinel, so the first c+1 hit is at
//!      least 4 letters after the last c hit (`diff >= wordsize`). With flag 0, either diff >=
//!      window, or `query_offset - diff` lies in context c, the context check fails
//!      (aa_ungapped.c:566-573) and `last_hit = s + D`. With flag 1, the c extension stopped at the
//!      end of context c (the query slice ends there; NCBI: the sentinel scores BLAST_SCORE_MIN),
//!      so `s_last_off - 2 < s` of the first c+1 hit and the hit resets the cell.
//!    - masked aliasing: raw diagonals d > d' share a reference cell when d - d' = m 2^k. Both
//!      exist only if s_len >= slack + 2. Then for any hit (q, s) on d and (q', s') on d',
//!      s' - s = m 2^k - (q - q') >= slack + 3: the cell sees all hits of d, then all of d', with
//!      a gap > slack >= window. At the first d' hit, flag 0 gives diff > slack >= window (too far);
//!      flag 1 comes from an extension of d whose right end is bounded by its context,
//!      `s_last_off <= s0 + len_c` for the triggering hit s0, so the first d' hit is masked only if
//!      slack + 3 < len_c - 2: impossible when every context is at most slack long.
//!    Hence the guard (`guard`): max context length <= slack, or longest subject chunk <= slack
//!    (no aliasing at all). Otherwise the reference path runs ("fallback: guard").
//! 4. Same list. Within a unit the saved InitHSPs are the reference's saved InitHSPs of that
//!    context, in the same order. The reference sorts the chunk list with the stable
//!    `score_compare_match` (score, s_start, length, absolute q_start; blast_extend.c:274-310,
//!    RI `sort_init_hsps_by_score_ncbi`). Comparator-equal records share s_start and the absolute
//!    q_start, so they lie on one raw diagonal of one context, in one unit, in reference order.
//!    A stable sort of the concatenation c = 0..n therefore returns the reference's sorted list
//!    element for element, and every later step (`get_ungapped_hsp_list`, chunk merge, frame
//!    append and sort) is the reference's code on the same input.
//!
//! Contexts of a multi-query batch are taken in at most `X_MAX_GROUPS` contiguous groups (one
//! context per group for one query). A group is one unit: inside it the diagonals of its contexts
//! share unit cells exactly as reference cells without aliasing, and items 1-4 apply to its
//! boundaries. Units depend only on the inputs, never on the thread count; with one thread they
//! run in order on the calling thread.
//!
//! `LOSAT_X_PAIRPARSHADOW=1` runs this path and the reference loop for every subject and stops
//! at the first difference of `combined_ungapped_hits`; the summary is printed by `main`.

use super::*;
use crate::algorithm::tblastx::lookup::{BackboneCell, BlastAaLookupTable, AA_HITS_PER_CELL};
use std::sync::OnceLock;

// No NCBI counterpart: counters for LOSAT_X_STATS, the shadow summary and LOSAT_TIMING; they do not
// change any value NCBI computes.
static BATCHES: AtomicU64 = AtomicU64::new(0);
static FALLBACK_GATE: AtomicU64 = AtomicU64::new(0);
static FALLBACK_GUARD: AtomicU64 = AtomicU64::new(0);
static SUBJECTS: AtomicU64 = AtomicU64::new(0);
static UNITS: AtomicU64 = AtomicU64::new(0);
static INIT_HSPS: AtomicU64 = AtomicU64::new(0);
static SHADOW_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static SHADOW_HSPS: AtomicU64 = AtomicU64::new(0);
static NS_SUBLOOKUP: AtomicU64 = AtomicU64::new(0);
static NS_UNITS: AtomicU64 = AtomicU64::new(0);
static NS_MERGE: AtomicU64 = AtomicU64::new(0);

/// Largest number of context groups (units per subject chunk).
const X_MAX_GROUPS: usize = 6;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:264-266
// ```c
// case eAaLookupTable:
//    offset_array_size = OFFSET_ARRAY_SIZE +
//       ((BlastAaLookupTable*)lookup->lut)->longest_chain;
// ```
// The offset-pair array of a unit is sized from its own lookup table, as NCBI sizes it from the
// table it scans with (lookup_wrap.c:255-266).
const OFFSET_ARRAY_SIZE: i32 = 4096;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-504
// ```c
// while (scan_range[1] <= scan_range[2]) {
//     /* scan the subject sequence for hits */
//     hits = scansub(lookup_wrap, subject,
//                               offset_pairs, array_size, scan_range);
// ```
// Dispatch switch: 0 = the reference loop, 1 = the unit path, 2 = both, compared.
/// 0 = off, 1 = LOSAT_X_PAIRPAR, 2 = LOSAT_X_PAIRPARSHADOW.
pub(super) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_PAIRPARSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_PAIRPAR").is_some() {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: LOSAT_TIMING read once; it does not change any value NCBI computes.
fn timing() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_TIMING").is_some())
}

// No NCBI counterpart: accumulates coarse wall time for the LOSAT_TIMING line.
fn add_ns(counter: &AtomicU64, t0: Option<Instant>) {
    if let Some(t0) = t0 {
        counter.fetch_add(t0.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
    }
}

/// The summary lines (printed by `main` at exit): the shadow summary, the LOSAT_X_STATS counters
/// and one LOSAT_TIMING line.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    let m = mode();
    if m == 0 {
        return;
    }
    let load = |c: &AtomicU64| c.load(AtomicOrdering::Relaxed);
    if m == 2 {
        eprintln!(
            "[X_SHADOW] PAIRPAR subjects_compared={} units={} hsps_compared={} fallback_gate={} fallback_guard={} (combined_ungapped_hits identical)",
            load(&SHADOW_SUBJECTS),
            load(&UNITS),
            load(&SHADOW_HSPS),
            load(&FALLBACK_GATE),
            load(&FALLBACK_GUARD)
        );
    }
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] PAIRPAR batches={} subjects={} units={} init_hsps={} fallback_gate={} fallback_guard={}",
            load(&BATCHES),
            load(&SUBJECTS),
            load(&UNITS),
            load(&INIT_HSPS),
            load(&FALLBACK_GATE),
            load(&FALLBACK_GUARD)
        );
    }
    if timing() {
        eprintln!(
            "[TIMING] x_pairpar: sublookup {:.3}s units {:.3}s merge {:.3}s (subjects={} units={})",
            load(&NS_SUBLOOKUP) as f64 / 1e9,
            load(&NS_UNITS) as f64 / 1e9,
            load(&NS_MERGE) as f64 / 1e9,
            load(&SUBJECTS),
            load(&UNITS)
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-63
// ```c
// diag_array_length = 1;
// /* What power of 2 is just longer than the query? */
// while (diag_array_length < (qlen+window_size))
// {
//         diag_array_length = diag_array_length << 1;
// }
// ...
// diag_table->diag_mask = diag_array_length-1;
// ```
// The reference table has diag_array_length cells and masks the diagonal with diag_mask. The units
// may use unmasked tables only when the masking cannot change a decision (module doc, item 3):
// every context is at most `slack` long, or no subject chunk is longer than `slack`.
/// The exactness guard: `max context length <= slack || longest subject chunk <= slack`, with
/// `slack = diag_array_size - query_length` (the reference's own table).
pub(super) fn guard(
    max_context_len: usize,
    longest_subject_chunk: usize,
    query_length: i32,
    diag_array_size: i32,
) -> bool {
    let slack = diag_array_size as i64 - query_length as i64;
    max_context_len as i64 <= slack || longest_subject_chunk as i64 <= slack
}

/// The per-batch inputs of `XPairPar::prepare`.
pub(super) struct XBatch<'a> {
    pub lookup: &'a BlastAaLookupTable,
    pub contexts: &'a [QueryContext],
    pub query_length: i32,
    pub diag_array_size: i32,
    pub subjects: &'a [FastaRecord],
    /// A trace, debug or diagnostics switch is on (they exist only in the reference loop).
    pub debugging: bool,
    /// LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS or the parallel-chunk path is on.
    pub chunk_switches: bool,
}

/// The per-subject inputs of `XPairPar::run_subject` (the values the reference loop reads).
pub(super) struct XSubject<'a> {
    pub s_frames_preliminary: &'a [QueryFrame],
    pub s_frames: &'a [QueryFrame],
    pub contexts: &'a [QueryContext],
    pub cutoff_scores: &'a [i32],
    pub x_dropoff_per_context: &'a [i32],
    pub window: i32,
    pub wordsize: i32,
    pub s_idx: usize,
    pub s_len: usize,
}

/// One unit-context group: its contexts' lookup entries and their query-offset bounds.
struct XGroup {
    lookup: BlastAaLookupTable,
    q_min: u32,
    q_max: u32,
}

/// One unit: subject frame `f`, its chunk number `k`, context group `g`.
#[derive(Clone, Copy)]
struct XUnit {
    f: usize,
    k: usize,
    chunk: SubjectChunk,
    g: usize,
}

/// Per-thread scratch of the units: the offset pairs and the unit diagonal table.
#[derive(Default)]
struct XScratch {
    offset_pairs: Vec<OffsetPair>,
    diag: Vec<DiagStruct>,
}

/// The read-only values of the per-hit body (`x_hit`).
struct XHitCtx<'a> {
    diag_offset: i32,
    window: i32,
    wordsize: i32,
    lookup: &'a BlastAaLookupTable,
    contexts: &'a [QueryContext],
    cutoff_scores: &'a [i32],
    x_dropoff_per_context: &'a [i32],
    subject: &'a [u8],
    s_f_idx: usize,
    s_frame: i8,
    s_idx: usize,
    s_len: usize,
}

/// Which branch of the two-hit test a hit took (tests read it; the unit loop ignores it).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum XHit {
    Masked,
    FlagReset,
    TooFar,
    Overlap,
    ContextBoundary,
    Extended,
    Saved,
}

/// The unit path of one batch: the sub-lookups of its context groups.
pub(super) struct XPairPar {
    groups: Vec<XGroup>,
    shadow: bool,
}

impl XPairPar {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:337-360
    // ```c
    // for (i = 0; i < lookup->backbone_size; i++) {
    //   /* if there are hits there, */
    //   if (lookup->thin_backbone[i] ) {
    //       Int4 * dest = NULL;
    //       /* set the corresponding bit in the pv_array */
    //       PV_SET(pv, i, PV_ARRAY_BTS);
    //       bbc[i].num_used = lookup->thin_backbone[i][1];
    // ```
    // Dispatch point (per batch, after the lookup build): checks the gates and the guard and builds
    // the sub-lookups. None means the reference loop runs for every subject of the batch.
    /// The unit path for this batch, or None (switch off, a gate or the guard fails).
    pub(super) fn prepare(b: XBatch<'_>) -> Option<Self> {
        let m = mode();
        if m == 0 {
            return None;
        }
        BATCHES.fetch_add(1, AtomicOrdering::Relaxed);
        if b.debugging || b.chunk_switches {
            FALLBACK_GATE.fetch_add(1, AtomicOrdering::Relaxed);
            return None;
        }
        let max_context_len = b.contexts.iter().map(|c| c.aa_len).max().unwrap_or(0);
        let longest_subject_chunk = b
            .subjects
            .iter()
            .map(|r| r.seq().len() / 3)
            .max()
            .unwrap_or(0)
            .min(tblastx_max_dbseq_len_for_run());
        if !guard(
            max_context_len,
            longest_subject_chunk,
            b.query_length,
            b.diag_array_size,
        ) {
            FALLBACK_GUARD.fetch_add(1, AtomicOrdering::Relaxed);
            return None;
        }
        let t0 = timing().then(Instant::now);
        let groups = build_groups(b.lookup, b.contexts);
        add_ns(&NS_SUBLOOKUP, t0);
        Some(Self {
            groups,
            shadow: m == 2,
        })
    }

    /// Shadow mode: the reference loop runs as well and the lists are compared.
    pub(super) fn shadow(&self) -> bool {
        self.shadow
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:479-491,561-584
    // ```c
    // status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
    //                                dbseq_chunk_overlap);
    // ...
    // BlastInitHitListReset(init_hitlist);
    // ...
    // BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
    //         hit_params->options, &hsp_list);
    // ...
    // Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
    // overlap = (backup.offset == backup.hard_ranges[backup.hm_index].left) ?
    //           0 : dbseq_chunk_overlap;
    // status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:234
    // ```c
    // Blast_InitHitListSortByScore(init_hitlist);
    // ```
    // The seed stage of one subject: the units of every (frame, chunk), then, per chunk in order,
    // the concatenated unit lists through the same sort, conversion, offset adjustment and chunk
    // merge as the reference (RI:3582-3598, 4533-4553), and per frame in s_f_idx order the append
    // and sort of RI:4556-4581.
    /// `combined_ungapped_hits` of one subject, from the units.
    pub(super) fn run_subject(
        &self,
        inp: &XSubject<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Vec<UngappedHit> {
        let max_dbseq_len = tblastx_max_dbseq_len_for_run();
        let mut frame_chunks: Vec<Vec<SubjectChunk>> =
            Vec::with_capacity(inp.s_frames_preliminary.len());
        let mut units: Vec<XUnit> = Vec::new();
        for (f, s_frame) in inp.s_frames_preliminary.iter().enumerate() {
            let mut split_state = SubjectSplitState::new(s_frame.aa_len);
            let mut chunks = Vec::new();
            loop {
                match split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP) {
                    SubjectChunkStatus::Done => break,
                    SubjectChunkStatus::Ok(chunk) => chunks.push(chunk),
                }
            }
            for (k, chunk) in chunks.iter().enumerate() {
                // The reference finds no hits in a chunk shorter than a word (RI:3070-3080).
                if chunk.length < inp.wordsize as usize {
                    continue;
                }
                for (g, group) in self.groups.iter().enumerate() {
                    if group.lookup.longest_chain > 0 {
                        units.push(XUnit {
                            f,
                            k,
                            chunk: *chunk,
                            g,
                        });
                    }
                }
            }
            frame_chunks.push(chunks);
        }

        let t_units = timing().then(Instant::now);
        let outputs = self.run_units(&units, inp, pool);
        add_ns(&NS_UNITS, t_units);
        let t_merge = timing().then(Instant::now);
        SUBJECTS.fetch_add(1, AtomicOrdering::Relaxed);
        UNITS.fetch_add(units.len() as u64, AtomicOrdering::Relaxed);

        let mut outputs = outputs.into_iter();
        let mut next_unit = 0usize;
        let mut combined_ungapped_hits: Vec<UngappedHit> = Vec::new();
        for (f, chunks) in frame_chunks.iter().enumerate() {
            let mut frame_ungapped_hits: Vec<UngappedHit> = Vec::new();
            for (k, chunk) in chunks.iter().enumerate() {
                let mut init_hsps: Vec<InitHSP> = Vec::new();
                while next_unit < units.len() && units[next_unit].f == f && units[next_unit].k == k
                {
                    init_hsps.extend(outputs.next().expect("one output per unit"));
                    next_unit += 1;
                }
                INIT_HSPS.fetch_add(init_hsps.len() as u64, AtomicOrdering::Relaxed);
                sort_init_hsps_by_score_ncbi(&mut init_hsps);
                let mut hits = if init_hsps.is_empty() {
                    Vec::new()
                } else {
                    get_ungapped_hsp_list(init_hsps, inp.contexts, inp.s_frames)
                };
                adjust_tblastx_chunk_subject_offsets(&mut hits, chunk.offset);
                if !hits.is_empty() {
                    merge_tblastx_subject_chunk_hits(
                        &mut frame_ungapped_hits,
                        hits,
                        chunk.offset,
                        chunk.overlap,
                    );
                }
            }
            if !frame_ungapped_hits.is_empty() {
                combined_ungapped_hits.extend(frame_ungapped_hits);
                if !ungapped_hits_is_sorted_by_score_ncbi(&combined_ungapped_hits) {
                    sort_ungapped_hits_by_score_ncbi(&mut combined_ungapped_hits);
                }
            }
        }
        add_ns(&NS_MERGE, t_merge);
        combined_ungapped_hits
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149
    // ```c
    // TBlastThreads the_threads(GetNumberOfThreads());
    // ```
    // Scheduling only: the units run on the search pool (an indexed map, results in unit order),
    // or in order on the calling thread when the search has one thread.
    fn run_units(
        &self,
        units: &[XUnit],
        inp: &XSubject<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Vec<Vec<InitHSP>> {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() {
            return pool.install(|| {
                units
                    .par_iter()
                    .map_init(XScratch::default, |scratch, unit| {
                        self.run_unit(unit, inp, scratch)
                    })
                    .collect()
            });
        }
        let _ = pool;
        let mut scratch = XScratch::default();
        units
            .iter()
            .map(|unit| self.run_unit(unit, inp, &mut scratch))
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:494-516
    // ```c
    // scan_range[0] = 0;
    // scan_range[1] = subject->seq_ranges[0].left;
    // scan_range[2] = subject->seq_ranges[0].right - wordsize;
    // ...
    // while (scan_range[1] <= scan_range[2]) {
    //     /* scan the subject sequence for hits */
    //     hits = scansub(lookup_wrap, subject,
    //                               offset_pairs, array_size, scan_range);
    // ...
    //     for (i = 0; i < hits; ++i) {
    //         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
    //         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
    // ...
    //         diag_coord = (query_offset - subject_offset) & diag_mask;
    // ```
    // The scan loop of the reference (RI:3082-3168) on one unit: the group's sub-lookup, the frame
    // chunk, a fresh table (RI:2930-2933: `DiagStruct::default()`, diag_offset = window) indexed by
    // the raw diagonal minus the group's lowest one, without the mask (module doc, item 3).
    /// The saved init HSPs of one unit, in its scan order.
    fn run_unit(&self, unit: &XUnit, inp: &XSubject<'_>, scratch: &mut XScratch) -> Vec<InitHSP> {
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
        let mut init_hsps: Vec<InitHSP> = Vec::new();
        if subject.len() < wordsize as usize {
            return init_hsps;
        }

        let offset_array_size: i32 = OFFSET_ARRAY_SIZE + sub.longest_chain.max(0);
        if scratch.offset_pairs.len() < offset_array_size as usize {
            scratch
                .offset_pairs
                .resize(offset_array_size as usize, OffsetPair::default());
        }
        // Unit cell of hit (q, s): q - s - lo, lo = q_min - (s_len - 1) (u32 wrapping), so the
        // cells run from 0 to (q_max - q_min) + (s_len - 1).
        let width = (group.q_max - group.q_min) as usize + chunk.length;
        let lo = group.q_min.wrapping_sub(chunk.length as u32 - 1);
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
            // [C] scan_range[0] = 0;
            // [C] scan_range[1] = subject->seq_ranges[0].left;
            // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
            let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];
            // [C] while (scan_range[1] <= scan_range[2])
            while scan_range[1] <= scan_range[2] {
                let prev_scan_left = scan_range[1];
                // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
                let hits = s_blast_aa_scan_subject(
                    sub,
                    subject,
                    &seq_ranges,
                    &mut scratch.offset_pairs[..],
                    offset_array_size,
                    &mut scan_range,
                );
                if hits == 0 && scan_range[1] == prev_scan_left {
                    // The reference's safety guard (RI:3157-3161).
                    break;
                }
                let diag_ptr = scratch.diag.as_mut_ptr();
                let offset_pairs_ptr = scratch.offset_pairs.as_ptr();
                // [C] for (i = 0; i < hits; ++i)
                for i in 0..hits as usize {
                    // SAFETY: i < hits <= offset_array_size <= offset_pairs.len().
                    let pair = unsafe { &*offset_pairs_ptr.add(i) };
                    let query_offset = pair.q_off;
                    let subject_offset = pair.s_off;
                    // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
                    // Unit cell: the raw diagonal minus the unit's lowest one (no mask).
                    let diag_coord =
                        query_offset.wrapping_sub(subject_offset).wrapping_sub(lo) as usize;
                    debug_assert!(diag_coord < width, "unit diagonal outside its table");
                    // SAFETY: q_min <= query_offset <= q_max (the group's entries) and
                    // subject_offset < s_len, so diag_coord < width = scratch.diag.len().
                    let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };
                    let _ = x_hit(diag_entry, query_offset, subject_offset, &h, &mut init_hsps);
                }
            }
        }
        init_hsps
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:518-606
// ```c
// /* If the reset bit is set, an extension just happened. */
// if (diag_array[diag_coord].flag) {
//     /* If we've already extended past this hit, skip it. */
//     if ((Int4) (subject_offset + diag_offset) <
//         diag_array[diag_coord].last_hit) {
//         continue;
//     }
//     /* Otherwise, start a new hit. */
//     else {
//         diag_array[diag_coord].last_hit =
//             subject_offset + diag_offset;
//         diag_array[diag_coord].flag = 0;
//     }
// }
// /* If the reset bit is cleared, try to start an extension. */
// else {
//     /* find the distance to the last hit on this diagonal */
//     last_hit = diag_array[diag_coord].last_hit - diag_offset;
//     diff = subject_offset - last_hit;
//
//     if (diff >= window) {
//         /* We are beyond the window for this diagonal; start a
//            new hit */
//         diag_array[diag_coord].last_hit =
//             subject_offset + diag_offset;
//         continue;
//     }
//
//     /* If the difference is less than the wordsize (i.e. last
//        hit and this hit overlap), give up */
//
//     if (diff < wordsize) {
//         continue;
//     }
// ...
//     curr_context = BSearchContextInfo(query_offset, query_info);
// ...
//     if (query_offset - diff <
//         query_info->contexts[curr_context].query_offset) {
//
//         /* there was no last hit for this diagnol; start a new hit */
//         diag_array[diag_coord].last_hit =
//             subject_offset + diag_offset;
//         continue;
//     }
//
//     cutoffs = word_params->cutoffs + curr_context;
//     score = s_BlastAaExtendTwoHit(matrix, subject, query,
//                                   last_hit + wordsize,
//                                   subject_offset, query_offset,
//                                   cutoffs->x_dropoff,
//                                   &hsp_q, &hsp_s,
//                                   &hsp_len, use_pssm,
//                                   wordsize, &right_extend,
//                                   &s_last_off);
// ...
//     /* if the hsp meets the score threshold, report it */
//     if (score >= cutoffs->cutoff_score)
//         BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
//                          query_offset, subject_offset, hsp_len,
//                          score);
//
//     /* If an extension to the right happened, reset the last hit
//        so that future hits to this diagonal must start over. */
//
//     if (right_extend) {
//         diag_array[diag_coord].flag = 1;
//         diag_array[diag_coord].last_hit =
//             s_last_off - (wordsize - 1) + diag_offset;
//     }
//     /* Otherwise, make the present hit into the previous hit for
//        this diagonal */
//     else {
//         diag_array[diag_coord].last_hit =
//             subject_offset + diag_offset;
//     }
// }
// ```
// The per-hit body of the reference loop (RI:3204-3546) statement for statement, without its
// trace, debug and diagnostics branches (the gates send those runs to the reference loop) and with
// `continue` written as `return`. As in the reference, the cell is updated before the init HSP is
// saved (both read only local values, so the order is immaterial).
/// The two-hit test and extension of one hit on its cell.
#[inline(always)]
fn x_hit(
    diag_entry: &mut DiagStruct,
    query_offset: u32,
    subject_offset: u32,
    h: &XHitCtx<'_>,
    init_hsps: &mut Vec<InitHSP>,
) -> XHit {
    let diag_offset = &h.diag_offset;
    let window = h.window;
    let wordsize = h.wordsize;
    // [C] if (diag_array[diag_coord].flag)
    if diag_entry.flag() != 0 {
        // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
        let subject_plus_offset = subject_offset.wrapping_add(*diag_offset as u32);
        if subject_plus_offset < diag_entry.last_hit() as u32 {
            return XHit::Masked;
        }
        // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
        // [C] diag_array[diag_coord].flag = 0;
        diag_entry.set_last_hit(subject_plus_offset as i32);
        diag_entry.set_flag(0);
        XHit::FlagReset
    }
    // [C] else
    else {
        // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
        let last_hit = diag_entry.last_hit() - *diag_offset;
        // [C] diff = subject_offset - last_hit;
        let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;

        // [C] if (diff >= window)
        if diff >= window {
            diag_entry.set_last_hit(subject_offset.wrapping_add(*diag_offset as u32) as i32);
            return XHit::TooFar;
        }

        // [C] if (diff < wordsize)
        if diff < wordsize {
            return XHit::Overlap;
        }

        // [C] curr_context = BSearchContextInfo(query_offset, query_info);
        let ctx_idx = h.lookup.get_context_idx(query_offset as i32);
        // SAFETY: get_context_idx returns an index below num_contexts == contexts.len().
        let ctx = unsafe { h.contexts.get_unchecked(ctx_idx) };
        let q_raw = query_offset.wrapping_sub(ctx.frame_base as u32) as usize;
        let query_full = &ctx.aa_seq;
        let query = &query_full[1..query_full.len() - 1];

        // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
        let q_minus_diff = query_offset.wrapping_sub(diff as u32);
        if q_minus_diff < ctx.frame_base as u32 {
            diag_entry.set_last_hit(subject_offset.wrapping_add(*diag_offset as u32) as i32);
            return XHit::ContextBoundary;
        }

        // [C] cutoffs = word_params->cutoffs + curr_context;
        // SAFETY: ctx_idx < contexts.len() == cutoff_scores.len() == x_dropoff_per_context.len().
        let cutoff = unsafe { *h.cutoff_scores.get_unchecked(ctx_idx) };
        let x_dropoff = unsafe { *h.x_dropoff_per_context.get_unchecked(ctx_idx) };

        // [C] score = s_BlastAaExtendTwoHit(matrix, subject, query,
        //                                   last_hit + wordsize, subject_offset, query_offset, ...)
        let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) =
            extend_hit_two_hit(
                query,
                h.subject,
                (last_hit + wordsize) as usize,
                subject_offset as usize,
                q_raw,
                x_dropoff,
                false,
            );

        let hsp_q: i32 = hsp_q_u as i32;
        let hsp_s: i32 = hsp_s_u as i32;
        let hsp_len: i32 = (hsp_qe_u - hsp_q_u) as i32;
        let s_last_off: i32 = s_last_off_u as i32;

        // [C] if (right_extend) { flag = 1; last_hit = s_last_off - (wordsize - 1) + diag_offset; }
        // [C] else { last_hit = subject_offset + diag_offset; }
        if right_extend {
            diag_entry.set_flag(1);
            diag_entry.set_last_hit(s_last_off - (wordsize - 1) + *diag_offset);
        } else {
            diag_entry.set_last_hit(subject_offset.wrapping_add(*diag_offset as u32) as i32);
        }

        // [C] if (score >= cutoffs->cutoff_score)
        // [C]     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s, query_offset, subject_offset, hsp_len, score);
        if score >= cutoff {
            let hsp_q_absolute = ctx.frame_base + hsp_q;
            let hsp_qe_absolute = ctx.frame_base + (hsp_q + hsp_len);
            init_hsps.push(InitHSP {
                q_start_absolute: hsp_q_absolute,
                q_end_absolute: hsp_qe_absolute,
                s_start: hsp_s,
                s_end: hsp_s + hsp_len,
                q_seed_absolute: query_offset as i32,
                s_seed: subject_offset as i32,
                score,
                ctx_idx,
                s_f_idx: h.s_f_idx,
                q_idx: ctx.q_idx,
                s_idx: h.s_idx as u32,
                q_frame: ctx.frame,
                s_frame: h.s_frame,
                q_orig_len: ctx.orig_len,
                s_orig_len: h.s_len,
            });
            XHit::Saved
        } else {
            XHit::Extended
        }
    }
}

// No NCBI counterpart: how the contexts are grouped into units (scheduling only; module doc).
/// The first context of each group: one group per context up to `X_MAX_GROUPS` contexts, else
/// `X_MAX_GROUPS` contiguous groups balanced by length.
fn group_starts(contexts: &[QueryContext]) -> Vec<usize> {
    let n = contexts.len();
    if n <= X_MAX_GROUPS {
        return (0..n).collect();
    }
    let total: usize = contexts.iter().map(|c| c.aa_len + 1).sum();
    let mut starts = vec![0usize];
    let mut cumulative = 0usize;
    for (c, ctx) in contexts.iter().enumerate() {
        if c > 0 && cumulative * X_MAX_GROUPS >= starts.len() * total {
            starts.push(c);
            if starts.len() == X_MAX_GROUPS {
                break;
            }
        }
        cumulative += ctx.aa_len + 1;
    }
    starts
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:282-299,337-360
// ```c
// for (i = 0; i < lookup->backbone_size; i++) {
//     if (lookup->thin_backbone[i]) {
// ...
//         if (lookup->thin_backbone[i][1] > AA_HITS_PER_CELL){
// ...
//             overflow_cells_needed += lookup->thin_backbone[i][1];
//         }
//         if (lookup->thin_backbone[i][1] > longest_chain)
//             longest_chain = lookup->thin_backbone[i][1];
// ...
//   if (lookup->thin_backbone[i] ) {
//       Int4 * dest = NULL;
//       /* set the corresponding bit in the pv_array */
//       PV_SET(pv, i, PV_ARRAY_BTS);
//       bbc[i].num_used = lookup->thin_backbone[i][1];
//       /* if there are three or fewer hits, */
//       if (lookup->thin_backbone[i][1] <= AA_HITS_PER_CELL)
//           /* copy them into the thick_backbone cell */
//           dest = bbc[i].payload.entries;
//       else /* more than three hits; copy to overflow array */
//       {
//           bbc[i].payload.overflow_cursor = overflow_cursor;
//           dest = (Int4 *) lookup->overflow;
//           dest += overflow_cursor;
//           overflow_cursor += lookup->thin_backbone[i][1];
//       }
// ```
// One table per context group in the layout above: each backbone chain of the batch table
// filtered to the group's contexts (`BSearchContextInfo` ranges of `frame_bases`), order kept;
// pv bit, inline entries or overflow cursor (cells in index order), and longest_chain rebuilt.
fn build_groups(lookup: &BlastAaLookupTable, contexts: &[QueryContext]) -> Vec<XGroup> {
    /// One group's table while it is filled.
    struct Build {
        backbone: Vec<BackboneCell>,
        overflow: Vec<i32>,
        pv: Vec<u32>,
        longest_chain: i32,
        q_min: u32,
        q_max: u32,
        chain: Vec<i32>,
    }
    let starts = group_starts(contexts);
    // The group of a query offset: the last group whose first context starts at or before it.
    let bounds: Vec<i32> = starts.iter().map(|&c| lookup.frame_bases[c]).collect();
    let group_of = |q: i32| -> usize {
        let mut g = bounds.len() - 1;
        while g > 0 && q < bounds[g] {
            g -= 1;
        }
        g
    };
    let mut builds: Vec<Build> = bounds
        .iter()
        .map(|_| Build {
            backbone: vec![BackboneCell::default(); lookup.backbone.len()],
            overflow: Vec::new(),
            pv: vec![0u32; lookup.pv.len()],
            longest_chain: 0,
            q_min: u32::MAX,
            q_max: 0,
            chain: Vec::new(),
        })
        .collect();
    for (idx, full_cell) in lookup.backbone.iter().enumerate() {
        if full_cell.num_used <= 0 {
            continue;
        }
        for b in builds.iter_mut() {
            b.chain.clear();
        }
        for &q in lookup.get_hits(idx) {
            builds[group_of(q)].chain.push(q);
        }
        for b in builds.iter_mut() {
            let count = b.chain.len();
            if count == 0 {
                continue;
            }
            crate::algorithm::tblastx::lookup::pv_set(&mut b.pv, idx);
            let cell = &mut b.backbone[idx];
            cell.num_used = count as i32;
            if count <= AA_HITS_PER_CELL {
                cell.entries[..count].copy_from_slice(&b.chain);
            } else {
                cell.entries[0] = b.overflow.len() as i32;
                b.overflow.extend_from_slice(&b.chain);
            }
            b.longest_chain = b.longest_chain.max(count as i32);
            for &q in &b.chain {
                b.q_min = b.q_min.min(q as u32);
                b.q_max = b.q_max.max(q as u32);
            }
        }
    }
    builds
        .into_iter()
        .map(|b| XGroup {
            q_min: if b.longest_chain > 0 { b.q_min } else { 0 },
            q_max: b.q_max,
            lookup: BlastAaLookupTable {
                backbone: b.backbone,
                overflow: b.overflow,
                pv: b.pv,
                frame_bases: lookup.frame_bases.clone(),
                num_contexts: lookup.num_contexts,
                query_length: lookup.query_length,
                word_length: lookup.word_length,
                alphabet_size: lookup.alphabet_size,
                charsize: lookup.charsize,
                mask: lookup.mask,
                longest_chain: b.longest_chain,
                threshold: lookup.threshold,
                row_max: lookup.row_max.clone(),
            },
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:844
// ```c
// if (Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax)) {
// ```
// No computation: LOSAT_X_PAIRPARSHADOW compares the subject's combined list of the unit path with
// the reference loop's, every field, and stops at the first difference.
/// Shadow check of one subject: the unit path's `combined_ungapped_hits` against the reference's.
pub(super) fn shadow_compare(s_idx: usize, unit: &[UngappedHit], reference: &[UngappedHit]) {
    fn same(a: &UngappedHit, b: &UngappedHit) -> bool {
        let UngappedHit {
            q_idx,
            s_idx,
            ctx_idx,
            s_f_idx,
            q_frame,
            s_frame,
            q_aa_start,
            q_aa_end,
            s_aa_start,
            s_aa_end,
            q_seed_off,
            s_seed_off,
            q_orig_len,
            s_orig_len,
            raw_score,
            e_value,
            num_ident,
            hsp_list_order,
            ordering_method,
            linked_set,
            start_of_chain,
            link_id,
            chain_next_link_id,
            hsp_link_num,
            num,
        } = a;
        *q_idx == b.q_idx
            && *s_idx == b.s_idx
            && *ctx_idx == b.ctx_idx
            && *s_f_idx == b.s_f_idx
            && *q_frame == b.q_frame
            && *s_frame == b.s_frame
            && *q_aa_start == b.q_aa_start
            && *q_aa_end == b.q_aa_end
            && *s_aa_start == b.s_aa_start
            && *s_aa_end == b.s_aa_end
            && *q_seed_off == b.q_seed_off
            && *s_seed_off == b.s_seed_off
            && *q_orig_len == b.q_orig_len
            && *s_orig_len == b.s_orig_len
            && *raw_score == b.raw_score
            && e_value.to_bits() == b.e_value.to_bits()
            && *num_ident == b.num_ident
            && *hsp_list_order == b.hsp_list_order
            && *ordering_method == b.ordering_method
            && *linked_set == b.linked_set
            && *start_of_chain == b.start_of_chain
            && *link_id == b.link_id
            && *chain_next_link_id == b.chain_next_link_id
            && *hsp_link_num == b.hsp_link_num
            && *num == b.num
    }
    let n = unit.len().max(reference.len());
    for k in 0..n {
        match (unit.get(k), reference.get(k)) {
            (Some(a), Some(b)) if same(a, b) => {}
            (a, b) => panic!(
                "LOSAT_X_PAIRPARSHADOW: subject {s_idx} frame {:?} index {k} differs (unit {} hits, reference {} hits): unit {a:?} vs reference {b:?}",
                a.or(b).map(|h| h.s_f_idx),
                unit.len(),
                reference.len()
            ),
        }
    }
    SHADOW_SUBJECTS.fetch_add(1, AtomicOrdering::Relaxed);
    SHADOW_HSPS.fetch_add(reference.len() as u64, AtomicOrdering::Relaxed);
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

    /// NCBISTDAA codes of a letter string.
    fn aa(s: &str) -> Vec<u8> {
        s.bytes()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect()
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

    /// Random residues from a small alphabet (many repeated words and extensions).
    fn random_residues(rng: &mut Rng, len: usize, alphabet: &[u8]) -> Vec<u8> {
        (0..len)
            .map(|_| alphabet[rng.below(alphabet.len() as u64) as usize])
            .collect()
    }

    fn lookup_for(frames: &[QueryFrame]) -> (BlastAaLookupTable, Vec<QueryContext>) {
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        build_ncbi_lookup(&[frames.to_vec()], 11, &params, true)
    }

    fn query_length_of(contexts: &[QueryContext]) -> i32 {
        contexts
            .last()
            .map(|c| c.frame_base + c.aa_len as i32)
            .unwrap_or(0)
    }

    fn diag_array_size_for(query_length: i32, window: i32) -> i32 {
        let mut size: i32 = 1;
        while size < query_length + window {
            size <<= 1;
        }
        size
    }

    /// The reference two-hit stage of one subject on one masked table shared by the frames
    /// (RI:3019-3572): the batch lookup, `(q - s) & mask`, the same per-hit body, the diagonal
    /// offset advanced after every frame (`Blast_ExtendWordExit`). Returns the sorted init-HSP list
    /// of every frame and how often each branch was taken.
    #[allow(clippy::too_many_arguments)]
    fn masked_reference(
        lookup: &BlastAaLookupTable,
        contexts: &[QueryContext],
        s_frames: &[QueryFrame],
        cutoff_scores: &[i32],
        x_dropoff_per_context: &[i32],
        window: i32,
        diag_array: &mut [DiagStruct],
        diag_offset: &mut i32,
    ) -> (Vec<Vec<InitHSP>>, Vec<XHit>) {
        let wordsize = 3;
        let mask = (diag_array.len() - 1) as u32;
        let mut per_frame = Vec::new();
        let mut outcomes = Vec::new();
        let array_size = OFFSET_ARRAY_SIZE + lookup.longest_chain.max(0);
        let mut offset_pairs = vec![OffsetPair::default(); array_size as usize];
        for (f, s_frame) in s_frames.iter().enumerate() {
            let subject = &s_frame.aa_seq[1..s_frame.aa_seq.len() - 1];
            let mut init_hsps = Vec::new();
            if subject.len() >= wordsize as usize {
                let h = XHitCtx {
                    diag_offset: *diag_offset,
                    window,
                    wordsize,
                    lookup,
                    contexts,
                    cutoff_scores,
                    x_dropoff_per_context,
                    subject,
                    s_f_idx: f,
                    s_frame: s_frame.frame,
                    s_idx: 0,
                    s_len: 999,
                };
                let seq_ranges = [(0, subject.len() as i32)];
                let mut scan_range = [0, 0, subject.len() as i32 - wordsize];
                while scan_range[1] <= scan_range[2] {
                    let hits = s_blast_aa_scan_subject(
                        lookup,
                        subject,
                        &seq_ranges,
                        &mut offset_pairs,
                        array_size,
                        &mut scan_range,
                    );
                    for pair in &offset_pairs[..hits as usize] {
                        let cell = (pair.q_off.wrapping_sub(pair.s_off) & mask) as usize;
                        outcomes.push(x_hit(
                            &mut diag_array[cell],
                            pair.q_off,
                            pair.s_off,
                            &h,
                            &mut init_hsps,
                        ));
                    }
                }
            }
            advance_tblastx_diag_offset(diag_offset, diag_array, window, subject.len());
            sort_init_hsps_by_score_ncbi(&mut init_hsps);
            per_frame.push(init_hsps);
        }
        (per_frame, outcomes)
    }

    /// The unit path's sorted init-HSP list of every frame (one chunk per frame here).
    fn unit_lists(pp: &XPairPar, inp: &XSubject<'_>) -> Vec<Vec<InitHSP>> {
        let mut out = Vec::new();
        let mut scratch = XScratch::default();
        for (f, s_frame) in inp.s_frames_preliminary.iter().enumerate() {
            let chunk = SubjectChunk {
                offset: 0,
                length: s_frame.aa_len,
                overlap: 0,
            };
            let mut init_hsps = Vec::new();
            if chunk.length >= 3 {
                for (g, group) in pp.groups.iter().enumerate() {
                    if group.lookup.longest_chain > 0 {
                        let unit = XUnit { f, k: 0, chunk, g };
                        init_hsps.extend(pp.run_unit(&unit, inp, &mut scratch));
                    }
                }
            }
            sort_init_hsps_by_score_ncbi(&mut init_hsps);
            out.push(init_hsps);
        }
        out
    }

    struct Case {
        frames: Vec<QueryFrame>,
        s_frames: Vec<QueryFrame>,
    }

    /// Compares the unit path with the masked reference started from `start` cells at
    /// `start_offset`; returns the reference's branch outcomes.
    fn check_case(case: &Case, start: DiagStruct, start_offset: i32, stale: bool) -> Vec<XHit> {
        let window = 40;
        let (lookup, contexts) = lookup_for(&case.frames);
        let query_length = query_length_of(&contexts);
        let size = diag_array_size_for(query_length, window);
        let cutoff_scores = vec![12; contexts.len()];
        let x_dropoff_per_context = vec![16; contexts.len()];
        let pp = XPairPar {
            groups: build_groups(&lookup, &contexts),
            shadow: false,
        };
        let mut diag_array = vec![start; size as usize];
        let mut diag_offset = start_offset;
        if stale {
            // Leave the values of an earlier run of the same subject in the table.
            let _ = masked_reference(
                &lookup,
                &contexts,
                &case.s_frames,
                &cutoff_scores,
                &x_dropoff_per_context,
                window,
                &mut diag_array,
                &mut diag_offset,
            );
        }
        let (reference, outcomes) = masked_reference(
            &lookup,
            &contexts,
            &case.s_frames,
            &cutoff_scores,
            &x_dropoff_per_context,
            window,
            &mut diag_array,
            &mut diag_offset,
        );
        let inp = XSubject {
            s_frames_preliminary: &case.s_frames,
            s_frames: &case.s_frames,
            contexts: &contexts,
            cutoff_scores: &cutoff_scores,
            x_dropoff_per_context: &x_dropoff_per_context,
            window,
            wordsize: 3,
            s_idx: 0,
            s_len: 999,
        };
        let units = unit_lists(&pp, &inp);
        assert_eq!(
            units, reference,
            "unit lists differ from the masked reference"
        );
        outcomes
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:337-360
    // ```c
    // PV_SET(pv, i, PV_ARRAY_BTS);
    // bbc[i].num_used = lookup->thin_backbone[i][1];
    // ```
    #[test]
    fn sub_lookup_chains_remerge_to_the_full_chain() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let alphabet = aa("ACDEKLMNPQRSTVWY");
        let frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
            .iter()
            .map(|&f| frame(f, &random_residues(&mut rng, 300, &alphabet)))
            .collect();
        let (lookup, contexts) = lookup_for(&frames);
        let groups = build_groups(&lookup, &contexts);
        assert_eq!(groups.len(), 6);
        for idx in 0..lookup.backbone.len() {
            let full = lookup.get_hits(idx);
            let rank = |q: i32| {
                full.iter()
                    .position(|&x| x == q)
                    .expect("entry of the chain")
            };
            let mut merged: Vec<i32> = Vec::new();
            for (g, group) in groups.iter().enumerate() {
                let sub = group.lookup.get_hits(idx);
                assert_eq!(
                    crate::algorithm::tblastx::lookup::pv_test(&group.lookup.pv, idx),
                    !sub.is_empty()
                );
                assert!(sub.len() as i32 <= group.lookup.longest_chain);
                // Order kept within the sub-chain, all entries of context g.
                assert!(sub.windows(2).all(|w| rank(w[0]) < rank(w[1])));
                assert!(sub.iter().all(|&q| lookup.get_context_idx(q) == g
                    && q as u32 >= group.q_min
                    && q as u32 <= group.q_max));
                merged.extend_from_slice(sub);
            }
            merged.sort_by_key(|&q| rank(q));
            assert_eq!(merged, full, "cell {idx}");
        }
        let total: usize = groups
            .iter()
            .map(|g| {
                g.lookup
                    .backbone
                    .iter()
                    .map(|c| c.num_used as usize)
                    .sum::<usize>()
            })
            .sum();
        let full_total: usize = lookup.backbone.iter().map(|c| c.num_used as usize).sum();
        assert_eq!(total, full_total);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-61
    // ```c
    // while (diag_array_length < (qlen+window_size))
    // {
    //         diag_array_length = diag_array_length << 1;
    // }
    // ```
    #[test]
    fn guard_edges() {
        let window = 40;
        // Query length just under a power of two: slack = window, contexts far longer, long
        // subject: no unit path.
        let ql = (1 << 20) - window;
        let size = diag_array_size_for(ql, window);
        assert_eq!(size, 1 << 20);
        assert!(!guard(ql as usize / 6, 1_000_000, ql, size));
        // Same query, subject chunk within the slack: no aliasing, units exact.
        assert!(guard(ql as usize / 6, window as usize, ql, size));
        assert!(!guard(ql as usize / 6, window as usize + 1, ql, size));
        // One letter longer: the table doubles, the slack covers every context.
        let size2 = diag_array_size_for(ql + 1, window);
        assert_eq!(size2, 1 << 21);
        assert!(guard((ql as usize + 1) / 6, 5_000_000, ql + 1, size2));
        // Context exactly the slack long passes, one more letter fails.
        let slack = (size2 - (ql + 1)) as usize;
        assert!(guard(slack, 5_000_000, ql + 1, size2));
        assert!(!guard(slack + 1, 5_000_000, ql + 1, size2));
        // p03 and d04 (map §5).
        assert!(guard(95_687, 98_048, 574_123, 1 << 20));
        assert!(guard(220_702, 202_064, 1_324_217, 1 << 21));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560-573
    // ```c
    // curr_context = BSearchContextInfo(query_offset, query_info);
    // ...
    // if (query_offset - diff <
    //     query_info->contexts[curr_context].query_offset) {
    // ```
    #[test]
    fn context_boundary_diagonal() {
        // Context 0 ends with WWW, context 1 starts with CCC: on the raw diagonal of the subject
        // "WWWGCCC..." the CCC hit of context 1 follows the WWW hit of context 0 by 4 letters
        // and fails the context check in the reference; in the unit it is the first hit of a
        // fresh cell. A second CCC word later on the same diagonal then extends in both.
        let mut c0 = aa("KLMNPQRSTVWW");
        c0.extend(aa("WWW"));
        let c1 = aa("CCCHHHIIIEEECCC");
        let frames = vec![frame(1, &c0), frame(2, &c1)];
        let s = aa("AAWWWWCCCHHHIIIEEECCCAA");
        let case = Case {
            frames,
            s_frames: vec![frame(1, &s)],
        };
        let outcomes = check_case(&case, DiagStruct::default(), 40, false);
        assert!(outcomes.contains(&XHit::ContextBoundary), "{outcomes:?}");
        assert!(outcomes.contains(&XHit::Saved), "{outcomes:?}");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-104
    // ```c
    // diag->offset = diag->window;
    // ...
    //     diag_struct_array[i].flag = 0;
    //     diag_struct_array[i].last_hit = -diag->window;
    // ```
    #[test]
    fn cleared_cells_hold_minus_window() {
        let window = 40;
        let cell = DiagStruct::clear(window);
        // Sign extension of the 31-bit field (round-1 bug: 2^31 - window).
        assert_eq!(cell.last_hit(), -window);
        assert_eq!(cell.flag(), 0);
        // A hit on a cleared cell at diag_offset = window is "too far", like a fresh unit cell.
        let mut rng = Rng(7);
        let alphabet = aa("ACDEKW");
        let q = random_residues(&mut rng, 60, &alphabet);
        let s = random_residues(&mut rng, 80, &alphabet);
        let case = Case {
            frames: vec![frame(1, &q), frame(-1, &q[10..])],
            s_frames: vec![frame(1, &s), frame(-1, &s[5..])],
        };
        let outcomes = check_case(&case, DiagStruct::clear(window), window, false);
        assert!(outcomes.contains(&XHit::Saved));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
    // ```c
    // if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //    ewp->diag_table->offset = ewp->diag_table->window;
    //    s_BlastDiagClear(ewp->diag_table);
    // } else {
    //    ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    // }
    // ```
    #[test]
    fn diag_offset_near_the_reset_limit() {
        let mut rng = Rng(11);
        let alphabet = aa("ACDEKW");
        let q = random_residues(&mut rng, 70, &alphabet);
        let s = random_residues(&mut rng, 90, &alphabet);
        let s_frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
            .iter()
            .map(|&f| frame(f, &s[(f.unsigned_abs() as usize)..]))
            .collect();
        let case = Case {
            frames: vec![frame(1, &q), frame(2, &q[3..]), frame(-1, &q[7..])],
            s_frames,
        };
        // Stale values below the limit, then frames that cross it (reset to window mid-subject).
        for start in [i32::MAX / 4 - 400, i32::MAX / 4 - 100, i32::MAX / 4] {
            check_case(&case, DiagStruct::default(), start, true);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
    // ```c
    // diag_coord = (query_offset - subject_offset) & diag_mask;
    // ...
    // if (diag_array[diag_coord].flag) {
    // ```
    #[test]
    fn random_cases_with_and_without_aliasing() {
        let mut rng = Rng(0x1234_5678_9abc_def1);
        let alphabets = [aa("ACDEKW"), aa("AW"), aa("ACDEFGHIKLMNPQRSTVWY")];
        let mut aliased = 0usize;
        let mut checked = 0usize;
        for case_no in 0..300 {
            let alphabet = &alphabets[case_no % alphabets.len()];
            let n_ctx = 1 + rng.below(6) as usize;
            let frames: Vec<QueryFrame> = (0..n_ctx)
                .map(|c| {
                    let len = 3 + rng.below(40) as usize;
                    frame(
                        [1i8, 2, 3, -1, -2, -3][c],
                        &random_residues(&mut rng, len, alphabet),
                    )
                })
                .collect();
            let s_len = 3 + rng.below(400) as usize;
            let s_frames: Vec<QueryFrame> = (0..1 + rng.below(3) as usize)
                .map(|f| {
                    frame(
                        [1i8, 2, 3][f],
                        &random_residues(&mut rng, s_len - f, alphabet),
                    )
                })
                .collect();
            let (_, contexts) = lookup_for(&frames);
            let window = 40;
            let query_length = query_length_of(&contexts);
            let size = diag_array_size_for(query_length, window);
            let max_ctx = contexts.iter().map(|c| c.aa_len).max().unwrap_or(0);
            if !guard(max_ctx, s_len, query_length, size) {
                continue;
            }
            if s_len as i64 > (size - query_length) as i64 {
                aliased += 1;
            }
            let case = Case { frames, s_frames };
            let start = match case_no % 3 {
                0 => (DiagStruct::default(), window),
                1 => (DiagStruct::clear(window), window),
                _ => (DiagStruct::default(), i32::MAX / 4 - 50),
            };
            check_case(&case, start.0, start.1, case_no % 2 == 1);
            checked += 1;
        }
        assert!(checked > 100, "checked {checked}");
        assert!(aliased > 20, "aliased {aliased}");
    }

    // No NCBI counterpart: grouping of contexts (scheduling only).
    #[test]
    fn many_contexts_make_six_contiguous_groups() {
        let mut rng = Rng(3);
        let alphabet = aa("ACDEKW");
        let queries: Vec<Vec<QueryFrame>> = (0..4)
            .map(|_| {
                [1i8, 2, 3, -1, -2, -3]
                    .iter()
                    .map(|&f| frame(f, &random_residues(&mut rng, 20, &alphabet)))
                    .collect()
            })
            .collect();
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (_lookup, contexts) = build_ncbi_lookup(&queries, 11, &params, true);
        let starts = group_starts(&contexts);
        assert_eq!(starts.len(), X_MAX_GROUPS);
        assert_eq!(starts[0], 0);
        assert!(starts.windows(2).all(|w| w[0] < w[1]));
        assert!(*starts.last().unwrap() < contexts.len());
    }
}
