//! EXPERIMENT (PAIRPAR, graduated on 2026-10-22: no variable; the rule is `decide` below): the seed
//! stage of one TBLASTX subject
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
//!    cells hold the same `last_hit - D` and flag, and `s + D` fits the 31-bit `last_hit` in the
//!    reference as it does in the unit (`last_hit_fits`; it fails only for a `-window_size` near
//!    2^29).
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
//!      a gap > slack >= window. At the first d' hit, flag 0 gives diff > slack >= window (too far).
//!      Flag 1 comes from the last right extension of d, started by a hit (q0, s0), which wrote
//!      `last_hit = s_last_off - (wordsize - 1) + D`. Its right walk ends at the end of q0's
//!      context at the latest (`extend_hit_two_hit` walks at most `len_c - q_raw` letters right of
//!      the word start, the query slice ends there), so `s_last_off - s0 <= QL - q0`. A hit
//!      (q', s') on d' has q' >= 0, so `s' - s0 = q' - q0 + m 2^k >= 2^k - q0 >= QL + window - q0`,
//!      and `s' >= s_last_off + window`: the first d' hit is never masked; it resets the cell.
//!    Hence no condition on the lengths: with LOSAT_X_TBXBAND (section below) the units take every
//!    batch the gates let through and `last_hit_fits` holds. Without the switch the guard of the
//!    earlier, weaker argument (`guard`: max context length <= slack, or longest subject chunk
//!    <= slack) sends the other batches to the reference path ("fallback: guard"). That argument
//!    bounded `s_last_off - s0` by `len_c` and `s' - s0` by `slack + 3` apart, and so needed
//!    `len_c <= slack`.
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
//! `LOSAT_X_TBXBAND` (in the strict set: on unless `LOSAT_X_TBXBAND` is `0`)
//! : `decide` does not apply `guard` (item 3) to a batch whose longest subject
//! chunk has at most as many letters as the reference table has cells (`x_band_fits`, a memory
//! rule: the units' tables grow with the subject chunk). The batches the guard let through change
//! nothing. The ones it rejected have a longest context and a longest subject chunk both longer than
//! the slack: for one query of N bases, N in (3/7 2^k, 2^(k-1)] for the table of 2^k cells (1.76-2.0
//! kb, 7.0-8.2 kb, ..., 56-65 kb, 112-131 kb, 225-262 kb, 0.45-0.52 Mb, 0.90-1.05 Mb, 1.8-2.1 Mb,
//! 3.6-4.2 Mb), and multi-query batches whose letters lie just below a power of two. With the rule,
//! those against subject frames of at most 2^k letters (genome pairs of similar size) go through
//! `LOSAT_X_PAIRPARGATE` like any other batch instead of always going to the reference frame loop,
//! which searches a subject on one thread; those against longer subjects (a gene against a genome)
//! keep the reference loop without LOSAT_X_TBXUNITMASK and take the units on masked tables with it
//! (section below). The `[X_STATS] PAIRPAR` line counts the band batches
//! the units took (`band`).
//!
//! PAIRLUT (graduated on 2026-10-22: no variable): for a batch that takes the unit path the group
//! tables are built directly from the
//! query (`lookup/x_pair_lut.rs`, exactness argument in its header) and the batch's full table is
//! replaced by a stub that nothing scans. Whether the unit path runs is decided before the lookup
//! build by `decide`, the function `prepare` decides with (`lut_build`); `prepare` stops the run if
//! direct groups exist for a batch it would leave to the reference loop, and the per-subject
//! dispatch stops it if the reference frame loop would run on the stub. If `lut_build` gives no
//! groups for a batch that takes the units, `prepare` builds them from the full table with
//! `build_groups` (the fallback).
//!
//! `LOSAT_X_SEEDCOMPACT` (in the strict set, `utils/xswitch.rs`: on unless `LOSAT_X_SEEDCOMPACT` is
//! `0`): a unit keeps its saved init HSPs
//! as 28-byte records, NCBI's own record shape, copied out of a per-thread vector at their exact
//! length; the merge expands them into today's `InitHSP`s, converts each chunk's list into a list
//! of its exact length and moves the first chunk's list into the empty frame list. The section
//! above `XUnitHsp` gives the exactness argument. The unit test `compact_check` compares, on the
//! unit path, the compact records' expansion with the full lists and the merge of the compact
//! records with the merge of the full lists (both sides are the unit path); the per-unit hits are
//! compared with `masked_reference`. The unit path against the frame loop is checked only through
//! the frozen NCBI bytes, at 1 thread (frame loop) and 4 threads (units).
//!
//! `LOSAT_X_SEEDUNMAP` (in the strict set: on unless `LOSAT_X_SEEDUNMAP` is `0`;
//! it acts on Linux with glibc only, where transparent huge
//! pages do not round mapped blocks up (`x_thp_rounds_up`), and only on the search pool, at two or
//! more threads): where a unit's query-group span and subject chunk are both at least
//! `X_MAPPED_MIN_SIDE` letters (`x_unit_mapped`), the unit diagonal table and the record vector of
//! a thread's split and the unit's copy of its records get a capacity that glibc always maps, so
//! they are returned to the system when dropped instead of staying in the worker threads' arenas
//! through the linking call. Only capacities change; the section above `X_MAPPED_BYTES` gives the
//! reason, the three rules (sides, huge pages, one thread) and why no value changes.
//!
//! `LOSAT_X_TBXPREFETCH` (in the strict set, a value switch: on unless it is `0`;
//! x86-64 only): the units of a subject whose largest unit
//! diagonal table has at least the switch's number of KiB (built-in 1,536, `X_PREFETCH_MIN_KIB`)
//! run the scan with a prefetch of the cell of the hit `X_PREFETCH_AHEAD` entries ahead in the
//! batch (`scan_unit` with `PF`). A prefetch reads nothing into the program and has no
//! architectural effect, so the cells, the hits and the init HSPs are those without it. The other
//! subjects run the scan without the hint (`PF` false), whose machine code is the one before the
//! switch.
//!
//! `LOSAT_X_TBXQSLICE` (in the strict set, a value switch: on unless it is `0`;
//! x86-64 only): unset (or `=1`), a context group of more than
//! 128 Ki letters whose span times the search's threads reaches 6 Mi letters is cut into slices of
//! at most 128 Ki consecutive query offsets; `=N` cuts every group of more than N letters into
//! slices of at most N (at most 64 slices a group). Each slice has its own lookup table, and the
//! units of a batch with a cut group scan their chunk once per slice on their one diagonal table
//! (x_tbx_qslice.rs, rule and exactness argument in its header): the table cells a slice's scan
//! touches at one subject position span the slice, not the group. Each cell sees its hits in
//! the same order, so the table and the chunk's sorted init HSPs are those of `scan_unit`. A batch
//! with no cut group runs `scan_unit`, whose machine code is the one before the switch.
//!
//! `LOSAT_X_TBXUNITMASK` (in the strict set, a value switch: on unless it is `0`;
//! unset or `=1` the built-in 256, `=N` N KiB): a unit takes
//! a table of P cells, P the smallest power of two of at least its group's extent G (the query
//! offsets from the group's first context to the next group's, the last group to the end of its
//! last context, `x_group_extents`) plus the window, and the cell of hit (q, s) is
//! (q - s - lo) & (P - 1) (`XMasked`), where its unmasked table (`x_unit_width`) has at least four
//! times P cells and P is at most the switch's KiB at eight threads or more, at most 4,096 cells at
//! two to seven, never at one (`x_mask_max_cells`), or, in a batch that `guard` rejects
//! (LOSAT_X_TBXBAND), where its subject chunk is longer than the reference table has cells
//! (`x_unit_mask_for`). Measured (workstation node 1): a unit's hit loop
//! does one more instruction per hit on a masked table and saves the clearing of the unmasked
//! table's every cell per unit, so masking pays where the group is short against the chunk (1.5 kb
//! against MG1655 and the fruit fly genome database x0.80-0.85 at 16 threads, max RSS 784 -> 309
//! MiB; Sakai 1-100,000 x MG1655 x0.987 at 16, 143 -> 61 MiB, where NCBI BLAST+ takes 130 MiB),
//! costs 0.8-1.7 % with groups of 100,000 letters or more (Sakai 1-300,000 x MG1655 with
//! `=512`). With fewer threads the clearing weighs less: at two and four threads tables of 1-4 Ki
//! cells gain (x0.895-0.971) and those of 32-64 Ki cells lose 0.2-0.9 %, and at one thread, where a
//! unit's table is the only one, Sakai 1-100,000 loses 0.6 %.
//! The units of a batch whose groups LOSAT_X_TBXQSLICE's rule cuts keep their unmasked tables
//! (`x_tbx_qslice::rule_cuts`, whatever the diagnostics switches say). With the switch, `decide`
//! does not apply `x_band_fits` to a batch whose groups the rule does not cut: every unit whose
//! chunk is longer than the reference table then takes a masked table, of at most the reference
//! table's cells. Why the units' init HSPs do not change, item 3's masked-aliasing argument for the
//! group: raw diagonals d > d' of the unit share a masked cell when d - d' = m P (m >= 1); every hit
//! of the group has q in [a, a + G), so for hits (q, s) on d and (q', s') on d',
//! s' - s = m P - (q - q') >= P - G + 1 > window, and the cell sees all hits of d before those of
//! d'. Flag 0: the cell holds `s_d + D` of a hit on d and diff = s' - s_d > window, too far, the
//! state of a fresh cell after its first hit. Flag 1: the last right extension of d, from (q0, s0),
//! stopped at the end of q0's context at the latest, at most a + G, so s_last_off - s0 <= a + G -
//! q0, while s' - s0 = q' - q0 + m P >= a - q0 + P: s' >= s_last_off + window, the hit is not masked,
//! it resets the cell to the fresh state. So the first hit of d' leaves its cell as on an unmasked
//! table, and by induction every hit takes the branch it takes there. The scan in slices cannot use
//! such a table: a slice's pass over the whole chunk brings hits of d' (low q') before the later
//! slices' hits of d. The masked loops are their own instances (`scan_unit_masked`,
//! `scan_unit_batched_masked`); the unmasked ones are the code before the switch. The
//! `[X_STATS] TBXUNITMASK` line counts the units with a masked table.

use super::*;
use crate::algorithm::tblastx::blast_gapalign::get_ungapped_hsp_list_into;
use crate::algorithm::tblastx::lookup::x_pair_lut::{XPairLutGroup, XPairLutQuery};
use crate::algorithm::tblastx::lookup::{BackboneCell, BlastAaLookupTable, AA_HITS_PER_CELL};
use std::sync::OnceLock;

// Diagnostic (LOSAT_TIMING, LOSAT_X_STATS): what the scan of a unit records besides its init
// HSPs (the seed / extension split of the units' time, the diagonal statistics). No NCBI
// counterpart.
#[path = "x_unit_probe.rs"]
mod x_unit_probe;

// EXPERIMENT (LOSAT_X_TBXKMERGE): the merge of a subject's unit lists
// without the full sorts (`merge_subject`). A child module of this file, as x_unit_probe is.
#[path = "x_kmerge.rs"]
mod x_kmerge;
use x_unit_probe::{XNoProbe, XProbe, XProbeSlot, XStatsProbe, XTimeProbe};

// EXPERIMENT (LOSAT_X_TBXEXTBATCH): the two-hit stage of a unit with its extensions held back and
// resolved in batches (x_tbx_extbatch.rs).
#[path = "x_tbx_extbatch.rs"]
mod x_tbx_extbatch;

// EXPERIMENT (LOSAT_X_TBXHITBF): the hit loop of a unit with the too-far
// and overlap outcomes of the two-hit test taken without a branch (x_tbx_hitbf.rs).
#[path = "x_tbx_hitbf.rs"]
mod x_tbx_hitbf;
// EXPERIMENT (LOSAT_X_TBXQSLICE): the scan of a unit in slices of its
// query group, one lookup table per slice, on the unit's one diagonal table (x_tbx_qslice.rs).
#[path = "x_tbx_qslice.rs"]
mod x_tbx_qslice;
use x_tbx_qslice::XSliced;

// No NCBI counterpart: counters for LOSAT_X_STATS and LOSAT_TIMING; they do not
// change any value NCBI computes.
static BATCHES: AtomicU64 = AtomicU64::new(0);
static FALLBACK_GATE: AtomicU64 = AtomicU64::new(0);
static FALLBACK_GUARD: AtomicU64 = AtomicU64::new(0);
static FALLBACK_COST: AtomicU64 = AtomicU64::new(0);
static SUBJECTS: AtomicU64 = AtomicU64::new(0);
static UNITS: AtomicU64 = AtomicU64::new(0);
static INIT_HSPS: AtomicU64 = AtomicU64::new(0);
static NS_SUBLOOKUP: AtomicU64 = AtomicU64::new(0);
static NS_UNITS: AtomicU64 = AtomicU64::new(0);
static NS_MERGE: AtomicU64 = AtomicU64::new(0);
// No NCBI counterpart: counters of PAIRLUT (batches whose groups were built directly, their groups,
// batches left to the reference builder, wall time of the direct build); they do not change any
// value NCBI computes.
static LUT_BATCHES: AtomicU64 = AtomicU64::new(0);
static LUT_GROUPS: AtomicU64 = AtomicU64::new(0);
static LUT_FALLBACK: AtomicU64 = AtomicU64::new(0);
static NS_LUT_DIRECT: AtomicU64 = AtomicU64::new(0);
// No NCBI counterpart: counters of LOSAT_X_SEEDCOMPACT (subjects merged from compact unit records,
// their units, the records expanded); they do not change any value NCBI computes.
static SC_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static SC_UNITS: AtomicU64 = AtomicU64::new(0);
static SC_RECORDS: AtomicU64 = AtomicU64::new(0);
// No NCBI counterpart: counters of LOSAT_X_SEEDUNMAP (blocks given the mapped capacity: unit
// diagonal tables, record vectors of a split, unit copies of the records); they do not change any
// value NCBI computes.
static SU_DIAG: AtomicU64 = AtomicU64::new(0);
static SU_RECORDS: AtomicU64 = AtomicU64::new(0);
static SU_COPIES: AtomicU64 = AtomicU64::new(0);
// No NCBI counterpart: the counter of LOSAT_X_TBXBAND (batches the units took that `guard`
// rejects); it does not change any value NCBI computes.
static BAND: AtomicU64 = AtomicU64::new(0);
// No NCBI counterpart: the counter of LOSAT_X_TBXUNITMASK (units that took a masked table); it does
// not change any value NCBI computes.
static UM_UNITS: AtomicU64 = AtomicU64::new(0);

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
// No NCBI counterpart: a LOSAT variable that chooses, per batch, between two LOSAT paths with the
// same result (the unit path and the reference frame loop, module doc items 1-4); NCBI has one path.
/// LOSAT_X_PAIRPARGATE: a strict-set switch (on unless it is `0`; `utils/xswitch.rs`).
pub(super) fn cost_gate() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_PAIRPARGATE"))
}

// No NCBI counterpart: a LOSAT variable that moves the line between the same two LOSAT paths at one
// thread (`units_gain`); NCBI has one path.
/// LOSAT_X_PAIRPARN1: a strict-set switch (on unless it is `0`; `utils/xswitch.rs`).
pub(super) fn n1_gate() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_PAIRPARN1"))
}

// No NCBI counterpart: a LOSAT variable that chooses, per batch, between two LOSAT paths with the
// same result (the unit path and the reference frame loop, module doc item 3); NCBI has one path.
// Read once per process, as the other modes of this file.
/// LOSAT_X_TBXBAND: a strict-set switch (on unless it is `0`;
/// `utils/xswitch.rs`). With it `decide` does not apply `guard` (module doc item 3).
fn tbx_band() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXBAND"))
}

// No NCBI counterpart: a LOSAT variable that chooses, per unit, between two tables with the same
// result (module doc, LOSAT_X_TBXUNITMASK); NCBI has one table per thread. Read once per process,
// as the other modes of this file.
/// LOSAT_X_TBXUNITMASK: a strict-set value switch (`utils/xswitch.rs`). `None` when it is off (`0`);
/// otherwise the largest masked table, in KiB, that a
/// unit takes for its size alone at `X_UNITMASK_MANY_THREADS` threads or more (`x_mask_max_cells`,
/// module doc).
fn unit_mask_kib() -> Option<usize> {
    static KIB: OnceLock<Option<usize>> = OnceLock::new();
    *KIB.get_or_init(|| {
        crate::utils::xswitch::value(
            "LOSAT_X_TBXUNITMASK",
            2..=X_UNITMASK_MAX_KIB,
            X_UNITMASK_KIB,
        )
    })
}

/// No NCBI counterpart: LOSAT_X_TBXUNITMASK's built-in, in KiB (module doc): masked tables of up to
/// 65,536 cells (group extents plus the window of up to 64 Ki letters).
const X_UNITMASK_KIB: usize = 256;

/// No NCBI counterpart: the largest value LOSAT_X_TBXUNITMASK takes (4 GiB); a value outside
/// 2..=X_UNITMASK_MAX_KIB other than `0`, or a non-numeric one, gives the built-in.
const X_UNITMASK_MAX_KIB: usize = 1 << 22;

/// No NCBI counterpart: LOSAT_X_TBXUNITMASK takes a masked table for its size only where the
/// unmasked one has at least this many times its cells (module doc).
const X_UNITMASK_RATIO: usize = 4;

/// No NCBI counterpart: from this many threads LOSAT_X_TBXUNITMASK's size rule takes masked tables
/// of up to the switch's KiB; below it, of up to `X_UNITMASK_FEW_THREADS_CELLS` (module doc).
const X_UNITMASK_MANY_THREADS: usize = 8;

/// No NCBI counterpart: the largest masked table, in cells (16 KiB), of LOSAT_X_TBXUNITMASK's size
/// rule at 2 to `X_UNITMASK_MANY_THREADS - 1` threads (module doc).
const X_UNITMASK_FEW_THREADS_CELLS: usize = 1 << 12;

// No NCBI counterpart: the size limit of LOSAT_X_TBXUNITMASK's rule (module doc).
/// The largest masked table, in cells, that a unit takes for its size at `threads` threads, for
/// the switch's `kib`: none at one thread, at most `X_UNITMASK_FEW_THREADS_CELLS` below
/// `X_UNITMASK_MANY_THREADS`, the switch's KiB from there.
fn x_mask_max_cells(kib: usize, threads: usize) -> usize {
    let cells = kib.saturating_mul(1024) / std::mem::size_of::<DiagStruct>();
    match threads {
        0 | 1 => 0,
        t if t < X_UNITMASK_MANY_THREADS => cells.min(X_UNITMASK_FEW_THREADS_CELLS),
        _ => cells,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:510
// ```c
//         diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// The cell of a hit in a unit's table: the raw diagonal minus the unit's lowest one, without a
// mask (`XUnmasked`, module doc item 3), or, with LOSAT_X_TBXUNITMASK, that value masked to the
// unit's table of `x_masked_cells` cells (`XMasked`, module doc). One type per form, so that the
// hit loops are compiled once for each and the unmasked loops keep their machine code.
/// How a unit maps a hit (q, s) to the cell of its table.
pub(super) trait XCells: Copy {
    /// Whether the table is masked (`XMasked`).
    const MASKED: bool;
    /// The cells of a unit whose lowest raw diagonal is `lo`; `mask` is the masked table's
    /// number of cells minus one (`XMasked`) or ignored (`XUnmasked`).
    fn of(lo: u32, mask: u32) -> Self;
    /// The cell of hit (q, s).
    fn cell(self, query_offset: u32, subject_offset: u32) -> usize;
}

/// The unmasked cells: q - s - lo (u32 wrapping).
#[derive(Clone, Copy)]
pub(super) struct XUnmasked {
    pub lo: u32,
}

impl XCells for XUnmasked {
    const MASKED: bool = false;
    #[inline(always)]
    fn of(lo: u32, _mask: u32) -> Self {
        Self { lo }
    }
    #[inline(always)]
    fn cell(self, query_offset: u32, subject_offset: u32) -> usize {
        query_offset
            .wrapping_sub(subject_offset)
            .wrapping_sub(self.lo) as usize
    }
}

/// LOSAT_X_TBXUNITMASK's cells: (q - s - lo) & mask (u32 wrapping), mask + 1 a power of two.
#[derive(Clone, Copy)]
pub(super) struct XMasked {
    pub lo: u32,
    pub mask: u32,
}

impl XCells for XMasked {
    const MASKED: bool = true;
    #[inline(always)]
    fn of(lo: u32, mask: u32) -> Self {
        Self { lo, mask }
    }
    #[inline(always)]
    fn cell(self, query_offset: u32, subject_offset: u32) -> usize {
        (query_offset
            .wrapping_sub(subject_offset)
            .wrapping_sub(self.lo)
            & self.mask) as usize
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-56
// ```c
// diag_array_length = 1;
// /* What power of 2 is just longer than the query? */
// while (diag_array_length < (qlen+window_size))
// {
//         diag_array_length = diag_array_length << 1;
// }
// ```
/// LOSAT_X_TBXUNITMASK: the cells of the masked table of a unit whose query group extends over
/// `extent` letters (`x_group_extents`): the smallest power of two of at least extent + window, as
/// NCBI sizes its table from the whole query.
fn x_masked_cells(extent: usize, window: i32) -> usize {
    (extent + window.max(0) as usize).next_power_of_two()
}

// No NCBI counterpart: the rule of LOSAT_X_TBXUNITMASK for one batch (module doc), apart from the
// switch so that a test can give it values.
/// Which units of a batch take a masked table (`x_unit_mask_for`).
#[derive(Clone, Copy)]
struct XMaskRule {
    /// A unit whose masked table has at most this many cells takes it where its unmasked table has
    /// at least `ratio` times as many (`X_UNITMASK_RATIO`)...
    max_cells: usize,
    ratio: usize,
    /// ... and, in a batch that `guard` rejects (LOSAT_X_TBXBAND), a unit whose subject chunk is
    /// longer than the reference table's cells (these) takes it whatever its size (`x_band_fits`).
    band_cells: Option<i32>,
}

// No NCBI counterpart: the rule of LOSAT_X_TBXUNITMASK (module doc).
/// The mask of the table of a unit whose group extends over `extent` letters, whose unmasked
/// table has `width` cells and whose subject chunk has `chunk_len` letters, under `rule`: the
/// masked table's cells minus one, or None (the unmasked table).
fn x_unit_mask_for(
    rule: &XMaskRule,
    extent: usize,
    window: i32,
    width: usize,
    chunk_len: usize,
) -> Option<u32> {
    let cells = x_masked_cells(extent, window);
    let small = cells <= rule.max_cells && width as u64 >= rule.ratio as u64 * cells as u64;
    let band = rule
        .band_cells
        .is_some_and(|diag_array_size| !x_band_fits(chunk_len, diag_array_size));
    (small || band).then(|| (cells - 1) as u32)
}

// No NCBI counterpart: the groups' query offsets as LOSAT_X_TBXQSLICE reads them (from the first
// context of a group to the first of the next, the last one to the end of its last context).
/// The number of query offsets each group of `contexts` extends over (`group_starts`).
fn x_group_extents(contexts: &[QueryContext]) -> Vec<usize> {
    x_tbx_qslice::x_group_spans(contexts, &group_starts(contexts))
        .into_iter()
        .map(|(a, b)| (b - a) as usize)
        .collect()
}

// No NCBI counterpart: chooses how LOSAT holds a unit's saved init HSPs between the units and the
// merge (as NCBI's own record, `XUnitHsp`, or as today's `InitHSP`); it does not change any value
// NCBI computes. Read once per process, as the other modes of this file.
/// A strict-set switch (on unless `LOSAT_X_SEEDCOMPACT` is `0`; `utils/xswitch.rs`).
fn seed_compact() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_SEEDCOMPACT"))
}

// No NCBI counterpart: chooses the capacity LOSAT gives three blocks of the unit path so that glibc
// maps them (section above `X_MAPPED_BYTES`); it does not change any value NCBI computes. Read once
// per process, as the other modes of this file.
/// A strict-set switch (on unless `LOSAT_X_SEEDUNMAP` is `0`;
/// `utils/xswitch.rs`). It joined the set on 2026-10-17 in the change that added it, on by
/// default under the Owner's rule of 2026-10-17 (a speed-up or memory reduction lands on by default).
fn seed_unmap() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_SEEDUNMAP"))
}

// No NCBI counterpart: whether LOSAT_X_SEEDUNMAP acts in this process (section above
// `X_MAPPED_BYTES`): the switch is on, the target is Linux with glibc, and transparent huge pages
// do not round mapped blocks up (`x_thp_rounds_up`; the two files are read only when the first two
// hold). Read once per process.
/// LOSAT_X_SEEDUNMAP as it acts here: the value the scratches of the search pool get.
fn seed_unmap_active() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        seed_unmap()
            && cfg!(all(target_os = "linux", target_env = "gnu"))
            && !x_thp_rounds_up(
                std::fs::read_to_string("/sys/kernel/mm/transparent_hugepage/enabled")
                    .ok()
                    .as_deref(),
                std::fs::read_to_string("/proc/self/status").ok().as_deref(),
            )
    })
}

// EXPERIMENT (LOSAT_X_TBXPREFETCH, module doc): the unit scan prefetches the cell of the hit this
// many entries ahead in its `scansub` batch. 16 hits are about 16 x 4 ns of scan on the measured
// CPUs, longer than a load from L3, where the tables that gained fit. No NCBI counterpart.
const X_PREFETCH_AHEAD: usize = 16;
// The built-in number of the switch: the smallest largest-unit table, in KiB, of a subject whose
// units prefetch. Measured with the hint on every table (A/B in docs/agents/diagnostics.md, entry
// `LOSAT_X_TBXPREFETCH`): tables of 1.65 MiB and more gained on both measured CPUs (L2 of one core
// 1.25 MiB on the Xeon 6338, 2 MiB on the i9-14900HX), 1.18 MiB lost 0.7 %; 1.5 MiB lies between.
const X_PREFETCH_MIN_KIB: usize = 1536;
// The largest number the switch takes (4 GiB of table); a value outside 2..=X_PREFETCH_MAX_KIB
// other than `0`, or a non-numeric one, gives the built-in (`xswitch::value`).
const X_PREFETCH_MAX_KIB: usize = 1 << 22;

// No NCBI counterpart: a cache hint in LOSAT's unit scan (module doc); NCBI's loop gives none, and
// a prefetch has no architectural effect, so no value NCBI computes changes. Read once per process,
// as the other modes of this file.
/// LOSAT_X_TBXPREFETCH: a strict-set value switch (`utils/xswitch.rs`). `None` when it is off (`0`)
/// and on targets other than x86-64; otherwise the
/// smallest largest-unit table, in KiB, of a subject whose units prefetch (2 to 4,194,304; the
/// built-in `X_PREFETCH_MIN_KIB` when the variable is not set or holds anything else).
fn tbx_prefetch_min_kib() -> Option<usize> {
    static MIN: OnceLock<Option<usize>> = OnceLock::new();
    *MIN.get_or_init(|| {
        if cfg!(target_arch = "x86_64") {
            crate::utils::xswitch::value(
                "LOSAT_X_TBXPREFETCH",
                2..=X_PREFETCH_MAX_KIB,
                X_PREFETCH_MIN_KIB,
            )
        } else {
            None
        }
    })
}

// No NCBI counterpart: the rule of LOSAT_X_TBXPREFETCH for one subject (module doc), apart from
// the switch so that a test can give it values.
/// Whether the units whose tables have `table_bytes` prefetch at the switch's value `min_kib`: one
/// of them has at least `min_kib` KiB.
fn x_prefetch_units(min_kib: Option<usize>, mut table_bytes: impl Iterator<Item = usize>) -> bool {
    let prefetch = min_kib.is_some_and(|kib| table_bytes.any(|bytes| bytes / 1024 >= kib));
    if min_kib.is_some() {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXPREFETCH", prefetch);
    }
    prefetch
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

/// The summary lines (printed by `main` at exit): the LOSAT_X_STATS counters and the LOSAT_TIMING
/// lines.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    x_kmerge::print_summary();
    // LOSAT_X_SEEDCOMPACT: printed whenever the switch is on, so also when the unit path does not
    // run (zero counters). `bytes_saved` is what the unit outputs no longer
    // request: records x (size of `InitHSP` - size of `XUnitHsp`).
    if seed_compact() && std::env::var_os("LOSAT_X_STATS").is_some() {
        let records = SC_RECORDS.load(AtomicOrdering::Relaxed);
        let saved_per_record =
            (std::mem::size_of::<InitHSP>() - std::mem::size_of::<XUnitHsp>()) as u64;
        eprintln!(
            "[X_STATS] SEEDCOMPACT subjects={} units={} records={} bytes_saved={}",
            SC_SUBJECTS.load(AtomicOrdering::Relaxed),
            SC_UNITS.load(AtomicOrdering::Relaxed),
            records,
            records * saved_per_record
        );
    }
    // LOSAT_X_SEEDUNMAP: printed whenever the switch is on; `active` is 1 where it acts
    // (`seed_unmap_active`). Zero counters at one thread, where no unit path runs, where the rule
    // never holds, and where it does not act.
    if seed_unmap() && std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] SEEDUNMAP active={} diag={} records={} copies={}",
            u8::from(seed_unmap_active()),
            SU_DIAG.load(AtomicOrdering::Relaxed),
            SU_RECORDS.load(AtomicOrdering::Relaxed),
            SU_COPIES.load(AtomicOrdering::Relaxed)
        );
    }
    let load = |c: &AtomicU64| c.load(AtomicOrdering::Relaxed);
    // LOSAT_X_TBXUNITMASK: printed whenever the switch is on and a TBLASTX batch reached the unit
    // path's decision (zero where no unit took a masked table); not for the other programs.
    if let Some(kib) = unit_mask_kib().filter(|_| load(&BATCHES) > 0) {
        if std::env::var_os("LOSAT_X_STATS").is_some() {
            eprintln!(
                "[X_STATS] TBXUNITMASK kib={kib} units_masked={}",
                load(&UM_UNITS)
            );
        }
    }
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] PAIRPAR batches={} subjects={} units={} init_hsps={} fallback_gate={} fallback_guard={} fallback_cost={} band={}",
            load(&BATCHES),
            load(&SUBJECTS),
            load(&UNITS),
            load(&INIT_HSPS),
            load(&FALLBACK_GATE),
            load(&FALLBACK_GUARD),
            load(&FALLBACK_COST),
            load(&BAND)
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
    // The units' seed / extension split (LOSAT_TIMING) or diagonal statistics (LOSAT_X_STATS).
    if let Some(line) = x_unit_probe::timing_line() {
        eprintln!("{line}");
    }
    for line in x_unit_probe::stats_lines() {
        eprintln!("{line}");
    }
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] PAIRLUT batches_direct={} groups={} fallback={}",
            load(&LUT_BATCHES),
            load(&LUT_GROUPS),
            load(&LUT_FALLBACK)
        );
    }
    if timing() {
        // The direct build (query preparation, group tables, stub); it is part of the search's
        // `[TIMING] build_lookup`, and `sublookup` above is then 0 (no `build_groups`).
        eprintln!(
            "[TIMING] x_pairlut: direct {:.3}s (batches={} groups={} fallback={})",
            load(&NS_LUT_DIRECT) as f64 / 1e9,
            load(&LUT_BATCHES),
            load(&LUT_GROUPS),
            load(&LUT_FALLBACK)
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
// use unmasked tables; the masking cannot change a decision (module doc, item 3). This test is the
// condition of the earlier, weaker argument of item 3: every context is at most `slack` long, or no
// subject chunk is longer than `slack`. `decide` applies it only without LOSAT_X_TBXBAND.
/// The first exactness guard of the earlier argument: `max context length <= slack || longest
/// subject chunk <= slack`, with `slack = diag_array_size - query_length` (the reference's own
/// table).
pub(super) fn guard(
    max_context_len: usize,
    longest_subject_chunk: usize,
    query_length: i32,
    diag_array_size: i32,
) -> bool {
    let slack = diag_array_size as i64 - query_length as i64;
    max_context_len as i64 <= slack || longest_subject_chunk as i64 <= slack
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
// ```c
// typedef struct DiagStruct {
//    signed int last_hit   : 31; /**< Offset of the last hit */
//    unsigned int flag      : 1 ; /**< Reset the next extension? */
// } DiagStruct;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:168-173
// ```c
//         if (ewp->diag_table->offset >= INT4_MAX / 4) {
//             ewp->diag_table->offset = ewp->diag_table->window;
//             s_BlastDiagClear(ewp->diag_table);
//         } else {
//             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
//         }
// ```
// The reference keeps `s + offset` in the 31-bit `last_hit`. Its offset is below INT4_MAX/4 before
// this step and so at most INT4_MAX/4 - 1 + subject_length + window after it; with an `s` of up to
// subject_length the stored value stays below 2^30, the first value the signed field cannot hold,
// while `window + 2 * subject_length < 2^29`. Past that the field wraps and the later hits of the
// diagonal are dropped. A unit starts every table at offset = window (module doc, item 1), so it
// equals the reference only where the reference does not wrap. The window is at most `slack`
// (the table has at least query_length + window cells), which keeps this test on the safe side
// without the window; `X_LAST_HIT_MARGIN` covers what a frame adds to a chunk's letters
// (sentinels, the chunk overlap).
/// The second exactness guard: the reference's `last_hit` cannot wrap,
/// `slack + 2 * longest subject chunk + margin < 2^29` with `slack = diag_array_size - query_length`.
pub(super) fn last_hit_fits(
    longest_subject_chunk: usize,
    query_length: i32,
    diag_array_size: i32,
) -> bool {
    let slack = diag_array_size as i64 - query_length as i64;
    slack + 2 * (longest_subject_chunk as i64) + X_LAST_HIT_MARGIN < (1 << 29)
}

/// No NCBI counterpart: room left below 2^29 in `last_hit_fits`.
const X_LAST_HIT_MARGIN: i64 = 1 << 16;

/// No NCBI counterpart: LOSAT_X_PAIRPARGATE keeps the unit path for a batch whose reference
/// diagonal table has at least this many cells (4 bytes each, 1 MiB; a query batch of more than
/// about 65 kb). Measured at one thread (session 13, DESIGN_s13_unitgate.md): from 2^18 cells the
/// small per-group tables of the units make them faster than the reference loop (x0.74-0.99); with
/// 2^17 cells and fewer the units' scan of every subject frame once per group makes them slower
/// (x1.03-1.72).
const X_GATE_DIAG_CELLS: i32 = 1 << 18;

/// No NCBI counterpart: with LOSAT_X_PAIRPARN1, LOSAT_X_PAIRPARGATE keeps the unit path at one
/// thread for a query batch of at least this many letters (`query_length`: the translated frames
/// with their sentinels; a nucleotide query of about 2.5 kb). The constant above was measured
/// before the units' branch-free hit loop and masked tables; measured again (workstation, one
/// thread, line A session 9, 2026-10-25) the units ran x0.82-0.94 of the reference loop for queries
/// of 3-40 kb against one genome (4.6 Mb), the same genome cut into 2-1,000 kb subjects and 3,000
/// zebrafish cDNAs, and x0.906 on 10 fly transcripts against the zebrafish cDNA database (the
/// search that set the constant above); x0.99-1.03 at 2-2.1 kb (about 4,000-4,200 letters) and
/// x1.15-1.96 at 1 kb and 300 bases, where the units' scan of every subject frame once per group
/// outweighs the hits.
const X_GATE_N1_LETTERS: i32 = 5_000;

/// No NCBI counterpart: LOSAT_X_PAIRPARGATE acts only for a lookup threshold of at least TBLASTX's
/// default (13), where it was measured. A lower threshold puts more words in every query position:
/// more hits share each scan position (the units' extra scans weigh less) and the full table grows
/// faster than the group tables (`-threshold +inf`, which NCBI converts to INT_MIN, every word: 14.2
/// GB for the reference loop, 3.3 GB for the units, 5.9 GB for NCBI on the option-sweep inputs).
const X_GATE_MIN_THRESHOLD: i32 = 13;

/// The per-batch inputs of `XPairPar::prepare`.
pub(super) struct XBatch<'a> {
    pub lookup: &'a BlastAaLookupTable,
    pub contexts: &'a [QueryContext],
    pub query_length: i32,
    pub diag_array_size: i32,
    pub subjects: crate::algorithm::tblastx::subjects::Subjects<'a>,
    /// A trace, debug or diagnostics switch is on (they exist only in the reference loop).
    pub debugging: bool,
    /// LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS or the parallel-chunk path is on.
    pub chunk_switches: bool,
    /// PAIRLUT: the groups built at the lookup build site (`lut_build`), or None.
    pub lut: Option<XLutBatch>,
    /// The threads of the batch's search pool (`SearchPool::threads`), as the build site reads them.
    pub threads: usize,
    /// The lookup threshold, as the build site reads it (`XLutEarly::threshold`).
    pub threshold: i32,
}

/// What the unit path does with a batch (`decide`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum XDecision {
    /// A trace, debug, diagnostics or chunk switch: the reference loop ("fallback: gate").
    Gate,
    /// An exactness guard fails (`last_hit_fits`, or `guard` without LOSAT_X_TBXBAND): the
    /// reference loop ("fallback: guard").
    Guard,
    /// LOSAT_X_PAIRPARGATE: the units would be slower than the reference loop, which runs
    /// ("fallback: cost"; `units_gain`).
    Cost,
    /// The units run.
    Units,
}

/// The inputs of `decide`: what the decision reads of a batch. None of them needs the lookup
/// table.
pub(super) struct XGates<'a> {
    pub contexts: &'a [QueryContext],
    pub query_length: i32,
    pub diag_array_size: i32,
    pub subjects: crate::algorithm::tblastx::subjects::Subjects<'a>,
    pub debugging: bool,
    pub chunk_switches: bool,
    /// The threads of the batch's search pool (`SearchPool::threads`).
    pub threads: usize,
    /// `cost_gate()`: LOSAT_X_PAIRPARGATE acts.
    pub cost_gate: bool,
    /// `n1_gate()`: LOSAT_X_PAIRPARN1 acts (one thread: the units from `X_GATE_N1_LETTERS`).
    pub n1_gate: bool,
    /// The lookup threshold (`lookup->threshold`, `ncbi_int4_from_double` of `-threshold`).
    pub threshold: i32,
    /// `tbx_band()`: LOSAT_X_TBXBAND, `guard` is not applied.
    pub band: bool,
    /// `unit_mask_kib().is_some()`: LOSAT_X_TBXUNITMASK, `x_band_fits` is not applied to a batch
    /// whose groups LOSAT_X_TBXQSLICE does not cut.
    pub unit_mask: bool,
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
// The decision of `prepare` (the gates, then the guard against the reference table size above,
// which LOSAT_X_TBXBAND drops, and the second guard, `last_hit_fits`) as one function, so that the
// lookup build site (`lut_build`, before the table exists) and `prepare` (after it) decide with the
// same code.
/// Whether the unit path runs for a batch.
pub(super) fn decide(g: &XGates<'_>) -> XDecision {
    if g.debugging || g.chunk_switches {
        return XDecision::Gate;
    }
    let (max_context_len, longest_subject_chunk) = x_guard_lengths(g.contexts, g.subjects);
    // EXPERIMENT (LOSAT_X_TBXBAND): with the switch the first guard is not applied (module doc
    // item 3) where the subject chunks fit the reference table (`x_band_fits`), or, with
    // LOSAT_X_TBXUNITMASK, where the units' tables are masked or small (`x_band_masked`).
    let band_takes =
        g.band && (x_band_fits(longest_subject_chunk, g.diag_array_size) || x_band_masked(g));
    let first_guard = band_takes
        || guard(
            max_context_len,
            longest_subject_chunk,
            g.query_length,
            g.diag_array_size,
        );
    // The old guard rejects, and the switch's memory rule does not take the batch: the reference
    // loop. The batches the switch takes are marked in `prepare` (`band`).
    if g.band && !first_guard {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXBAND", false);
    }
    if !first_guard || !last_hit_fits(longest_subject_chunk, g.query_length, g.diag_array_size) {
        return XDecision::Guard;
    }
    if g.cost_gate {
        let gain = units_gain(g);
        // The gate's rule sends a batch without gain to the reference loop (false).
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARGATE", gain);
        if !gain {
            return XDecision::Cost;
        }
    }
    XDecision::Units
}

// No NCBI counterpart: a memory rule of LOSAT_X_TBXBAND (module doc). The units need no length
// condition for exactness (item 3), but a unit's unmasked diagonal table has a cell for every raw
// diagonal of its query group and subject chunk, one table per thread, where the reference loop
// has one table of `diag_array_size` cells per subject in flight (LOSAT_X_TBXUNITMASK masks the
// tables of the units this rule would stop, `x_band_masked`). With a subject chunk of at most that
// many letters a unit's table is less than twice the reference's. Measured (workstation, 8 and 16
// threads): a 1 Mb and a 3.9 Mb query against 1 and 4.6 Mb genomes (chunks of 0.16 and 0.18 of the
// table) ran x0.08-0.13 with a lower max RSS; 3 fruit fly transcripts (2,005 bases, 4,096 cells)
// against the fruit fly genome database and a 1.9 kb query against MG1655 (chunks of 2,000 and
// 380 times the table) ran x0.29-0.74, but their max RSS grew from 318 to 600-873 MiB (NCBI 368
// MiB at 8 threads) and from 21 to 68-111 MiB, so such batches keep the reference loop without
// LOSAT_X_TBXUNITMASK.
/// LOSAT_X_TBXBAND: the longest subject chunk has at most as many letters as the reference's
/// diagonal table has cells.
fn x_band_fits(longest_subject_chunk: usize, diag_array_size: i32) -> bool {
    longest_subject_chunk as i64 <= diag_array_size as i64
}

// No NCBI counterpart: the memory rule of LOSAT_X_TBXBAND with LOSAT_X_TBXUNITMASK (module doc).
// A unit whose subject chunk is longer than the reference's table has cells then takes a masked
// table (`x_unit_mask_for`: `x_masked_cells`, at most the reference's), so `x_band_fits` holds for
// every unit that keeps an unmasked one. The scan in slices of LOSAT_X_TBXQSLICE never takes a
// masked table (module doc), so a batch whose groups its rule cuts keeps `x_band_fits`.
// `rule_cuts` reads only the contexts, the thread count and process-wide switches (not the
// diagnostics ones), so the lookup build site and `prepare` agree, with and without diagnostics.
/// LOSAT_X_TBXBAND with LOSAT_X_TBXUNITMASK: the switch is on and LOSAT_X_TBXQSLICE's rule cuts
/// no group of the batch.
fn x_band_masked(g: &XGates<'_>) -> bool {
    g.unit_mask && !x_tbx_qslice::rule_cuts(g.contexts, &group_starts(g.contexts), g.threads)
}

// No NCBI counterpart: the lengths the guards compare, from the batch's contexts and subjects (a
// subject's frame length, at most one chunk).
/// The longest context and the longest subject chunk of a batch, in letters.
fn x_guard_lengths(
    contexts: &[QueryContext],
    subjects: crate::algorithm::tblastx::subjects::Subjects<'_>,
) -> (usize, usize) {
    let max_context_len = contexts.iter().map(|c| c.aa_len).max().unwrap_or(0);
    let longest_subject_chunk = subjects
        .lengths()
        .map(|length| length / 3)
        .max()
        .unwrap_or(0)
        .min(tblastx_max_dbseq_len_for_run());
    (max_context_len, longest_subject_chunk)
}

// No NCBI counterpart: a choice between two LOSAT paths with the same result (module doc items
// 1-4), from the batch's sizes and the thread count only.
/// LOSAT_X_PAIRPARGATE: whether the unit path is expected to be faster than the reference loop for
/// a batch. Either the lookup threshold is below `X_GATE_MIN_THRESHOLD` (not measured; the units
/// stay), or the reference's diagonal table is large (`X_GATE_DIAG_CELLS`), or, at one thread with
/// LOSAT_X_PAIRPARN1, the batch has at least `X_GATE_N1_LETTERS` letters, or the search has several
/// threads and one subject holds more than one and a half thread's share of the bases
/// (`2 * longest * threads > 3 * total`): the reference loop searches a subject on one thread, the
/// units split it over the pool. Several threads without such a subject keep the reference loop
/// below `X_GATE_DIAG_CELLS`: with the subjects spread over the pool the units used less CPU there
/// but lost wall time on many short subjects (x1.02-1.06 at eight threads, session 9 of line A).
fn units_gain(g: &XGates<'_>) -> bool {
    if g.threshold < X_GATE_MIN_THRESHOLD || g.diag_array_size >= X_GATE_DIAG_CELLS {
        return true;
    }
    if g.threads <= 1 {
        if !g.n1_gate {
            return false;
        }
        let gain = g.query_length >= X_GATE_N1_LETTERS;
        // The one-thread rule sends a batch below the line to the reference loop (false).
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARN1", gain);
        return gain;
    }
    let (longest, total) = g
        .subjects
        .lengths()
        .fold((0u128, 0u128), |(longest, total), length| {
            (longest.max(length as u128), total + length as u128)
        });
    2 * longest * g.threads as u128 > 3 * total
}

/// The inputs of `lut_build`: what the lookup build site of the batch has (RI:2311).
pub(super) struct XLutEarly<'a> {
    pub query_frames: &'a [Vec<QueryFrame>],
    /// The lookup threshold (`lookup->threshold`).
    pub threshold: i32,
    pub window: i32,
    pub subjects: crate::algorithm::tblastx::subjects::Subjects<'a>,
    /// As `XBatch::debugging`.
    pub debugging: bool,
    /// `requested_parallel` and `use_serial_scan_chunks` of the batch (RI:2130, 2139).
    pub requested_parallel: bool,
    pub use_serial_scan_chunks: bool,
    pub pool: &'a crate::utils::threading::SearchPool<'a>,
}

/// The direct groups of a batch, made at the lookup build site and handed to `prepare`.
pub(super) struct XLutBatch {
    groups: Vec<XGroup>,
    /// EXPERIMENT (LOSAT_X_TBXQSLICE): the slices of the batch's groups when one is cut
    /// (x_tbx_qslice.rs); `groups` is then empty.
    sliced: Option<Vec<XSliced>>,
    /// The stub of the full table and the batch's contexts, until the build site takes them
    /// (`take_stub`).
    stub: Option<(BlastAaLookupTable, Vec<QueryContext>)>,
}

impl XLutBatch {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:264-266
    // ```c
    //    case eAaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastAaLookupTable*)lookup->lut)->longest_chain;
    // ```
    // The reference code reads `longest_chain` and the header of the batch table; the stub carries
    // the full table's values (`XPairLutQuery::x_into_stub`).
    /// The table and contexts that take the place of `build_ncbi_lookup`'s for this batch: the
    /// stub.
    pub(super) fn take_stub(&mut self) -> Option<(BlastAaLookupTable, Vec<QueryContext>)> {
        self.stub.take()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-250
// ```c
//     if (backup->offset + MAX_DBSEQ_LEN <
//         backup->hard_ranges[backup->hm_index].right) {
//
//         subject->length = MAX_DBSEQ_LEN;
//         backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
// ```
// `chunk_switches` of the batch before the lookup build: the expression of RI:2385-2394
// (`use_parallel_chunks`, which reads the thread request, two switches and the subject lengths
// only) or'ed with `use_serial_scan_chunks`, as RI:2562 passes it to `prepare`. `prepare` decides
// again with RI's own value and stops the run if the two decisions differ for a batch with direct
// groups.
fn early_chunk_switches(
    requested_parallel: bool,
    use_serial_scan_chunks: bool,
    subjects: crate::algorithm::tblastx::subjects::Subjects<'_>,
) -> bool {
    let use_parallel_chunks = requested_parallel
        && !use_serial_scan_chunks
        && (std::env::var_os("LOSAT_TBLASTX_PARALLEL_CHUNKS").is_some()
            || (cfg!(all(target_arch = "wasm32", feature = "wasm-threads"))
                && subjects.lengths().any(|length| {
                    tblastx_estimated_subject_chunk_work_items(
                        length,
                        tblastx_max_dbseq_len_for_run(),
                    ) > 1
                })));
    use_serial_scan_chunks || use_parallel_chunks
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-97
// ```c
// BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
//                       &lookup_wrap->lut);
// ((BlastAaLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
// BlastAaLookupIndexQuery( (BlastAaLookupTable*) lookup_wrap->lut, matrix,
//                           query, lookup_segments, 0);
// ```
// Dispatch point 0 (per batch, at the lookup build): when the unit path will run for the batch,
// the tables its units scan with are built here, one per context group, from the reference query
// preparation; the full table NCBI builds above is then not built. None leaves the batch to the
// reference builder and `build_groups`: the direct builder's preconditions
// (`XPairLutQuery::x_prepare`), a gate or the guard (`decide`), or a lookup location outside one
// group.
/// The direct groups of the batch, or None (the reference builder runs).
pub(super) fn lut_build(e: XLutEarly<'_>) -> Option<XLutBatch> {
    let t0 = timing().then(Instant::now);
    let built = lut_direct(&e);
    add_ns(&NS_LUT_DIRECT, t0);
    crate::utils::xreach::x_reach!("LOSAT_X_PAIRLUT", built.is_some());
    match &built {
        Some(batch) => {
            LUT_BATCHES.fetch_add(1, AtomicOrdering::Relaxed);
            LUT_GROUPS.fetch_add(batch.groups.len() as u64, AtomicOrdering::Relaxed);
        }
        None => {
            LUT_FALLBACK.fetch_add(1, AtomicOrdering::Relaxed);
        }
    }
    built
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-57
// ```c
// diag_array_length = 1;
// /* What power of 2 is just longer than the query? */
// while (diag_array_length < (qlen+window_size))
// {
//         diag_array_length = diag_array_length << 1;
// }
// ```
// The early decision: the inputs of `decide` from the prepared query, by the expressions RI uses
// after the lookup build (RI:2502-2510 `query_length`, `check_diag_table_window`,
// `diag_array_size`), then the groups and the stub.
fn lut_direct(e: &XLutEarly<'_>) -> Option<XLutBatch> {
    // The reference query preparation (contexts, concatenated query, lookup locations).
    let q = XPairLutQuery::x_prepare(e.query_frames, e.threshold, true)?;
    let contexts = q.contexts();
    let query_length: i32 = contexts
        .last()
        .map(|c| c.frame_base + c.aa_len as i32)
        .unwrap_or(0);
    // RI:2506 reports this error after the reference build.
    if crate::blastinput::app::check_diag_table_window(query_length, e.window, "TBLASTX").is_err() {
        return None;
    }
    let mut diag_array_size: i32 = 1;
    while diag_array_size < (query_length + e.window) {
        diag_array_size <<= 1;
    }
    let gates = XGates {
        contexts,
        query_length,
        diag_array_size,
        subjects: e.subjects,
        debugging: e.debugging,
        chunk_switches: early_chunk_switches(
            e.requested_parallel,
            e.use_serial_scan_chunks,
            e.subjects,
        ),
        threads: e.pool.threads(),
        cost_gate: cost_gate(),
        n1_gate: n1_gate(),
        threshold: e.threshold,
        band: tbx_band(),
        unit_mask: unit_mask_kib().is_some(),
    };
    if decide(&gates) != XDecision::Units {
        return None;
    }
    // EXPERIMENT (LOSAT_X_TBXQSLICE): the tables of the slices of the groups the switch cuts,
    // built directly as the groups are (x_tbx_qslice.rs).
    let starts = group_starts(contexts);
    let (groups, sliced) = match x_tbx_qslice::partition(contexts, &starts, e.pool.threads()) {
        None => (lut_direct_groups(&q, e.pool, lut_on_pool())?, None),
        Some(parts) => {
            let sliced =
                x_tbx_qslice::lut_direct_sliced(&q, &starts, &parts, e.pool, lut_on_pool())?;
            (Vec::new(), Some(sliced))
        }
    };
    // The stub's `longest_chain` sums every cell over a partition of the entries: the whole
    // groups, or the slices (x_tbx_qslice.rs, item 5).
    let stub = match &sliced {
        None => q.x_into_stub(groups.iter().map(|g| &g.lookup)),
        Some(sliced) => q.x_into_stub(x_tbx_qslice::tables(sliced)),
    };
    Some(XLutBatch {
        groups,
        sliced,
        stub: Some(stub),
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:454-456,461-465
// ```c
//     BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
//                                       lookup->charsize, lookup->word_length,
//                                       query, location);
// ...
//     for (i = 0; i < lookup->backbone_size; i++) {
//         if (exact_backbone[i] != NULL) {
//             s_AddWordHits(lookup, matrix, query->sequence,
//                           exact_backbone[i], query_bias, row_max);
//             sfree(exact_backbone[i]);
// ```
// One table per context group (`group_starts`, the groups of `build_groups`), each from the
// lookup locations of its contexts (`XPairLutQuery::x_build_group`). Scheduling only: the groups
// are independent. They are built in order on the calling thread; with `on_pool` (the tuning knob
// LOSAT_X_PAIRLUT_POOL, `lut_on_pool`) and a search pool they are built on the pool (an indexed
// map, results in group order; each thread has its own neighbour-list memo). The calling thread
// is the default because the pool gained at most 17 ms per batch (AP027131, 6 x 220,702 aa) and
// cost memory: the memo of every pool thread, and tables allocated in the pool threads' arenas
// and freed by the search thread (+17 MiB on a 50-batch run at 8 threads).
/// The groups `build_groups` would filter out of the full table of `q`, built directly. None
/// when a lookup location does not lie inside one group.
fn lut_direct_groups(
    q: &XPairLutQuery,
    pool: &crate::utils::threading::SearchPool<'_>,
    on_pool: bool,
) -> Option<Vec<XGroup>> {
    let starts = group_starts(q.contexts());
    let locations = q.x_group_locations(&starts)?;
    Some(lut_direct_tables(q, &locations, pool, on_pool))
}

// No NCBI counterpart: scheduling of the builds (the NCBI reference above `lut_direct_groups`;
// each table is built by `XPairLutQuery::x_build_group`).
/// One table per list of lookup locations (`XPairLutQuery::x_build_group`), in list order: the
/// groups of `lut_direct_groups`, or the slices of LOSAT_X_TBXQSLICE (x_tbx_qslice.rs). On the
/// search pool with `on_pool` (LOSAT_X_PAIRLUT_POOL), as `lut_direct_groups` says.
fn lut_direct_tables(
    q: &XPairLutQuery,
    locations: &[Vec<(i32, i32)>],
    pool: &crate::utils::threading::SearchPool<'_>,
    on_pool: bool,
) -> Vec<XGroup> {
    let to_group = |g: XPairLutGroup| XGroup {
        scan_bf: crate::algorithm::tblastx::x_scan_bf::x_scan_bf_for(&g.table),
        lookup: g.table,
        q_min: g.q_min,
        q_max: g.q_max,
    };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if on_pool && pool.enabled() {
        return pool.install(|| {
            locations
                .par_iter()
                .map(|group_locations| to_group(q.x_build_group(group_locations)))
                .collect()
        });
    }
    let _ = (pool, on_pool);
    locations
        .iter()
        .map(|group_locations| to_group(q.x_build_group(group_locations)))
        .collect()
}

// No NCBI counterpart: a tuning knob for where the groups are built (the calling thread or the
// search pool); the tables are the same either way.
/// LOSAT_X_PAIRLUT_POOL: the groups on the search pool instead of the calling thread.
fn lut_on_pool() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_PAIRLUT_POOL").is_some())
}

// No NCBI counterpart: NCBI's word finder hands its init-HSP list straight to
// `BLAST_GetUngappedHSPList` (blast_engine.c:561-562). Here the lists of a subject's units are
// carried from `XPairPar::run_subject_units` to `XPairPar::merge_subject` as a value, so that
// nothing of the group tables is needed between the two.
/// The units of one subject after they ran: the chunks of each frame, the units in their order,
/// and each unit's saved init HSPs. It owns its lists and holds no reference to a group table.
pub(super) struct XSubjectUnits {
    frame_chunks: Vec<Vec<SubjectChunk>>,
    units: Vec<XUnit>,
    outputs: XUnitOutputs,
}

// No NCBI counterpart: the two forms of the units' saved init HSPs (section above `XUnitHsp`),
// chosen once per subject by `run_subject_units`; `merge_subject` reads the form once per subject.
/// The saved init HSPs of each unit, in unit order.
enum XUnitOutputs {
    /// Switch off: today's lists.
    Full(Vec<Vec<InitHSP>>),
    /// LOSAT_X_SEEDCOMPACT: the compact records, each list at its exact length.
    Compact(Vec<Vec<XUnitHsp>>),
}

/// The per-subject inputs of `XPairPar::run_subject_units` and `XPairPar::merge_subject` (the
/// values the reference loop reads).
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
    /// EXPERIMENT (LOSAT_X_TBXSCANBF): `x_scan_bf_for(&lookup)`, decided once per table.
    scan_bf: bool,
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
    /// LOSAT_X_SEEDCOMPACT: the records of the unit that runs (`run_unit_compact`), cleared at
    /// the start of every unit and kept across units; empty while the switch is off.
    records: Vec<XUnitHsp>,
    /// LOSAT_X_SEEDUNMAP as it acts in this process (`seed_unmap_active`) when the scratch was made
    /// (`x_unit_mapped`); always false for the serial path of `run_units`.
    unmap: bool,
    /// LOSAT_TIMING: the clocks of the unit that runs (x_unit_probe.rs).
    time_probe: XTimeProbe,
    /// LOSAT_X_STATS: the diagonal statistics of the unit that runs (x_unit_probe.rs); empty
    /// without the variable.
    stats_probe: XStatsProbe,
}

impl XScratch {
    /// An empty scratch; `unmap` is LOSAT_X_SEEDUNMAP as it acts here on the search pool
    /// (`seed_unmap_active`), false on the serial path, or a test's choice.
    fn new(unmap: bool) -> Self {
        Self {
            unmap,
            ..Self::default()
        }
    }
}

// ---------------------------------------------------------------------------------------------
// EXPERIMENT (LOSAT_X_SEEDUNMAP): the unit path's large per-split blocks mapped by glibc.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
// ```c
//         diag_table->hit_level_array = (DiagStruct *)
//             calloc(diag_table->diag_array_length, sizeof(DiagStruct));
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:336-341
// ```c
//     if (num >= num_avail) {
//         if (init_hitlist->do_not_reallocate)
//             return FALSE;
//         num_avail *= 2;
//         match_array = (BlastInitHSP *)
//             realloc(match_array, num_avail * sizeof(BlastInitHSP));
// ```
// No NCBI counterpart for the capacities: NCBI allocates one diagonal table and one init hit list
// per search thread (`BlastExtendWordNew`, `BLAST_InitHitListNew`, blast_engine.c:1043) and keeps
// them for the thread's subjects. The unit path holds, per rayon split of a subject's units
// (`run_units`, `map_init`), a unit diagonal table (`XScratch::diag`) and a record vector
// (`XScratch::records`, LOSAT_X_SEEDCOMPACT), and per unit a copy of its records that lives until
// the subject's merge. glibc maps a block of at least its mmap threshold and unmaps it when it is
// freed, but it raises the threshold to the size of each mapped block that is freed (up to
// 32 MiB; in a TBLASTX search SEG's buffer and the pair lookup's group tables are such blocks) and
// its trim threshold to twice that. These blocks then come from the worker threads' arenas, and
// when they are freed the arenas keep their freed tops through the linking call, where the search
// peaks (AP027131 x AP027133, genetic code 4, 8 threads: 10-14 % of the max RSS). With the switch,
// for a unit where `x_unit_mapped` holds, the three blocks get a capacity above glibc's largest
// dynamic threshold (`X_MAPPED_BYTES`), so that glibc maps them and unmaps them when they are
// dropped (at the end of the split, and in the merge); the pages beyond their length are never
// touched and are not resident. Other allocators return large blocks on their own.
//
// Where mapping costs more than it returns, three rules:
// 1. Long on both sides. A mapped table is new zeroed pages at every split, about one page fault
//    per 1,024 cells of its width (query-group span + chunk), where the arena handed the same block
//    back; the unit's scan work grows with span x chunk, so the cost per work is about
//    1 / span + 1 / chunk and small only when both are long. Measured at 8 threads with a rule that
//    mapped every unit whose chunk was at least a quarter of its table: 200 queries of about 3 kb
//    from AP027131 against AP027133 (genetic code 4; 50 query batches, group spans of about 4,000
//    letters, chunks of 202,000) mapped 1,748 tables, sys 0.05-0.2 -> 0.4-0.7 s and x1.05-1.07 of
//    the wall time; MelaMJNV against AP027133 in six pieces (span about 96,000, chunks about
//    33,000) mapped 159 tables, sys 0.034 -> 0.070 s; AP027131 against 200 contigs of about 3 kb
//    x1.07 with every table mapped. So a unit is mapped only when its query-group span and its
//    chunk are both at least `X_MAPPED_MIN_SIDE` letters (the table is then at least 512 KiB, above
//    `X_MAPPED_MIN`); AP027131 x AP027133 (span about 220,000, chunk 202,000), NZ_CP006932 against
//    itself (219,000 / 219,000) and MelaMJNV x PemoMJNVA (96,000 / 98,000) are mapped, the three
//    searches above are not. A copy is mapped only when it has at least `X_MAPPED_MIN` bytes.
// 2. Transparent huge pages that round mapped blocks up. Where the kernel gives anonymous mappings
//    huge pages by default (`/sys/kernel/mm/transparent_hugepage/enabled` "[always]": WSL2's
//    kernel, the RHEL family's default), the 32 MiB regions get 2 MiB pages on first touch, so a
//    mapped block of 0.1-0.5 MiB is resident as 2 MiB (8 threads, 5 pairs: the 200 queries above
//    53-56 -> 60-67 MiB of max RSS, AP027131 x AP027133 137-140 -> 124-145, higher in 4 of 5
//    pairs; the same binary with THP disabled for the process by prctl PR_SET_THP_DISABLE: 52 ->
//    44, 135 -> 122 and NZ_CP006932 185 -> 168 MiB, lower in 5 of 5). So the switch acts only where
//    that file selects "[madvise]" or "[never]", or `/proc/self/status` says `THP_enabled: 0`
//    (`x_thp_rounds_up`, read once per process with `std::fs`, `seed_unmap_active`); with the file
//    unreadable it does nothing.
// 3. One thread. `run_units` then takes its serial path on the calling thread, whose blocks come
//    from glibc's main arena, which glibc trims: nothing is retained there (AP027131 x AP027133 at
//    one thread: 78 MiB of max RSS with and without MALLOC_ARENA_MAX or a fixed threshold), and
//    mapping only costs (77 -> 78 MiB and 2,800 more page faults with the blocks mapped). So only
//    the scratches of the search pool (`pool.enabled()`, two or more threads) carry the switch; the
//    serial scratch keeps today's capacities.
//
// Why no value changes: (1) `scan_unit` still clears the table and resizes it to the unit's width,
// writing every cell it reads, and reads no capacity; (2) the record vector is still cleared at the
// start of the unit and filled by the same pushes; (3) the copy is the records in their order
// (`extend_from_slice` instead of `to_vec`), and the merge (`merge_subject_compact`,
// `x_expand_unit`) reads the lengths and the elements of the copies, never their capacity.
// ---------------------------------------------------------------------------------------------

/// No NCBI counterpart: glibc's largest dynamic mmap threshold (DEFAULT_MMAP_THRESHOLD_MAX on
/// 64-bit) and one page more; a block of at least this capacity is mapped.
const X_MAPPED_BYTES: usize = (32 << 20) + 4096;

/// No NCBI counterpart: glibc's default mmap threshold; a smaller copy is left to the arenas.
const X_MAPPED_MIN: usize = 128 << 10;

/// No NCBI counterpart: a unit's blocks are mapped only when its query-group span and its subject
/// chunk both have at least this many letters (section above, rule 1).
const X_MAPPED_MIN_SIDE: usize = 1 << 16;

/// No NCBI counterpart: the capacity of a mapped block of `len` elements of `T`.
fn x_mapped_len<T>(len: usize) -> usize {
    len.max(X_MAPPED_BYTES / std::mem::size_of::<T>() + 1)
}

/// No NCBI counterpart: the rule of LOSAT_X_SEEDUNMAP for a unit whose query group spans `span`
/// letters (`x_unit_span`) and whose subject chunk has `chunk_len` letters: the switch acts (`on`,
/// `XScratch::unmap`), the target is Linux with glibc, and both sides have at least
/// `X_MAPPED_MIN_SIDE` letters (section above, rule 1; the table is then at least 512 KiB, so the
/// table's own minimum of `X_MAPPED_MIN` is implied and not tested).
fn x_unit_mapped(on: bool, span: usize, chunk_len: usize) -> bool {
    on && cfg!(all(target_os = "linux", target_env = "gnu"))
        && span >= X_MAPPED_MIN_SIDE
        && chunk_len >= X_MAPPED_MIN_SIDE
}

/// No NCBI counterpart (LOSAT_X_SEEDUNMAP, section above `X_MAPPED_BYTES`): a unit's copy of its
/// records, exact (`to_vec`), or, for a unit where `x_unit_mapped` holds (`mapped`) whose records
/// have at least `X_MAPPED_MIN` bytes, the records in their order in a vector of the mapped capacity.
fn x_unit_copy(records: &[XUnitHsp], mapped: bool) -> Vec<XUnitHsp> {
    if mapped && std::mem::size_of_val(records) >= X_MAPPED_MIN {
        let mut copy = Vec::with_capacity(x_mapped_len::<XUnitHsp>(records.len()));
        copy.extend_from_slice(records);
        SU_COPIES.fetch_add(1, AtomicOrdering::Relaxed);
        copy
    } else {
        records.to_vec()
    }
}

/// No NCBI counterpart (LOSAT_X_SEEDUNMAP, section above): `x_unit_mapped` for a unit whose query
/// group spans `span` letters; where it holds, `records` (empty) gets the mapped capacity unless it
/// has that much room. The statements of `XPairPar::x_unit_records`, for LOSAT_X_TBXQSLICE's units
/// (x_tbx_qslice.rs), whose groups have no table in `groups`. `x_unit_records` keeps its own copy:
/// calling this function from it changed the machine code of the default unit scan (one-thread
/// instruction counts -2.0 % to +0.1 %, v4-A DECISIONS 19).
fn x_records_for(
    span: usize,
    unit: &XUnit,
    scratch: &XScratch,
    records: &mut Vec<XUnitHsp>,
) -> bool {
    let mapped = x_unit_mapped(scratch.unmap, span, unit.chunk.length);
    if scratch.unmap {
        crate::utils::xreach::x_reach!("LOSAT_X_SEEDUNMAP", mapped);
    }
    if mapped && records.capacity() * std::mem::size_of::<XUnitHsp>() < X_MAPPED_BYTES {
        *records = Vec::with_capacity(x_mapped_len::<XUnitHsp>(0));
        SU_RECORDS.fetch_add(1, AtomicOrdering::Relaxed);
    }
    mapped
}

/// No NCBI counterpart: whether transparent huge pages round the mapped blocks of LOSAT_X_SEEDUNMAP
/// up (section above, rule 2), from the text of `/sys/kernel/mm/transparent_hugepage/enabled`
/// (`enabled`, the selected value in brackets) and of `/proc/self/status` (`status`). False only
/// when the process has huge pages disabled (`THP_enabled: 0`, prctl PR_SET_THP_DISABLE) or the
/// kernel gives them only on request ("[madvise]") or never ("[never]"); true for "[always]", an
/// unreadable file or text it does not recognise.
fn x_thp_rounds_up(enabled: Option<&str>, status: Option<&str>) -> bool {
    let disabled_here = status.is_some_and(|text| {
        text.lines().any(|line| {
            line.strip_prefix("THP_enabled:")
                .is_some_and(|value| value.trim() == "0")
        })
    });
    let selected = enabled.and_then(|text| {
        let start = text.find('[')? + 1;
        let end = start + text[start..].find(']')?;
        Some(&text[start..end])
    });
    !(disabled_here || matches!(selected, Some("madvise") | Some("never")))
}

// No NCBI counterpart: the unit cells of `scan_unit` (the cell of hit (q, s) is q - s - lo,
// lo = q_min - (s_len - 1)), so they run from 0 to (q_max - q_min) + (s_len - 1). One pair of
// helpers for the table of `scan_unit` and for LOSAT_X_SEEDUNMAP's rule (in `scan_unit` and
// `run_unit_compact`), so that they agree.
/// The query-group span of a unit of `group`: its lookup entries' highest query offset minus the
/// lowest.
fn x_unit_span(group: &XGroup) -> usize {
    (group.q_max - group.q_min) as usize
}

/// The number of cells of the unit diagonal table of a unit of `group` on `chunk`.
fn x_unit_width(group: &XGroup, chunk: &SubjectChunk) -> usize {
    x_unit_span(group) + chunk.length
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
    /// EXPERIMENT (LOSAT_X_TBXUNGAPSENT): with `SENT` (a unit `x_unit_sent` accepted: its chunk is
    /// the whole frame), the frame's buffer `[0, subject.., 0]`; else empty.
    subject_sent: &'a [u8],
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

// ---------------------------------------------------------------------------------------------
// EXPERIMENT (LOSAT_X_SEEDCOMPACT): compact, exactly sized seed-stage lists of the unit path.
//
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:142-155
// ```c
// typedef struct BlastUngappedData {
//    Int4 q_start; /**< Start of the ungapped alignment in query */
//    Int4 s_start; /**< Start of the ungapped alignment in subject */
//    Int4 length;  /**< Length of the ungapped alignment */
//    Int4 score;   /**< Score of the ungapped alignment */
// } BlastUngappedData;
//
// /** Structure to hold the initial HSP information */
// typedef struct BlastInitHSP {
//     BlastOffsetPair offsets; /**< Offsets in query and subject, or, in PHI
//                                 BLAST, start and end of pattern in subject. */
//     BlastUngappedData* ungapped_data; /**< Pointer to a structure holding
//                                          ungapped alignment information */
// } BlastInitHSP;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:360-372
// ```c
// BlastSaveInitHsp(BlastInitHitList * ungapped_hsps, Int4 q_start, Int4 s_start,
//                  Int4 q_off, Int4 s_off, Int4 len, Int4 score)
// {
//     BlastUngappedData *ungapped_data = NULL;
//
//     ungapped_data = (BlastUngappedData *) malloc(sizeof(BlastUngappedData));
//
//     ungapped_data->q_start = q_start;
//     ungapped_data->s_start = s_start;
//     ungapped_data->length = len;
//     ungapped_data->score = score;
//
//     BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2371-2375
// ```c
// s_GetUngappedHSPContext(const BlastQueryInfo* query_info,
//                         const BlastInitHSP* init_hsp)
// {
//     return BSearchContextInfo(init_hsp->offsets.qs_offsets.q_off, query_info);
// }
// ```
// NCBI saves an initial HSP as six values: the query start in the concatenated query, the subject
// start, the length, the score and the two seed offsets. Its context is searched again from the
// query offset when the list is converted, and the frames are those of the call. `XUnitHsp` is that
// record plus the context index `x_hit` has already found: 28 bytes, against the 72 of an `InitHSP`
// on 64-bit targets, whose other fields are the two ends (start plus length) and seven values that
// are the same for every init HSP of a unit or of a context.
//
// What the switch changes (read once per process, `seed_compact`; not in the strict set):
// 1. `x_hit` pushes into an `XInitSink`. For `Vec<InitHSP>` (switch off) that is today's push. For
//    `Vec<XUnitHsp>` the sink keeps `q_start_absolute`, `s_start`, `s_end - s_start`,
//    `q_seed_absolute`, `s_seed`, `score` and `ctx_idx` of the `InitHSP` that `x_hit` built.
// 2. The records of a unit go into its thread's scratch vector (`XScratch::records`, cleared at the
//    start of every unit and kept across units, as NCBI keeps one init hit list per thread and
//    resets it per chunk, blast_engine.c:491,1043), and the unit's output is an exact copy of the
//    records (`to_vec`, `run_unit_compact`; with LOSAT_X_SEEDUNMAP a large copy has the mapped
//    capacity, section above `X_MAPPED_BYTES`, and the same elements).
// 3. `merge_subject_compact` fills each chunk's init list by expanding the records of the chunk's
//    units (`x_expand_unit`), converts it with `get_ungapped_hsp_list_into` into a Vec reserved at
//    the list's length, and a chunk list that meets an empty frame list becomes the frame list.
//
// Why every chunk's init list, and so `combined_ungapped_hits`, is today's:
// a. The expansion gives every field of the `InitHSP` that `x_hit` pushes. `q_start_absolute`,
//    `s_start`, `q_seed_absolute`, `s_seed`, `score` and `ctx_idx` are kept. `x_hit` sets
//    `s_end = hsp_s + hsp_len`, `q_start_absolute = ctx.frame_base + hsp_q` and
//    `q_end_absolute = ctx.frame_base + (hsp_q + hsp_len)`, so the record's `len` is `hsp_len`,
//    `s_start + len` is `s_end` and `q_start + len` is `q_end_absolute` (the terms are not negative,
//    so one sum overflows exactly when the other does). `scan_unit` builds the `XHitCtx` with
//    `s_f_idx: unit.f`, `s_frame: inp.s_frames_preliminary[unit.f].frame`, `s_idx: inp.s_idx`,
//    `s_len: inp.s_len` and `contexts: inp.contexts`, and `x_hit` takes `q_idx`, `q_frame` and
//    `q_orig_len` from `ctx = contexts[ctx_idx]`; the expansion takes the same values from the
//    unit's `f`, the same `XSubject` (dispatch point 3 of run_impl.rs passes one `x_subject` to
//    `run_subject_units` and `merge_subject`) and `inp.contexts[ctx_idx]`. `ctx_idx` is below the
//    number of contexts; `u32::try_from` stops the run rather than truncate it.
// b. The expansion maps a unit's records one to one, in their order, and the chunk loop takes the
//    units of a chunk in unit order (today's loop), so each chunk's init list is today's element for
//    element.
// c. No element changes after that: the exact copy holds the scratch's records in their order; the
//    conversion's body reads the length and the elements of the Vec it fills, never its capacity
//    (the Vec starts empty: `get_ungapped_hsp_list_into` asserts it); and the move gives the empty
//    frame list the chunk list's elements in their order, which is what the empty-list branch of
//    `merge_tblastx_subject_chunk_hits` computes (`combined.append(&mut incoming)`).
// d. So the sort, conversion, offset adjustment, chunk merge and frame append and sort get today's
//    inputs and run today's statements: `combined_ungapped_hits` is today's list. The unit test
//    `compact_check` compares the compact merge with the full-list merge (both are the unit path);
//    unit path against frame loop is checked through the frozen NCBI bytes at 1 and 4 threads.
// ---------------------------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:142-147
// ```c
// typedef struct BlastUngappedData {
//    Int4 q_start; /**< Start of the ungapped alignment in query */
//    Int4 s_start; /**< Start of the ungapped alignment in subject */
//    Int4 length;  /**< Length of the ungapped alignment */
//    Int4 score;   /**< Score of the ungapped alignment */
// } BlastUngappedData;
// ```
// The values of a saved init HSP that vary inside a unit: NCBI's `BlastUngappedData` and
// `BlastInitHSP.offsets`, and the context `BSearchContextInfo` gives for the query offset.
/// One saved init HSP of a unit (LOSAT_X_SEEDCOMPACT): 28 bytes.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct XUnitHsp {
    /// `InitHSP::q_start_absolute` (NCBI `ungapped_data->q_start`, concatenated query).
    q_start: i32,
    /// `InitHSP::s_start` (NCBI `ungapped_data->s_start`).
    s_start: i32,
    /// `InitHSP::s_end - InitHSP::s_start` (NCBI `ungapped_data->length`).
    len: i32,
    /// `InitHSP::q_seed_absolute` (NCBI `offsets.qs_offsets.q_off`).
    q_seed: i32,
    /// `InitHSP::s_seed` (NCBI `offsets.qs_offsets.s_off`).
    s_seed: i32,
    /// `InitHSP::score` (NCBI `ungapped_data->score`).
    score: i32,
    /// `InitHSP::ctx_idx`.
    ctx_idx: u32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:360-372
// ```c
// BlastSaveInitHsp(BlastInitHitList * ungapped_hsps, Int4 q_start, Int4 s_start,
//                  Int4 q_off, Int4 s_off, Int4 len, Int4 score)
// ...
//     BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
// ```
// Where `x_hit` saves an init HSP (`BlastSaveInitHsp`). `x_hit` and the scan of a unit
// (`scan_unit`) are generic over the sink, so each is compiled once per sink and no test of the
// switch runs per hit; the push statement of `x_hit` is today's.
/// The list a unit saves its init HSPs into.
trait XInitSink {
    /// Saves one init HSP as `x_hit` built it.
    fn push(&mut self, hsp: InitHSP);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:351-355
// ```c
//     match_array[num].offsets.qs_offsets.q_off = q_off;
//     match_array[num].offsets.qs_offsets.s_off = s_off;
//     match_array[num].ungapped_data = ungapped_data;
//
//     init_hitlist->total++;
// ```
// Switch off: today's list, today's push.
impl XInitSink for Vec<InitHSP> {
    #[inline(always)]
    fn push(&mut self, hsp: InitHSP) {
        Vec::push(self, hsp);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:367-372
// ```c
//     ungapped_data->q_start = q_start;
//     ungapped_data->s_start = s_start;
//     ungapped_data->length = len;
//     ungapped_data->score = score;
//
//     BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
// ```
// LOSAT_X_SEEDCOMPACT: the record of the section above (item 1). `len` is `hsp_len` of `x_hit`
// (`s_end: hsp_s + hsp_len`). After inlining, the fields of `hsp` that the record drops are not
// computed.
impl XInitSink for Vec<XUnitHsp> {
    #[inline(always)]
    fn push(&mut self, hsp: InitHSP) {
        Vec::push(
            self,
            XUnitHsp {
                q_start: hsp.q_start_absolute,
                s_start: hsp.s_start,
                len: hsp.s_end - hsp.s_start,
                q_seed: hsp.q_seed_absolute,
                s_seed: hsp.s_seed,
                score: hsp.score,
                ctx_idx: u32::try_from(hsp.ctx_idx)
                    .expect("LOSAT_X_SEEDCOMPACT: a context index fits 32 bits"),
            },
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4756-4767
// ```c
//       context = s_GetUngappedHSPContext(query_info, init_hsp);
//       s_AdjustInitialHSPOffsets(init_hsp,
//                                 query_info->contexts[context].query_offset);
//       ungapped_data = init_hsp->ungapped_data;
//       Blast_HSPInit(ungapped_data->q_start,
//                     ungapped_data->length+ungapped_data->q_start,
//                     ungapped_data->s_start,
//                     ungapped_data->length+ungapped_data->s_start,
//                     init_hsp->offsets.qs_offsets.q_off,
//                     init_hsp->offsets.qs_offsets.s_off,
//                     context, query_info->contexts[context].frame,
//                     subject->frame, ungapped_data->score, NULL, &new_hsp);
// ```
// LOSAT_X_SEEDCOMPACT: the `InitHSP`s of one unit's records, appended to a chunk's init list in
// the records' order, with the values `x_hit` takes them from (section above, item a): the ends
// as NCBI's conversion computes them (`length + start`), the subject frame of the unit's frame `f`
// (`XUnit::f`), the context's fields from `inp.contexts[ctx_idx]`, the subject's index and length
// from `inp`.
/// Appends the `InitHSP`s of the records of a unit of subject frame `f` to `init_hsps`.
fn x_expand_unit(records: &[XUnitHsp], f: usize, inp: &XSubject<'_>, init_hsps: &mut Vec<InitHSP>) {
    let unit = XExpand::new(f, inp);
    init_hsps.extend(records.iter().map(|r| unit.init_hsp(r)));
}

/// The values a record of a unit of subject frame `f` is expanded with (`x_expand_unit`).
#[derive(Clone, Copy)]
struct XExpand<'a, 'b> {
    f: usize,
    /// `XHitCtx::s_frame` of the unit (`scan_unit`).
    s_frame: i8,
    /// `XHitCtx::s_idx` of the unit (`scan_unit`).
    s_idx: u32,
    inp: &'a XSubject<'b>,
}

impl<'a, 'b> XExpand<'a, 'b> {
    fn new(f: usize, inp: &'a XSubject<'b>) -> Self {
        Self {
            f,
            s_frame: inp.s_frames_preliminary[f].frame,
            s_idx: inp.s_idx as u32,
            inp,
        }
    }

    /// The `InitHSP` of one record.
    #[inline]
    fn init_hsp(&self, r: &XUnitHsp) -> InitHSP {
        let ctx_idx = r.ctx_idx as usize;
        let ctx = &self.inp.contexts[ctx_idx];
        InitHSP {
            q_start_absolute: r.q_start,
            q_end_absolute: r.q_start + r.len,
            s_start: r.s_start,
            s_end: r.s_start + r.len,
            q_seed_absolute: r.q_seed,
            s_seed: r.s_seed,
            score: r.score,
            ctx_idx,
            s_f_idx: self.f,
            q_idx: ctx.q_idx,
            s_idx: self.s_idx,
            q_frame: ctx.frame,
            s_frame: self.s_frame,
            q_orig_len: ctx.orig_len,
            s_orig_len: self.inp.s_len,
        }
    }
}

/// The unit path of one batch: the sub-lookups of its context groups.
pub(super) struct XPairPar {
    /// The groups' tables. Empty when LOSAT_X_TBXQSLICE cuts a group of the batch (`sliced`).
    groups: Vec<XGroup>,
    /// PAIRLUT: the groups were built directly and the batch's table is a stub without
    /// entries, which the reference frame loop must not scan.
    stub: bool,
    /// EXPERIMENT (LOSAT_X_TBXQSLICE): the batch's groups as slices, when the switch cuts one of
    /// them (x_tbx_qslice.rs); the units then index this list and scan with its tables.
    sliced: Option<Vec<XSliced>>,
    /// EXPERIMENT (LOSAT_X_TBXUNITMASK): the batch's rule, None when the switch is off or
    /// LOSAT_X_TBXQSLICE's rule cuts a group of the batch (module doc).
    mask_rule: Option<XMaskRule>,
    /// LOSAT_X_TBXUNITMASK: each group's extent (`x_group_extents`) while `mask_rule` is Some,
    /// else empty.
    extents: Vec<usize>,
    /// EXPERIMENT (LOSAT_X_TBXUNGAPSENT): the batch's extensions may take the walks the
    /// sentinels stop (`x_tbx_sent_contexts`); a unit takes them where its chunk is the whole frame.
    sent: bool,
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
    // Dispatch point (per batch, after the lookup build): checks the gates and the guard (`decide`)
    // and builds the sub-lookups, or takes the ones PAIRLUT built at the lookup build site
    // (`lut_build`). None means the reference loop runs for every subject of the batch.
    /// The unit path for this batch, or None (a gate, the guard or the cost rule sends it to the
    /// reference loop).
    pub(super) fn prepare(b: XBatch<'_>) -> Option<Self> {
        let decision = decide(&XGates {
            contexts: b.contexts,
            query_length: b.query_length,
            diag_array_size: b.diag_array_size,
            subjects: b.subjects,
            debugging: b.debugging,
            chunk_switches: b.chunk_switches,
            threads: b.threads,
            cost_gate: cost_gate(),
            n1_gate: n1_gate(),
            threshold: b.threshold,
            band: tbx_band(),
            unit_mask: unit_mask_kib().is_some(),
        });
        BATCHES.fetch_add(1, AtomicOrdering::Relaxed);
        match decision {
            XDecision::Gate => {
                FALLBACK_GATE.fetch_add(1, AtomicOrdering::Relaxed);
            }
            XDecision::Guard => {
                FALLBACK_GUARD.fetch_add(1, AtomicOrdering::Relaxed);
            }
            XDecision::Cost => {
                FALLBACK_COST.fetch_add(1, AtomicOrdering::Relaxed);
            }
            XDecision::Units => {}
        }
        // The unit path runs for the batch (`Units`), or the switches that act only on it are left
        // with the reference loop.
        let units = decision == XDecision::Units;
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPAR", units);
        if seed_compact() {
            crate::utils::xreach::x_reach!("LOSAT_X_SEEDCOMPACT", units);
        }
        if x_kmerge::on() {
            crate::utils::xreach::x_reach!("LOSAT_X_TBXKMERGE", units);
        }
        if !units {
            x_tbx_hitbf::x_reach_not_taken();
        }
        if decision != XDecision::Units {
            // PAIRLUT: direct groups mean that the build site decided for the unit path (and that
            // the batch's table is a stub); the reference loop must not get this batch.
            assert!(
                b.lut.is_none(),
                "PAIRLUT: direct groups were built for a batch the unit path does not take ({decision:?})"
            );
            return None;
        }
        let (groups, sliced, stub) = match b.lut {
            // PAIRLUT: the groups built at the lookup build site (`lut_build`); the batch's table
            // is the stub the build site took.
            Some(lut) => {
                let stub = lut.stub.is_none();
                (lut.groups, lut.sliced, stub)
            }
            None => {
                let t0 = timing().then(Instant::now);
                // EXPERIMENT (LOSAT_X_TBXQSLICE): the slices filtered out of the batch's table, as
                // the groups are; the whole groups only where no group is cut.
                let starts = group_starts(b.contexts);
                let sliced = x_tbx_qslice::partition(b.contexts, &starts, b.threads)
                    .map(|parts| x_tbx_qslice::build_sliced(b.lookup, &parts));
                let groups = if sliced.is_none() {
                    build_groups(b.lookup, b.contexts)
                } else {
                    Vec::new()
                };
                add_ns(&NS_SUBLOOKUP, t0);
                (groups, sliced, false)
            }
        };
        // LOSAT_X_TBXBAND: a batch the first guard rejects, counted for the `[X_STATS] PAIRPAR`
        // line (without the switch no such batch gets here).
        let band = tbx_band() && {
            let (max_context_len, longest_subject_chunk) = x_guard_lengths(b.contexts, b.subjects);
            !guard(
                max_context_len,
                longest_subject_chunk,
                b.query_length,
                b.diag_array_size,
            )
        };
        if band {
            BAND.fetch_add(1, AtomicOrdering::Relaxed);
            crate::utils::xreach::x_reach!("LOSAT_X_TBXBAND", true);
        }
        // EXPERIMENT (LOSAT_X_TBXUNITMASK): the batch's rule and the groups' extents for the units'
        // masked tables; none where LOSAT_X_TBXQSLICE's rule cuts a group (module doc;
        // `x_band_masked`), so a batch takes the same tables with and without the diagnostics that
        // stop the slicing.
        let mask_rule = unit_mask_kib()
            .filter(|_| {
                let ungrouped = sliced.is_none()
                    && !x_tbx_qslice::rule_cuts(b.contexts, &group_starts(b.contexts), b.threads);
                // A batch with cut groups keeps the unmasked tables.
                if !ungrouped {
                    crate::utils::xreach::x_reach!("LOSAT_X_TBXUNITMASK", false);
                }
                ungrouped
            })
            .map(|kib| XMaskRule {
                max_cells: x_mask_max_cells(kib, b.threads),
                ratio: X_UNITMASK_RATIO,
                band_cells: band.then_some(b.diag_array_size),
            });
        let extents = if mask_rule.is_some() {
            x_group_extents(b.contexts)
        } else {
            Vec::new()
        };
        debug_assert!(extents.is_empty() || extents.len() == groups.len());
        // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): checked once for the batch's contexts.
        let sent = x_tbx_sent_contexts(b.contexts.iter().map(|c| c.aa_seq.as_slice()));
        Some(Self {
            groups,
            stub,
            sliced,
            mask_rule,
            extents,
            sent,
        })
    }

    /// PAIRLUT: the batch's lookup table is a stub (the units scan with the direct group
    /// tables; the reference frame loop must run zero times).
    pub(super) fn stub(&self) -> bool {
        self.stub
    }

    // No NCBI counterpart: LOSAT_X_TBXUNITMASK's choice of a unit's table (module doc).
    /// The mask of `unit`'s table (`x_unit_mask_for`), or None for the unmasked table.
    fn x_unit_mask(&self, unit: &XUnit, window: i32) -> Option<u32> {
        let rule = self.mask_rule.as_ref()?;
        x_unit_mask_for(
            rule,
            self.extents[unit.g],
            window,
            x_unit_width(&self.groups[unit.g], &unit.chunk),
            unit.chunk.length,
        )
    }

    // No NCBI counterpart: the size of a unit's table, for LOSAT_X_TBXPREFETCH's rule.
    /// The number of cells of `unit`'s table: the masked one's, or `x_unit_width`.
    fn x_unit_cells(&self, unit: &XUnit, window: i32) -> usize {
        match self.x_unit_mask(unit, window) {
            Some(mask) => mask as usize + 1,
            None => x_unit_width(&self.groups[unit.g], &unit.chunk),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:479-491
    // ```c
    // status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
    //                                dbseq_chunk_overlap);
    // ...
    // BlastInitHitListReset(init_hitlist);
    // ```
    // The seed stage of one subject, first half: the chunks of every frame and the units of every
    // (frame, chunk), run. This is the only part of the seed stage that reads the group tables
    // (`self.groups`: `longest_chain` while the units are listed, the tables in `scan_unit`).
    // `merge_subject` is the second half; the two were one function (`run_subject`) and are
    // called one after the other at dispatch point 2 of run_impl.rs, where the seed release
    // ends its hold on the group tables between them (x_seed_release.rs).
    /// The units of one subject, run: what `merge_subject` builds `combined_ungapped_hits` from.
    pub(super) fn run_subject_units(
        &self,
        inp: &XSubject<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> XSubjectUnits {
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
                for g in 0..self.x_group_count() {
                    if self.x_group_has_entries(g) {
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
        let outputs = if let Some(sliced) = &self.sliced {
            // EXPERIMENT (LOSAT_X_TBXQSLICE): the units of a batch with a cut group scan in
            // slices (x_tbx_qslice.rs); the units and their order are the same.
            self.run_subject_unit_fns_sliced(sliced, &units, inp, pool)
        } else {
            // EXPERIMENT (LOSAT_X_TBXPREFETCH, module doc): all units of the subject scan with the
            // prefetch when one of their tables is large enough; the units and their order are
            // the same either way.
            // LOSAT_X_TBXUNITMASK: the tables' own sizes (`x_unit_cells`).
            let prefetch = x_prefetch_units(
                tbx_prefetch_min_kib(),
                units.iter().map(|unit| {
                    self.x_unit_cells(unit, inp.window) * std::mem::size_of::<DiagStruct>()
                }),
            );
            if self.mask_rule.is_some() {
                let masked = units
                    .iter()
                    .filter(|unit| self.x_unit_mask(unit, inp.window).is_some())
                    .count();
                UM_UNITS.fetch_add(masked as u64, AtomicOrdering::Relaxed);
                if masked > 0 {
                    crate::utils::xreach::x_reach!("LOSAT_X_TBXUNITMASK", true);
                }
                if masked < units.len() {
                    crate::utils::xreach::x_reach!("LOSAT_X_TBXUNITMASK", false);
                }
            }
            if prefetch {
                self.run_subject_unit_fns::<true>(&units, inp, pool)
            } else {
                self.run_subject_unit_fns::<false>(&units, inp, pool)
            }
        };
        add_ns(&NS_UNITS, t_units);
        XSubjectUnits {
            frame_chunks,
            units,
            outputs,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:561-562,578-581
    // ```c
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
    // The seed stage of one subject, second half: the list the reference builds per chunk in order
    // (sort, conversion, offset adjustment and chunk merge, RI:3582-3598, 4533-4553) and per frame
    // in s_f_idx order (append and sort, RI:4556-4581); `merge_subject_listed` runs those
    // statements, LOSAT_X_TBXKMERGE (x_kmerge.rs) computes the same list without the full sorts.
    // It has no `self`: it reads the units' own lists and the subject's inputs, never a group
    // table, so the tables may already be freed when it runs.
    /// `combined_ungapped_hits` of one subject, from its units.
    pub(super) fn merge_subject(x_units: XSubjectUnits, inp: &XSubject<'_>) -> Vec<UngappedHit> {
        // EXPERIMENT (LOSAT_X_TBXKMERGE), dispatch: the same list without the full sorts
        // (x_kmerge.rs); with the switch off the body below runs.
        if x_kmerge::on() {
            return x_kmerge::merge_subject(x_units, inp);
        }
        Self::merge_subject_listed(x_units, inp)
    }

    /// `merge_subject` with the switch LOSAT_X_TBXKMERGE off: the units' lists of a chunk
    /// concatenated in unit order and sorted, converted, merged into the frame list, and the frame
    /// lists appended and sorted (the statements of the reference).
    fn merge_subject_listed(x_units: XSubjectUnits, inp: &XSubject<'_>) -> Vec<UngappedHit> {
        let XSubjectUnits {
            frame_chunks,
            units,
            outputs,
        } = x_units;
        // EXPERIMENT (LOSAT_X_SEEDCOMPACT): the form `run_subject_units` chose, read once per
        // subject; today's lists continue into today's statements below.
        let outputs = match outputs {
            XUnitOutputs::Full(outputs) => outputs,
            XUnitOutputs::Compact(records) => {
                return Self::merge_subject_compact(frame_chunks, units, records, inp);
            }
        };
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:561-562,578-581
    // ```c
    // BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
    //         hit_params->options, &hsp_list);
    // ...
    // Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
    // overlap = (backup.offset == backup.hard_ranges[backup.hm_index].left) ?
    //           0 : dbseq_chunk_overlap;
    // status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2876-2881
    // ```c
    //    /* If no previous HSP list, just return a copy of the new one. */
    //    if (!combined_hsp_list) {
    //       *combined_hsp_list_ptr = hsp_list;
    //       *hsp_list_ptr = NULL;
    //       return 0;
    //    }
    // ```
    // EXPERIMENT (LOSAT_X_SEEDCOMPACT): `merge_subject` on the compact records (section above
    // `XUnitHsp`, items 3 and a-d). The statements are those of `merge_subject_listed` except three: a
    // unit's list is appended to the chunk's init list by `x_expand_unit` instead of `extend`, the
    // conversion fills a Vec reserved at the init list's length (`get_ungapped_hsp_list_into`), and
    // a chunk list that meets an empty frame list becomes the frame list, as NCBI hands the first
    // list on (above) and as the empty-list branch of `merge_tblastx_subject_chunk_hits`
    // (`combined.append(&mut incoming)`) computes, without its copy. Each unit's records are
    // dropped once they are expanded.
    /// `combined_ungapped_hits` of one subject, from its units' compact records.
    fn merge_subject_compact(
        frame_chunks: Vec<Vec<SubjectChunk>>,
        units: Vec<XUnit>,
        records: Vec<Vec<XUnitHsp>>,
        inp: &XSubject<'_>,
    ) -> Vec<UngappedHit> {
        let t_merge = timing().then(Instant::now);
        SUBJECTS.fetch_add(1, AtomicOrdering::Relaxed);
        UNITS.fetch_add(units.len() as u64, AtomicOrdering::Relaxed);
        SC_SUBJECTS.fetch_add(1, AtomicOrdering::Relaxed);
        SC_UNITS.fetch_add(units.len() as u64, AtomicOrdering::Relaxed);
        SC_RECORDS.fetch_add(
            records.iter().map(|r| r.len() as u64).sum::<u64>(),
            AtomicOrdering::Relaxed,
        );

        let mut records = records.into_iter();
        let mut next_unit = 0usize;
        let mut combined_ungapped_hits: Vec<UngappedHit> = Vec::new();
        for (f, chunks) in frame_chunks.iter().enumerate() {
            let mut frame_ungapped_hits: Vec<UngappedHit> = Vec::new();
            for (k, chunk) in chunks.iter().enumerate() {
                let mut init_hsps: Vec<InitHSP> = Vec::new();
                while next_unit < units.len() && units[next_unit].f == f && units[next_unit].k == k
                {
                    let unit_records = records.next().expect("one output per unit");
                    x_expand_unit(&unit_records, units[next_unit].f, inp, &mut init_hsps);
                    next_unit += 1;
                }
                INIT_HSPS.fetch_add(init_hsps.len() as u64, AtomicOrdering::Relaxed);
                sort_init_hsps_by_score_ncbi(&mut init_hsps);
                let mut hits = if init_hsps.is_empty() {
                    Vec::new()
                } else {
                    let exact = Vec::with_capacity(init_hsps.len());
                    get_ungapped_hsp_list_into(init_hsps, inp.contexts, inp.s_frames, exact)
                };
                adjust_tblastx_chunk_subject_offsets(&mut hits, chunk.offset);
                if !hits.is_empty() && frame_ungapped_hits.is_empty() {
                    frame_ungapped_hits = hits;
                } else if !hits.is_empty() {
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

    // EXPERIMENT (LOSAT_X_SEEDCOMPACT): one read of the switch per subject chooses the form of the
    // units' lists (section above `XUnitHsp`); the units, their order and the scan are the same in
    // both.
    // Diagnostic (LOSAT_TIMING, LOSAT_X_STATS): with either variable the same units run with a
    // probe (x_unit_probe.rs); without them, the functions without probes.
    // EXPERIMENT (LOSAT_X_TBXPREFETCH): `PF` picks the unit functions with the prefetch (module
    // doc); each (form, probe, `PF`) is its own function item, so `run_units` and the scan are
    // compiled once for each and the scan without the hint keeps its own code.
    // EXPERIMENT (LOSAT_X_TBXKMERGE): the units sort their lists.
    // No NCBI counterpart: dispatch only.
    /// The units of one subject, run by the unit function of the switches' choice.
    fn run_subject_unit_fns<const PF: bool>(
        &self,
        units: &[XUnit],
        inp: &XSubject<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> XUnitOutputs {
        let sort = x_kmerge::on();
        // EXPERIMENT (LOSAT_X_TBXEXTBATCH): the batched extensions (x_tbx_extbatch.rs), without
        // the probes, for the subjects of its rule (`x_tbx_extbatch::active`: the window leaves
        // room for its pending marks and a unit's query-group span reaches the switch's value).
        if x_unit_probe::mode() == 0
            && x_tbx_extbatch::active(
                inp.window,
                inp.wordsize,
                units.iter().map(|unit| x_unit_span(&self.groups[unit.g])),
            )
        {
            return if seed_compact() {
                XUnitOutputs::Compact(self.run_units(
                    units,
                    inp,
                    pool,
                    Self::run_unit_compact_batched::<PF>,
                    sort,
                ))
            } else {
                XUnitOutputs::Full(self.run_units(
                    units,
                    inp,
                    pool,
                    Self::run_unit_batched::<PF>,
                    sort,
                ))
            };
        }
        match (seed_compact(), x_unit_probe::mode()) {
            (true, 0) => XUnitOutputs::Compact(self.run_units(
                units,
                inp,
                pool,
                Self::run_unit_compact::<PF>,
                sort,
            )),
            (false, 0) => {
                XUnitOutputs::Full(self.run_units(units, inp, pool, Self::run_unit::<PF>, sort))
            }
            (true, 1) => XUnitOutputs::Compact(self.run_units(
                units,
                inp,
                pool,
                Self::run_unit_compact_probed::<XTimeProbe, PF>,
                sort,
            )),
            (false, 1) => XUnitOutputs::Full(self.run_units(
                units,
                inp,
                pool,
                Self::run_unit_probed::<XTimeProbe, PF>,
                sort,
            )),
            (true, _) => XUnitOutputs::Compact(self.run_units(
                units,
                inp,
                pool,
                Self::run_unit_compact_probed::<XStatsProbe, PF>,
                sort,
            )),
            (false, _) => XUnitOutputs::Full(self.run_units(
                units,
                inp,
                pool,
                Self::run_unit_probed::<XStatsProbe, PF>,
                sort,
            )),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149
    // ```c
    // TBlastThreads the_threads(GetNumberOfThreads());
    // ```
    // Scheduling only: the units run on the search pool (an indexed map, results in unit order),
    // or in order on the calling thread when the search has one thread. `unit_fn` is
    // `Self::run_unit` (today's lists) or `Self::run_unit_compact` (LOSAT_X_SEEDCOMPACT), each
    // with or without LOSAT_X_TBXPREFETCH's hint (`PF`), a function item, so this function is
    // compiled once for each. The scratches of the pool carry LOSAT_X_SEEDUNMAP as it acts in this
    // process (`seed_unmap_active`, `XScratch::unmap`); the serial scratch never maps (section
    // above `X_MAPPED_BYTES`, rule 3: at one thread nothing is retained, mapping only costs).
    // EXPERIMENT (LOSAT_X_TBXKMERGE, x_kmerge.rs item 1): with `sort` each unit's list is sorted
    // in the unit's own task, after its scan, by the comparator of the chunk's init sort.
    fn run_units<T: Send + x_kmerge::XUnitSort>(
        &self,
        units: &[XUnit],
        inp: &XSubject<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
        unit_fn: impl Fn(&Self, &XUnit, &XSubject<'_>, &mut XScratch) -> T + Sync,
        sort: bool,
    ) -> Vec<T> {
        let run = |unit: &XUnit, scratch: &mut XScratch| {
            let t0 = std::time::Instant::now();
            let mut list = unit_fn(self, unit, inp, scratch);
            if sort {
                list.x_sort_in(scratch);
            }
            crate::utils::xtline::rec("unit", inp.s_idx as u32, t0, std::time::Instant::now());
            list
        };
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() {
            let unmap = seed_unmap_active();
            return pool.install(|| {
                units
                    .par_iter()
                    .map_init(|| XScratch::new(unmap), |scratch, unit| run(unit, scratch))
                    .collect()
            });
        }
        let _ = pool;
        let mut scratch = XScratch::new(false);
        units.iter().map(|unit| run(unit, &mut scratch)).collect()
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491
    // ```c
    //         BlastInitHitListReset(init_hitlist);
    // ```
    // The saved init HSPs of one unit as today's list, which starts empty as NCBI's list of a
    // chunk does (`scan_unit` fills it).
    /// The saved init HSPs of one unit, in its scan order.
    fn run_unit<const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<InitHSP> {
        let mut init_hsps: Vec<InitHSP> = Vec::new();
        // EXPERIMENT (LOSAT_X_TBXHITBF): `scan_unit` or its branch-free hit loop (x_tbx_hitbf.rs).
        self.scan_unit_or_bf::<_, PF>(unit, inp, scratch, &mut init_hsps);
        init_hsps
    }

    // Diagnostic (LOSAT_TIMING, LOSAT_X_STATS): `run_unit` with a probe (x_unit_probe.rs). Its
    // own function, so that `run_unit` keeps the scan without probes. No NCBI counterpart.
    /// `run_unit`, recording into the thread's probe `P`.
    #[inline(never)]
    fn run_unit_probed<P: XProbeSlot, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<InitHSP> {
        let mut probe = std::mem::take(P::slot(scratch));
        let mut init_hsps: Vec<InitHSP> = Vec::new();
        self.scan_unit::<_, _, PF>(unit, inp, scratch, &mut init_hsps, &mut probe);
        *P::slot(scratch) = probe;
        init_hsps
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1043
    // ```c
    //     aux_struct->init_hitlist = BLAST_InitHitListNew();
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:229-235
    // ```c
    // void BlastInitHitListReset(BlastInitHitList * init_hitlist)
    // {
    //     Int4 index;
    //
    //     for (index = 0; index < init_hitlist->total; ++index)
    //         sfree(init_hitlist->init_hsp_array[index].ungapped_data);
    //     init_hitlist->total = 0;
    // ```
    // EXPERIMENT (LOSAT_X_SEEDUNMAP, section above `X_MAPPED_BYTES`): whether the blocks of `unit`
    // are mapped (`x_unit_mapped` on its query-group span and chunk), and for such a unit the
    // cleared record vector replaced by an empty one of the mapped capacity unless it has that much
    // room. Shared by `run_unit_compact` and `run_unit_compact_probed`. No NCBI counterpart.
    /// `x_unit_mapped` for `unit`; where it holds, `records` (empty) gets the mapped capacity.
    fn x_unit_records(
        &self,
        unit: &XUnit,
        scratch: &XScratch,
        records: &mut Vec<XUnitHsp>,
    ) -> bool {
        let span = x_unit_span(&self.groups[unit.g]);
        let mapped = x_unit_mapped(scratch.unmap, span, unit.chunk.length);
        if scratch.unmap {
            crate::utils::xreach::x_reach!("LOSAT_X_SEEDUNMAP", mapped);
        }
        if mapped && records.capacity() * std::mem::size_of::<XUnitHsp>() < X_MAPPED_BYTES {
            *records = Vec::with_capacity(x_mapped_len::<XUnitHsp>(0));
            SU_RECORDS.fetch_add(1, AtomicOrdering::Relaxed);
        }
        mapped
    }

    // No NCBI counterpart: the unit listing of the unit path (scheduling only).
    /// The number of context groups of the batch (the `g` of its units).
    fn x_group_count(&self) -> usize {
        match &self.sliced {
            Some(sliced) => sliced.len(),
            None => self.groups.len(),
        }
    }

    // No NCBI counterpart: the unit listing of the unit path (scheduling only).
    /// Whether group `g` has lookup entries (its units are listed): its table, or LOSAT_X_TBXQSLICE's
    /// slices of it.
    fn x_group_has_entries(&self, g: usize) -> bool {
        match &self.sliced {
            Some(sliced) => sliced[g].has_entries(),
            None => self.groups[g].lookup.longest_chain > 0,
        }
    }

    // EXPERIMENT (LOSAT_X_SEEDCOMPACT): the records of one unit (section above `XUnitHsp`, item
    // 2). NCBI keeps one init hit list per thread and empties it for every chunk, keeping its
    // array; the unit fills its thread's scratch vector, emptied here and kept for the thread's
    // next unit, and returns an exact copy of it (`to_vec`: the records in their order, with
    // no spare room), since the outputs of a subject's units live until its merge.
    // EXPERIMENT (LOSAT_X_SEEDUNMAP, section above `X_MAPPED_BYTES`): for a unit where
    // `x_unit_mapped` holds, the scratch vector is replaced, before the scan, by an empty one of the
    // mapped capacity unless it has that much room, and a copy of at least `X_MAPPED_MIN` bytes is
    // the records in their order in a vector of the mapped capacity (its pages beyond the records
    // are never touched).
    /// The saved init HSPs of one unit as compact records, in its scan order, at their exact length
    /// (LOSAT_X_SEEDUNMAP: a large copy has the mapped capacity instead).
    fn run_unit_compact<const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<XUnitHsp> {
        let mut records = std::mem::take(&mut scratch.records);
        records.clear();
        let mapped = self.x_unit_records(unit, scratch, &mut records);
        // EXPERIMENT (LOSAT_X_TBXHITBF): `scan_unit` or its branch-free hit loop (x_tbx_hitbf.rs).
        self.scan_unit_or_bf::<_, PF>(unit, inp, scratch, &mut records);
        let exact = x_unit_copy(&records, mapped);
        scratch.records = records;
        exact
    }

    // Diagnostic (LOSAT_TIMING, LOSAT_X_STATS): `run_unit_compact` with a probe
    // (x_unit_probe.rs), in its own function as `run_unit_probed`; LOSAT_X_SEEDUNMAP acts here as
    // in `run_unit_compact` (`x_unit_records`, `x_unit_copy`). No NCBI counterpart.
    /// `run_unit_compact`, recording into the thread's probe `P`.
    #[inline(never)]
    fn run_unit_compact_probed<P: XProbeSlot, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<XUnitHsp> {
        let mut probe = std::mem::take(P::slot(scratch));
        let mut records = std::mem::take(&mut scratch.records);
        records.clear();
        let mapped = self.x_unit_records(unit, scratch, &mut records);
        self.scan_unit::<_, _, PF>(unit, inp, scratch, &mut records, &mut probe);
        let exact = x_unit_copy(&records, mapped);
        scratch.records = records;
        *P::slot(scratch) = probe;
        exact
    }

    // The scan loop of the reference (RI:3082-3168) on one unit: the group's sub-lookup, the frame
    // chunk, a fresh table (RI:2930-2933: `DiagStruct::default()`, diag_offset = window) indexed by
    // the raw diagonal minus the group's lowest one, without the mask (module doc, item 3). The
    // body of what was `run_unit`; the list it saves into is the caller's (`XInitSink`).
    // EXPERIMENT (LOSAT_X_TBXPREFETCH, module doc): with `PF`, on x86-64, each hit of a batch
    // first prefetches the cell of the hit X_PREFETCH_AHEAD entries after it; without `PF` the
    // loop is the one before the switch.
    /// Scans one unit and saves its init HSPs, in scan order, into `init_hsps`; `probe` records
    /// the diagnostics of LOSAT_TIMING / LOSAT_X_STATS (`XNoProbe`: nothing); `PF` adds the
    /// prefetch of LOSAT_X_TBXPREFETCH.
    #[inline(always)]
    fn scan_unit<S: XInitSink, P: XProbe, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
        probe: &mut P,
    ) {
        self.scan_unit_hits::<S, P, PF, false, false>(unit, inp, scratch, init_hsps, probe);
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF, x_tbx_hitbf.rs): with `BF` the hits of a batch go through
    // `x_hits_bf`, which takes the too-far and overlap outcomes of `x_hit` without a branch and
    // hands the others to `x_hit`; without it this is the loop before the switch. `BF` runs only
    // with `XNoProbe` (the probes count every outcome) and where `x_tbx_hitbf::fits` holds.
    /// `scan_unit`, with `BF` the branch-free hit loop of LOSAT_X_TBXHITBF and with `SENT` the
    /// extensions of `extend_hit_two_hit_sent` (LOSAT_X_TBXUNGAPSENT).
    fn scan_unit_hits<S: XInitSink, P: XProbe, const PF: bool, const BF: bool, const SENT: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
        probe: &mut P,
    ) {
        // EXPERIMENT (LOSAT_X_TBXUNITMASK): the masked table where the rule holds (module doc).
        match self.x_unit_mask(unit, inp.window) {
            None => self.scan_unit_cells::<S, P, XUnmasked, PF, BF, SENT>(
                unit, inp, scratch, init_hsps, probe, 0,
            ),
            Some(mask) => self
                .scan_unit_masked::<S, P, PF, BF, SENT>(unit, inp, scratch, init_hsps, probe, mask),
        }
    }

    // EXPERIMENT (LOSAT_X_TBXUNITMASK): the masked form in a function of its own, so that the
    // unmasked form inlined into `scan_unit_hits` keeps the code it had before the switch.
    /// `scan_unit_cells` with `XMasked`.
    #[inline(never)]
    fn scan_unit_masked<
        S: XInitSink,
        P: XProbe,
        const PF: bool,
        const BF: bool,
        const SENT: bool,
    >(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
        probe: &mut P,
        mask: u32,
    ) {
        self.scan_unit_cells::<S, P, XMasked, PF, BF, SENT>(
            unit, inp, scratch, init_hsps, probe, mask,
        );
    }

    // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): never on wasm32 (line W keeps its walks there).
    /// Whether the extensions of `unit` take the walks the sentinels stop: the batch's contexts
    /// allow them (`sent`) and the unit's chunk is its whole subject frame.
    #[inline(always)]
    pub(super) fn x_unit_sent(&self, unit: &XUnit, inp: &XSubject<'_>) -> bool {
        cfg!(not(target_arch = "wasm32"))
            && x_tbx_sent_subject(
                self.sent,
                &inp.s_frames_preliminary[unit.f].aa_seq,
                unit.chunk.offset,
                unit.chunk.length,
            )
    }

    // The statements of `scan_unit_hits` for one form of the unit's cells (`XCells`): with
    // `XUnmasked` the table and the loop before LOSAT_X_TBXUNITMASK, with `XMasked` a table of
    // `mask + 1` cells and the masked cell of every hit.
    /// `scan_unit_hits` with the cells `C` (`mask`: the masked table's cells minus one).
    #[inline(always)]
    fn scan_unit_cells<
        S: XInitSink,
        P: XProbe,
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
        probe: &mut P,
        mask: u32,
    ) {
        // LOSAT_X_TBXUNGAPSENT: `SENT` only on the units `x_unit_sent` accepts (`x_hit`'s SAFETY;
        // once per unit).
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
            return;
        }

        let offset_array_size: i32 = OFFSET_ARRAY_SIZE + sub.longest_chain.max(0);
        // EXPERIMENT (LOSAT_X_TBXSCANBF): with the switch, X_SCAN_BF_SLACK entries more.
        // EXPERIMENT (LOSAT_X_TBXHITBF): with `BF` and `PF`, X_PREFETCH_AHEAD entries more
        // (`x_hits_bf` reads the entry that far ahead without a bound test).
        let offset_pairs_len =
            crate::algorithm::tblastx::x_scan_bf::x_offset_pairs_len(offset_array_size)
                + if BF && PF { X_PREFETCH_AHEAD } else { 0 };
        if scratch.offset_pairs.len() < offset_pairs_len {
            scratch
                .offset_pairs
                .resize(offset_pairs_len, OffsetPair::default());
        }
        // Unit cell of hit (q, s): q - s - lo, lo = q_min - (s_len - 1) (u32 wrapping), so the
        // cells run from 0 to (q_max - q_min) + (s_len - 1). LOSAT_X_TBXUNITMASK: that value
        // masked to the table's mask + 1 cells.
        let lo = group.q_min.wrapping_sub(chunk.length as u32 - 1);
        let cells = C::of(lo, mask);
        let width = if C::MASKED {
            mask as usize + 1
        } else {
            x_unit_width(group, &chunk)
        };
        probe.unit_begin(width);
        // EXPERIMENT (LOSAT_X_SEEDUNMAP, section above `X_MAPPED_BYTES`): a table too small for a
        // unit where `x_unit_mapped` holds is replaced by an empty one of the mapped capacity; the
        // clear and resize below write every cell of the unit's width as before.
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
                // One batch per call (`x_unit_batch`: its own function on wasm32 only).
                if !x_unit_batch::<S, P, C, PF, BF, SENT>(
                    group,
                    subject,
                    &seq_ranges,
                    offset_array_size,
                    &mut scan_range,
                    scratch,
                    cells,
                    width,
                    &h,
                    init_hsps,
                    probe,
                ) {
                    // The reference's safety guard (RI:3157-3161).
                    break;
                }
            }
        }
        probe.unit_end();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:501-508
// ```c
// while (scan_range[1] <= scan_range[2]) {
//     /* scan the subject sequence for hits */
//     hits = scansub(lookup_wrap, subject,
//                               offset_pairs, array_size, scan_range);
// ...
//     for (i = 0; i < hits; ++i) {
// ```
// One pass of the scan loop of `scan_unit` (one `scansub` call and the two-hit test of its hits),
// statement for statement the loop body it was. No NCBI counterpart for the function boundary: on
// wasm32 (LOSAT Web) it is a function of its own, so that the later batches of a unit run the
// optimized code. V8 runs the first calls of a freshly compiled function in its baseline code
// (Liftoff) and switches to its optimized code (TurboFan) only at a later call, never inside a
// running call (no on-stack replacement for Wasm), and a unit was one long call. Elsewhere it is
// inlined into `scan_unit`.
/// One batch of `scan_unit`'s loop: scans the next hits of `subject` and runs the two-hit test of
/// each on its cell of the unit table; false when the loop stops (no hit and no progress).
#[cfg_attr(target_arch = "wasm32", inline(never))]
#[cfg_attr(not(target_arch = "wasm32"), inline(always))]
#[allow(clippy::too_many_arguments)]
fn x_unit_batch<
    S: XInitSink,
    P: XProbe,
    C: XCells,
    const PF: bool,
    const BF: bool,
    const SENT: bool,
>(
    group: &XGroup,
    subject: &[u8],
    seq_ranges: &[(i32, i32)],
    offset_array_size: i32,
    scan_range: &mut [i32; 3],
    scratch: &mut XScratch,
    cells: C,
    width: usize,
    h: &XHitCtx<'_>,
    init_hsps: &mut S,
    probe: &mut P,
) -> bool {
    let prev_scan_left = scan_range[1];
    // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
    // EXPERIMENT (LOSAT_X_TBXSCANBF): the branch-free scan with the switch.
    let hits = crate::algorithm::tblastx::x_scan_bf::x_scan_subject(
        group.scan_bf,
        &group.lookup,
        subject,
        seq_ranges,
        &mut scratch.offset_pairs[..],
        offset_array_size,
        scan_range,
    );
    if hits == 0 && scan_range[1] == prev_scan_left {
        return false;
    }
    probe.batch(hits as usize);
    let diag_ptr = scratch.diag.as_mut_ptr();
    let offset_pairs_ptr = scratch.offset_pairs.as_ptr();
    // EXPERIMENT (LOSAT_X_TBXHITBF): the loop below, with the too-far and overlap outcomes of
    // `x_hit` taken without a branch (x_tbx_hitbf.rs).
    if BF {
        // SAFETY: hits <= offset_array_size, and with `PF` the X_PREFETCH_AHEAD entries after them
        // are in offset_pairs too (`scan_unit_hits` sizes it); every cell is below width =
        // scratch.diag.len() as below; `fits` holds (the caller's rule). `x_hit` gets the cell of
        // its hit, as below.
        unsafe {
            x_tbx_hitbf::x_hits_bf::<C, PF>(
                diag_ptr,
                offset_pairs_ptr,
                hits as usize,
                cells,
                width,
                h.window,
                h.wordsize,
                |cell, query_offset, subject_offset| {
                    x_hit::<S, P, SENT>(
                        &mut *diag_ptr.add(cell),
                        query_offset,
                        subject_offset,
                        h,
                        init_hsps,
                        probe,
                    );
                },
            )
        };
        probe.batch_end();
        return true;
    }
    // [C] for (i = 0; i < hits; ++i)
    for i in 0..hits as usize {
        // EXPERIMENT (LOSAT_X_TBXPREFETCH): the cell of the hit X_PREFETCH_AHEAD
        // entries ahead, computed as below. A prefetch never faults and has no
        // architectural effect.
        #[cfg(target_arch = "x86_64")]
        if PF && i + X_PREFETCH_AHEAD < hits as usize {
            // SAFETY: i + X_PREFETCH_AHEAD < hits <= offset_pairs.len(), an entry
            // `scansub` wrote in this batch. The address is only prefetched
            // (`wrapping_add`: no bound needed; it lies in the table as below).
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
        // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
        // Unit cell: the raw diagonal minus the unit's lowest one (`XUnmasked`: no mask;
        // LOSAT_X_TBXUNITMASK's `XMasked`: masked to the table).
        let diag_coord = cells.cell(query_offset, subject_offset);
        debug_assert!(diag_coord < width, "unit diagonal outside its table");
        // SAFETY: q_min <= query_offset <= q_max (the group's entries) and
        // subject_offset < s_len, so an unmasked diag_coord is below width = scratch.diag.len();
        // a masked one is at most the mask, width - 1.
        let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };
        probe.at(diag_coord, subject_offset);
        let outcome = x_hit::<S, P, SENT>(
            diag_entry,
            query_offset,
            subject_offset,
            h,
            init_hsps,
            probe,
        );
        probe.hit(outcome, h.window);
    }
    probe.batch_end();
    true
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
// saved (both read only local values, so the order is immaterial). The init HSP goes to the
// caller's sink (today's `Vec<InitHSP>`, or the records of LOSAT_X_SEEDCOMPACT, section above
// `XUnitHsp`), by the same push statement.
/// The two-hit test and extension of one hit on its cell (`SENT`: LOSAT_X_TBXUNGAPSENT, the
/// extension of `extend_hit_two_hit_sent`).
#[inline(always)]
fn x_hit<S: XInitSink, P: XProbe, const SENT: bool>(
    diag_entry: &mut DiagStruct,
    query_offset: u32,
    subject_offset: u32,
    h: &XHitCtx<'_>,
    init_hsps: &mut S,
    probe: &mut P,
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
        // Diagnostic (LOSAT_TIMING, LOSAT_X_STATS): the probe around the extension
        // (x_unit_probe.rs; `XNoProbe` does nothing). No NCBI counterpart.
        let mark = probe.ext_begin();
        // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): with `SENT` the walks the sentinels stop.
        let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) = if !SENT {
            extend_hit_two_hit(
                query,
                h.subject,
                (last_hit + wordsize) as usize,
                subject_offset as usize,
                q_raw,
                x_dropoff,
                false,
            )
        } else {
            // SAFETY: `SENT` runs on the units `x_unit_sent` accepts (`scan_unit_cells`): the
            // batch's contexts passed `x_tbx_sent_contexts` (XPairPar::sent), and `h.subject_sent`
            // is the buffer of a frame `x_tbx_sent_subject` accepted.
            unsafe {
                extend_hit_two_hit_sent(
                    query_full,
                    h.subject_sent,
                    (last_hit + wordsize) as usize,
                    subject_offset as usize,
                    q_raw,
                    x_dropoff,
                )
            }
        };
        probe.ext_end(
            mark,
            query,
            h.subject,
            (last_hit + wordsize) as usize,
            subject_offset as usize,
            q_raw,
            x_dropoff,
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
    let starts = group_starts(contexts);
    let bounds: Vec<i32> = starts.iter().map(|&c| lookup.frame_bases[c]).collect();
    build_tables(lookup, &bounds)
}

// No NCBI counterpart beyond `build_groups`'s: its filter for any partition of the query offsets
// into intervals (LOSAT_X_TBXQSLICE's slices, x_tbx_qslice.rs, or the groups).
/// One table per interval of query offsets, each chain of `lookup` filtered to the interval with
/// its order kept: interval i holds the offsets q with `bounds[i] <= q < bounds[i + 1]` (the first
/// also those below `bounds[0]`, the last those above). `bounds` ascends.
fn build_tables(lookup: &BlastAaLookupTable, bounds: &[i32]) -> Vec<XGroup> {
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
    // The group of a query offset: the last group whose first context starts at or before it
    // (the first group when none does).
    debug_assert!(bounds.windows(2).all(|w| w[0] <= w[1]), "bounds ascend");
    let group_of = |q: i32| -> usize { bounds.partition_point(|&b| b <= q).saturating_sub(1) };
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
        .map(|b| {
            let lookup = BlastAaLookupTable {
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
            };
            XGroup {
                q_min: if b.longest_chain > 0 { b.q_min } else { 0 },
                q_max: b.q_max,
                scan_bf: crate::algorithm::tblastx::x_scan_bf::x_scan_bf_for(&lookup),
                lookup,
            }
        })
        .collect()
}

/// Every field of `a` equals the one of `b` (the e-value by its bits); the comparison of the unit
/// tests.
#[cfg(test)]
fn x_same_hit(a: &UngappedHit, b: &UngappedHit) -> bool {
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

#[cfg(test)]
mod tests {
    use super::x_kmerge::XUnitSort;
    use super::*;
    use crate::algorithm::tblastx::lookup::x_pair_lut::{
        x_assert_same_contexts, x_assert_same_group, x_assert_stub_scalars,
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:343-346
    // ```c
    //             PV_SET(pv, i, PV_ARRAY_BTS);
    //             bbc[i].num_used = lookup->thin_backbone[i][1];
    //             /* if there are three or fewer hits, */
    //             if (lookup->thin_backbone[i][1] <= AA_HITS_PER_CELL)
    // ```
    // No computation: compares the direct groups with the groups `build_groups` filtered out of the
    // reference builder's table (every field of every table, `q_min`, `q_max`), the stub with that
    // table, and the contexts; the first difference stops the test.
    fn lut_shadow_compare(
        lut: &XLutBatch,
        groups: &[XGroup],
        sliced: Option<&[XSliced]>,
        lookup: &BlastAaLookupTable,
        contexts: &[QueryContext],
    ) {
        assert_eq!(
            lut.groups.len(),
            groups.len(),
            "number of groups differs (direct vs reference)"
        );
        for (g, (direct, reference)) in lut.groups.iter().zip(groups).enumerate() {
            x_assert_same_group(
                g,
                &direct.lookup,
                (direct.q_min, direct.q_max),
                &reference.lookup,
                (reference.q_min, reference.q_max),
            );
        }
        // EXPERIMENT (LOSAT_X_TBXQSLICE): the slices, table by table, and each group's bounds.
        x_tbx_qslice::lut_shadow_compare(lut.sliced.as_deref(), sliced);
        let (stub, stub_contexts) = lut
            .stub
            .as_ref()
            .expect("the batch keeps the stub for the comparison");
        x_assert_stub_scalars(stub, lookup);
        x_assert_same_contexts(stub_contexts, contexts);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:844
    // ```c
    // if (Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax)) {
    // ```
    // No computation: compares the subject's combined list of the unit path with the reference
    // loop's, every field, and stops at the first difference.
    /// Check of one subject: the unit path's `combined_ungapped_hits` against the reference's.
    fn shadow_compare(s_idx: usize, unit: &[UngappedHit], reference: &[UngappedHit]) {
        let n = unit.len().max(reference.len());
        for k in 0..n {
            match (unit.get(k), reference.get(k)) {
                (Some(a), Some(b)) if x_same_hit(a, b) => {}
                (a, b) => panic!(
                    "subject {s_idx} frame {:?} index {k} differs (unit {} hits, reference {} hits): unit {a:?} vs reference {b:?}",
                    a.or(b).map(|h| h.s_f_idx),
                    unit.len(),
                    reference.len()
                ),
            }
        }
    }

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
    /// of every frame, how often each branch was taken, and how many hits met a cell whose flag 1
    /// was written by a right extension on another raw diagonal (the masked-aliasing case of
    /// module doc item 3).
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
    ) -> (Vec<Vec<InitHSP>>, Vec<XHit>, usize) {
        let wordsize = 3;
        let mask = (diag_array.len() - 1) as u32;
        let mut per_frame = Vec::new();
        let mut outcomes = Vec::new();
        // The raw diagonal of the right extension that last set each cell's flag.
        let mut flag_writer = vec![i64::MIN; diag_array.len()];
        let mut aliased = 0usize;
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
                    subject_sent: &[],
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
                        let raw = pair.q_off as i64 - pair.s_off as i64;
                        if diag_array[cell].flag() != 0
                            && flag_writer[cell] != i64::MIN
                            && flag_writer[cell] != raw
                        {
                            aliased += 1;
                        }
                        let outcome = x_hit::<_, _, false>(
                            &mut diag_array[cell],
                            pair.q_off,
                            pair.s_off,
                            &h,
                            &mut init_hsps,
                            &mut XNoProbe,
                        );
                        if matches!(outcome, XHit::Extended | XHit::Saved)
                            && diag_array[cell].flag() != 0
                        {
                            flag_writer[cell] = raw;
                        }
                        outcomes.push(outcome);
                    }
                }
            }
            advance_tblastx_diag_offset(diag_offset, diag_array, window, subject.len());
            sort_init_hsps_by_score_ncbi(&mut init_hsps);
            per_frame.push(init_hsps);
        }
        (per_frame, outcomes, aliased)
    }

    /// The unit path's sorted init-HSP list of every frame (one chunk per frame here).
    fn unit_lists(pp: &XPairPar, inp: &XSubject<'_>) -> Vec<Vec<InitHSP>> {
        let mut out = Vec::new();
        let mut scratch = XScratch::new(false);
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
                        init_hsps.extend(pp.run_unit::<false>(&unit, inp, &mut scratch));
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
            stub: false,
            sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
            sliced: None,
            mask_rule: None,
            extents: Vec::new(),
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
        let (reference, outcomes, _) = masked_reference(
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
    fn last_hit_fits_edges() {
        // Default window, a query of 20 M letters, the longest chunk NCBI makes: far inside.
        let ql = 20_000_000;
        let size = diag_array_size_for(ql, 40);
        assert!(last_hit_fits(5_000_000 / 3, ql, size));
        // A window just under 2^29 (the table has 2^30 cells): the reference's second scan
        // stores values of 2^30 and more, so no unit path.
        let window = 536_870_000;
        let size = diag_array_size_for(2_000, window);
        assert_eq!(size, 1 << 30);
        assert!(!last_hit_fits(20_000, 2_000, size));
        // The bound itself: slack + 2 * chunk + margin one below and at 2^29.
        let size = 1 << 29;
        let ql = (1 << 29) - (1 << 28);
        let half_margin = (X_LAST_HIT_MARGIN / 2) as usize;
        assert!(last_hit_fits((1 << 27) - half_margin - 1, ql, size));
        assert!(!last_hit_fits((1 << 27) - half_margin, ql, size));
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

    // No NCBI counterpart: LOSAT_X_TBXBAND (module doc item 3). Batches the first guard rejects: a
    // context and the subject longer than the slack, the subject made of stretches copied from the
    // long context so that extensions run past `slack + 3 + (wordsize - 1)` letters, the bound the
    // earlier argument could not exclude. The units' lists are the masked reference's in every
    // case, also with the batched extensions of LOSAT_X_TBXEXTBATCH.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:518-528,598-604
    // ```c
    // if (diag_array[diag_coord].flag) {
    //     if ((Int4) (subject_offset + diag_offset) <
    //         diag_array[diag_coord].last_hit) {
    //         continue;
    //     }
    // ...
    // if (right_extend) {
    //     diag_array[diag_coord].flag = 1;
    //     diag_array[diag_coord].last_hit =
    //         s_last_off - (wordsize - 1) + diag_offset;
    // ```
    #[test]
    fn units_equal_the_masked_reference_where_the_first_guard_fails() {
        let mut rng = Rng(0x7e57_ab1e_5eed_0001);
        let alphabets = [aa("AW"), aa("ACW"), aa("ACDEKW")];
        let window = 40;
        let (mut checked, mut long_reach, mut aliased_cases) = (0usize, 0usize, 0usize);
        for case_no in 0..400 {
            let alphabet = &alphabets[case_no % alphabets.len()];
            // One long context (longer than the slack) and up to two short ones.
            let long = 90 + rng.below(120) as usize;
            let mut frames = vec![frame(1, &random_residues(&mut rng, long, alphabet))];
            for f in 0..rng.below(3) as usize {
                let len = 3 + rng.below(20) as usize;
                frames.push(frame(
                    [2i8, 3][f],
                    &random_residues(&mut rng, len, alphabet),
                ));
            }
            let s_len = 150 + rng.below(350) as usize;
            let mut residues = Vec::with_capacity(s_len);
            while residues.len() < s_len {
                if rng.below(2) == 0 {
                    let from = rng.below(long as u64 / 2) as usize;
                    let take = (long - from).min(s_len - residues.len());
                    residues.extend_from_slice(&frames[0].aa_seq[1 + from..1 + from + take]);
                } else {
                    let take = (5 + rng.below(60) as usize).min(s_len - residues.len());
                    residues.extend(random_residues(&mut rng, take, alphabet));
                }
            }
            let s_frames = vec![frame(1, &residues)];
            let (lookup, contexts) = lookup_for(&frames);
            let query_length = query_length_of(&contexts);
            let size = diag_array_size_for(query_length, window);
            let max_ctx = contexts.iter().map(|c| c.aa_len).max().unwrap_or(0);
            if guard(max_ctx, s_len, query_length, size) {
                continue;
            }
            assert!(last_hit_fits(s_len, query_length, size));
            let slack = size - query_length;
            let cutoff_scores = vec![12; contexts.len()];
            let x_dropoff_per_context = vec![16; contexts.len()];
            let mut diag_array = vec![DiagStruct::default(); size as usize];
            let mut diag_offset = window;
            let (reference, outcomes, aliased) = masked_reference(
                &lookup,
                &contexts,
                &s_frames,
                &cutoff_scores,
                &x_dropoff_per_context,
                window,
                &mut diag_array,
                &mut diag_offset,
            );
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
            let pp = XPairPar {
                groups: build_groups(&lookup, &contexts),
                stub: false,
                sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
                sliced: None,
                mask_rule: None,
                extents: Vec::new(),
            };
            assert_eq!(
                unit_lists(&pp, &inp),
                reference,
                "case {case_no}: unit lists differ from the masked reference"
            );
            assert_eq!(
                unit_lists_batched(&pp, &inp),
                reference,
                "case {case_no}: batched unit lists differ from the masked reference"
            );
            // A saved extension that ran further right of its seed than the earlier argument's
            // bound allowed: the reach `s_last_off - s0` is at least `s_end - s_seed - 1`
            // (`s_end` is one past the best end, which is at most the last letter scanned), and the
            // earlier argument fails only for a reach above `slack + 3 + 2`.
            if reference
                .iter()
                .flatten()
                .any(|h| h.s_end - h.s_seed - 1 > slack + 3 + 2)
            {
                long_reach += 1;
            }
            // A hit of the reference met a flag set by a right extension on another raw diagonal
            // sharing its cell (the masked-aliasing case), and took the unmasked branch.
            if aliased > 0 {
                aliased_cases += 1;
            }
            assert!(outcomes.len() >= aliased);
            checked += 1;
        }
        assert!(checked > 150, "checked {checked}");
        assert!(long_reach > 50, "long_reach {long_reach}");
        assert!(aliased_cases > 20, "aliased_cases {aliased_cases}");
    }

    // No NCBI counterpart: LOSAT_X_TBXUNITMASK (module doc). Units on masked tables against the
    // same units on unmasked tables, in every hit loop (with and without LOSAT_X_TBXPREFETCH's hint
    // and LOSAT_X_TBXHITBF, and LOSAT_X_TBXEXTBATCH's batched extensions), and the sorted lists
    // against the reference loop on its own masked table. Query groups far shorter than the
    // subject, so that many raw diagonals share each masked cell; subjects made of pieces of the
    // query, so that long right extensions set flags that hits of aliasing diagonals then meet.
    #[test]
    fn masked_unit_tables_give_the_unmasked_units() {
        let mut rng = Rng(0x7e57_ab1e_5eed_0b1b);
        let alphabets = [aa("AW"), aa("ACW"), aa("ACDEKW")];
        let frame_numbers = [1i8, 2, 3, -1, -2, -3];
        let (mut masked_units, mut aliased_cases, mut saved) = (0usize, 0usize, 0usize);
        let (mut exact_fits, mut wide_groups, mut control_differs) = (0usize, 0usize, 0usize);
        let mut sent_units = 0usize;
        // LOSAT_X_TBXUNGAPSENT is on (the batches' contexts have their sentinels).
        let on = x_tbx_sent_contexts(std::iter::empty());
        for case_no in 0..300 {
            let alphabet = &alphabets[case_no % alphabets.len()];
            let window = [40, 3 + x_tbx_extbatch::BATCH as i32, 64][case_no % 3];
            // Every fourth case two queries (7 to 12 contexts: groups of several contexts); every
            // fifth one context lengths whose group extent plus the window is a power of two (a
            // group's extent is its length and the sentinel after it, the last group's its length).
            let queries = if case_no % 4 == 3 { 2 } else { 1 };
            let n_frames = if queries == 2 {
                7 + rng.below(6) as usize
            } else {
                1 + rng.below(6) as usize
            };
            let exact = queries == 1 && case_no % 5 == 0;
            let mut frames = Vec::with_capacity(n_frames);
            for i in 0..n_frames {
                let len = if exact {
                    let cells =
                        [64usize, 128][i % 2].max((window as usize + 5).next_power_of_two());
                    cells - window as usize - usize::from(i + 1 < n_frames)
                } else {
                    4 + rng.below(80) as usize
                };
                frames.push(frame(
                    frame_numbers[i % 6],
                    &random_residues(&mut rng, len, alphabet),
                ));
            }
            let s_len = 1500 + rng.below(2500) as usize;
            let mut residues = Vec::with_capacity(s_len);
            while residues.len() < s_len {
                if rng.below(3) != 0 {
                    let q = &frames[rng.below(n_frames as u64) as usize];
                    let len = q.aa_seq.len() - 2;
                    let from = rng.below(len as u64) as usize;
                    let take = (len - from).min(s_len - residues.len());
                    residues.extend_from_slice(&q.aa_seq[1 + from..1 + from + take]);
                } else {
                    let take = (3 + rng.below(40) as usize).min(s_len - residues.len());
                    residues.extend(random_residues(&mut rng, take, alphabet));
                }
            }
            let s_frames = vec![frame(1, &residues)];
            let (lookup, contexts) = if queries == 2 {
                let (a, b) = frames.split_at(n_frames / 2);
                let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
                build_ncbi_lookup(&[a.to_vec(), b.to_vec()], 11, &params, true)
            } else {
                lookup_for(&frames)
            };
            let query_length = query_length_of(&contexts);
            let size = diag_array_size_for(query_length, window);
            let cutoff_scores = vec![12; contexts.len()];
            let x_dropoff_per_context = vec![16; contexts.len()];
            let mut diag_array = vec![DiagStruct::default(); size as usize];
            let mut diag_offset = window;
            let (reference, _, _) = masked_reference(
                &lookup,
                &contexts,
                &s_frames,
                &cutoff_scores,
                &x_dropoff_per_context,
                window,
                &mut diag_array,
                &mut diag_offset,
            );
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
            let plain = XPairPar {
                groups: build_groups(&lookup, &contexts),
                stub: false,
                sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
                sliced: None,
                mask_rule: None,
                extents: Vec::new(),
            };
            // Ratio 1: every unit whose unmasked table is at least the masked one's size.
            let pp = XPairPar {
                groups: build_groups(&lookup, &contexts),
                stub: false,
                sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
                sliced: None,
                mask_rule: Some(XMaskRule {
                    max_cells: usize::MAX,
                    ratio: 1,
                    band_cells: None,
                }),
                extents: x_group_extents(&contexts),
            };
            assert_eq!(pp.extents.len(), pp.groups.len());
            // Negative control: tables of half the masked size (below extent + window, against
            // the module doc's condition), where the units may differ.
            let control = XPairPar {
                groups: build_groups(&lookup, &contexts),
                stub: false,
                sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
                sliced: None,
                mask_rule: pp.mask_rule,
                extents: pp
                    .extents
                    .iter()
                    .map(|&e| (x_masked_cells(e, window) / 2).saturating_sub(window as usize))
                    .collect(),
            };
            let chunk = SubjectChunk {
                offset: 0,
                length: s_len,
                overlap: 0,
            };
            let mut scratch = XScratch::new(false);
            let mut aliased = false;
            for g in 0..pp.groups.len() {
                if pp.groups[g].lookup.longest_chain <= 0 {
                    continue;
                }
                let unit = XUnit {
                    f: 0,
                    k: 0,
                    chunk,
                    g,
                };
                let mask = pp.x_unit_mask(&unit, window).expect("a masked table");
                assert!(plain.x_unit_mask(&unit, window).is_none());
                masked_units += 1;
                // The masked table covers the group's extent plus the window (module doc).
                assert!(mask as usize + 1 >= pp.extents[g] + window as usize);
                exact_fits += usize::from(mask as usize + 1 == pp.extents[g] + window as usize);
                wide_groups += usize::from(group_starts(&contexts).len() < contexts.len());
                let mut want: Vec<InitHSP> = Vec::new();
                plain.scan_unit_hits::<_, XNoProbe, false, false, false>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut want,
                    &mut XNoProbe,
                );
                saved += want.len();
                // Two saved extensions of the unit on raw diagonals that share a masked cell.
                let diagonals: Vec<i64> = want
                    .iter()
                    .map(|h| h.q_seed_absolute as i64 - h.s_seed as i64)
                    .collect();
                let cells = mask as i64 + 1;
                aliased |= diagonals.iter().any(|&a| {
                    diagonals
                        .iter()
                        .any(|&b| a != b && (a - b).rem_euclid(cells) == 0)
                });
                // LOSAT_X_TBXUNGAPSENT: the subject is one whole frame with its sentinels, so its
                // units take the sentinel walks whenever the batch does.
                let sent = pp.x_unit_sent(&unit, &inp);
                assert_eq!(sent, pp.sent);
                sent_units += usize::from(sent);
                let mut lists: Vec<(&str, Vec<InitHSP>)> = Vec::new();
                for name in [
                    "scan",
                    "scan pf",
                    "bf",
                    "bf pf",
                    "scan sent",
                    "bf sent",
                    "bf pf sent",
                ] {
                    if name.ends_with("sent") && !sent {
                        continue;
                    }
                    let mut got: Vec<InitHSP> = Vec::new();
                    match name {
                        "scan" => pp.scan_unit_hits::<_, XNoProbe, false, false, false>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            &mut XNoProbe,
                        ),
                        "scan pf" => pp.scan_unit_hits::<_, XNoProbe, true, false, false>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            &mut XNoProbe,
                        ),
                        "bf" => pp.scan_unit_hits::<_, XNoProbe, false, true, false>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            &mut XNoProbe,
                        ),
                        "bf pf" => pp.scan_unit_hits::<_, XNoProbe, true, true, false>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            &mut XNoProbe,
                        ),
                        "scan sent" => pp.scan_unit_hits::<_, XNoProbe, false, false, true>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            &mut XNoProbe,
                        ),
                        "bf sent" => pp.scan_unit_hits::<_, XNoProbe, false, true, true>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            &mut XNoProbe,
                        ),
                        _ => pp.scan_unit_hits::<_, XNoProbe, true, true, true>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            &mut XNoProbe,
                        ),
                    }
                    assert_eq!(scratch.diag.len(), mask as usize + 1);
                    lists.push((name, got));
                }
                for (name, len) in [("batched", x_tbx_extbatch::BATCH), ("batched 1", 1)] {
                    let mut got: Vec<InitHSP> = Vec::new();
                    pp.scan_unit_batched_hits::<_, false, false, false>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        len,
                    );
                    lists.push((name, got));
                    let mut got: Vec<InitHSP> = Vec::new();
                    pp.scan_unit_batched_hits::<_, true, true, false>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        len,
                    );
                    lists.push((name, got));
                    if sent {
                        let mut got: Vec<InitHSP> = Vec::new();
                        pp.scan_unit_batched_hits::<_, true, true, true>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            len,
                        );
                        lists.push((name, got));
                        let mut got: Vec<InitHSP> = Vec::new();
                        pp.scan_unit_batched_hits::<_, false, false, true>(
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut got,
                            len,
                        );
                        lists.push((name, got));
                    }
                }
                for (name, got) in &lists {
                    assert_eq!(
                        got, &want,
                        "case {case_no} group {g}: {name} on the masked table differs"
                    );
                }
                if control.x_unit_mask(&unit, window).map(|m| m as usize + 1)
                    == Some((mask as usize + 1) / 2)
                {
                    let mut got: Vec<InitHSP> = Vec::new();
                    control.scan_unit_hits::<_, XNoProbe, false, false, false>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut got,
                        &mut XNoProbe,
                    );
                    control_differs += usize::from(got != want);
                }
            }
            assert_eq!(
                unit_lists(&pp, &inp),
                reference,
                "case {case_no}: masked unit lists differ from the reference"
            );
            assert_eq!(
                unit_lists_batched(&pp, &inp),
                reference,
                "case {case_no}: masked batched unit lists differ from the reference"
            );
            aliased_cases += usize::from(aliased);
        }
        assert!(masked_units > 600, "masked_units {masked_units}");
        assert_eq!(sent_units, if on { masked_units } else { 0 });
        assert!(aliased_cases > 100, "aliased_cases {aliased_cases}");
        assert!(saved > 10_000, "saved {saved}");
        assert!(exact_fits > 50, "exact_fits {exact_fits}");
        assert!(wide_groups > 100, "wide_groups {wide_groups}");
        assert!(control_differs > 10, "control_differs {control_differs}");
    }

    // No NCBI counterpart: the sizes of LOSAT_X_TBXUNITMASK (module doc).
    #[test]
    fn masked_table_sizes_and_rule() {
        assert_eq!(x_masked_cells(0, 40), 64);
        assert_eq!(x_masked_cells(24, 40), 64);
        assert_eq!(x_masked_cells(25, 40), 128);
        assert_eq!(x_masked_cells(1 << 20, 0), 1 << 20);
        assert_eq!(x_masked_cells((1 << 20) - 40, 40), 1 << 20);
        assert_eq!(x_masked_cells((1 << 20) - 39, 40), 1 << 21);
        // Small masked tables where the unmasked one has at least `ratio` times their cells; in a
        // band batch, every unit whose chunk is longer than the reference table.
        let rule = |max_cells: usize, band_cells: Option<i32>| XMaskRule {
            max_cells,
            ratio: 4,
            band_cells,
        };
        assert_eq!(x_unit_mask_for(&rule(64, None), 24, 40, 256, 200), Some(63));
        assert_eq!(x_unit_mask_for(&rule(64, None), 24, 40, 255, 200), None);
        assert_eq!(x_unit_mask_for(&rule(63, None), 24, 40, 256, 200), None);
        assert_eq!(
            x_unit_mask_for(&rule(128, None), 25, 40, 512, 200),
            Some(127)
        );
        assert_eq!(x_unit_mask_for(&rule(128, None), 25, 40, 511, 200), None);
        assert_eq!(x_unit_mask_for(&rule(0, Some(256)), 24, 40, 300, 256), None);
        assert_eq!(
            x_unit_mask_for(&rule(0, Some(256)), 24, 40, 300, 257),
            Some(63)
        );
        assert_eq!(
            x_unit_mask_for(&rule(0, None), 24, 40, 1 << 20, 1 << 20),
            None
        );
        // The size limit by thread count: none at one, 4,096 cells at 2-7, the switch's KiB at 8+.
        assert_eq!(x_mask_max_cells(256, 1), 0);
        assert_eq!(x_mask_max_cells(256, 2), 4096);
        assert_eq!(x_mask_max_cells(256, 7), 4096);
        assert_eq!(x_mask_max_cells(256, 8), 65536);
        assert_eq!(x_mask_max_cells(2, 4), 512);
        assert_eq!(x_mask_max_cells(2, 16), 512);
        // The groups: one per context up to X_MAX_GROUPS, from a context's first offset to the
        // next group's (sentinels included), the last to the end of its last context.
        let frames: Vec<QueryFrame> = [5usize, 7, 9]
            .iter()
            .zip([1i8, 2, 3])
            .map(|(&len, f)| frame(f, &vec![aa("A")[0]; len]))
            .collect();
        let (_, contexts) = lookup_for(&frames);
        let extents = x_group_extents(&contexts);
        assert_eq!(extents.len(), 3);
        assert_eq!(
            extents[0],
            (contexts[1].frame_base - contexts[0].frame_base) as usize
        );
        assert_eq!(
            extents[1],
            (contexts[2].frame_base - contexts[1].frame_base) as usize
        );
        assert_eq!(extents[2], 9);
        assert!(extents[0] > 5 && extents[1] > 7);
    }

    /// The unit lists of `inp` with the batched extensions of LOSAT_X_TBXEXTBATCH (`unit_lists`
    /// with `run_unit_batched`).
    fn unit_lists_batched(pp: &XPairPar, inp: &XSubject<'_>) -> Vec<Vec<InitHSP>> {
        let mut out = Vec::new();
        let mut scratch = XScratch::new(false);
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
                        init_hsps.extend(pp.run_unit_batched::<false>(&unit, inp, &mut scratch));
                    }
                }
            }
            sort_init_hsps_by_score_ncbi(&mut init_hsps);
            out.push(init_hsps);
        }
        out
    }

    /// A frame with mask intervals (start, end exclusive, as SEG gives them); the lookup locations
    /// are their complement.
    fn masked_frame(f: i8, residues: &[u8], masks: &[(usize, usize)]) -> QueryFrame {
        let mut masked = frame(f, residues);
        masked.seg_masks = masks.to_vec();
        masked
    }

    /// A random frame: now and then empty or shorter than a word, with zero to two mask intervals
    /// inside it, now and then masked from end to end.
    fn random_masked_frame(rng: &mut Rng, f: i8, max_len: u64, alphabet: &[u8]) -> QueryFrame {
        let len = match rng.below(8) {
            0 => 0,
            1 => 2,
            _ => 3 + rng.below(max_len) as usize,
        };
        let residues = random_residues(rng, len, alphabet);
        let mut masks: Vec<(usize, usize)> = Vec::new();
        if len > 0 {
            match rng.below(5) {
                0 => masks.push((0, len)),
                1 => {}
                _ => {
                    let a = rng.below(len as u64) as usize;
                    let b = a + 1 + rng.below((len - a) as u64) as usize;
                    masks.push((a, b));
                    if b + 2 < len && rng.below(2) == 0 {
                        let c = b + 1 + rng.below((len - b - 1) as u64) as usize;
                        let d = c + 1 + rng.below((len - c) as u64) as usize;
                        masks.push((c, d));
                    }
                }
            }
        }
        masked_frame(f, &residues, &masks)
    }

    /// The direct groups of `queries` (`lut_direct_groups`, then the stub) against `build_groups`
    /// of the reference builder's table, with `lut_shadow_compare` (every
    /// field of every group table, `q_min`, `q_max`, the stub's header fields and `longest_chain`,
    /// the contexts). Returns the number of groups and of groups with entries.
    fn lut_check(
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> (usize, usize) {
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (full, contexts) = build_ncbi_lookup(queries, threshold, &params, true);
        let reference = build_groups(&full, &contexts);
        let q = XPairLutQuery::x_prepare(queries, threshold, true).expect("direct builder applies");
        // On the pool when the caller has one, else in order on this thread.
        let direct = lut_direct_groups(&q, pool, true).expect("every location inside one group");
        let stub = q.x_into_stub(direct.iter().map(|g| &g.lookup));
        let lut = XLutBatch {
            groups: direct,
            sliced: None,
            stub: Some(stub),
        };
        lut_shadow_compare(&lut, &reference, None, &full, &contexts);
        (
            reference.len(),
            reference
                .iter()
                .filter(|g| g.lookup.longest_chain > 0)
                .count(),
        )
    }

    const FRAMES: [i8; 6] = [1, 2, 3, -1, -2, -3];

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:454-456,461-465
    // ```c
    //     BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
    //                                       lookup->charsize, lookup->word_length,
    //                                       query, location);
    // ...
    //     for (i = 0; i < lookup->backbone_size; i++) {
    //         if (exact_backbone[i] != NULL) {
    //             s_AddWordHits(lookup, matrix, query->sequence,
    //                           exact_backbone[i], query_bias, row_max);
    //             sfree(exact_backbone[i]);
    // ```
    #[test]
    fn direct_groups_equal_build_groups_for_one_query() {
        let mut rng = Rng(0x0dd5_eed5_1234_5678);
        let alphabets = [aa("ACDEKW"), aa("AW"), aa("ACDEFGHIKLMNPQRSTVWY")];
        let serial = crate::utils::threading::SearchPool::x_serial();
        let mut with_entries = 0usize;
        let mut empty_groups = 0usize;
        for case_no in 0..90 {
            let alphabet = &alphabets[case_no % alphabets.len()];
            // 1 to 6 contexts: one group per context.
            let n_ctx = 1 + case_no % 6;
            let frames: Vec<QueryFrame> = (0..n_ctx)
                .map(|c| random_masked_frame(&mut rng, FRAMES[c], 300, alphabet))
                .collect();
            for threshold in [0, 11, 13] {
                let (groups, filled) = lut_check(std::slice::from_ref(&frames), threshold, &serial);
                assert_eq!(groups, n_ctx);
                with_entries += filled;
                empty_groups += groups - filled;
            }
        }
        assert!(with_entries > 300, "groups with entries {with_entries}");
        assert!(empty_groups > 100, "empty groups {empty_groups}");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:580-590
    // ```c
    //         for (i = 0; i < alphabet_size; i++) {
    //             if (score + row[i] >= threshold) {
    //                 subject_word[current_pos] = i;
    //                 for (j = 0; j < offset_list[1]; j++) {
    //                     BlastLookupAddWordHit(lookup->thin_backbone, wordsize,
    //                                           charsize, subject_word,
    //                                           query_bias + offset_list[j + 2]);
    //                 }
    // ```
    // A threshold of INT_MIN (the conversion of a threshold beyond Int4): every word is a
    // neighbour of every query word, so every cell of a group holds every offset of the group.
    #[test]
    fn direct_groups_equal_build_groups_at_threshold_int_min() {
        let mut rng = Rng(0x7357_0000_0000_0001);
        let alphabet = aa("ACDEFGHIKLMNPQRSTVWY");
        let serial = crate::utils::threading::SearchPool::x_serial();
        for n_ctx in [1usize, 4, 6] {
            let frames: Vec<QueryFrame> = (0..n_ctx)
                .map(|c| random_masked_frame(&mut rng, FRAMES[c], 9, &alphabet))
                .collect();
            let (groups, _) = lut_check(std::slice::from_ref(&frames), i32::MIN, &serial);
            assert_eq!(groups, n_ctx);
        }
        // Two queries: twelve contexts in six groups.
        let queries: Vec<Vec<QueryFrame>> = (0..2)
            .map(|_| {
                FRAMES
                    .iter()
                    .map(|&f| random_masked_frame(&mut rng, f, 6, &alphabet))
                    .collect()
            })
            .collect();
        let (groups, _) = lut_check(&queries, i32::MIN, &serial);
        assert!((2..=X_MAX_GROUPS).contains(&groups), "groups {groups}");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1045-1047
    // ```c
    //       start_offset = query_info->contexts[context].query_offset;
    //       end_offset = query_info->contexts[context].query_length
    //           + start_offset - 1;
    // ```
    // Several queries make more than X_MAX_GROUPS contexts: a group then holds several contexts,
    // and its direct table is built from the lookup locations of all of them.
    #[test]
    fn direct_groups_equal_build_groups_for_several_queries() {
        let mut rng = Rng(0x0bad_cafe_dead_beef);
        let alphabets = [aa("ACDEKW"), aa("ACDEFGHIKLMNPQRSTVWY")];
        let serial = crate::utils::threading::SearchPool::x_serial();
        let mut six_groups = 0usize;
        for case_no in 0..24 {
            let alphabet = &alphabets[case_no % alphabets.len()];
            let n_queries = 2 + case_no % 5;
            let queries: Vec<Vec<QueryFrame>> = (0..n_queries)
                .map(|_| {
                    FRAMES
                        .iter()
                        .map(|&f| random_masked_frame(&mut rng, f, 120, alphabet))
                        .collect()
                })
                .collect();
            for threshold in [0, 11, 13] {
                // 12 to 36 contexts in at most X_MAX_GROUPS groups (`group_starts` balances them
                // by length and can make fewer), so some group holds several contexts.
                let (groups, _) = lut_check(&queries, threshold, &serial);
                assert!((1..=X_MAX_GROUPS).contains(&groups), "groups {groups}");
                six_groups += usize::from(groups == X_MAX_GROUPS);
            }
        }
        assert!(six_groups > 0, "cases with six groups {six_groups}");
        // Four queries with frames of one length: six groups of four contexts.
        let even: Vec<Vec<QueryFrame>> = (0..4)
            .map(|_| {
                FRAMES
                    .iter()
                    .map(|&f| masked_frame(f, &random_residues(&mut rng, 60, &alphabets[1]), &[]))
                    .collect()
            })
            .collect();
        assert_eq!(lut_check(&even, 13, &serial), (X_MAX_GROUPS, X_MAX_GROUPS));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:96-100
    // ```c
    //         /* if this location is too small to fit a complete word, skip the
    //            location */
    //
    //         if (word_length > to - from + 1)
    //             continue;
    // ```
    // A context masked from end to end has no lookup location, an empty frame has no letter and a
    // two-letter frame no word: their groups have no entry, `q_min` = `q_max` = 0, and the units
    // skip them (`longest_chain` 0).
    #[test]
    fn direct_groups_of_contexts_without_words_are_empty() {
        let mut rng = Rng(0x00c0_ffee_0000_0042);
        let alphabet = aa("ACDEKLMNPQRSTVWY");
        let serial = crate::utils::threading::SearchPool::x_serial();
        let frames = vec![
            masked_frame(1, &random_residues(&mut rng, 80, &alphabet), &[(10, 20)]),
            masked_frame(2, &random_residues(&mut rng, 70, &alphabet), &[(0, 70)]),
            masked_frame(3, &[], &[]),
            masked_frame(-1, &random_residues(&mut rng, 2, &alphabet), &[]),
            masked_frame(-2, &random_residues(&mut rng, 60, &alphabet), &[]),
            masked_frame(
                -3,
                &random_residues(&mut rng, 50, &alphabet),
                &[(0, 3), (47, 50)],
            ),
        ];
        for threshold in [0, 11, 13] {
            assert_eq!(
                lut_check(std::slice::from_ref(&frames), threshold, &serial),
                (6, 3)
            );
        }
        let q = XPairLutQuery::x_prepare(std::slice::from_ref(&frames), 13, true)
            .expect("direct builder applies");
        let groups =
            lut_direct_groups(&q, &serial, false).expect("every location inside one group");
        for g in [1usize, 2, 3] {
            assert_eq!(groups[g].lookup.longest_chain, 0);
            assert_eq!((groups[g].q_min, groups[g].q_max), (0, 0));
            assert!(groups[g].lookup.overflow.is_empty());
        }
        // The stub has the full table's longest chain although it has no entry.
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (full, _) = build_ncbi_lookup(std::slice::from_ref(&frames), 13, &params, true);
        let (stub, contexts) = q.x_into_stub(groups.iter().map(|g| &g.lookup));
        assert!(full.longest_chain > 0);
        assert_eq!(stub.longest_chain, full.longest_chain);
        assert_eq!(contexts.len(), 6);
        assert!(stub.pv.iter().all(|&word| word == 0) && stub.overflow.is_empty());
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149
    // ```c
    // TBlastThreads the_threads(GetNumberOfThreads());
    // ```
    // The groups built on a search pool (every thread with its own neighbour-list memo) are the
    // groups built in order on one thread.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    #[test]
    fn direct_groups_on_a_search_pool_equal_build_groups() {
        let mut rng = Rng(0x9001_9002_9003_9004);
        let alphabet = aa("ACDEFGHIKLMNPQRSTVWY");
        let one: Vec<Vec<QueryFrame>> = vec![FRAMES
            .iter()
            .map(|&f| random_masked_frame(&mut rng, f, 2000, &alphabet))
            .collect()];
        let many: Vec<Vec<QueryFrame>> = (0..4)
            .map(|_| {
                FRAMES
                    .iter()
                    .map(|&f| random_masked_frame(&mut rng, f, 400, &alphabet))
                    .collect()
            })
            .collect();
        crate::utils::threading::with_search_pool(4, "tblastx", |pool| {
            assert!(pool.enabled());
            for threshold in [0, 11, 13] {
                // The second round reads the lists the threads memoised in the first.
                for _round in 0..2 {
                    assert_eq!(lut_check(&one, threshold, pool).0, 6);
                    let (groups, _) = lut_check(&many, threshold, pool);
                    assert!((2..=X_MAX_GROUPS).contains(&groups), "groups {groups}");
                }
            }
            Ok(())
        })
        .expect("search pool");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-57
    // ```c
    // diag_array_length = 1;
    // /* What power of 2 is just longer than the query? */
    // while (diag_array_length < (qlen+window_size))
    // {
    //         diag_array_length = diag_array_length << 1;
    // }
    // ```
    // `decide` is the decision of `prepare`: the gates, then the guard.
    #[test]
    fn decide_checks_the_gates_and_the_guard() {
        let mut rng = Rng(5);
        let alphabet = aa("ACDEKW");
        let window = 40;
        // One context of 216 letters: the table has 256 cells and the slack is 40.
        let (_, contexts) = lookup_for(&[frame(1, &random_residues(&mut rng, 216, &alphabet))]);
        let query_length = query_length_of(&contexts);
        let diag_array_size = diag_array_size_for(query_length, window);
        assert_eq!((query_length, diag_array_size), (216, 256));
        let short = [FastaRecord::new("", b"s", &[b'A'; 120])];
        let long = [
            FastaRecord::new("", b"s", &[b'A'; 30]),
            FastaRecord::new("", b"l", &[b'A'; 123]),
        ];
        let none: [FastaRecord; 0] = [];
        let decision_band =
            |subjects: &[FastaRecord], debugging: bool, chunk_switches: bool, band: bool| {
                decide(&XGates {
                    contexts: &contexts,
                    query_length,
                    diag_array_size,
                    subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(subjects),
                    debugging,
                    chunk_switches,
                    threads: 1,
                    cost_gate: false,
                    n1_gate: false,
                    threshold: 13,
                    band,
                    unit_mask: false,
                })
            };
        let decision = |subjects: &[FastaRecord], debugging: bool, chunk_switches: bool| {
            decision_band(subjects, debugging, chunk_switches, false)
        };
        assert_eq!(decision(&short, false, false), XDecision::Units);
        assert_eq!(decision(&none, false, false), XDecision::Units);
        assert_eq!(decision(&long, false, false), XDecision::Guard);
        assert_eq!(decision(&short, true, false), XDecision::Gate);
        assert_eq!(decision(&short, false, true), XDecision::Gate);
        // The gates come before the guard, as in the counters of LOSAT_X_STATS.
        assert_eq!(decision(&long, true, false), XDecision::Gate);
        // LOSAT_X_TBXBAND: the first guard is not applied; the gates still are.
        assert_eq!(decision_band(&long, false, false, true), XDecision::Units);
        assert_eq!(decision_band(&short, false, false, true), XDecision::Units);
        assert_eq!(decision_band(&long, true, false, true), XDecision::Gate);
        assert_eq!(decision_band(&short, false, true, true), XDecision::Gate);
        // The memory rule: a subject chunk longer than the table's 256 cells keeps the guard.
        let longer = [FastaRecord::new("", b"l", &[b'A'; 3 * 257])];
        let fits = [FastaRecord::new("", b"l", &[b'A'; 3 * 256])];
        assert_eq!(decision_band(&longer, false, false, true), XDecision::Guard);
        assert_eq!(decision_band(&fits, false, false, true), XDecision::Units);
        assert!(x_band_fits(256, 256) && !x_band_fits(257, 256));
        // LOSAT_X_TBXUNITMASK: the memory rule is not applied where LOSAT_X_TBXQSLICE cuts no group
        // (small groups at one thread); the gates and the second guard still are.
        let decision_mask = |subjects: &[FastaRecord], debugging: bool| {
            decide(&XGates {
                contexts: &contexts,
                query_length,
                diag_array_size,
                subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(subjects),
                debugging,
                chunk_switches: false,
                threads: 1,
                cost_gate: false,
                n1_gate: false,
                threshold: 13,
                band: true,
                unit_mask: true,
            })
        };
        assert_eq!(decision_mask(&longer, false), XDecision::Units);
        assert_eq!(decision_mask(&fits, false), XDecision::Units);
        assert_eq!(decision_mask(&longer, true), XDecision::Gate);
        // LOSAT_X_PAIRPARGATE still decides after the guard: this table of 256 cells at one thread.
        let gated = |band: bool| {
            decide(&XGates {
                contexts: &contexts,
                query_length,
                diag_array_size,
                subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(&long),
                debugging: false,
                chunk_switches: false,
                threads: 1,
                cost_gate: true,
                n1_gate: false,
                threshold: 13,
                band,
                unit_mask: false,
            })
        };
        assert_eq!(gated(true), XDecision::Cost);
        assert_eq!(gated(false), XDecision::Guard);
        // The second guard (`last_hit_fits`): with a window just under 2^29 the slack is almost
        // 2^29, so the first guard holds for the short subject and the units still do not run.
        let wide = diag_array_size_for(query_length, 536_870_000);
        assert_eq!(wide, 1 << 29);
        assert!(guard(216, 40, query_length, wide));
        let gates = XGates {
            contexts: &contexts,
            query_length,
            diag_array_size: wide,
            subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(&short),
            debugging: false,
            chunk_switches: false,
            threads: 1,
            cost_gate: false,
            n1_gate: false,
            threshold: 13,
            band: false,
            unit_mask: false,
        };
        assert_eq!(decide(&gates), XDecision::Guard);
        // LOSAT_X_TBXBAND does not lift the second guard.
        assert_eq!(
            decide(&XGates {
                band: true,
                ..gates
            }),
            XDecision::Guard
        );
    }

    // No NCBI counterpart: LOSAT_X_PAIRPARGATE chooses between two LOSAT paths with the same result.
    // The cost gate comes after the gates and the guard, keeps the units below
    // the default threshold, for a reference table of 2^18 cells or more, and at several threads
    // for a subject that holds more than one and a half thread's share of the bases.
    #[test]
    fn the_cost_gate_keeps_the_units_where_they_gain() {
        let mut rng = Rng(7);
        let alphabet = aa("ACDEKW");
        let window = 40;
        let (_, contexts) = lookup_for(&[frame(1, &random_residues(&mut rng, 216, &alphabet))]);
        let query_length = query_length_of(&contexts);
        let diag_array_size = diag_array_size_for(query_length, window);
        assert_eq!(diag_array_size, 256);
        let one = [FastaRecord::new("", b"s", &[b'A'; 120])];
        // 120 of 200 bases in one subject: 2 * 120 * t > 3 * 200 from t = 3.
        let mixed = [
            FastaRecord::new("", b"s", &[b'A'; 120]),
            FastaRecord::new("", b"t", &[b'A'; 40]),
            FastaRecord::new("", b"u", &[b'A'; 40]),
        ];
        let none: [FastaRecord; 0] = [];
        let long = [
            FastaRecord::new("", b"s", &[b'A'; 30]),
            FastaRecord::new("", b"l", &[b'A'; 123]),
        ];
        let decision_at = |threshold: i32,
                           diag_array_size: i32,
                           subjects: &[FastaRecord],
                           debugging: bool,
                           threads: usize,
                           cost_gate: bool| {
            decide(&XGates {
                contexts: &contexts,
                query_length,
                diag_array_size,
                subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(subjects),
                debugging,
                chunk_switches: false,
                threads,
                cost_gate,
                n1_gate: false,
                threshold,
                band: false,
                unit_mask: false,
            })
        };
        let decision = |diag_array_size: i32,
                        subjects: &[FastaRecord],
                        debugging: bool,
                        threads: usize,
                        cost_gate: bool| {
            decision_at(13, diag_array_size, subjects, debugging, threads, cost_gate)
        };
        // Without the gate nothing changes.
        for t in [1, 2, 8] {
            assert_eq!(decision(256, &one, false, t, false), XDecision::Units);
        }
        // One thread: the table size alone.
        assert_eq!(decision(256, &one, false, 1, true), XDecision::Cost);
        assert_eq!(decision(1 << 17, &one, false, 1, true), XDecision::Cost);
        assert_eq!(decision(1 << 18, &one, false, 1, true), XDecision::Units);
        assert_eq!(decision(256, &none, false, 1, true), XDecision::Cost);
        // Several threads: one subject is all of the bases.
        assert_eq!(decision(256, &one, false, 2, true), XDecision::Units);
        assert_eq!(decision(256, &none, false, 2, true), XDecision::Cost);
        assert_eq!(decision(256, &mixed, false, 2, true), XDecision::Cost);
        assert_eq!(decision(256, &mixed, false, 3, true), XDecision::Units);
        assert_eq!(decision(256, &mixed, false, 8, true), XDecision::Units);
        assert_eq!(decision(1 << 18, &mixed, false, 2, true), XDecision::Units);
        // Below the default threshold (and INT_MIN, `-threshold +inf`) the units stay.
        assert_eq!(decision_at(12, 256, &one, false, 1, true), XDecision::Units);
        assert_eq!(
            decision_at(i32::MIN, 256, &one, false, 1, true),
            XDecision::Units
        );
        assert_eq!(decision_at(14, 256, &one, false, 1, true), XDecision::Cost);
        // The gates and the guard come first.
        assert_eq!(decision(256, &one, true, 1, true), XDecision::Gate);
        assert_eq!(decision(256, &long, false, 1, true), XDecision::Guard);
    }

    // No NCBI counterpart: LOSAT_X_PAIRPARN1 moves the one-thread line of LOSAT_X_PAIRPARGATE from
    // 2^18 cells to batches of 5,000 letters; the rule of several threads, the threshold rule and
    // the gates before them are unchanged.
    #[test]
    fn the_one_thread_line_keeps_the_units_from_5000_letters() {
        let mut rng = Rng(11);
        let alphabet = aa("ACDEKW");
        let (_, contexts) = lookup_for(&[frame(1, &random_residues(&mut rng, 216, &alphabet))]);
        let one = [FastaRecord::new("", b"s", &[b'A'; 120])];
        // 10 subjects of 40 bases: no subject holds one and a half thread's share at 2 or 8 threads.
        let spread: Vec<FastaRecord> = (0..10)
            .map(|_| FastaRecord::new("", b"s", &[b'A'; 40]))
            .collect();
        // The batch's length in letters (the table follows it, as at the build site).
        let decision = |query_length: i32,
                        subjects: &[FastaRecord],
                        threads: usize,
                        threshold: i32,
                        n1_gate: bool,
                        debugging: bool| {
            decide(&XGates {
                contexts: &contexts,
                query_length,
                diag_array_size: diag_array_size_for(query_length, 40),
                subjects: crate::algorithm::tblastx::subjects::Subjects::Fasta(subjects),
                debugging,
                chunk_switches: false,
                threads,
                cost_gate: true,
                n1_gate,
                threshold,
                band: false,
                unit_mask: false,
            })
        };
        assert_eq!(diag_array_size_for(140_000, 40), 1 << 18);
        // One thread: the units from 5,000 letters with the switch, from 2^18 cells without it.
        assert_eq!(decision(4_999, &one, 1, 13, true, false), XDecision::Cost);
        assert_eq!(decision(5_000, &one, 1, 13, true, false), XDecision::Units);
        assert_eq!(
            decision(100_000, &spread, 1, 13, true, false),
            XDecision::Units
        );
        assert_eq!(decision(5_000, &one, 1, 13, false, false), XDecision::Cost);
        assert_eq!(
            decision(100_000, &one, 1, 13, false, false),
            XDecision::Cost
        );
        assert_eq!(
            decision(140_000, &one, 1, 13, false, false),
            XDecision::Units
        );
        // Several threads keep their rule: spread subjects below 2^18 cells take the reference loop.
        for threads in [2, 8] {
            assert_eq!(
                decision(5_000, &spread, threads, 13, true, false),
                XDecision::Cost
            );
            assert_eq!(
                decision(100_000, &spread, threads, 13, true, false),
                XDecision::Cost
            );
            assert_eq!(
                decision(5_000, &one, threads, 13, true, false),
                XDecision::Units
            );
        }
        // The threshold rule and the gates come first.
        assert_eq!(decision(4_999, &one, 1, 12, true, false), XDecision::Units);
        assert_eq!(decision(5_000, &one, 1, 13, true, true), XDecision::Gate);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:264-266
    // ```c
    //    case eAaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastAaLookupTable*)lookup->lut)->longest_chain;
    // ```
    // The build site gets the stub once.
    #[test]
    fn the_stub_is_taken_once() {
        let mut rng = Rng(9);
        let alphabet = aa("ACDEKW");
        let serial = crate::utils::threading::SearchPool::x_serial();
        let queries = vec![vec![
            frame(1, &random_residues(&mut rng, 50, &alphabet)),
            frame(2, &random_residues(&mut rng, 40, &alphabet)),
        ]];
        let batch = || {
            let q = XPairLutQuery::x_prepare(&queries, 11, true).expect("direct builder applies");
            let groups =
                lut_direct_groups(&q, &serial, false).expect("every location inside one group");
            let stub = q.x_into_stub(groups.iter().map(|g| &g.lookup));
            XLutBatch {
                groups,
                sliced: None,
                stub: Some(stub),
            }
        };
        let mut direct = batch();
        let (stub, contexts) = direct.take_stub().expect("the stub of the batch");
        assert_eq!(contexts.len(), 2);
        assert_eq!(stub.num_contexts, 2);
        assert!(direct.take_stub().is_none());
        assert!(direct.stub.is_none());
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

    /// One random case of the LOSAT_X_SEEDCOMPACT tests: one query of one to six contexts (now
    /// and then two or three queries, so that a group holds several contexts), its groups, and one
    /// to six subject frames of `s_min` to `s_min + s_spread - 1` letters (now and then shorter
    /// than a word).
    struct CompactCase {
        pp: XPairPar,
        contexts: Vec<QueryContext>,
        s_frames: Vec<QueryFrame>,
        cutoff_scores: Vec<i32>,
        x_dropoff_per_context: Vec<i32>,
    }

    impl CompactCase {
        fn random(rng: &mut Rng, alphabet: &[u8], s_min: usize, s_spread: u64) -> Self {
            let n_queries = if rng.below(4) == 0 {
                2 + rng.below(2) as usize
            } else {
                1
            };
            let mut queries: Vec<Vec<QueryFrame>> = Vec::new();
            for _ in 0..n_queries {
                let n_ctx = 1 + rng.below(6) as usize;
                let mut frames = Vec::new();
                for &f in &FRAMES[..n_ctx] {
                    let len = 3 + rng.below(60) as usize;
                    frames.push(frame(f, &random_residues(rng, len, alphabet)));
                }
                queries.push(frames);
            }
            let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
            let (lookup, contexts) = build_ncbi_lookup(&queries, 11, &params, true);
            let n_s = 1 + rng.below(6) as usize;
            let mut s_frames = Vec::new();
            for &f in &FRAMES[..n_s] {
                let len = if rng.below(10) == 0 {
                    rng.below(3) as usize
                } else {
                    s_min + rng.below(s_spread) as usize
                };
                s_frames.push(frame(f, &random_residues(rng, len, alphabet)));
            }
            Self {
                pp: XPairPar {
                    groups: build_groups(&lookup, &contexts),
                    stub: false,
                    sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
                    sliced: None,
                    mask_rule: None,
                    extents: Vec::new(),
                },
                cutoff_scores: vec![12; contexts.len()],
                x_dropoff_per_context: vec![16; contexts.len()],
                contexts,
                s_frames,
            }
        }

        fn inp(&self) -> XSubject<'_> {
            XSubject {
                s_frames_preliminary: &self.s_frames,
                s_frames: &self.s_frames,
                contexts: &self.contexts,
                cutoff_scores: &self.cutoff_scores,
                x_dropoff_per_context: &self.x_dropoff_per_context,
                window: 40,
                wordsize: 3,
                s_idx: 7,
                s_len: 4_321,
            }
        }
    }

    /// The chunks of every subject frame and the units of every (frame, chunk), listed with the
    /// statements of `run_subject_units` for chunks of at most `max_dbseq_len` letters
    /// (`seed_compact_test_listing_is_run_subject_units` checks it at MAX_DBSEQ_LEN).
    fn list_units(
        pp: &XPairPar,
        inp: &XSubject<'_>,
        max_dbseq_len: usize,
    ) -> (Vec<Vec<SubjectChunk>>, Vec<XUnit>) {
        let mut frame_chunks = Vec::new();
        let mut units = Vec::new();
        for (f, s_frame) in inp.s_frames_preliminary.iter().enumerate() {
            let mut split_state = SubjectSplitState::new(s_frame.aa_len);
            let mut chunks = Vec::new();
            while let SubjectChunkStatus::Ok(chunk) =
                split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP)
            {
                chunks.push(chunk);
            }
            for (k, chunk) in chunks.iter().enumerate() {
                if chunk.length < inp.wordsize as usize {
                    continue;
                }
                for (g, group) in pp.groups.iter().enumerate() {
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
        (frame_chunks, units)
    }

    /// What `compact_check` saw: units, units without a saved init HSP, compact records, frames
    /// whose list got hits from two or more chunks (the chunk merge ran), and hits of the subject.
    #[derive(Default)]
    struct CompactSeen {
        units: usize,
        empty_units: usize,
        records: usize,
        merged_frames: usize,
        hits: usize,
    }

    impl CompactSeen {
        fn add(&mut self, other: CompactSeen) {
            self.units += other.units;
            self.empty_units += other.empty_units;
            self.records += other.records;
            self.merged_frames += other.merged_frames;
            self.hits += other.hits;
        }
    }

    /// The units of `case` (chunks of at most `max_dbseq_len` letters) run on `pool` both ways:
    /// every unit's compact records have no spare room and expand to today's list of that unit,
    /// element for element, and `merge_subject_listed` of the compact records gives today's
    /// `combined_ungapped_hits`, every field (`shadow_compare`), as does
    /// the merge of LOSAT_X_TBXKMERGE of either form, sorted as the units' tasks sort them.
    fn compact_check(
        case: &CompactCase,
        pool: &crate::utils::threading::SearchPool<'_>,
        max_dbseq_len: usize,
    ) -> CompactSeen {
        let inp = case.inp();
        let (frame_chunks, units) = list_units(&case.pp, &inp, max_dbseq_len);
        let full = case
            .pp
            .run_units(&units, &inp, pool, XPairPar::run_unit::<false>, false);
        let compact = case.pp.run_units(
            &units,
            &inp,
            pool,
            XPairPar::run_unit_compact::<false>,
            false,
        );
        assert_eq!((full.len(), compact.len()), (units.len(), units.len()));
        let mut seen = CompactSeen {
            units: units.len(),
            ..CompactSeen::default()
        };
        for (unit, (list, records)) in units.iter().zip(full.iter().zip(&compact)) {
            let at = format!("unit f {} k {} g {}", unit.f, unit.k, unit.g);
            // No spare room, except in a copy that LOSAT_X_SEEDUNMAP gave the mapped capacity.
            let mapped_copy = records.len() * std::mem::size_of::<XUnitHsp>() >= X_MAPPED_MIN
                && records.capacity() == x_mapped_len::<XUnitHsp>(records.len());
            assert!(
                records.capacity() == records.len() || mapped_copy,
                "{at}: spare room"
            );
            let mut expanded = Vec::new();
            x_expand_unit(records, unit.f, &inp, &mut expanded);
            assert_eq!(&expanded, list, "{at}");
            seen.empty_units += usize::from(list.is_empty());
            seen.records += records.len();
        }
        for (f, chunks) in frame_chunks.iter().enumerate() {
            let chunks_with_hits = (0..chunks.len())
                .filter(|&k| {
                    units
                        .iter()
                        .zip(&full)
                        .any(|(u, list)| u.f == f && u.k == k && !list.is_empty())
                })
                .count();
            seen.merged_frames += usize::from(chunks_with_hits >= 2);
        }
        // LOSAT_X_TBXKMERGE (x_kmerge.rs): the units' lists sorted as their tasks sort them, merged
        // by the switch's merge; both forms give today's list.
        let mut scratch = XScratch::new(false);
        let mut sorted_full = full.clone();
        sorted_full
            .iter_mut()
            .for_each(|list| list.x_sort_in(&mut scratch));
        let mut sorted_compact = compact.clone();
        sorted_compact
            .iter_mut()
            .for_each(|records| records.x_sort_in(&mut scratch));
        let today = XPairPar::merge_subject_listed(
            XSubjectUnits {
                frame_chunks: frame_chunks.clone(),
                units: units.clone(),
                outputs: XUnitOutputs::Full(full),
            },
            &inp,
        );
        let with_switch = XPairPar::merge_subject_listed(
            XSubjectUnits {
                frame_chunks: frame_chunks.clone(),
                units: units.clone(),
                outputs: XUnitOutputs::Compact(compact),
            },
            &inp,
        );
        shadow_compare(inp.s_idx, &with_switch, &today);
        for outputs in [
            XUnitOutputs::Full(sorted_full),
            XUnitOutputs::Compact(sorted_compact),
        ] {
            let merged = x_kmerge::merge_subject(
                XSubjectUnits {
                    frame_chunks: frame_chunks.clone(),
                    units: units.clone(),
                    outputs,
                },
                &inp,
            );
            shadow_compare(inp.s_idx, &merged, &today);
        }
        seen.hits = today.len();
        seen
    }

    const COMPACT_ALPHABETS: [&str; 3] = ["ACDEKW", "AW", "ACDEFGHIKLMNPQRSTVWY"];

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:367-372
    // ```c
    //     ungapped_data->q_start = q_start;
    //     ungapped_data->s_start = s_start;
    //     ungapped_data->length = len;
    //     ungapped_data->score = score;
    //
    //     BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
    // ```
    // LOSAT_X_SEEDCOMPACT on random units, one chunk per frame, serially (one scratch for every
    // unit of a case): each unit's records are its init HSPs of today, and the merge is today's.
    #[test]
    fn seed_compact_unit_records_expand_to_todays_lists() {
        let mut rng = Rng(0x5eed_c0de_0000_0001);
        let serial = crate::utils::threading::SearchPool::x_serial();
        let mut seen = CompactSeen::default();
        for case_no in 0..150 {
            let alphabet = aa(COMPACT_ALPHABETS[case_no % COMPACT_ALPHABETS.len()]);
            let case = CompactCase::random(&mut rng, &alphabet, 3, 400);
            seen.add(compact_check(&case, &serial, MAX_DBSEQ_LEN));
        }
        assert!(seen.units > 1_000, "units {}", seen.units);
        assert!(seen.empty_units > 50, "empty units {}", seen.empty_units);
        assert!(seen.records > 2_000, "records {}", seen.records);
        assert!(seen.hits > 1_000, "hits {}", seen.hits);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2876-2881
    // ```c
    //    /* If no previous HSP list, just return a copy of the new one. */
    //    if (!combined_hsp_list) {
    //       *combined_hsp_list_ptr = hsp_list;
    //       *hsp_list_ptr = NULL;
    //       return 0;
    //    }
    // ```
    // Subjects of several chunks per frame (chunks of 103, 160 and 300 letters with the overlap
    // of 100: the first chunk's list moves into the empty frame list, the later ones go through
    // `merge_tblastx_subject_chunk_hits`): the compact merge gives today's list, every field.
    #[test]
    fn seed_compact_merge_with_several_chunks_per_frame_is_todays_merge() {
        let mut rng = Rng(0x5eed_c0de_0000_0002);
        let serial = crate::utils::threading::SearchPool::x_serial();
        let mut seen = CompactSeen::default();
        for case_no in 0..36 {
            let alphabet = aa(COMPACT_ALPHABETS[case_no % COMPACT_ALPHABETS.len()]);
            let case = CompactCase::random(&mut rng, &alphabet, 200, 500);
            let max_dbseq_len = [DBSEQ_CHUNK_OVERLAP + 3, 160, 300][case_no % 3];
            seen.add(compact_check(&case, &serial, max_dbseq_len));
        }
        assert!(
            seen.merged_frames > 20,
            "merged frames {}",
            seen.merged_frames
        );
        assert!(seen.hits > 1_000, "hits {}", seen.hits);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1043
    // ```c
    //     aux_struct->init_hitlist = BLAST_InitHitListNew();
    // ```
    // On search pools of 1, 2 and 8 threads (one scratch vector per rayon job, kept across its
    // units): the same records, lists and merge as today, with one and with several chunks per
    // frame.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    #[test]
    fn seed_compact_lists_on_search_pools_are_todays() {
        for threads in [1usize, 2, 8] {
            crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
                assert_eq!(pool.enabled(), threads > 1);
                let mut rng = Rng(0x5eed_c0de_0000_0003);
                let mut seen = CompactSeen::default();
                for case_no in 0..30 {
                    let alphabet = aa(COMPACT_ALPHABETS[case_no % COMPACT_ALPHABETS.len()]);
                    let (case, max_dbseq_len) = if case_no % 2 == 0 {
                        (
                            CompactCase::random(&mut rng, &alphabet, 3, 400),
                            MAX_DBSEQ_LEN,
                        )
                    } else {
                        (CompactCase::random(&mut rng, &alphabet, 200, 400), 160)
                    };
                    seen.add(compact_check(&case, pool, max_dbseq_len));
                }
                assert!(
                    seen.records > 500,
                    "threads {threads}: records {}",
                    seen.records
                );
                assert!(
                    seen.merged_frames > 5,
                    "threads {threads}: merged frames {}",
                    seen.merged_frames
                );
                Ok(())
            })
            .expect("search pool");
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2371-2375
    // ```c
    // s_GetUngappedHSPContext(const BlastQueryInfo* query_info,
    //                         const BlastInitHSP* init_hsp)
    // {
    //     return BSearchContextInfo(init_hsp->offsets.qs_offsets.q_off, query_info);
    // }
    // ```
    // One record, through the sink and back, for context indices below and above 16 bits (a
    // batch of 70,001 contexts with distinct query indices, frames, lengths and offsets); no
    // record expands to nothing. The record is 28 bytes.
    #[test]
    fn seed_compact_one_record_and_context_indices_above_16_bits_round_trip() {
        assert_eq!(std::mem::size_of::<XUnitHsp>(), 28);
        let (_, template) = lookup_for(&[frame(1, &aa("ACDEKWACDEKW"))]);
        let contexts: Vec<QueryContext> = (0..70_001usize)
            .map(|c| QueryContext {
                q_idx: (c / 6) as u32,
                frame: FRAMES[c % 6],
                orig_len: 100 + c,
                frame_base: (c * 13) as i32,
                ..template[0].clone()
            })
            .collect();
        let s_frames = vec![
            frame(2, &aa("ACDEKWACDEKW")),
            frame(-3, &aa("WKEDCAWKEDCA")),
        ];
        let inp = XSubject {
            s_frames_preliminary: &s_frames,
            s_frames: &s_frames,
            contexts: &contexts,
            cutoff_scores: &[],
            x_dropoff_per_context: &[],
            window: 40,
            wordsize: 3,
            s_idx: 3,
            s_len: 38,
        };
        for (ctx_idx, f) in [(0usize, 0usize), (65_535, 1), (65_536, 0), (70_000, 1)] {
            let ctx = &contexts[ctx_idx];
            // The InitHSP `x_hit` pushes for hsp_q 2, hsp_s 1, hsp_len 5 (its statements).
            let (hsp_q, hsp_s, hsp_len) = (2, 1, 5);
            let hsp = InitHSP {
                q_start_absolute: ctx.frame_base + hsp_q,
                q_end_absolute: ctx.frame_base + (hsp_q + hsp_len),
                s_start: hsp_s,
                s_end: hsp_s + hsp_len,
                q_seed_absolute: ctx.frame_base + 4,
                s_seed: 3,
                score: 31,
                ctx_idx,
                s_f_idx: f,
                q_idx: ctx.q_idx,
                s_idx: inp.s_idx as u32,
                q_frame: ctx.frame,
                s_frame: s_frames[f].frame,
                q_orig_len: ctx.orig_len,
                s_orig_len: inp.s_len,
            };
            let mut records: Vec<XUnitHsp> = Vec::new();
            XInitSink::push(&mut records, hsp);
            assert_eq!(records.len(), 1);
            assert_eq!(records[0].ctx_idx as usize, ctx_idx);
            let mut expanded = Vec::new();
            x_expand_unit(&records, f, &inp, &mut expanded);
            assert_eq!(expanded, vec![hsp], "context {ctx_idx}");
            // Today's sink keeps the InitHSP itself.
            let mut today: Vec<InitHSP> = Vec::new();
            XInitSink::push(&mut today, hsp);
            assert_eq!(today, vec![hsp]);
        }
        let mut untouched: Vec<InitHSP> = Vec::new();
        x_expand_unit(&[], 1, &inp, &mut untouched);
        assert!(untouched.is_empty());
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:479-491
    // ```c
    // status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
    //                                dbseq_chunk_overlap);
    // ...
    // BlastInitHitListReset(init_hitlist);
    // ```
    // The listing of the tests (`list_units`) is the listing of `run_subject_units` at
    // MAX_DBSEQ_LEN, and `run_subject_units` gives the form the switch of this process selects
    // (today's lists without LOSAT_X_SEEDCOMPACT), with the lists of `run_units` (sorted by the
    // units where LOSAT_X_TBXKMERGE is on, x_kmerge.rs item 1).
    #[test]
    fn seed_compact_test_listing_is_run_subject_units() {
        let mut rng = Rng(0x5eed_c0de_0000_0004);
        let serial = crate::utils::threading::SearchPool::x_serial();
        for case_no in 0..20 {
            let alphabet = aa(COMPACT_ALPHABETS[case_no % COMPACT_ALPHABETS.len()]);
            let case = CompactCase::random(&mut rng, &alphabet, 3, 300);
            let inp = case.inp();
            let x_units = case.pp.run_subject_units(&inp, &serial);
            let (frame_chunks, units) = list_units(&case.pp, &inp, tblastx_max_dbseq_len_for_run());
            assert_eq!(x_units.frame_chunks, frame_chunks);
            let key = |u: &XUnit| (u.f, u.k, u.chunk, u.g);
            assert_eq!(
                x_units.units.iter().map(key).collect::<Vec<_>>(),
                units.iter().map(key).collect::<Vec<_>>()
            );
            match (seed_compact(), x_units.outputs) {
                (false, XUnitOutputs::Full(outputs)) => {
                    assert_eq!(
                        outputs,
                        case.pp.run_units(
                            &units,
                            &inp,
                            &serial,
                            XPairPar::run_unit::<false>,
                            x_kmerge::on()
                        )
                    );
                }
                (true, XUnitOutputs::Compact(outputs)) => {
                    assert_eq!(
                        outputs,
                        case.pp.run_units(
                            &units,
                            &inp,
                            &serial,
                            XPairPar::run_unit_compact::<false>,
                            x_kmerge::on()
                        )
                    );
                }
                (on, _) => panic!("run_subject_units chose the other form (switch {on})"),
            }
        }
    }

    // No NCBI counterpart: LOSAT_X_SEEDUNMAP's rule at its boundaries: a query-group span or a
    // chunk of 65,535 and 65,536 letters, a long side with a short one (the 200 short queries
    // against AP027133, and the reverse), no overflow at the largest sizes; never with the switch
    // off, and only on Linux with glibc. Two minimum sides make a table of at least 512 KiB, above
    // `X_MAPPED_MIN`. The mapped capacities pass glibc's largest dynamic mmap threshold for both
    // element types.
    #[test]
    fn seed_unmap_rule_at_its_boundaries() {
        let glibc = cfg!(all(target_os = "linux", target_env = "gnu"));
        let side = X_MAPPED_MIN_SIDE;
        assert_eq!(side, 65_536);
        assert_eq!(std::mem::size_of::<DiagStruct>(), 4);
        assert!(2 * side * std::mem::size_of::<DiagStruct>() >= X_MAPPED_MIN);
        for (span, chunk_len, holds) in [
            (side - 1, side, false),
            (side, side - 1, false),
            (side, side, true),
            (side - 1, side - 1, false),
            (1 << 22, side, true),
            (side, 1 << 22, true),
            (4_000, 202_000, false),
            (220_000, 1_000, false),
            (220_000, 202_000, true),
            (0, 0, false),
            (usize::MAX, usize::MAX, true),
        ] {
            let at = format!("span {span} chunk {chunk_len}");
            assert_eq!(x_unit_mapped(true, span, chunk_len), glibc && holds, "{at}");
            assert!(!x_unit_mapped(false, span, chunk_len), "{at}: switch off");
        }
        assert_eq!(X_MAPPED_BYTES, (32 << 20) + 4096);
        assert!(x_mapped_len::<DiagStruct>(0) * std::mem::size_of::<DiagStruct>() > X_MAPPED_BYTES);
        assert!(x_mapped_len::<XUnitHsp>(0) * std::mem::size_of::<XUnitHsp>() > X_MAPPED_BYTES);
        assert_eq!(x_mapped_len::<DiagStruct>(1 << 24), 1 << 24);
        assert_eq!(x_mapped_len::<XUnitHsp>(1 << 21), 1 << 21);
    }

    // No NCBI counterpart: LOSAT_X_SEEDUNMAP's THP gate on the texts of
    // /sys/kernel/mm/transparent_hugepage/enabled and /proc/self/status: huge pages round the
    // mapped blocks up with "[always]", an unreadable or unrecognised file, unless the process has
    // them disabled (`THP_enabled: 0`); not with "[madvise]" or "[never]".
    #[test]
    fn seed_unmap_thp_gate_reads_the_two_files() {
        let status_1 = "Name:\tLOSAT\nUmask:\t0022\nTHP_enabled:\t1\nThreads:\t8\n";
        let status_0 = "Name:\tLOSAT\nUmask:\t0022\nTHP_enabled:\t0\nThreads:\t8\n";
        let status_without = "Name:\tLOSAT\nUmask:\t0022\nThreads:\t8\n";
        for (enabled, status, rounds_up) in [
            (Some("always [madvise] never\n"), Some(status_1), false),
            (Some("[always] madvise never\n"), Some(status_1), true),
            (Some("always madvise [never]\n"), Some(status_1), false),
            (None, Some(status_1), true),
            (Some("[always] madvise never\n"), Some(status_0), false),
            (None, Some(status_0), false),
            (Some("[always] madvise never\n"), Some(status_without), true),
            (
                Some("always [madvise] never\n"),
                Some(status_without),
                false,
            ),
            (Some("[always] madvise never\n"), None, true),
            (Some("always madvise [never]\n"), None, false),
            (None, None, true),
            (Some(""), None, true),
            (Some("always madvise never\n"), None, true),
            (Some("[madvise"), None, true),
            (
                Some("[always] madvise never\n"),
                Some("THP_enabled:\t01\n"),
                true,
            ),
        ] {
            assert_eq!(
                x_thp_rounds_up(enabled, status),
                rounds_up,
                "{enabled:?} {status:?}"
            );
        }
    }

    /// One LOSAT_X_SEEDUNMAP case: one query context of `W` blocks of 100 letters around `q_fill`
    /// letters of `D`, and two subject frames of `s_len` letters in blocks of 20 `W` and 380 `P`
    /// (one chunk each). Every query position has lookup entries (the context's span is about
    /// `q_fill + 200`), the subject hits only the `W` blocks, and each block saves about 220 init
    /// HSPs, so long sides stay cheap and a unit's copy passes `X_MAPPED_MIN`.
    fn unmap_case(q_fill: usize, s_len: usize) -> CompactCase {
        let query = format!("{0}{1}{0}", "W".repeat(100), "D".repeat(q_fill));
        let frames = vec![frame(1, &aa(&query))];
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (lookup, contexts) = build_ncbi_lookup(&[frames], 11, &params, true);
        let block = format!("{}{}", "W".repeat(20), "P".repeat(380));
        let subject: String = block.chars().cycle().take(s_len).collect();
        let shifted: String = block.chars().cycle().skip(200).take(s_len).collect();
        CompactCase {
            pp: XPairPar {
                groups: build_groups(&lookup, &contexts),
                stub: false,
                sent: x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice())),
                sliced: None,
                mask_rule: None,
                extents: Vec::new(),
            },
            cutoff_scores: vec![12; contexts.len()],
            x_dropoff_per_context: vec![16; contexts.len()],
            contexts,
            s_frames: vec![frame(2, &aa(&subject)), frame(-1, &aa(&shifted))],
        }
    }

    /// The units of `case`, run in order on a scratch with LOSAT_X_SEEDUNMAP and on one without:
    /// the same records, the same copies (with the mapped capacity where a large copy of a unit of
    /// the rule gets it, on Linux with glibc), the same lists of `run_unit` and the same merge.
    /// Returns the units where the rule holds and the copies given the mapped capacity.
    fn unmap_check(case: &CompactCase) -> (usize, usize) {
        let inp = case.inp();
        let (frame_chunks, units) = list_units(&case.pp, &inp, MAX_DBSEQ_LEN);
        assert!(!units.is_empty());
        let mut today = XScratch::new(false);
        let mut mapped = XScratch::new(true);
        let (mut held, mut copies) = (0usize, 0usize);
        let mut today_records = Vec::new();
        let mut mapped_records = Vec::new();
        for unit in &units {
            let at = format!("unit f {} k {} g {}", unit.f, unit.k, unit.g);
            let span = x_unit_span(&case.pp.groups[unit.g]);
            let rule = x_unit_mapped(true, span, unit.chunk.length);
            held += usize::from(rule);
            let exact = case.pp.run_unit_compact::<false>(unit, &inp, &mut today);
            let copy = case.pp.run_unit_compact::<false>(unit, &inp, &mut mapped);
            assert_eq!(copy, exact, "{at}");
            assert_eq!(
                exact.capacity(),
                exact.len(),
                "{at}: spare room without the switch"
            );
            let large = rule && copy.len() * std::mem::size_of::<XUnitHsp>() >= X_MAPPED_MIN;
            let capacity = if large {
                x_mapped_len::<XUnitHsp>(copy.len())
            } else {
                copy.len()
            };
            assert_eq!(copy.capacity(), capacity, "{at}");
            copies += usize::from(large);
            // Today's lists on the scratch with the mapped table, and the records expand to them.
            let list = case.pp.run_unit::<false>(unit, &inp, &mut today);
            assert_eq!(
                case.pp.run_unit::<false>(unit, &inp, &mut mapped),
                list,
                "{at}"
            );
            let mut expanded = Vec::new();
            x_expand_unit(&copy, unit.f, &inp, &mut expanded);
            assert_eq!(expanded, list, "{at}");
            // LOSAT_X_TBXKMERGE: the task's sort of the copy, with the scratch it was copied from as
            // the buffer (the mapped capacity where the rule holds), is `slice::sort_by` and
            // allocates nothing.
            let buffer = mapped.records.capacity();
            let mut in_task = copy.clone();
            in_task.x_sort_in(&mut mapped);
            let mut sorted = copy.clone();
            sorted.x_sort();
            assert_eq!(in_task, sorted, "{at}");
            assert_eq!(
                mapped.records.capacity(),
                buffer,
                "{at}: the sort allocated"
            );
            today_records.push(exact);
            mapped_records.push(copy);
        }
        if held > 0 {
            assert!(mapped.diag.capacity() >= x_mapped_len::<DiagStruct>(0));
            assert!(mapped.records.capacity() >= x_mapped_len::<XUnitHsp>(0));
        } else {
            assert!(mapped.diag.capacity() < x_mapped_len::<DiagStruct>(0));
            assert!(mapped.records.capacity() < x_mapped_len::<XUnitHsp>(0));
        }
        assert!(today.diag.capacity() < x_mapped_len::<DiagStruct>(0));
        assert!(today.records.capacity() < x_mapped_len::<XUnitHsp>(0));
        let today_hits = XPairPar::merge_subject_listed(
            XSubjectUnits {
                frame_chunks: frame_chunks.clone(),
                units: units.clone(),
                outputs: XUnitOutputs::Compact(today_records),
            },
            &inp,
        );
        let mapped_hits = XPairPar::merge_subject_listed(
            XSubjectUnits {
                frame_chunks,
                units,
                outputs: XUnitOutputs::Compact(mapped_records),
            },
            &inp,
        );
        shadow_compare(inp.s_idx, &mapped_hits, &today_hits);
        assert!(!today_hits.is_empty());
        (held, copies)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
    // ```c
    //         diag_table->hit_level_array = (DiagStruct *)
    //             calloc(diag_table->diag_array_length, sizeof(DiagStruct));
    // ```
    // LOSAT_X_SEEDUNMAP on units where the rule holds (a context of about 66,200 letters against
    // subject frames of 66,000: both sides long; every unit's copy passes 128 KiB), where the
    // query side is short (a context of about 1,200 letters; the copies pass 128 KiB too) and where
    // the chunk is short (2,000 letters): the same records, lists and merge, and mapped capacities
    // only where the rule holds.
    #[test]
    fn seed_unmap_units_with_the_mapped_capacities_are_todays() {
        let glibc = cfg!(all(target_os = "linux", target_env = "gnu"));
        for (q_fill, s_len, rule) in [
            (66_000, 66_000, true),
            (1_000, 66_000, false),
            (66_000, 2_000, false),
        ] {
            let at = format!("query fill {q_fill} subjects of {s_len}");
            let case = unmap_case(q_fill, s_len);
            let (held, copies) = unmap_check(&case);
            let expected = if glibc && rule { 2 } else { 0 };
            assert_eq!((held, copies), (expected, expected), "{at}");
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1043
    // ```c
    //     aux_struct->init_hitlist = BLAST_InitHitListNew();
    // ```
    // On search pools of 1, 2 and 8 threads (at one thread the serial scratch, which never maps;
    // at 2 and 8 one scratch per rayon job, with the switch as it acts in this process): the
    // records of units where the rule holds are those of a scratch without it.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    #[test]
    fn seed_unmap_units_on_search_pools_are_todays() {
        let case = unmap_case(66_000, 66_000);
        let inp = case.inp();
        let (_, units) = list_units(&case.pp, &inp, MAX_DBSEQ_LEN);
        let mut today = XScratch::new(false);
        let expected: Vec<Vec<XUnitHsp>> = units
            .iter()
            .map(|unit| case.pp.run_unit_compact::<false>(unit, &inp, &mut today))
            .collect();
        assert!(expected.iter().all(|records| !records.is_empty()));
        for threads in [1usize, 2, 8] {
            crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
                assert_eq!(pool.enabled(), threads > 1);
                let got = case.pp.run_units(
                    &units,
                    &inp,
                    pool,
                    XPairPar::run_unit_compact::<false>,
                    false,
                );
                assert_eq!(got, expected, "threads {threads}");
                Ok(())
            })
            .expect("search pool");
        }
    }

    // LOSAT_X_TBXPREFETCH's rule at its boundaries: off, never; on, only where one of the tables
    // has at least the switch's number of KiB (the built-in 1,536 and the smallest value, 2).
    #[test]
    fn tbx_prefetch_rule_at_its_boundaries() {
        let built_in = X_PREFETCH_MIN_KIB * 1024;
        assert!(!x_prefetch_units(None, [usize::MAX].into_iter()));
        assert!(!x_prefetch_units(
            Some(X_PREFETCH_MIN_KIB),
            std::iter::empty()
        ));
        assert!(!x_prefetch_units(
            Some(X_PREFETCH_MIN_KIB),
            [0, built_in - 1].into_iter()
        ));
        assert!(x_prefetch_units(
            Some(X_PREFETCH_MIN_KIB),
            [0, built_in, 1].into_iter()
        ));
        assert!(x_prefetch_units(Some(2), [2048].into_iter()));
        assert!(!x_prefetch_units(Some(2), [2047].into_iter()));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-516
    // ```c
    //     for (i = 0; i < hits; ++i) {
    //         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
    //         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
    // ```
    // The unit functions with the prefetch of LOSAT_X_TBXPREFETCH (`PF`) on search pools of 1, 2
    // and 8 threads, in both list forms: the random cases of the compact tests (one and several
    // chunks per frame) and the large case of the unmap tests (batches far longer than
    // X_PREFETCH_AHEAD) give the lists of the functions without it.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    #[test]
    fn tbx_prefetch_units_are_todays() {
        let mut rng = Rng(0x5eed_c0de_0000_0005);
        let mut cases = vec![(unmap_case(66_000, 66_000), MAX_DBSEQ_LEN)];
        for case_no in 0..12 {
            let alphabet = aa(COMPACT_ALPHABETS[case_no % COMPACT_ALPHABETS.len()]);
            cases.push(if case_no % 2 == 0 {
                (
                    CompactCase::random(&mut rng, &alphabet, 3, 400),
                    MAX_DBSEQ_LEN,
                )
            } else {
                (CompactCase::random(&mut rng, &alphabet, 200, 400), 160)
            });
        }
        for threads in [1usize, 2, 8] {
            crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
                let mut records = 0usize;
                for (case_no, (case, max_dbseq_len)) in cases.iter().enumerate() {
                    let inp = case.inp();
                    let (_, units) = list_units(&case.pp, &inp, *max_dbseq_len);
                    let at = format!("threads {threads} case {case_no}");
                    let compact = case.pp.run_units(
                        &units,
                        &inp,
                        pool,
                        XPairPar::run_unit_compact::<false>,
                        false,
                    );
                    records += compact.iter().map(Vec::len).sum::<usize>();
                    assert_eq!(
                        case.pp.run_units(
                            &units,
                            &inp,
                            pool,
                            XPairPar::run_unit_compact::<true>,
                            false
                        ),
                        compact,
                        "{at}"
                    );
                    assert_eq!(
                        case.pp
                            .run_units(&units, &inp, pool, XPairPar::run_unit::<true>, false),
                        case.pp
                            .run_units(&units, &inp, pool, XPairPar::run_unit::<false>, false),
                        "{at}"
                    );
                }
                assert!(records > 500, "threads {threads}: records {records}");
                Ok(())
            })
            .expect("search pool");
        }
    }
}
