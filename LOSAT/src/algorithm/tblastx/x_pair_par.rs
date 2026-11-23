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
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805,835-844
//! NCBI runs the two-hit word finder once per subject frame (and per chunk of a frame longer than
//! MAX_DBSEQ_LEN), with one scan stream that emits the hits of every query context, on one diagonal
//! table shared by the frames. This path runs the same two-hit logic (`x_hit`, the unit path's own
//! form of `word_finder::two_hit_one_hit`: its two-hit statements on the unit's unmasked cell,
//! without its trace, debug, diagnostics and HSP-saving statistics, which the gates send to the
//! reference loop; the extension, sink and probe are its own) on units: for each (frame, chunk),
//! one scan per query context with a lookup table that holds only that context's entries
//! (`XGroup`), on a fresh diagonal table without masking. The
//! unit lists are concatenated in context order and then go through the reference's own steps in
//! the reference's order: the stable init-HSP score sort (`Blast_InitHitListSortByScore`),
//! `BLAST_GetUngappedHSPList`, the chunk offset adjustment and merge, and the per-frame append and
//! sort (`frame_chunks::append_frame_hits`). Only the placement of work changes.
//!
//! # Why the units give the reference's `combined_ungapped_hits`
//!
//! Notation: QL = `query_length` (concatenated contexts with their shared sentinels), 2^k the
//! reference table (`diag_array_size`, next power of two >= QL + window), slack = 2^k - QL
//! (>= window), D the reference `diag_offset` of a chunk, len_c a context length, s_len the
//! chunk length. A hit is (q, s); its raw diagonal is d = q - s; the reference cell is
//! d & (2^k - 1), the unit cell is d itself (offset by the unit's lowest diagonal).
//!
//! 1. Cell-local state. The two-hit test of a hit reads and writes only its own cell (`last_hit`
//!    relative to D, and `flag`); cutoffs, x-drop and D are fixed during a chunk; the extension
//!    reads only the two sequences (aa_ungapped.c:516-606; the doc of
//!    `word_finder::two_hit_one_hit`). So a hit's decision depends only on the earlier hits of its
//!    cell. All comparisons are of `last_hit - D` with `s`, so a unit starting at D = window and
//!    the reference at any D agree whenever their cells hold the same `last_hit - D` and flag, and
//!    `s + D` fits the 31-bit `last_hit` in the reference as it does in the unit (`last_hit_fits`;
//!    it fails only for a `-window_size` near 2^29).
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
//!      `D` at or above `INT4_MAX` divided by 4, every cell is `-window` (`s_BlastDiagClear`, blast_extend.c:87-104) with
//!      D = window: `last_hit - D = -2 window`, too far. The value is read through the
//!      sign-extending `DiagStruct::last_hit` (round-1 bug: `-window` read as `2^31 - window`).
//!    - cross-sentinel diagonals: one raw diagonal can run from context c into c+1 (q grows with
//!      s, so all c hits come first). Words do not span the sentinel, so the first c+1 hit is at
//!      least `wordsize` + 1 letters after the last c hit (`diff >= wordsize`). With flag 0,
//!      either diff >= window, or `query_offset - diff` lies in context c, the context check fails
//!      (aa_ungapped.c:566-573) and `last_hit = s + D`. With flag 1, the c extension's right walk
//!      stopped at the NULLB after context c (NCBI's walk runs on the batch, and so does the units'
//!      through `XHitCtx::query_batch`, extension/query_batch.rs; the NULLB's BLAST_SCORE_MIN ends
//!      every walk of an X-drop of at most 32,768), so `s_last_off - (wordsize - 1) < s` of the
//!      first c+1 hit and the hit resets the cell. A walk of a larger X-drop that passes the NULLB
//!      without a new best (its stop is past it; `QueryBatch` notes the pass) writes a `last_hit`
//!      past it, which masks the hits under it of the later contexts in the reference's cell, and
//!      in the unit's cell where those contexts are in the unit's group. A unit of another group
//!      starts that cell fresh, so where a unit of a batch of several groups notes a pass,
//!      `run_subject_units` runs the units of every chunk of the subject again, a chunk's units in
//!      group order on one unmasked table of all the groups' diagonals (`x_units_in_order`): a
//!      cell then sees the hits of every context in ascending s, as the reference's cell does
//!      (items 1, 2), and starts fresh. The units' refusals (extension/query_batch.rs) wait for
//!      that choice, since a unit that missed a pass may refuse a hit the reference masks. With
//!      one group the unit's cell is the reference's (and a masked table's aliasing stays harmless:
//!      a walk stops at the batch's end, within the group's extent, below).
//!    - masked aliasing: raw diagonals d > d' share a reference cell when d - d' = m 2^k. Both
//!      exist only if s_len >= slack + 2. Then for any hit (q, s) on d and (q', s') on d',
//!      s' - s = m 2^k - (q - q') >= slack + 3: the cell sees all hits of d, then all of d', with
//!      a gap > slack >= window. At the first d' hit, flag 0 gives diff > slack >= window (too far).
//!      Flag 1 comes from the last right extension of d, started by a hit (q0, s0), which wrote
//!      `last_hit = s_last_off - (wordsize - 1) + D`. Its right walk ends at the end of the batch
//!      at the latest (`extend_hit_two_hit_fitting` walks the context's letters right of the word
//!      start, and `QueryBatch::walk_right` goes on at most to the end of the last context), so
//!      `s_last_off - s0 <= QL - q0`. A hit
//!      (q', s') on d' has q' >= 0, so `s' - s0 = q' - q0 + m 2^k >= 2^k - q0 >= QL + window - q0`,
//!      and `s' >= s_last_off + window`: the first d' hit is never masked; it resets the cell.
//!
//!    Hence no condition on the lengths: with LOSAT_X_TBXBAND (section below) the units take every
//!    batch the gates let through and `last_hit_fits` holds. Without the switch the guard of the
//!    earlier, weaker argument (`guard`: max context length <= slack, or longest subject chunk
//!    <= slack) sends the other batches to the reference path ("fallback: guard"). That argument
//!    bounded `s_last_off - s0` by `len_c` and `s' - s0` by `slack + 3` apart, and so needed
//!    that `len_c` is at most `slack`.
//! 4. Same list. Within a unit the saved InitHSPs are the reference's saved InitHSPs of that
//!    context, in the same order. The reference sorts the chunk list with the stable
//!    `score_compare_match` (score, s_start, length, absolute q_start; blast_extend.c:274-310,
//!    `sort_init_hsps_by_score_ncbi` in `word_finder::finish_chunk_scan`). Comparator-equal records
//!    share s_start and the absolute q_start, so they lie on one raw diagonal of one context, in
//!    one unit, in reference order. A stable sort of the concatenation c = 0..n therefore returns
//!    the reference's sorted list element for element, and every later step
//!    (`get_ungapped_hsp_list`, chunk merge, frame append and sort) is the reference's code on the
//!    same input.
//!
//! Contexts of a multi-query batch are taken in at most `X_MAX_GROUPS` contiguous groups (without
//! LOSAT_X_TBXFEWGROUPS, one context per group for one query). A group is one unit: inside it the diagonals of its contexts
//! share unit cells exactly as reference cells without aliasing, and items 1-4 apply to its
//! boundaries. The groups depend on the inputs and, with LOSAT_X_TBXFEWGROUPS, on the units the pool
//! can use (`x_groups_wanted`, `group_starts`), which changes no output: a group boundary only decides
//! where the concatenation above is cut; with one thread they run in order on the calling thread.
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
//! `LOSAT_X_TBXCTXTABLE` (in the strict set: on unless it is `0`): the units' two-hit test finds
//! the context of a hit that passes it (NCBI's `BSearchContextInfo`) in the `XContextTable` of
//! blastp's and blastx's fused loops (`common/x_context_table.rs`),
//! made once per batch from the frame bases (`x_ctx_table`): a coarse table of at most 16,384
//! entries gives the context at the start of the offset's block and a walk over the frame bases
//! moves on to the last context that starts at or before it, the index of the binary search
//! (`BlastAaLookupTable::get_context_idx`). The binary search took every context of the batch
//! (NCBI's narrowing by the shortest and longest context is not ported), and its comparison
//! depends on the hit's query offset: on 40 phage genomes of 37 to 46 kb all against all (six
//! contexts a batch) it was mispredicted about once per search, 27 % of the unit scan's
//! mispredictions. With `=0` the table covers no offset and every hit takes the binary search, the
//! path before the switch. Same index, so the same init HSPs.
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
//! two to seven, at most 8,192 at one (`x_mask_max_cells`), or, in a batch that `guard` rejects
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
//! d'. Flag 0: the cell holds `s_d` plus `D` of a hit on d and diff = s' - s_d > window, too far, the
//! state of a fresh cell after its first hit. Flag 1: the last right extension of d, from (q0, s0),
//! stopped at the end of q0's context at the latest, at most a + G, so s_last_off - s0 <= a + G -
//! q0, while s' - s0 = q' - q0 + m P >= a - q0 + P: s' >= s_last_off + window, the hit is not masked,
//! it resets the cell to the fresh state. (A walk past a NULLB stops later, but in a batch of one
//! group still at the batch's end, at most a + G; with several groups such a pass makes the
//! subject's units run again on one unmasked table, item 3.) So the first hit of d' leaves its cell as on an unmasked
//! table, and by induction every hit takes the branch it takes there. The scan in slices cannot use
//! such a table: a slice's pass over the whole chunk brings hits of d' (low q') before the later
//! slices' hits of d. The masked loops are their own instances (`scan_unit_masked`,
//! `scan_unit_batched_masked`); the unmasked ones are the code before the switch. The
//! `[X_STATS] TBXUNITMASK` line counts the units with a masked table.

use super::*;
use crate::algorithm::common::x_context_table::XContextTable;
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
// The offset-pair array of a unit is sized from its own lookup table, as NCBI sizes it from the
// table it scans with (lookup_wrap.c:255-266).
const OFFSET_ARRAY_SIZE: i32 = 4096;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-505
// No NCBI counterpart: how many units a batch takes (`group_starts`, scheduling only).
/// LOSAT_X_TBXFEWGROUPS: a strict-set switch (on unless it is `0`; `utils/xswitch.rs`).
fn few_groups() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXFEWGROUPS"))
}

/// No NCBI counterpart: LOSAT_X_TBXFEWGROUPS asks for this many units per thread of the pool: the
/// groups serve only parallelism, and each scans every subject frame and has its tables.
const X_UNITS_PER_THREAD: usize = 2;

/// No NCBI counterpart: the frames of a subject (`XUnit::f`), which are the units of one group.
const X_SUBJECT_FRAMES: usize = 6;

/// No NCBI counterpart: the groups a batch needs for `threads`: enough for `X_UNITS_PER_THREAD`
/// units per thread on the longest subject's share of the bases (its frames are all the units a
/// chromosome among plasmids can split into), so one group at one thread or with subjects enough.
fn x_groups_wanted(
    threads: usize,
    subjects: crate::algorithm::tblastx::subjects::Subjects<'_>,
) -> usize {
    if !few_groups() {
        return X_MAX_GROUPS;
    }
    let (longest, total) = subjects.longest_and_total();
    let (longest, total) = (longest as u128, u128::from(total));
    let shares = usize::try_from(total / longest.max(1))
        .unwrap_or(usize::MAX)
        .max(1);
    (threads * X_UNITS_PER_THREAD)
        .div_ceil(shares.saturating_mul(X_SUBJECT_FRAMES))
        .clamp(1, X_MAX_GROUPS)
}

/// No NCBI counterpart: LOSAT_X_TBXFEWGROUPS takes fewer groups up to this many letters of a batch
/// (a nucleotide query of about 65 kb). Measured (workstation, one thread, 2026-10-31): 60 kb x 1 Mb
/// x0.94; 40 x 40 phages of 43 kb x0.81; 50 kb x 4.6 Mb x0.87; 50 cDNAs x 1 Mb x0.77. Above it the
/// group's diagonal table outgrows the cache: 100-300 kb queries x1.03-1.05 in one or three groups,
/// x1.05-1.11 for 291 kb x 590 kb (BLOSUM45), so those batches keep their six groups.
const X_FEW_GROUPS_LETTERS: usize = 1 << 17;

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

/// No NCBI counterpart: the largest masked table, in cells (32 KiB), of LOSAT_X_TBXUNITMASK's size
/// rule at one thread (module doc).
const X_UNITMASK_ONE_THREAD_CELLS: usize = 1 << 13;

// No NCBI counterpart: the size limit of LOSAT_X_TBXUNITMASK's rule (module doc).
/// The largest masked table, in cells, that a unit takes for its size at `threads` threads, for
/// the switch's `kib`: at most `X_UNITMASK_ONE_THREAD_CELLS` at one thread,
/// `X_UNITMASK_FEW_THREADS_CELLS` below `X_UNITMASK_MANY_THREADS`, the switch's KiB from there.
fn x_mask_max_cells(kib: usize, threads: usize) -> usize {
    let cells = kib.saturating_mul(1024) / std::mem::size_of::<DiagStruct>();
    match threads {
        0 | 1 => cells.min(X_UNITMASK_ONE_THREAD_CELLS),
        t if t < X_UNITMASK_MANY_THREADS => cells.min(X_UNITMASK_FEW_THREADS_CELLS),
        _ => cells,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:510
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
    /// PROBE (exp/fx-tbxhit, line FX s4): `XUnmasked`'s lowest raw diagonal as an exact integer
    /// (q_min - (s_len - 1) lies inside i32: query offsets are NCBI's Int4), so that the cell of
    /// (q, s) is exactly q - s - lo_exact; 0 for `XMasked` (unused).
    #[inline(always)]
    fn lo_exact(self) -> isize {
        0
    }
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
    #[inline(always)]
    fn lo_exact(self) -> isize {
        self.lo as i32 as isize
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
fn x_group_extents(contexts: &[QueryContext], wanted: usize) -> Vec<usize> {
    x_tbx_qslice::x_group_spans(contexts, &group_starts(contexts, wanted))
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

// No NCBI counterpart: which LOSAT form finds a hit's context in the units' two-hit test (module
// doc); both give `BSearchContextInfo`'s index. Read once per process, as the other modes of this
// file.
/// LOSAT_X_TBXCTXTABLE: a strict-set switch (`utils/xswitch.rs`).
fn ctx_table_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXCTXTABLE"))
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:168-173
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
pub(super) struct XBatch<'a, M: TbxScore = Blosum62Score> {
    pub lookup: &'a BlastAaLookupTable,
    pub contexts: &'a [QueryContext],
    pub query_length: i32,
    pub diag_array_size: i32,
    pub subjects: crate::algorithm::tblastx::subjects::Subjects<'a>,
    /// A trace, debug or diagnostics switch is on, the subject has lower-case letters with
    /// `-lcase_masking`, or the word finder is the one-hit one (`-window_size 0`): they exist only
    /// in the reference loop.
    pub reference_loop: bool,
    /// LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS or the parallel-chunk path is on.
    pub chunk_switches: bool,
    /// PAIRLUT: the groups built at the lookup build site (`lut_build`), or None.
    pub lut: Option<XLutBatch>,
    /// The threads of the batch's search pool (`SearchPool::threads`), as the build site reads them.
    pub threads: usize,
    /// The lookup threshold, as the build site reads it (`XLutEarly::threshold`).
    pub threshold: i32,
    /// The search's matrix (`-matrix`; no data for BLOSUM62's instance).
    pub sc: M,
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
    pub reference_loop: bool,
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
// The decision of `prepare` (the gates, then the guard against the reference table size above,
// which LOSAT_X_TBXBAND drops, and the second guard, `last_hit_fits`) as one function, so that the
// lookup build site (`lut_build`, before the table exists) and `prepare` (after it) decide with the
// same code.
/// Whether the unit path runs for a batch.
pub(super) fn decide(g: &XGates<'_>) -> XDecision {
    if g.reference_loop || g.chunk_switches {
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

// No NCBI counterpart: `decide` for a batch whose groups are never cut (SC10b PR-3: words of 4
// letters, `x_tbx_qslice::x_slices_fit`). With LOSAT_X_TBXBAND and LOSAT_X_TBXUNITMASK the memory
// rule of such a batch holds whatever LOSAT_X_TBXQSLICE's rule says (`x_band_masked` without
// `rule_cuts`), as `prepare` keeps its masked tables: the first guard holds. `decide` reads the
// contexts for that guard only (the longest context, `x_guard_lengths`, and `x_band_masked`); with
// no context the guard holds (a slack of `diag_array_size - query_length > 0` letters), and the
// gates before it and the second guard and the cost gate after it read the same values. So this is
// `decide` on the gates without their contexts.
/// `decide` for a batch whose groups are never cut.
#[inline(never)]
pub(super) fn decide_uncut(g: &XGates<'_>) -> XDecision {
    if g.band && g.unit_mask {
        return decide(&XGates {
            contexts: &[],
            query_length: g.query_length,
            diag_array_size: g.diag_array_size,
            subjects: g.subjects,
            reference_loop: g.reference_loop,
            chunk_switches: g.chunk_switches,
            threads: g.threads,
            cost_gate: g.cost_gate,
            n1_gate: g.n1_gate,
            threshold: g.threshold,
            band: g.band,
            unit_mask: g.unit_mask,
        });
    }
    decide(g)
}

// No NCBI counterpart: the choice between `decide` and `decide_uncut` by the search's word length
// (`x_tbx_qslice::x_slices_fit`; constant for BLOSUM62's instance, which calls `decide`).
/// The decision for a batch of `word`-letter words.
#[inline(always)]
pub(super) fn decide_for(word: usize, g: &XGates<'_>) -> XDecision {
    if x_tbx_qslice::x_slices_fit(word) {
        decide(g)
    } else {
        decide_uncut(g)
    }
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
// `rule_cuts` reads only the contexts, the thread count, the subjects' lengths (through the
// groups, `x_groups_wanted`) and process-wide switches (not the
// diagnostics ones), so the lookup build site and `prepare` agree, with and without diagnostics.
/// LOSAT_X_TBXBAND with LOSAT_X_TBXUNITMASK: the switch is on and LOSAT_X_TBXQSLICE's rule cuts
/// no group of the batch.
fn x_band_masked(g: &XGates<'_>) -> bool {
    g.unit_mask
        && !x_tbx_qslice::rule_cuts(
            g.contexts,
            &group_starts(g.contexts, x_groups_wanted(g.threads, g.subjects)),
            g.threads,
        )
}

// No NCBI counterpart: the lengths the guards compare, from the batch's contexts and subjects (a
// subject's frame length, at most one chunk).
/// The longest context and the longest subject chunk of a batch, in letters.
fn x_guard_lengths(
    contexts: &[QueryContext],
    subjects: crate::algorithm::tblastx::subjects::Subjects<'_>,
) -> (usize, usize) {
    let max_context_len = contexts.iter().map(|c| c.aa_len).max().unwrap_or(0);
    let longest_subject_chunk =
        (subjects.longest_and_total().0 / 3).min(tblastx_max_dbseq_len_for_run());
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
/// units split it over the pool. With LOSAT_X_TBXFEWGROUPS several threads take the one-thread line
/// too: one group for spread subjects, as at one thread. Without it they keep the reference loop
/// below `X_GATE_DIAG_CELLS`: the six groups used less CPU there but lost wall time on many short
/// subjects (x1.02-1.06 at eight threads, session 9 of line A).
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
    // LOSAT_X_TBXFEWGROUPS: the one-thread line, as the batch takes the groups the pool needs
    // (one with subjects enough) and no more scans of each frame than at one thread.
    if few_groups() && g.query_length >= X_GATE_N1_LETTERS {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXFEWGROUPS", true);
        return true;
    }
    let (longest, total) = g.subjects.longest_and_total();
    let (longest, total) = (longest as u128, u128::from(total));
    2 * longest * g.threads as u128 > 3 * total
}

/// The inputs of `lut_build`: what the lookup build site of the batch has
/// (`batch_setup::batch_lookup`).
pub(super) struct XLutEarly<'a, M: TbxScore = Blosum62Score> {
    pub query_frames: &'a [Vec<QueryFrame>],
    /// The lookup threshold (`lookup->threshold`).
    pub threshold: i32,
    pub window: i32,
    pub subjects: crate::algorithm::tblastx::subjects::Subjects<'a>,
    /// As `XBatch::reference_loop`.
    pub reference_loop: bool,
    /// `requested_parallel` and `use_serial_scan_chunks` of the batch (`batch_setup::batch_flags`).
    pub requested_parallel: bool,
    pub use_serial_scan_chunks: bool,
    pub pool: &'a crate::utils::threading::SearchPool<'a>,
    /// The search's matrix (`-matrix`; no data for BLOSUM62's instance): the groups' tables hold
    /// its neighbour words.
    pub sc: M,
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
    // The reference code reads `longest_chain` and the header of the batch table; the stub carries
    // the full table's values (`XPairLutQuery::x_into_stub`).
    /// The table and contexts that take the place of `build_ncbi_lookup`'s for this batch: the
    /// stub.
    pub(super) fn take_stub(&mut self) -> Option<(BlastAaLookupTable, Vec<QueryContext>)> {
        self.stub.take()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-250
// `chunk_switches` of the batch before the lookup build: the expression in `search_query_batch`
// (`use_parallel_chunks`, which reads the thread request, two switches and the subject lengths
// only) or'ed with `use_serial_scan_chunks`, as `seed_tables::scan_tables` passes it to `prepare`.
// `prepare` decides again with `search_query_batch`'s own value and stops the run if the two
// decisions differ for a batch with direct groups.
fn early_chunk_switches(
    requested_parallel: bool,
    use_serial_scan_chunks: bool,
    subjects: crate::algorithm::tblastx::subjects::Subjects<'_>,
) -> bool {
    let use_parallel_chunks = requested_parallel
        && !use_serial_scan_chunks
        && (std::env::var_os("LOSAT_TBLASTX_PARALLEL_CHUNKS").is_some()
            || (cfg!(all(target_arch = "wasm32", feature = "wasm-threads"))
                && tblastx_estimated_subject_chunk_work_items(
                    subjects.longest_and_total().0,
                    tblastx_max_dbseq_len_for_run(),
                ) > 1));
    use_serial_scan_chunks || use_parallel_chunks
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-97
// Dispatch point 0 (per batch, at the lookup build): when the unit path will run for the batch,
// the tables its units scan with are built here, one per context group, from the reference query
// preparation; the full table NCBI builds above is then not built. None leaves the batch to the
// reference builder and `build_groups`: the direct builder's preconditions
// (`XPairLutQuery::x_prepare`), a gate or the guard (`decide`), or a lookup location outside one
// group.
/// The direct groups of the batch, or None (the reference builder runs).
pub(super) fn lut_build<M: TbxScore>(e: XLutEarly<'_, M>) -> Option<XLutBatch> {
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
// The early decision: the inputs of `decide` from the prepared query, by the expressions
// `seed_tables::scan_tables` uses after the lookup build (`query_length`,
// `check_diag_table_window`, `diag_array_size`), then the groups and the stub.
fn lut_direct<M: TbxScore>(e: &XLutEarly<'_, M>) -> Option<XLutBatch> {
    // The reference query preparation (contexts, concatenated query, lookup locations), with
    // another matrix's scores or another word length that of `build_ncbi_lookup_for_profile`
    // (SC7b; `TbxScore::word`, 3 in BLOSUM62's instance). A word length neither builder covers
    // gives None here (`x_tables_cover`: TBLASTX takes words of 2 to 4 letters, all covered).
    let q = if M::BLOSUM62 {
        XPairLutQuery::x_prepare(e.query_frames, e.threshold, true)?
    } else {
        XPairLutQuery::x_prepare_for(
            e.query_frames,
            e.threshold,
            true,
            e.sc.matrix(),
            e.sc.word(),
        )?
    };
    let contexts = q.contexts();
    let query_length: i32 = contexts
        .last()
        .map(|c| c.frame_base + c.aa_len as i32)
        .unwrap_or(0);
    // `seed_tables::scan_tables` reports this error after the reference build.
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
        reference_loop: e.reference_loop,
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
    if decide_for(e.sc.word(), &gates) != XDecision::Units {
        return None;
    }
    // EXPERIMENT (LOSAT_X_TBXQSLICE): the tables of the slices of the groups the switch cuts,
    // built directly as the groups are (x_tbx_qslice.rs).
    let threads = e.pool.threads();
    let wanted = x_groups_wanted(threads, e.subjects);
    let starts = group_starts(contexts, wanted);
    // SC10b PR-3: words whose backbone index the 16-bit direct builder does not cover (4 letters)
    // take the joint builder of 32-bit indices on the pool, and their groups are not cut
    // (`lut_wide`; `x_tbx_qslice::x_slices_fit`). Constant for BLOSUM62's instance (3 letters).
    if !crate::algorithm::tblastx::lookup::x_pair_lut::x_direct_covers(e.sc.word()) {
        return lut_wide(q, &starts, e.pool);
    }
    let (groups, sliced) = match x_tbx_qslice::partition(contexts, &starts, threads) {
        None => (lut_direct_groups(&q, wanted, e.pool, lut_on_pool())?, None),
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
// One table per context group (`group_starts`, the groups of `build_groups`), each from the
// lookup locations of its contexts (`XPairLutQuery::x_build_group`). Scheduling only: the groups
// are independent. They are built in order on the calling thread; with `on_pool` (the tuning knob
// LOSAT_X_PAIRLUT_POOL, `lut_on_pool`) and a search pool they are built on the pool (an indexed
// map, results in group order; each build takes a memo of neighbour lists, `x_with_tbx_memo`). The
// calling thread is the default because the pool gained at most 17 ms per batch (AP027131, 6 x
// 220,702 aa) and cost memory: a memo for each group built at once (kept for the next builds), and
// tables allocated in the pool threads' arenas and freed by the search thread (+17 MiB on a
// 50-batch run at 8 threads).
/// The groups `build_groups` would filter out of the full table of `q`, built directly. None
/// when a lookup location does not lie inside one group.
fn lut_direct_groups(
    q: &XPairLutQuery,
    wanted: usize,
    pool: &crate::utils::threading::SearchPool<'_>,
    on_pool: bool,
) -> Option<Vec<XGroup>> {
    let starts = group_starts(q.contexts(), wanted);
    let locations = q.x_group_locations(&starts)?;
    Some(lut_direct_tables(q, &locations, pool, on_pool))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:454-456,461-465 (cited at
// `lut_direct_groups`). SC10b PR-3: the groups of `lut_direct_groups` for words whose backbone index
// does not fit 16 bits, built together by the joint builder of 32-bit indices on the search pool
// (`XPairLutQuery::x_build_tables_wide`, the tables `x_build_group` would build; lookup/
// x_pair_lut_wide.rs). Scheduling only: the pool's threads build the neighbour lists and fill the
// tables; the tables are the same at any thread count.
/// The groups of `starts` (`group_starts`) built by the joint builder. None when a lookup location
/// does not lie inside one group.
fn lut_wide_groups(
    q: &XPairLutQuery,
    starts: &[usize],
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Option<Vec<XGroup>> {
    let locations = q.x_group_locations(starts)?;
    let mapped = locations.len() <= X_MAX_GROUPS;
    Some(
        q.x_build_tables_wide(&locations, pool, mapped)
            .into_iter()
            .map(XGroup::from_direct)
            .collect(),
    )
}

// No NCBI counterpart beyond `lut_wide_groups`' and `x_into_stub`'s: the end of `lut_direct` for
// words of 4 letters, whose groups are never cut (`x_tbx_qslice::x_slices_fit`).
/// The batch of `lut_direct` for words the 16-bit direct builder does not cover: the groups of
/// `starts` by the joint builder and the stub. None when a lookup location does not lie inside one
/// group.
#[inline(never)]
fn lut_wide(
    q: XPairLutQuery,
    starts: &[usize],
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Option<XLutBatch> {
    let groups = lut_wide_groups(&q, starts, pool)?;
    let stub = q.x_into_stub(groups.iter().map(|g| &g.lookup));
    Some(XLutBatch {
        groups,
        sliced: None,
        stub: Some(stub),
    })
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
    let to_group = XGroup::from_direct;
    // A cut group's slices are many tables (x_tbx_qslice.rs), each a region of address space.
    let mapped = locations.len() <= X_MAX_GROUPS;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if on_pool && pool.enabled() {
        return pool.install(|| {
            locations
                .par_iter()
                .map(|group_locations| to_group(q.x_build_group(group_locations, mapped)))
                .collect()
        });
    }
    let _ = (pool, on_pool);
    locations
        .iter()
        .map(|group_locations| to_group(q.x_build_group(group_locations, mapped)))
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

impl XSubjectUnits {
    // EXPERIMENT (LOSAT_X_TBXBOARD): the records the units returned (a bound on the HSPs the
    // subject's merge keeps; tblastx/x_batch_pipe.rs `XGate`).
    pub(super) fn records(&self) -> usize {
        match &self.outputs {
            XUnitOutputs::Full(lists) => lists.iter().map(Vec::len).sum(),
            XUnitOutputs::Compact(lists) => lists.iter().map(Vec::len).sum(),
        }
    }
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
    /// The search's failure: an extension that LOSAT refuses fails the search at this subject
    /// (extension/query_batch.rs).
    pub failure: &'a crate::algorithm::tblastx::subjects::XReadFailure,
    /// Where the units note a right walk that went on past a NULLB (`run_subject_units` gives
    /// its units one; None elsewhere).
    pub passed: Option<&'a std::sync::atomic::AtomicBool>,
    /// `-subject_besthit`: each frame's list is filtered after each chunk's merge
    /// (`frame_subject_besthit`).
    pub subject_besthit: bool,
}

impl<'a> XSubject<'a> {
    /// The query batch as this subject's extensions read it (extension/query_batch.rs).
    #[inline(always)]
    fn query_batch(&self) -> QueryBatch<'a> {
        QueryBatch::noting(self.contexts, self.failure, self.passed)
    }
}

/// One unit-context group: its contexts' lookup entries and their query-offset bounds.
struct XGroup {
    /// The group's table. With `rank` not empty its backbone holds only the cells with entries
    /// (`x_scan_bf::x_compact`): it is read through `scan` and `hits`, never by index.
    lookup: BlastAaLookupTable,
    q_min: u32,
    q_max: u32,
    /// EXPERIMENT (LOSAT_X_TBXSCANBF): `x_scan_bf_for(&lookup)`, decided once per table.
    scan_bf: bool,
    /// The presence-vector rank of a compacted table, empty for a table of the full backbone.
    rank: Vec<u32>,
    /// EXPERIMENT (LOSAT_X_TBXSCANBYTE): the presence bytes of a table of the full backbone that
    /// NCBI's loop scanned and the scan on presence bytes takes (`x_scan_bf::x_tbx_scans_byte`),
    /// `None` for the others. `XPairPar::prepare` makes them (`x_make_presence`) when every table of
    /// the batch is built; they are freed with the table.
    ///
    /// Bytes made where each table was built (#1228: 32 KiB a table) sat between the group tables,
    /// which a batch builds one after the other on one thread, in that thread's heap. On Sakai x
    /// MG1655 `-subject` at 16 threads (84 slice tables of LOSAT_X_TBXQSLICE, 2.6 MiB of bytes)
    /// glibc then kept the freed tables in the main heap after the seed release (`[heap]` 624 MiB
    /// to the end of the run, 58 MiB with the switch off), and max RssAnon was 2,225 MiB against
    /// 1,984 MiB with the switch off. The scan reads this field as before: a lazy cell read in
    /// `scan` (made by its first call) gave the same memory as the switch off but changed the code
    /// of the units' scan loops (`scan_unit_bf`: two stack reloads per hit, 1.6 % more
    /// instructions on ViPTreeGen's Begomovirus.40).
    presence: Option<crate::algorithm::tblastx::x_scan_byte::XPresence>,
}

impl XGroup {
    // No NCBI counterpart: the layout of a group table (`x_scan_bf::x_compact`); its cells,
    // entries and presence bits are those of `lookup`.
    /// The group of `lookup`. A table the branch-free scan takes (`x_scan_bf_for`: fewer than
    /// half of its cells have entries) keeps only the cells with entries; the others get their
    /// presence bytes when LOSAT_X_TBXSCANBYTE scans them (`x_make_presence`).
    fn new(mut lookup: BlastAaLookupTable, q_min: u32, q_max: u32) -> Self {
        let scan_bf = crate::algorithm::tblastx::x_scan_bf::x_scan_bf_for(&lookup);
        let rank = if scan_bf {
            crate::algorithm::tblastx::x_scan_bf::x_compact(&mut lookup)
        } else {
            Vec::new()
        };
        if scan_bf {
            crate::algorithm::tblastx::x_scan_bf::x_reach_scan_byte_not_taken();
        }
        Self {
            lookup,
            q_min,
            q_max,
            scan_bf,
            rank,
            presence: None,
        }
    }

    // No NCBI counterpart: as `new`, for a table `x_build_group` laid out (compacted there where
    // `new` would compact it).
    /// The group of a direct table (`XPairLutQuery::x_build_group`).
    fn from_direct(g: XPairLutGroup) -> Self {
        let scan_bf = if g.rank.is_empty() {
            crate::algorithm::tblastx::x_scan_bf::x_scan_bf_for(&g.table)
        } else {
            crate::algorithm::tblastx::x_scan_bf::x_reach_scan_bf(true);
            true
        };
        if scan_bf {
            crate::algorithm::tblastx::x_scan_bf::x_reach_scan_byte_not_taken();
        }
        Self {
            lookup: g.table,
            q_min: g.q_min,
            q_max: g.q_max,
            scan_bf,
            rank: g.rank,
            presence: None,
        }
    }

    // No NCBI counterpart: the presence vector in another layout, made once per table.
    /// EXPERIMENT (LOSAT_X_TBXSCANBYTE): makes the presence bytes of this table when the scan on
    /// presence bytes takes it (a table of the full backbone, `x_tbx_scans_byte`); called once per
    /// table by `XPairPar::prepare` after the batch's tables are built.
    fn x_make_presence(&mut self) {
        if !self.scan_bf && crate::algorithm::tblastx::x_scan_bf::x_tbx_scans_byte(&self.lookup) {
            self.presence = Some(crate::algorithm::tblastx::x_scan_bf::x_tbx_presence(
                &self.lookup,
            ));
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-511
    /// One scan call of this group's table: `x_scan_subject_byte` on a table with presence bytes
    /// (LOSAT_X_TBXSCANBYTE), `x_scan_subject` (LOSAT_X_TBXSCANBF's choice), or
    /// `x_scan_subject_compact` on a compacted table. The offset array holds
    /// `x_offset_pairs_len(array_size)` entries.
    #[inline(always)]
    fn scan(
        &self,
        subject: &[u8],
        seq_ranges: &[(i32, i32)],
        offset_pairs: &mut [OffsetPair],
        array_size: i32,
        scan_range: &mut [i32; 3],
    ) -> i32 {
        if let Some(presence) = &self.presence {
            crate::algorithm::tblastx::x_scan_bf::x_scan_subject_byte(
                &self.lookup,
                presence,
                subject,
                seq_ranges,
                offset_pairs,
                array_size,
                scan_range,
            )
        } else if self.rank.is_empty() {
            crate::algorithm::tblastx::x_scan_bf::x_scan_subject(
                self.scan_bf,
                &self.lookup,
                subject,
                seq_ranges,
                offset_pairs,
                array_size,
                scan_range,
            )
        } else {
            crate::algorithm::tblastx::x_scan_bf::x_scan_subject_compact(
                &self.lookup,
                &self.rank,
                subject,
                seq_ranges,
                offset_pairs,
                array_size,
                scan_range,
            )
        }
    }

    /// The entries of cell `index` (none for a cell without its presence bit).
    #[cfg(test)]
    fn hits(&self, index: usize) -> &[i32] {
        if self.rank.is_empty() {
            return self.lookup.get_hits(index);
        }
        let word = index >> crate::algorithm::tblastx::lookup::PV_ARRAY_BTS;
        let bit = index & crate::algorithm::tblastx::lookup::PV_ARRAY_MASK;
        if (self.lookup.pv[word] >> bit) & 1 == 0 {
            return &[];
        }
        let below = self.lookup.pv[word] & ((1u32 << bit) - 1);
        let cell = &self.lookup.backbone[self.rank[word] as usize + below.count_ones() as usize];
        let num = cell.num_used as usize;
        if num <= AA_HITS_PER_CELL {
            &cell.entries[..num]
        } else {
            let cursor = cell.cursor();
            &self.lookup.overflow[cursor..cursor + num]
        }
    }
}

/// One unit: subject frame `f`, its chunk number `k`, context group `g` (`XPairPar::x_unit`).
#[derive(Clone, Copy)]
pub(super) struct XUnit {
    f: usize,
    k: usize,
    chunk: SubjectChunk,
    g: usize,
    /// EXPERIMENT (LOSAT_X_TBXUNGAPSENT): whether the frame `f` may take the walks the sentinels
    /// stop (`x_tbx_sent_frame`: the batch's contexts allow them and every byte of the frame's
    /// buffer is a key's). That check reads the whole buffer, so `x_unit` finds it once for each
    /// chunk where the units are listed (`run_subject_units`), not in every unit's
    /// `x_unit_sent`.
    sent_frame: bool,
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:336-341
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
/// kernel gives them only on request (`[madvise]`) or never (`[never]`); true for `[always]`, an
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
struct XHitCtx<'a, M: TbxScore = Blosum62Score> {
    diag_offset: i32,
    window: i32,
    wordsize: i32,
    lookup: &'a BlastAaLookupTable,
    /// LOSAT_X_TBXCTXTABLE: the batch's `x_ctx_table`, the context of a query offset without the
    /// binary search (`context_of`).
    ctx_of: &'a XContextTable,
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
    /// The search's matrix (`-matrix`; no data for BLOSUM62's instance): another matrix's
    /// extensions take its table walks (`extend_hit_two_hit_with`).
    sc: M,
    /// The batch's contexts as the extensions read them past a context's ends, and the search's
    /// failure (extension/query_batch.rs).
    query_batch: QueryBatch<'a>,
}

impl<M: TbxScore> XHitCtx<'_, M> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243
    /// `BSearchContextInfo` of a query offset (`x_context_of`).
    #[inline(always)]
    fn context_of(&self, query_offset: u32) -> usize {
        x_context_of(self.ctx_of, self.lookup, query_offset)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243 (as `context_of`)
/// LOSAT_X_TBXCTXTABLE: the context of `query_offset` from `table` (the binary search's index,
/// common/x_context_table.rs), or from the binary search over `lookup`'s frame bases for an offset
/// outside it (every offset with `=0`, whose table covers none).
#[inline(always)]
fn x_context_of(table: &XContextTable, lookup: &BlastAaLookupTable, query_offset: u32) -> usize {
    match table.get(query_offset as i32) {
        Some(idx) => idx,
        None => lookup.get_context_idx(query_offset as i32),
    }
}

// No NCBI counterpart: the table behind `XHitCtx::context_of` (the one blastp's and blastx's fused
// loops take), built once per batch from the frame bases `get_context_idx` searches.
/// LOSAT_X_TBXCTXTABLE: the `XContextTable` of a batch's contexts (`lookup.frame_bases`,
/// `num_contexts` of them) over its `query_length` offsets with the switch `on`, over none without
/// it (every offset then takes the binary search).
fn x_ctx_table(lookup: &BlastAaLookupTable, on: bool) -> XContextTable {
    let bases = &lookup.frame_bases[..lookup.num_contexts.min(lookup.frame_bases.len())];
    let span = if on {
        lookup.query_length.max(0) as usize
    } else {
        0
    };
    XContextTable::new(bases, span)
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:360-372
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2371-2375
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
//    unit's `f`, the same `XSubject` (dispatch point 3, `subject_prep::run_unit_path`, passes one
//    `x_subject` to `run_subject_units` and `merge_subject`) and `inp.contexts[ctx_idx]`.
//    `ctx_idx` is below the number of contexts; `u32::try_from` stops the run rather than
//    truncate it.
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
// Where `x_hit` saves an init HSP (`BlastSaveInitHsp`). `x_hit` and the scan of a unit
// (`scan_unit`) are generic over the sink, so each is compiled once per sink and no test of the
// switch runs per hit; the push statement of `x_hit` is today's.
/// The list a unit saves its init HSPs into.
trait XInitSink {
    /// Saves one init HSP as `x_hit` built it.
    fn push(&mut self, hsp: InitHSP);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:351-355
// Switch off: today's list, today's push.
impl XInitSink for Vec<InitHSP> {
    #[inline(always)]
    fn push(&mut self, hsp: InitHSP) {
        Vec::push(self, hsp);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:367-372
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
pub(super) struct XPairPar<M: TbxScore = Blosum62Score> {
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
    /// Another matrix's needs its byte table (`TbxScore::x_sent_walks`); only for words of 2 and 3
    /// letters, whose walks the units have (`scan_unit_or_bf`, `TbxScore::sent_two_hit_of`).
    sent: bool,
    /// The search's matrix (`-matrix`; no data for BLOSUM62's instance).
    sc: M,
    /// LOSAT_X_TBXCTXTABLE: the context of a query offset of the batch (`XContextTable`, as blastp's
    /// and blastx's, over the batch's frame bases, `x_ctx_table`), which the units' two-hit test reads for every
    /// hit that passes it.
    ctx_of: XContextTable,
}

impl<M: TbxScore> XPairPar<M> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:337-360
    // Dispatch point (per batch, after the lookup build): checks the gates and the guard (`decide`)
    // and builds the sub-lookups, or takes the ones PAIRLUT built at the lookup build site
    // (`lut_build`). None means the reference loop runs for every subject of the batch.
    /// The unit path for this batch, or None (a gate, the guard or the cost rule sends it to the
    /// reference loop).
    pub(super) fn prepare(b: XBatch<'_, M>) -> Option<Self> {
        let wanted = x_groups_wanted(b.threads, b.subjects);
        let decision = decide_for(
            b.sc.word(),
            &XGates {
                contexts: b.contexts,
                query_length: b.query_length,
                diag_array_size: b.diag_array_size,
                subjects: b.subjects,
                reference_loop: b.reference_loop,
                chunk_switches: b.chunk_switches,
                threads: b.threads,
                cost_gate: cost_gate(),
                n1_gate: n1_gate(),
                threshold: b.threshold,
                band: tbx_band(),
                unit_mask: unit_mask_kib().is_some(),
            },
        );
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
        if units {
            if let Some(few) = x_takes_few_groups(b.contexts, wanted) {
                crate::utils::xreach::x_reach!("LOSAT_X_TBXFEWGROUPS", few);
            }
        }
        if seed_compact() {
            crate::utils::xreach::x_reach!("LOSAT_X_SEEDCOMPACT", units);
        }
        if x_kmerge::on() {
            crate::utils::xreach::x_reach!("LOSAT_X_TBXKMERGE", units);
        }
        if ctx_table_on() {
            crate::utils::xreach::x_reach!("LOSAT_X_TBXCTXTABLE", units);
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
        let (mut groups, mut sliced, stub) = match b.lut {
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
                let starts = group_starts(b.contexts, wanted);
                // As at the build site (`lut_direct`): no cut where a slice's table does not fit
                // beside its window (`x_slices_fit`; constant for BLOSUM62's instance).
                let sliced = if x_tbx_qslice::x_slices_fit(b.sc.word()) {
                    x_tbx_qslice::partition(b.contexts, &starts, b.threads)
                        .map(|parts| x_tbx_qslice::build_sliced(b.lookup, &parts))
                } else {
                    None
                };
                let groups = if sliced.is_none() {
                    build_groups(b.lookup, b.contexts, wanted)
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
        // A batch whose groups are never cut (`x_slices_fit`; constant for BLOSUM62's instance)
        // keeps the rule's masked tables whatever the slicing rule says.
        let mask_rule = if x_tbx_qslice::x_slices_fit(b.sc.word()) {
            unit_mask_kib()
                .filter(|_| {
                    let ungrouped = sliced.is_none()
                        && !x_tbx_qslice::rule_cuts(
                            b.contexts,
                            &group_starts(b.contexts, wanted),
                            b.threads,
                        );
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
                })
        } else {
            unit_mask_kib().map(|kib| XMaskRule {
                max_cells: x_mask_max_cells(kib, b.threads),
                ratio: X_UNITMASK_RATIO,
                band_cells: band.then_some(b.diag_array_size),
            })
        };
        let extents = if mask_rule.is_some() {
            x_group_extents(b.contexts, wanted)
        } else {
            Vec::new()
        };
        debug_assert!(extents.is_empty() || extents.len() == groups.len());
        // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): checked once for the batch's contexts. The units'
        // extensions take the walks of the search's word length, 2 or 3 letters
        // (`scan_unit_or_bf`; `TbxScore::word`, 3 for BLOSUM62's instance).
        let sent = matches!(b.sc.word(), 2 | 3)
            && b.sc.x_sent_walks()
            && x_tbx_sent_contexts(b.contexts.iter().map(|c| c.aa_seq.as_slice()));
        let ctx_of = x_ctx_table(b.lookup, ctx_table_on());
        // EXPERIMENT (LOSAT_X_TBXSCANBYTE): the presence bytes of the tables the scan on presence
        // bytes takes, made now that every table of the batch is built, not between the tables
        // (`XGroup::presence`).
        for group in groups.iter_mut().chain(
            sliced
                .iter_mut()
                .flatten()
                .flat_map(x_tbx_qslice::XSliced::slices_mut),
        ) {
            group.x_make_presence();
        }
        Some(Self {
            groups,
            stub,
            sliced,
            mask_rule,
            extents,
            sent,
            sc: b.sc,
            ctx_of,
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
    // The seed stage of one subject, first half: the chunks of every frame and the units of every
    // (frame, chunk), run. This is the only part of the seed stage that reads the group tables
    // (`self.groups`: `longest_chain` while the units are listed, the tables in `scan_unit`).
    // `merge_subject` is the second half; the two were one function (`run_subject`) and are
    // called one after the other at dispatch point 2 (`subject_prep::run_unit_path`), where the
    // seed release ends its hold on the group tables between them (x_seed_release.rs).
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
                // The reference finds no hits in a chunk shorter than a word
                // (`word_finder::begin_chunk_scan`).
                if chunk.length < inp.wordsize as usize {
                    continue;
                }
                let first = self.x_unit(inp, f, k, *chunk, 0);
                for g in 0..self.x_group_count() {
                    if self.x_group_has_entries(g) {
                        units.push(XUnit { g, ..first });
                    }
                }
            }
            frame_chunks.push(chunks);
        }

        let t_units = timing().then(Instant::now);
        // The units refuse into a failure of their own and note a right walk past a NULLB
        // (`x_units_in_order`).
        let failure = crate::algorithm::tblastx::subjects::XReadFailure::new();
        let passed = std::sync::atomic::AtomicBool::new(false);
        let units_inp = XSubject {
            failure: &failure,
            passed: Some(&passed),
            ..*inp
        };
        let inp_all = inp;
        let inp = &units_inp;
        let mut outputs = if let Some(sliced) = &self.sliced {
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
        // A walk past a NULLB in a batch of several groups may have masked hits of another
        // group's unit (module doc, item 3): the units of each chunk again, in group order on
        // one table, with the search's failure. Else the units' refusal is the search's.
        if passed.load(AtomicOrdering::Relaxed) && self.x_group_count() > 1 {
            self.x_units_in_order(&units, inp_all, &mut outputs);
        } else if failure.failed() {
            if let Some(error) = failure.take() {
                inp_all.failure.fail(inp_all.s_idx, error);
            }
        }
        add_ns(&NS_UNITS, t_units);
        XSubjectUnits {
            frame_chunks,
            units,
            outputs,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:595-599
    // NCBI's one table holds the cell of a diagonal for every context of the batch, so a right
    // walk that went on past a NULLB (extension/query_batch.rs) masks the hits of the next
    // contexts under its stop. A unit's table holds only its group's contexts; where the stop lies
    // in another group's contexts, that group's unit does not see it. This runs the units of
    // every chunk of the subject again, a chunk's units one after the other in group order, on
    // one table of all the batch's diagonals (module doc, item 3), and replaces their lists.
    /// The units of `units` again, each chunk's in group order on one table of the batch's
    /// diagonals, into `outputs`; their refusal goes to `inp.failure`.
    #[cold]
    #[inline(never)]
    fn x_units_in_order(&self, units: &[XUnit], inp: &XSubject<'_>, outputs: &mut XUnitOutputs) {
        let mut scratch = XScratch::new(false);
        let sort = x_kmerge::on();
        let mut start = 0;
        while start < units.len() {
            let (f, k) = (units[start].f, units[start].k);
            let end = start
                + units[start..]
                    .iter()
                    .take_while(|unit| unit.f == f && unit.k == k)
                    .count();
            let chunk_units = &units[start..end];
            match outputs {
                XUnitOutputs::Full(lists) => {
                    let mut again: Vec<Vec<InitHSP>> =
                        chunk_units.iter().map(|_| Vec::new()).collect();
                    self.x_scan_in_order(chunk_units, inp, &mut scratch, &mut again);
                    for (slot, mut list) in lists[start..end].iter_mut().zip(again) {
                        if sort {
                            x_kmerge::XUnitSort::x_sort_in(&mut list, &mut scratch);
                        }
                        *slot = list;
                    }
                }
                XUnitOutputs::Compact(lists) => {
                    let mut again: Vec<Vec<XUnitHsp>> =
                        chunk_units.iter().map(|_| Vec::new()).collect();
                    self.x_scan_in_order(chunk_units, inp, &mut scratch, &mut again);
                    for (slot, records) in lists[start..end].iter_mut().zip(again) {
                        let mut list = x_unit_copy(&records, false);
                        if sort {
                            x_kmerge::XUnitSort::x_sort_in(&mut list, &mut scratch);
                        }
                        *slot = list;
                    }
                }
            }
            start = end;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-516
    // The scan loop of `scan_unit` (no prefetch, no branch-free loop, no sentinel walks, no probe)
    // for the units of one chunk (`chunk_units`, in group order), each group's tables (its whole
    // table, or LOSAT_X_TBXQSLICE's slices in slice order) one after the other, on one fresh table
    // whose cells are the raw diagonals of all the groups (`XUnmasked` from the lowest group's
    // bound, the width of `x_unit_width` over the groups' span). The hits of one diagonal come
    // from the groups in query order, so in ascending subject offset: each cell sees the
    // reference's hits in the reference's order (module doc, items 1-3; x_tbx_qslice.rs item 2),
    // and each unit's list holds its group's init HSPs in its scan order (item 4).
    /// Scans the units of one chunk in group order on one table, each into its sink.
    fn x_scan_in_order<S: XInitSink>(
        &self,
        chunk_units: &[XUnit],
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        sinks: &mut [S],
    ) {
        let tables = |g: usize| -> &[XGroup] {
            match &self.sliced {
                Some(sliced) => sliced[g].slices(),
                None => std::slice::from_ref(&self.groups[g]),
            }
        };
        let with_entries = |g: usize| {
            tables(g)
                .iter()
                .filter(|table| table.lookup.longest_chain > 0)
        };
        let bounds = chunk_units
            .iter()
            .flat_map(|unit| with_entries(unit.g))
            .map(|table| (table.q_min, table.q_max))
            .reduce(|(lo, hi), (q_min, q_max)| (lo.min(q_min), hi.max(q_max)));
        let Some((q_min, q_max)) = bounds else {
            return;
        };
        let wordsize = inp.wordsize;
        let window = inp.window;
        let unit0 = &chunk_units[0];
        let s_frame = &inp.s_frames_preliminary[unit0.f];
        let subject_full = &s_frame.aa_seq;
        let subject_all = &subject_full[1..subject_full.len() - 1];
        let chunk = unit0.chunk;
        let chunk_end = chunk.offset.saturating_add(chunk.length);
        let subject = &subject_all[chunk.offset..chunk_end];
        if subject.len() < wordsize as usize {
            return;
        }
        let width = (q_max - q_min) as usize + chunk.length;
        let cells = XUnmasked::of(q_min.wrapping_sub(chunk.length as u32 - 1), 0);
        scratch.diag.clear();
        scratch.diag.resize(width, DiagStruct::default());
        let h = XHitCtx {
            diag_offset: window,
            window,
            wordsize,
            // `get_context_idx` reads only the batch's context bounds, which every group and
            // slice table carries.
            lookup: &tables(unit0.g)[0].lookup,
            ctx_of: &self.ctx_of,
            contexts: inp.contexts,
            cutoff_scores: inp.cutoff_scores,
            x_dropoff_per_context: inp.x_dropoff_per_context,
            subject,
            subject_sent: &[],
            s_f_idx: unit0.f,
            s_frame: s_frame.frame,
            s_idx: inp.s_idx,
            s_len: inp.s_len,
            sc: self.sc,
            query_batch: inp.query_batch(),
        };
        let base_seq_ranges: [(i32, i32); 1] = [(0, chunk.length as i32)];
        let scan_interiors = tblastx_scan_interiors(chunk.length, None);
        for (unit, sink) in chunk_units.iter().zip(sinks.iter_mut()) {
            for table in with_entries(unit.g) {
                // NCBI sizes the offset array from the table it scans (c++/src/algo/blast/core/lookup_wrap.c:264-266).
                let offset_array_size: i32 = OFFSET_ARRAY_SIZE + table.lookup.longest_chain;
                let offset_pairs_len =
                    crate::algorithm::tblastx::x_scan_bf::x_offset_pairs_len(offset_array_size);
                if scratch.offset_pairs.len() < offset_pairs_len {
                    scratch
                        .offset_pairs
                        .resize(offset_pairs_len, OffsetPair::default());
                }
                for &(interior_start, interior_end) in &scan_interiors {
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
                    // Initialize `scan_range` to range index 0, start offset, and last word position.
                    let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];
                    // Scan batches until the current subject offset advances past the end bound.
                    while scan_range[1] <= scan_range[2] {
                        if !x_unit_batch::<S, XNoProbe, XUnmasked, _, false, false, 0>(
                            table,
                            subject,
                            &seq_ranges,
                            offset_array_size,
                            &mut scan_range,
                            scratch,
                            cells,
                            width,
                            &h,
                            sink,
                            &mut XNoProbe,
                        ) {
                            // The reference's safety guard (`word_finder::scan_next_hits`).
                            break;
                        }
                    }
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:561-562,578-581
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:234
    // The seed stage of one subject, second half: the list the reference builds per chunk in order
    // (sort, conversion and offset adjustment, `word_finder::finish_chunk_scan`; chunk merge,
    // `frame_chunks::add_chunk_result`) and per frame in s_f_idx order (append and sort,
    // `frame_chunks::append_frame_hits`); `merge_subject_listed` runs those statements,
    // LOSAT_X_TBXKMERGE (x_kmerge.rs) computes the same list without the full sorts. It has no
    // `self`: it reads the units' own lists and the subject's inputs, never a group table, so the
    // tables may already be freed when it runs.
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
                    )
                    .expect(TWO_HIT_CHUNK_MERGE_IS_DEFINED);
                    if inp.subject_besthit {
                        frame_subject_besthit(&mut frame_ungapped_hits);
                    }
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2876-2881
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
                let merged = !hits.is_empty();
                if !hits.is_empty() && frame_ungapped_hits.is_empty() {
                    frame_ungapped_hits = hits;
                } else if !hits.is_empty() {
                    merge_tblastx_subject_chunk_hits(
                        &mut frame_ungapped_hits,
                        hits,
                        chunk.offset,
                        chunk.overlap,
                    )
                    .expect(TWO_HIT_CHUNK_MERGE_IS_DEFINED);
                }
                if merged && inp.subject_besthit {
                    frame_subject_besthit(&mut frame_ungapped_hits);
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
    // Scheduling only: the units run on the search pool (an indexed map, results in unit order),
    // or in order on the calling thread when the search has one thread. `unit_fn` is
    // `Self::run_unit` (today's lists) or `Self::run_unit_compact` (LOSAT_X_SEEDCOMPACT), each
    // with or without LOSAT_X_TBXPREFETCH's hint (`PF`), a function item, so this function is
    // compiled once for each. The scratches of the pool carry LOSAT_X_SEEDUNMAP as it acts in this
    // process (`seed_unmap_active`, `XScratch::unmap`); the serial scratch never maps (section
    // above `X_MAPPED_BYTES`, rule 3: at one thread nothing is retained, mapping only costs).
    // EXPERIMENT (LOSAT_X_TBXKMERGE, x_kmerge.rs item 1): with `sort` each unit's list is sorted
    // in the unit's own task, after its scan, by the comparator of the chunk's init sort.
    fn run_units<T: Send + Sync + x_kmerge::XUnitSort>(
        &self,
        units: &[XUnit],
        inp: &XSubject<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
        unit_fn: impl Fn(&Self, &XUnit, &XSubject<'_>, &mut XScratch) -> T + Sync,
        sort: bool,
    ) -> Vec<T> {
        let run = |unit: &XUnit, scratch: &mut XScratch| {
            let mut list = unit_fn(self, unit, inp, scratch);
            if sort {
                list.x_sort_in(scratch);
            }
            list
        };
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() {
            let unmap = seed_unmap_active();
            // EXPERIMENT (LOSAT_X_TBXBOARD): the units on the search's board, by batch and subject
            // order (tblastx/x_work_board.rs); the same unit function on each unit, so the same
            // lists in unit order.
            if let Some((board, batch)) = crate::algorithm::tblastx::x_work_board::current() {
                let results: Vec<OnceLock<T>> = (0..units.len()).map(|_| OnceLock::new()).collect();
                let item =
                    |i: usize, taker: &mut crate::algorithm::tblastx::x_work_board::XTaker| {
                        let scratch = taker
                            .scratch
                            .get_or_insert_with(|| Box::new(XScratch::new(unmap)))
                            .downcast_mut::<XScratch>()
                            .expect("a board taker's scratch is the unit path's");
                        let list = run(&units[i], scratch);
                        if results[i].set(list).is_err() {
                            unreachable!("a unit is taken once");
                        }
                    };
                board.run_units(batch, inp.s_idx, units.len(), &item);
                return results
                    .into_iter()
                    .map(|cell| cell.into_inner().expect("every unit ran"))
                    .collect();
            }
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:229-235
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

    // The scan loop of the reference (`scan_chunk_plain` with `word_finder::scan_next_hits`) on one
    // unit: the group's sub-lookup, the frame chunk, a fresh table (`DiagStruct::default()`,
    // diag_offset = window, as `process_subject` resets its table per subject) indexed by the raw
    // diagonal minus the group's lowest one, without the mask (module doc, item 3). The body of
    // what was `run_unit`; the list it saves into is the caller's (`XInitSink`).
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
        self.scan_unit_hits::<S, P, PF, false, 0>(unit, inp, scratch, init_hsps, probe);
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF, x_tbx_hitbf.rs): with `BF` the hits of a batch go through
    // `x_hits_bf`, which takes the too-far and overlap outcomes of `x_hit` without a branch and
    // hands the others to `x_hit`; without it this is the loop before the switch. `BF` runs only
    // with `XNoProbe` (the probes count every outcome) and where `x_tbx_hitbf::fits` holds.
    /// `scan_unit`, with `BF` the branch-free hit loop of LOSAT_X_TBXHITBF and with `SENT` (a word
    /// length, 0: none) the extensions of `TbxScore::sent_two_hit_of` (LOSAT_X_TBXUNGAPSENT).
    fn scan_unit_hits<
        S: XInitSink,
        P: XProbe,
        const PF: bool,
        const BF: bool,
        const SENT: usize,
    >(
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
        const SENT: usize,
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

    /// The unit of the context group `g` on the chunk `k` (`chunk`) of the subject frame `f`.
    pub(super) fn x_unit(
        &self,
        inp: &XSubject<'_>,
        f: usize,
        k: usize,
        chunk: SubjectChunk,
        g: usize,
    ) -> XUnit {
        let sent_frame = x_tbx_sent_frame(self.sent, &inp.s_frames_preliminary[f].aa_seq);
        XUnit {
            f,
            k,
            chunk,
            g,
            sent_frame,
        }
    }

    // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): never on wasm32 (line W keeps its walks there).
    /// Whether the extensions of `unit` take the walks the sentinels stop: its frame may
    /// (`XUnit::sent_frame`) and the unit's chunk is its whole subject frame.
    #[inline(always)]
    pub(super) fn x_unit_sent(&self, unit: &XUnit, inp: &XSubject<'_>) -> bool {
        cfg!(not(target_arch = "wasm32"))
            && x_tbx_sent_subject(
                unit.sent_frame,
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
        const SENT: usize,
    >(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
        probe: &mut P,
        mask: u32,
    ) {
        // LOSAT_X_TBXUNGAPSENT: the walks of `SENT`-letter words only on the units `x_unit_sent`
        // accepts (`x_hit`'s SAFETY; once per unit), of the search's word length (the dispatches of
        // x_tbx_hitbf.rs, `scan_unit_or_bf`, choose `SENT` and check it, outside this function so
        // that its loop keeps its registers).
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
            // Initialize `scan_range` to range index 0, start offset, and last word position.
            let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];
            // Scan batches until the current subject offset advances past the end bound.
            while scan_range[1] <= scan_range[2] {
                // One batch per call (`x_unit_batch`: its own function on wasm32 only).
                if !x_unit_batch::<S, P, C, _, PF, BF, SENT>(
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
                    // The reference's safety guard (`word_finder::scan_next_hits`).
                    break;
                }
            }
        }
        probe.unit_end();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:501-509
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
    M: TbxScore,
    const PF: bool,
    const BF: bool,
    const SENT: usize,
>(
    group: &XGroup,
    subject: &[u8],
    seq_ranges: &[(i32, i32)],
    offset_array_size: i32,
    scan_range: &mut [i32; 3],
    scratch: &mut XScratch,
    cells: C,
    width: usize,
    h: &XHitCtx<'_, M>,
    init_hsps: &mut S,
    probe: &mut P,
) -> bool {
    let prev_scan_left = scan_range[1];
    // Scan the next batch of subject word hits via `scansub` into `offset_pairs`.
    // EXPERIMENT (LOSAT_X_TBXSCANBF): the branch-free scan with the switch.
    let hits = group.scan(
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
                    x_hit::<S, P, _, SENT>(
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
    // Process each hit returned in this batch in scan order.
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
        // Map the hit diagonal (`query_offset` minus `subject_offset`) into the unit table.
        // Unit cell: the raw diagonal minus the unit's lowest one (`XUnmasked`: no mask;
        // LOSAT_X_TBXUNITMASK's `XMasked`: masked to the table).
        let diag_coord = cells.cell(query_offset, subject_offset);
        debug_assert!(diag_coord < width, "unit diagonal outside its table");
        // SAFETY: q_min <= query_offset <= q_max (the group's entries) and
        // subject_offset < s_len, so an unmasked diag_coord is below width = scratch.diag.len();
        // a masked one is at most the mask, width - 1.
        let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };
        probe.at(diag_coord, subject_offset);
        let outcome = x_hit::<S, P, _, SENT>(
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
// The unit path's own form of the per-hit body of the reference loop
// (`word_finder::two_hit_one_hit`): its two-hit statements on the unit's unmasked cell, without
// its trace, debug, diagnostics and HSP-saving statistics (the gates send those runs to the
// reference loop), with `continue` written as `return`; the extension, the sink and the probe are
// this path's own. As in the reference, the cell is updated before
// the init HSP is saved (both read only local values, so the order is immaterial). The init HSP
// goes to the caller's sink (today's `Vec<InitHSP>`, or the records of LOSAT_X_SEEDCOMPACT, section
// above `XUnitHsp`), by the same push statement.
/// The two-hit test and extension of one hit on its cell (`SENT`: LOSAT_X_TBXUNGAPSENT, the
/// sentinel walks of `SENT`-letter words, `TbxScore::sent_two_hit_of`; 0: none).
#[inline(always)]
fn x_hit<S: XInitSink, P: XProbe, M: TbxScore, const SENT: usize>(
    diag_entry: &mut DiagStruct,
    query_offset: u32,
    subject_offset: u32,
    h: &XHitCtx<'_, M>,
    init_hsps: &mut S,
    probe: &mut P,
) -> XHit {
    let diag_offset = &h.diag_offset;
    let window = h.window;
    let wordsize = h.wordsize;
    // When the diagonal's `flag` is set, skip hits covered by a prior right extension.
    if diag_entry.flag() != 0 {
        // Skip the hit if the offset-adjusted subject position precedes `last_hit`
        // (signed `Int4` comparison in NCBI; Rust compares unsigned).
        let subject_plus_offset = subject_offset.wrapping_add(*diag_offset as u32);
        if subject_plus_offset < diag_entry.last_hit() as u32 {
            return XHit::Masked;
        }
        // Otherwise record the new offset-adjusted `last_hit` and clear `flag`.
        diag_entry.set_last_hit(subject_plus_offset as i32);
        diag_entry.set_flag(0);
        XHit::FlagReset
    }
    // Diagonal flag is clear: test distance from the previous hit.
    else {
        // Recover the previous hit's subject offset by subtracting `diag_offset`.
        let last_hit = diag_entry.last_hit() - *diag_offset;
        // Compute the subject distance from `last_hit` to `subject_offset`.
        let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;

        // Outside the two-hit window: record this hit as the new `last_hit` and return.
        if diff >= window {
            diag_entry.set_last_hit(subject_offset.wrapping_add(*diag_offset as u32) as i32);
            return XHit::TooFar;
        }

        // Overlapping word hit: ignore without updating `last_hit`.
        if diff < wordsize {
            return XHit::Overlap;
        }

        // Look up the query context containing `query_offset` (as `BSearchContextInfo`).
        let ctx_idx = h.context_of(query_offset);
        // SAFETY: context_of returns an index below num_contexts == contexts.len().
        let ctx = unsafe { h.contexts.get_unchecked(ctx_idx) };
        let q_raw = query_offset.wrapping_sub(ctx.frame_base as u32) as usize;
        let query_full = &ctx.aa_seq;
        let query = &query_full[1..query_full.len() - 1];

        // Reject hits whose earlier partner lies before the start of the current query context.
        let q_minus_diff = query_offset.wrapping_sub(diff as u32);
        if q_minus_diff < ctx.frame_base as u32 {
            diag_entry.set_last_hit(subject_offset.wrapping_add(*diag_offset as u32) as i32);
            return XHit::ContextBoundary;
        }

        // Load the score cutoff and X-dropoff parameters for the current query context.
        // SAFETY: ctx_idx < contexts.len() == cutoff_scores.len() == x_dropoff_per_context.len().
        let cutoff = unsafe { *h.cutoff_scores.get_unchecked(ctx_idx) };
        let x_dropoff = unsafe { *h.x_dropoff_per_context.get_unchecked(ctx_idx) };

        // Run two-hit ungapped extension (as `s_BlastAaExtendTwoHit`) between the prior hit end
        // and the current hit.
        // Diagnostic (LOSAT_TIMING, LOSAT_X_STATS): the probe around the extension
        // (x_unit_probe.rs; `XNoProbe` does nothing). No NCBI counterpart.
        let mark = probe.ext_begin();
        // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): with `SENT` the walks the sentinels stop, of
        // `SENT`-letter words.
        // Another matrix's extension takes its table walks (`extend_hit_two_hit_with`).
        // NCBI's walks run on the query batch: past the context's ends they go on into it
        // (extension/query_batch.rs; aa_ungapped.c:841, 897).
        let in_batch = h.query_batch.at(ctx_idx, h.s_idx);
        let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) =
            if !M::BLOSUM62 && SENT == 0 {
                extend_hit_two_hit_with(
                    h.sc,
                    query,
                    in_batch,
                    h.subject,
                    (last_hit + wordsize) as usize,
                    subject_offset as usize,
                    q_raw,
                    x_dropoff,
                    false,
                )
            } else if SENT == 0 {
                extend_hit_two_hit_fitting(
                    query,
                    in_batch,
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
                // is the buffer of a frame `x_tbx_sent_subject` accepted; `SENT` is the search's
                // word length (the dispatch `scan_unit_or_bf` chooses and checks it).
                unsafe {
                    h.sc.sent_two_hit_of::<SENT>(
                        query_full,
                        in_batch,
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
            h.sc,
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

        // On a right extension, set `flag` to 1 and `last_hit` to the extension end `s_last_off`
        // less `wordsize - 1`, biased by `diag_offset`; otherwise `last_hit` is the subject offset
        // biased by `diag_offset` and `flag` stays 0.
        if right_extend {
            diag_entry.set_flag(1);
            diag_entry.set_last_hit(s_last_off - (wordsize - 1) + *diag_offset);
        } else {
            diag_entry.set_last_hit(subject_offset.wrapping_add(*diag_offset as u32) as i32);
        }

        // Save the initial HSP (as `BlastSaveInitHsp`) whenever `score` meets or exceeds the context cutoff.
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

// No NCBI counterpart: LOSAT_X_TBXFEWGROUPS's decision for a batch (scheduling only; module doc).
/// Whether a batch that `wanted` would cut less than `X_MAX_GROUPS` does extends over at most
/// `X_FEW_GROUPS_LETTERS` letters; None where `wanted` changes nothing.
fn x_takes_few_groups(contexts: &[QueryContext], wanted: usize) -> Option<bool> {
    let n = contexts.len();
    if wanted >= X_MAX_GROUPS || n <= wanted {
        return None;
    }
    let letters =
        contexts[n - 1].frame_base + contexts[n - 1].aa_len as i32 - contexts[0].frame_base;
    Some(letters as usize <= X_FEW_GROUPS_LETTERS)
}

// No NCBI counterpart: how the contexts are grouped into units (scheduling only; module doc).
/// The first context of each group: one group per context up to `X_MAX_GROUPS` contexts, else
/// `X_MAX_GROUPS` contiguous groups balanced by length; `wanted` groups in place of `X_MAX_GROUPS`
/// for a batch that extends over at most `X_FEW_GROUPS_LETTERS` letters (LOSAT_X_TBXFEWGROUPS,
/// `x_groups_wanted`).
fn group_starts(contexts: &[QueryContext], wanted: usize) -> Vec<usize> {
    let n = contexts.len();
    let max_groups = if x_takes_few_groups(contexts, wanted) == Some(true) {
        wanted
    } else {
        X_MAX_GROUPS
    };
    if n <= max_groups {
        return (0..n).collect();
    }
    let total: usize = contexts.iter().map(|c| c.aa_len + 1).sum();
    let mut starts = vec![0usize];
    let mut cumulative = 0usize;
    for (c, ctx) in contexts.iter().enumerate() {
        if c > 0 && cumulative * max_groups >= starts.len() * total {
            starts.push(c);
            if starts.len() == max_groups {
                break;
            }
        }
        cumulative += ctx.aa_len + 1;
    }
    starts
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:282-299,337-360
// One table per context group in the layout above: each backbone chain of the batch table
// filtered to the group's contexts (`BSearchContextInfo` ranges of `frame_bases`), order kept;
// pv bit, inline entries or overflow cursor (cells in index order), and longest_chain rebuilt.
fn build_groups(
    lookup: &BlastAaLookupTable,
    contexts: &[QueryContext],
    wanted: usize,
) -> Vec<XGroup> {
    let starts = group_starts(contexts, wanted);
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
                cell.set_cursor(b.overflow.len());
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
            XGroup::new(
                lookup,
                if b.longest_chain > 0 { b.q_min } else { 0 },
                b.q_max,
            )
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
#[path = "x_pair_par_tests.rs"]
mod tests;
