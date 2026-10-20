//! EXPERIMENT (LOSAT_X_TBXQSLICE / LOSAT_X_TBXQSLICESHADOW): the scan of a unit of the unit path
//! (x_pair_par.rs) in slices of its query group, each slice with its own lookup table, all of them
//! on the unit's one diagonal table. A child module of x_pair_par.rs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-516
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
//! ```
//! A unit scans its subject chunk once, and at each subject position `s` the lookup gives the hits
//! of every query offset of the group, so the cells it touches, `q - s`, lie anywhere in the
//! table (the group's span plus the chunk, 4 bytes a cell). Where the group is a long query frame
//! (Sakai, 1.83 M letters a frame) one table is about 7 MB: it fits the L3 at one thread, but not
//! 8 or 16 of them, and the units' cycles grow by 39 % (8 threads) and 82 % (16 threads) for the
//! same instructions (v4 line A, llc-sakai-20261018/DESIGN_qslice.md §1). Here such a group is
//! cut into slices of consecutive query offsets, one lookup table each (the group's table with
//! every chain filtered to the slice), and the unit scans its chunk once per slice, in slice order,
//! on the same table, which is not cleared between slices. While slice j `[a, b)` is scanned, the
//! cells touched at subject position `s` are `[a - s, b - s)`: a window of the slice's width that
//! moves one cell per subject position. With slices of at most the switch's letters (built-in
//! `X_QSLICE_LETTERS`) the window fits the core's L2.
//!
//! # Why the result is that of `scan_unit`
//!
//! Following the items of x_pair_par.rs's proof:
//! 1. A hit's decision reads and writes only its own cell, and depends on the earlier hits of that
//!    cell (item 1). The extension reads only the two sequences and the matrix.
//! 2. Same hits, same order on each cell. A slice's chain is the group's chain filtered to the
//!    slice, order kept, so the slices' scans emit exactly the group's hits, each once (a hit
//!    (q, s) comes from the slice that holds q). The hits of one cell (diagonal d) are (q, q - d)
//!    with distinct s; those from slice j have s in `[a_j - d, b_j - d)`, below every s of slice
//!    j + 1, and each slice's scan emits its hits in ascending s. So each cell sees its hits in
//!    ascending s, the order of `scan_unit`. Only hits of different cells change their order.
//! 3. Same start. The table is set up once per unit, as `scan_unit` sets it up (the group's bounds,
//!    hence the same `lo` and width); every cell starts fresh, and by 2 it then takes the same
//!    updates. The table at the end of the unit is `scan_unit`'s, cell for cell.
//! 4. Same list after the chunk's sort. The unit saves the same init HSPs (1-3), each where
//!    `scan_unit` saves it relative to the other init HSPs of its cell; init HSPs of different
//!    cells may change order. The chunk's list is concatenated in unit order and sorted with the
//!    stable `score_compare_match` (LOSAT_X_TBXKMERGE: each unit's list sorted alone, then merged
//!    with the same comparator, x_kmerge.rs); records the comparator calls equal share `s_start`
//!    and the absolute `q_start`, so they lie on one cell of one unit, in the same relative order.
//!    The stable sorts therefore return `scan_unit`'s sorted lists element for element.
//! 5. The tables. `XPairLutQuery::x_build_group` builds a table from any list of lookup locations
//!    (x_pair_lut.rs, items 1-5: the full chains filtered to the offsets of the words of those
//!    locations, order kept). A slice's locations are the group's, each cut at the slice bounds by
//!    word start (`x_slice_locations`): a word lies in one location of the list and starts in one
//!    slice, and whether a word is added depends only on its own letters
//!    (`blast_lookup_index_query_exact_matches`), so the slice's exact backbone is the group's
//!    filtered to the slice, order kept, and so is its table. Where the reference builder runs
//!    (`build_groups`), the slices are its filter at the slice bounds (`build_tables`). The stub's
//!    `longest_chain` sums each cell over a partition of the entries, the slices as well as the
//!    groups. Each slice's scan sizes its offset array from its own table, as NCBI sizes it from
//!    the table it scans (lookup_wrap.c:264-266); where a scan call ends does not matter (item 2 of
//!    x_pair_par.rs).
//!
//! # The rule
//!
//! `LOSAT_X_TBXQSLICE` is a value switch. Unset (or `=1`), a group is cut into slices of at most
//! `X_QSLICE_LETTERS` letters where its query offsets (its contexts' extent, sentinels included)
//! span more letters and, times the search's threads, at least `X_QSLICE_MIN_SPAN_THREADS`: there
//! the diagonal tables of all threads at once (4 bytes a cell) crowd the L3; below it the whole
//! group's scan (with LOSAT_X_TBXPREFETCH where the table is large) is as fast or faster (Sakai at
//! one thread: the slices x1.015). `=N` (2 to `X_QSLICE_MAX_LETTERS`) cuts every group of more
//! than N letters into slices of at most N, whatever the threads (tuning and tests; with at most
//! `X_QSLICE_MAX_SLICES` slices a group, wider where N is small). A group is cut into the fewest
//! equal slices (`X_QSLICE_MAX_SLICES` at most; `x_partition`); the cut does not depend on the
//! subjects. A batch in which no group is cut runs as before (`scan_unit`, whose machine code does
//! not change); in a batch with a cut group every group's units take this path (a group that is not
//! cut is one slice). The size of a slice changes no output. x86-64 only (where it was measured).
//! The probes of LOSAT_TIMING / LOSAT_X_STATS (x_unit_probe.rs) and the shadow of
//! LOSAT_X_TBXEXTBATCH (which compares `scan_unit` with the batched scan of the whole group,
//! unit by unit) need the group's table: with them nothing is cut.
//!
//! Interactions, read by the window (the widest slice's span of entries) in place of the group:
//! LOSAT_X_TBXPREFETCH prefetches where the window reaches its size, and LOSAT_X_TBXEXTBATCH holds
//! the extensions back where the window reaches its span (both rules unchanged otherwise). The
//! pending marks of LOSAT_X_TBXEXTBATCH stay valid across slices (a mark is resolved before the
//! next hit on its cell is tested, from whichever slice); the batch is flushed at the end of the
//! unit.
//!
//! `LOSAT_X_TBXQSLICESHADOW=1` builds the whole groups as well and runs `scan_unit` with them on
//! every unit of a cut batch, on its own table: it stops at the first difference of the final
//! tables or of the unit's init HSPs after the chunk's stable sort (`XUnitSort::x_sort`).

use super::x_kmerge::XUnitSort;
use super::x_tbx_extbatch::{self, XPendBatch, XShadowList, XShadowRef, BATCH};
use super::*;

// No NCBI counterpart: counters of LOSAT_X_TBXQSLICE (batches with a cut group, groups cut, slice
// tables; the units and slice scans in shadow mode only) and of its shadow; they do not change any
// value NCBI computes.
static QS_BATCHES: AtomicU64 = AtomicU64::new(0);
static QS_GROUPS_CUT: AtomicU64 = AtomicU64::new(0);
static QS_SLICES: AtomicU64 = AtomicU64::new(0);
static QS_UNITS: AtomicU64 = AtomicU64::new(0);
static QS_SCANS: AtomicU64 = AtomicU64::new(0);
static QS_SHADOW_UNITS: AtomicU64 = AtomicU64::new(0);
static QS_SHADOW_RECORDS: AtomicU64 = AtomicU64::new(0);

/// The most letters of a slice when the switch is not given a number (128 Ki: a window of
/// 512 KiB, which fits the core's L2 with the slice's backbone).
const X_QSLICE_LETTERS: usize = 1 << 17;
/// The rule's bound when the switch is not given a number: a group is cut only where its span times
/// the search's threads reaches this many letters (6 Mi: 24 MiB of diagonal tables at once; on the
/// workstation 4.3 M lost 0.5 % and 7.3 M gained 28 %).
const X_QSLICE_MIN_SPAN_THREADS: u64 = 6 << 20;
/// The largest number the switch takes; a value outside 2..=X_QSLICE_MAX_LETTERS other than `0`,
/// or a non-numeric one, gives the built-in rule (`xswitch::value`).
const X_QSLICE_MAX_LETTERS: usize = 1 << 30;
/// The most slices of one group: each slice table has its own backbone (about 460 KB), so a small
/// value of the switch on a long group widens the slices instead.
const X_QSLICE_MAX_SLICES: u64 = 64;

/// The switch's rule: the most letters of a slice, and the least span times threads of a group that
/// is cut (0: every group longer than a slice).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct XSliceRule {
    letters: usize,
    min_span_threads: u64,
}

// No NCBI counterpart: which scan the units of the unit path use (scheduling only). Read once per
// process, as the other modes of x_pair_par.rs.
/// LOSAT_X_TBXQSLICE: a strict-set value switch (`utils/xswitch.rs`). `None` when it is off (`0`,
/// or `LOSAT_X_DEFAULT=0` and it is not set) and on targets other than x86-64; the built-in rule
/// when it is not given a number (the built-in passed to `xswitch::value` is 0, outside the range);
/// otherwise slices of at most the number given, on every group.
fn slice_rule() -> Option<XSliceRule> {
    static RULE: OnceLock<Option<XSliceRule>> = OnceLock::new();
    *RULE.get_or_init(|| {
        if !cfg!(target_arch = "x86_64") {
            return None;
        }
        match crate::utils::xswitch::value("LOSAT_X_TBXQSLICE", 2..=X_QSLICE_MAX_LETTERS, 0)? {
            0 => Some(XSliceRule {
                letters: X_QSLICE_LETTERS,
                min_span_threads: X_QSLICE_MIN_SPAN_THREADS,
            }),
            letters => Some(XSliceRule {
                letters,
                min_span_threads: 0,
            }),
        }
    })
}

// No NCBI counterpart: a comparison of two LOSAT code paths.
/// LOSAT_X_TBXQSLICESHADOW (set to any value): the whole groups are built as well and every unit of
/// a cut batch runs `scan_unit` with them, compared.
pub(super) fn shadow() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_TBXQSLICESHADOW").is_some())
}

// No NCBI counterpart: the rule of the switch for one batch (module doc). The groups are those of
// `group_starts`; group g holds the query offsets from its first context's `frame_base` to the
// next group's (the last group: to the end of its last context). `threads` are the search pool's.
/// For each group, the first query offset of each of its slices (the group's own first offset
/// first), or None when no group is cut: the switch is off, a probe runs, or no group passes the
/// rule.
pub(super) fn partition(
    contexts: &[QueryContext],
    starts: &[usize],
    threads: usize,
) -> Option<Vec<Vec<i32>>> {
    let rule = slice_rule()?;
    if x_unit_probe::mode() != 0 || x_tbx_extbatch::shadow() {
        return None;
    }
    let min_len = rule.min_span_threads.div_ceil(threads.max(1) as u64);
    x_partition(&x_group_spans(contexts, starts), rule.letters, min_len)
}

// No NCBI counterpart: the rule's input (module doc).
/// The query offsets `[a, b)` of each group (`group_starts`'s first contexts `starts`).
fn x_group_spans(contexts: &[QueryContext], starts: &[usize]) -> Vec<(i32, i32)> {
    let end = contexts
        .last()
        .map(|c| c.frame_base + c.aa_len as i32)
        .unwrap_or(0);
    starts
        .iter()
        .enumerate()
        .map(|(g, &c)| {
            let b = starts
                .get(g + 1)
                .map(|&next| contexts[next].frame_base)
                .unwrap_or(end);
            (contexts[c].frame_base, b)
        })
        .collect()
}

// No NCBI counterpart: the rule's arithmetic (module doc).
/// `partition` for groups of query offsets `spans` and slices of at most `width` letters: a group
/// of `len > width` letters and at least `min_len` is cut into `k = ceil(len / width)` slices (at
/// most `X_QSLICE_MAX_SLICES`) starting at `a + floor(j len / k)`, j = 0..k; None when no group is
/// cut.
fn x_partition(spans: &[(i32, i32)], width: usize, min_len: u64) -> Option<Vec<Vec<i32>>> {
    let mut cut = false;
    let parts = spans
        .iter()
        .map(|&(a, b)| {
            let len = (i64::from(b) - i64::from(a)).max(0) as u64;
            let k = if len >= min_len {
                len.div_ceil(width.max(1) as u64)
                    .clamp(1, X_QSLICE_MAX_SLICES)
            } else {
                1
            };
            cut |= k > 1;
            (0..k).map(|j| a + (j * len / k) as i32).collect()
        })
        .collect();
    cut.then_some(parts)
}

/// One group of a batch with a cut group: its slices, in query order, each an `XGroup` (the slice's
/// table, and its entries' bounds), and the bounds of the group's entries, those of its whole
/// table.
pub(super) struct XSliced {
    q_min: u32,
    q_max: u32,
    /// The largest span of a slice's entries (`x_unit_span`), the window the rules read.
    window: usize,
    slices: Vec<XGroup>,
}

// No NCBI counterpart: the bounds of a cut group's tables (scheduling and the rules only).
impl XSliced {
    /// Whether the group has lookup entries (`longest_chain > 0` of its whole table).
    pub(super) fn has_entries(&self) -> bool {
        self.slices.iter().any(|s| s.lookup.longest_chain > 0)
    }

    /// The group's span, `x_unit_span` of its whole table.
    fn span(&self) -> usize {
        (self.q_max - self.q_min) as usize
    }

    /// The window: the largest span of a slice's entries.
    fn window(&self) -> usize {
        self.window
    }
}

// No NCBI counterpart: the stub's `longest_chain` reads every table (module doc, item 5).
/// The tables of every slice of every group, in order (the stub's partition of the entries).
pub(super) fn tables(sliced: &[XSliced]) -> impl Iterator<Item = &BlastAaLookupTable> {
    sliced
        .iter()
        .flat_map(|g| g.slices.iter().map(|s| &s.lookup))
}

// No NCBI counterpart: the slice tables arranged by group (scheduling only).
/// The groups of `parts` from their slice `tables`, in (group, slice) order. A group's bounds are
/// those of its slices with entries (0 and 0 without any, as for a whole table without entries).
fn x_assemble(parts: &[Vec<i32>], tables: Vec<XGroup>) -> Vec<XSliced> {
    let mut tables = tables.into_iter();
    let sliced: Vec<XSliced> = parts
        .iter()
        .map(|p| {
            let slices: Vec<XGroup> = tables.by_ref().take(p.len()).collect();
            assert_eq!(
                slices.len(),
                p.len(),
                "LOSAT_X_TBXQSLICE: a table per slice"
            );
            let full = || slices.iter().filter(|s| s.lookup.longest_chain > 0);
            XSliced {
                q_min: full().map(|s| s.q_min).min().unwrap_or(0),
                q_max: full().map(|s| s.q_max).max().unwrap_or(0),
                window: full().map(x_unit_span).max().unwrap_or(0),
                slices,
            }
        })
        .collect();
    assert!(
        tables.next().is_none(),
        "LOSAT_X_TBXQSLICE: a slice per table"
    );
    sliced
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:105-115
// ```c
//         seq = query->sequence + from;
//         word_target = seq + lut_word_length;
//
//         for (offset = from; offset <= to; offset++, seq++) {
//
//             if (seq >= word_target) {
//                 BlastLookupAddWordHit(backbone,
//                                       lut_word_length, charsize,
//                                       seq - lut_word_length,
//                                       offset - lut_word_length);
//             }
// ```
// A location [from, to] adds the words that start at from..=to - word_length + 1 and hold no
// invalid letter (module doc, item 5). Each location is cut at the slice starts by word start: the
// piece of slice j is [max(from, a_j), min(last word start, b_j - 1) + word_length - 1], which adds
// exactly the location's words that start in the slice; locations without a word and empty pieces
// are left out (the indexer skips them).
/// The lookup locations of each slice of a group whose locations are `locations` (list order) and
/// whose slices start at `parts`: every location cut at the slice bounds by word start.
fn x_slice_locations(
    locations: &[(i32, i32)],
    parts: &[i32],
    word_length: i32,
) -> Vec<Vec<(i32, i32)>> {
    let mut out: Vec<Vec<(i32, i32)>> = vec![Vec::new(); parts.len()];
    for &(from, to) in locations {
        let last_start = to - (word_length - 1);
        for (j, pieces) in out.iter_mut().enumerate() {
            let first = if j == 0 { from } else { from.max(parts[j]) };
            let last = match parts.get(j + 1) {
                Some(&next) => last_start.min(next - 1),
                None => last_start,
            };
            if first <= last {
                pieces.push((first, last + (word_length - 1)));
            }
        }
    }
    out
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:454-456
// ```c
//     BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
//                                       lookup->charsize, lookup->word_length,
//                                       query, location);
// ```
// The slices' tables at the lookup build site (LOSAT_X_PAIRLUT), each by `x_build_group` from the
// slice's locations (module doc, item 5), on the calling thread or the pool as the groups are.
/// The groups of `starts` cut at `parts`, built directly from the query. None when a lookup
/// location does not lie inside one group (the reference builder must run).
pub(super) fn lut_direct_sliced(
    q: &XPairLutQuery,
    starts: &[usize],
    parts: &[Vec<i32>],
    pool: &crate::utils::threading::SearchPool<'_>,
    on_pool: bool,
) -> Option<Vec<XSliced>> {
    let locations = q.x_group_locations(starts)?;
    let word_length = crate::algorithm::tblastx::lookup::LOOKUP_WORD_LENGTH as i32;
    let pieces: Vec<Vec<(i32, i32)>> = locations
        .iter()
        .zip(parts)
        .flat_map(|(group_locations, p)| x_slice_locations(group_locations, p, word_length))
        .collect();
    Some(x_assemble(
        parts,
        lut_direct_tables(q, &pieces, pool, on_pool),
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:337-343
// ```c
// for (i = 0; i < lookup->backbone_size; i++) {
//   /* if there are hits there, */
//   if (lookup->thin_backbone[i] ) {
//       Int4 * dest = NULL;
//       /* set the corresponding bit in the pv_array */
//       PV_SET(pv, i, PV_ARRAY_BTS);
// ```
// Where the reference builder made the batch's table: the slices are its chains filtered at the
// slice bounds, as `build_groups` filters them at the group bounds (`build_tables`).
/// The groups cut at `parts`, filtered out of the batch's full table.
pub(super) fn build_sliced(lookup: &BlastAaLookupTable, parts: &[Vec<i32>]) -> Vec<XSliced> {
    let bounds: Vec<i32> = parts.iter().flatten().copied().collect();
    x_assemble(parts, build_tables(lookup, &bounds))
}

// No NCBI counterpart: LOSAT_X_PAIRLUTSHADOW's comparison (x_pair_par.rs) of the slices built
// directly with those filtered out of the reference builder's table.
/// Asserts that both builders cut the same groups into the same slices with the same tables.
pub(super) fn lut_shadow_compare(direct: Option<&[XSliced]>, reference: Option<&[XSliced]>) {
    match (direct, reference) {
        (None, None) => {}
        (Some(direct), Some(reference)) => {
            assert_eq!(
                direct.len(),
                reference.len(),
                "LOSAT_X_PAIRLUTSHADOW: number of cut groups differs (direct vs reference)"
            );
            let mut t = 0usize;
            for (g, (a, b)) in direct.iter().zip(reference).enumerate() {
                assert_eq!(
                    (a.q_min, a.q_max, a.slices.len()),
                    (b.q_min, b.q_max, b.slices.len()),
                    "LOSAT_X_PAIRLUTSHADOW: group {g}: bounds or slices differ (direct vs reference)"
                );
                // `t` numbers the slice tables of the batch ("group" in the messages).
                for (x, y) in a.slices.iter().zip(&b.slices) {
                    x_assert_same_group(
                        t,
                        &x.lookup,
                        (x.q_min, x.q_max),
                        &y.lookup,
                        (y.q_min, y.q_max),
                    );
                    t += 1;
                }
            }
            // The slice tables count as the groups compared (`[X_SHADOW] PAIRLUT groups_compared=`).
            LUT_SHADOW_GROUPS.fetch_add(t as u64, AtomicOrdering::Relaxed);
        }
        _ => panic!("LOSAT_X_PAIRLUTSHADOW: one builder cut the groups and the other did not"),
    }
}

// No NCBI counterpart: counters of the summary.
/// Counts a batch with a cut group (LOSAT_X_TBXQSLICE's summary).
pub(super) fn count_batch(sliced: &[XSliced]) {
    QS_BATCHES.fetch_add(1, AtomicOrdering::Relaxed);
    let cut = sliced.iter().filter(|g| g.slices.len() > 1).count();
    QS_GROUPS_CUT.fetch_add(cut as u64, AtomicOrdering::Relaxed);
    let slices: usize = sliced.iter().map(|g| g.slices.len()).sum();
    QS_SLICES.fetch_add(slices as u64, AtomicOrdering::Relaxed);
}

/// The summary lines (printed by x_pair_par.rs's `print_summary`): the shadow's and, in shadow
/// mode, the counters.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub(super) fn print_summary() {
    if !shadow() {
        return;
    }
    let load = |c: &AtomicU64| c.load(AtomicOrdering::Relaxed);
    eprintln!(
        "[X_SHADOW] TBXQSLICE units_compared={} records_compared={} (records after the chunk sort and unit tables identical)",
        load(&QS_SHADOW_UNITS),
        load(&QS_SHADOW_RECORDS)
    );
    eprintln!(
        "[X_STATS] TBXQSLICE batches={} groups_cut={} slices={} units={} slice_scans={} rule={:?}",
        load(&QS_BATCHES),
        load(&QS_GROUPS_CUT),
        load(&QS_SLICES),
        load(&QS_UNITS),
        load(&QS_SCANS),
        slice_rule()
    );
}

impl XPairPar {
    // EXPERIMENT (LOSAT_X_TBXQSLICE): `run_subject_unit_fns` for a batch with a cut group. The rules
    // of LOSAT_X_TBXPREFETCH and LOSAT_X_TBXEXTBATCH read each unit's window (module doc). No NCBI
    // counterpart: dispatch only.
    /// The units of one subject of a batch with a cut group (`sliced`), run.
    pub(super) fn run_subject_unit_fns_sliced(
        &self,
        sliced: &[XSliced],
        units: &[XUnit],
        inp: &XSubject<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> XUnitOutputs {
        let windows = || units.iter().map(|unit| sliced[unit.g].window());
        let prefetch = x_prefetch_units(
            tbx_prefetch_min_kib(),
            windows().map(|w| w * std::mem::size_of::<DiagStruct>()),
        );
        let batched = x_tbx_extbatch::active(inp.window, inp.wordsize, windows());
        match (batched, prefetch) {
            (false, false) => self.run_units_sliced::<false, false>(units, inp, pool),
            (false, true) => self.run_units_sliced::<false, true>(units, inp, pool),
            (true, false) => self.run_units_sliced::<true, false>(units, inp, pool),
            (true, true) => self.run_units_sliced::<true, true>(units, inp, pool),
        }
    }

    // No NCBI counterpart: dispatch only.
    /// `run_units` with the sliced unit functions: batched extensions (`B`) and the prefetch
    /// (`PF`), today's lists or the compact records (LOSAT_X_SEEDCOMPACT).
    fn run_units_sliced<const B: bool, const PF: bool>(
        &self,
        units: &[XUnit],
        inp: &XSubject<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> XUnitOutputs {
        let sort = x_kmerge::mode() == 1;
        if seed_compact() {
            XUnitOutputs::Compact(self.run_units(
                units,
                inp,
                pool,
                Self::run_unit_compact_sliced::<B, PF>,
                sort,
            ))
        } else {
            XUnitOutputs::Full(self.run_units(
                units,
                inp,
                pool,
                Self::run_unit_sliced::<B, PF>,
                sort,
            ))
        }
    }

    // No NCBI counterpart: the unit's group (scheduling only).
    /// The group of `unit` in a batch with a cut group.
    fn x_sliced(&self, unit: &XUnit) -> &XSliced {
        &self
            .sliced
            .as_ref()
            .expect("LOSAT_X_TBXQSLICE: a sliced unit of a batch without slices")[unit.g]
    }

    // EXPERIMENT (LOSAT_X_TBXQSLICE): `run_unit` in slices.
    /// The saved init HSPs of one unit.
    #[inline(never)]
    fn run_unit_sliced<const B: bool, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<InitHSP> {
        let mut init_hsps: Vec<InitHSP> = Vec::new();
        self.scan_unit_sliced_or_shadow::<_, B, PF>(unit, inp, scratch, &mut init_hsps);
        init_hsps
    }

    // EXPERIMENT (LOSAT_X_TBXQSLICE): `run_unit_compact` in slices; LOSAT_X_SEEDUNMAP acts as there
    // (`x_records_for` on the group's span, `x_unit_copy`).
    /// The saved init HSPs of one unit as compact records.
    #[inline(never)]
    fn run_unit_compact_sliced<const B: bool, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<XUnitHsp> {
        let mut records = std::mem::take(&mut scratch.records);
        records.clear();
        let mapped = x_records_for(self.x_sliced(unit).span(), unit, scratch, &mut records);
        self.scan_unit_sliced_or_shadow::<_, B, PF>(unit, inp, scratch, &mut records);
        let exact = x_unit_copy(&records, mapped);
        scratch.records = records;
        exact
    }

    // EXPERIMENT (LOSAT_X_TBXQSLICESHADOW): `scan_unit` with the whole group (on its own table, into
    // its own list), then the scan in slices, compared (module doc).
    /// The scan of one unit in slices, compared with `scan_unit` in shadow mode.
    #[inline(always)]
    fn scan_unit_sliced_or_shadow<
        S: XShadowList + XUnitSort + Clone,
        const B: bool,
        const PF: bool,
    >(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        let group = self.x_sliced(unit);
        // EXPERIMENT (LOSAT_X_TBXHITBF): the scan in slices or its branch-free hit loop
        // (x_tbx_hitbf.rs).
        if !shadow() {
            self.scan_unit_sliced_or_bf::<S, B, PF>(group, unit, inp, scratch, init_hsps);
            return;
        }
        let mut reference = S::default();
        self.scan_unit::<S, XShadowRef, PF>(unit, inp, scratch, &mut reference, &mut XShadowRef);
        let reference_diag = std::mem::take(&mut scratch.diag);
        let scans = self.scan_unit_sliced_or_bf::<S, B, PF>(group, unit, inp, scratch, init_hsps);
        let mut sorted = init_hsps.clone();
        sorted.x_sort();
        reference.x_sort();
        if sorted != reference || scratch.diag != reference_diag {
            let first = sorted.first_difference(&reference);
            let cell = scratch
                .diag
                .iter()
                .zip(&reference_diag)
                .position(|(a, b)| a != b);
            panic!(
                "LOSAT_X_TBXQSLICESHADOW: subject {} frame {} chunk {} group {} ({} slices): sorted records {} vs {} (first difference at {}), tables {} vs {} cells (first difference at {:?})",
                inp.s_idx,
                unit.f,
                unit.k,
                unit.g,
                group.slices.len(),
                sorted.len_of(),
                reference.len_of(),
                first,
                scratch.diag.len(),
                reference_diag.len(),
                cell
            );
        }
        QS_UNITS.fetch_add(1, AtomicOrdering::Relaxed);
        QS_SCANS.fetch_add(scans, AtomicOrdering::Relaxed);
        QS_SHADOW_UNITS.fetch_add(1, AtomicOrdering::Relaxed);
        QS_SHADOW_RECORDS.fetch_add(sorted.len_of() as u64, AtomicOrdering::Relaxed);
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
    // ```
    // The scan loop of `scan_unit` (RI:3082-3168 on one unit) once per slice with the slice's table,
    // in slice order, on the unit's one table (module doc, items 1-3), with `x_hit`, or with
    // `x_hit_batched` of LOSAT_X_TBXEXTBATCH (`B`; the batch flushed at the end of the unit), and
    // LOSAT_X_TBXPREFETCH's hint (`PF`) as in `scan_unit`.
    /// Scans one unit of the cut batch's group `group` in slices and saves its init HSPs into
    /// `init_hsps`. Returns the number of slices scanned.
    #[inline(always)]
    pub(super) fn scan_unit_sliced<S: XInitSink, const B: bool, const PF: bool>(
        &self,
        group: &XSliced,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) -> u64 {
        self.scan_unit_sliced_hits::<S, B, PF, false>(group, unit, inp, scratch, init_hsps)
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF, x_tbx_hitbf.rs): with `BF` the hits of a batch go through
    // `x_hits_bf`, which takes the too-far and overlap outcomes without a branch and hands the
    // others to `x_hit` (or `x_hit_batched` with `B`); without it this is the loop before the
    // switch. `BF` runs only where `x_tbx_hitbf::fits` holds.
    /// `scan_unit_sliced`, with `BF` the branch-free hit loop of LOSAT_X_TBXHITBF.
    #[inline(always)]
    pub(super) fn scan_unit_sliced_hits<
        S: XInitSink,
        const B: bool,
        const PF: bool,
        const BF: bool,
    >(
        &self,
        group: &XSliced,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) -> u64 {
        let wordsize = inp.wordsize;
        let window = inp.window;
        let s_frame = &inp.s_frames_preliminary[unit.f];
        let subject_full = &s_frame.aa_seq;
        let subject_all = &subject_full[1..subject_full.len() - 1];
        let chunk = unit.chunk;
        let chunk_end = chunk.offset.saturating_add(chunk.length);
        let subject = &subject_all[chunk.offset..chunk_end];
        if subject.len() < wordsize as usize {
            return 0;
        }

        // The unit's table, as `scan_unit` sets it up from the whole group's bounds (module doc,
        // item 3); LOSAT_X_SEEDUNMAP as there.
        let span = group.span();
        let width = span + chunk.length;
        let lo = group.q_min.wrapping_sub(chunk.length as u32 - 1);
        if scratch.diag.capacity() < width && x_unit_mapped(scratch.unmap, span, chunk.length) {
            scratch.diag = Vec::with_capacity(x_mapped_len::<DiagStruct>(width));
            SU_DIAG.fetch_add(1, AtomicOrdering::Relaxed);
        }
        scratch.diag.clear();
        scratch.diag.resize(width, DiagStruct::default());

        let h = XHitCtx {
            diag_offset: window,
            window,
            wordsize,
            // `get_context_idx` reads only the batch's context bounds (`frame_bases`,
            // `num_contexts`), which every slice table carries.
            lookup: &group.slices[0].lookup,
            contexts: inp.contexts,
            cutoff_scores: inp.cutoff_scores,
            x_dropoff_per_context: inp.x_dropoff_per_context,
            subject,
            s_f_idx: unit.f,
            s_frame: s_frame.frame,
            s_idx: inp.s_idx,
            s_len: inp.s_len,
        };
        let mut b = XPendBatch::new(BATCH);
        let mut scans = 0u64;

        let base_seq_ranges: [(i32, i32); 1] = [(0, chunk.length as i32)];
        let scan_interiors = tblastx_scan_interiors(chunk.length, None);
        for slice in &group.slices {
            let sub = &slice.lookup;
            if sub.longest_chain <= 0 {
                continue;
            }
            scans += 1;
            // NCBI sizes the offset array from the table it scans (lookup_wrap.c:264-266).
            let offset_array_size: i32 = OFFSET_ARRAY_SIZE + sub.longest_chain;
            // EXPERIMENT (LOSAT_X_TBXHITBF): with `BF` and `PF`, X_PREFETCH_AHEAD entries more, as
            // in `scan_unit_hits`.
            let offset_pairs_len =
                crate::algorithm::tblastx::x_scan_bf::x_offset_pairs_len(offset_array_size)
                    + if BF && PF { X_PREFETCH_AHEAD } else { 0 };
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
                // [C] scan_range[0] = 0;
                // [C] scan_range[1] = subject->seq_ranges[0].left;
                // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
                let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];
                // [C] while (scan_range[1] <= scan_range[2])
                while scan_range[1] <= scan_range[2] {
                    let prev_scan_left = scan_range[1];
                    // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
                    let hits = crate::algorithm::tblastx::x_scan_bf::x_scan_subject(
                        slice.scan_bf,
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
                    // EXPERIMENT (LOSAT_X_TBXHITBF): the loop below, with the too-far and overlap
                    // outcomes taken without a branch (x_tbx_hitbf.rs).
                    if BF {
                        // SAFETY: hits <= offset_array_size, and with `PF` the X_PREFETCH_AHEAD
                        // entries after them are in offset_pairs too (sized above); every cell is
                        // below width = scratch.diag.len() as below; `fits` holds (the caller's
                        // rule); every candidate's cell is a cell of this unit's table.
                        unsafe {
                            x_tbx_hitbf::x_hits_bf::<PF>(
                                diag_ptr,
                                offset_pairs_ptr,
                                hits as usize,
                                lo,
                                width,
                                window,
                                wordsize,
                                |cell, query_offset, subject_offset| {
                                    if B {
                                        x_tbx_extbatch::x_hit_batched(
                                            diag_ptr,
                                            cell,
                                            query_offset,
                                            subject_offset,
                                            &h,
                                            &mut b,
                                            init_hsps,
                                        );
                                    } else {
                                        x_hit(
                                            &mut *diag_ptr.add(cell),
                                            query_offset,
                                            subject_offset,
                                            &h,
                                            init_hsps,
                                            &mut XNoProbe,
                                        );
                                    }
                                },
                            )
                        };
                        continue;
                    }
                    // [C] for (i = 0; i < hits; ++i)
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
                        // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
                        // Unit cell: the raw diagonal minus the unit's lowest one (no mask).
                        let diag_coord =
                            query_offset.wrapping_sub(subject_offset).wrapping_sub(lo) as usize;
                        debug_assert!(diag_coord < width, "unit diagonal outside its table");
                        if B {
                            // SAFETY: the slice's entries lie between the group's q_min and q_max
                            // and subject_offset < s_len, so diag_coord < width =
                            // scratch.diag.len(); every candidate's cell is such a cell.
                            unsafe {
                                x_tbx_extbatch::x_hit_batched(
                                    diag_ptr,
                                    diag_coord,
                                    query_offset,
                                    subject_offset,
                                    &h,
                                    &mut b,
                                    init_hsps,
                                )
                            };
                        } else {
                            // SAFETY: as above.
                            let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };
                            x_hit(
                                diag_entry,
                                query_offset,
                                subject_offset,
                                &h,
                                init_hsps,
                                &mut XNoProbe,
                            );
                        }
                    }
                }
            }
        }
        if B && b.pending() > 0 {
            // SAFETY: the candidates' cells are cells of this unit's table.
            unsafe { x_tbx_extbatch::flush(&mut b, scratch.diag.as_mut_ptr(), &h, init_hsps) };
        }
        scans
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

    fn aa(s: &str) -> Vec<u8> {
        s.bytes()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect()
    }

    // No NCBI counterpart: the rule's arithmetic (module doc).
    #[test]
    fn the_partition_cuts_only_groups_wider_than_a_slice() {
        assert_eq!(x_partition(&[(1, 101), (102, 150)], 100, 0), None);
        let parts = x_partition(&[(1, 102), (103, 150)], 100, 0).expect("a cut");
        assert_eq!(parts, vec![vec![1, 51], vec![103]]);
        // Equal slices, the last ending at the group's end.
        let parts = x_partition(&[(10, 1010)], 300, 0).expect("a cut");
        assert_eq!(parts, vec![vec![10, 260, 510, 760]]);
        // At most X_QSLICE_MAX_SLICES slices; strictly ascending starts.
        let parts = x_partition(&[(0, 1_000_000)], 2, 0).expect("a cut");
        assert_eq!(parts[0].len() as u64, X_QSLICE_MAX_SLICES);
        assert!(parts[0].windows(2).all(|w| w[0] < w[1]));
        // A group shorter than the cap of slices: one letter each.
        let parts = x_partition(&[(5, 9)], 1, 0).expect("a cut");
        assert_eq!(parts, vec![vec![5, 6, 7, 8]]);
        // An empty group stays one slice.
        assert_eq!(x_partition(&[(7, 7)], 2, 0), None);
        // The built-in rule's bound: only the groups of at least `min_len` letters are cut.
        assert_eq!(x_partition(&[(0, 1000), (1000, 1500)], 100, 1001), None);
        let parts = x_partition(&[(0, 1000), (1000, 1500)], 100, 600).expect("a cut");
        assert_eq!(parts[0].len(), 10);
        assert_eq!(parts[1], vec![1000]);
        assert_eq!(
            X_QSLICE_MIN_SPAN_THREADS.div_ceil(8),
            768 << 10,
            "8 threads: groups of 768 Ki letters or more"
        );
    }

    // No NCBI counterpart: the cut of the lookup locations (module doc, item 5): every word start of
    // every location goes to exactly one piece, in the slice that holds it, in list order.
    #[test]
    fn slice_locations_keep_every_word_once_in_its_slice() {
        let mut rng = Rng(0x51ce_0000_1234_abcd);
        for _ in 0..2000 {
            let wl = 3i32;
            let a = rng.below(50) as i32;
            let b = a + 1 + rng.below(400) as i32;
            let k = 1 + rng.below(8) as usize;
            let mut parts: Vec<i32> = vec![a];
            let mut cuts: Vec<i32> = (1..k)
                .map(|_| a + 1 + rng.below((b - a) as u64) as i32)
                .collect();
            cuts.sort_unstable();
            cuts.dedup();
            cuts.retain(|&c| c > a && c < b);
            parts.extend(cuts);
            let mut locations = Vec::new();
            let mut at = a;
            while at < b {
                let from = at + rng.below(4) as i32;
                let to = (from + rng.below(30) as i32).min(b - 1);
                if from <= to {
                    locations.push((from, to));
                }
                at = to + 1 + rng.below(3) as i32;
            }
            let pieces = x_slice_locations(&locations, &parts, wl);
            let starts = |list: &[(i32, i32)]| -> Vec<i32> {
                list.iter()
                    .flat_map(|&(from, to)| from..=(to - wl + 1))
                    .collect()
            };
            let mut all = Vec::new();
            for (j, list) in pieces.iter().enumerate() {
                let lo = parts[j];
                let hi = parts.get(j + 1).copied().unwrap_or(i32::MAX);
                for &(from, to) in list {
                    assert!(to - from + 1 >= wl, "a piece holds a word");
                    assert!(
                        locations.iter().any(|&(f, t)| f <= from && to <= t),
                        "a piece lies in a location"
                    );
                    for s in from..=(to - wl + 1) {
                        assert!(
                            j == 0 && s < lo || (lo <= s && s < hi),
                            "word {s} in slice {j}"
                        );
                    }
                }
                let mut ordered = starts(list);
                let sorted = ordered.windows(2).all(|w| w[0] < w[1]);
                assert!(sorted, "slice {j}: list order is offset order");
                all.append(&mut ordered);
            }
            all.sort_unstable();
            assert_eq!(all, starts(&locations));
        }
    }

    // No NCBI counterpart: one instance of the scan in slices by the bits of `variant`: batched
    // extensions (2), the prefetch (1), the branch-free hit loop of LOSAT_X_TBXHITBF (4).
    fn scan_sliced_variant<S: XInitSink>(
        pp: &XPairPar,
        variant: usize,
        sliced: &XSliced,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        got: &mut S,
    ) -> u64 {
        match variant {
            0 => {
                pp.scan_unit_sliced_hits::<_, false, false, false>(sliced, unit, inp, scratch, got)
            }
            1 => pp.scan_unit_sliced_hits::<_, false, true, false>(sliced, unit, inp, scratch, got),
            2 => pp.scan_unit_sliced_hits::<_, true, false, false>(sliced, unit, inp, scratch, got),
            3 => pp.scan_unit_sliced_hits::<_, true, true, false>(sliced, unit, inp, scratch, got),
            4 => pp.scan_unit_sliced_hits::<_, false, false, true>(sliced, unit, inp, scratch, got),
            5 => pp.scan_unit_sliced_hits::<_, false, true, true>(sliced, unit, inp, scratch, got),
            6 => pp.scan_unit_sliced_hits::<_, true, false, true>(sliced, unit, inp, scratch, got),
            _ => pp.scan_unit_sliced_hits::<_, true, true, true>(sliced, unit, inp, scratch, got),
        }
    }

    // No NCBI counterpart: compares the scan in slices with `scan_unit` on the same units (module
    // doc, items 1-4), for slices built directly and filtered (item 5), with and without the
    // batched extensions, the prefetch and the branch-free hit loop (LOSAT_X_TBXHITBF).
    #[test]
    fn slices_equal_scan_unit() {
        let alphabet = aa("ACDEKLMNPQRSTVWY");
        let narrow = aa("AKLE");
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let serial = crate::utils::threading::SearchPool::x_serial();
        let mut rng = Rng(0x5eed_51ce_0123_4567);
        let mut records = 0usize;
        let mut cut_groups = 0usize;
        for case in 0..12 {
            let abc = if case % 3 == 2 { &narrow } else { &alphabet };
            let q_len = 200 + rng.below(500) as usize;
            let frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
                .iter()
                .map(|&f| frame(f, &residues(&mut rng, q_len, abc, &[])))
                .collect();
            let source: Vec<u8> = frames
                .iter()
                .flat_map(|f| f.aa_seq[1..f.aa_len + 1].to_vec())
                .collect();
            let s_len = 300 + rng.below(600) as usize;
            let s_frames: Vec<QueryFrame> = if case % 4 == 3 {
                frames.clone()
            } else {
                [1i8, 2, 3, -1, -2, -3]
                    .iter()
                    .map(|&f| frame(f, &residues(&mut rng, s_len, abc, &source)))
                    .collect()
            };
            let queries = [frames.clone()];
            let threshold = 11;
            let (lookup, contexts) = build_ncbi_lookup(&queries, threshold, &params, true);
            let starts = group_starts(&contexts);
            // Every fourth case cuts no group (one slice each, as a group of a cut batch that is
            // not cut).
            let spans = x_group_spans(&contexts, &starts);
            let width = [7usize, 60, 150, 1000][case % 4];
            let parts = x_partition(&spans, width, 0)
                .unwrap_or_else(|| spans.iter().map(|&(a, _)| vec![a]).collect());
            let filtered = build_sliced(&lookup, &parts);
            let q = XPairLutQuery::x_prepare(&queries, threshold, true).expect("direct builder");
            let direct = lut_direct_sliced(&q, &starts, &parts, &serial, false).expect("groups");
            lut_shadow_compare(Some(&direct), Some(&filtered));
            cut_groups += direct.iter().filter(|g| g.slices.len() > 1).count();
            let groups = build_groups(&lookup, &contexts);
            // The stub's longest_chain over the slices is the full table's.
            let (stub, _) = q.x_into_stub(tables(&direct));
            assert_eq!(stub.longest_chain, lookup.longest_chain);
            for (g, s) in groups.iter().zip(&direct) {
                assert_eq!((g.q_min, g.q_max), (s.q_min, s.q_max), "group bounds");
                assert_eq!(g.lookup.longest_chain > 0, s.has_entries());
            }
            let cutoff_scores: Vec<i32> = (0..contexts.len())
                .map(|i| 10 + (i as i32 % 4) * 3)
                .collect();
            let x_dropoff_per_context: Vec<i32> = (0..contexts.len())
                .map(|i| [7, 16, 30, 200][i % 4])
                .collect();
            let pp = XPairPar {
                groups,
                shadow: false,
                stub: false,
                sliced: Some(direct),
            };
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
                let offset = if case % 2 == 1 { s_frame.aa_len / 3 } else { 0 };
                let chunk = SubjectChunk {
                    offset,
                    length: s_frame.aa_len - offset,
                    overlap: 0,
                };
                for g in 0..pp.groups.len() {
                    if pp.groups[g].lookup.longest_chain <= 0 {
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
                    reference.x_sort();
                    for sliced in [&pp.sliced.as_ref().unwrap()[g], &filtered[g]] {
                        let mut plain: Vec<Vec<InitHSP>> = Vec::new();
                        for variant in 0..8 {
                            let mut got: Vec<InitHSP> = Vec::new();
                            scan_sliced_variant(
                                &pp,
                                variant,
                                sliced,
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut got,
                            );
                            assert!(
                                scratch.diag == reference_diag,
                                "case {case} frame {f} group {g} variant {variant}: table"
                            );
                            if variant < 4 {
                                plain.push(got.clone());
                            } else {
                                // LOSAT_X_TBXHITBF: the same records in the same order.
                                assert_eq!(
                                    got,
                                    plain[variant - 4],
                                    "case {case} frame {f} group {g} variant {variant}: order"
                                );
                            }
                            got.x_sort();
                            assert_eq!(
                                got, reference,
                                "case {case} frame {f} group {g} variant {variant}: records"
                            );
                        }
                        // LOSAT_X_TBXHITBFSHADOW: both scans compared, the branch-free one's
                        // records kept.
                        let mut shadowed: Vec<InitHSP> = Vec::new();
                        pp.scan_unit_sliced_mode::<_, true, false>(
                            2,
                            sliced,
                            &unit,
                            &inp,
                            &mut scratch,
                            &mut shadowed,
                        );
                        assert_eq!(
                            shadowed, plain[2],
                            "case {case} frame {f} group {g}: shadow"
                        );
                    }
                    records += reference.len();
                    // The compact records, sorted, are today's too.
                    let mut compact_ref: Vec<XUnitHsp> = Vec::new();
                    pp.scan_unit::<_, _, false>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut compact_ref,
                        &mut XNoProbe,
                    );
                    let mut compact: Vec<XUnitHsp> = Vec::new();
                    pp.scan_unit_sliced::<_, true, false>(
                        &pp.sliced.as_ref().unwrap()[g],
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut compact,
                    );
                    let mut compact_bf: Vec<XUnitHsp> = Vec::new();
                    pp.scan_unit_sliced_hits::<_, true, false, true>(
                        &pp.sliced.as_ref().unwrap()[g],
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut compact_bf,
                    );
                    assert_eq!(
                        compact_bf, compact,
                        "case {case} frame {f} group {g}: compact, branch-free"
                    );
                    compact.x_sort();
                    compact_ref.x_sort();
                    assert_eq!(
                        compact, compact_ref,
                        "case {case} frame {f} group {g}: compact"
                    );
                }
            }
        }
        assert!(cut_groups > 40, "cut groups {cut_groups}");
        assert!(records > 500, "records {records}");
    }
}
