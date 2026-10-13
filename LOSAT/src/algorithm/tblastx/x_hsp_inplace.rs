//! EXPERIMENT (LOSAT_X_HSPINPLACE / LOSAT_X_HSPINPLACESHADOW): the TBLASTX HSP list is sorted,
//! cut into linking groups, put together again and replayed into NCBI's output order without a
//! second whole copy of the list alive. Only where the records lie changes; every comparison,
//! every per-group linking call and every value written is the reference's.
//!
//! This file is a child module of `sum_stats_linking/linking.rs` (declared there with `#[path]`,
//! as `lookup/x_lut_direct.rs` is for `backbone.rs`), so it uses that file's comparators,
//! per-group functions and replay without changing them. `ungapped_hit_sort.rs` reaches the sort
//! through `sum_stats_linking::x_hsp_inplace`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1383
//! ```c
//! void Blast_HSPListSortByScore(BlastHSPList* hsp_list)
//! {
//!     if (!hsp_list || hsp_list->hspcnt <= 1)
//!         return;
//!
//!     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
//!         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//!               ScoreCompareHSPs);
//!     }
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:476-486
//! ```c
//! link_hsp_array =
//!    (LinkHSPStruct**) malloc(total_number_of_hsps*sizeof(LinkHSPStruct*));
//! for (index = 0; index < total_number_of_hsps; ++index) {
//!    link_hsp_array[index] = (LinkHSPStruct*) calloc(1, sizeof(LinkHSPStruct));
//!    link_hsp_array[index]->hsp = hsp_array[index];
//! }
//!
//! /* Sort by (reverse) position. */
//! if (kTranslatedQuery) {
//!    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//!          s_RevCompareHSPsTbx);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:513-534
//! ```c
//! Int4 strand_factor = (kTranslatedQuery ? 3 : 1);
//! for (index=0;index<number_of_hsps;index++)
//! {
//!   H=link_hsp_array[index];
//! ...
//!   if (H->prev != NULL &&
//!       ((H->hsp->context/strand_factor) !=
//!        (H->prev->hsp->context/strand_factor) ||
//!        (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
//!   { /* If frame switches, then start new list. */
//!      hp_frame_number[cur_frame]--;
//!      hp_frame_number[++cur_frame]++;
//!      hp_frame_start[cur_frame] = H;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
//! ```c
//! for (frame_index=0; frame_index<num_query_frames; frame_index++)
//! {
//!    hp_start = s_LinkHSPStructReset(hp_start);
//!    hp_start->next = hp_frame_start[frame_index];
//!    hp_frame_start[frame_index]->prev = hp_start;
//!    number_of_hsps = hp_frame_number[frame_index];
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994,1012-1020,1044-1061,1080-1082
//! ```c
//! if (kTranslatedQuery) {
//!    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//!          s_RevCompareHSPsTransl);
//!    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//!          s_FwdCompareHSPsTransl);
//! ...
//! /* hook up the HSP's. */
//! first_hsp = NULL;
//! for (index=0, last_hsp=NULL;index<total_number_of_hsps; index++)
//! {
//!    H = link_hsp_array[index];
//!
//!    /* If this is not a single piece or the start of a chain, then Skip it. */
//!    if (H->linked_set == TRUE && H->start_of_chain == FALSE)
//!       continue;
//! ...
//!    /* The first one has the number of links correct. */
//!    num_links = H->hsp_link.num[ordering_method];
//!    link = H->hsp_link.link[ordering_method];
//!    while (link)
//!    {
//!       H->hsp->num = num_links;
//! ...
//!    /* Set these for last link in chain. */
//!    H->hsp->num = num_links;
//! ...
//! for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
//!    hsp_list->hsp_array[index] = H->hsp;
//! ```
//! NCBI sorts arrays of pointers in place, runs its frame groups as sublists of the one sorted
//! array and writes the final order back into the same `hsp_array`: it never holds a second copy
//! of the HSPs. LOSAT's reference path owns complete `UngappedHit` records (152 bytes) and copies
//! the whole list at each of these steps (`hits.to_vec()` in every sort, groups grown by `push`
//! while the sorted list is drained, a result grown by `extend`, a cloned replay). This path keeps
//! one list.
//!
//! # Why the result is the reference's (each point checked against the code it replaces)
//!
//! 1. Sort (`sort_in_place`; reference `sort_ungapped_hits_by_index_replay`,
//!    `ungapped_hit_sort.rs`). The reference sorts the indices `0..n` with `slice::sort_by` and a
//!    comparator that reads only the two records, taken from an untouched copy of the list. This
//!    path makes the same `sort_by` call on the same initial indices and reads the records from
//!    the list itself, which is not written during the sort: the same records at the same
//!    indices, so the same comparison results in the same sequence and the same index order.
//!    The reference then writes `hits[dst] = original[indices[dst]]`; `apply_permutation_in_place`
//!    produces exactly that arrangement by walking each cycle of the permutation with swaps
//!    (after the walk of a cycle every position `dst` of it holds the record that was at
//!    `indices[dst]`, and the cycle's first record has travelled to the position that asks for
//!    it). No record is cloned, so no field can change.
//! 2. Groups (`groups_off_the_end`; reference `linking.rs`, the `current_group.push(hit)` loop).
//!    The reference starts a new group at a record whose `translated_context_group` or sign of
//!    `s_frame` differs from the record before it, so its groups are the maximal runs of
//!    consecutive records of the sorted list, in list order. This path finds the same run starts
//!    with the same two tests on the same neighbours, cuts the runs off the end of the list
//!    (`Vec::split_off` copies the tail records, in order, into a Vec of exactly that size), keeps
//!    the remaining head as the first group, and reverses the group list: the same groups with
//!    the same records in the same order. Only the capacity of the Vecs differs, which neither
//!    kernel reads (they index, read and write the records and return the same Vec).
//! 3. Per-group linking. The same functions (`link_hsp_group_incremental`, or
//!    `link_hsp_group_ncbi` when the literal port is asked for) with the same arguments, serially
//!    in group order or on the pool exactly as in the reference (`into_par_iter().enumerate()`,
//!    indexed collect, sort by index). The scratch buffers handed to the kernels start empty
//!    unless the literal port can run (`LOSAT_LINKING_LEGACY`, a linking trace, or
//!    `LOSAT_LINK_FAST_SHADOW`): the incremental kernel never touches them, and the literal port
//!    clears and reserves them itself for every group (`pool_hsp_links.clear(); .reserve(n)`,
//!    `pool_lh_helpers.clear(); .resize(..)`), so their initial capacity changes no value.
//! 4. Assembly (`assemble`). The reference concatenates the linked groups in group order
//!    (`ordered.extend(processed)`, or `flat_map` over the index-sorted pool results). This path
//!    extends one Vec of exact capacity with the groups in the same order.
//! 5. Output order (`replay_in_place`; reference: two `sort_linked_hsp_references` calls and
//!    `replay_ncbi_output_list`). The reference stable-sorts a Vec of references that starts in
//!    list order, twice; this path stable-sorts the indices `0..n` with the same two comparators
//!    applied to the same records, so `order[i]` is the list position of the reference's
//!    `sorted_hits[i]`. The walk is the reference's statement for statement with
//!    `results[order[i]]` in place of `sorted_hits[i]`: the same `id_to_index`, the same skip of
//!    chain members, the same `num` writes, the same stops, hence the same `replayed_indices`.
//!    The reference returns clones of `sorted_hits[replayed_indices[k]]`, taken after all `num`
//!    writes. When `replayed_indices` is a permutation of `0..n` (every record exactly once,
//!    which is what linked chains give), moving the records by that permutation gives the same
//!    list; otherwise (a record left out or reached twice) this path returns the reference's own
//!    expression, the clones. Both outcomes are counted.
//!
//! The reference functions keep their bodies; each has one dispatch at its head. A linking call
//! with `LOSAT_DEBUG_CUTOFFS` or `LOSAT_DIAGNOSTICS` set runs the reference body (its debug output
//! lives there).
//!
//! `LOSAT_X_HSPINPLACESHADOW=1` runs every sort call and every linking call both ways (the
//! reference on a clone of the input), compares the two results record by record, every field,
//! `e_value` by bits, stops at the first difference and continues with the reference's result;
//! the summary is printed by `main`.
//!
//! EXPERIMENT (LOSAT_X_LINKSORT / LOSAT_X_LINKSORTSHADOW, on by default as a strict-set switch):
//! the two index sort sites of this file (`sort_in_place`, `replay_in_place_counted`) then call
//! `x_link_sort.rs`, which returns the same index order from one final sort instead of two and
//! from a stable parallel sort for long lists (its header has the argument). With the switch off
//! the statements described in items 1 and 5 run.

use super::*;
// Diagnostic (LOSAT_TIMING): the clock of a linking call (tblastx/x_link_timing.rs). It reads
// the clock only in a call the batch function armed, and the linking reads nothing from it.
use crate::algorithm::tblastx::x_link_timing::{KernelTimes, LinkLap, LinkStage};
use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

// No NCBI counterpart: counters for LOSAT_X_STATS and the shadow summary; they do not change any
// value NCBI computes.
static SORTS: AtomicU64 = AtomicU64::new(0);
static LINKINGS: AtomicU64 = AtomicU64::new(0);
static HITS: AtomicU64 = AtomicU64::new(0);
static GROUPS: AtomicU64 = AtomicU64::new(0);
static REPLAY_INPLACE: AtomicU64 = AtomicU64::new(0);
static REPLAY_FALLBACKS: AtomicU64 = AtomicU64::new(0);
static GATE_REFERENCE: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: how the shadow reaches a reference body that has the dispatch at its head.
// The shadow sets the flag of the calling thread and calls the reference function; the dispatch at
// the head of that call takes the flag and lets the body run. Nothing between the two can hand the
// thread to other work, and the two reference functions have separate flags, so the sort inside a
// reference linking call is shadowed as well. The unit tests reach the reference bodies the same
// way (`tests::reference_sort`, `tests::reference_linking`), with the switch on as well; nothing
// else sets a flag, so a search with the switch on only reads them.
thread_local! {
    static SORT_REFERENCE_CALL: Cell<bool> = const { Cell::new(false) };
    static LINKING_REFERENCE_CALL: Cell<bool> = const { Cell::new(false) };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
// ```c
// /* Sort by (reverse) position. */
// if (kTranslatedQuery) {
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ```
// Dispatch switch: 0 = the reference (copying) path, 1 = this path, 2 = both, compared.
/// 0 = off, 1 = LOSAT_X_HSPINPLACE, 2 = LOSAT_X_HSPINPLACESHADOW.
pub(crate) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_HSPINPLACESHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_HSPINPLACE") {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: the reference linking body prints the LOSAT_DEBUG_CUTOFFS block and passes
// the LOSAT_DIAGNOSTICS flag to the literal port; with either set the reference body runs.
fn gate() -> bool {
    static GATE: OnceLock<bool> = OnceLock::new();
    *GATE.get_or_init(|| std::env::var_os("LOSAT_DEBUG_CUTOFFS").is_some() || diagnostics_enabled())
}

/// The summary lines (printed by `main` at exit): the shadow summary and the LOSAT_X_STATS
/// counters.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    let m = mode();
    if m == 0 {
        return;
    }
    let load = |c: &AtomicU64| c.load(AtomicOrdering::Relaxed);
    if m == 2 {
        eprintln!(
            "[X_SHADOW] HSPINPLACE sorts={} linkings={} hits_compared={} replay_fallbacks={}",
            load(&SORTS),
            load(&LINKINGS),
            load(&HITS),
            load(&REPLAY_FALLBACKS)
        );
    }
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] HSPINPLACE sorts={} linkings={} hits={} groups={} replay_inplace={} replay_fallbacks={} gate_reference={}",
            load(&SORTS),
            load(&LINKINGS),
            load(&HITS),
            load(&GROUPS),
            load(&REPLAY_INPLACE),
            load(&REPLAY_FALLBACKS),
            load(&GATE_REFERENCE)
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
// ```c
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
// ```
// NCBI's qsort moves the pointers inside the one array. `source[dst]` names the record that
// belongs at `dst`; each cycle of that permutation is walked once with swaps (module doc, item 1).
/// Rearranges `items` so that `items[dst]` is the element that was at `source[dst]`. `source` must
/// be a permutation of `0..items.len()`; it is left as the identity.
fn apply_permutation_in_place<T>(items: &mut [T], source: &mut [usize]) {
    assert_eq!(
        items.len(),
        source.len(),
        "LOSAT_X_HSPINPLACE: permutation length"
    );
    for start in 0..source.len() {
        let mut dst = start;
        loop {
            let src = source[dst];
            // `dst` is settled by the swap below (or needs none): mark it.
            source[dst] = dst;
            if src == start {
                // The cycle is closed: the record that started at `start` has been carried here.
                break;
            }
            items.swap(dst, src);
            dst = src;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
// ```c
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
// ```c
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ```
// The index sort of `sort_ungapped_hits_by_index_replay` (same `sort_by`, same comparator, same
// records at the same indices), then the records moved instead of copied and cloned back.
// `site` names the caller for LOSAT_X_LINKSORT and is not read with that switch off.
/// The result of `sort_ungapped_hits_by_index_replay` without a copy of the list.
fn sort_in_place(
    hits: &mut [UngappedHit],
    compare: fn(&UngappedHit, &UngappedHit) -> Ordering,
    site: x_link_sort::Site,
) {
    if hits.len() <= 1 {
        return;
    }
    let mut indices: Vec<usize> = (0..hits.len()).collect();
    {
        let original: &[UngappedHit] = hits;
        // EXPERIMENT (LOSAT_X_LINKSORT / LOSAT_X_LINKSORTSHADOW), dispatch: the same index order
        // from a stable parallel sort when the list is long and the call may use the pool
        // (tblastx/x_link_sort.rs, item 2 of its header). With the switch off the statement
        // below runs.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
        // ```c
        //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
        //          s_RevCompareHSPsTbx);
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
        // ```c
        //         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
        //               ScoreCompareHSPs);
        // ```
        if x_link_sort::mode() != 0 {
            x_link_sort::sort_indices(&mut indices, original, compare, site);
        } else {
            indices.sort_by(|&lhs_index, &rhs_index| {
                compare(&original[lhs_index], &original[rhs_index])
            });
        }
    }
    apply_permutation_in_place(hits, &mut indices);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1082
// ```c
// for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
//    hsp_list->hsp_array[index] = H->hsp;
// ```
// No computation: the shadow compares the two lists, every field of every HSP, and stops at the
// first difference. A new field of `UngappedHit` does not compile here until it is compared.
/// Shadow check: `actual` (this path) against `reference`, record by record.
fn shadow_compare(what: &str, actual: &[UngappedHit], reference: &[UngappedHit]) {
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
    let n = actual.len().max(reference.len());
    for k in 0..n {
        match (actual.get(k), reference.get(k)) {
            (Some(a), Some(b)) if same(a, b) => {}
            (a, b) => panic!(
                "LOSAT_X_HSPINPLACESHADOW: {what}, index {k} differs (in place {} hits, reference {} hits): in place {a:?} vs reference {b:?}",
                actual.len(),
                reference.len()
            ),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
// ```c
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
// ```
// The dispatch of `sort_ungapped_hits_by_index_replay` (called after its `len <= 1` return, only
// when the switch or its shadow is on).
/// `true`: `hits` is sorted (in place, or by both ways and compared). `false`: the caller's
/// reference body sorts it.
pub(crate) fn sort_dispatch(
    hits: &mut [UngappedHit],
    compare: fn(&UngappedHit, &UngappedHit) -> Ordering,
) -> bool {
    let m = mode();
    if m == 0 {
        return false;
    }
    if SORT_REFERENCE_CALL.with(|flag| flag.replace(false)) {
        // A call for the reference body: the shadow's own (below) or a unit test's.
        return false;
    }
    match m {
        1 => {
            sort_in_place(hits, compare, x_link_sort::Site::Dispatch);
            SORTS.fetch_add(1, AtomicOrdering::Relaxed);
            HITS.fetch_add(hits.len() as u64, AtomicOrdering::Relaxed);
            true
        }
        _ => {
            let mut reference = hits.to_vec();
            SORT_REFERENCE_CALL.with(|flag| flag.set(true));
            crate::algorithm::tblastx::ungapped_hit_sort::sort_ungapped_hits_by_index_replay(
                &mut reference,
                compare,
            );
            assert!(
                !SORT_REFERENCE_CALL.with(|flag| flag.get()),
                "LOSAT_X_HSPINPLACESHADOW: the reference sort did not reach its dispatch"
            );
            sort_in_place(hits, compare, x_link_sort::Site::Dispatch);
            shadow_compare("sort", hits, &reference);
            SORTS.fetch_add(1, AtomicOrdering::Relaxed);
            HITS.fetch_add(reference.len() as u64, AtomicOrdering::Relaxed);
            // The search continues with the reference's result.
            hits.clone_from_slice(&reference);
            true
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:513-534
// ```c
// Int4 strand_factor = (kTranslatedQuery ? 3 : 1);
// for (index=0;index<number_of_hsps;index++)
// {
//   H=link_hsp_array[index];
// ...
//   if (H->prev != NULL &&
//       ((H->hsp->context/strand_factor) !=
//        (H->prev->hsp->context/strand_factor) ||
//        (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
//   { /* If frame switches, then start new list. */
//      hp_frame_number[cur_frame]--;
//      hp_frame_number[++cur_frame]++;
//      hp_frame_start[cur_frame] = H;
// ```
// NCBI's groups are sublists of the sorted array (`hp_frame_start`, `hp_frame_number`). The same
// run starts as the reference's push loop (the test of a record against the one before it), then
// the runs are cut off the end of the list so that the list and the groups together never hold
// more than the list plus the one group being cut (module doc, item 2).
/// The reference's groups (maximal runs of one `translated_context_group` and one sign of
/// `s_frame`, in list order) as exactly sized Vecs.
fn groups_off_the_end(mut hits: Vec<UngappedHit>) -> Vec<Vec<UngappedHit>> {
    if hits.is_empty() {
        return Vec::new();
    }
    let mut starts: Vec<usize> = Vec::new();
    for index in 1..hits.len() {
        let prev = &hits[index - 1];
        let hit = &hits[index];
        let ctx_group = translated_context_group(hit);
        let s_sign = hit.s_frame.signum();
        let prev_ctx_group = translated_context_group(prev);
        let prev_s_sign = prev.s_frame.signum();
        if prev_ctx_group != ctx_group || prev_s_sign != s_sign {
            // Frame switch - a new group starts here
            starts.push(index);
        }
    }
    // The list may arrive with spare capacity (it was grown by `extend`): give it back first.
    hits.shrink_to_fit();
    let mut frame_groups: Vec<Vec<UngappedHit>> = Vec::with_capacity(starts.len() + 1);
    for &start in starts.iter().rev() {
        // `start >= 1`: the tail records move, in order, into a Vec of exactly their number.
        let group = hits.split_off(start);
        hits.shrink_to(start);
        frame_groups.push(group);
    }
    // What is left is the first group.
    frame_groups.push(hits);
    frame_groups.reverse();
    frame_groups
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:448-454
// ```c
// hsp_array = hsp_list->hsp_array;
//
// lh_helper_size = MAX(1024,hsp_list->hspcnt+5);
// lh_helper = (LinkHelpStruct *)
//    calloc(lh_helper_size, sizeof(LinkHelpStruct));
// ```
// The scratch of the literal port `link_hsp_group_ncbi`, sized as in the reference when that port
// can run. The incremental kernel does not use it, so it then starts empty (module doc, item 3):
// the reference's capacity is 152 bytes per HSP of the whole list for every linking thread.
fn buffer_pools(initial_lh_size: usize, total_hits: usize, literal_port: bool) -> BufferPools {
    if literal_port {
        BufferPools {
            lh_helpers: Vec::with_capacity(initial_lh_size),
            hsp_links: Vec::with_capacity(total_hits),
        }
    } else {
        BufferPools {
            lh_helpers: Vec::new(),
            hsp_links: Vec::new(),
        }
    }
}

/// The values every per-group linking call of one linking call shares.
struct GroupArgs<'a> {
    use_legacy: bool,
    diag_enabled: bool,
    params: &'a KarlinParams,
    cutoffs: &'a LinkHspCutoffs,
    linking_params: &'a LinkingParams,
    query_contexts: &'a [QueryContext],
    subject_frame_bases: &'a [i32],
    length_adj_per_context: &'a [i64],
    eff_searchsp_per_context: &'a [i64],
    log_k_by_ctx: &'a [f64],
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// ```c
// for (frame_index=0; frame_index<num_query_frames; frame_index++)
// {
//    hp_start = s_LinkHSPStructReset(hp_start);
//    hp_start->next = hp_frame_start[frame_index];
//    hp_frame_start[frame_index]->prev = hp_start;
//    number_of_hsps = hp_frame_number[frame_index];
// ```
// One group through the reference's per-group function with the reference's arguments.
fn link_one_group(
    group_hits: Vec<UngappedHit>,
    args: &GroupArgs<'_>,
    pool: &mut BufferPools,
) -> Vec<UngappedHit> {
    if args.use_legacy {
        link_hsp_group_ncbi(
            group_hits,
            args.params,
            args.cutoffs,
            args.linking_params.gap_decay_rate,
            args.diag_enabled,
            args.linking_params.subject_len_nucl,
            args.query_contexts,
            args.subject_frame_bases,
            args.length_adj_per_context,
            args.eff_searchsp_per_context,
            args.log_k_by_ctx,
            &mut pool.lh_helpers,
            &mut pool.hsp_links,
        )
    } else {
        link_hsp_group_incremental(
            group_hits,
            args.params,
            args.cutoffs,
            args.linking_params.gap_decay_rate,
            args.linking_params.subject_len_nucl,
            args.query_contexts,
            args.subject_frame_bases,
            args.length_adj_per_context,
            args.eff_searchsp_per_context,
            args.log_k_by_ctx,
            &mut pool.lh_helpers,
            &mut pool.hsp_links,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1082
// ```c
// for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
//    hsp_list->hsp_array[index] = H->hsp;
// ```
// NCBI's groups stay in the one array. Here the linked groups go back into one list in group
// order, the reference's concatenation, with one allocation of the exact size; each group's
// memory is released as soon as its records are in the list (module doc, item 4).
/// The linked groups, in group order, as one list of capacity `total`.
fn assemble(groups: impl Iterator<Item = Vec<UngappedHit>>, total: usize) -> Vec<UngappedHit> {
    let mut ordered: Vec<UngappedHit> = Vec::with_capacity(total);
    for group in groups {
        ordered.extend(group);
    }
    ordered
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
// ```c
// if (kTranslatedQuery) {
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTransl);
//    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_FwdCompareHSPsTransl);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1012-1020
// ```c
// /* hook up the HSP's. */
// first_hsp = NULL;
// for (index=0, last_hsp=NULL;index<total_number_of_hsps; index++)
// {
//    H = link_hsp_array[index];
//
//    /* If this is not a single piece or the start of a chain, then Skip it. */
//    if (H->linked_set == TRUE && H->start_of_chain == FALSE)
//       continue;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1044-1061
// ```c
// /* The first one has the number of links correct. */
// num_links = H->hsp_link.num[ordering_method];
// link = H->hsp_link.link[ordering_method];
// while (link)
// {
//    H->hsp->num = num_links;
//    H->xsum = H->hsp_link.xsum[ordering_method];
//    H->next = (LinkHSPStruct*) link;
//    H->prev = last_hsp;
//    last_hsp = H;
//    H = H->next;
//    if (H != NULL)
//        link = H->hsp_link.link[ordering_method];
//    else
//        break;
// }
// /* Set these for last link in chain. */
// H->hsp->num = num_links;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1082
// ```c
// for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
//    hsp_list->hsp_array[index] = H->hsp;
// ```
// The two final sorts on indices instead of references, the walk of `replay_ncbi_output_list`
// with `results[order[i]]` for its `sorted_hits[i]`, and the records moved into the replayed
// order when it holds every record exactly once (module doc, item 5).
/// The result of the reference's two final sorts and `replay_ncbi_output_list` for a list of at
/// least two HSPs, without a second list unless the replay leaves a record out or repeats one.
fn replay_in_place(
    results: Vec<UngappedHit>,
    lap: &mut LinkLap,
    allow_parallel: bool,
) -> Vec<UngappedHit> {
    let (replayed, in_place) = replay_in_place_counted(results, lap, allow_parallel);
    if in_place {
        REPLAY_INPLACE.fetch_add(1, AtomicOrdering::Relaxed);
    } else {
        REPLAY_FALLBACKS.fetch_add(1, AtomicOrdering::Relaxed);
    }
    replayed
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1012-1020,1080-1082
// ```c
// /* hook up the HSP's. */
// first_hsp = NULL;
// for (index=0, last_hsp=NULL;index<total_number_of_hsps; index++)
// {
//    H = link_hsp_array[index];
//
//    /* If this is not a single piece or the start of a chain, then Skip it. */
//    if (H->linked_set == TRUE && H->start_of_chain == FALSE)
//       continue;
// ...
// for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
//    hsp_list->hsp_array[index] = H->hsp;
// ```
// The body of `replay_in_place`; the flag says which of the two outcomes produced the list.
// `allow_parallel` is the linking call's, for LOSAT_X_LINKSORT; it is not read with that switch
// off.
/// The replayed list and whether it was made in place (`true`) or by the reference's clones.
fn replay_in_place_counted(
    mut results: Vec<UngappedHit>,
    lap: &mut LinkLap,
    allow_parallel: bool,
) -> (Vec<UngappedHit>, bool) {
    let n = results.len();
    let mut order: Vec<usize> = (0..n).collect();
    // EXPERIMENT (LOSAT_X_LINKSORT / LOSAT_X_LINKSORTSHADOW), dispatch: the index order of the
    // two sorts below from one stable sort with `s_FwdCompareHSPsTransl`, ties by
    // `s_RevCompareHSPsTransl`, on the pool when the list is long and the call may use it
    // (tblastx/x_link_sort.rs, items 1 and 2 of its header). With the switch off the two sorts
    // below run.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
    // ```c
    // if (kTranslatedQuery) {
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTransl);
    //    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_FwdCompareHSPsTransl);
    // ```
    if x_link_sort::mode() != 0 {
        x_link_sort::final_sort(
            &mut order,
            &results,
            x_link_sort::Site::Final(allow_parallel),
        );
    } else {
        for compare in [
            rev_compare_hsps_transl as fn(&UngappedHit, &UngappedHit) -> Ordering,
            fwd_compare_hsps_transl,
        ] {
            let records: &[UngappedHit] = &results;
            // Stable, as `sort_linked_hsp_references`: comparator-equal records keep their order.
            order.sort_by(|&lhs, &rhs| compare(&records[lhs], &records[rhs]));
        }
    }
    lap.lap(LinkStage::FinalSorts);

    let max_link_id = results.iter().map(|h| h.link_id).max().unwrap_or(0);
    let mut id_to_index = vec![SENTINEL_IDX; max_link_id + 1];
    for (index, &position) in order.iter().enumerate() {
        let hit = &results[position];
        if hit.link_id < id_to_index.len() {
            id_to_index[hit.link_id] = index;
        }
    }

    let mut replayed_indices = Vec::with_capacity(n);
    for index in 0..n {
        let hit = &results[order[index]];
        if hit.linked_set && !hit.start_of_chain {
            continue;
        }

        let mut current = index;
        let mut hops = 0usize;
        // Every HSP of a chain gets the head's count; an HSP without a link keeps 1.
        let num_links = results[order[index]]
            .chain_next_link_id
            .map(|_| i32::from(results[order[index]].hsp_link_num));
        loop {
            replayed_indices.push(current);
            if let Some(num_links) = num_links {
                results[order[current]].num = num_links;
            }
            let Some(next_id) = results[order[current]].chain_next_link_id else {
                break;
            };
            if next_id >= id_to_index.len() {
                break;
            }
            let next_index = id_to_index[next_id];
            if next_index == SENTINEL_IDX {
                break;
            }
            current = next_index;
            hops += 1;
            if hops > n {
                break;
            }
        }
    }
    drop(id_to_index);

    // Is every sorted index replayed exactly once?
    let mut is_permutation = replayed_indices.len() == n;
    if is_permutation {
        let mut seen = vec![false; n];
        for &index in &replayed_indices {
            if std::mem::replace(&mut seen[index], true) {
                is_permutation = false;
                break;
            }
        }
    }

    if is_permutation {
        // Output position k takes `sorted_hits[replayed_indices[k]]`, the record at list position
        // `order[replayed_indices[k]]`.
        for slot in replayed_indices.iter_mut() {
            *slot = order[*slot];
        }
        drop(order);
        apply_permutation_in_place(&mut results, &mut replayed_indices);
        (results, true)
    } else {
        // A record left out or reached twice: the reference's own result, its clones.
        let replayed = replayed_indices
            .into_iter()
            .map(|index| results[order[index]].clone())
            .collect();
        (replayed, false)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:476-486
// ```c
// link_hsp_array =
//    (LinkHSPStruct**) malloc(total_number_of_hsps*sizeof(LinkHSPStruct*));
// for (index = 0; index < total_number_of_hsps; ++index) {
//    link_hsp_array[index] = (LinkHSPStruct*) calloc(1, sizeof(LinkHSPStruct));
//    link_hsp_array[index]->hsp = hsp_array[index];
// }
//
// /* Sort by (reverse) position. */
// if (kTranslatedQuery) {
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// ```c
// for (frame_index=0; frame_index<num_query_frames; frame_index++)
// {
//    hp_start = s_LinkHSPStructReset(hp_start);
//    hp_start->next = hp_frame_start[frame_index];
//    hp_frame_start[frame_index]->prev = hp_start;
//    number_of_hsps = hp_frame_number[frame_index];
// ```
// The steps of `apply_sum_stats_even_gap_linking_with_parallel` in its order, each with the
// reference's functions and values, on one list (module doc, items 1-5).
/// The result of `apply_sum_stats_even_gap_linking_with_parallel` without a second whole list.
#[allow(clippy::too_many_arguments)]
fn link_in_place(
    mut hits: Vec<UngappedHit>,
    params: &KarlinParams,
    linking_params: &LinkingParams,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    allow_parallel: bool,
) -> Vec<UngappedHit> {
    // No NCBI counterpart: the clock of this call for LOSAT_TIMING. The first local, so that it
    // is dropped last.
    let mut lap = LinkLap::start();
    if hits.is_empty() {
        return hits;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1773-1776
    // ```c
    //     /* Remove any information on number of linked HSPs from previous
    //        linking. */
    //     for (index = 0; index < hsp_list->hspcnt; ++index)
    //         hsp_list->hsp_array[index]->num = 1;
    // ```
    for (link_id, hit) in hits.iter_mut().enumerate() {
        hit.link_id = link_id;
        hit.chain_next_link_id = None;
        hit.num = 1;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:492-493
    // ```c
    // cutoff[0] = link_hsp_params->cutoff_small_gap;
    // cutoff[1] = link_hsp_params->cutoff_big_gap;
    // ```
    let cutoffs = linking_params.cutoffs;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:751-752
    // ```c
    //                        H_hsp_xsum + score*kbp[H->hsp->context]->Lambda -
    //                        kbp[H->hsp->context]->logK;
    // ```
    // Per-context logK, as in the reference.
    let log_k_by_ctx: Vec<f64> = query_contexts
        .iter()
        .map(|ctx| ctx.karlin_params.k.ln())
        .collect();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
    // ```c
    // /* Sort by (reverse) position. */
    // if (kTranslatedQuery) {
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTbx);
    // ```
    // The reference's `sort_hsps_by_ncbi_link_order(&mut hits, rev_compare_hsps_tbx)`.
    // `allow_parallel` goes to the sort for LOSAT_X_LINKSORT (tblastx/x_link_sort.rs).
    lap.lap(LinkStage::Reset);
    sort_in_place(
        &mut hits,
        rev_compare_hsps_tbx,
        x_link_sort::Site::First(allow_parallel),
    );
    lap.lap(LinkStage::Sort);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:522-531
    // ```c
    // if (H->prev != NULL &&
    //     ((H->hsp->context/strand_factor) !=
    //      (H->prev->hsp->context/strand_factor) ||
    //      (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
    // { /* If frame switches, then start new list. */
    //    hp_frame_number[cur_frame]--;
    //    hp_frame_number[++cur_frame]++;
    //    hp_frame_start[cur_frame] = H;
    //    H->prev->next = NULL;
    //    H->prev = NULL;
    // ```
    let total_hits = hits.len(); // NCBI: hsp_list->hspcnt
    let frame_groups = groups_off_the_end(hits);
    lap.lap(LinkStage::Groups);
    lap.sizes(total_hits, frame_groups.iter().map(Vec::len));
    let kernel_times = lap.kernels();
    GROUPS.fetch_add(frame_groups.len() as u64, AtomicOrdering::Relaxed);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
    // ```c
    // for (frame_index=0; frame_index<num_query_frames; frame_index++)
    // {
    //    hp_start = s_LinkHSPStructReset(hp_start);
    //    hp_start->next = hp_frame_start[frame_index];
    //    hp_frame_start[frame_index]->prev = hp_start;
    //    number_of_hsps = hp_frame_number[frame_index];
    // ```
    // The groups are linked in frame_index order or, as in the reference, on the pool with the
    // results put back into that order.
    crate::utils::threading::report_stage(
        "tblastx",
        "linking",
        frame_groups.len(),
        allow_parallel
            && frame_groups.len() > 1
            && cfg!(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            )),
    );
    let initial_lh_size = (total_hits + 5).max(1024); // NCBI: MAX(1024, hspcnt+5)
    let use_legacy = legacy_linking_requested();
    let literal_port = use_legacy || link_shadow_enabled();
    let args = GroupArgs {
        use_legacy,
        diag_enabled: diagnostics_enabled(),
        params,
        cutoffs: &cutoffs,
        linking_params,
        query_contexts,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        log_k_by_ctx: &log_k_by_ctx,
    };

    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let results: Vec<UngappedHit> = if allow_parallel && frame_groups.len() > 1 {
        thread_local! {
            static POOLS: RefCell<Option<BufferPools>> = const { RefCell::new(None) };
        }

        let mut indexed_results: Vec<(usize, Vec<UngappedHit>)> = frame_groups
            .into_par_iter()
            .enumerate()
            .map(|(group_idx, group_hits)| {
                POOLS.with(|pools| {
                    let mut pools_opt = pools.borrow_mut();
                    let pool = pools_opt.get_or_insert_with(|| {
                        buffer_pools(initial_lh_size, total_hits, literal_port)
                    });
                    (
                        group_idx,
                        kernel_times.time(|| link_one_group(group_hits, &args, pool)),
                    )
                })
            })
            .collect();

        // Sort by group index to preserve NCBI processing order
        indexed_results.sort_by_key(|(idx, _)| *idx);
        lap.lap(LinkStage::Kernels);

        let ordered = assemble(
            indexed_results.into_iter().map(|(_, hits)| hits),
            total_hits,
        );
        lap.lap(LinkStage::Assemble);
        ordered
    } else {
        link_groups_serially(
            frame_groups,
            &args,
            initial_lh_size,
            total_hits,
            literal_port,
            &mut lap,
            &kernel_times,
        )
    };

    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let results: Vec<UngappedHit> = link_groups_serially(
        frame_groups,
        &args,
        initial_lh_size,
        total_hits,
        literal_port,
        &mut lap,
        &kernel_times,
    );

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
    // ```c
    // if (kTranslatedQuery) {
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTransl);
    //    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_FwdCompareHSPsTransl);
    // ```
    // The reference's empty/singleton shortcut, then its final sorts and replay.
    lap.add_kernels(&kernel_times);
    if results.len() <= 1 {
        return results;
    }
    let replayed = replay_in_place(results, &mut lap, allow_parallel);
    lap.lap(LinkStage::Replay);
    replayed
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// ```c
// for (frame_index=0; frame_index<num_query_frames; frame_index++)
// {
//    hp_start = s_LinkHSPStructReset(hp_start);
//    hp_start->next = hp_frame_start[frame_index];
//    hp_frame_start[frame_index]->prev = hp_start;
//    number_of_hsps = hp_frame_number[frame_index];
// ```
// The reference's serial loop: the groups in frame_index order with one scratch. Each kernel
// returns the Vec it was given, which stays where it is until all groups are linked; the list
// that takes them is allocated only then (no kernel runs while that list is being filled).
fn link_groups_serially(
    mut frame_groups: Vec<Vec<UngappedHit>>,
    args: &GroupArgs<'_>,
    initial_lh_size: usize,
    total_hits: usize,
    literal_port: bool,
    lap: &mut LinkLap,
    kernel_times: &KernelTimes,
) -> Vec<UngappedHit> {
    let mut pools = buffer_pools(initial_lh_size, total_hits, literal_port);
    for slot in frame_groups.iter_mut() {
        let group_hits = std::mem::take(slot);
        *slot = kernel_times.time(|| link_one_group(group_hits, args, &mut pools));
    }
    drop(pools);
    lap.lap(LinkStage::Kernels);
    let ordered = assemble(frame_groups.into_iter(), total_hits);
    lap.lap(LinkStage::Assemble);
    ordered
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1773-1776
// ```c
//     /* Remove any information on number of linked HSPs from previous
//        linking. */
//     for (index = 0; index < hsp_list->hspcnt; ++index)
//         hsp_list->hsp_array[index]->num = 1;
// ```
// The dispatch of `apply_sum_stats_even_gap_linking_with_parallel` (called after its empty-list
// return and before it touches the list, only when the switch or its shadow is on).
/// `Ok`: the linked list (this path, or both ways compared and the reference's result). `Err`:
/// the untouched input, for the caller's reference body.
#[allow(clippy::too_many_arguments)]
pub(super) fn linking_dispatch(
    hits: Vec<UngappedHit>,
    params: &KarlinParams,
    linking_params: &LinkingParams,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    allow_parallel: bool,
) -> Result<Vec<UngappedHit>, Vec<UngappedHit>> {
    let m = mode();
    if m == 0 {
        return Err(hits);
    }
    if LINKING_REFERENCE_CALL.with(|flag| flag.replace(false)) {
        // A call for the reference body: the shadow's own (below) or a unit test's.
        return Err(hits);
    }
    if gate() {
        GATE_REFERENCE.fetch_add(1, AtomicOrdering::Relaxed);
        return Err(hits);
    }
    if m == 1 {
        LINKINGS.fetch_add(1, AtomicOrdering::Relaxed);
        HITS.fetch_add(hits.len() as u64, AtomicOrdering::Relaxed);
        return Ok(link_in_place(
            hits,
            params,
            linking_params,
            query_contexts,
            subject_frame_bases,
            length_adj_per_context,
            eff_searchsp_per_context,
            allow_parallel,
        ));
    }

    let reference_input = hits.clone();
    let actual = link_in_place(
        hits,
        params,
        linking_params,
        query_contexts,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        allow_parallel,
    );
    LINKING_REFERENCE_CALL.with(|flag| flag.set(true));
    let reference = apply_sum_stats_even_gap_linking_with_parallel(
        reference_input,
        params,
        linking_params,
        query_contexts,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        allow_parallel,
    );
    assert!(
        !LINKING_REFERENCE_CALL.with(|flag| flag.get()),
        "LOSAT_X_HSPINPLACESHADOW: the reference linking call did not reach its dispatch"
    );
    shadow_compare("linking", &actual, &reference);
    LINKINGS.fetch_add(1, AtomicOrdering::Relaxed);
    HITS.fetch_add(reference.len() as u64, AtomicOrdering::Relaxed);
    // The search continues with the reference's result.
    Ok(reference)
}

#[cfg(test)]
mod tests {
    use super::super::tests::{fast_test_contexts, fast_test_cutoffs, mock_hit};
    use super::*;
    use crate::algorithm::tblastx::blast_gapalign::score_compare_ungapped_hits_ncbi;
    use crate::algorithm::tblastx::ungapped_hit_sort::sort_ungapped_hits_by_index_replay;

    /// xorshift64, for reproducible random lists (as in the linking fuzz test).
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
        /// A random permutation of `0..n` (Fisher-Yates).
        fn shuffled(&mut self, n: usize) -> Vec<usize> {
            let mut values: Vec<usize> = (0..n).collect();
            for index in (1..n).rev() {
                let other = self.below(index as u64 + 1) as usize;
                values.swap(index, other);
            }
            values
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
    // ```c
    //         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //               ScoreCompareHSPs);
    // ```
    // No computation: the tests compare with the reference bodies, and the switch is on by
    // default, so the call goes past the dispatch as the shadow's own call does.
    /// The body of `sort_ungapped_hits_by_index_replay`, whatever the switch says. The dispatch
    /// takes the flag when the call reaches it: not with the switch off, and not for a list of
    /// at most one HSP, which returns before it.
    fn reference_sort(
        hits: &mut [UngappedHit],
        compare: fn(&UngappedHit, &UngappedHit) -> Ordering,
    ) {
        SORT_REFERENCE_CALL.with(|flag| flag.set(true));
        sort_ungapped_hits_by_index_replay(hits, compare);
        let left = SORT_REFERENCE_CALL.with(|flag| flag.replace(false));
        assert_eq!(
            left,
            mode() == 0 || hits.len() <= 1,
            "the dispatch of the reference sort takes the flag when it is reached"
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
    // ```c
    // /* Sort by (reverse) position. */
    // if (kTranslatedQuery) {
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTbx);
    // ```
    // No computation: as `reference_sort`, for the linking call. The reference linking body
    // sorts the list once with the index-replay sort, on the calling thread, before anything
    // else can sort; the second flag sends that sort to its reference body as well.
    /// The body of `apply_sum_stats_even_gap_linking_with_parallel` with the reference sort
    /// inside, whatever the switch says. The linking dispatch takes its flag unless the list is
    /// empty, the sort dispatch unless the list has at most one HSP; neither with the switch off.
    #[allow(clippy::too_many_arguments)]
    fn reference_linking(
        hits: Vec<UngappedHit>,
        params: &KarlinParams,
        linking_params: &LinkingParams,
        query_contexts: &[QueryContext],
        subject_frame_bases: &[i32],
        length_adj_per_context: &[i64],
        eff_searchsp_per_context: &[i64],
        allow_parallel: bool,
    ) -> Vec<UngappedHit> {
        let n = hits.len();
        LINKING_REFERENCE_CALL.with(|flag| flag.set(true));
        SORT_REFERENCE_CALL.with(|flag| flag.set(true));
        let linked = apply_sum_stats_even_gap_linking_with_parallel(
            hits,
            params,
            linking_params,
            query_contexts,
            subject_frame_bases,
            length_adj_per_context,
            eff_searchsp_per_context,
            allow_parallel,
        );
        let linking_left = LINKING_REFERENCE_CALL.with(|flag| flag.replace(false));
        let sort_left = SORT_REFERENCE_CALL.with(|flag| flag.replace(false));
        assert_eq!(
            linking_left,
            mode() == 0 || n == 0,
            "the dispatch of the reference linking call takes the flag when it is reached"
        );
        assert_eq!(
            sort_left,
            mode() == 0 || n <= 1,
            "the sort of the reference linking call takes the flag when it reaches its dispatch"
        );
        linked
    }

    /// The same HSPs in the same order, every field, E-values bit for bit.
    fn assert_same_list(label: &str, actual: &[UngappedHit], expected: &[UngappedHit]) {
        assert_eq!(actual.len(), expected.len(), "{label}: length");
        for (index, (a, b)) in actual.iter().zip(expected).enumerate() {
            assert_eq!(
                a.e_value.to_bits(),
                b.e_value.to_bits(),
                "{label}: e_value of HSP {index}"
            );
            assert_eq!(format!("{a:?}"), format!("{b:?}"), "{label}: HSP {index}");
        }
        // The shadow's comparison must agree with the one above.
        shadow_compare(label, actual, expected);
    }

    const FRAMES: [i8; 6] = [1, 2, 3, -1, -2, -3];

    /// A random HSP of one of `queries` queries: few distinct coordinates and scores when
    /// `span` is small (many comparator-equal HSPs), every record distinct in `hsp_list_order`,
    /// `link_id` and the bits of `e_value`.
    fn random_hit(rng: &mut Rng, index: usize, span: u64, queries: usize) -> UngappedHit {
        let len = 4 + rng.below(40) as usize;
        let q0 = rng.below(span) as usize;
        let s0 = if rng.below(3) == 0 {
            rng.below(span) as usize
        } else {
            q0 + rng.below(25) as usize
        };
        let mut hit = mock_hit(q0, q0 + len, s0, s0 + len, 38 + rng.below(40) as i32);
        hit.ctx_idx = rng.below(6 * queries as u64) as usize;
        hit.q_idx = (hit.ctx_idx / 6) as u32;
        hit.q_frame = FRAMES[hit.ctx_idx % 6];
        hit.s_frame = FRAMES[rng.below(6) as usize];
        hit.s_f_idx = rng.below(6) as usize;
        hit.q_orig_len = 300_000;
        hit.s_orig_len = 300_000;
        hit.hsp_list_order = index;
        hit.link_id = index;
        hit.e_value = f64::from_bits(0x3ff0_0000_0000_0000 + index as u64);
        hit
    }

    fn random_hits(rng: &mut Rng, n: usize, span: u64, queries: usize) -> Vec<UngappedHit> {
        (0..n)
            .map(|index| random_hit(rng, index, span, queries))
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
    // ```c
    //         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //               ScoreCompareHSPs);
    // ```
    // The arrangement the reference builds from a copy (`hits[dst] = original[src].clone()`)
    // against the same arrangement made by swaps: n = 0, 1, 2, identity, reversal, one cycle
    // through every position (both directions), pair swaps and random permutations.
    #[test]
    fn permutation_in_place_equals_copy_and_replay() {
        fn check(source: Vec<usize>, label: &str) {
            let items: Vec<(usize, String)> = (0..source.len())
                .map(|index| (index, format!("record {index}")))
                .collect();
            let expected: Vec<(usize, String)> =
                source.iter().map(|&src| items[src].clone()).collect();
            let mut actual = items.clone();
            let mut work = source.clone();
            apply_permutation_in_place(&mut actual, &mut work);
            assert_eq!(actual, expected, "{label}: {source:?}");
            assert!(
                work.iter().enumerate().all(|(index, &src)| index == src),
                "{label}: the permutation is left as the identity"
            );
        }

        check(Vec::new(), "n = 0");
        check(vec![0], "n = 1");
        check(vec![0, 1], "n = 2 identity");
        check(vec![1, 0], "n = 2 swap");
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        for n in [3usize, 4, 7, 64, 1000, 4097] {
            check((0..n).collect(), "identity");
            check((0..n).rev().collect(), "reversal");
            check((0..n).map(|i| (i + 1) % n).collect(), "one cycle, forward");
            check(
                (0..n).map(|i| (i + n - 1) % n).collect(),
                "one cycle, backward",
            );
            check(
                (0..n)
                    .map(|i| if (i ^ 1) < n { i ^ 1 } else { i })
                    .collect(),
                "pair swaps",
            );
            for round in 0..20 {
                check(rng.shuffled(n), &format!("random {round}"));
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
    // ```c
    //         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //               ScoreCompareHSPs);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486,991-994
    // ```c
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTbx);
    // ...
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTransl);
    //    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_FwdCompareHSPsTransl);
    // ```
    // The in-place sort against the reference index-replay sort, with each of the four NCBI
    // comparators LOSAT sorts HSP lists with, on lists with many comparator-equal HSPs.
    #[test]
    fn sort_in_place_equals_the_index_replay_sort() {
        let comparators: [(&str, fn(&UngappedHit, &UngappedHit) -> Ordering); 4] = [
            ("ScoreCompareHSPs", score_compare_ungapped_hits_ncbi),
            ("s_RevCompareHSPsTbx", rev_compare_hsps_tbx),
            ("s_RevCompareHSPsTransl", rev_compare_hsps_transl),
            ("s_FwdCompareHSPsTransl", fwd_compare_hsps_transl),
        ];
        let mut rng = Rng(0x0123_4567_89AB_CDEF);
        let mut equal_neighbours = 0usize;
        for n in [0usize, 1, 2, 3, 17, 256, 5000] {
            for span in [3u64, 40, 100_000] {
                for (name, compare) in comparators {
                    let hits = random_hits(&mut rng, n, span, 2);
                    let mut expected = hits.clone();
                    reference_sort(&mut expected, compare);
                    let mut actual = hits;
                    sort_in_place(&mut actual, compare, x_link_sort::Site::Dispatch);
                    assert_same_list(&format!("{name} n {n} span {span}"), &actual, &expected);
                    equal_neighbours += actual
                        .windows(2)
                        .filter(|pair| compare(&pair[0], &pair[1]) == Ordering::Equal)
                        .count();
                }
            }
        }
        // The lists must hold comparator-equal HSPs, whose order only stability decides.
        assert!(equal_neighbours > 1000, "{equal_neighbours}");
    }

    /// The grouping loop of `apply_sum_stats_even_gap_linking_with_parallel`, statement for
    /// statement.
    fn push_grouping(hits: Vec<UngappedHit>) -> Vec<Vec<UngappedHit>> {
        let mut frame_groups: Vec<Vec<UngappedHit>> = Vec::new();
        let mut current_group: Vec<UngappedHit> = Vec::new();

        for hit in hits {
            let ctx_group = translated_context_group(&hit);
            let s_sign = hit.s_frame.signum();

            if let Some(prev) = current_group.last() {
                let prev_ctx_group = translated_context_group(prev);
                let prev_s_sign = prev.s_frame.signum();

                if prev_ctx_group != ctx_group || prev_s_sign != s_sign {
                    // Frame switch - start new group
                    if !current_group.is_empty() {
                        frame_groups.push(std::mem::take(&mut current_group));
                    }
                }
            }
            current_group.push(hit);
        }
        if !current_group.is_empty() {
            frame_groups.push(current_group);
        }
        frame_groups
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:522-531
    // ```c
    // if (H->prev != NULL &&
    //     ((H->hsp->context/strand_factor) !=
    //      (H->prev->hsp->context/strand_factor) ||
    //      (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
    // { /* If frame switches, then start new list. */
    //    hp_frame_number[cur_frame]--;
    //    hp_frame_number[++cur_frame]++;
    //    hp_frame_start[cur_frame] = H;
    //    H->prev->next = NULL;
    //    H->prev = NULL;
    // ```
    // The groups cut off the end of the list against the reference's push loop: no HSP, one HSP,
    // one group, many groups (sorted lists of one and of three queries), groups of one HSP each,
    // and a list that arrives with spare capacity. Every group Vec must be exactly sized.
    #[test]
    fn groups_cut_off_the_end_equal_the_push_grouping() {
        fn check(hits: Vec<UngappedHit>, label: &str) -> usize {
            let expected = push_grouping(hits.clone());
            let actual = groups_off_the_end(hits);
            assert_eq!(actual.len(), expected.len(), "{label}: number of groups");
            for (index, (a, b)) in actual.iter().zip(&expected).enumerate() {
                assert_same_list(&format!("{label}: group {index}"), a, b);
                assert_eq!(
                    a.capacity(),
                    a.len(),
                    "{label}: group {index} is exactly sized"
                );
            }
            actual.len()
        }

        let mut rng = Rng(0xD1B5_4A32_D192_ED03);
        assert_eq!(check(Vec::new(), "no HSP"), 0);
        assert_eq!(check(random_hits(&mut rng, 1, 50, 1), "one HSP"), 1);

        // One group: one query strand and one subject strand.
        let mut one_group = random_hits(&mut rng, 300, 50, 1);
        for hit in one_group.iter_mut() {
            hit.q_frame = 2;
            hit.s_frame = -1;
        }
        assert_eq!(check(one_group, "one group"), 1);

        // Sorted lists, as the linking call makes them: 4 groups per query.
        for queries in [1usize, 3] {
            for n in [2usize, 5, 64, 3000] {
                let mut hits = random_hits(&mut rng, n, 60, queries);
                sort_in_place(
                    &mut hits,
                    rev_compare_hsps_tbx,
                    x_link_sort::Site::First(false),
                );
                let groups = check(hits, &format!("sorted, {queries} queries, n {n}"));
                assert!(groups <= 4 * queries);
                if n == 3000 {
                    assert_eq!(groups, 4 * queries);
                }
            }
        }

        // Groups of one HSP each: every neighbour differs in the subject strand or the query.
        let mut singles = random_hits(&mut rng, 97, 50, 1);
        for (index, hit) in singles.iter_mut().enumerate() {
            hit.q_idx = (index / 2) as u32;
            hit.q_frame = 1;
            hit.s_frame = if index % 2 == 0 { 1 } else { -1 };
        }
        assert_eq!(check(singles, "groups of one HSP"), 97);

        // Runs of random lengths, in a list with spare capacity (grown by push).
        let mut runs: Vec<UngappedHit> = Vec::new();
        let mut expected_groups = 0usize;
        let mut index = 0usize;
        for run in 0..40 {
            expected_groups += 1;
            for _ in 0..1 + rng.below(if run % 5 == 0 { 400 } else { 6 }) {
                let mut hit = random_hit(&mut rng, index, 50, 1);
                hit.q_idx = (run / 2) as u32;
                hit.q_frame = 3;
                hit.s_frame = if run % 2 == 0 { -2 } else { 3 };
                runs.push(hit);
                index += 1;
            }
        }
        runs.reserve(runs.len() + 3);
        assert!(runs.capacity() > runs.len());
        assert_eq!(check(runs, "runs with spare capacity"), expected_groups);
    }

    /// The end of the reference linking call: its two final sorts and its replay.
    fn reference_replay(mut results: Vec<UngappedHit>) -> Vec<UngappedHit> {
        let mut references: Vec<&mut UngappedHit> = results.iter_mut().collect();
        sort_linked_hsp_references(&mut references, rev_compare_hsps_transl);
        sort_linked_hsp_references(&mut references, fwd_compare_hsps_transl);
        replay_ncbi_output_list(references)
    }

    /// Random HSPs linked into random chains the way the kernels leave them: unique `link_id`s
    /// in an order unrelated to the list order, heads with `start_of_chain` and the chain's
    /// count, members with `linked_set`, single HSPs with neither link nor `linked_set`.
    fn chained_hits(rng: &mut Rng, n: usize, max_chain: u64) -> Vec<UngappedHit> {
        let mut hits = random_hits(rng, n, 40, 2);
        let ids = rng.shuffled(n);
        for (hit, &id) in hits.iter_mut().zip(&ids) {
            hit.link_id = id;
            hit.num = 1;
            hit.start_of_chain = rng.below(2) == 0;
        }
        let positions = rng.shuffled(n);
        let mut next = 0usize;
        while next < n {
            let len = (1 + rng.below(max_chain) as usize).min(n - next);
            let chain = &positions[next..next + len];
            next += len;
            if len == 1 {
                continue;
            }
            for (k, &position) in chain.iter().enumerate() {
                let next_id = chain.get(k + 1).map(|&p| hits[p].link_id);
                let hit = &mut hits[position];
                hit.linked_set = true;
                hit.start_of_chain = k == 0;
                hit.chain_next_link_id = next_id;
                hit.hsp_link_num = if k == 0 { len as i16 } else { 1 };
                hit.ordering_method = (len % 2) as u8;
            }
        }
        hits
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1018-1020,1044-1061,1080-1082
    // ```c
    //    /* If this is not a single piece or the start of a chain, then Skip it. */
    //    if (H->linked_set == TRUE && H->start_of_chain == FALSE)
    //       continue;
    // ...
    //    /* The first one has the number of links correct. */
    //    num_links = H->hsp_link.num[ordering_method];
    //    link = H->hsp_link.link[ordering_method];
    //    while (link)
    //    {
    //       H->hsp->num = num_links;
    // ...
    //    /* Set these for last link in chain. */
    //    H->hsp->num = num_links;
    // ...
    // for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
    //    hsp_list->hsp_array[index] = H->hsp;
    // ```
    // Lists whose chains hold every HSP exactly once (what linking produces): the in-place
    // replay must give the reference's list, in the list's own allocation.
    #[test]
    fn replay_in_place_equals_the_reference_replay_on_chains() {
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        let mut chained = 0usize;
        for case in 0..400 {
            let n = match case % 8 {
                0 => 2,
                1 => 3,
                2 => 2000,
                _ => 2 + rng.below(300) as usize,
            };
            let max_chain = [1u64, 2, 6, 40][case % 4];
            let hits = chained_hits(&mut rng, n, max_chain);
            let expected = reference_replay(hits.clone());
            let allocation = hits.as_ptr();
            let (actual, in_place) = replay_in_place_counted(hits, &mut LinkLap::off(), false);
            let label = format!("case {case} n {n} max_chain {max_chain}");
            assert!(in_place, "{label}: every HSP is replayed once");
            assert_eq!(actual.as_ptr(), allocation, "{label}: no second list");
            assert_same_list(&label, &actual, &expected);
            chained += actual.iter().filter(|hit| hit.num > 1).count();
        }
        assert!(chained > 10_000, "{chained}");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1018-1020,1046-1056
    // ```c
    //    if (H->linked_set == TRUE && H->start_of_chain == FALSE)
    //       continue;
    // ...
    //    link = H->hsp_link.link[ordering_method];
    //    while (link)
    //    {
    // ...
    //       H->next = (LinkHSPStruct*) link;
    // ...
    //       H = H->next;
    // ```
    // Link fields NCBI's linking does not produce, which the reference replay still answers: a
    // chain member no head reaches (left out), a member two heads reach (replayed twice), both
    // at once (as many entries as HSPs, but not every HSP), a loop, a link to an unknown HSP and
    // a repeated `link_id`. The replay must then return the reference's list (its clones).
    #[test]
    fn replay_falls_back_when_the_replay_is_not_a_permutation() {
        fn base(n: usize) -> Vec<UngappedHit> {
            (0..n)
                .map(|index| {
                    let at = 10 + 30 * index;
                    let mut hit = mock_hit(at, at + 20, at, at + 20, 50 + index as i32);
                    hit.link_id = index;
                    hit.hsp_list_order = index;
                    hit.num = 1;
                    hit.start_of_chain = true;
                    hit
                })
                .collect()
        }
        fn head(hit: &mut UngappedHit, next: usize, count: i16) {
            hit.linked_set = true;
            hit.start_of_chain = true;
            hit.chain_next_link_id = Some(next);
            hit.hsp_link_num = count;
        }
        fn member(hit: &mut UngappedHit) {
            hit.linked_set = true;
            hit.start_of_chain = false;
        }
        fn check(hits: Vec<UngappedHit>, expected_len: usize, expect_in_place: bool, label: &str) {
            let expected = reference_replay(hits.clone());
            assert_eq!(expected.len(), expected_len, "{label}: reference length");
            let (actual, in_place) = replay_in_place_counted(hits, &mut LinkLap::off(), false);
            assert_eq!(in_place, expect_in_place, "{label}");
            assert_same_list(label, &actual, &expected);
        }

        // A member no head reaches is left out.
        let mut hits = base(5);
        member(&mut hits[2]);
        check(hits, 4, false, "orphan member");

        // A member two heads reach is replayed twice, with the count of the later head.
        let mut hits = base(5);
        head(&mut hits[0], 3, 2);
        head(&mut hits[4], 3, 7);
        member(&mut hits[3]);
        check(hits, 6, false, "shared member");

        // Both: as many entries as HSPs, one HSP twice and one missing.
        let mut hits = base(5);
        head(&mut hits[0], 3, 2);
        head(&mut hits[4], 3, 7);
        member(&mut hits[3]);
        member(&mut hits[1]);
        check(hits, 5, false, "shared member and orphan member");

        // A loop: the walk stops after more hops than HSPs (the two single HSPs, then head,
        // member, head, member, head).
        let mut hits = base(4);
        head(&mut hits[1], 2, 2);
        member(&mut hits[2]);
        hits[2].chain_next_link_id = Some(1);
        check(hits, 2 + 5, false, "loop");

        // A link to an HSP that is not in the list ends the chain: still every HSP once.
        let mut hits = base(4);
        head(&mut hits[1], 99, 2);
        check(hits, 4, true, "link to an unknown HSP");

        // A repeated link_id: the later HSP of the sorted order answers for the id.
        let mut hits = base(5);
        head(&mut hits[0], 3, 2);
        member(&mut hits[3]);
        member(&mut hits[4]);
        hits[4].link_id = 3;
        check(hits, 4, false, "repeated link_id");
    }

    /// Two queries of six contexts each, with the Karlin block of the linking tests.
    fn two_query_contexts() -> (KarlinParams, Vec<QueryContext>) {
        let (params, contexts, _) = fast_test_contexts();
        let mut all = contexts.clone();
        for (index, context) in contexts.iter().enumerate() {
            let mut second = context.clone();
            second.q_idx = 1;
            second.frame_base = ((6 + index) * 100_001) as i32;
            all.push(second);
        }
        (params, all)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,483-486,990-994
    // ```c
    // s_BlastEvenGapLinkHSPs(EBlastProgramType program_number, BlastHSPList* hsp_list,
    // ...
    // /* Sort by (reverse) position. */
    // if (kTranslatedQuery) {
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTbx);
    // ...
    // if (kTranslatedQuery) {
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTransl);
    // ```
    // The whole linking call: this path against `apply_sum_stats_even_gap_linking_with_parallel`
    // on random HSP lists of one and two queries (up to eight groups, chains, comparator-equal
    // HSPs), serially and on the pool, with the cutoff sets of the kernel fuzz test.
    #[test]
    fn linking_in_place_equals_the_reference_call() {
        let (params, contexts) = two_query_contexts();
        let cutoff_sets = [
            fast_test_cutoffs(),
            LinkHspCutoffs {
                cutoff_small_gap: 20,
                cutoff_big_gap: 20,
                gap_prob: 0.5,
                ignore_small_gaps: false,
            },
            LinkHspCutoffs {
                cutoff_small_gap: 0,
                cutoff_big_gap: 30,
                gap_prob: 0.0,
                ignore_small_gaps: true,
            },
        ];
        let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(400);
        let mut rng = Rng(0xA076_1D64_78BD_642F);
        let mut total_hsps = 0usize;
        let mut total_chained = 0usize;
        let mut most_groups = 0usize;
        for case in 0..cases {
            let n = match case % 10 {
                0 => case / 10 % 4,
                1 => 2500,
                _ => 1 + rng.below(400) as usize,
            };
            let queries = 1 + case % 2;
            let span = [60u64, 400, 3000, 20_000][case / 2 % 4];
            let mut hits = random_hits(&mut rng, n, span, queries);
            if case % 3 == 0 {
                // A list that arrives with spare capacity, as the search's lists do.
                hits.reserve(n + 7);
            }
            let linking = LinkingParams {
                subject_len_nucl: 300_000,
                gap_decay_rate: 0.5,
                cutoffs: cutoff_sets[case % cutoff_sets.len()],
            };
            for allow_parallel in [false, true] {
                let label = format!(
                    "case {case} n {n} queries {queries} span {span} parallel {allow_parallel}"
                );
                let expected = reference_linking(
                    hits.clone(),
                    &params,
                    &linking,
                    &contexts,
                    &[0; 6],
                    &[49; 12],
                    &[9_372_428_362; 12],
                    allow_parallel,
                );
                let actual = link_in_place(
                    hits.clone(),
                    &params,
                    &linking,
                    &contexts,
                    &[0; 6],
                    &[49; 12],
                    &[9_372_428_362; 12],
                    allow_parallel,
                );
                assert_same_list(&label, &actual, &expected);
                assert_eq!(actual.len(), n, "{label}: every HSP is in the result");
                if allow_parallel {
                    total_hsps += actual.len();
                    total_chained += actual.iter().filter(|hit| hit.linked_set).count();
                }
            }
            let mut sorted = hits;
            sort_in_place(
                &mut sorted,
                rev_compare_hsps_tbx,
                x_link_sort::Site::First(false),
            );
            most_groups = most_groups.max(groups_off_the_end(sorted).len());
        }
        if cases >= 400 {
            // The lists must exercise chains and every group of two queries.
            assert!(
                total_chained * 10 > total_hsps,
                "{total_chained} of {total_hsps}"
            );
            assert_eq!(most_groups, 8);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1082
    // ```c
    // for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
    //    hsp_list->hsp_array[index] = H->hsp;
    // ```
    // The shadow's comparison must see a difference in any field, also one of E-value bits only
    // and one of length.
    #[test]
    fn shadow_compare_sees_every_kind_of_difference() {
        let mut rng = Rng(0x5851_F42D_4C95_7F2D);
        let reference = random_hits(&mut rng, 6, 50, 1);
        shadow_compare("equal lists", &reference, &reference.clone());

        let differs = |change: &dyn Fn(&mut Vec<UngappedHit>)| {
            let mut actual = reference.clone();
            change(&mut actual);
            std::panic::catch_unwind(|| shadow_compare("test", &actual, &reference)).is_err()
        };
        assert!(differs(&|hits| hits[3].num += 1));
        assert!(differs(&|hits| hits[0].chain_next_link_id = Some(0)));
        assert!(differs(&|hits| hits[5].hsp_link_num = 9));
        assert!(differs(
            &|hits| hits[2].start_of_chain = !hits[2].start_of_chain
        ));
        assert!(differs(&|hits| hits.swap(1, 4)));
        assert!(differs(&|hits| {
            hits.pop();
        }));
        // 0.0 and -0.0 are equal as numbers; the comparison is of bits.
        let mut zero = reference.clone();
        zero[1].e_value = 0.0;
        let mut negative_zero = zero.clone();
        negative_zero[1].e_value = -0.0;
        assert!(
            std::panic::catch_unwind(|| shadow_compare("zero", &negative_zero, &zero)).is_err()
        );
    }
}
