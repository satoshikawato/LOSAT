//! The TBLASTX HSP list is sorted, linked as sub-slices and replayed into NCBI's output order
//! without a second whole copy of the list alive. Only where the records lie changes; every
//! comparison, every per-group linking call and every value written is the reference's.
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
//! of the HSPs. LOSAT owns complete `UngappedHit` records (88 bytes); the first LOSAT form copied
//! the whole list at each of these steps (`hits.to_vec()` in every sort, groups grown by `push`
//! while the sorted list is drained, a result grown by `extend`, a cloned replay). This path
//! keeps one list. The statements of that first form are described below as "the reference". Its
//! linking body in `linking.rs` (the push grouping, the clones, the two reference sorts and the
//! replay) remains: `LOSAT_DEBUG_CUTOFFS` and `LOSAT_DIAGNOSTICS` still take it (`gate()`), and
//! the unit tests compare this path with it. Its index-replay sort copy is gone.
//!
//! # Why the result is the reference's (each point checked against the code it replaced)
//!
//! 1. Sort (`sort_in_place`). The reference sorted the indices `0..n` with `slice::sort_by` and
//!    a comparator that reads only the two records, taken from an untouched copy of the list.
//!    This path sorts the same initial indices (`x_link_sort::sort_indices`, which returns the
//!    same order as that `sort_by`) and reads the records from the list itself, which is not
//!    written during the sort: the same records at the same indices, so the same comparison
//!    results in the same sequence and the same index order. The reference then wrote
//!    `hits[dst] = original[indices[dst]]`; `apply_permutation_in_place` produces exactly that
//!    arrangement by walking each cycle of the permutation with swaps (after the walk of a cycle
//!    every position `dst` of it holds the record that was at `indices[dst]`, and the cycle's
//!    first record has travelled to the position that asks for it). No record is cloned, so no
//!    field can change.
//! 2. Groups (`group_starts`, `x_link_slices.rs`). The reference started a new group at a record
//!    whose `translated_context_group` or sign of `s_frame` differs from the record before it
//!    (`current_group.push(hit)` loop), so its groups are the maximal runs of consecutive records
//!    of the sorted list, in list order. `group_starts` finds the same run starts with the same
//!    two tests on the same neighbours, and `x_link_slices.rs` links the runs where they lie, as
//!    disjoint sub-slices of the one list: the same groups with the same records in the same
//!    order, with no cut and no second list.
//! 3. Per-group linking. The same functions (`link_hsp_group_incremental`, or
//!    `link_hsp_group_ncbi` when the literal port is asked for) with the same arguments, serially
//!    in group order or on the pool. The scratch buffers handed to the kernels start empty
//!    unless the literal port can run (`LOSAT_LINKING_LEGACY`, a linking trace, or
//!    `LOSAT_LINK_FAST_SHADOW`): the incremental kernel never touches them, and the literal port
//!    clears and reserves them itself for every group (`pool_hsp_links.clear(); .reserve(n)`,
//!    `pool_lh_helpers.clear(); .resize(..)`), so their initial capacity changes no value.
//! 4. Assembly. The reference concatenated the linked groups in group order. The groups lie in
//!    group order inside the one list, so there is nothing to put together.
//! 5. Output order (`replay_in_place`; reference: two `sort_linked_hsp_references` calls and
//!    `replay_ncbi_output_list`, which the `LOSAT_DEBUG_CUTOFFS` / `LOSAT_DIAGNOSTICS` body still
//!    runs). The reference stable-sorts a Vec of references that starts in
//!    list order, twice; this path stable-sorts the indices `0..n` with the same two comparators
//!    applied to the same records (`x_link_sort::final_sort`, one sort with the same result), so
//!    `order[i]` is the list position of the reference's `sorted_hits[i]`. The walk is the
//!    reference's statement for statement with `results[order[i]]` in place of `sorted_hits[i]`:
//!    the same `id_to_index`, the same skip of chain members, the same `num` writes, the same
//!    stops, hence the same `replayed_indices`. The reference returns clones of
//!    `sorted_hits[replayed_indices[k]]`, taken after all `num` writes. When `replayed_indices`
//!    is a permutation of `0..n` (every record exactly once, which is what linked chains give),
//!    moving the records by that permutation gives the same list; otherwise (a record left out
//!    or reached twice) this path returns the reference's own expression, the clones. Both
//!    outcomes are counted.
//!
//! A linking call with `LOSAT_DEBUG_CUTOFFS` or `LOSAT_DIAGNOSTICS` set runs the reference body
//! in `linking.rs` (its debug output lives there); its first sort is the sort of item 1.
//!
//! The index sort sites of this file (`sort_in_place`, `replay_in_place_counted`) call
//! `x_link_sort.rs`, which returns the index order of the reference's sorts from one final sort
//! instead of two and from a stable parallel sort for long lists (its header has the argument).
//! The group step of `link_in_place` is `x_link_slices.rs`, which links the groups as disjoint
//! sub-slices of the one sorted list (its header has the argument).

use super::*;
// Diagnostic (LOSAT_TIMING): the clock of a linking call (tblastx/x_link_timing.rs). It reads
// the clock only in a call the batch function armed, and the linking reads nothing from it.
use crate::algorithm::tblastx::x_link_timing::{LinkLap, LinkStage};
#[cfg(test)]
use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

// No NCBI counterpart: counters for LOSAT_X_STATS; they do not change any value NCBI computes.
static SORTS: AtomicU64 = AtomicU64::new(0);
static LINKINGS: AtomicU64 = AtomicU64::new(0);
static HITS: AtomicU64 = AtomicU64::new(0);
static GROUPS: AtomicU64 = AtomicU64::new(0);
static REPLAY_INPLACE: AtomicU64 = AtomicU64::new(0);
static REPLAY_FALLBACKS: AtomicU64 = AtomicU64::new(0);
static GATE_REFERENCE: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: how a unit test reaches the reference linking body, which has the dispatch
// at its head. The test sets the flag of the calling thread and calls the reference function; the
// dispatch at the head of that call takes the flag and lets the body run. Nothing else sets the
// flag (`tests::reference_linking`), so a search never reads it.
#[cfg(test)]
thread_local! {
    static LINKING_REFERENCE_CALL: Cell<bool> = const { Cell::new(false) };
}

// No NCBI counterpart: the reference linking body prints the LOSAT_DEBUG_CUTOFFS block and passes
// the LOSAT_DIAGNOSTICS flag to the literal port; with either set the reference body runs.
fn gate() -> bool {
    static GATE: OnceLock<bool> = OnceLock::new();
    *GATE.get_or_init(|| std::env::var_os("LOSAT_DEBUG_CUTOFFS").is_some() || diagnostics_enabled())
}

/// The summary line (printed by `main` at exit): the LOSAT_X_STATS counters.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_none() {
        return;
    }
    let load = |c: &AtomicU64| c.load(AtomicOrdering::Relaxed);
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
    assert_eq!(items.len(), source.len(), "permutation length");
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
// The index sort of `sort_ungapped_hits_by_index_replay` (same comparator, same records at the
// same indices), then the records moved instead of copied and cloned back. `site` names the
// caller for `x_link_sort`.
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
        // The same index order from a stable parallel sort when the list is long and the call may
        // use the pool (tblastx/x_link_sort.rs, item 2 of its header).
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
        x_link_sort::sort_indices(&mut indices, original, compare, site);
    }
    apply_permutation_in_place(hits, &mut indices);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1082
// ```c
// for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
//    hsp_list->hsp_array[index] = H->hsp;
// ```
// No computation: the unit tests of this file and of `x_link_slices.rs` compare two lists record
// by record, every field, `e_value` by bits, and stop at the first difference. A new field of
// `UngappedHit` does not compile here until it is compared.
/// Test check: `actual` (the in-place or sliced path) against `reference`, record by record.
#[cfg(test)]
pub(super) fn x_assert_same_hits(what: &str, actual: &[UngappedHit], reference: &[UngappedHit]) {
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
                "{what}, index {k} differs (in place {} hits, reference {} hits): in place {a:?} vs reference {b:?}",
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
// The sort of `sort_ungapped_hits_by_index_replay` (called after its `len <= 1` return).
/// Sorts `hits` in place.
pub(crate) fn sort_dispatch(
    hits: &mut [UngappedHit],
    compare: fn(&UngappedHit, &UngappedHit) -> Ordering,
) {
    sort_in_place(hits, compare, x_link_sort::Site::Dispatch);
    SORTS.fetch_add(1, AtomicOrdering::Relaxed);
    HITS.fetch_add(hits.len() as u64, AtomicOrdering::Relaxed);
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
// run starts as the reference's push loop: the test of a record against the one before it. This
// loop is the one the reference's push loop ran; `x_link_slices.rs` finds its groups with this
// function.
/// The list positions where a group of the reference's grouping starts, after the first: every
/// `index >= 1` whose record differs from the record before it in `translated_context_group` or
/// in the sign of `s_frame`, in increasing order.
pub(super) fn group_starts(hits: &[UngappedHit]) -> Vec<usize> {
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
    starts
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
// the reference's capacity is 88 bytes per HSP of the whole list for every linking thread.
// `x_link_slices.rs` makes its scratch with this function.
pub(super) fn buffer_pools(
    initial_lh_size: usize,
    total_hits: usize,
    literal_port: bool,
) -> BufferPools {
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

/// The values every per-group linking call of one linking call shares (`x_link_slices.rs` hands
/// them to the per-group functions).
pub(super) struct GroupArgs<'a> {
    pub(super) use_legacy: bool,
    pub(super) diag_enabled: bool,
    pub(super) params: &'a KarlinParams,
    pub(super) cutoffs: &'a LinkHspCutoffs,
    pub(super) linking_params: &'a LinkingParams,
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) subject_frame_bases: &'a [i32],
    pub(super) length_adj_per_context: &'a [i64],
    pub(super) eff_searchsp_per_context: &'a [i64],
    pub(super) log_k_by_ctx: &'a [f64],
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
// The reference's two final sorts as one composite sort on indices (`x_link_sort::final_sort`),
// the walk of `replay_ncbi_output_list` with `results[order[i]]` for its `sorted_hits[i]`, and
// the records moved into the replayed order when it holds every record exactly once (module
// doc, item 5).
/// The result of the reference's two final sorts (one composite sort here) and
/// `replay_ncbi_output_list` for a list of at least two HSPs, without a second list unless the
/// replay leaves a record out or repeats one.
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
// `allow_parallel` is the linking call's, for `x_link_sort`.
/// The replayed list and whether it was made in place (`true`) or by the reference's clones.
fn replay_in_place_counted(
    mut results: Vec<UngappedHit>,
    lap: &mut LinkLap,
    allow_parallel: bool,
) -> (Vec<UngappedHit>, bool) {
    let n = results.len();
    let mut order: Vec<usize> = (0..n).collect();
    // The index order of the two reference sorts from one stable sort with
    // `s_FwdCompareHSPsTransl`, ties by `s_RevCompareHSPsTransl`, on the pool when the list is
    // long and the call may use it (tblastx/x_link_sort.rs, items 1 and 2 of its header).
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
    // ```c
    // if (kTranslatedQuery) {
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTransl);
    //    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_FwdCompareHSPsTransl);
    // ```
    x_link_sort::final_sort(
        &mut order,
        &results,
        x_link_sort::Site::Final(allow_parallel),
    );
    lap.lap(LinkStage::FinalSorts);

    let max_link_id = results
        .iter()
        .map(|h| h.link_id as usize)
        .max()
        .unwrap_or(0);
    let mut id_to_index = vec![SENTINEL_IDX; max_link_id + 1];
    for (index, &position) in order.iter().enumerate() {
        let hit = &results[position];
        if (hit.link_id as usize) < id_to_index.len() {
            id_to_index[hit.link_id as usize] = index;
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
            let next_id = next_id as usize;
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
        hit.link_id = u32::try_from(link_id)
            .expect("link_id: a position in the linking list fits u32 (NCBI's hspcnt is Int4)");
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
    // `allow_parallel` goes to the sort (tblastx/x_link_sort.rs).
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
    // The groups linked as disjoint sub-slices of the one sorted list, with no cut and no second
    // list (tblastx/x_link_slices.rs, exactness argument in its header).
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
    // ```c
    // for (frame_index=0; frame_index<num_query_frames; frame_index++)
    // {
    //    hp_start = s_LinkHSPStructReset(hp_start);
    //    hp_start->next = hp_frame_start[frame_index];
    //    hp_frame_start[frame_index]->prev = hp_start;
    //    number_of_hsps = hp_frame_number[frame_index];
    // ```
    // NCBI links its groups as sublists of the one sorted array.
    let (results, groups) = x_link_slices::link_groups(
        hits,
        &x_link_slices::Call {
            params,
            cutoffs: &cutoffs,
            linking_params,
            query_contexts,
            subject_frame_bases,
            length_adj_per_context,
            eff_searchsp_per_context,
            log_k_by_ctx: &log_k_by_ctx,
            allow_parallel,
        },
        &mut lap,
    );
    GROUPS.fetch_add(groups as u64, AtomicOrdering::Relaxed);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
    // ```c
    // if (kTranslatedQuery) {
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTransl);
    //    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_FwdCompareHSPsTransl);
    // ```
    // The reference's empty/singleton shortcut, then its final sorts and replay.
    if results.len() <= 1 {
        return results;
    }
    let replayed = replay_in_place(results, &mut lap, allow_parallel);
    lap.lap(LinkStage::Replay);
    replayed
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1773-1776
// ```c
//     /* Remove any information on number of linked HSPs from previous
//        linking. */
//     for (index = 0; index < hsp_list->hspcnt; ++index)
//         hsp_list->hsp_array[index]->num = 1;
// ```
// The dispatch of `apply_sum_stats_even_gap_linking_with_parallel` (called after its empty-list
// return and before it touches the list).
/// `Ok`: the linked list (this path). `Err`: the untouched input, for the caller's reference body
/// (`LOSAT_DEBUG_CUTOFFS` or `LOSAT_DIAGNOSTICS`, or a unit test's reference call).
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
    // A unit test's call for the reference body (`tests::reference_linking`).
    #[cfg(test)]
    if LINKING_REFERENCE_CALL.with(|flag| flag.replace(false)) {
        return Err(hits);
    }
    if gate() {
        GATE_REFERENCE.fetch_add(1, AtomicOrdering::Relaxed);
        return Err(hits);
    }
    LINKINGS.fetch_add(1, AtomicOrdering::Relaxed);
    HITS.fetch_add(hits.len() as u64, AtomicOrdering::Relaxed);
    Ok(link_in_place(
        hits,
        params,
        linking_params,
        query_contexts,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        allow_parallel,
    ))
}

#[cfg(test)]
mod tests {
    use super::super::tests::{fast_test_contexts, fast_test_cutoffs, mock_hit};
    use super::*;

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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
    // ```c
    // /* Sort by (reverse) position. */
    // if (kTranslatedQuery) {
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTbx);
    // ```
    // No computation: the test calls the reference linking body past the dispatch. That body
    // sorts the list once with the index-replay sort, which is the in-place sort here as well.
    /// The body of `apply_sum_stats_even_gap_linking_with_parallel`. The linking dispatch takes
    /// the flag unless the list is empty.
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
        assert_eq!(
            linking_left,
            n == 0,
            "the dispatch of the reference linking call takes the flag when it is reached"
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
        // The comparison of `x_assert_same_hits` must agree with the one above.
        x_assert_same_hits(label, actual, expected);
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
        hit.ctx_idx = u32::try_from(rng.below(6 * queries as u64)).unwrap();
        hit.q_idx = (hit.ctx_idx as usize / 6) as u32;
        hit.q_frame = FRAMES[hit.ctx_idx as usize % 6];
        hit.s_frame = FRAMES[rng.below(6) as usize];
        hit.s_f_idx = u8::try_from(rng.below(6)).unwrap();
        hit.q_orig_len = 300_000;
        hit.s_orig_len = 300_000;
        hit.hsp_list_order = u32::try_from(index).unwrap();
        hit.link_id = u32::try_from(index).unwrap();
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
            hit.link_id = u32::try_from(id).unwrap();
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
                    hit.link_id = u32::try_from(index).unwrap();
                    hit.hsp_list_order = u32::try_from(index).unwrap();
                    hit.num = 1;
                    hit.start_of_chain = true;
                    hit
                })
                .collect()
        }
        fn head(hit: &mut UngappedHit, next: usize, count: i16) {
            hit.linked_set = true;
            hit.start_of_chain = true;
            hit.chain_next_link_id = Some(u32::try_from(next).unwrap());
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
            most_groups =
                most_groups.max(group_starts(&sorted).len() + usize::from(!sorted.is_empty()));
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
}
