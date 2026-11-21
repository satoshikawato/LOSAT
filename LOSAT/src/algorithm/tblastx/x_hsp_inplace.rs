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
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:476-486
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:513-534
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994,1012-1020,1044-1061,1080-1082
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
use crate::utils::glibc_fma_libm::GlibcLibm;
// Diagnostic (LOSAT_TIMING): the clock of a linking call (tblastx/x_link_timing.rs). It reads
// the clock only in a call the batch function armed, and the linking reads nothing from it.
use crate::algorithm::tblastx::x_link_timing::{LinkLap, LinkStage};
use crate::utils::x_order::apply_permutation_in_place;
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
// The index sort of `sort_ungapped_hits_by_index_replay` (same comparator, same records at the
// same indices), then the records moved instead of copied and cloned back. `site` names the
// caller for `x_link_sort`.
/// The result of `sort_ungapped_hits_by_index_replay` without a copy of the list.
fn sort_in_place<C>(hits: &mut [UngappedHit], compare: C, site: x_link_sort::Site)
where
    C: Fn(&UngappedHit, &UngappedHit) -> Ordering + Sync,
{
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
        x_link_sort::sort_indices(&mut indices, original, compare, site);
    }
    apply_permutation_in_place(hits, &mut indices);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1082
// No computation: the unit tests of this file and of `x_link_slices.rs` compare two lists record
// by record, every field, `e_value` by bits, and stop at the first difference. A new field of
// `UngappedHit` does not compile here until it is compared.
/// Test check: `actual` (the in-place or sliced path) against `reference`, record by record.
#[cfg(test)]
pub(crate) fn x_assert_same_hits(what: &str, actual: &[UngappedHit], reference: &[UngappedHit]) {
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
// The sort of `sort_ungapped_hits_by_index_replay` (called after its `len <= 1` return).
/// Sorts `hits` in place.
pub(crate) fn sort_dispatch<C>(hits: &mut [UngappedHit], compare: C)
where
    C: Fn(&UngappedHit, &UngappedHit) -> Ordering + Sync,
{
    sort_in_place(hits, compare, x_link_sort::Site::Dispatch);
    SORTS.fetch_add(1, AtomicOrdering::Relaxed);
    HITS.fetch_add(hits.len() as u64, AtomicOrdering::Relaxed);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:513-534
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1012-1020
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1044-1061
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1082
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// The steps of `apply_sum_stats_even_gap_linking_with_parallel` in its order, each with the
// reference's functions and values, on one list (module doc, items 1-5). `compare` is the
// comparator of the sort (`rev_compare_hsps_tbx`; `rev_compare_hsps_tbx_int4` for a list of the
// one-hit word finder that holds an HSP of length 0 or less, `link_in_place_int4`), and
// `link_modes` gives the group step's `use_legacy` and `fast_shadow` (`x_link_slices::link_modes`
// reads them from the environment; the one-hit form asks for the literal port). Each caller's
// instance is its own, so the default one is the code of this function with those two fixed.
/// The result of `apply_sum_stats_even_gap_linking_with_parallel` without a second whole list.
#[allow(clippy::too_many_arguments)]
fn link_in_place<C, M>(
    mut hits: Vec<UngappedHit>,
    params: &KarlinParams,
    linking_params: &LinkingParams,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    allow_parallel: bool,
    compare: C,
    link_modes: M,
) -> Vec<UngappedHit>
where
    C: Fn(&UngappedHit, &UngappedHit) -> Ordering + Sync,
    M: FnOnce() -> (bool, bool),
{
    // No NCBI counterpart: the clock of this call for LOSAT_TIMING. The first local, so that it
    // is dropped last.
    let mut lap = LinkLap::start();
    if hits.is_empty() {
        return hits;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1773-1776
    for (link_id, hit) in hits.iter_mut().enumerate() {
        hit.link_id = u32::try_from(link_id)
            .expect("link_id: a position in the linking list fits u32 (NCBI's hspcnt is Int4)");
        hit.chain_next_link_id = None;
        hit.num = 1;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:492-493
    let cutoffs = linking_params.cutoffs;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:751-752
    // Per-context logK, as in the reference.
    let log_k_by_ctx: Vec<f64> = query_contexts
        .iter()
        .map(|ctx| ctx.karlin_params.k.glibc_ln())
        .collect();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
    // The reference's `sort_hsps_by_ncbi_link_order(&mut hits, rev_compare_hsps_tbx)`.
    // `allow_parallel` goes to the sort (tblastx/x_link_sort.rs).
    lap.lap(LinkStage::Reset);
    sort_in_place(&mut hits, compare, x_link_sort::Site::First(allow_parallel));
    lap.lap(LinkStage::Sort);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:522-531
    // The groups linked as disjoint sub-slices of the one sorted list, with no cut and no second
    // list (tblastx/x_link_slices.rs, exactness argument in its header).
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
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
        link_modes,
    );
    GROUPS.fetch_add(groups as u64, AtomicOrdering::Relaxed);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
    // The reference's empty/singleton shortcut, then its final sorts and replay.
    if results.len() <= 1 {
        return results;
    }
    let replayed = replay_in_place(results, &mut lap, allow_parallel);
    lap.lap(LinkStage::Replay);
    replayed
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486 and :540-551
// `link_in_place` for a list of the one-hit word finder that holds an HSP of length 0 or less
// (`linking.rs`, `link_one_hit_with_int4_ends`): the sort's comparator reads the ends as NCBI's
// `Int4` (`rev_compare_hsps_tbx_int4`), and every group goes to the literal port
// `link_hsp_group_ncbi`, which reads the four coordinates as `Int4` (`frame_relative_coords`) and
// so computes NCBI's trims, windows and sums on them. The replay's comparators read the starts
// only. The incremental kernel is not used: it visits the HSPs by decreasing end and assumes every
// end is after its start (`linking_incr.rs`). The reference body of LOSAT_DEBUG_CUTOFFS and
// LOSAT_DIAGNOSTICS is not used either (`gate`): it may take that kernel, and its prints are not
// values of the search. The LOSAT_X_STATS counters count the call as `linking_dispatch` does.
/// `link_in_place` with `Int4` ends and the literal port, for a one-hit list.
#[allow(clippy::too_many_arguments)]
pub(super) fn link_in_place_int4(
    hits: Vec<UngappedHit>,
    params: &KarlinParams,
    linking_params: &LinkingParams,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    allow_parallel: bool,
) -> Vec<UngappedHit> {
    LINKINGS.fetch_add(1, AtomicOrdering::Relaxed);
    HITS.fetch_add(hits.len() as u64, AtomicOrdering::Relaxed);
    link_in_place(
        hits,
        params,
        linking_params,
        query_contexts,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        allow_parallel,
        rev_compare_hsps_tbx_int4,
        || (true, false),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1773-1776
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
        rev_compare_hsps_tbx,
        x_link_slices::link_modes,
    ))
}

#[cfg(test)]
#[path = "x_hsp_inplace_tests.rs"]
mod tests;
