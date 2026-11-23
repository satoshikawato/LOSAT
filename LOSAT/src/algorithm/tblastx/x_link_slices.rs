//! The groups of a TBLASTX linking call linked as disjoint sub-slices of the one sorted HSP
//! list. The first form of the in-place path cut the sorted list into one Vec per group, linked
//! each, and copied the linked groups into a new list of the whole size, so that the list existed
//! twice while the second one was filled (AP027131 x AP027133, genetic code 4: 49 MiB twice in
//! the first linking call, the heap peak of the search). Here the groups are linked where they
//! lie: no cut, no assembly, no second list. The list after the per-group linking is the
//! reference's list, record for record, at any number of threads.
//!
//! This file is a child module of `sum_stats_linking/linking.rs` (declared there with `#[path]`,
//! beside `x_hsp_inplace.rs` and `x_link_sort.rs`), so it calls that file's per-group functions
//! without changing them. `x_hsp_inplace::link_in_place` calls it at its group step; nothing else
//! does. The reference body (`LOSAT_DIAGNOSTICS`, `LOSAT_DEBUG_CUTOFFS`) has its own group loop
//! (the `current_group.push(hit)` grouping and the `extend` of the linked groups) and is not
//! touched; the unit tests below compare this file with a copy of that loop.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:522-531
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1082
//! NCBI's groups are sublists of the one sorted array (`hp_frame_start[frame_index]` and
//! `hp_frame_number[frame_index]`), each linked in place, and the array is read back in chain
//! order: the slices are closer to NCBI than the cut Vecs. NCBI has no second array of the HSPs.
//!
//! # Why the list after the per-group linking is the reference's (each point checked)
//!
//! 1. The slices are the reference's groups (same test, same order). The reference body starts
//!    a group at a record whose `translated_context_group` or sign of `s_frame` differs from the
//!    record before it, so its group k holds the records at positions [s_k, s_{k+1}) of the
//!    sorted list, in list order, with s_0 = 0 and the last end at n. `link_groups_as` takes the
//!    starts s_k from `x_hsp_inplace::group_starts`, the same two tests on the same neighbours of
//!    the same sorted list (the reset and the sort before it are `link_in_place`'s), and
//!    `group_slices` splits with `split_at_mut(s)` from the last start down, the head last,
//!    reversed: slice k is positions [s_k, s_{k+1}) of the same list. The same records in the same
//!    order; only where they lie differs (the reference moves them into Vecs of their own, a
//!    bitwise move; the slices leave them in the list).
//! 2. A kernel's result depends only on its group's records, the shared arguments and its
//!    scratch, all the same as in the reference. The arguments: `link_groups_as` builds
//!    `GroupArgs` with the expressions of the reference body (`legacy_linking_requested()`,
//!    `diagnostics_enabled()`, the same references, the larger of 1024 and `hspcnt` plus 5), and
//!    the scratch with the same `buffer_pools`. The kernels: `link_one_slice` takes the branch of
//!    the reference's per-group call on `use_legacy`. The default branch calls
//!    `link_hsp_group_incremental_slice` (`linking.rs`), the statements of
//!    `link_hsp_group_incremental` on the slice, which calls `link_hsp_group_incr_slice`
//!    (`linking_incr.rs`); `link_hsp_group_incr` only passes its Vec to that function as a slice
//!    and returns the Vec, and both reach the one `link_group`. The body of `link_group` takes a
//!    slice, and every operation of the body on `group_hits` (`len`, indexing, `iter`) is a slice
//!    operation, which a `Vec` reaches through `Deref` on the same records. With
//!    `LOSAT_LINKING_LEGACY` (or a linking trace) `link_one_slice` passes `to_vec()` of the slice,
//!    a clone of every record (`UngappedHit` derives `Clone` and every field is a plain value:
//!    the same bits), to `link_hsp_group_ncbi` and writes its result back with
//!    `clone_from_slice`, record for record.
//! 3. A kernel writes only fields of its own records and changes no order, so the list after the
//!    kernels equals the reference's concatenated list record for record. `link_group` reads
//!    `group_hits.len()`, `group_hits[0].ctx_idx` and `group_hits.iter()`, writes
//!    `group_hits[best_i].hsp_link_num` and five fields of `&mut group_hits[cur]`, and returns
//!    the slice it was given; it has no push, removal, swap or sort. The literal port
//!    `link_hsp_group_ncbi` writes `start_of_chain`, `hsp_link_num`, `ordering_method`,
//!    `e_value`, `linked_set` and `chain_next_link_id` of `group_hits[i]` by index, otherwise only
//!    reads the records (indexing, `iter`, `&group_hits` to its scan), and returns the Vec it was
//!    given: its result has the group's length and order. This is the finding on the literal
//!    port; `clone_from_slice`, which panics on unequal lengths, puts record i back at position
//!    i. The reference concatenates the linked groups in group order (on the pool after
//!    `sort_by_key` by group index); here the groups are already in group order inside the
//!    list, and no kernel writes outside its slice (the borrow checker gives each kernel its
//!    slice alone).
//! 4. The per-thread scratch is cleared or sized per group as in the reference. It is made by
//!    `buffer_pools` with the reference's arguments, once per call serially and once per worker
//!    thread on the pool (a thread-local of this file, made on first use, as the reference's
//!    `POOLS`). The incremental kernel does not read it; it reaches the literal port only under
//!    pool_hsp_links.reserve(n)` and `pool_lh_helpers.clear(); pool_lh_helpers.resize(..)` for
//!    every group before it reads them, so neither their contents nor their capacity at entry
//!    changes a value.
//! 5. The final sorts and the replay then get the same list: `link_in_place` passes the result to
//!    the `len <= 1` return and `replay_in_place`, as the reference body passes its list to its
//!    own. The list's capacity differs from the reference's (`shrink_to_fit` here, an exact
//!    `with_capacity` there); nothing reads a capacity.
//!
//! The kernels of different groups touch disjoint records and their own scratch, so the order in
//! which the pool runs them changes nothing: the result does not depend on the number of threads.
//! The pool is used under the reference's condition (`allow_parallel && groups > 1` and the `cfg`
//! of its per-group linking), so the search uses no more threads than `-num_threads`.
//!
//! # Counters
//!
//! With `LOSAT_X_STATS`, `main` prints `[X_STATS] LINKSLICES linkings= groups= parallel_linkings=
//! legacy_fallbacks=` (linking calls linked as slices, their groups, those linked on the pool,
//! and calls handed back to a cut, which the literal port never needs: always 0, see item 3).
//!
//! # Memory
//!
//! In the cut form, at the end of the kernels the linked groups (one Vec each, exactly sized) and
//! the new list of capacity n are alive together while the list is filled. Here the one list stays,
//! and the only allocation of the group step is the vector of `starts.len() + 1` slices.
//!
//! Measured once on 2026-10-13 (one binary; a loaded machine). massif, AP027131 x AP027133
//! (genetic code 4, one thread), heap requested: with the cut the peak is 107.3 MiB in the
//! first linking call (the new list 49.1 MiB and the linked groups 36.8 + 12.3 MiB); with the
//! slices the linking call reaches 77.3 MiB at most and the peak of the search is 102.4 MiB at
//! the end of the seed stage, as with the cut. Max RSS (medians of three alternating runs): that
//! pair 110 -> 111 MiB at one thread and 162 -> 163 at eight (its RSS peak is in the seed stage);
//! NZ_CP006932 against itself (516,395 HSPs in the first linking call) 142 -> 132 MiB at one
//! thread and 242 -> 222 at eight. `groups` of the first linking call 0.011 -> 0.002 s and
//! `assemble` 0.010 -> 0.000 s (the pair, one thread).

#[cfg(test)]
use super::x_hsp_inplace::x_assert_same_hits;
use super::x_hsp_inplace::{buffer_pools, group_starts, GroupArgs};
use super::{
    diagnostics_enabled, legacy_linking_requested, link_hsp_group_incremental_slice,
    link_hsp_group_ncbi, link_shadow_enabled, BufferPools, KarlinParams, LinkHspCutoffs,
    LinkingParams, QueryContext, UngappedHit,
};
// Diagnostic (LOSAT_TIMING): the clock of the linking call (tblastx/x_link_timing.rs); the
// linking reads nothing from it.
use crate::algorithm::tblastx::x_link_timing::{KernelTimes, LinkLap, LinkStage};
// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1290-1297
// As in `linking.rs`: a build without threads has no pool.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use std::cell::RefCell;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

// No NCBI counterpart: counters for LOSAT_X_STATS; they do not change any value NCBI computes.
static LINKINGS: AtomicU64 = AtomicU64::new(0);

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
thread_local! {
    // The kernel scratch of a pool thread (module doc, item 4): one per worker thread, made on
    // first use, shared by the pool path and the helped path.
    static POOLS: RefCell<Option<BufferPools>> = const { RefCell::new(None) };
}

// EXPERIMENT (LOSAT_X_TBXBOARD): whether this thread's linking calls are a TBXPIPE chain that
// must not wait in the pool (tblastx/x_batch_pipe.rs, item 6). No NCBI counterpart (scheduling
// only).
thread_local! {
    static HELPED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// The calling thread's helped mode until the guard drops (it restores the mode before).
pub(crate) struct XHelped(bool);

impl Drop for XHelped {
    fn drop(&mut self) {
        HELPED.with(|helped| helped.set(self.0));
    }
}

// EXPERIMENT (LOSAT_X_TBXBOARD): the linking calls and score sorts that follow on this thread
// are a chain of the batch pipe that must not wait in the pool: its groups are linked on the
// search's board (`x_link_board`), and its sorts are sorted there or stay on the thread
// (x_link_sort.rs `x_board_sort`, `pool_allowed`). No NCBI counterpart.
/// Sets this thread's helped mode to `on` until the guard drops.
pub(crate) fn x_helped(on: bool) -> XHelped {
    XHelped(HELPED.with(|helped| helped.replace(on)))
}

/// Whether this thread's linking calls are in the helped mode.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn x_helped_now() -> bool {
    HELPED.with(std::cell::Cell::get)
}
static GROUPS: AtomicU64 = AtomicU64::new(0);
static PARALLEL_LINKINGS: AtomicU64 = AtomicU64::new(0);

/// The summary line (printed by `main` at exit): the LOSAT_X_STATS counters.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_none() {
        return;
    }
    let load = |c: &AtomicU64| c.load(AtomicOrdering::Relaxed);
    // `legacy_fallbacks` stays 0: the literal port keeps the order and the number of the
    // records of its group (module doc, item 3), so no call goes back to a cut.
    eprintln!(
        "[X_STATS] LINKSLICES linkings={} groups={} parallel_linkings={} legacy_fallbacks=0",
        load(&LINKINGS),
        load(&GROUPS),
        load(&PARALLEL_LINKINGS)
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,553-558
// No computation: the arguments of the linking call that its group step reads, the values
// `link_in_place` builds `GroupArgs` from.
/// What `link_in_place` hands the group step besides the sorted list.
pub(super) struct Call<'a> {
    pub(super) params: &'a KarlinParams,
    pub(super) cutoffs: &'a LinkHspCutoffs,
    pub(super) linking_params: &'a LinkingParams,
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) subject_frame_bases: &'a [i32],
    pub(super) length_adj_per_context: &'a [i64],
    pub(super) eff_searchsp_per_context: &'a [i64],
    pub(super) log_k_by_ctx: &'a [f64],
    pub(super) allow_parallel: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:522-531
// NCBI's groups are sublists of the sorted array: `hp_frame_start[k]` and `hp_frame_number[k]`.
// The runs between the starts of `group_starts`, as views of the list: disjoint sub-slices of the
// one sorted list, not copies (module doc, item 1).
/// The groups of `hits` (maximal runs of one `translated_context_group` and one sign of
/// `s_frame`, whose starts after the first are `starts`) as disjoint sub-slices, in list order.
fn group_slices<'h>(hits: &'h mut [UngappedHit], starts: &[usize]) -> Vec<&'h mut [UngappedHit]> {
    if hits.is_empty() {
        return Vec::new();
    }
    let mut groups: Vec<&'h mut [UngappedHit]> = Vec::with_capacity(starts.len() + 1);
    let mut rest: &'h mut [UngappedHit] = hits;
    for &start in starts.iter().rev() {
        // `start >= 1` and below the length of what is left (the starts increase): the tail is
        // the run from `start` to the end of what is left.
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(start);
        groups.push(tail);
        rest = head;
    }
    // What is left is the first group.
    groups.push(rest);
    groups.reverse();
    groups
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// One group through the reference's per-group function with the reference's arguments, on its
// records where they lie (module doc, item 2): the reference body's per-group call, on a slice.
/// Links one group in place.
fn link_one_slice(group: &mut [UngappedHit], args: &GroupArgs<'_>, pool: &mut BufferPools) {
    if args.use_legacy {
        // The literal port takes and returns a Vec with the group's records in their order
        // (module doc, item 3): a copy in, the linked records back into their places.
        let linked = link_hsp_group_ncbi(
            group.to_vec(),
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
        );
        group.clone_from_slice(&linked);
    } else {
        link_hsp_group_incremental_slice(
            group,
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
        );
    }
}

/// The sizes the scratch is made with (`buffer_pools`).
#[derive(Clone, Copy)]
struct Scratch {
    initial_lh_size: usize,
    total_hits: usize,
    literal_port: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// The groups in frame_index order with one scratch, or on the pool with one scratch per worker
// thread, as the reference body does (module doc, items 2 and 4). Each kernel runs once, on its
// own slice, timed by the kernel clock.
/// Links every group; `true` when the pool ran them.
fn link_slices(
    groups: Vec<&mut [UngappedHit]>,
    args: &GroupArgs<'_>,
    scratch: Scratch,
    kernel_times: &KernelTimes,
    on_pool: bool,
) -> bool {
    // EXPERIMENT (LOSAT_X_TBXBOARD): the groups on the search's board (tblastx/x_work_board.rs),
    // taken before any unit, for a call on the pool and for a pipe chain (`x_helped_now`).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if (on_pool || x_helped_now())
        && groups.len() > 1
        && rayon::current_thread_index().is_some()
        && rayon::current_num_threads() > 1
    {
        if let Some((board, _)) = crate::algorithm::tblastx::x_work_board::current() {
            x_link_board(board, groups, args, &scratch, kernel_times);
            return on_pool;
        }
    }
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if on_pool {
        groups.into_par_iter().for_each(|group| {
            POOLS.with(|pools| {
                let mut pools_opt = pools.borrow_mut();
                let pool = pools_opt.get_or_insert_with(|| {
                    buffer_pools(
                        scratch.initial_lh_size,
                        scratch.total_hits,
                        scratch.literal_port,
                    )
                });
                kernel_times.time(|| link_one_slice(group, args, pool));
            })
        });
        return true;
    }
    let _ = on_pool;
    let mut pools = buffer_pools(
        scratch.initial_lh_size,
        scratch.total_hits,
        scratch.literal_port,
    );
    for group in groups {
        kernel_times.time(|| link_one_slice(group, args, &mut pools));
    }
    false
}

// EXPERIMENT (LOSAT_X_TBXBOARD): the groups of a linking call on the search's board, the largest
// first (tblastx/x_work_board.rs). No NCBI counterpart: the groups are disjoint and each kernel's
// result depends on its group alone (module doc, items 2-4), so which thread links which group,
// and when, changes nothing. No thread waits in the pool: the board's takers and this call's
// helper jobs take groups, and this thread blocks until the groups others took are linked (the
// kernels make no pool call).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_link_board(
    board: &crate::algorithm::tblastx::x_work_board::XBoard,
    mut groups: Vec<&mut [UngappedHit]>,
    args: &GroupArgs<'_>,
    scratch: &Scratch,
    kernel_times: &KernelTimes,
) {
    // The largest group first (equal sizes in group order): the call ends when its largest group
    // is linked.
    groups.sort_by_key(|group| std::cmp::Reverse(group.len()));
    let n = groups.len();
    let slices: Vec<(usize, usize)> = groups
        .iter_mut()
        .map(|group| (group.as_mut_ptr() as usize, group.len()))
        .collect();
    let run = |i: usize, _: &mut crate::algorithm::tblastx::x_work_board::XTaker| {
        let (ptr, len) = slices[i];
        // SAFETY: the slices are disjoint parts of the list this call borrows mutably (`groups`
        // came from `split_at_mut`); the board runs index `i` once, and `run_groups` returns only
        // after every index has run.
        let group = unsafe { std::slice::from_raw_parts_mut(ptr as *mut UngappedHit, len) };
        POOLS.with(|pools| {
            let mut pools_opt = pools.borrow_mut();
            let pool = pools_opt.get_or_insert_with(|| {
                buffer_pools(
                    scratch.initial_lh_size,
                    scratch.total_hits,
                    scratch.literal_port,
                )
            });
            kernel_times.time(|| link_one_slice(group, args, pool));
        });
    };
    board.run_groups(n, &run);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:522-531
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419
// No NCBI counterpart for the two switches: LOSAT_LINKING_LEGACY links every group with the literal
// port, LOSAT_LINK_FAST_SHADOW also runs it and compares (`link_hsp_group_incremental`). They do
// not change any value NCBI computes.
/// `link_groups_as`'s `use_legacy` and `fast_shadow` from the environment: the group step of
/// `x_hsp_inplace::link_in_place` for every list but a one-hit list with an HSP of length 0 or
/// less.
pub(super) fn link_modes() -> (bool, bool) {
    (legacy_linking_requested(), link_shadow_enabled())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// The group step of `x_hsp_inplace::link_in_place` on the sorted list (module doc), with the two
// variables of `link_groups_as` that `link_modes` gives (`link_modes` itself for the default list).
// Returns the list and its number of groups.
/// The sorted list with every group linked in place (the reference's list after its concatenation).
pub(super) fn link_groups<M>(
    hits: Vec<UngappedHit>,
    call: &Call<'_>,
    lap: &mut LinkLap,
    link_modes: M,
) -> (Vec<UngappedHit>, usize)
where
    M: FnOnce() -> (bool, bool),
{
    let (use_legacy, fast_shadow) = link_modes();
    link_groups_as(hits, call, lap, use_legacy, fast_shadow)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// `link_groups` with the two variables that choose the literal port and size the scratch given
// (`link_modes`; the tests run both ways, and a one-hit list with an HSP of length 0 or less asks
// for the literal port). The statements of the reference's group step in their order, with the
// cut replaced by the slices and the assembly by nothing.
/// `link_groups` with `use_legacy` and `fast_shadow` given.
fn link_groups_as(
    mut hits: Vec<UngappedHit>,
    call: &Call<'_>,
    lap: &mut LinkLap,
    use_legacy: bool,
    fast_shadow: bool,
) -> (Vec<UngappedHit>, usize) {
    // The list may arrive with spare capacity (it was grown by `extend`); give it back first, so
    // the slices below and the replay at the end work on a list of exactly its length.
    hits.shrink_to_fit();
    let total_hits = hits.len(); // Total HSP count (NCBI `hsp_list->hspcnt`)
    let starts = group_starts(&hits);
    let groups = group_slices(&mut hits, &starts);
    lap.lap(LinkStage::Groups);
    lap.sizes(total_hits, groups.iter().map(|group| group.len()));
    let kernel_times = lap.kernels();
    let group_count = groups.len();

    crate::utils::threading::report_stage(
        "tblastx",
        "linking",
        group_count,
        call.allow_parallel
            && group_count > 1
            && cfg!(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            )),
    );
    let initial_lh_size = (total_hits + 5).max(1024); // NCBI: the larger of 1024 and `hspcnt + 5`
    let literal_port = use_legacy || fast_shadow;
    let args = GroupArgs {
        use_legacy,
        diag_enabled: diagnostics_enabled(),
        params: call.params,
        cutoffs: call.cutoffs,
        linking_params: call.linking_params,
        query_contexts: call.query_contexts,
        subject_frame_bases: call.subject_frame_bases,
        length_adj_per_context: call.length_adj_per_context,
        eff_searchsp_per_context: call.eff_searchsp_per_context,
        log_k_by_ctx: call.log_k_by_ctx,
    };
    let scratch = Scratch {
        initial_lh_size,
        total_hits,
        literal_port,
    };
    let on_pool = link_slices(
        groups,
        &args,
        scratch,
        &kernel_times,
        call.allow_parallel && group_count > 1,
    );
    lap.lap(LinkStage::Kernels);
    // No assembly: the groups are in group order inside the list (module doc, item 3).
    lap.lap(LinkStage::Assemble);
    lap.add_kernels(&kernel_times);

    LINKINGS.fetch_add(1, AtomicOrdering::Relaxed);
    GROUPS.fetch_add(group_count as u64, AtomicOrdering::Relaxed);
    if on_pool {
        PARALLEL_LINKINGS.fetch_add(1, AtomicOrdering::Relaxed);
    }
    (hits, group_count)
}

#[cfg(test)]
#[path = "x_link_slices_tests.rs"]
mod tests;
