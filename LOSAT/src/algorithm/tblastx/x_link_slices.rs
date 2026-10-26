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
//! ```c
//! if (H->prev != NULL &&
//!     ((H->hsp->context/strand_factor) !=
//!      (H->prev->hsp->context/strand_factor) ||
//!      (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
//! { /* If frame switches, then start new list. */
//!    hp_frame_number[cur_frame]--;
//!    hp_frame_number[++cur_frame]++;
//!    hp_frame_start[cur_frame] = H;
//!    H->prev->next = NULL;
//!    H->prev = NULL;
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
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1082
//! ```c
//! for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
//!    hsp_list->hsp_array[index] = H->hsp;
//! ```
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
//!    `diagnostics_enabled()`, the same references, `MAX(1024, hspcnt+5)`), and the scratch with
//!    the same `buffer_pools`. The kernels: `link_one_slice` takes the branch of the reference's
//!    per-group call on `use_legacy`. The default branch calls
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
//!    `LOSAT_LINK_FAST_SHADOW`, and the literal port runs `pool_hsp_links.clear();
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
// NCBI reference: ncbi-blast/c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
// ```c
// CMTArgs(...)
// {
// #ifdef NCBI_NO_THREADS
//     m_NumThreads = CThreadable::kMinNumThreads;
//     m_MTMode = eNotSupported;
// #endif
// }
// ```
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
static HELPED_LINKINGS: AtomicU64 = AtomicU64::new(0);

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
thread_local! {
    // The kernel scratch of a pool thread (module doc, item 4): one per worker thread, made on
    // first use, shared by the pool path and the helped path.
    static POOLS: RefCell<Option<BufferPools>> = const { RefCell::new(None) };
}

// EXPERIMENT (LOSAT_X_TBXLINKHELP): whether this thread's linking calls are a TBXPIPE chain that
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

// EXPERIMENT (LOSAT_X_TBXLINKHELP): the linking calls and score sorts that follow on this thread
// are a chain of the batch pipe that must not wait in the pool: its groups are linked by
// the search's board, and its sorts stay on the thread (`x_link_board`). No NCBI
// counterpart.
/// Sets this thread's helped mode to `on` until the guard drops.
pub(crate) fn x_helped(on: bool) -> XHelped {
    XHelped(HELPED.with(|helped| helped.replace(on)))
}

/// Whether this thread's linking calls are in the helped mode.
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
// ```c
// s_BlastEvenGapLinkHSPs(EBlastProgramType program_number, BlastHSPList* hsp_list,
// ...
// for (frame_index=0; frame_index<num_query_frames; frame_index++)
// {
//    hp_start = s_LinkHSPStructReset(hp_start);
// ```
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
// ```c
// if (H->prev != NULL &&
//     ((H->hsp->context/strand_factor) !=
//      (H->prev->hsp->context/strand_factor) ||
//      (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
// { /* If frame switches, then start new list. */
//    hp_frame_number[cur_frame]--;
//    hp_frame_number[++cur_frame]++;
//    hp_frame_start[cur_frame] = H;
// ```
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
// ```c
// for (frame_index=0; frame_index<num_query_frames; frame_index++)
// {
//    hp_start = s_LinkHSPStructReset(hp_start);
//    hp_start->next = hp_frame_start[frame_index];
//    hp_frame_start[frame_index]->prev = hp_start;
//    number_of_hsps = hp_frame_number[frame_index];
// ```
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
// ```c
// for (frame_index=0; frame_index<num_query_frames; frame_index++)
// {
//    hp_start = s_LinkHSPStructReset(hp_start);
//    hp_start->next = hp_frame_start[frame_index];
//    hp_frame_start[frame_index]->prev = hp_start;
//    number_of_hsps = hp_frame_number[frame_index];
// ```
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
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
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
            HELPED_LINKINGS.fetch_add(1, AtomicOrdering::Relaxed);
            return on_pool;
        }
    }
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
                let n = group.len() as u32 | group.first().map_or(0, |h| h.s_idx << 24);
                let t0 = std::time::Instant::now();
                kernel_times.time(|| link_one_slice(group, args, pool));
                crate::utils::xtline::rec("kernel", n, t0, std::time::Instant::now());
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
        let n = group.len() as u32 | group.first().map_or(0, |h| h.s_idx << 24);
        let t0 = std::time::Instant::now();
        kernel_times.time(|| link_one_slice(group, args, &mut pools));
        crate::utils::xtline::rec("kernel", n, t0, std::time::Instant::now());
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
        let tag = group.len() as u32 | group.first().map_or(0, |h| h.s_idx << 24);
        let t0 = std::time::Instant::now();
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
        crate::utils::xtline::rec("kernel", tag, t0, std::time::Instant::now());
    };
    board.run_groups(n, &run, rayon::current_num_threads() - 1);
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
// The group step of `link_in_place` on the sorted list (module doc). Returns
// the list and its number of groups.
/// The sorted list with every group linked in place (the reference's list after its concatenation).
pub(super) fn link_groups(
    hits: Vec<UngappedHit>,
    call: &Call<'_>,
    lap: &mut LinkLap,
) -> (Vec<UngappedHit>, usize) {
    link_groups_as(
        hits,
        call,
        lap,
        legacy_linking_requested(),
        link_shadow_enabled(),
    )
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
// `link_groups` with the two variables that choose the literal port and size the scratch given
// (the tests run both ways). The statements of the reference's group step in their order, with the
// cut replaced by the slices and the assembly by nothing.
/// `link_groups` with `legacy_linking_requested()` and `link_shadow_enabled()` given.
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
    let total_hits = hits.len(); // NCBI: hsp_list->hspcnt
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
    let initial_lh_size = (total_hits + 5).max(1024); // NCBI: MAX(1024, hspcnt+5)
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
mod tests {
    use super::super::tests::{fast_test_contexts, fast_test_cutoffs, mock_hit};
    use super::super::{
        apply_sum_stats_even_gap_linking_with_parallel, fwd_compare_hsps_transl,
        link_hsp_group_incremental, replay_ncbi_output_list, rev_compare_hsps_tbx,
        rev_compare_hsps_transl, sort_hsps_by_ncbi_link_order, sort_linked_hsp_references,
        translated_context_group,
    };
    use super::*;

    /// xorshift64, for reproducible random lists (as in `x_hsp_inplace.rs`).
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
    }

    const FRAMES: [i8; 6] = [1, 2, 3, -1, -2, -3];

    /// A random HSP of one of `queries` queries (the generator of `x_hsp_inplace.rs`'s tests):
    /// few distinct coordinates and scores when `span` is small, every record distinct in
    /// `hsp_list_order` and in the bits of `e_value`.
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

    /// The list as the group step of a linking call receives it: link ids, no chain link and
    /// `num = 1` set as the call does, then sorted by `s_RevCompareHSPsTbx`.
    fn sorted(mut hits: Vec<UngappedHit>) -> Vec<UngappedHit> {
        for (link_id, hit) in hits.iter_mut().enumerate() {
            hit.link_id = u32::try_from(link_id).unwrap();
            hit.chain_next_link_id = None;
            hit.num = 1;
        }
        sort_hsps_by_ncbi_link_order(&mut hits, rev_compare_hsps_tbx);
        hits
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

    /// The group step of the reference body (`linking.rs`): the push grouping, each group
    /// through the per-group function the body chooses with one scratch sized as there, and the
    /// linked groups concatenated in group order.
    fn reference_groups(
        sorted: &[UngappedHit],
        call: &Call<'_>,
        use_legacy: bool,
    ) -> Vec<UngappedHit> {
        let total_hits = sorted.len();
        let mut pools = BufferPools {
            lh_helpers: Vec::with_capacity((total_hits + 5).max(1024)),
            hsp_links: Vec::with_capacity(total_hits),
        };
        let mut ordered = Vec::new();
        for group_hits in push_grouping(sorted.to_vec()) {
            let processed = if use_legacy {
                link_hsp_group_ncbi(
                    group_hits,
                    call.params,
                    call.cutoffs,
                    call.linking_params.gap_decay_rate,
                    false,
                    call.linking_params.subject_len_nucl,
                    call.query_contexts,
                    call.subject_frame_bases,
                    call.length_adj_per_context,
                    call.eff_searchsp_per_context,
                    call.log_k_by_ctx,
                    &mut pools.lh_helpers,
                    &mut pools.hsp_links,
                )
            } else {
                link_hsp_group_incremental(
                    group_hits,
                    call.params,
                    call.cutoffs,
                    call.linking_params.gap_decay_rate,
                    call.linking_params.subject_len_nucl,
                    call.query_contexts,
                    call.subject_frame_bases,
                    call.length_adj_per_context,
                    call.eff_searchsp_per_context,
                    call.log_k_by_ctx,
                    &mut pools.lh_helpers,
                    &mut pools.hsp_links,
                )
            };
            ordered.extend(processed);
        }
        ordered
    }

    /// The end of the linking call: its two final sorts and its replay (the reference's).
    fn reference_replay(mut results: Vec<UngappedHit>) -> Vec<UngappedHit> {
        if results.len() <= 1 {
            return results;
        }
        let mut references: Vec<&mut UngappedHit> = results.iter_mut().collect();
        sort_linked_hsp_references(&mut references, rev_compare_hsps_transl);
        sort_linked_hsp_references(&mut references, fwd_compare_hsps_transl);
        replay_ncbi_output_list(references)
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

    // No NCBI counterpart: `body` once on the test thread and once inside pools of 1, 2 and 8
    // threads (a build without the parallel feature has no pool and runs it on the test thread
    // only), as in `x_link_sort.rs`'s tests. The flag tells `body` whether it runs in a pool.
    fn on_and_off_pools(body: impl Fn(&str, bool) + Sync) {
        body("no pool", false);
        #[cfg(feature = "parallel")]
        for threads in [1usize, 2, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("test pool");
            pool.install(|| body(&format!("pool of {threads}"), true));
        }
    }

    /// Does this build link groups on the pool?
    const PARALLEL_BUILD: bool = cfg!(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ));

    /// The cutoff sets of the kernel fuzz test (as in `x_hsp_inplace.rs`'s tests).
    fn cutoff_sets() -> [LinkHspCutoffs; 3] {
        [
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
        ]
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
    // ```
    // Module doc, item 1: the slices are the reference's groups (the push loop), the same records
    // in the same order, and each one is a view of the list at its run's position: no HSP, one
    // HSP, one group, sorted lists of one and three queries, and groups of one HSP each.
    #[test]
    fn the_slices_are_the_groups_of_the_reference() {
        fn check(mut hits: Vec<UngappedHit>, label: &str) -> usize {
            let expected = push_grouping(hits.clone());
            let starts = group_starts(&hits);
            let base = hits.as_ptr() as usize;
            let size = std::mem::size_of::<UngappedHit>();
            let groups = group_slices(&mut hits, &starts);
            assert_eq!(groups.len(), expected.len(), "{label}: number of groups");
            let mut position = 0usize;
            for (index, (slice, group)) in groups.iter().zip(&expected).enumerate() {
                assert_eq!(
                    slice.as_ptr() as usize,
                    base + position * size,
                    "{label}: group {index} lies in the list at its run's position"
                );
                position += slice.len();
                assert_same_list(&format!("{label}: group {index}"), slice, group);
            }
            groups.len()
        }

        let mut rng = Rng(0xD1B5_4A32_D192_ED03);
        assert_eq!(check(Vec::new(), "no HSP"), 0);
        assert_eq!(check(random_hits(&mut rng, 1, 50, 1), "one HSP"), 1);
        let mut one_group = random_hits(&mut rng, 300, 50, 1);
        for hit in one_group.iter_mut() {
            hit.q_frame = 2;
            hit.s_frame = -1;
        }
        assert_eq!(check(one_group, "one group"), 1);
        for queries in [1usize, 3] {
            for n in [2usize, 5, 64, 3000] {
                let hits = sorted(random_hits(&mut rng, n, 60, queries));
                let groups = check(hits, &format!("sorted, {queries} queries, n {n}"));
                assert!(groups <= 4 * queries);
                if n == 3000 {
                    assert_eq!(groups, 4 * queries);
                }
            }
        }
        let mut singles = random_hits(&mut rng, 97, 50, 1);
        for (index, hit) in singles.iter_mut().enumerate() {
            hit.q_idx = (index / 2) as u32;
            hit.q_frame = 1;
            hit.s_frame = if index % 2 == 0 { 1 } else { -1 };
        }
        assert_eq!(check(singles, "groups of one HSP"), 97);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558,1080-1082
    // ```c
    // for (frame_index=0; frame_index<num_query_frames; frame_index++)
    // {
    //    hp_start = s_LinkHSPStructReset(hp_start);
    //    hp_start->next = hp_frame_start[frame_index];
    // ...
    // for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
    //    hsp_list->hsp_array[index] = H->hsp;
    // ```
    // Module doc, items 1-5: on random lists of one and two queries (up to eight groups, chains,
    // comparator-equal HSPs) the slices give the reference's group step, every field, with the
    // incremental kernel and with the literal port, serially and on pools of 1, 2 and 8 threads,
    // in the list's own allocation; after the final sorts and the replay the result is the whole
    // linking call's.
    #[test]
    fn the_slices_give_the_reference_list_serially_and_on_pools() {
        let (params, contexts) = two_query_contexts();
        let log_k: Vec<f64> = contexts.iter().map(|c| c.karlin_params.k.ln()).collect();
        let cutoffs = cutoff_sets();
        let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(120);
        let mut rng = Rng(0xA076_1D64_78BD_642F);
        let (mut most_groups, mut chained, mut total) = (0usize, 0usize, 0usize);
        for case in 0..cases {
            let n = match case % 10 {
                0 => case / 10 % 4,
                1 => 2500,
                _ => 1 + rng.below(400) as usize,
            };
            let queries = 1 + case % 2;
            let span = [60u64, 400, 3000, 20_000][case / 2 % 4];
            let hits = random_hits(&mut rng, n, span, queries);
            let linking = LinkingParams {
                subject_len_nucl: 300_000,
                gap_decay_rate: 0.5,
                cutoffs: cutoffs[case % cutoffs.len()],
            };
            let list = sorted(hits.clone());
            let expected_groups = push_grouping(list.clone()).len();
            most_groups = most_groups.max(expected_groups);
            for use_legacy in [false, true] {
                let call = |allow_parallel| Call {
                    params: &params,
                    cutoffs: &linking.cutoffs,
                    linking_params: &linking,
                    query_contexts: &contexts,
                    subject_frame_bases: &[0; 6],
                    length_adj_per_context: &[49; 12],
                    eff_searchsp_per_context: &[9_372_428_362; 12],
                    log_k_by_ctx: &log_k,
                    allow_parallel,
                };
                let expected = reference_groups(&list, &call(false), use_legacy);
                if !use_legacy {
                    total += expected.len();
                    chained += expected.iter().filter(|hit| hit.linked_set).count();
                }
                on_and_off_pools(|where_, _| {
                    for allow_parallel in [false, true] {
                        let label = format!(
                            "case {case} n {n} queries {queries} span {span} legacy {use_legacy} {where_} parallel {allow_parallel}"
                        );
                        let input = list.clone();
                        let allocation = input.as_ptr();
                        let before = PARALLEL_LINKINGS.load(AtomicOrdering::Relaxed);
                        let (actual, groups) = link_groups_as(
                            input,
                            &call(allow_parallel),
                            &mut LinkLap::off(),
                            use_legacy,
                            use_legacy,
                        );
                        assert_eq!(actual.as_ptr(), allocation, "{label}: no second list");
                        assert_eq!(groups, expected_groups, "{label}");
                        assert_same_list(&label, &actual, &expected);
                        if PARALLEL_BUILD && allow_parallel && groups > 1 {
                            assert!(
                                PARALLEL_LINKINGS.load(AtomicOrdering::Relaxed) > before,
                                "{label}: the groups ran on the pool"
                            );
                        }
                    }
                });
            }
            // The whole linking call: the reference body against the slices
            // followed by the reference's final sorts and replay.
            for allow_parallel in [false, true] {
                let call = Call {
                    params: &params,
                    cutoffs: &linking.cutoffs,
                    linking_params: &linking,
                    query_contexts: &contexts,
                    subject_frame_bases: &[0; 6],
                    length_adj_per_context: &[49; 12],
                    eff_searchsp_per_context: &[9_372_428_362; 12],
                    log_k_by_ctx: &log_k,
                    allow_parallel,
                };
                let whole_call = apply_sum_stats_even_gap_linking_with_parallel(
                    hits.clone(),
                    &params,
                    &linking,
                    &contexts,
                    &[0; 6],
                    &[49; 12],
                    &[9_372_428_362; 12],
                    allow_parallel,
                );
                let (linked, _) =
                    link_groups_as(list.clone(), &call, &mut LinkLap::off(), false, false);
                assert_same_list(
                    &format!("case {case} whole call, parallel {allow_parallel}"),
                    &reference_replay(linked),
                    &whole_call,
                );
            }
        }
        if cases >= 120 {
            // The lists must exercise chains and every group of two queries.
            assert!(chained * 10 > total, "{chained} of {total}");
            assert_eq!(most_groups, 8);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
    // ```c
    // for (frame_index=0; frame_index<num_query_frames; frame_index++)
    // {
    //    hp_start = s_LinkHSPStructReset(hp_start);
    //    hp_start->next = hp_frame_start[frame_index];
    // ```
    // The edges: an empty list (no group, no kernel), one HSP, one group of many HSPs with the
    // pool allowed, and a list that arrives with spare capacity (given back, as the reference's cut
    // does).
    #[test]
    fn empty_single_and_one_group_lists() {
        let (params, contexts) = two_query_contexts();
        let log_k: Vec<f64> = contexts.iter().map(|c| c.karlin_params.k.ln()).collect();
        let linking = LinkingParams {
            subject_len_nucl: 300_000,
            gap_decay_rate: 0.5,
            cutoffs: fast_test_cutoffs(),
        };
        let call = Call {
            params: &params,
            cutoffs: &linking.cutoffs,
            linking_params: &linking,
            query_contexts: &contexts,
            subject_frame_bases: &[0; 6],
            length_adj_per_context: &[49; 12],
            eff_searchsp_per_context: &[9_372_428_362; 12],
            log_k_by_ctx: &log_k,
            allow_parallel: true,
        };
        let (linked, groups) = link_groups_as(Vec::new(), &call, &mut LinkLap::off(), false, false);
        assert!(linked.is_empty());
        assert_eq!(groups, 0);

        let mut rng = Rng(0x5851_F42D_4C95_7F2D);
        let one = sorted(random_hits(&mut rng, 1, 50, 1));
        let (linked, groups) =
            link_groups_as(one.clone(), &call, &mut LinkLap::off(), false, false);
        assert_eq!(groups, 1);
        assert_same_list("one HSP", &linked, &reference_groups(&one, &call, false));

        let mut one_group = random_hits(&mut rng, 500, 300, 1);
        for hit in one_group.iter_mut() {
            hit.q_frame = -3;
            hit.ctx_idx = 5;
            hit.s_frame = 2;
        }
        let one_group = sorted(one_group);
        let (linked, groups) =
            link_groups_as(one_group.clone(), &call, &mut LinkLap::off(), false, false);
        assert_eq!(groups, 1);
        assert_same_list(
            "one group",
            &linked,
            &reference_groups(&one_group, &call, false),
        );
        assert!(
            linked.iter().any(|hit| hit.linked_set),
            "the group has chains"
        );

        let mut spare = sorted(random_hits(&mut rng, 300, 400, 2));
        let expected = reference_groups(&spare, &call, false);
        spare.reserve(1000);
        assert!(spare.capacity() > spare.len());
        let (linked, _) = link_groups_as(spare, &call, &mut LinkLap::off(), false, false);
        assert_same_list("spare capacity", &linked, &expected);
    }
}
