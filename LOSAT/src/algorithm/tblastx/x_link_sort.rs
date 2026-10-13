//! EXPERIMENT (LOSAT_X_LINKSORT / LOSAT_X_LINKSORTSHADOW): the index sorts of the TBLASTX HSP
//! list (`x_hsp_inplace.rs`). With the switch on, the linking call sorts its output order once
//! where it sorts twice, and an index sort of `PARALLEL_MIN` indices or more runs as a stable
//! parallel sort when its call may use the search pool. Every sort returns the index order of
//! today's statements, element for element, at any number of threads.
//!
//! This file is a child module of `sum_stats_linking/linking.rs` (declared there with `#[path]`,
//! beside `x_hsp_inplace.rs`), so it names that file's comparators without changing them.
//! `x_hsp_inplace.rs` calls it at its two index sort sites (`sort_in_place`,
//! `replay_in_place_counted`); nothing else does.
//!
//! The switch is not in the strict set (`utils/xswitch.rs`): it is on when `LOSAT_X_LINKSORT` is
//! set, to any value, and off otherwise, and `LOSAT_X_DEFAULT` does not touch it. With it off the
//! statements of `x_hsp_inplace.rs` run as before, after one read of `mode()` per index sort. It
//! acts only where `LOSAT_X_HSPINPLACE` does: the reference bodies (`LOSAT_X_HSPINPLACE=0`,
//! `LOSAT_X_DEFAULT=0`) have their own sorts and are not touched.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
//! ```c
//! /* Sort by (reverse) position. */
//! if (kTranslatedQuery) {
//!    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//!          s_RevCompareHSPsTbx);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
//! ```c
//! if (kTranslatedQuery) {
//!    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//!          s_RevCompareHSPsTransl);
//!    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//!          s_FwdCompareHSPsTransl);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
//! ```c
//!         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//!               ScoreCompareHSPs);
//! ```
//! NCBI sorts with the C library's `qsort`, which is not stable: where its comparator returns 0
//! the order is the library's. LOSAT's defined behaviour is the stable sorts of its reference
//! bodies (`ungapped_hit_sort.rs`, `sort_linked_hsp_references` in `linking.rs`), which the
//! fixtures compare with NCBI's reports. This file reproduces the index order of those stable
//! sorts and changes no comparator. No NCBI counterpart for a parallel sort or for one sort in
//! place of two: both are other ways to compute the same order.
//!
//! # Why the index order is today's
//!
//! Terms. A comparator is a *total preorder* when, for all records a, b, c, `compare(a, b)` is
//! the reverse of `compare(b, a)` (so `compare(a, a)` is `Equal`), and `compare(a, b) != Greater`
//! with `compare(b, c) != Greater` gives `compare(a, c) != Greater`, `Equal` only when both are.
//! A comparator that takes a tuple of integers from each record and compares the two tuples
//! component by component with `Ord::cmp`, each component in a fixed direction, is one: it is
//! the lexicographic order of integer tuples, a total order, applied to the keys.
//!
//! (S) The result of a stable sort is unique. Let the comparator be a total preorder and the
//! input the sequence x_0 .. x_{n-1}. A sorted output never has x_j before x_i when
//! `compare(x_i, x_j)` is `Less`, and a stable sort keeps x_i before x_j when the two are
//! `Equal` and i < j. So in the output x_i is before x_j exactly when `compare(x_i, x_j)` is
//! `Less`, or `Equal` with i < j. That fixes the relative order of every pair, hence the whole
//! sequence. The two sorts used here are documented as stable:
//!
//! - `slice::sort_by` (Rust 1.92, alloc/src/slice.rs:136-144): "Sorts the slice in ascending
//!   order with a comparison function, preserving initial order of equal elements. This sort is
//!   stable (i.e., does not reorder equal elements) and *O*(*n* \* log(*n*)) worst-case. If the
//!   comparison function `compare` does not implement a [total order], the function may panic;
//!   even if the function exits normally, the resulting order of elements in the slice is
//!   unspecified."
//! - `ParallelSliceMut::par_sort_by` (rayon 1.12.0, src/slice/mod.rs:399-404): "Sorts the slice
//!   in parallel with a comparator function. This sort is stable (i.e., does not reorder equal
//!   elements) and *O*(*n* \* log(*n*)) worst-case. The comparator function must define a total
//!   ordering for the elements in the slice. If the ordering is not total, the order of the
//!   elements is unspecified."
//!
//! Both documents say "total order" of a comparator that may return `Equal` for different
//! elements, which is what "does not reorder equal elements" is about: the total preorder above.
//! The elements sorted here are indices, all different, and the records they name are not
//! written while a sort runs.
//!
//! 1. One final sort instead of two (`final_sort`). Today `replay_in_place_counted` sorts
//!    `order = 0..n` stably with B = `rev_compare_hsps_transl`, then stably with A =
//!    `fwd_compare_hsps_transl`. By (S), after the first sort i is before j exactly when B is
//!    `Less`, or `Equal` with i < j. By (S) again, after the second sort i is before j exactly
//!    when A is `Less`, or A is `Equal` and i was before j: that is, when A is `Less`, or A is
//!    `Equal` and B is `Less`, or both are `Equal` and i < j. The composite
//!    `A(a, b).then_with(|| B(a, b))` is `Less` exactly in the first two cases and `Equal` exactly
//!    when both are (`Ordering::then_with` returns `self` when it is not `Equal`, otherwise the
//!    result of the function), and it is a total preorder (the key of A followed by the key of
//!    B). By (S) one stable sort of `0..n` with it puts i before j in the same three cases: the
//!    same sequence. For these two comparators A is `Equal` exactly when B is, because both read
//!    the keys (query strand group, `q_aa_start`, `s_aa_start`); the argument does not use that,
//!    and the unit tests also run the composite on pairs of comparators without that property.
//! 2. Large index sorts on the pool (`sort_stable`). The same comparator on the same indices and
//!    records; by (S) `par_sort_by` returns the sequence `sort_by` returns, and its
//!    documentation makes the result depend on nothing else, so not on the number of threads.
//!    (rayon 1.12.0 sorts fixed chunks of 2,000 elements, reverses a chunk only when it is
//!    strictly descending, and merges: src/slice/sort.rs:1510-1608.)
//! 3. The comparators these sorts receive, each a lexicographic comparison of integer fields of
//!    `UngappedHit` (`chaining.rs`: `q_idx: u32`, `q_frame`/`s_frame: i8`, the four offsets
//!    `usize`, `raw_score: i32`); none reads `e_value` or any other float:
//!    - `rev_compare_hsps_tbx` (`linking.rs`; first sort of the linking call, and the sort of the
//!      reference linking body when it reaches `sort_dispatch`): `q_idx * 2 + (q_frame < 0)`
//!      ascending, `s_frame.signum()` ascending, then `q_aa_start`, `q_aa_end`, `s_aa_start`,
//!      `s_aa_end` descending.
//!    - `rev_compare_hsps_transl` (`linking.rs`; B): the group ascending, `q_aa_start` and
//!      `s_aa_start` descending.
//!    - `fwd_compare_hsps_transl` (`linking.rs`; A): the group, `q_aa_start` and `s_aa_start`
//!      ascending.
//!    - `score_compare_ungapped_hits_ncbi` (`blast_gapalign.rs`; every score sort, the other
//!      caller of `sort_dispatch`): `raw_score` descending, `s_aa_start` ascending, `s_aa_end`
//!      descending, `q_aa_start` ascending, `q_aa_end` descending.
//!
//!    `sort_dispatch` takes its comparator as an argument. A comparator that is not a total
//!    preorder already has an unspecified result, and may panic, in today's `sort_by` (quoted
//!    above), so this path asks nothing of a comparator that today's statement does not.
//!
//! # The pool and `-num_threads`
//!
//! A serial build (no `parallel` feature, wasm32 without `wasm-threads`) has no parallel sort:
//! the `cfg` of `sort_stable` is the `cfg` of the per-group linking in `link_in_place`.
//!
//! The two sorts of a linking call (`Site::First`, `Site::Final`) use the pool when the call's
//! `allow_parallel` is true: the flag, and the `cfg`, under which the per-group linking of the
//! same call uses it. The flag is `use_parallel` of the batch function (`run_impl.rs`), which is
//! `SearchPool::enabled()`: the search has a pool, so `-num_threads` is 2 or more.
//!
//! No flag reaches `sort_dispatch` (`Site::Dispatch`). It uses the pool when the calling thread
//! is a worker of a rayon pool, `rayon::current_thread_index().is_some()`:
//!
//! - `-num_threads` N > 1. `with_search_pool` (`utils/threading.rs`) builds a pool of N threads,
//!   checks `pool.current_num_threads() == requested`, and runs the whole search inside
//!   `ThreadPool::install`, on the pool's slot zero. rayon-core 1.13.0, `install`: "Executes `op`
//!   within the thread pool. Any attempts to use `join`, `scope`, or parallel iterators will then
//!   operate within that thread pool." `par_sort_by` is made of a parallel iterator over the
//!   chunks and of `join`, so every part of the sort runs on those N threads, whether the caller
//!   is slot zero (one subject) or a worker that runs a subject (several subjects). No other
//!   pool exists outside tests, and `with_search_pool` refuses to start on a rayon worker.
//! - `-num_threads 1`. `with_search_pool` builds no pool and runs the search on the thread that
//!   called it. The searches this repository starts are called from `main` and from the Web
//!   adapter (a crate without rayon), neither a rayon worker. `current_thread_index`: "If called
//!   from a Rayon worker thread, returns the index of that thread within its current pool; if
//!   not called from a Rayon thread, returns `None`." So the sort is the serial statement.
//!   `rayon::current_num_threads()` is not the test: off a pool it is "the number of threads for
//!   the global thread pool", and asking starts that pool.
//! - Not shown by this repository's code: a program that uses LOSAT as a library and calls a
//!   `-num_threads 1` search from a thread of its own rayon pool. There the test is true and a
//!   large dispatch sort runs on that program's pool. `pool_threads` of blastn
//!   (`x_mb_lookup_par.rs`, `x_pair_seed.rs`) reads the same property of the thread.
//!
//! A parallel sort waits for its parts, and rayon may run another job of the pool on the
//! waiting thread: another subject of a search with several subjects, a group of another
//! linking call. The per-group linking and the unit path already wait this way. The sort sites
//! hold nothing such a job takes: the thread-local flags of `x_hsp_inplace.rs` are taken before
//! a sort starts, the timer's flag (`x_link_timing.rs`) at the start of the linking call it
//! times, the scratch `POOLS` is borrowed only inside a per-group call, which sorts through
//! neither site and uses no pool, and no lock is held at a sort site (the files that sort have
//! no lock and no other thread-local). In a search with several subjects the `sort` and
//! `final_sorts` stages of `LOSAT_TIMING`, and `reap` and `merge` of `x_subject`, can therefore
//! hold time of another subject, as `kernels` can.
//!
//! `PARALLEL_MIN` is 32,768. Measured on 2026-10-13 at 8 threads (medians of three alternating
//! runs without and with the switch, a loaded machine): MelaMJNV x PemoMJNVA, whose four longest
//! sorts are of about 34,400 indices, `sort` 0.004 s both ways, so the threshold is not clearly
//! wrong; AP027131 x AP027133 (338,859 HSPs) `sort` 0.085 to 0.039 s, `final_sorts` 0.067 to
//! 0.009 s, `reap` 0.073 to 0.045 s; at one thread `final_sorts` 0.060 to 0.035 s (one sort
//! instead of two; the other sorts are today's statements there). What is left of `sort` and
//! `reap` includes the move of the records, which stays on one thread.
//!
//! # Shadow and counters
//!
//! `LOSAT_X_LINKSORTSHADOW=1` turns the path on and runs every index sort both ways: this
//! file's sort on a copy of the index vector and today's statement or statements on the vector
//! itself. The two must be equal element for element; the first difference stops the run with
//! the switch, the sort (first, final or dispatch), the length and the position. The search goes
//! on with today's result. `main` prints `[X_SHADOW] LINKSORT sorts= indices_compared=` at exit,
//! and with `LOSAT_X_STATS` (switch or shadow on)
//! `[X_STATS] LINKSORT sorts= indices= parallel_sorts= parallel_indices= composite_sorts=`.

use super::{fwd_compare_hsps_transl, rev_compare_hsps_transl, UngappedHit};
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
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use std::sync::OnceLock;

// No NCBI counterpart: counters for LOSAT_X_STATS and the shadow summary; they do not change any
// value NCBI computes.
static SORTS: AtomicU64 = AtomicU64::new(0);
static INDICES: AtomicU64 = AtomicU64::new(0);
static PARALLEL_SORTS: AtomicU64 = AtomicU64::new(0);
static PARALLEL_INDICES: AtomicU64 = AtomicU64::new(0);
static COMPOSITE_SORTS: AtomicU64 = AtomicU64::new(0);
static SHADOW_SORTS: AtomicU64 = AtomicU64::new(0);
static SHADOW_INDICES: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: NCBI's qsort calls run on the thread of their search. Below this length an
// index sort stays the serial statement: rayon sorts chunks of 2,000 elements as separate jobs
// and merges them, which costs more than it saves on a short list.
/// The smallest index sort that goes to the pool.
pub(super) const PARALLEL_MIN: usize = 32_768;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486,990-994
// ```c
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ...
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTransl);
//    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_FwdCompareHSPsTransl);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
// ```c
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
// ```
// No computation: names the index sort for the shadow's message and carries what its caller
// allows.
/// Which index sort of `x_hsp_inplace.rs`, with the linking call's `allow_parallel` where there
/// is one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Site {
    /// `link_in_place`: the sort before the groups (`s_RevCompareHSPsTbx`).
    First(bool),
    /// `replay_in_place_counted`: the sorts after the groups (`s_RevCompareHSPsTransl`, then
    /// `s_FwdCompareHSPsTransl`).
    Final(bool),
    /// `sort_dispatch`: every other HSP sort (the score sorts, and the first sort of a
    /// reference linking body).
    Dispatch,
}

impl Site {
    fn name(self) -> &'static str {
        match self {
            Site::First(_) => "first",
            Site::Final(_) => "final",
            Site::Dispatch => "dispatch",
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
// ```c
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ```
// Dispatch switch, read once: 0 = today's statements, 1 = this path, 2 = both, compared. Not in
// `xswitch::STRICT_SET`, so not read with `xswitch::on`.
/// 0 = off, 1 = LOSAT_X_LINKSORT, 2 = LOSAT_X_LINKSORTSHADOW.
pub(super) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_LINKSORTSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_LINKSORT").is_some() {
            1
        } else {
            0
        }
    })
}

/// The summary lines (printed by `main` at exit): the shadow summary and the LOSAT_X_STATS
/// counters, only with the switch or its shadow on.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    let m = mode();
    if m == 0 {
        return;
    }
    let load = |c: &AtomicU64| c.load(AtomicOrdering::Relaxed);
    if m == 2 {
        eprintln!(
            "[X_SHADOW] LINKSORT sorts={} indices_compared={}",
            load(&SHADOW_SORTS),
            load(&SHADOW_INDICES)
        );
    }
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] LINKSORT sorts={} indices={} parallel_sorts={} parallel_indices={} composite_sorts={}",
            load(&SORTS),
            load(&INDICES),
            load(&PARALLEL_SORTS),
            load(&PARALLEL_INDICES),
            load(&COMPOSITE_SORTS)
        );
    }
}

// No NCBI counterpart: counts one index sort of this path for LOSAT_X_STATS.
fn count(indices: usize, parallel: bool, composite: bool) {
    SORTS.fetch_add(1, AtomicOrdering::Relaxed);
    INDICES.fetch_add(indices as u64, AtomicOrdering::Relaxed);
    if parallel {
        PARALLEL_SORTS.fetch_add(1, AtomicOrdering::Relaxed);
        PARALLEL_INDICES.fetch_add(indices as u64, AtomicOrdering::Relaxed);
    }
    if composite {
        COMPOSITE_SORTS.fetch_add(1, AtomicOrdering::Relaxed);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149
// ```c
// TBlastThreads the_threads(GetNumberOfThreads());
// ```
// Scheduling only: whether an index sort may run on the search pool (module doc, "The pool and
// `-num_threads`"). A linking call's sorts follow its `allow_parallel`, as its per-group linking
// does; a sort that no flag reaches follows the calling thread, which is a rayon worker exactly
// when the search runs inside its pool. `rayon::current_thread_index` reads a thread-local and
// starts no pool.
/// May the sort of `site` use the pool?
fn pool_allowed(site: Site) -> bool {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    {
        match site {
            Site::First(allow_parallel) | Site::Final(allow_parallel) => allow_parallel,
            Site::Dispatch => rayon::current_thread_index().is_some(),
        }
    }
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    {
        let _ = site;
        false
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
// ```c
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
// ```
// The statement of `x_hsp_inplace.rs` (and of the reference bodies): a stable sort of the indices
// with the comparator applied to the records they name.
/// Today's index sort.
fn sort_today<C>(indices: &mut [usize], records: &[UngappedHit], compare: C)
where
    C: Fn(&UngappedHit, &UngappedHit) -> Ordering,
{
    indices.sort_by(|&lhs, &rhs| compare(&records[lhs], &records[rhs]));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
// ```c
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ```
// The same stable sort of the same indices with the same comparator: rayon's stable parallel sort
// when the call may use the pool and the list is long, today's statement otherwise. For a total
// preorder both return the one stable order (module doc, (S) and item 2).
/// The index order of `sort_today`; `true` when the sort ran on the pool.
fn sort_stable<C>(
    indices: &mut [usize],
    records: &[UngappedHit],
    compare: C,
    on_pool: bool,
    parallel_min: usize,
) -> bool
where
    C: Fn(&UngappedHit, &UngappedHit) -> Ordering + Sync,
{
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if on_pool && indices.len() >= parallel_min {
        indices.par_sort_by(|&lhs, &rhs| compare(&records[lhs], &records[rhs]));
        return true;
    }
    let _ = (on_pool, parallel_min);
    sort_today(indices, records, compare);
    false
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
// ```c
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ```
// No computation: the shadow compares two index vectors and stops at the first difference.
/// Shadow check: `new` (this path) against `today`, element for element.
fn shadow_check(site: Site, new: &[usize], today: &[usize]) {
    let differs = new
        .iter()
        .zip(today)
        .position(|(lhs, rhs)| lhs != rhs)
        .or_else(|| (new.len() != today.len()).then(|| new.len().min(today.len())));
    if let Some(position) = differs {
        panic!(
            "LOSAT_X_LINKSORTSHADOW: the {} index sort of {} indices differs at position {position}: this path {:?}, today {:?}",
            site.name(),
            today.len(),
            new.get(position),
            today.get(position)
        );
    }
    SHADOW_SORTS.fetch_add(1, AtomicOrdering::Relaxed);
    SHADOW_INDICES.fetch_add(today.len() as u64, AtomicOrdering::Relaxed);
}

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
// The index sort of `x_hsp_inplace::sort_in_place` with the switch or its shadow on.
/// Sorts `indices` into the order of today's `sort_by` with `compare` over `records`.
pub(super) fn sort_indices(
    indices: &mut [usize],
    records: &[UngappedHit],
    compare: fn(&UngappedHit, &UngappedHit) -> Ordering,
    site: Site,
) {
    sort_indices_as(
        indices,
        records,
        compare,
        site,
        mode() == 2,
        pool_allowed(site),
        PARALLEL_MIN,
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
// ```c
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ```
// The body of `sort_indices` with its three decisions as arguments (the unit tests set them).
// With `shadow`, this path sorts a copy and today's statement sorts the vector the caller goes
// on with.
fn sort_indices_as(
    indices: &mut [usize],
    records: &[UngappedHit],
    compare: fn(&UngappedHit, &UngappedHit) -> Ordering,
    site: Site,
    shadow: bool,
    on_pool: bool,
    parallel_min: usize,
) {
    if shadow {
        let mut new = indices.to_vec();
        let parallel = sort_stable(&mut new, records, compare, on_pool, parallel_min);
        count(new.len(), parallel, false);
        sort_today(indices, records, compare);
        shadow_check(site, &new, indices);
    } else {
        let parallel = sort_stable(indices, records, compare, on_pool, parallel_min);
        count(indices.len(), parallel, false);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
// ```c
// if (kTranslatedQuery) {
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTransl);
//    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_FwdCompareHSPsTransl);
// ```
// The two final sorts of `x_hsp_inplace::replay_in_place_counted` with the switch or its shadow
// on: one stable sort with `s_FwdCompareHSPsTransl`, ties by `s_RevCompareHSPsTransl` (module doc,
// item 1).
/// Sorts `order` into the order of today's two final sorts over `records`.
pub(super) fn final_sort(order: &mut [usize], records: &[UngappedHit], site: Site) {
    final_sort_as(
        order,
        records,
        rev_compare_hsps_transl,
        fwd_compare_hsps_transl,
        site,
        mode() == 2,
        pool_allowed(site),
        PARALLEL_MIN,
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
// ```c
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTransl);
//    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_FwdCompareHSPsTransl);
// ```
// The body of `final_sort` for any two comparators: `first` and `second` are today's two sorts
// in the order they run, so `second` decides and `first` breaks its ties. With `shadow`, this
// path sorts a copy and today's two statements sort the vector the caller goes on with.
#[allow(clippy::too_many_arguments)]
fn final_sort_as<B, A>(
    order: &mut [usize],
    records: &[UngappedHit],
    first: B,
    second: A,
    site: Site,
    shadow: bool,
    on_pool: bool,
    parallel_min: usize,
) where
    B: Fn(&UngappedHit, &UngappedHit) -> Ordering + Sync,
    A: Fn(&UngappedHit, &UngappedHit) -> Ordering + Sync,
{
    let composite = |a: &UngappedHit, b: &UngappedHit| second(a, b).then_with(|| first(a, b));
    if shadow {
        let mut new = order.to_vec();
        let parallel = sort_stable(&mut new, records, composite, on_pool, parallel_min);
        count(new.len(), parallel, true);
        sort_today(order, records, &first);
        sort_today(order, records, &second);
        shadow_check(site, &new, order);
    } else {
        let parallel = sort_stable(order, records, composite, on_pool, parallel_min);
        count(order.len(), parallel, true);
    }
}

#[cfg(test)]
mod tests {
    use super::super::rev_compare_hsps_tbx;
    use super::super::tests::mock_hit;
    use super::*;
    use crate::algorithm::tblastx::blast_gapalign::score_compare_ungapped_hits_ncbi;

    type Compare = fn(&UngappedHit, &UngappedHit) -> Ordering;

    /// The four comparators the index sorts receive (module doc, item 3).
    const COMPARATORS: [(&str, Compare); 4] = [
        ("ScoreCompareHSPs", score_compare_ungapped_hits_ncbi),
        ("s_RevCompareHSPsTbx", rev_compare_hsps_tbx),
        ("s_RevCompareHSPsTransl", rev_compare_hsps_transl),
        ("s_FwdCompareHSPsTransl", fwd_compare_hsps_transl),
    ];

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

    /// `n` random HSPs of two queries. A small `span` gives few distinct offsets and scores, so
    /// many HSPs that a comparator calls equal; every record differs in `hsp_list_order`.
    fn random_hits(rng: &mut Rng, n: usize, span: u64) -> Vec<UngappedHit> {
        (0..n)
            .map(|index| {
                let q0 = rng.below(span) as usize;
                let s0 = rng.below(span) as usize;
                let q_len = 4 + rng.below(3) as usize;
                let s_len = 4 + rng.below(3) as usize;
                let score = 38 + rng.below(span.min(40)) as i32;
                let mut hit = mock_hit(q0, q0 + q_len, s0, s0 + s_len, score);
                hit.ctx_idx = rng.below(12) as usize;
                hit.q_idx = (hit.ctx_idx / 6) as u32;
                hit.q_frame = FRAMES[hit.ctx_idx % 6];
                hit.s_frame = FRAMES[rng.below(6) as usize];
                hit.hsp_list_order = index;
                hit
            })
            .collect()
    }

    /// The lists the sorts are tried on: random, already in the comparator's order, and in the
    /// reverse of it (whole runs that descend, with equal neighbours inside them).
    fn arrangements(rng: &mut Rng, n: usize, span: u64, compare: Compare) -> Vec<Vec<UngappedHit>> {
        let random = random_hits(rng, n, span);
        let mut ascending = random.clone();
        ascending.sort_by(compare);
        let mut descending = ascending.clone();
        descending.reverse();
        vec![random, ascending, descending]
    }

    // No NCBI counterpart: `body` once on the test thread and once inside pools of 1, 2 and 8
    // threads (a serial build has no pool and runs it on the test thread only). The flag tells
    // `body` whether it runs on a pool worker.
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

    /// Does this build have the parallel sort?
    const PARALLEL_BUILD: bool = cfg!(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ));

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:359-366
    // ```c
    // 	if (h1->query.offset < h2->query.offset)
    // 		return  1;
    // 	if (h1->query.offset > h2->query.offset)
    // 		return -1;
    // 	if (h1->query.end < h2->query.end)
    // 		return  1;
    // 	if (h1->query.end > h2->query.end)
    // 		return -1;
    // ```
    // Module doc, "Terms" and item 3: each comparator, and the composite of the final sort, is a
    // total preorder on random records (reversal, and transitivity of `Less` and of `Equal`).
    #[test]
    fn the_comparators_are_total_preorders() {
        let composite: Compare =
            |a, b| fwd_compare_hsps_transl(a, b).then_with(|| rev_compare_hsps_transl(a, b));
        let mut all = COMPARATORS.to_vec();
        all.push(("composite of the final sort", composite));
        let mut rng = Rng(0x243F_6A88_85A3_08D3);
        for (name, compare) in all {
            let hits = random_hits(&mut rng, 60, 2);
            let mut equal_pairs = 0usize;
            for a in &hits {
                assert_eq!(compare(a, a), Ordering::Equal, "{name}");
                for b in &hits {
                    let ab = compare(a, b);
                    assert_eq!(compare(b, a), ab.reverse(), "{name}");
                    equal_pairs += usize::from(ab == Ordering::Equal);
                    for c in &hits {
                        let bc = compare(b, c);
                        if ab != Ordering::Greater && bc != Ordering::Greater {
                            let expected = if ab == Ordering::Equal && bc == Ordering::Equal {
                                Ordering::Equal
                            } else {
                                Ordering::Less
                            };
                            assert_eq!(compare(a, c), expected, "{name}");
                        }
                    }
                }
            }
            // Different records that the comparator calls equal must be among them.
            assert!(equal_pairs > hits.len(), "{name}: {equal_pairs}");
        }
    }

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
    // Module doc, item 2: this path's index order equals today's for each comparator, on and off
    // pools of 1, 2 and 8 threads, for lengths around rayon's insertion-sort limit (20) and chunk
    // length (2,000) with every list sent to the parallel sort, and for lengths below, at and
    // above `PARALLEL_MIN` with the real threshold. Random lists, lists already in order, and
    // lists in the reverse order (rayon reverses a chunk that descends strictly).
    #[test]
    fn the_index_order_is_todays_at_any_number_of_threads() {
        let short: [usize; 11] = [0, 1, 2, 3, 20, 21, 1999, 2000, 2001, 3999, 4001];
        let mut rng = Rng(0x0123_4567_89AB_CDEF);
        for (name, compare) in COMPARATORS {
            // (threshold, list, today's index order)
            let mut cases: Vec<(usize, Vec<UngappedHit>, Vec<usize>)> = Vec::new();
            let mut add = |parallel_min: usize, hits: Vec<UngappedHit>| {
                let mut today: Vec<usize> = (0..hits.len()).collect();
                sort_today(&mut today, &hits, compare);
                cases.push((parallel_min, hits, today));
            };
            for n in short {
                for span in [3u64, 100_000] {
                    for hits in arrangements(&mut rng, n, span, compare) {
                        add(0, hits);
                    }
                }
            }
            add(PARALLEL_MIN, random_hits(&mut rng, PARALLEL_MIN - 1, 60));
            for hits in arrangements(&mut rng, PARALLEL_MIN, 60, compare) {
                add(PARALLEL_MIN, hits);
            }
            add(PARALLEL_MIN, random_hits(&mut rng, 70_001, 60));
            let equal_neighbours: usize = cases
                .iter()
                .map(|(_, hits, today)| {
                    today
                        .windows(2)
                        .filter(|pair| compare(&hits[pair[0]], &hits[pair[1]]) == Ordering::Equal)
                        .count()
                })
                .sum();
            // The lists must hold comparator-equal HSPs, whose order only stability decides.
            assert!(equal_neighbours > 10_000, "{name}: {equal_neighbours}");
            on_and_off_pools(|where_, in_pool| {
                for (case, (parallel_min, hits, today)) in cases.iter().enumerate() {
                    let n = hits.len();
                    // Not allowed to use the pool: today's statement, once is enough.
                    for on_pool in [false, true] {
                        if !on_pool && in_pool {
                            continue;
                        }
                        let label = format!("{name} case {case} n {n} {where_} on_pool {on_pool}");
                        let mut new: Vec<usize> = (0..n).collect();
                        let parallel = sort_stable(&mut new, hits, compare, on_pool, *parallel_min);
                        assert_eq!(&new, today, "{label}");
                        assert_eq!(
                            parallel,
                            PARALLEL_BUILD && on_pool && n >= *parallel_min,
                            "{label}"
                        );
                    }
                }
            });
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
    // ```c
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTbx);
    // ```
    // (S) of the module doc does not need the indices to start as `0..n`: the same equality from
    // a shuffled index vector, where the stable order is "comparator, then position in the
    // vector".
    #[test]
    fn the_index_order_is_todays_from_any_starting_order() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        for (name, compare) in COMPARATORS {
            for n in [50usize, 2500, 12_000] {
                let hits = random_hits(&mut rng, n, 6);
                let mut start: Vec<usize> = (0..n).collect();
                for index in (1..n).rev() {
                    start.swap(index, rng.below(index as u64 + 1) as usize);
                }
                let mut today = start.clone();
                sort_today(&mut today, &hits, compare);
                on_and_off_pools(|where_, _| {
                    let mut new = start.clone();
                    sort_stable(&mut new, &hits, compare, true, 0);
                    assert_eq!(new, today, "{name} n {n} {where_}");
                });
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
    // ```c
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTransl);
    //    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_FwdCompareHSPsTransl);
    // ```
    // Module doc, item 1: one stable sort with "second, ties by first" equals the two stable
    // sorts in sequence. The pair of the linking call, and pairs of comparators whose ties do
    // not coincide; the lists hold pairs of HSPs on which the two comparators disagree, on which
    // only the second ties, and on which both tie.
    #[test]
    fn one_composite_sort_equals_the_two_sorts_in_sequence() {
        let pairs: [(&str, Compare, Compare); 5] = [
            (
                "the linking call: s_RevCompareHSPsTransl, then s_FwdCompareHSPsTransl",
                rev_compare_hsps_transl,
                fwd_compare_hsps_transl,
            ),
            (
                "ScoreCompareHSPs, then s_FwdCompareHSPsTransl",
                score_compare_ungapped_hits_ncbi,
                fwd_compare_hsps_transl,
            ),
            (
                "s_RevCompareHSPsTbx, then ScoreCompareHSPs",
                rev_compare_hsps_tbx,
                score_compare_ungapped_hits_ncbi,
            ),
            (
                "s_RevCompareHSPsTbx, then s_FwdCompareHSPsTransl",
                rev_compare_hsps_tbx,
                fwd_compare_hsps_transl,
            ),
            (
                "ScoreCompareHSPs, then s_RevCompareHSPsTbx",
                score_compare_ungapped_hits_ncbi,
                rev_compare_hsps_tbx,
            ),
        ];
        let mut rng = Rng(0xD1B5_4A32_D192_ED03);
        for (pair_index, (name, first, second)) in pairs.into_iter().enumerate() {
            let mut disagree = 0usize;
            let mut second_ties_only = 0usize;
            let mut both_tie = 0usize;
            let mut cases: Vec<(usize, Vec<UngappedHit>)> = Vec::new();
            for (n, span, parallel_min) in [
                (0usize, 2u64, 0usize),
                (1, 2, 0),
                (2, 2, 0),
                (300, 2, 0),
                (300, 40, 0),
                (5000, 2, 0),
                (5000, 100_000, 0),
                (20_000, 1, 0),
                (PARALLEL_MIN - 1, 3, PARALLEL_MIN),
                (PARALLEL_MIN + 3, 3, PARALLEL_MIN),
            ] {
                let hits = random_hits(&mut rng, n, span);
                for pair in hits.windows(2) {
                    let (f, s) = (first(&pair[0], &pair[1]), second(&pair[0], &pair[1]));
                    disagree += usize::from(f != Ordering::Equal && f == s.reverse());
                    second_ties_only += usize::from(s == Ordering::Equal && f != Ordering::Equal);
                    both_tie += usize::from(s == Ordering::Equal && f == Ordering::Equal);
                }
                cases.push((parallel_min, hits));
            }
            assert!(disagree > 100, "{name}: {disagree}");
            assert!(both_tie > 10, "{name}: {both_tie}");
            if pair_index == 0 {
                // The pair of the linking call reads the same keys: a tie of one is a tie of both.
                assert_eq!(second_ties_only, 0, "{name}");
            } else {
                assert!(second_ties_only > 10, "{name}: {second_ties_only}");
            }
            on_and_off_pools(|where_, _| {
                for (parallel_min, hits) in &cases {
                    let n = hits.len();
                    let mut today: Vec<usize> = (0..n).collect();
                    sort_today(&mut today, hits, first);
                    sort_today(&mut today, hits, second);
                    for (shadow, on_pool) in [(false, false), (false, true), (true, true)] {
                        let mut order: Vec<usize> = (0..n).collect();
                        final_sort_as(
                            &mut order,
                            hits,
                            first,
                            second,
                            Site::Final(on_pool),
                            shadow,
                            on_pool,
                            *parallel_min,
                        );
                        assert_eq!(
                            order, today,
                            "{name} n {n} {where_} shadow {shadow} on_pool {on_pool}"
                        );
                    }
                }
            });
        }
    }

    /// The text of a panic raised by `body`, if it panics.
    fn panic_text(body: impl FnOnce() + std::panic::UnwindSafe) -> Option<String> {
        std::panic::catch_unwind(body).err().map(|payload| {
            payload
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| payload.downcast_ref::<&str>().map(|text| text.to_string()))
                .unwrap_or_default()
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
    // ```c
    //    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //          s_RevCompareHSPsTbx);
    // ```
    // The shadow's comparison passes equal index vectors and stops on a planted difference with
    // the switch, the sort, the length and the first differing position in its message.
    #[test]
    fn the_shadow_check_reports_a_planted_difference() {
        let today: Vec<usize> = vec![4, 0, 3, 1, 2, 5];
        for site in [Site::First(true), Site::Final(false), Site::Dispatch] {
            assert_eq!(panic_text(|| shadow_check(site, &today, &today)), None);
        }
        assert_eq!(panic_text(|| shadow_check(Site::Dispatch, &[], &[])), None);

        // Two neighbours exchanged: what a sort that is not stable does to equal HSPs.
        let mut planted = today.clone();
        planted.swap(2, 3);
        for (site, name) in [
            (Site::First(false), "first"),
            (Site::Final(true), "final"),
            (Site::Dispatch, "dispatch"),
        ] {
            let text = panic_text(|| shadow_check(site, &planted, &today)).expect("a difference");
            assert!(text.starts_with("LOSAT_X_LINKSORTSHADOW: "), "{text}");
            assert!(
                text.contains(&format!("the {name} index sort of 6 indices")),
                "{text}"
            );
            assert!(text.contains("differs at position 2:"), "{text}");
            assert!(text.contains("this path Some(1), today Some(3)"), "{text}");
        }

        // A difference in the last element, and one of length.
        let mut last = today.clone();
        last[5] = 6;
        let text = panic_text(|| shadow_check(Site::Dispatch, &last, &today)).expect("last");
        assert!(text.contains("differs at position 5:"), "{text}");
        let text =
            panic_text(|| shadow_check(Site::Dispatch, &today[..4], &today)).expect("length");
        assert!(text.contains("differs at position 4:"), "{text}");
        assert!(text.contains("this path None, today Some(2)"), "{text}");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
    // ```c
    //         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //               ScoreCompareHSPs);
    // ```
    // With the shadow on, the caller's vector holds today's order (today's statement sorted it),
    // on and off the pool, and the sort is counted as compared. The counters are shared by the
    // tests of this process, so only their growth is checked.
    #[test]
    fn a_shadowed_sort_leaves_todays_order_and_counts() {
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        for (name, compare) in COMPARATORS {
            let hits = random_hits(&mut rng, 6000, 5);
            let mut today: Vec<usize> = (0..hits.len()).collect();
            sort_today(&mut today, &hits, compare);
            on_and_off_pools(|where_, _| {
                for on_pool in [false, true] {
                    let before = (
                        SORTS.load(AtomicOrdering::Relaxed),
                        SHADOW_SORTS.load(AtomicOrdering::Relaxed),
                        SHADOW_INDICES.load(AtomicOrdering::Relaxed),
                        PARALLEL_SORTS.load(AtomicOrdering::Relaxed),
                    );
                    let mut indices: Vec<usize> = (0..hits.len()).collect();
                    sort_indices_as(
                        &mut indices,
                        &hits,
                        compare,
                        Site::Dispatch,
                        true,
                        on_pool,
                        0,
                    );
                    assert_eq!(indices, today, "{name} {where_} on_pool {on_pool}");
                    assert!(SORTS.load(AtomicOrdering::Relaxed) > before.0);
                    assert!(SHADOW_SORTS.load(AtomicOrdering::Relaxed) > before.1);
                    assert!(SHADOW_INDICES.load(AtomicOrdering::Relaxed) >= before.2 + 6000);
                    if PARALLEL_BUILD && on_pool {
                        assert!(PARALLEL_SORTS.load(AtomicOrdering::Relaxed) > before.3);
                    }
                }
            });
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149
    // ```c
    // TBlastThreads the_threads(GetNumberOfThreads());
    // ```
    // Module doc, "The pool and `-num_threads`": a linking call's sorts follow its flag; a
    // dispatch sort follows the calling thread, which is not a rayon worker off a pool (the
    // thread of a `-num_threads 1` search) and is one inside `ThreadPool::install`, where the
    // current pool has the threads it was built with.
    #[test]
    fn the_pool_decision_follows_the_flag_or_the_calling_thread() {
        on_and_off_pools(|where_, in_pool| {
            for site in [Site::First(false), Site::Final(false)] {
                assert!(!pool_allowed(site), "{site:?} {where_}");
            }
            for site in [Site::First(true), Site::Final(true)] {
                assert_eq!(pool_allowed(site), PARALLEL_BUILD, "{site:?} {where_}");
            }
            assert_eq!(
                pool_allowed(Site::Dispatch),
                PARALLEL_BUILD && in_pool,
                "{where_}"
            );
        });
        #[cfg(feature = "parallel")]
        {
            assert_eq!(rayon::current_thread_index(), None);
            for threads in [1usize, 2, 8] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("test pool");
                let (index, size) =
                    pool.install(|| (rayon::current_thread_index(), rayon::current_num_threads()));
                assert!(index.is_some_and(|index| index < threads));
                assert_eq!(size, threads);
            }
        }
    }
}
