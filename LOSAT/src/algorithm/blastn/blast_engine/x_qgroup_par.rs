//! EXPERIMENT (LOSAT_X_QGROUPPAR, shadow LOSAT_X_QGROUPPARSHADOW): the gapped stages of a
//! blastn query batch against one subject run per query on the search pool: the preliminary
//! gapped loop with its purge, and the traceback loop with its two purges and the
//! re-evaluation between them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3917-3922
//! ```c
//!    for (index=0; index<init_hitlist->total; index++)
//! ...
//!       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
//!                                         hit_options->min_diag_separation))
//!       {
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,404-405
//! ```c
//!    for (index=0; index < num_initial_hsps; index++) {
//! ...
//!           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
//!                                hit_options->min_diag_separation)) {
//! ```
//! NCBI runs both loops for one subject on one thread, over the HSPs of all queries of the batch
//! in one order (score order for the preliminary loop, the hit lists' order for the traceback),
//! with one interval tree. `search_query_batch` does the same on the thread that owns the subject;
//! with one subject (`-subject` with one record) the other threads of the pool only compute
//! alignments ahead of that loop (LOSAT_X_AHEAD), which leaves most of them waiting: at 16 threads
//! the two loops, the purges and the sorts around them were 57 % of a blastn batch's wall
//! (31,168 cDNAs against 3 Mb; v4-D session 8 probe).
//!
//! Why a query's results do not depend on the other queries: the tree's answers are per context.
//! `containing_hsp` returns false across contexts, and `add_hsp`'s walk over the HSPs that share
//! an endpoint with the new one (which rejects the new HSP or unlinks an old one) also compares
//! only HSPs of the same context; both walks visit every node that can hold such an HSP, so the
//! tree's shape (which depends on what else it holds) cannot change an answer, and the
//! LOSAT_X_ITREEFAST grid gives the tree's yes / no (`interval_tree.rs`). Only the trace reads
//! which HSP contains another, and the trace keeps the owner loop. The purges sort by context
//! first and compare only within a context (`purge_prelim_hits_with_common_endpoints`,
//! `purge_hsps_with_common_endpoints_ex`), the re-evaluation and the identity test read one HSP
//! and its context's cutoff and statistics, the alignment scratch carries nothing from one call
//! to the next that changes a result (LOSAT_X_AHEAD's helpers already compute with scratch of
//! their own), and the sorts are stable. So running the loop of one query over that query's HSPs,
//! in the relative order they have in the batch's list, with a tree of its own, gives the HSPs the
//! batch's loop gives for that query, in the same relative order.
//!
//! With the switch on, a search pool of two or more threads and one subject on the owner thread
//! (`!use_parallel`; with several subjects the pool already searches subjects in parallel):
//! - Preliminary stage (`collect_prelim_hits_for_chunk`, one chunk at a time on the owner): the
//!   chunk's ungapped hits are split by query in the order the seed stage saved them, each query's
//!   hits are sorted with the batch's comparator (stable, so in the batch's relative order), looped
//!   and purged; the purged lists are concatenated in query order. The serial purge's output is
//!   ordered by context first, so this concatenation is the serial purge's output itself, and the
//!   stable score sort after it (unchanged) gives the serial `chunk_prelim_hits` element for
//!   element: the merge of a subject's chunks and `-subject_besthit` see the same list.
//! - Traceback (`process_subject` with the kept HSPs): the kept HSPs are the hit lists of the
//!   queries one after another (`prelim_hits_by_subject`), so each query's HSPs are one slice. Each
//!   slice runs the traceback loop, the first purge (trimmed HSPs at its end), the re-evaluation of
//!   the trimmed HSPs, the second purge, the score sort, the final tree and the e-values; the final
//!   lists are concatenated in query order and sorted by score with the same stable sort. The
//!   first purge leaves the HSPs it keeps in context order, then the HSPs its second pass trimmed,
//!   then those its first pass trimmed, each group in the reverse order of trimming; one query's
//!   trims happen in the same relative order in both forms, so its HSPs come in the same relative
//!   order, and the re-evaluation takes the same HSPs (the trimmed ones). The second purge orders
//!   by context first and keeps, within a context, the input's relative order for its ties, so a
//!   query's list after it is the serial list's part for that query; the score sort, the tree and
//!   the e-values keep that per query; and the serial score sort leaves ties between queries in
//!   the order of its input, which is query order, as the stable sort of the concatenation does.
//!   So the subject's final list is the serial one element for element.
//!
//! The queries are handed out to the threads from an atomic counter, the query with the most HSPs
//! first (the counts are known before each stage), each thread with its own tree and alignment
//! scratch; measured on 5,000 cDNAs, that order keeps 16 threads 0.99 busy where contiguous groups
//! of queries kept them 0.66-0.79 busy. The owner loop with LOSAT_X_AHEAD stays where a batch has
//! fewer than two queries with HSPs or one query holds more than a quarter of the stage's HSPs (a
//! long genome query: the loop's helpers then do better than the queries), for a subject searched
//! in chunks in parallel, and on the diagnostic paths (`-verbose`, LOSAT_DEBUG_BLASTN,
//! BLEMIR_DEBUG, LOSAT_DEBUG_COORDS, the trace of `utils/blastn_trace`), which print in the loop's
//! order.
//!
//! The shadow runs both forms of each stage and compares the lists element by element (every
//! field, through `Debug`, which prints floats exactly); on a difference it prints
//! `[X_SHADOW] QGROUPPAR ... DIFFERENT` and panics, otherwise `main` prints
//! `[X_SHADOW] QGROUPPAR prelim= traceback= hsps= diffs=0 identical` at exit (`prelim=` and
//! `traceback=` count the stages compared, `hsps=` the HSPs compared). The search keeps the owner
//! loop's lists. With `LOSAT_X_STATS`, `main` prints
//! `[X_STATS] QGROUPPAR prelim= prelim_queries= traceback= traceback_queries= owner=` when a stage
//! was planned: the stages run per query, the query tasks they ran and the stages with HSPs that
//! kept the owner loop by the rule. With `LOSAT_WASI_THREADS_DEBUG=1` each stage run per query
//! prints `[losat-thread-stage] program=blastn stage=prelim_queries|traceback_queries`. In the
//! per-query form the `[TIMING]` steps inside the stages (alignment, purge passes, re-evaluation,
//! final tree) are sums over the threads, while `gapped_extend` and `traceback_prune` stay the
//! stages' wall times.
//!
//! The chunk's ungapped hits are grouped by query in place (a counting sort's destinations, 4
//! bytes per hit, applied by swaps), so the stage holds no second copy of them; the shadow keeps
//! one for the owner loop.
//!
//! Programs that reach it: BLASTN only (megablast, dc-megablast, blastn, blastn-short), from
//! `search_query_batch` (the loop, a LOSAT_X_QBATCHPAR slot of two or more threads, and a `-db`
//! search of a database of one sequence).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use crate::utils::threading::SearchPool;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

// No NCBI counterpart: reads the switches once (LOSAT_X_QGROUPPAR is a strict-set switch: on
// unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set; the shadow, off unless set to any
// value, implies the switch). They choose between two LOSAT code paths with the same result and
// change no value NCBI computes.
/// 0 = off, 1 = LOSAT_X_QGROUPPAR, 2 = LOSAT_X_QGROUPPARSHADOW.
fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_QGROUPPARSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_QGROUPPAR") {
            1
        } else {
            0
        }
    })
}

/// Whether the stages may run per query (the switch or its shadow).
pub(super) fn on() -> bool {
    mode() != 0
}

/// Whether both forms run and are compared.
pub(super) fn shadow() -> bool {
    mode() == 2
}

// No NCBI counterpart: counters of the summary lines; they change no value NCBI computes.
static PRELIM: AtomicU64 = AtomicU64::new(0);
static PRELIM_QUERIES: AtomicU64 = AtomicU64::new(0);
static TRACEBACK: AtomicU64 = AtomicU64::new(0);
static TRACEBACK_QUERIES: AtomicU64 = AtomicU64::new(0);
static OWNER: AtomicU64 = AtomicU64::new(0);
static COMPARED_PRELIM: AtomicU64 = AtomicU64::new(0);
static COMPARED_TRACEBACK: AtomicU64 = AtomicU64::new(0);
static COMPARED_HSPS: AtomicU64 = AtomicU64::new(0);
static DIFFS: AtomicU64 = AtomicU64::new(0);

/// Which stage a plan or a comparison is for.
#[derive(Clone, Copy)]
pub(super) enum Stage {
    Prelim,
    Traceback,
}

/// The number of HSPs of each of `num_queries` queries in `query_of` (query indices).
pub(super) fn counts(query_of: impl Iterator<Item = u32>, num_queries: usize) -> Vec<u32> {
    let mut counts = vec![0u32; num_queries];
    for query in query_of {
        counts[query as usize] += 1;
    }
    counts
}

/// What a stage does with its HSPs per query.
#[derive(Debug, PartialEq, Eq)]
enum Choice {
    /// No HSPs: nothing to run.
    Empty,
    /// The owner loop: fewer than two queries with HSPs, or one query with more than a quarter.
    Owner,
    /// The queries as tasks on the pool.
    Queries,
}

fn choose(counts: &[u32]) -> Choice {
    let total: u64 = counts.iter().map(|&count| u64::from(count)).sum();
    let largest = counts.iter().copied().max().unwrap_or(0);
    let with_hsps = counts.iter().filter(|&&count| count > 0).count();
    if total == 0 {
        Choice::Empty
    } else if with_hsps < 2 || u64::from(largest) * 4 > total {
        Choice::Owner
    } else {
        Choice::Queries
    }
}

/// The order in which the queries with HSPs are handed out (most HSPs first, then query order),
/// or `None` when the owner loop stays: no HSPs, fewer than two queries with HSPs, or one query
/// with more than a quarter of the HSPs (counted as `owner=` unless there are no HSPs).
pub(super) fn plan(counts: &[u32], stage: Stage) -> Option<Vec<u32>> {
    match choose(counts) {
        Choice::Empty => return None,
        Choice::Owner => {
            OWNER.fetch_add(1, Ordering::Relaxed);
            return None;
        }
        Choice::Queries => {}
    }
    let mut order: Vec<u32> = (0..counts.len() as u32)
        .filter(|&query| counts[query as usize] > 0)
        .collect();
    order.sort_by(|&a, &b| counts[b as usize].cmp(&counts[a as usize]).then(a.cmp(&b)));
    let (stages, queries) = match stage {
        Stage::Prelim => (&PRELIM, &PRELIM_QUERIES),
        Stage::Traceback => (&TRACEBACK, &TRACEBACK_QUERIES),
    };
    stages.fetch_add(1, Ordering::Relaxed);
    queries.fetch_add(order.len() as u64, Ordering::Relaxed);
    Some(order)
}

/// Where each query's HSPs start in a list ordered by query whose per-query counts are `counts`
/// (one more entry than `counts`: the last is the list's length).
pub(super) fn starts(counts: &[u32]) -> Vec<usize> {
    let mut starts = Vec::with_capacity(counts.len() + 1);
    let mut at = 0usize;
    starts.push(at);
    for &count in counts {
        at += count as usize;
        starts.push(at);
    }
    starts
}

/// Groups `items` by query (`query_of`) in place, each query's items in their order in `items`:
/// the destinations of a counting sort (4 bytes per item) applied by swaps, so no second copy of
/// the items is made. With `starts(counts)` the slice of query `q` is then
/// `starts[q]..starts[q + 1]`.
pub(super) fn group_in_place<T>(items: &mut [T], query_of: impl Fn(&T) -> u32, counts: &[u32]) {
    let mut next = starts(counts);
    next.pop();
    let mut destination: Vec<u32> = items
        .iter()
        .map(|item| {
            let query = query_of(item) as usize;
            let at = next[query];
            next[query] += 1;
            u32::try_from(at).expect("fewer than 2^32 hits in a chunk")
        })
        .collect();
    // Each swap puts the item at `at` where it belongs, so at most one swap per item.
    for at in 0..items.len() {
        while destination[at] as usize != at {
            let to = destination[at] as usize;
            items.swap(at, to);
            destination.swap(at, to);
        }
    }
}

/// Whether the query indices of `query_of` never decrease (each query's items are one slice).
pub(super) fn in_query_order(mut query_of: impl Iterator<Item = u32>) -> bool {
    let Some(mut last) = query_of.next() else {
        return true;
    };
    for query in query_of {
        if query < last {
            return false;
        }
        last = query;
    }
    true
}

/// Runs `task` once for every query of `order` on `pool`: the threads take the queries from an
/// atomic counter in the order given, each with a state of its own made by `new_state`. Returns
/// each query's result at its query index (`None` for the queries not in `order`).
pub(super) fn run<S, R>(
    pool: &SearchPool<'_>,
    stage: Stage,
    order: &[u32],
    num_queries: usize,
    new_state: impl Fn() -> S + Sync,
    task: impl Fn(&mut S, u32) -> R + Sync,
) -> Vec<Option<R>>
where
    R: Send,
{
    let mut results: Vec<Option<R>> = Vec::with_capacity(num_queries);
    results.resize_with(num_queries, || None);
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    {
        let workers = pool.threads().min(order.len());
        crate::utils::threading::report_stage(
            "blastn",
            match stage {
                Stage::Prelim => "prelim_queries",
                Stage::Traceback => "traceback_queries",
            },
            order.len(),
            pool.enabled() && workers > 1,
        );
        if pool.enabled() && workers > 1 {
            let next = std::sync::atomic::AtomicUsize::new(0);
            let parts: Vec<Vec<(u32, R)>> = pool.install(|| {
                (0..workers)
                    .into_par_iter()
                    .map(|_| {
                        let mut state = new_state();
                        let mut done = Vec::new();
                        loop {
                            let at = next.fetch_add(1, Ordering::Relaxed);
                            let Some(&query) = order.get(at) else {
                                break;
                            };
                            done.push((query, task(&mut state, query)));
                        }
                        done
                    })
                    .collect()
            });
            for (query, result) in parts.into_iter().flatten() {
                results[query as usize] = Some(result);
            }
            return results;
        }
    }
    let _ = stage;
    let mut state = new_state();
    for &query in order {
        results[query as usize] = Some(task(&mut state, query));
    }
    results
}

/// Sorts `items` stably by `compare`: on `pool` with rayon's parallel merge sort, which is stable
/// as `sort_by` is, so both give the same order; serially without a pool of two or more threads.
pub(super) fn sort_stable<T: Send>(
    pool: &SearchPool<'_>,
    items: &mut [T],
    compare: impl Fn(&T, &T) -> std::cmp::Ordering + Send + Sync,
) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool.enabled() && pool.threads() > 1 {
        pool.install(|| items.par_sort_by(compare));
        return;
    }
    let _ = pool;
    items.sort_by(compare);
}

/// The shadow's comparison of the owner loop's list (`slow`) with the per-query form's (`fast`):
/// counts the comparison, and on a difference prints the `DIFFERENT` line and panics.
pub(super) fn compare<T: std::fmt::Debug>(stage: Stage, slow: &[T], fast: &[T]) {
    let (compared, name) = match stage {
        Stage::Prelim => (&COMPARED_PRELIM, "preliminary"),
        Stage::Traceback => (&COMPARED_TRACEBACK, "traceback"),
    };
    compared.fetch_add(1, Ordering::Relaxed);
    COMPARED_HSPS.fetch_add(slow.len() as u64, Ordering::Relaxed);
    let first = (0..slow.len().max(fast.len()))
        .find(|&at| format!("{:?}", slow.get(at)) != format!("{:?}", fast.get(at)));
    if let Some(at) = first {
        DIFFS.fetch_add(1, Ordering::Relaxed);
        eprintln!(
            "[X_SHADOW] QGROUPPAR {} HSPs {} / {} DIFFERENT (first at {}: {:?} / {:?})",
            name,
            slow.len(),
            fast.len(),
            at,
            slow.get(at),
            fast.get(at)
        );
        panic!("LOSAT_X_QGROUPPARSHADOW: the {name} HSPs of the query tasks differ from the loop's at {at}");
    }
}

/// The lines `main` prints at exit: the shadow's summary (with the shadow) and the counters (with
/// `LOSAT_X_STATS`, when a stage was planned). No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub(super) fn print_summary() {
    if shadow() {
        let diffs = DIFFS.load(Ordering::Relaxed);
        eprintln!(
            "[X_SHADOW] QGROUPPAR prelim={} traceback={} hsps={} diffs={diffs} {}",
            COMPARED_PRELIM.load(Ordering::Relaxed),
            COMPARED_TRACEBACK.load(Ordering::Relaxed),
            COMPARED_HSPS.load(Ordering::Relaxed),
            if diffs == 0 { "identical" } else { "DIFFERENT" }
        );
    }
    let planned = PRELIM.load(Ordering::Relaxed)
        + TRACEBACK.load(Ordering::Relaxed)
        + OWNER.load(Ordering::Relaxed);
    if std::env::var_os("LOSAT_X_STATS").is_some() && planned > 0 {
        eprintln!(
            "[X_STATS] QGROUPPAR prelim={} prelim_queries={} traceback={} traceback_queries={} owner={}",
            PRELIM.load(Ordering::Relaxed),
            PRELIM_QUERIES.load(Ordering::Relaxed),
            TRACEBACK.load(Ordering::Relaxed),
            TRACEBACK_QUERIES.load(Ordering::Relaxed),
            OWNER.load(Ordering::Relaxed)
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_plan_hands_out_the_largest_queries_first() {
        let order = plan(&[3, 0, 3, 3, 1, 2], Stage::Prelim).expect("no query above a quarter");
        assert_eq!(order, vec![0, 2, 3, 5, 4]);
    }

    #[test]
    fn the_owner_loop_stays_for_one_query_or_a_dominant_one() {
        assert!(plan(&[0, 7, 0], Stage::Traceback).is_none());
        assert!(plan(&[], Stage::Traceback).is_none());
        // 5 of 16 is above a quarter, 4 of 16 is not.
        assert!(plan(&[5, 4, 4, 3], Stage::Traceback).is_none());
        assert!(plan(&[4, 4, 4, 4], Stage::Traceback).is_some());
    }

    #[test]
    fn group_in_place_keeps_each_querys_order() {
        let mut items = [(2u32, 'a'), (0, 'b'), (2, 'c'), (1, 'd'), (0, 'e')];
        let counts = counts(items.iter().map(|item| item.0), 4);
        assert_eq!(counts, vec![2, 1, 2, 0]);
        group_in_place(&mut items, |item| item.0, &counts);
        assert_eq!(items, [(0, 'b'), (0, 'e'), (1, 'd'), (2, 'a'), (2, 'c')]);
        assert_eq!(starts(&counts), vec![0, 2, 3, 5, 5]);
        // Longer cycles: a reversed list of 3 queries x 50 items, each query's items numbered.
        let mut items: Vec<(u32, u32)> = (0..150u32).rev().map(|i| (i % 3, i)).collect();
        let three = super::counts(items.iter().map(|item| item.0), 3);
        group_in_place(&mut items, |item| item.0, &three);
        let mut want: Vec<(u32, u32)> = (0..150u32).rev().map(|i| (i % 3, i)).collect();
        want.sort_by_key(|item| item.0);
        assert_eq!(items, want);
    }

    #[test]
    fn in_query_order_reads_slices() {
        assert!(in_query_order([0u32, 0, 1, 3, 3].into_iter()));
        assert!(in_query_order(std::iter::empty()));
        assert!(!in_query_order([0u32, 2, 1].into_iter()));
    }

    #[test]
    fn the_rule_names_empty_stages_and_dominant_queries() {
        assert_eq!(choose(&[0, 0, 0]), Choice::Empty);
        assert_eq!(choose(&[]), Choice::Empty);
        assert_eq!(choose(&[0, 7, 0]), Choice::Owner);
        assert_eq!(choose(&[5, 4, 4, 3]), Choice::Owner);
        assert_eq!(choose(&[4, 4, 4, 4]), Choice::Queries);
    }

    #[test]
    fn sort_stable_keeps_ties_in_order() {
        let items: Vec<(u32, u32)> = (0..2_000u32).map(|i| ((i * 7) % 13, i)).collect();
        let mut want = items.clone();
        want.sort_by_key(|item| item.0);
        let mut serial = items.clone();
        sort_stable(&SearchPool::x_serial(), &mut serial, |a, b| a.0.cmp(&b.0));
        assert_eq!(serial, want);
        #[cfg(feature = "parallel")]
        crate::utils::threading::with_search_pool(4, "blastn", |pool| {
            let mut pooled = items.clone();
            sort_stable(pool, &mut pooled, |a, b| a.0.cmp(&b.0));
            assert_eq!(pooled, want);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    #[should_panic(expected = "differ from the loop's at 1")]
    fn compare_panics_on_a_difference() {
        compare(Stage::Prelim, &[1, 2, 3], &[1, 5, 3]);
    }

    #[test]
    #[should_panic(expected = "differ from the loop's at 2")]
    fn compare_panics_on_a_length_mismatch() {
        compare(Stage::Traceback, &[1, 2], &[1, 2, 3]);
    }

    /// On a pool of four threads every query of the order runs exactly once and its result lands
    /// at its index, also with more threads than queries.
    #[cfg(feature = "parallel")]
    #[test]
    fn run_on_a_pool_runs_each_query_once() {
        use std::sync::atomic::AtomicUsize;
        crate::utils::threading::with_search_pool(4, "blastn", |pool| {
            for order in [(0..40u32).rev().collect::<Vec<_>>(), vec![3u32, 1]] {
                let calls: Vec<AtomicUsize> = (0..41).map(|_| AtomicUsize::new(0)).collect();
                let results = run(
                    pool,
                    Stage::Prelim,
                    &order,
                    41,
                    || (),
                    |_, query| {
                        calls[query as usize].fetch_add(1, Ordering::SeqCst);
                        query * 3
                    },
                );
                for (query, result) in results.iter().enumerate() {
                    let ran = order.contains(&(query as u32));
                    assert_eq!(*result, ran.then_some(query as u32 * 3), "query {query}");
                    assert_eq!(calls[query].load(Ordering::SeqCst), usize::from(ran));
                }
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn run_gives_every_query_its_own_result() {
        let order = [4u32, 1, 2];
        let serial = run(
            &SearchPool::x_serial(),
            Stage::Traceback,
            &order,
            6,
            || 0usize,
            |calls, query| {
                *calls += 1;
                query * 10
            },
        );
        assert_eq!(serial, vec![None, Some(10), Some(20), None, Some(40), None]);
    }
}
