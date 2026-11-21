//! EXPERIMENT (LOSAT_X_QGROUPPAR): the gapped stages of a
//! blastn query batch against one subject run per query on the search pool: the preliminary
//! gapped loop with its purge, and the traceback loop with its two purges and the
//! re-evaluation between them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3917-3922
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,404-405
//! NCBI runs both loops for one subject on one thread, over the HSPs of all queries of the batch
//! in one order (score order for the preliminary loop, the hit lists' order for the traceback),
//! with one interval tree. `search_query_batch` does the same on the thread that owns the subject;
//! with one subject (`-subject` with one record) the other threads of the pool only compute
//! alignments ahead of that loop (the look-ahead of `utils/xahead.rs`), which leaves most of them waiting: at 16 threads
//! the two loops, the purges and the sorts around them were 57 % of a blastn batch's wall
//! (31,168 cDNAs against 3 Mb; v4-D session 8 probe).
//!
//! Why a query's results do not depend on the other queries: the tree's answers are per context.
//! `containing_hsp` returns false across contexts, and `add_hsp`'s walk over the HSPs that share
//! an endpoint with the new one (which rejects the new HSP or unlinks an old one) also compares
//! only HSPs of the same context; both walks visit every node that can hold such an HSP, so the
//! tree's shape (which depends on what else it holds) cannot change an answer, and the
//! interval-tree grid gives the tree's yes / no (`interval_tree.rs`). Only the trace reads
//! which HSP contains another, and the trace keeps the owner loop. The purges sort by context
//! first and compare only within a context (`purge_prelim_hits_with_common_endpoints`,
//! `purge_hsps_with_common_endpoints_ex`), the re-evaluation and the identity test read one HSP
//! and its context's cutoff and statistics, the alignment scratch carries nothing from one call
//! to the next that changes a result (the look-ahead helpers already compute with scratch of
//! their own), and the sorts are stable. So running the loop of one query over that query's HSPs,
//! in the relative order they have in the batch's list, with a tree of its own, gives the HSPs the
//! batch's loop gives for that query, in the same relative order.
//!
//! With the switch on, a search pool of two or more threads and one subject at a time on the owner
//! thread (`!use_parallel`: one subject, or several one at a time under LOSAT_X_BNFEWSUBJ; when the
//! pool searches subjects in parallel it is already busy):
//! - Preliminary stage (`prelim_hsp_list::gapped_hsp_list`, which
//!   `prelim_subject::prelim_chunk_hits` calls for one chunk at a time on the owner): the chunk's
//!   ungapped hits are split by query in the order the seed stage saved them, each query's hits are
//!   sorted with the batch's comparator (stable, so in the batch's relative order), looped and
//!   purged; the purged lists are concatenated in query order. The serial purge's output is ordered
//!   by context first, so this concatenation is the serial purge's output itself, and the stable
//!   score sort after it (unchanged) gives the serial `chunk_prelim_hits` element for element: the
//!   merge of a subject's chunks and `-subject_besthit` see the same list.
//! - Traceback (`traceback_list::traceback_subject_hits`, which
//!   `traceback_subject::traceback_subject` calls with the kept HSPs): the kept HSPs are the hit
//!   lists of the queries one after another (`prelim_hits_by_subject`), so each query's HSPs are
//!   one slice. Each slice runs the traceback loop, the first purge (trimmed HSPs at its end), the
//!   re-evaluation of the trimmed HSPs, the second purge, the score sort, the final tree and the
//!   e-values; the final lists are concatenated in query order and sorted by score with the same
//!   stable sort. The first purge leaves the HSPs it keeps in context order, then the HSPs its
//!   second pass trimmed, then those its first pass trimmed, each group in the reverse order of
//!   trimming; one query's trims happen in the same relative order in both forms, so its HSPs come
//!   in the same relative order, and the re-evaluation takes the same HSPs (the trimmed ones). The
//!   second purge orders by context first and keeps, within a context, the input's relative order
//!   for its ties, so a query's list after it is the serial list's part for that query; the score
//!   sort, the tree and the e-values keep that per query; and the serial score sort leaves ties
//!   between queries in the order of its input, which is query order, as the stable sort of the
//!   concatenation does. So the subject's final list is the serial one element for element.
//!
//! The queries are handed out to the threads from an atomic counter, the query with the most HSPs
//! first (the counts are known before each stage), each thread with its own tree and alignment
//! scratch; measured on 5,000 cDNAs, that order keeps 16 threads 0.99 busy where contiguous groups
//! of queries kept them 0.66-0.79 busy. The owner loop with the look-ahead stays where a batch has
//! fewer than two queries with HSPs or one query holds more than a quarter of the stage's HSPs (a
//! long genome query: the loop's helpers then do better than the queries), for a subject searched
//! in chunks in parallel, and on the diagnostic paths (`-verbose`, LOSAT_DEBUG_BLASTN,
//! BLEMIR_DEBUG, LOSAT_DEBUG_COORDS, the trace of `utils/blastn_trace`), which print in the loop's
//! order.
//!
//! With `LOSAT_X_STATS`, `main` prints
//! `[X_STATS] QGROUPPAR prelim= prelim_queries= traceback= traceback_queries= owner=` when a stage
//! was planned: the stages run per query, the query tasks they ran and the stages with HSPs that
//! kept the owner loop by the rule. With `LOSAT_WASI_THREADS_DEBUG=1` each stage run per query
//! prints `[losat-thread-stage] program=blastn stage=prelim_queries|traceback_queries`. In the
//! per-query form the `[TIMING]` steps inside the stages (alignment, purge passes, re-evaluation,
//! final tree) are sums over the threads, while `gapped_extend` and `traceback_prune` stay the
//! stages' wall times.
//!
//! The chunk's ungapped hits are grouped by query in place (a counting sort's destinations, 4
//! bytes per hit, applied by swaps), so the stage holds no second copy of them.
//!
//! Programs that reach it: BLASTN only (megablast, dc-megablast, blastn, blastn-short), from
//! `search_query_batch` (the loop, and a `-db` search of a database of one sequence).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use crate::utils::threading::SearchPool;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

// No NCBI counterpart: reads the switch once (LOSAT_X_QGROUPPAR is a strict-set switch: on
// unless it is `0`). It chooses between two LOSAT code paths with the same result and changes no
// value NCBI computes.
/// Whether the stages may run per query (the switch).
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_QGROUPPAR"))
}

// No NCBI counterpart: counters of the summary lines; they change no value NCBI computes.
static PRELIM: AtomicU64 = AtomicU64::new(0);
static PRELIM_QUERIES: AtomicU64 = AtomicU64::new(0);
static TRACEBACK: AtomicU64 = AtomicU64::new(0);
static TRACEBACK_QUERIES: AtomicU64 = AtomicU64::new(0);
static OWNER: AtomicU64 = AtomicU64::new(0);

/// Which stage a plan is for.
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
            crate::utils::xreach::x_reach!("LOSAT_X_QGROUPPAR", false);
            OWNER.fetch_add(1, Ordering::Relaxed);
            return None;
        }
        Choice::Queries => {
            crate::utils::xreach::x_reach!("LOSAT_X_QGROUPPAR", true);
        }
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
    let _ = pool;
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

/// The line `main` prints at exit: the counters (with `LOSAT_X_STATS`, when a stage was planned).
/// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub(super) fn print_summary() {
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
#[path = "x_qgroup_par_tests.rs"]
mod tests;
