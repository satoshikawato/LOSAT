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
//! Why a query's results do not depend on the other queries: the tree answers per context
//! (`containing_hsp` returns false across contexts, `interval_tree.rs`), the purges sort by
//! context first and compare only within a context (`purge_prelim_hits_with_common_endpoints`,
//! `purge_hsps_with_common_endpoints_ex`), the re-evaluation and the identity test read one HSP
//! and its context's cutoff and statistics, and the sorts are stable. So running the loop of one
//! query over that query's HSPs, in the relative order they have in the batch's list, with a tree
//! of its own, gives the HSPs the batch's loop gives for that query, in the same relative order.
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
//!   the trimmed HSPs and the second purge; the lists are concatenated in query order. The second
//!   purge orders by context first and keeps, within a context, the input's relative order for its
//!   ties, which is the same in both forms (kept HSPs before trimmed ones, each group in the
//!   order of the first purge), so the concatenation is the serial second purge's output; the
//!   score sort, the final tree and the e-values run on it as before.
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
//! `[X_STATS] QGROUPPAR prelim= prelim_queries= traceback= traceback_queries= owner=`, the stages
//! run per query, the query tasks they ran and the stages that kept the owner loop by the rule.
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

/// The order in which the queries with HSPs are handed out (most HSPs first, then query order),
/// or `None` when the owner loop stays: fewer than two queries with HSPs, or one query with more
/// than a quarter of the HSPs.
pub(super) fn plan(counts: &[u32], stage: Stage) -> Option<Vec<u32>> {
    let total: u64 = counts.iter().map(|&count| u64::from(count)).sum();
    let largest = counts.iter().copied().max().unwrap_or(0);
    let with_hsps = counts.iter().filter(|&&count| count > 0).count();
    if with_hsps < 2 || u64::from(largest) * 4 > total {
        OWNER.fetch_add(1, Ordering::Relaxed);
        return None;
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

/// The items of `items` grouped by query (`query_of`), each query's items in their order in
/// `items` (a counting sort); with `starts(counts)` the slice of query `q` is
/// `starts[q]..starts[q + 1]`.
pub(super) fn by_query<T: Clone>(
    items: &[T],
    query_of: impl Fn(&T) -> u32,
    counts: &[u32],
) -> Vec<T> {
    let starts = starts(counts);
    let mut next = starts[..counts.len()].to_vec();
    let mut slots: Vec<Option<T>> = vec![None; items.len()];
    for item in items {
        let query = query_of(item) as usize;
        slots[next[query]] = Some(item.clone());
        next[query] += 1;
    }
    slots
        .into_iter()
        .map(|slot| slot.expect("every slot of a counted query is filled"))
        .collect()
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
    let mut state = new_state();
    for &query in order {
        results[query as usize] = Some(task(&mut state, query));
    }
    results
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
/// `LOSAT_X_STATS`). No NCBI counterpart: it prints counters and changes no value NCBI computes.
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
    if std::env::var_os("LOSAT_X_STATS").is_some() {
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
        let order = plan(&[3, 0, 5, 3, 1], Stage::Prelim).expect("no query above a quarter");
        assert_eq!(order, vec![2, 0, 3, 4]);
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
    fn by_query_keeps_each_querys_order() {
        let items = [(2u32, 'a'), (0, 'b'), (2, 'c'), (1, 'd'), (0, 'e')];
        let counts = counts(items.iter().map(|item| item.0), 4);
        assert_eq!(counts, vec![2, 1, 2, 0]);
        let grouped = by_query(&items, |item| item.0, &counts);
        assert_eq!(
            grouped,
            vec![(0, 'b'), (0, 'e'), (1, 'd'), (2, 'a'), (2, 'c')]
        );
        assert_eq!(starts(&counts), vec![0, 2, 3, 5, 5]);
    }

    #[test]
    fn in_query_order_reads_slices() {
        assert!(in_query_order([0u32, 0, 1, 3, 3].into_iter()));
        assert!(in_query_order(std::iter::empty()));
        assert!(!in_query_order([0u32, 2, 1].into_iter()));
    }

    #[test]
    fn run_gives_every_query_its_own_result() {
        let order = [4u32, 1, 2];
        let serial = run(
            &SearchPool::x_serial(),
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
