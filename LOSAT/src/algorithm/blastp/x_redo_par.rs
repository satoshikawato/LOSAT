//! EXPERIMENT (LOSAT_X_BPREDOPAR): blastp's composition-based redo with jobs of a few matches on
//! the search's pool.
//!
//! NCBI's redo (`Blast_RedoAlignmentCore_MT`, below) reads the matches of all queries from the
//! preliminary stream into one array, in the stream's order (the queries in order, and the HSP
//! lists of a query from best to worst), and splits one loop over that array between its threads.
//! Each thread keeps its own heaps there, so NCBI's output with several threads depends on the
//! thread count; LOSAT gives the one-thread result at every thread count. LOSAT's two parallel
//! paths (`blast_engine.rs`, `run_resolved_in_pool`) run whole queries as jobs, or the matches of
//! the only query that has matches, and a search waits for its largest query: in the blocks of
//! the human x yeast blastp-fast job (`LOSAT_X_BPQBLOCK=20000`) the largest query is 33 % of its
//! block's redo work at the median and 54 % at p90.
//!
//! With the switch on, more than one thread and at least two matches, the matches whose list has
//! HSPs are split between the threads as in NCBI's loop, in jobs of up to four consecutive matches
//! of one query, and the steps whose result depends on order run afterwards on the calling thread,
//! in the order of the one-thread loop:
//! 1. the matches in stream order: for each query, each position of its preliminary hit list whose
//!    list has HSPs (the one-thread loop skips the others before any test); a job is up to four
//!    consecutive of them of one query (one vector of the Newton lanes);
//! 2. windows of consecutive jobs, each computed on the pool (one `pool.install`) and then
//!    replayed. A window holds at most 64 jobs per thread (`LOSAT_X_BPREDOPARWINDOW` sets another
//!    number, for tests) and at most a cap of jobs at risk per query (below), and it may end inside
//!    a query, so it holds at most 4 x 64 x threads results before their replay. A window boundary
//!    changes nothing: the replay is in stream order and the heaps carry over;
//! 3. a job builds its query's workspace (`build_query_workspace`) if no other job of the window has
//!    (one `OnceLock` per query and window, read only afterwards; a query that spans two windows
//!    has it built again from the same inputs, a small cost that grows with the query's length),
//!    and its own parameters
//!    (`build_redo_align_params`), which it uses for its matches one after another as the one-thread
//!    loop uses a query's parameters for all its matches: the `Cell` they hold is set during a call
//!    and set back to its previous value at the end of it;
//! 4. the job prefetches the first Newton problems of its own matches
//!    against the heaps of its query's batch as the window's start left them
//!    (`x_newton_lanes_prefetch_matches`, as the one-thread loop calls it for the next matches);
//!    then each match runs `postprocess_preliminary_hits` with the arguments of the one-thread loop,
//!    with per-worker scratch (composition workspace, gapped-alignment scratch, and the
//!    subject-range SEG cache backed by one shared store for the search; with `LOSAT_X_BPSEGWIN`
//!    the store is emptied after a window once it holds 4,096 entries or more,
//!    `kappa.rs` `x_bp_seg_win`);
//! 5. the replay, per match in order: the early termination test (`BlastCompo_EarlyTermination`)
//!    on the heaps of the query's batch; then the match's error, if it has one (the first error of
//!    a match that the test does not end is the one the one-thread loop reports); then
//!    `BlastCompo_HeapWouldInsert` / `BlastCompo_HeapInsert` on the query's heap. A query's
//!    preliminary list is released as the one-thread loop releases it (`x_prelim_free.rs`), after the replay of the window that holds its last match.
//!
//! The matches the test ends. A heap's `ecutoff` is fixed when it is made; `n` only grows (the redo
//! only inserts; `pop` comes after it); and a heap filled to its cutoff (`n` at least
//! `heapThreshold` and `worstEvalue` at most `ecutoff`) stays so: an insert with an E-value within
//! the cutoff keeps the worst one kept within it, and an insert with a larger E-value does not
//! replace the root, whose E-value is within it (`BlastCompoHeap::insert`). So for a match's fixed
//! `best_evalue` the test is monotone over the replay: once true, true at every later point. Hence:
//! - a match whose test is already true on the heaps at a window's start is not computed (nor
//!   prefetched): its replay finds the test true and drops it (asserted);
//! - the test can be true for a match of query q only for a `best_evalue` above `EVALUE_STRETCH`
//!   times `ecutoff`, and only if every heap of q's batch can be filled by then: the heaps of the
//!   queries before q and q's own hold at most their matches (`n` <= inserts <= matches), and those
//!   of the queries after q are untouched when q's matches are tested (stream order). With the
//!   heaps as the search makes them (`heapThreshold` = `max_target_seqs` >= 1) that is the last
//!   query of a batch whose queries all have at least `max_target_seqs` matches. At a window's
//!   start, also: a query of the batch whose matches are all replayed keeps its heap as it is, so
//!   if that heap is not filled no match of the batch can be ended any more; and q's heap gains at
//!   most one list per match of q, so q's next `heapThreshold` minus `n` matches cannot be ended. A
//!   job with a match that these conditions leave open, not ended at its window's start, is at
//!   risk; the other jobs are never ended and are computed ahead freely;
//! - a computed result is dropped only in a job at risk, in the window during which its batch's
//!   heaps become filled; at the start of any later window the test is already true for the
//!   query's matches that it can end, so at most one window drops results of a query. A window
//!   holds at most cap(q) jobs at risk of a query q, and may hold jobs at risk of several
//!   queries. cap(q) is `threads` at first and after a window that dropped a result of q;
//!   otherwise it grows by q's jobs at risk that the window computed (it doubles when it closed
//!   the window), up to the job budget. So that window drops at most 4 x cap(q) results of q,
//!   with cap(q) <= threads + q's jobs at risk computed, and all used, before it; and where q's
//!   batch's heaps never fill, q's windows are back to their job budget after a number of windows
//!   logarithmic in q's work at risk.
//!
//! Why the output is the one-thread loop's: a match's result is a function of that match's inputs
//! (the same query workspace values, parameter values, subject and matrix), computed by the same
//! sequential function; matches share no reduction, so no floating-point operation changes order.
//! The scratch objects carry no value from one call to the next (NCBI reuses them per thread
//! across any split of the matches, `NRrecord_tld[tid]`). The memos return the values a
//! computation would: the SEG store keeps the masked range a computation gives, and a Newton result
//! prefetched by the job (solved four per vector with the exact solver's operation sequence,
//! `x_newton_lanes.rs`) is kept in the worker thread's store and used only by a call on that
//! thread whose input has the same bits; a call without a stored result solves its problem
//! itself. Which matches are prefetched together changes only which results are stored. Every
//! step that depends on order (the test, the heaps' admission and their internal layout, and which
//! error is reported) is replayed on the same heaps in the same sequence; a match is left out
//! before its replay only where that replay's test is true.
//!
//! `LOSAT_X_BPREDOPAR` is a switch of the strict set (on by default since 2026-10-17): on unless it
//! is `0`. With `LOSAT_TIMING=1` each search the path runs
//! prints
//! `[TIMING] x_bpredopar: matches=<n> jobs=<j> windows=<w> split=<s> ended=<e> skipped=<k> wasted=<x> par_s=<s> replay_s=<s>`
//! (`split`: windows that end inside a query; `ended`: the matches the early termination test
//! ended, as many as the one-thread loop ends; `skipped`: those of them not computed; `wasted`:
//! those computed and dropped).
//!
//! `LOSAT_X_BPREDOLPT` (a switch of the strict set, on by default since 2026-10-17; it acts only
//! where this path runs) computes a window of two or more jobs in one `rayon::scope` of
//! `min(threads, jobs)` tasks instead of a `par_iter` over its jobs: each task holds one set of
//! scratch and takes the window's jobs one at a time from one counter, in the order of their
//! estimated cost (the sum over a job's computed matches of query letters times subject letters,
//! largest first, ties by job index), and stores each job's results in the job's slot
//! (`x_compute_window_lpt`). The window's results are read by job index as without it. With it on,
//! the `[TIMING] x_bpredopar:` line ends with `lpt_windows=<windows computed so>`.
//!
//! This path (graduated 2026-10-27; the switch was `LOSAT_X_BPREDOWINEXPAND`) leaves each compact
//! preliminary list (`x_prelim_compact.rs`) as it is until
//! a window reaches it, instead of expanding the lists of all queries before the redo: before each
//! window the lists from the first one not expanded are expanded in query order, by the same
//! `x_expand_hit_list`, until the queries expanded from the window's start hold more jobs than a
//! window holds, or every list is expanded (`XExpandAhead`). A query's jobs are its matches in
//! jobs of up to four, counted from its compact list as from its expansion (the lists with HSPs;
//! `x_end_above` and the path's choice read the same counts), and a window holds at most its job
//! budget and then reads the next job (one query's), so the window and the first match after it
//! lie within expanded lists (asserted after each window). The windows, the replay and the lists released are
//! therefore those of the full expansion, and the expanded lists held at a time are those of the
//! queries of about one window instead of all queries.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3398-3399,3421,3447-3449

use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};
use std::time::Instant;

use anyhow::Result;
use rayon::prelude::*;

use crate::algorithm::tblastx::lookup::QueryContext;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::BlastCompositionWorkspace;
use crate::stats::search_space::SearchSpace;
use crate::stats::{BlastGumbelBlk, KarlinParams};
use crate::utils::threading::SearchPool;

use super::blast_engine::preliminary_hit_from_local_match_hsp;
use super::gapalign::{BlastpPreliminaryHsp, GapAlignScratch};
use super::hsp::{
    blast_compo_early_termination, BlastCompoHeap, BlastpHitList, BlastpHspList,
    COMPO_HEAP_EVALUE_STRETCH,
};
#[cfg(test)]
use super::hsp::{fill_results_from_compo_heaps, BlastpHsp};
use super::kappa::{
    build_query_workspace, postprocess_preliminary_hits, x_newton_lanes_prefetch_matches,
    BlastRedoAlignParams, BlastpKappaQueryWorkspace, BlastpKappaSubjectRangeCache,
    BlastpPostprocessResult, XSharedSubjectRanges,
};
use super::x_prelim_compact::{x_expand_hit_list, XPrelimHsp};

// No NCBI counterpart: reads the LOSAT_X_BPREDOPAR switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPREDOPAR): a strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_redo_par() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPREDOPAR"))
}

// No NCBI counterpart: reads the LOSAT_X_BPREDOLPT switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPREDOLPT): a strict-set switch, on unless it is `0`.
/// It acts only where the redo runs
/// in windows (`LOSAT_X_BPREDOPAR`).
pub(crate) fn x_bp_redo_lpt() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPREDOLPT"))
}

/// Jobs per thread that a window holds at most, unless `LOSAT_X_BPREDOPARWINDOW` gives another
/// number.
const X_WINDOW_JOBS_PER_THREAD: usize = 64;

// No NCBI counterpart: reads the LOSAT_X_BPREDOPARWINDOW tuning variable once, when the path
// first runs; it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDOPARWINDOW): the jobs per thread of a window when the variable is a
/// positive integer (tests make small windows with it), `X_WINDOW_JOBS_PER_THREAD` when it is
/// not set or holds anything else (0, a negative number, not a number), as a value switch of
/// `utils/xswitch.rs` gives its built-in value.
fn x_window_jobs_per_thread() -> usize {
    static JOBS: OnceLock<usize> = OnceLock::new();
    *JOBS.get_or_init(|| {
        std::env::var("LOSAT_X_BPREDOPARWINDOW")
            .ok()
            .and_then(|value| value.trim().parse::<usize>().ok())
            .filter(|&jobs| jobs > 0)
            .unwrap_or(X_WINDOW_JOBS_PER_THREAD)
    })
}

// No NCBI counterpart: matches per job, one vector of the Newton lanes (four problems); it
// changes no value NCBI computes.
/// Matches of one query per job at most.
const X_JOB_MATCHES: usize = 4;

/// What the redo of a search reads besides the preliminary hit lists and the heaps: the values
/// that `redo_run::one_thread_redo` reads (made of `redo_run::RedoInputs` in `x_windows_redo`).
pub(crate) struct XRedoParSearch<'a, P> {
    pub pool: &'a SearchPool<'a>,
    pub threads: usize,
    pub contexts: &'a [QueryContext],
    /// `blast_engine::query_nomask_sequence`.
    pub query_nomask: fn(&QueryContext) -> &[u8],
    pub search_spaces: &'a [SearchSpace],
    /// The encoded subject set, or a database's listed subjects, read per match.
    pub subjects: super::x_redo_subjects::XRedoSubjects<'a>,
    /// Each query's batch (`kappa_query_batches`): the heaps its early termination test reads.
    pub query_batches: &'a [Range<usize>],
    /// Each query's context of the E-value adjustment with a composition P-value
    /// (`kappa::evalue_adjust_contexts`).
    pub evalue_adjust_contexts: &'a [super::kappa::EvalueAdjustContext],
    pub matrix: ScoringMatrix,
    pub gapped_params: &'a KarlinParams,
    pub gapped_gumbel: Option<&'a BlastGumbelBlk>,
    pub evalue: f64,
    /// `build_redo_align_params` of `run_call` (`RedoInputs`).
    pub build_redo_align_params: &'a P,
    pub timing: bool,
    /// `postprocess_preliminary_hits`' `engine_positives` (the Seq-aligns are requested).
    pub engine_positives: bool,
}

/// One match: its query and the position of its list in the query's preliminary hit list.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct XMatch {
    q_idx: usize,
    position: usize,
}

/// The scratch of one worker.
struct XScratch {
    composition_workspace: BlastCompositionWorkspace,
    gap_scratch: GapAlignScratch,
    hits: Vec<BlastpPreliminaryHsp>,
    range_cache: BlastpKappaSubjectRangeCache,
}

impl XScratch {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334
    fn new(range_cache: BlastpKappaSubjectRangeCache, matrix: ScoringMatrix) -> Self {
        let mut scratch = Self {
            composition_workspace: BlastCompositionWorkspace::new_blosum62(),
            gap_scratch: GapAlignScratch::new(),
            hits: Vec::new(),
            range_cache,
        };
        // Another matrix than BLOSUM62: its workspace (item SC2b).
        if matrix != ScoringMatrix::Blosum62 {
            super::blast_engine::blastp_use_matrix_workspace(
                &mut scratch.composition_workspace,
                matrix,
            );
        }
        scratch
    }
}

#[derive(Debug, Default)]
struct XCounts {
    matches: usize,
    jobs: usize,
    windows: usize,
    /// Windows that end inside a query (the next window continues it).
    split: usize,
    /// Matches the early termination test ended that were not computed (ended before their
    /// window).
    skipped: usize,
    /// Matches the early termination test ended that were computed.
    wasted: usize,
    released: usize,
    par_s: f64,
    replay_s: f64,
    /// With `LOSAT_X_BPREDOLPT`: the windows computed in one scope with jobs largest first.
    lpt_windows: Option<usize>,
}

/// The matches of one query's preliminary hit list: its lists with HSPs. A compact list
/// (`x_prelim_compact.rs`) has the lists, and the HSPs per list, of its expansion.
fn x_list_matches<H>(list: &BlastpHitList<H>) -> usize {
    list.hsplist_array
        .iter()
        .take(list.hsplist_count)
        .filter(|local_match| !local_match.hsps.is_empty())
        .count()
}

/// The matches of each query: those of its compact list while it is not expanded
/// (`compact` is empty or has one entry per query), else those of its
/// preliminary hit list.
fn x_query_matches(
    lists: &[Option<BlastpHitList>],
    compact: &[Option<BlastpHitList<XPrelimHsp>>],
) -> Vec<usize> {
    (0..lists.len())
        .map(|q_idx| match compact.get(q_idx).and_then(Option::as_ref) {
            Some(list) => x_list_matches(list),
            None => lists[q_idx].as_ref().map_or(0, x_list_matches),
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3398-3399,3511
/// The matches that are redone: those whose list has HSPs, counted from a query's compact list
/// while it is not expanded (`x_query_matches`).
pub(crate) fn x_matches(
    lists: &[Option<BlastpHitList>],
    compact: &[Option<BlastpHitList<XPrelimHsp>>],
) -> usize {
    x_query_matches(lists, compact).into_iter().sum()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1576
// No NCBI counterpart for the bound: it decides only which results are computed ahead; the
// replay runs the test for every match.
/// For each query, the E-value above which the early termination test may end one of its matches
/// during the replay into `heaps` (the heaps before the replay), or `None` when it ends none of
/// them; `matches` holds each query's matches (`x_query_matches`). The test is true only when
/// every heap of the batch is filled to its cutoff and the E-value exceeds `EVALUE_STRETCH` times
/// each heap's cutoff. In stream order a query's matches are tested after those of the queries
/// before it in its batch and before those of the queries after it, so when they are tested the
/// heaps up to the query hold at most `n` plus their queries' matches, and those after it are as
/// they are now.
fn x_end_above(
    matches: &[usize],
    heaps: &[BlastCompoHeap],
    query_batches: &[Range<usize>],
) -> Vec<Option<f64>> {
    (0..matches.len())
        .map(|q_idx| {
            let batch = query_batches[q_idx].clone();
            // The heaps after the query first: they are not filled unless the query is the
            // last one of its batch, so this ends the test at once for the other queries.
            let may_end = (q_idx + 1..batch.end).all(|other| heaps[other].filled_to_cutoff())
                && (batch.start..=q_idx).all(|other| heaps[other].may_fill(matches[other]));
            may_end.then(|| {
                batch
                    .map(|other| COMPO_HEAP_EVALUE_STRETCH * heaps[other].ecutoff)
                    .fold(f64::NEG_INFINITY, f64::max)
            })
        })
        .collect()
}

/// The bounds of the windows of a search: their jobs, and the jobs at risk of a query that a
/// window holds at first (those with a match that the early termination test may end and had not
/// ended at the window's start; `x_window`), which `x_redo_windows` adapts per query.
#[derive(Clone, Copy, Debug)]
struct XWindowLimits {
    jobs: usize,
    at_risk_jobs: usize,
}

/// A window: consecutive jobs of the stream, computed together and then replayed.
#[derive(Debug, PartialEq)]
struct XWindow {
    /// Its matches in stream order.
    matches: Vec<XMatch>,
    /// Per match: the early termination test was true for it on the heaps at the window's start,
    /// so it is not computed.
    ended: Vec<bool>,
    /// Its jobs as ranges of `matches`.
    jobs: Vec<Range<usize>>,
    /// The first match after the window, where the next one starts (`None`: no match is left).
    next: Option<XMatch>,
    /// Its jobs at risk: (query, jobs) for each query that has any, in stream order.
    at_risk: Vec<(usize, usize)>,
}

/// The results of a window's matches in order; `None` for a match not computed.
type XWindowResults = Vec<Option<Result<BlastpPostprocessResult>>>;

/// What the early termination test ended in a window's replay.
#[derive(Debug, Default, PartialEq, Eq)]
struct XReplayed {
    /// Matches without a result (not computed).
    skipped: usize,
    /// Matches with a result (computed and dropped).
    wasted: usize,
    /// The queries of the dropped results, each once.
    dropped: Vec<usize>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3398-3399,3421
/// The next job from `cursor` (the next list position to read) into `job`: up to
/// `X_JOB_MATCHES` consecutive matches of one query in stream order, none when no match is
/// left. `cursor` moves past them.
fn x_next_job(lists: &[Option<BlastpHitList>], cursor: &mut XMatch, job: &mut Vec<XMatch>) {
    job.clear();
    while cursor.q_idx < lists.len() {
        if let Some(list) = lists[cursor.q_idx].as_ref() {
            while cursor.position < list.hsplist_count {
                if !list.hsplist_array[cursor.position].hsps.is_empty() {
                    if job.len() == X_JOB_MATCHES {
                        return;
                    }
                    job.push(*cursor);
                }
                cursor.position += 1;
            }
        }
        if !job.is_empty() {
            return;
        }
        *cursor = XMatch {
            q_idx: cursor.q_idx + 1,
            position: 0,
        };
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529
/// The window that starts at `start`, on the heaps as the replay of the windows before left them:
/// the jobs from there while it holds fewer than `jobs` jobs and, for a job at risk of query q,
/// fewer than `caps[q]` jobs at risk of q. A match whose test is already true is marked ended;
/// the test stays true until its replay. A match may be ended later only if all of these hold
/// (each is necessary for the test to be true at its replay; a job with such a match, not ended,
/// is at risk):
/// - its E-value is above `end_above` of its query (`x_end_above`, from the heaps before the
///   replay);
/// - no query of its batch whose matches are all replayed (the queries before the window's first)
///   has a heap that is not filled to its cutoff: such a heap receives no insert any more;
/// - its query's heap can be filled by then: it holds `n` now and gains at most one list per match
///   of the query that the replay redoes before this one, so the query's next
///   `heapThreshold` minus `n` matches of the window cannot be ended.
fn x_window(
    lists: &[Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    query_batches: &[Range<usize>],
    end_above: &[Option<f64>],
    caps: &[usize],
    start: XMatch,
    jobs: usize,
) -> XWindow {
    let mut window = XWindow {
        matches: Vec::new(),
        ended: Vec::new(),
        jobs: Vec::new(),
        next: None,
        at_risk: Vec::new(),
    };
    let mut cursor = start;
    let mut job = Vec::with_capacity(X_JOB_MATCHES);
    // The query of the last job: whether its batch can still end matches, how many of its
    // matches the window holds, and how many of its jobs at risk.
    let mut query: Option<(usize, bool, usize, usize)> = None;
    loop {
        x_next_job(lists, &mut cursor, &mut job);
        let Some(&first) = job.first() else {
            return window;
        };
        let q_idx = first.q_idx;
        let (batch_open, in_window, at_risk_jobs) = match query {
            Some((current, batch_open, in_window, at_risk_jobs)) if current == q_idx => {
                (batch_open, in_window, at_risk_jobs)
            }
            _ => {
                let batch = query_batches[q_idx].clone();
                let open = end_above[q_idx].is_some()
                    && (batch.start..batch.end.min(start.q_idx))
                        .all(|done| heaps[done].filled_to_cutoff());
                (open, 0, 0)
            }
        };
        let heap = &heaps[q_idx];
        let mark = window.matches.len();
        let mut at_risk = false;
        for (index, &matched) in job.iter().enumerate() {
            let best_evalue = x_local_match(lists, matched).best_evalue;
            let may_end = batch_open
                && end_above[q_idx].is_some_and(|above| best_evalue > above)
                && heap.may_fill(in_window + index);
            let ended = may_end
                && blast_compo_early_termination(best_evalue, &heaps[query_batches[q_idx].clone()]);
            at_risk |= may_end && !ended;
            window.ended.push(ended);
        }
        // A window holds at least one job.
        if !window.jobs.is_empty()
            && (window.jobs.len() >= jobs || (at_risk && at_risk_jobs >= caps[q_idx]))
        {
            window.ended.truncate(mark);
            window.next = Some(first);
            return window;
        }
        if at_risk {
            match window.at_risk.last_mut() {
                Some((current, count)) if *current == q_idx => *count += 1,
                _ => window.at_risk.push((q_idx, 1)),
            }
        }
        query = Some((
            q_idx,
            batch_open,
            in_window + job.len(),
            at_risk_jobs + usize::from(at_risk),
        ));
        window.matches.extend_from_slice(&job);
        window.jobs.push(mark..window.matches.len());
    }
}

fn x_local_match(lists: &[Option<BlastpHitList>], matched: XMatch) -> &BlastpHspList {
    &lists[matched.q_idx]
        .as_ref()
        .expect("a match's query has its preliminary hit list")
        .hsplist_array[matched.position]
}

/// The compact preliminary lists (`x_prelim_compact.rs`) that
/// the windows expand as they reach them, in query order.
pub(crate) struct XExpandAhead<'a> {
    /// A query's compact list until it is expanded into its preliminary hit list.
    compact: &'a mut [Option<BlastpHitList<XPrelimHsp>>],
    /// Each query's jobs: its matches in jobs of up to `X_JOB_MATCHES` (`x_next_job`).
    jobs: Vec<usize>,
    /// The queries before it have their lists expanded.
    expanded_to: usize,
    /// The jobs from the next window's start to `expanded_to`.
    ahead: usize,
}

impl<'a> XExpandAhead<'a> {
    /// The expansion of `compact` ahead of the windows of the redo of `lists`.
    pub(crate) fn new(
        lists: &[Option<BlastpHitList>],
        compact: &'a mut [Option<BlastpHitList<XPrelimHsp>>],
    ) -> Self {
        debug_assert_eq!(lists.len(), compact.len());
        Self {
            compact,
            jobs: Vec::new(),
            expanded_to: 0,
            ahead: 0,
        }
    }

    /// Before the windows: the jobs of each query from its matches (`x_query_matches`).
    fn start(&mut self, matches: &[usize]) {
        self.jobs = matches
            .iter()
            .map(|&count| count.div_ceil(X_JOB_MATCHES))
            .collect();
        self.expanded_to = 0;
        self.ahead = 0;
    }

    /// Before a window: expands the lists from the first one not expanded until the queries
    /// expanded from the window's start hold more than `jobs` jobs (the window's budget and the
    /// first job after it) or every list is expanded.
    fn expand(&mut self, lists: &mut [Option<BlastpHitList>], jobs: usize) {
        while self.expanded_to < lists.len() && self.ahead <= jobs {
            let q_idx = self.expanded_to;
            // As the full expansion: the compact list replaces the query's list.
            if let Some(list) = self.compact[q_idx].take() {
                lists[q_idx] = Some(x_expand_hit_list(list));
            }
            self.ahead += self.jobs[q_idx];
            self.expanded_to += 1;
        }
    }

    // No NCBI counterpart: a check of the windows' bounds; it changes no value NCBI computes.
    /// After `window` (from `start`) was chosen: it and the first match after it lie within the
    /// expanded lists.
    fn check(&mut self, lists: &[Option<BlastpHitList>], start: XMatch, window: &XWindow) {
        assert!(
            window.next.map_or(self.expanded_to == lists.len(), |next| {
                next.q_idx < self.expanded_to
            }),
            "a window from {start:?} reads past the expanded lists (next {:?}, {} of {} expanded)",
            window.next,
            self.expanded_to,
            lists.len()
        );
        self.ahead -= window.jobs.len();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2308-2336
/// The query workspace of the one-thread loop for query `q_idx`.
fn x_query_workspace<P>(search: &XRedoParSearch<'_, P>, q_idx: usize) -> BlastpKappaQueryWorkspace {
    let ctx = &search.contexts[q_idx];
    let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
    build_query_workspace(
        query_raw,
        ctx.aa_len as i32,
        search.search_spaces[q_idx].length_adjustment as i32,
        search.search_spaces[q_idx].effective_space,
        search.evalue_adjust_contexts[q_idx],
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3542,3624-3628
/// The redo of one match: the call of the one-thread loop, with the same arguments.
fn x_redo_match<P>(
    search: &XRedoParSearch<'_, P>,
    q_idx: usize,
    local_match: &BlastpHspList,
    query_workspace: &BlastpKappaQueryWorkspace,
    redo_align_params: &BlastRedoAlignParams,
    scratch: &mut XScratch,
) -> Result<BlastpPostprocessResult> {
    scratch.hits.clear();
    scratch.hits.extend(
        local_match
            .hsps
            .iter()
            .map(preliminary_hit_from_local_match_hsp),
    );
    let subject = search.subjects.read(local_match.oid)?;
    let subject_raw = subject.residues();
    postprocess_preliminary_hits(
        &scratch.hits,
        query_workspace,
        &mut scratch.composition_workspace,
        (search.query_nomask)(&search.contexts[q_idx]),
        subject_raw,
        search.matrix,
        search.gapped_params,
        search.gapped_gumbel,
        redo_align_params,
        search.evalue,
        &mut scratch.gap_scratch,
        &mut scratch.range_cache,
        search.engine_positives,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3447-3449,3493-3511
/// The redo of a window's matches on the pool, in jobs of up to `X_JOB_MATCHES` consecutive
/// matches of one query; the results in match order, `None` for a match ended before the window
/// (not computed). `heaps` are the heaps as the window's start left them; nothing changes them
/// while the jobs run. NCBI's loop also runs the early termination test and the heap insertion of
/// each match on its thread; here they are replayed afterwards (`x_replay`). With `lpt`
/// (`LOSAT_X_BPREDOLPT`), more than one thread on the search's pool and two or more jobs, the
/// jobs run as `x_compute_window_lpt` runs them; the second value says whether they did.
fn x_compute_window<P>(
    search: &XRedoParSearch<'_, P>,
    lists: &[Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    window: &XWindow,
    shared: &Arc<XSharedSubjectRanges>,
    lpt: bool,
) -> (XWindowResults, bool)
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let (Some(first), Some(last)) = (window.matches.first(), window.matches.last()) else {
        return (Vec::new(), false);
    };
    if window.ended.iter().all(|&ended| ended) {
        return (window.ended.iter().map(|_| None).collect(), false);
    }
    let queries = first.q_idx..last.q_idx + 1;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3244-3248
    // Built by the first job of the query that needs it; the others wait for the same value.
    let workspaces: Vec<OnceLock<BlastpKappaQueryWorkspace>> =
        queries.clone().map(|_| OnceLock::new()).collect();
    // One job's results in match order, with one set of scratch; both ways of running the jobs
    // below call it.
    let run_job = |scratch: &mut XScratch, job: &Range<usize>| -> XWindowResults {
        let job_matches = &window.matches[job.clone()];
        let job_ended = &window.ended[job.clone()];
        if job_ended.iter().all(|&ended| ended) {
            return job_ended.iter().map(|_| None).collect();
        }
        let q_idx = job_matches[0].q_idx;
        let query_workspace =
            workspaces[q_idx - queries.start].get_or_init(|| x_query_workspace(search, q_idx));
        // The parameters of the one-thread loop for this query, used for the job's matches one
        // after another as that loop uses them for all of the query's.
        let redo_align_params = match (search.build_redo_align_params)(q_idx) {
            Ok(params) => params,
            Err(error) => {
                return job_ended
                    .iter()
                    .map(|&ended| (!ended).then(|| Err(anyhow::anyhow!("{error:#}"))))
                    .collect();
            }
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529
        // No NCBI counterpart: the Newton lanes solve the first Newton problems of the job's
        // matches, four per vector (kappa::x_newton_lanes_prefetch_matches); a result is used
        // only by a call on this thread with the same input bits. It does not change any value
        // NCBI computes. Given the heaps of the batch as the window's start left them, as the
        // one-thread loop gives its heaps, it leaves out the matches the test ends on them: those
        // not computed here.
        {
            let list = lists[q_idx]
                .as_ref()
                .expect("a match's query has its preliminary hit list");
            x_newton_lanes_prefetch_matches(
                &list.hsplist_array
                    [job_matches[0].position..=job_matches[job_matches.len() - 1].position],
                search.subjects,
                &heaps[search.query_batches[q_idx].clone()],
                query_workspace,
                &scratch.composition_workspace,
                &redo_align_params,
                &mut scratch.range_cache,
            );
        }
        job_matches
            .iter()
            .zip(job_ended)
            .map(|(&matched, &ended)| {
                (!ended).then(|| {
                    x_redo_match(
                        search,
                        q_idx,
                        x_local_match(lists, matched),
                        query_workspace,
                        &redo_align_params,
                        scratch,
                    )
                })
            })
            .collect()
    };
    let new_scratch = || {
        XScratch::new(
            BlastpKappaSubjectRangeCache::new(Some(shared.clone())),
            search.matrix,
        )
    };
    let lpt = lpt && {
        let acts = search.pool.enabled() && search.threads > 1 && window.jobs.len() > 1;
        crate::utils::xreach::x_reach!("LOSAT_X_BPREDOLPT", acts);
        acts
    };
    let per_job: Vec<XWindowResults> = if lpt {
        let order = x_lpt_order(&x_job_costs(search, lists, window));
        x_compute_window_lpt(
            search.pool,
            search.threads.min(window.jobs.len()),
            &order,
            new_scratch,
            |scratch, job| run_job(scratch, &window.jobs[job]),
        )
    } else {
        search.pool.install(|| {
            window
                .jobs
                .par_iter()
                // Each rayon work unit owns its scratch; how many are made depends on the
                // scheduler, and no value depends on it.
                .map_init(new_scratch, run_job)
                .collect()
        })
    };
    (per_job.into_iter().flatten().collect(), lpt)
}

// No NCBI counterpart: an estimate that decides only the order in which the threads take a
// window's jobs (`x_lpt_order`); it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDOLPT): the estimated cost of each job of `window`: the sum over its
/// matches that it computes (those not ended before the window) of the query's letters times the
/// subject's letters.
fn x_job_costs<P>(
    search: &XRedoParSearch<'_, P>,
    lists: &[Option<BlastpHitList>],
    window: &XWindow,
) -> Vec<u64> {
    window
        .jobs
        .iter()
        .map(|job| {
            window.matches[job.clone()]
                .iter()
                .zip(&window.ended[job.clone()])
                .filter(|&(_, &ended)| !ended)
                .map(|(&matched, _)| {
                    let query = search.contexts[matched.q_idx].aa_len as u64;
                    let subject = search
                        .subjects
                        .aa_len(x_local_match(lists, matched).oid as usize)
                        .map_or(0, |length| length as u64);
                    query.saturating_mul(subject)
                })
                .fold(0, u64::saturating_add)
        })
        .collect()
}

// No NCBI counterpart: the order in which the threads take a window's jobs; it changes no value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDOLPT): the job indices by `costs`, largest first; jobs of equal cost
/// in index order.
fn x_lpt_order(costs: &[u64]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..costs.len()).collect();
    order.sort_by_key(|&job| (std::cmp::Reverse(costs[job]), job));
    order
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3447-3449,3493-3497
// No NCBI counterpart for the order: NCBI's threads take fixed shares of its loop; here the
// threads take the jobs one at a time, largest first. Scheduling only.
/// EXPERIMENT (LOSAT_X_BPREDOLPT): `run_job(scratch, job)` for every job index of `order`, run
/// in one `rayon::scope` on `pool` by `tasks` tasks. Each task takes the next index of `order`
/// from one counter until none is left, with one scratch (`new_scratch`, made for its first job);
/// each job's result goes to the job's own slot. Returns the results by job index.
///
/// Why the window's results are the `par_iter` path's: a job's results are a function of its
/// matches (`run_job` is the same function with the same arguments); the query workspaces it
/// reads are per-query `OnceLock`s of the same values, whichever job builds one; every scratch
/// carries no value from one call to the next, so it does not matter how many jobs a task runs
/// with it, and the memos (the SEG store, a thread's prefetched Newton results) return the values
/// a computation would (module header; the `par_iter` path also runs several jobs with one
/// scratch, as many as the scheduler gives a work unit); the results are read by job index, so
/// `x_replay` sees the same sequence, errors included (each job's error is its own value in its
/// slot; the replay reports the first one in stream order, as from the `par_iter` path's
/// results). The order only changes which thread computes what, and when. No floating-point
/// operation changes order: jobs share no reduction. A task waits for nothing but a workspace
/// another job is building, and that building runs no rayon work, so no task depends on a task
/// that could be stolen below it.
fn x_compute_window_lpt<S, R>(
    pool: &SearchPool<'_>,
    tasks: usize,
    order: &[usize],
    new_scratch: impl Fn() -> S + Sync,
    run_job: impl Fn(&mut S, usize) -> R + Sync,
) -> Vec<R>
where
    R: Send + Sync,
{
    // A disabled pool would run the scope on the global pool, not the search pool.
    debug_assert!(pool.enabled());
    let slots: Vec<OnceLock<R>> = order.iter().map(|_| OnceLock::new()).collect();
    let next = AtomicUsize::new(0);
    pool.install(|| {
        rayon::scope(|scope| {
            for _ in 0..tasks {
                scope.spawn(|_| {
                    let mut scratch = None;
                    while let Some(&job) = order.get(next.fetch_add(1, Ordering::Relaxed)) {
                        let result = run_job(scratch.get_or_insert_with(&new_scratch), job);
                        assert!(slots[job].set(result).is_ok(), "job {job} ran twice");
                    }
                });
            }
        });
    });
    slots
        .into_iter()
        .enumerate()
        .map(|(job, slot)| {
            slot.into_inner()
                .unwrap_or_else(|| panic!("job {job} did not run"))
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529,3705-3719
/// The steps of the one-thread loop that depend on order, for a window's matches in stream order,
/// with the results computed for them: the early termination test on the heaps of the query's
/// batch, the match's error, and the heap insertion.
fn x_replay(
    matches: &[XMatch],
    results: XWindowResults,
    lists: &[Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    query_batches: &[Range<usize>],
) -> Result<XReplayed> {
    assert_eq!(matches.len(), results.len(), "one result per match");
    let mut replayed = XReplayed::default();
    for (&matched, postprocessed) in matches.iter().zip(results) {
        let local_match = x_local_match(lists, matched);
        if blast_compo_early_termination(
            local_match.best_evalue,
            &heaps[query_batches[matched.q_idx].clone()],
        ) {
            if postprocessed.is_some() {
                replayed.wasted += 1;
                if replayed.dropped.last() != Some(&matched.q_idx) {
                    replayed.dropped.push(matched.q_idx);
                }
            } else {
                replayed.skipped += 1;
            }
            continue;
        }
        // A match without a result had the test true at its window's start; the test is
        // monotone over the replay (`x_end_above`, header), so it is true here.
        let Some(postprocessed) = postprocessed else {
            unreachable!(
                "query {} match {}: ended before its window, not ended at its replay",
                matched.q_idx, matched.position
            );
        };
        let postprocessed = postprocessed?;
        let best_score = postprocessed.best_score;
        let best_evalue = postprocessed.best_evalue;
        let Some(hsp_list) = postprocessed.hsp_list else {
            continue;
        };
        let redone_match = &mut heaps[matched.q_idx];
        if redone_match.would_insert(best_evalue, best_score, local_match.oid) {
            redone_match.insert(hsp_list, best_evalue, best_score, local_match.oid);
        }
    }
    Ok(replayed)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3398-3399,3447-3449,3735-3738
/// The windows of a search, each computed (`compute`, with the lists and the heaps as the replay
/// of the windows before left them) and then replayed in stream order into `heaps`. A query's
/// preliminary list is released after the replay of the window that holds its last match. With
/// `expand`, the compact lists are expanded ahead of each window.
#[allow(clippy::too_many_arguments)]
fn x_redo_windows(
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    query_batches: &[Range<usize>],
    limits: XWindowLimits,
    timing: bool,
    mut expand: Option<&mut XExpandAhead<'_>>,
    compute: &mut dyn FnMut(
        &[Option<BlastpHitList>],
        &[BlastCompoHeap],
        &XWindow,
    ) -> XWindowResults,
) -> Result<XCounts> {
    let matches = x_query_matches(
        lists,
        expand.as_deref().map_or(&[][..], |expand| &*expand.compact),
    );
    let end_above = x_end_above(&matches, heaps, query_batches);
    if let Some(expand) = expand.as_deref_mut() {
        expand.start(&matches);
    }
    let mut counts = XCounts::default();
    let mut start = Some(XMatch {
        q_idx: 0,
        position: 0,
    });
    let mut released_to = 0;
    // The cap on the jobs at risk of query q in a window, `caps[q]`: `limits.at_risk_jobs` at
    // first and after a window that dropped a computed result of q; otherwise grown by q's jobs
    // at risk that the window computed (all used; it doubles when it closed the window), up to
    // the job budget. So it is at most `limits.at_risk_jobs` plus q's jobs at risk computed and
    // used before.
    let budget = limits.jobs.max(limits.at_risk_jobs);
    let mut caps = vec![limits.at_risk_jobs; lists.len()];
    while let Some(cursor) = start {
        if let Some(expand) = expand.as_deref_mut() {
            expand.expand(lists, limits.jobs);
        }
        let window = x_window(
            lists,
            heaps,
            query_batches,
            &end_above,
            &caps,
            cursor,
            limits.jobs,
        );
        if let Some(expand) = expand.as_deref_mut() {
            expand.check(lists, cursor, &window);
        }
        if !window.jobs.is_empty() {
            counts.windows += 1;
            counts.matches += window.matches.len();
            counts.jobs += window.jobs.len();
            if window
                .next
                .zip(window.matches.last())
                .is_some_and(|(next, last)| next.q_idx == last.q_idx)
            {
                counts.split += 1;
            }
            let par_start = timing.then(Instant::now);
            let results = compute(lists, heaps, &window);
            if let Some(par_start) = par_start {
                counts.par_s += par_start.elapsed().as_secs_f64();
            }
            let replay_start = timing.then(Instant::now);
            let replayed = x_replay(&window.matches, results, lists, heaps, query_batches)?;
            counts.skipped += replayed.skipped;
            counts.wasted += replayed.wasted;
            if let Some(replay_start) = replay_start {
                counts.replay_s += replay_start.elapsed().as_secs_f64();
            }
            for &(q_idx, at_risk_jobs) in &window.at_risk {
                caps[q_idx] = caps[q_idx].saturating_add(at_risk_jobs).min(budget);
            }
            // A dropped result is of a job at risk; after its window the pre-check leaves out
            // every match of the query that the test ends.
            for &q_idx in &replayed.dropped {
                caps[q_idx] = limits.at_risk_jobs;
            }
        }
        // Nothing reads a query's preliminary
        // list after its redo; the one-thread loop releases each one after its query
        // (`x_prelim_free.rs`). The queries before the next window's first match are done.
        let done = window.next.map_or(lists.len(), |next| next.q_idx);
        for list in &mut lists[released_to..done] {
            if list.take().is_some() {
                counts.released += 1;
            }
        }
        released_to = done;
        start = window.next;
    }
    Ok(counts)
}

/// The redo of a search: windows of jobs computed on the pool and replayed in stream order into
/// `heaps`; each query's preliminary list is released after its redo; with `expand`, the compact
/// lists are expanded ahead of the windows.
fn x_redo<P>(
    search: &XRedoParSearch<'_, P>,
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    expand: Option<&mut XExpandAhead<'_>>,
) -> Result<XCounts>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    // The one-thread loop builds the parameters of query 0 before any match. Their only error
    // (the initial matrix, kept once for the search) is the same for every query, so it is
    // reported here, before any match, as there.
    (search.build_redo_align_params)(0)?;
    // One store of masked subject ranges for the search, as in the query-parallel redo
    // (`XSharedSubjectRanges`); in a block search of the subjects encoded for all blocks the store
    // of all blocks (`x_query_blocks::x_block_seg_store`): an entry depends on its key only (the
    // subject's index in that encoding, the range), the same in every block.
    let block_store = super::x_query_blocks::x_block_seg_store(search.pool, search.subjects);
    let shared = block_store
        .clone()
        .unwrap_or_else(|| Arc::new(XSharedSubjectRanges::new()));
    let threads = search.threads.max(1);
    let limits = XWindowLimits {
        jobs: x_window_jobs_per_thread().saturating_mul(threads),
        at_risk_jobs: threads,
    };
    // EXPERIMENT (LOSAT_X_BPREDOLPT): the windows' jobs in one scope, largest first.
    let lpt = x_bp_redo_lpt();
    // EXPERIMENT (LOSAT_X_BPSEGWIN): a store of `seg_win` entries or more after a window is
    // emptied.
    let seg_win = super::kappa::x_bp_seg_win();
    let mut seg_emptied = false;
    let mut lpt_windows = 0;
    let counts = x_redo_windows(
        lists,
        heaps,
        search.query_batches,
        limits,
        search.timing,
        expand,
        &mut |lists, heaps, window| {
            let (results, lpt_ran) = x_compute_window(search, lists, heaps, window, &shared, lpt);
            if seg_win.is_some_and(|min| shared.x_len() >= min) {
                shared.x_clear();
                seg_emptied = true;
            }
            lpt_windows += usize::from(lpt_ran);
            results
        },
    );
    if seg_win.is_some() {
        crate::utils::xreach::x_reach!("LOSAT_X_BPSEGWIN", seg_emptied);
    }
    // The block search keeps the store for the next block only as LOSAT_X_BPSEGWIN keeps it for
    // the next window (fewer entries than its bound); without the switch a search's store ends with
    // the search, as outside the block search.
    if block_store.is_some() && seg_win.is_none_or(|min| shared.x_len() >= min) {
        shared.x_clear();
    }
    let mut counts = counts?;
    counts.lpt_windows = lpt.then_some(lpt_windows);
    Ok(counts)
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
fn x_print_timing(counts: &XCounts) {
    eprintln!(
        "[TIMING] x_bpredopar: matches={} jobs={} windows={} split={} ended={} skipped={} wasted={} par_s={:.6} replay_s={:.6}{}",
        counts.matches,
        counts.jobs,
        counts.windows,
        counts.split,
        counts.skipped + counts.wasted,
        counts.skipped,
        counts.wasted,
        counts.par_s,
        counts.replay_s,
        counts
            .lpt_windows
            .map_or(String::new(), |windows| format!(" lpt_windows={windows}"))
    );
}

/// EXPERIMENT (LOSAT_X_BPREDOPAR): the redo of a search in jobs of a few matches, into `heaps`
/// (the redone heaps of the search). The preliminary lists are released after their redo, as the
/// one-thread loop releases them. With `compact` (one
/// entry per query), a query's compact list stands for its list until the windows expand it; none
/// is left afterwards. Returns the number of lists released.
pub(crate) fn x_redo_par_run<P>(
    search: &XRedoParSearch<'_, P>,
    lists: &mut [Option<BlastpHitList>],
    compact: Option<&mut [Option<BlastpHitList<XPrelimHsp>>]>,
    heaps: &mut [BlastCompoHeap],
) -> Result<usize>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let mut expand = compact.map(|compact| XExpandAhead::new(lists, compact));
    let counts = x_redo(search, lists, heaps, expand.as_mut())?;
    if search.timing {
        x_print_timing(&counts);
    }
    Ok(counts.released)
}

// No NCBI counterpart: the comparison of two results of the same match (floating-point values by
// their bits), for the unit test that compares the windows' computations.
#[cfg(test)]
fn x_hsp_difference(a: &BlastpHsp, b: &BlastpHsp) -> bool {
    !(a.identity.to_bits() == b.identity.to_bits()
        && a.length == b.length
        && a.mismatch == b.mismatch
        && a.gapopen == b.gapopen
        && a.q_start == b.q_start
        && a.q_end == b.q_end
        && a.s_start == b.s_start
        && a.s_end == b.s_end
        && a.gapped_q_start == b.gapped_q_start
        && a.gapped_s_start == b.gapped_s_start
        && a.e_value.to_bits() == b.e_value.to_bits()
        && a.bit_score.to_bits() == b.bit_score.to_bits()
        && a.num_ident == b.num_ident
        && a.query_context == b.query_context
        && a.query_frame == b.query_frame
        && a.subject_frame == b.subject_frame
        && a.query_length == b.query_length
        && a.q_idx == b.q_idx
        && a.s_idx == b.s_idx
        && a.raw_score == b.raw_score
        && a.gap_info == b.gap_info
        && a.num_positives == b.num_positives
        && a.comp_adjustment_method == b.comp_adjustment_method)
}

// No NCBI counterpart: the comparison of two results of the same match.
#[cfg(test)]
fn x_list_difference(a: &BlastpHspList, b: &BlastpHspList) -> Option<String> {
    if a.oid != b.oid
        || a.query_index != b.query_index
        || a.best_evalue.to_bits() != b.best_evalue.to_bits()
        || a.hsps.len() != b.hsps.len()
    {
        return Some(format!(
            "list oid {} / {}, query {} / {}, best E {:e} / {:e}, HSPs {} / {}",
            a.oid,
            b.oid,
            a.query_index,
            b.query_index,
            a.best_evalue,
            b.best_evalue,
            a.hsps.len(),
            b.hsps.len()
        ));
    }
    a.hsps
        .iter()
        .zip(&b.hsps)
        .position(|(x, y)| x_hsp_difference(x, y))
        .map(|index| {
            format!(
                "list oid {} HSP {index}: {:?} / {:?}",
                a.oid, a.hsps[index], b.hsps[index]
            )
        })
}

// No NCBI counterpart: the comparison of a computed result (first) with the one-thread
// computation's (second).
#[cfg(test)]
fn x_result_difference(
    computed: &Result<BlastpPostprocessResult>,
    reference: &std::result::Result<BlastpPostprocessResult, String>,
) -> Option<String> {
    match (computed, reference) {
        (Err(computed), Err(reference)) => {
            let computed = format!("{computed:#}");
            (computed != *reference).then(|| format!("errors {computed} / {reference}"))
        }
        (Ok(_), Err(reference)) => Some(format!("a result / the error {reference}")),
        (Err(computed), Ok(_)) => Some(format!("the error {computed:#} / a result")),
        (Ok(a), Ok(b)) => {
            if a.best_score != b.best_score || a.best_evalue.to_bits() != b.best_evalue.to_bits() {
                return Some(format!(
                    "best score {} / {}, best E {:e} / {:e}",
                    a.best_score, b.best_score, a.best_evalue, b.best_evalue
                ));
            }
            match (&a.hsp_list, &b.hsp_list) {
                (None, None) => None,
                (Some(a), Some(b)) => x_list_difference(a, b),
                (a, b) => Some(format!(
                    "HSP list {} / {}",
                    if a.is_some() { "present" } else { "absent" },
                    if b.is_some() { "present" } else { "absent" }
                )),
            }
        }
    }
}

#[cfg(test)]
#[path = "x_redo_par_tests.rs"]
mod tests;
