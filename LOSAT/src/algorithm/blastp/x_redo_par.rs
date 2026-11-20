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
//! only inserts; `pop` comes after it); and a heap filled to its cutoff (`n >= heapThreshold` and
//! `worstEvalue <= ecutoff`) stays so: an insert with an E-value within the cutoff keeps the worst
//! one kept within it, and an insert with a larger E-value does not replace the root, whose E-value
//! is within it (`BlastCompoHeap::insert`). So for a match's fixed `best_evalue` the test is
//! monotone over the replay: once true, true at every later point. Hence:
//! - a match whose test is already true on the heaps at a window's start is not computed (nor
//!   prefetched): its replay finds the test true and drops it (asserted);
//! - the test can be true for a match of query q only for `best_evalue > EVALUE_STRETCH x ecutoff`,
//!   and only if every heap of q's batch can be filled by then: the heaps of the queries before q
//!   and q's own hold at most their matches (`n` <= inserts <= matches), and those of the queries
//!   after q are untouched when q's matches are tested (stream order). With the heaps as the search
//!   makes them (`heapThreshold` = `max_target_seqs` >= 1) that is the last query of a batch whose
//!   queries all have at least `max_target_seqs` matches. At a window's start, also: a query of
//!   the batch whose matches are all replayed keeps its heap as it is, so if that heap is not
//!   filled no match of the batch can be ended any more; and q's heap gains at most one list per
//!   match of q, so q's next `heapThreshold - n` matches cannot be ended. A job with a match that
//!   these conditions leave open, not ended at its window's start, is at risk; the other jobs are
//!   never ended and are computed ahead freely;
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
//! With `-use_sw_tback` (item X4-mt) the Smith-Waterman redo runs on this path as well, with one
//! more step that depends on order. A BLASTP match has one window, and its redo reads the heaps
//! once: its query's `BlastCompo_HeapWouldInsert` at the window's first round, before any
//! alignment, with that round's E-value and score, which no other match changes. Every later step
//! of the match (its alignments, a failed or an endless redo) follows from its own inputs and that
//! answer; on FALSE the match ends there with no alignment and status 0. A job answers the test
//! from the heap as the window's start left it where the answer is certain (`x_sw_heap_test`):
//! TRUE when the E-value is within the heap's cutoff, or the heap holds fewer lists than its
//! threshold minus the query's matches before this one in the window (each inserts at most one);
//! FALSE when the test without its side effect says so (`BlastCompoHeap::would_insert_view`), as it
//! stays at every later point of the replay (the proof is at `x_sw_heap_test`). Otherwise the job
//! answers FALSE, so the match's redo stops at the test as NCBI's does on FALSE, and keeps the
//! E-value and score; the replay asks the query's heap with them at the point where the one-thread
//! loop asks it (before the match's insert, after every match before it). On FALSE the job's result
//! is NCBI's (no alignment, status 0) and nothing was computed in vain. On TRUE the replay stops at
//! the match: the match and the window's matches after it, with their results, go to the next
//! window, which finishes the match (a continuation: the same redo with the test answered TRUE, as
//! the replay answered it) and, in the same way, each carried match whose deferred test is TRUE for
//! certain on the heaps where the replay stopped (`x_sw_continued`: the carried matches before it
//! that keep no list do not count), while its other jobs compute the matches after the carried
//! ones. The carried results stay valid: a certain answer stays certain, and a deferred one is
//! asked at its replay. The replay also keeps the batch's status of NCBI's one-thread redo (`XSwBatch`): a
//! match that the redo reaches sets it, failed or not; a match the early termination test ends
//! leaves it; a batch whose status says failed at its end fails the search there, and an endless
//! redo fails it at that match. A continuation is replayed before every match after it, so a batch
//! ends after its continuations. A Smith-Waterman job is at risk only as an X-drop job is (the early
//! termination test may end its match). With `LOSAT_TIMING=1` the `[TIMING] x_bpredopar:` line of
//! such a search ends with `turned=<t> continued=<c>` (the deferred answers that the replay answered
//! FALSE, and TRUE).
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
//! ```c
//!         while (BlastHSPStreamRead(hsp_stream, &localMatch)
//!                 != kBlastHSPStream_Eof) {
//! ...
//!             theseMatches[i] = head->match;
//! ...
//!         int b;
//! #pragma omp for schedule(static)
//!         for (b = 0; b < numMatches; ++b) {
//! ```

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
use super::kappa::smith_waterman::{SwMatchScratch, SwRedoEndless, SwRedoFailed};
use super::kappa::{
    build_query_workspace, postprocess_preliminary_hits,
    postprocess_preliminary_hits_smith_waterman, x_newton_lanes_prefetch_matches,
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
    /// The errors of the Smith-Waterman redo's query batches (`-use_sw_tback`).
    pub sw_errors: XSwErrors,
}

/// The errors of a query batch of the Smith-Waterman redo, those of the one-thread loop
/// (`redo_run::sw_batch_failed`, `redo_run::sw_endless`): the batch that starts at query
/// `batch_start` fails at its end because its status says failed (given its first query's heap),
/// or at a match whose redo does not end.
#[derive(Clone, Copy)]
pub(crate) struct XSwErrors {
    pub batch_failed: fn(usize, &BlastCompoHeap) -> anyhow::Error,
    pub endless: fn(usize) -> anyhow::Error,
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
    /// The Smith-Waterman redo's memory, made for the worker's first such match.
    sw: Option<SwMatchScratch>,
}

impl XScratch {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334
    // ```c
    //         if ((int) compo_adjust_mode > 1 && !positionBased) {
    //             NRrecord_tld[i] = Blast_CompositionWorkspaceNew();
    //             status_code = Blast_CompositionWorkspaceInit(
    //                     NRrecord_tld[i],
    //                     scoringParams->options->matrix
    //             );
    // ```
    fn new(range_cache: BlastpKappaSubjectRangeCache, matrix: ScoringMatrix) -> Self {
        let mut scratch = Self {
            composition_workspace: BlastCompositionWorkspace::new_blosum62(),
            gap_scratch: GapAlignScratch::new(),
            hits: Vec::new(),
            range_cache,
            sw: None,
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
    /// The Smith-Waterman redo: the heap tests that a job deferred and its replay answered FALSE
    /// (`XReplayed::turned`).
    turned: Option<usize>,
    /// The Smith-Waterman redo: the heap tests that a job deferred and its replay answered TRUE,
    /// each finished by a continuation in the next window (`XCarry`).
    continued: usize,
    /// The seconds of the continued matches' jobs, which their continuations repeat up to the
    /// heap test (`XCarry::first_round_s`).
    rerun_s: f64,
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
// ```c
//         while (BlastHSPStreamRead(hsp_stream, &localMatch)
//                 != kBlastHSPStream_Eof) {
// ...
//                 BlastHSPList* localMatch = theseMatches[b];
// ```
/// The matches that are redone: those whose list has HSPs, counted from a query's compact list
/// while it is not expanded (`x_query_matches`).
pub(crate) fn x_matches(
    lists: &[Option<BlastpHitList>],
    compact: &[Option<BlastpHitList<XPrelimHsp>>],
) -> usize {
    x_query_matches(lists, compact).into_iter().sum()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1561-1576
// ```c
// int
// BlastCompo_EarlyTermination(double evalue,
//                             BlastCompo_Heap significantMatches[],
//                             int numQueries)
// {
//     ...
// }
// ```
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
    /// With `SW`, after a replay that stopped at a deferred heap test answered TRUE (`XCarry`): the
    /// carried matches whose test is TRUE (`x_sw_continued`), each the match of a job of its own (a
    /// continuation), with the test's E-value and score; by match index.
    continued: Vec<(usize, (f64, i32))>,
    /// The matches the replay before left, `matches[..carried]`: those not continued are in no job,
    /// their results come with the carry.
    carried: usize,
}

impl XWindow {
    /// The jobs of the stream that the window holds (not its continuation).
    fn stream_jobs(&self) -> usize {
        self.jobs.len() - self.continued.len()
    }
}

/// What a Smith-Waterman window's replay leaves to the next window when it stops at a match whose
/// deferred heap test it answered TRUE: that match and the window's matches after it, with their
/// results (`None`: ended before their window); the first one's result is its job's, which stopped
/// at the test (the replay of a match depends on the inserts of the matches before it).
#[derive(Debug)]
struct XCarry {
    matches: Vec<XMatch>,
    results: XWindowResults,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:251-276,329-400
/// The carried matches whose heap test, deferred by their job, is TRUE at its replay, with the
/// test's E-value and score. The first one's test the replay answered TRUE. A query's heap changes
/// only by the inserts of its own matches, so the tests after it are known query by query:
/// - in a query after the stopped match's, the query's carried matches are replayed ahead as
///   `x_replay` will replay them, on a copy of its heap without the alignments
///   (`BlastCompoHeap::keys`): a result with alignments is inserted where the heap's test says so,
///   and a deferred test is asked of the copy (FALSE: the match inserts nothing). This holds up to
///   the query's first TRUE, its first match that the early termination test may end at its
///   replay (`end_above`), or a redo that would end the search there;
/// - after that, and in the stopped match's query, a test is TRUE for certain where the heap (the
///   copy, or the heap as the replay left it) holds fewer lists than its threshold minus the
///   matches since that may insert one: those with alignments, those continued, and those whose
///   deferred test the heap without its side effect does not answer FALSE (a FALSE stays FALSE,
///   `x_sw_heap_test`). A match ended before its window inserts none.
///
/// A match that the early termination test may end is not continued (its continuation could be
/// computed in vain); nor is any match after a redo that would end the search.
fn x_sw_continued(
    carry: &XCarry,
    lists: &[Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    end_above: &[Option<f64>],
) -> Vec<(usize, (f64, i32))> {
    let mut continued = Vec::new();
    let mut query = None;
    // The query's heap replayed ahead (copied at its first match replayed so), exact while
    // `exact`; after that the count of the matches that may insert a list.
    let mut ahead: Option<BlastCompoHeap> = None;
    let mut exact = false;
    let mut may_insert = 0usize;
    for (index, (&matched, computed)) in carry.matches.iter().zip(&carry.results).enumerate() {
        let q_idx = matched.q_idx;
        if query != Some(q_idx) {
            query = Some(q_idx);
            ahead = None;
            exact = index > 0;
            may_insert = 0;
        }
        let Some(computed) = computed else {
            continue;
        };
        if index == 0 {
            let test = computed
                .deferred
                .expect("the stopped match's test was deferred");
            continued.push((index, test));
            may_insert = 1;
            continue;
        }
        let local_match = x_local_match(lists, matched);
        let inserts = match &computed.result {
            Ok(postprocessed) => postprocessed.hsp_list.is_some(),
            Err(error) if error.is::<SwRedoFailed>() => false,
            // The replay ends the search here (or the early termination test ends the match).
            Err(_) => break,
        };
        let may_end = end_above[q_idx].is_some_and(|above| local_match.best_evalue > above);
        if exact && !may_end {
            let heap = ahead.get_or_insert_with(|| heaps[q_idx].keys());
            match (computed.deferred, &computed.result) {
                (Some((evalue, score)), _) => {
                    if heap.would_insert(evalue, score, local_match.oid) {
                        continued.push((index, (evalue, score)));
                        exact = false;
                        may_insert = 1;
                    }
                }
                (None, Ok(postprocessed)) if inserts => {
                    let (evalue, score) = (postprocessed.best_evalue, postprocessed.best_score);
                    if heap.would_insert(evalue, score, local_match.oid) {
                        let keys = BlastpHspList {
                            oid: local_match.oid,
                            query_index: local_match.query_index,
                            hsps: Vec::new(),
                            best_evalue: evalue,
                        };
                        heap.insert(keys, evalue, score, local_match.oid);
                    }
                }
                (None, _) => {}
            }
            continue;
        }
        exact = false;
        let heap = ahead.as_ref().unwrap_or(&heaps[q_idx]);
        match computed.deferred {
            Some(test) if !may_end && heap.n.saturating_add(may_insert) < heap.heap_threshold => {
                continued.push((index, test));
                may_insert += 1;
            }
            Some((evalue, score)) => {
                may_insert += usize::from(heap.would_insert_view(evalue, score, local_match.oid));
            }
            None => may_insert += usize::from(inserts),
        }
    }
    continued
}

/// A match computed by its job: the result of its redo and, for a Smith-Waterman match whose heap
/// test the job deferred (`x_sw_heap_test`: the redo stopped at the test as on FALSE), that test's
/// E-value and score, which the replay asks the query's heap with.
#[derive(Debug)]
struct XComputed {
    result: Result<BlastpPostprocessResult>,
    deferred: Option<(f64, i32)>,
    /// With `LOSAT_TIMING=1` and `deferred`: the seconds of the job's redo of the match, which
    /// stopped at the test (what a continuation of the match does again); else 0.
    first_round_s: f64,
}

impl From<Result<BlastpPostprocessResult>> for XComputed {
    fn from(result: Result<BlastpPostprocessResult>) -> Self {
        Self {
            result,
            deferred: None,
            first_round_s: 0.0,
        }
    }
}

/// The results of a window's matches in order; `None` for a match not computed.
type XWindowResults = Vec<Option<XComputed>>;

/// What the early termination test ended in a window's replay, and the Smith-Waterman heap tests
/// that turned.
#[derive(Debug, Default, PartialEq, Eq)]
struct XReplayed {
    /// Matches without a result (not computed).
    skipped: usize,
    /// Matches with a result (computed and dropped).
    wasted: usize,
    /// Smith-Waterman matches whose heap test the job deferred and the replay answered FALSE:
    /// replayed without alignment, with the job's result.
    turned: usize,
    /// The queries of the dropped results, each once.
    dropped: Vec<usize>,
}

impl XReplayed {
    fn drop_result_of(&mut self, q_idx: usize) {
        if self.dropped.last() != Some(&q_idx) {
            self.dropped.push(q_idx);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3398-3399,3421
// ```c
//         while (BlastHSPStreamRead(hsp_stream, &localMatch)
//                 != kBlastHSPStream_Eof) {
// ...
//             theseMatches[i] = head->match;
// ```
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
// ```c
//                 if (BlastCompo_EarlyTermination(
//                         localMatch->best_evalue,
//                         redoneMatches,
//                         numQueries
//                 )) {
// ```
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
///   `heapThreshold - n` matches of the window cannot be ended.
///
/// With `SW` (the Smith-Waterman redo) and `carry` (the replay before stopped at a deferred heap
/// test answered TRUE), the window starts with the carried matches: those whose test is TRUE for
/// certain (`x_sw_continued`) are each a job of its own (a continuation), the others are in no job
/// (their results come with the carry); its other jobs start at `start`.
/// The carried matches count as the window's matches of their query (each inserts at most one
/// list). A deferred heap test wastes nothing (`x_sw_heap_test`), so a Smith-Waterman job is at risk
/// only as an X-drop job is.
fn x_window<const SW: bool>(
    lists: &[Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    query_batches: &[Range<usize>],
    end_above: &[Option<f64>],
    caps: &[usize],
    carry: Option<&XCarry>,
    start: XMatch,
    jobs: usize,
) -> XWindow {
    let mut window = XWindow {
        matches: Vec::new(),
        ended: Vec::new(),
        jobs: Vec::new(),
        next: None,
        at_risk: Vec::new(),
        continued: Vec::new(),
        carried: 0,
    };
    // The queries before the window's first match are replayed.
    let mut first_query = start.q_idx;
    // The query of the last job: whether its batch can still end matches, how many of its
    // matches the window holds, and how many of its jobs at risk.
    let mut query: Option<(usize, bool, usize, usize)> = None;
    let batch_open = |q_idx: usize, first_query: usize| {
        let batch = query_batches[q_idx].clone();
        end_above[q_idx].is_some()
            && (batch.start..batch.end.min(first_query)).all(|done| heaps[done].filled_to_cutoff())
    };
    if let (true, Some(carry)) = (SW, carry) {
        first_query = carry.matches[0].q_idx;
        window.matches.extend_from_slice(&carry.matches);
        window
            .ended
            .extend(carry.results.iter().map(|computed| computed.is_none()));
        window.carried = carry.matches.len();
        window.continued = x_sw_continued(carry, lists, heaps, end_above);
        window
            .jobs
            .extend(window.continued.iter().map(|&(index, _)| index..index + 1));
        let last = window.matches[window.matches.len() - 1].q_idx;
        let in_window = window.matches.len()
            - window
                .matches
                .partition_point(|matched| matched.q_idx < last);
        query = Some((last, batch_open(last, first_query), in_window, 0));
    }
    let mut cursor = start;
    let mut job = Vec::with_capacity(X_JOB_MATCHES);
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
            _ => (batch_open(q_idx, first_query), 0, 0),
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
        self.ahead -= window.stream_jobs();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2308-2336
// ```c
// static BlastCompo_QueryInfo *
// s_GetQueryInfo(Uint1 * query_data, const BlastQueryInfo * blast_query_info, Boolean skip)
// {
//     ...
//     query_info->eff_search_space = (double) query_context->eff_searchsp;
//     query_info->origin = query_context->query_offset;
//     query_info->seq.data = &query_data[query_info->origin];
//     query_info->seq.length = query_context->query_length;
// ```
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
// ```c
//                 query_index = localMatch->query_index;
// ...
//                         *pStatusCode =
//                                 Blast_RedoOneMatch(
//                                         alignments,             // thread-local
//                                         redo_align_params,      // thread-local
//                                         incoming_aligns,        // thread-local
// ```
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1463-1473
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:251-276
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:329-400
/// The answer that a job gives to the Smith-Waterman heap test of a match (NCBI asks
/// `BlastCompo_HeapWouldInsert` of the query's heap at the first round of the match's window, with
/// that round's E-value and score, while the query has no alignment): `heap` is the query's heap as
/// the window's start left it, `before` the query's matches in the window before this one. TRUE
/// when the answer at the replay is TRUE for certain, FALSE when it is FALSE for certain, else
/// FALSE with `deferred` set to the E-value and score that the replay asks the heap with: the redo
/// stops at the test as NCBI's does on FALSE, which is the match's result if the replay answers
/// FALSE; on TRUE a continuation in the next window finishes it (`XCarry`).
///
/// Certain TRUE: the E-value is within the heap's cutoff (the test's second condition), or `n +
/// before` is below the threshold (the first condition at the replay: only this query's matches
/// insert into its heap, each at most one list, and an insert at or past the threshold leaves `n`).
///
/// Certain FALSE: the test without its side effect says FALSE on the window's start, so `n >=
/// heapThreshold`, `e > ecutoff`, `e >= worstEvalue` and the candidate is not better than the root
/// (`compo_heap_record_compare`, E-value, then score, then subject index). A later insert keeps `n`
/// at or past the threshold. It either adds a list while every kept E-value and the new one are
/// within the cutoff (so the root's E-value, the new worst, is within it and below `e`), or replaces
/// the root by a better record and restores the heap (the new root is not worse than the old one, so
/// the worst E-value does not grow and the candidate is not better than the new root), or changes
/// nothing. Each case keeps all four conditions, so the test stays FALSE up to the replay
/// (`hsp.rs`, `test_compo_heap_would_insert_view_equals_would_insert_and_a_false_stays_false`).
fn x_sw_heap_test(
    heap: &BlastCompoHeap,
    before: usize,
    (evalue, score, subject_index): (f64, i32, u32),
    deferred: &mut Option<(f64, i32)>,
) -> bool {
    if evalue <= heap.ecutoff || heap.n.saturating_add(before) < heap.heap_threshold {
        return true;
    }
    if !heap.would_insert_view(evalue, score, subject_index) {
        return false;
    }
    // BLASTP's match has one window, whose test comes before its first alignment: one test.
    debug_assert!(deferred.is_none(), "a second heap test of one match");
    *deferred = Some((evalue, score));
    false
}

/// How a Smith-Waterman job answers its match's heap test.
#[derive(Clone, Copy)]
enum XSwTest<'h> {
    /// From the query's heap as the window's start left it, with the query's matches in the window
    /// before this one (`x_sw_heap_test`).
    Heap(&'h BlastCompoHeap, usize),
    /// A continuation (`XCarry`): TRUE, as the replay answered the test with this E-value and
    /// score.
    Continued((f64, i32)),
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3602-3621
/// The Smith-Waterman redo of one match (`-use_sw_tback`): the call of the one-thread loop, with its
/// heap test answered as `test` says.
fn x_redo_match_sw<P>(
    search: &XRedoParSearch<'_, P>,
    q_idx: usize,
    local_match: &BlastpHspList,
    query_workspace: &BlastpKappaQueryWorkspace,
    redo_align_params: &BlastRedoAlignParams,
    scratch: &mut XScratch,
    test: XSwTest<'_>,
) -> XComputed {
    scratch.hits.clear();
    scratch.hits.extend(
        local_match
            .hsps
            .iter()
            .map(preliminary_hit_from_local_match_hsp),
    );
    let subject = match search.subjects.read(local_match.oid) {
        Ok(subject) => subject,
        Err(error) => {
            return XComputed {
                result: Err(error),
                deferred: None,
                first_round_s: 0.0,
            }
        }
    };
    let subject_raw = subject.residues();
    let mut deferred = None;
    let started = search.timing.then(Instant::now);
    let result = postprocess_preliminary_hits_smith_waterman(
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
        scratch.sw.get_or_insert_with(SwMatchScratch::new),
        &mut |evalue, score| match test {
            XSwTest::Heap(heap, before) => x_sw_heap_test(
                heap,
                before,
                (evalue, score, local_match.oid),
                &mut deferred,
            ),
            XSwTest::Continued(asked) => {
                // The same inputs give the replay's test the same E-value and score.
                debug_assert!(
                    asked.0.to_bits() == evalue.to_bits() && asked.1 == score,
                    "a continuation's heap test differs from the deferred one"
                );
                true
            }
        },
    );
    let first_round_s = match (started, deferred) {
        (Some(started), Some(_)) => started.elapsed().as_secs_f64(),
        _ => 0.0,
    };
    XComputed {
        result,
        deferred,
        first_round_s,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3447-3449,3493-3511
// ```c
//         int b;
// #pragma omp for schedule(static)
//         for (b = 0; b < numMatches; ++b) {
// ...
//                 seqSrc               = seqsrc_tld[tid];
//                 scoringParams        = score_params_tld[tid];
//                 hitParams            = hit_params_tld[tid];
//                 redoneMatches        = redoneMatches_tld[tid];
//                 alignments           = alignments_tld[tid];
//                 incoming_align_set   = incoming_align_set_tld[tid];
//                 NRrecord             = NRrecord_tld[tid];
// ...
//                 BlastHSPList* localMatch = theseMatches[b];
// ```
/// The redo of a window's matches on the pool, in jobs of up to `X_JOB_MATCHES` consecutive
/// matches of one query; the results of the jobs' matches in order (all of the window's but the
/// carried ones), `None` for a match ended before the window (not computed). `heaps` are the heaps as the window's start left them; nothing changes them
/// while the jobs run. NCBI's loop also runs the early termination test and the heap insertion of
/// each match on its thread; here they are replayed afterwards (`x_replay`). With `lpt`
/// (`LOSAT_X_BPREDOLPT`), more than one thread on the search's pool and two or more jobs, the
/// jobs run as `x_compute_window_lpt` runs them; the second value says whether they did. With `SW`
/// each match runs the Smith-Waterman redo (`x_redo_match_sw`), its heap test answered from
/// `heaps`, or TRUE for a continuation (a match in `window.continued`).
fn x_compute_window<P, const SW: bool>(
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
    // ```c
    //         query_info_tld[i] = s_GetQueryInfo(
    //                 queryBlk->sequence,
    //                 queryInfo,
    //                 (program_number == eBlastTypeBlastx)
    //         );
    // ```
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
                    .map(|&ended| {
                        (!ended).then(|| XComputed {
                            result: Err(anyhow::anyhow!("{error:#}")),
                            deferred: None,
                            first_round_s: 0.0,
                        })
                    })
                    .collect();
            }
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529
        // ```c
        //                 if (BlastCompo_EarlyTermination(
        //                         localMatch->best_evalue,
        //                         redoneMatches,
        //                         numQueries
        //                 )) {
        // ```
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
        // With `SW`: the query's matches in the window before the job's first (they are
        // consecutive in stream order).
        let before_job = if SW {
            job.start
                - window
                    .matches
                    .partition_point(|matched| matched.q_idx < q_idx)
        } else {
            0
        };
        job_matches
            .iter()
            .zip(job_ended)
            .enumerate()
            .map(|(index, (&matched, &ended))| {
                (!ended).then(|| {
                    if SW {
                        let test = match window
                            .continued
                            .binary_search_by_key(&(job.start + index), |&(at, _)| at)
                        {
                            Ok(found) => XSwTest::Continued(window.continued[found].1),
                            Err(_) => XSwTest::Heap(&heaps[q_idx], before_job + index),
                        };
                        x_redo_match_sw(
                            search,
                            q_idx,
                            x_local_match(lists, matched),
                            query_workspace,
                            &redo_align_params,
                            scratch,
                            test,
                        )
                    } else {
                        x_redo_match(
                            search,
                            q_idx,
                            x_local_match(lists, matched),
                            query_workspace,
                            &redo_align_params,
                            scratch,
                        )
                        .into()
                    }
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
// ```c
//         int b;
// #pragma omp for schedule(static)
//         for (b = 0; b < numMatches; ++b) {
// ...
//                 seqSrc               = seqsrc_tld[tid];
//                 scoringParams        = score_params_tld[tid];
//                 hitParams            = hit_params_tld[tid];
// ```
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
// ```c
//                 if (BlastCompo_EarlyTermination(
//                         localMatch->best_evalue,
//                         redoneMatches,
//                         numQueries
//                 )) {
// ...
//                     if (BlastCompo_HeapWouldInsert(
//                             &redoneMatches[query_index],
//                             best_evalue,
//                             best_score,
//                             localMatch->oid
//                     )) {
//                         *pStatusCode =
//                                 BlastCompo_HeapInsert(
//                                         &redoneMatches[query_index],
//                                         hsp_list,
//                                         best_evalue,
//                                         best_score,
//                                         localMatch->oid,
//                                         &discarded_aligns
//                                 );
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3644-3646,3735-3737,3774-3778
// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:250-252
/// The status of a query batch in NCBI's one-thread Smith-Waterman redo (`status_code_tld[0]`),
/// replayed in stream order (`redo_run::one_thread_redo`'s `sw_last_failed`): a match that the redo
/// reaches sets it (failed, or 0 for a redo that ends without failing, with alignments or
/// without, as when its heap test is FALSE); a match the early termination test ends leaves it; at
/// the batch's end a failed status fails the search there (`XSwErrors::batch_failed`). It is kept
/// across windows: a window may end inside a batch.
struct XSwBatch {
    errors: XSwErrors,
    /// The first query of the batch of the last match replayed, if any.
    batch_start: Option<usize>,
    /// The status of that batch says failed.
    failed: bool,
}

impl XSwBatch {
    fn new(errors: XSwErrors) -> Self {
        Self {
            errors,
            batch_start: None,
            failed: false,
        }
    }

    /// Before the replay of a match of the batch that starts at `batch_start`: ends the batch of
    /// the match before if this one is in a later batch (the one-thread loop checks the status at
    /// each batch's last query; the batches between them have no match the redo reaches, so their
    /// status stays 0).
    fn enter(&mut self, batch_start: usize, heaps: &[BlastCompoHeap]) -> Result<()> {
        if self.batch_start != Some(batch_start) {
            self.end(heaps)?;
            self.batch_start = Some(batch_start);
        }
        Ok(())
    }

    /// A match that the redo reaches sets the batch's status: `failed` or 0.
    fn reached(&mut self, failed: bool) {
        self.failed = failed;
    }

    /// The end of the batch of the last match replayed: the error of a failed status.
    fn end(&mut self, heaps: &[BlastCompoHeap]) -> Result<()> {
        match self.batch_start {
            Some(start) if self.failed => Err((self.errors.batch_failed)(start, &heaps[start])),
            _ => {
                self.failed = false;
                Ok(())
            }
        }
    }
}

/// The steps of the one-thread loop that depend on order, for a window's matches in stream order,
/// with the results computed for them: the early termination test on the heaps of the query's
/// batch, the match's error, and the heap insertion. With `SW` (the Smith-Waterman redo) also each
/// deferred heap test on the query's heap (a FALSE replays the match with the job's result, a redo
/// that stopped at the test; a TRUE stops the replay there and leaves the match and the matches
/// after it to the next window in `carry`), and the batch's status (`batch`): a failed redo sets
/// it, an endless one fails the search at that match.
fn x_replay<const SW: bool>(
    matches: &[XMatch],
    results: XWindowResults,
    lists: &[Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    query_batches: &[Range<usize>],
    batch: &mut XSwBatch,
    carry: &mut Option<XCarry>,
) -> Result<XReplayed> {
    assert_eq!(matches.len(), results.len(), "one result per match");
    let mut replayed = XReplayed::default();
    let mut results = results.into_iter();
    for (index, &matched) in matches.iter().enumerate() {
        let computed = results.next().expect("one result per match");
        let local_match = x_local_match(lists, matched);
        let batch_start = query_batches[matched.q_idx].start;
        if SW {
            batch.enter(batch_start, heaps)?;
        }
        if blast_compo_early_termination(
            local_match.best_evalue,
            &heaps[query_batches[matched.q_idx].clone()],
        ) {
            if computed.is_some() {
                replayed.wasted += 1;
                replayed.drop_result_of(matched.q_idx);
            } else {
                replayed.skipped += 1;
            }
            continue;
        }
        // A match without a result had the test true at its window's start; the test is
        // monotone over the replay (`x_end_above`, header), so it is true here.
        let Some(XComputed {
            result,
            deferred,
            first_round_s,
        }) = computed
        else {
            unreachable!(
                "query {} match {}: ended before its window, not ended at its replay",
                matched.q_idx, matched.position
            );
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1463-1476
        // The heap test that the job left to the replay, on the heap as the matches before this
        // one left it. On FALSE NCBI's redo of the match ends at its first round with no
        // alignment and status 0: the job's result, a redo that stopped there. On TRUE the
        // match's redo goes on past the test: a continuation in the next window, and the replay
        // of the matches after it waits for it.
        if let (true, Some(test)) = (SW, deferred) {
            if heaps[matched.q_idx].would_insert(test.0, test.1, local_match.oid) {
                let stopped = XComputed {
                    result,
                    deferred,
                    first_round_s,
                };
                *carry = Some(XCarry {
                    matches: matches[index..].to_vec(),
                    results: std::iter::once(Some(stopped)).chain(results).collect(),
                });
                return Ok(replayed);
            }
            replayed.turned += 1;
        }
        let postprocessed = match result {
            Ok(postprocessed) => {
                if SW {
                    batch.reached(false);
                }
                postprocessed
            }
            // A failed redo leaves the match without alignments and sets the batch's status.
            Err(error) if SW && error.is::<SwRedoFailed>() => {
                batch.reached(true);
                continue;
            }
            Err(error) if SW && error.is::<SwRedoEndless>() => {
                return Err((batch.errors.endless)(batch_start));
            }
            Err(error) => return Err(error),
        };
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
// ```c
//         while (BlastHSPStreamRead(hsp_stream, &localMatch)
//                 != kBlastHSPStream_Eof) {
// ...
// #pragma omp for schedule(static)
//         for (b = 0; b < numMatches; ++b) {
// ...
// match_loop_cleanup:
//                 if (seqSrc) {
//                     localMatch = Blast_HSPListFree(localMatch);
//                 } else {
// ```
/// The windows of a search, each computed (`compute`, with the lists and the heaps as the replay
/// of the windows before left them) and then replayed in stream order into `heaps`. A query's
/// preliminary list is released after the replay of the window that holds its last match. With
/// `expand`, the compact lists are expanded ahead of each window. With `SW` (the Smith-Waterman
/// redo) the replay keeps the query batches' status across the windows (`XSwBatch`, its errors
/// `sw_errors`), and the batch of the last match ends after the last window.
#[allow(clippy::too_many_arguments)]
fn x_redo_windows<const SW: bool>(
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    query_batches: &[Range<usize>],
    limits: XWindowLimits,
    timing: bool,
    mut expand: Option<&mut XExpandAhead<'_>>,
    sw_errors: XSwErrors,
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
    if SW {
        counts.turned = Some(0);
    }
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
    let mut sw_batch = XSwBatch::new(sw_errors);
    // With `SW`: what the last replay left to the next window.
    let mut carry: Option<XCarry> = None;
    // With a carry and no match of the stream left, the window holds the carry only.
    let stream_end = XMatch {
        q_idx: lists.len(),
        position: 0,
    };
    while start.is_some() || carry.is_some() {
        let cursor = start.unwrap_or(stream_end);
        if let Some(expand) = expand.as_deref_mut() {
            expand.expand(lists, limits.jobs);
        }
        let window = x_window::<SW>(
            lists,
            heaps,
            query_batches,
            &end_above,
            &caps,
            carry.as_ref(),
            cursor,
            limits.jobs,
        );
        if let Some(expand) = expand.as_deref_mut() {
            expand.check(lists, cursor, &window);
        }
        if !window.jobs.is_empty() {
            counts.windows += 1;
            counts.matches += window.matches.len() - window.carried;
            counts.jobs += window.stream_jobs();
            if window
                .next
                .zip(window.matches.last())
                .is_some_and(|(next, last)| next.q_idx == last.q_idx)
            {
                counts.split += 1;
            }
            let par_start = timing.then(Instant::now);
            let mut results = compute(lists, heaps, &window);
            if let Some(par_start) = par_start {
                counts.par_s += par_start.elapsed().as_secs_f64();
            }
            // The carried results, the continued ones replaced by their continuations', then the
            // results of the jobs after them.
            if let Some(carry) = carry.take() {
                let mut computed = results.into_iter();
                results = carry.results;
                for index in window.jobs.iter().flat_map(|job| job.clone()) {
                    let result = computed.next().expect("one result per match of a job");
                    if index < window.carried {
                        let stopped = std::mem::replace(&mut results[index], result);
                        counts.continued += 1;
                        counts.rerun_s += stopped.map_or(0.0, |stopped| stopped.first_round_s);
                    } else {
                        results.push(result);
                    }
                }
            }
            let replay_start = timing.then(Instant::now);
            let replayed = match x_replay::<SW>(
                &window.matches,
                results,
                lists,
                heaps,
                query_batches,
                &mut sw_batch,
                &mut carry,
            ) {
                Ok(replayed) => replayed,
                Err(error) => return Err(x_timing_at_error(timing, &counts, error)),
            };
            counts.skipped += replayed.skipped;
            counts.wasted += replayed.wasted;
            if SW {
                *counts.turned.get_or_insert(0) += replayed.turned;
            }
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
        // A carried match keeps its query's list until its replay.
        let done = match &carry {
            Some(carry) => carry.matches[0].q_idx,
            None => window.next.map_or(lists.len(), |next| next.q_idx),
        };
        for list in &mut lists[released_to..done] {
            if list.take().is_some() {
                counts.released += 1;
            }
        }
        released_to = done;
        start = window.next;
    }
    // The end of the last batch that the replay reached (the one-thread loop checks its status
    // at the batch's last query; nothing after its last match changes it).
    if SW {
        if let Err(error) = sw_batch.end(heaps) {
            return Err(x_timing_at_error(timing, &counts, error));
        }
    }
    Ok(counts)
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
/// A redo that ends in `error` prints its counters up to there (with `timing`): the search
/// prints them only after a redo that ends.
fn x_timing_at_error(timing: bool, counts: &XCounts, error: anyhow::Error) -> anyhow::Error {
    if timing {
        x_print_timing(counts);
    }
    error
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
    // reported here, before any match, as there. They say which redo the search runs, the same for
    // every query (`-use_sw_tback`): chosen once here, so the X-drop redo's windows are compiled
    // without the Smith-Waterman steps.
    if (search.build_redo_align_params)(0)?.smith_waterman {
        x_redo_with::<P, true>(search, lists, heaps, expand)
    } else {
        x_redo_with::<P, false>(search, lists, heaps, expand)
    }
}

/// `x_redo` for the redo `SW` (the Smith-Waterman redo or the X-drop redo).
fn x_redo_with<P, const SW: bool>(
    search: &XRedoParSearch<'_, P>,
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    expand: Option<&mut XExpandAhead<'_>>,
) -> Result<XCounts>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
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
    let counts = x_redo_windows::<SW>(
        lists,
        heaps,
        search.query_batches,
        limits,
        search.timing,
        expand,
        search.sw_errors,
        &mut |lists, heaps, window| {
            let (results, lpt_ran) =
                x_compute_window::<P, SW>(search, lists, heaps, window, &shared, lpt);
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
        "[TIMING] x_bpredopar: matches={} jobs={} windows={} split={} ended={} skipped={} wasted={} par_s={:.6} replay_s={:.6}{}{}",
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
            .map_or(String::new(), |windows| format!(" lpt_windows={windows}")),
        counts
            .turned
            .map_or(String::new(), |turned| format!(
                " turned={turned} continued={} rerun_s={:.6}",
                counts.continued, counts.rerun_s
            ))
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
mod tests {
    use super::*;
    use crate::algorithm::blastp::encoding::EncodedProtein;
    use crate::utils::glibc_fma_libm::GlibcLibm;

    const ECUTOFF: f64 = 0.002;

    fn hsp(oid: u32, e_value: f64) -> BlastpHsp {
        BlastpHsp {
            num: 1,
            identity: 50.0,
            length: 10,
            mismatch: 5,
            gapopen: 0,
            q_start: 1,
            q_end: 10,
            s_start: 1,
            s_end: 10,
            gapped_q_start: 0,
            gapped_s_start: 0,
            e_value,
            bit_score: 20.0,
            num_ident: 5,
            query_context: 0,
            query_frame: 0,
            subject_frame: 0,
            query_length: 10,
            q_idx: 0,
            s_idx: oid,
            raw_score: 40,
            gap_info: None,
            num_positives: 6,
            comp_adjustment_method: 2,
        }
    }

    /// A list of subject `oid` with best E-value `best_evalue`; without HSPs when `empty`.
    fn list(oid: u32, best_evalue: f64, empty: bool) -> BlastpHspList {
        BlastpHspList {
            oid,
            query_index: 0,
            hsps: if empty {
                Vec::new()
            } else {
                vec![hsp(oid, best_evalue)]
            },
            best_evalue,
        }
    }

    fn hit_list(lists: Vec<BlastpHspList>) -> Option<BlastpHitList> {
        Some(BlastpHitList {
            hsplist_count: lists.len(),
            hsplist_max: 1050,
            worst_evalue: 0.0,
            low_score: 0,
            heapified: false,
            hsplist_array: lists,
            hsplist_current: 0,
            num_hits: 0,
        })
    }

    fn redone(oid: u32, best_evalue: f64) -> Result<BlastpPostprocessResult> {
        Ok(BlastpPostprocessResult {
            hsp_list: Some(list(oid, best_evalue, false)),
            best_score: 40,
            best_evalue,
        })
    }

    /// A computed match without a deferred heap test.
    fn done(result: Result<BlastpPostprocessResult>) -> Option<XComputed> {
        Some(result.into())
    }

    /// A computed match whose redo failed with `message`.
    fn failed(message: &str) -> Option<XComputed> {
        done(Err(anyhow::anyhow!("{message}")))
    }

    // No NCBI counterpart: the tests' stand-ins for `redo_run::sw_batch_failed` and
    // `redo_run::sw_endless`, which name the batch and what its first query's heap holds.
    fn test_batch_failed(batch_start: usize, first_query: &BlastCompoHeap) -> anyhow::Error {
        anyhow::anyhow!(
            "batch {batch_start} failed, first query holds {}",
            first_query.n
        )
    }

    fn test_endless(batch_start: usize) -> anyhow::Error {
        anyhow::anyhow!("endless in batch {batch_start}")
    }

    const TEST_SW_ERRORS: XSwErrors = XSwErrors {
        batch_failed: test_batch_failed,
        endless: test_endless,
    };

    /// The replay of the X-drop redo's window.
    fn replay(
        matches: &[XMatch],
        results: XWindowResults,
        lists: &[Option<BlastpHitList>],
        heaps: &mut [BlastCompoHeap],
        batches: &[Range<usize>],
    ) -> Result<XReplayed> {
        x_replay::<false>(
            matches,
            results,
            lists,
            heaps,
            batches,
            &mut XSwBatch::new(TEST_SW_ERRORS),
            &mut None,
        )
    }

    fn heaps_of(count: usize, heap_threshold: usize) -> Vec<BlastCompoHeap> {
        (0..count)
            .map(|_| BlastCompoHeap::new(heap_threshold, ECUTOFF))
            .collect()
    }

    /// Each query alone in its batch.
    fn alone(count: usize) -> Vec<Range<usize>> {
        (0..count).map(|q_idx| q_idx..q_idx + 1).collect()
    }

    fn at(q_idx: usize, position: usize) -> XMatch {
        XMatch { q_idx, position }
    }

    /// The window from `start` of at most `jobs` jobs, without a cap on jobs at risk.
    fn window_from(
        lists: &[Option<BlastpHitList>],
        heaps: &[BlastCompoHeap],
        batches: &[Range<usize>],
        end_above: &[Option<f64>],
        start: XMatch,
        jobs: usize,
    ) -> XWindow {
        let caps = vec![usize::MAX; lists.len()];
        x_window::<false>(lists, heaps, batches, end_above, &caps, None, start, jobs)
    }

    impl XWindow {
        fn at_risk_jobs(&self) -> usize {
            self.at_risk.iter().map(|&(_, jobs)| jobs).sum()
        }
    }

    fn replayed(skipped: usize, wasted: usize, dropped: &[usize]) -> XReplayed {
        XReplayed {
            skipped,
            wasted,
            turned: 0,
            dropped: dropped.to_vec(),
        }
    }

    /// All matches of `lists` in stream order.
    fn stream(lists: &[Option<BlastpHitList>]) -> Vec<XMatch> {
        let window = window_from(
            lists,
            &heaps_of(lists.len(), 1),
            &alone(lists.len()),
            &vec![None; lists.len()],
            at(0, 0),
            usize::MAX,
        );
        assert_eq!(window.next, None);
        window.matches
    }

    #[test]
    fn jobs_hold_up_to_four_matches_of_one_query_and_windows_end_at_their_job_budget() {
        let lists = vec![
            hit_list(vec![
                list(1, 1e-9, false),
                list(2, 1e-8, true),
                list(3, 1e-7, false),
                list(7, 1e-7, false),
                list(8, 1e-7, false),
                list(9, 1e-6, false),
                list(10, 1e-6, false),
            ]),
            None,
            hit_list(vec![list(4, 1e-9, false)]),
            hit_list(vec![list(5, 1e-9, false), list(6, 1e-8, false)]),
            None,
        ];
        assert_eq!(x_matches(&lists, &[]), 9);
        let (heaps, batches, none) = (heaps_of(5, 1), alone(5), vec![None; 5]);
        let window = |start, jobs| window_from(&lists, &heaps, &batches, &none, start, jobs);
        // The list without HSPs (query 0, position 1) is no match.
        let all = window(at(0, 0), usize::MAX);
        assert_eq!(
            all.matches,
            [
                at(0, 0),
                at(0, 2),
                at(0, 3),
                at(0, 4),
                at(0, 5),
                at(0, 6),
                at(2, 0),
                at(3, 0),
                at(3, 1)
            ]
        );
        assert_eq!(all.jobs, [0..4, 4..6, 6..7, 7..9]);
        assert_eq!((all.next, all.ended), (None, vec![false; 9]));
        // One job: the window ends inside query 0, and the next one starts at its fifth match
        // with the same job as above.
        let first = window(at(0, 0), 1);
        assert_eq!((first.jobs, first.next), (vec![0..4], Some(at(0, 5))));
        let second = window(at(0, 5), 1);
        assert_eq!(second.matches, [at(0, 5), at(0, 6)]);
        assert_eq!(second.next, Some(at(2, 0)));
        // The queries without matches after the last match end the stream.
        let last = window(at(2, 0), 2);
        assert_eq!((last.jobs, last.next), (vec![0..1, 1..3], None));
        assert!(window(at(4, 0), 1).jobs.is_empty());
    }

    /// Two queries of one batch, one subject each to keep (`-max_target_seqs 1`).
    fn two_query_batch() -> (
        Vec<Option<BlastpHitList>>,
        Vec<BlastCompoHeap>,
        Vec<Range<usize>>,
    ) {
        let lists = vec![
            hit_list(vec![list(10, 1e-20, false)]),
            hit_list(vec![
                list(20, 1e-20, false),
                list(21, 1.0, false),
                list(22, 0.005, false),
            ]),
        ];
        (lists, heaps_of(2, 1), vec![0..2, 0..2])
    }

    #[test]
    fn the_replay_drops_what_the_early_termination_test_ends_and_reports_the_next_error() {
        let (lists, mut heaps, batches) = two_query_batch();
        let matches = stream(&lists);
        assert_eq!(matches.len(), 4);
        // The match of best E 1.0 comes after both heaps are filled to the cutoff: the test
        // ends it, so its error is dropped. The next one (0.005 <= 5 * 0.002) is redone, and its
        // error is the one reported.
        let results = vec![
            done(redone(10, 1e-20)),
            done(redone(20, 1e-20)),
            failed("ended"),
            failed("reported"),
        ];
        let error = replay(&matches, results, &lists, &mut heaps, &batches)
            .expect_err("the redone match's error");
        assert_eq!(error.to_string(), "reported");
        assert_eq!((heaps[0].n, heaps[1].n), (1, 1));
    }

    #[test]
    fn the_replay_inserts_in_stream_order_and_counts_ended_results() {
        let (lists, heaps, batches) = two_query_batch();
        let matches = stream(&lists);
        // The third match is ended with a result (computed) or without one (not computed).
        for third in [failed("ended"), None] {
            let computed = third.is_some();
            let mut heaps = heaps.clone();
            let results = vec![
                done(redone(10, 1e-20)),
                done(redone(20, 1e-20)),
                third,
                done(redone(22, 0.005)),
            ];
            let ended = replay(&matches, results, &lists, &mut heaps, &batches)
                .expect("no error is redone");
            assert_eq!(
                ended,
                if computed {
                    replayed(0, 1, &[1])
                } else {
                    replayed(1, 0, &[])
                }
            );
            // Query 1 keeps its best subject: 22 (E 0.005, above the cutoff) does not replace 20.
            let hit_lists = fill_results_from_compo_heaps(1, &mut heaps);
            let oids: Vec<Vec<u32>> = hit_lists
                .iter()
                .map(|hit_list| {
                    let hit_list = hit_list.as_ref().expect("a hit list per query");
                    hit_list.hsplist_array[..hit_list.hsplist_count]
                        .iter()
                        .map(|list| list.oid)
                        .collect()
                })
                .collect();
            assert_eq!(oids, [vec![10], vec![20]]);
        }
    }

    #[test]
    #[should_panic(expected = "ended before its window, not ended at its replay")]
    fn a_match_without_a_result_is_ended_at_its_replay() {
        let (lists, mut heaps, batches) = two_query_batch();
        let matches = stream(&lists);
        // The first match cannot be ended (no heap is filled): it needs its result.
        let results = vec![None, None, None, None];
        let _ = replay(&matches, results, &lists, &mut heaps, &batches);
    }

    #[test]
    fn a_query_alone_in_its_batch_reads_only_its_own_heap() {
        // Query 0 has no match, so in a batch with query 1 the test is never true for query 1;
        // alone in its batch, query 1's second match is ended once its own heap is filled.
        let lists = vec![
            None,
            hit_list(vec![list(20, 1e-20, false), list(21, 1.0, false)]),
        ];
        for (batches, ended) in [
            (vec![0..2, 0..2], replayed(0, 0, &[])),
            (vec![0..1, 1..2], replayed(0, 1, &[1])),
        ] {
            let mut heaps = heaps_of(2, 1);
            let matches = stream(&lists);
            let results = vec![done(redone(20, 1e-20)), done(redone(21, 1.0))];
            assert_eq!(
                replay(&matches, results, &lists, &mut heaps, &batches).expect("no error"),
                ended,
                "{batches:?}"
            );
        }
    }

    /// The pre-check of a window rests on this: for a fixed E-value, once the early termination
    /// test is true on a batch's heaps it stays true after any later `would_insert` / `insert`
    /// (the replay's steps) on them.
    #[test]
    fn the_early_termination_test_stays_true_after_later_inserts() {
        // A fixed xorshift sequence.
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = move |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound as u64) as usize
        };
        // E-values around the cutoff (0.002), its stretch (0.01) and the epsilon of the heap's
        // comparison (1e-180).
        let evalues = [
            0.0,
            1e-200,
            1e-181,
            1e-5,
            ECUTOFF,
            0.004,
            0.01,
            0.010_000_1,
            0.05,
            1.0,
            10.0,
        ];
        let mut true_then_tested = 0;
        for heap_threshold in [1, 2, 3, 5] {
            for _ in 0..200 {
                let mut heaps = heaps_of(3, heap_threshold);
                let mut ended = [false; 11];
                for _ in 0..40 {
                    let q_idx = next(3);
                    // Inserts within the cutoff twice as often, so that heaps fill.
                    let evalue = evalues[if next(3) == 0 {
                        next(evalues.len())
                    } else {
                        next(5)
                    }];
                    let (score, oid) = (next(4) as i32, next(6) as u32);
                    if heaps[q_idx].would_insert(evalue, score, oid) {
                        heaps[q_idx].insert(list(oid, evalue, false), evalue, score, oid);
                    }
                    for (index, &tested) in evalues.iter().enumerate() {
                        let now = blast_compo_early_termination(tested, &heaps);
                        assert!(
                            now || !ended[index],
                            "E {tested}: true before an insert of E {evalue}, false after"
                        );
                        true_then_tested += usize::from(ended[index]);
                        ended[index] = now;
                    }
                }
            }
        }
        assert!(true_then_tested > 1000, "{true_then_tested}");
    }

    #[test]
    fn matches_may_be_ended_for_a_query_whose_batch_heaps_can_all_be_filled() {
        // Batch 0..3: queries 0 and 2 have two matches, query 1 one; query 3 has a list
        // without HSPs only, query 4 none.
        let lists = vec![
            hit_list(vec![list(1, 1e-9, false), list(2, 1.0, false)]),
            hit_list(vec![list(1, 1e-9, false)]),
            hit_list(vec![
                list(1, 1e-9, false),
                list(2, 1.0, false),
                list(3, 1.0, true),
            ]),
            hit_list(vec![list(1, 1.0, true)]),
            None,
        ];
        let batches = vec![0..3, 0..3, 0..3, 3..4, 4..5];
        let stretch = Some(COMPO_HEAP_EVALUE_STRETCH * ECUTOFF);
        // One subject to keep: only the last query of the batch.
        assert_eq!(
            x_end_above(&x_query_matches(&lists, &[]), &heaps_of(5, 1), &batches),
            [None, None, stretch, None, None]
        );
        // Two: query 1 cannot fill its heap.
        assert_eq!(
            x_end_above(&x_query_matches(&lists, &[]), &heaps_of(5, 2), &batches),
            [None; 5]
        );
        // A heap after a query that is filled already counts as such.
        let mut heaps = heaps_of(5, 1);
        heaps[2].insert(list(1, 1e-9, false), 1e-9, 40, 1);
        assert_eq!(
            x_end_above(&x_query_matches(&lists, &[]), &heaps, &batches),
            [None, stretch, stretch, None, None]
        );
    }

    /// The driver with results made from the preliminary lists (`redone` of each list's best E)
    /// instead of the redo; returns the counts and how many matches were computed.
    fn drive(
        lists: &mut [Option<BlastpHitList>],
        heaps: &mut [BlastCompoHeap],
        batches: &[Range<usize>],
        limits: XWindowLimits,
        window_seen: impl FnMut(&[Option<BlastpHitList>], &XWindow),
    ) -> (XCounts, usize) {
        drive_expanding(lists, heaps, batches, limits, None, window_seen)
    }

    /// `drive` with the compact lists of `expand` expanded ahead of the windows.
    fn drive_expanding(
        lists: &mut [Option<BlastpHitList>],
        heaps: &mut [BlastCompoHeap],
        batches: &[Range<usize>],
        limits: XWindowLimits,
        expand: Option<&mut XExpandAhead<'_>>,
        mut window_seen: impl FnMut(&[Option<BlastpHitList>], &XWindow),
    ) -> (XCounts, usize) {
        let mut computed = 0;
        let counts = x_redo_windows::<false>(
            lists,
            heaps,
            batches,
            limits,
            false,
            expand,
            TEST_SW_ERRORS,
            &mut |lists, _, window| {
                window_seen(lists, window);
                window
                    .matches
                    .iter()
                    .zip(&window.ended)
                    .map(|(&matched, &ended)| {
                        (!ended).then(|| {
                            computed += 1;
                            let local_match = x_local_match(lists, matched);
                            redone(local_match.oid, local_match.best_evalue).into()
                        })
                    })
                    .collect()
            },
        )
        .expect("no error");
        (counts, computed)
    }

    #[test]
    fn a_query_alone_in_its_batch_computes_one_window_of_jobs_at_risk_once_its_heap_fills() {
        // `-max_target_seqs 1`: the first match fills the heap, and the test ends every later
        // one (E 1.0 > 5 x 0.002). The one-thread loop redoes only the first.
        const MATCHES: u32 = 2000;
        const THREADS: usize = 4;
        let mut lists = vec![hit_list(
            (0..MATCHES)
                .map(|oid| list(oid, if oid == 0 { 1e-20 } else { 1.0 }, false))
                .collect(),
        )];
        let mut heaps = heaps_of(1, 1);
        let limits = XWindowLimits {
            jobs: X_WINDOW_JOBS_PER_THREAD * THREADS,
            at_risk_jobs: THREADS,
        };
        let (counts, computed) = drive(&mut lists, &mut heaps, &alone(1), limits, |_, _| {});
        // The first window holds the jobs at risk that the cap allows (4 x 4 matches; only the
        // first match fits in the empty heap, and it shares its job with three at risk); the test
        // ends 15 of them at the replay, which returns the cap to four, and every later match
        // before its window, which leaves no job at risk.
        assert_eq!(computed, X_JOB_MATCHES * THREADS);
        assert_eq!(counts.wasted, X_JOB_MATCHES * THREADS - 1);
        assert_eq!(counts.skipped, MATCHES as usize - X_JOB_MATCHES * THREADS);
        assert_eq!(counts.skipped + counts.wasted, MATCHES as usize - 1);
        // The later windows hold their whole job budget (1,024 matches).
        assert_eq!((counts.windows, counts.split), (3, 2));
        assert_eq!((counts.matches, counts.released), (MATCHES as usize, 1));
        assert_eq!(heaps[0].n, 1);
    }

    /// The windows of a search where the heaps can be filled but are not: the cap on jobs at risk
    /// grows while the test ends nothing.
    #[test]
    fn the_cap_on_jobs_at_risk_grows_while_the_test_ends_nothing() {
        // `-max_target_seqs 5`, a query alone in its batch with 4,000 matches whose results stay
        // above the cutoff (E 1.0): the heap keeps five lists but is never filled to its cutoff,
        // so the test ends nothing, though every match after the fifth may be ended.
        const MATCHES: u32 = 4000;
        const THREADS: usize = 4;
        let mut lists = vec![hit_list(
            (0..MATCHES).map(|oid| list(oid, 1.0, false)).collect(),
        )];
        let mut heaps = heaps_of(1, 5);
        let limits = XWindowLimits {
            jobs: X_WINDOW_JOBS_PER_THREAD * THREADS,
            at_risk_jobs: THREADS,
        };
        let mut at_risk = Vec::new();
        let (counts, computed) = drive(&mut lists, &mut heaps, &alone(1), limits, |_, window| {
            at_risk.push(window.at_risk_jobs())
        });
        // The first window: the first job (its matches fit in the empty heap) and four at risk;
        // then the cap doubles with each window up to the job budget (256): 1,000 jobs in nine
        // windows, where a cap fixed at four takes 250.
        assert_eq!(at_risk, [4, 8, 16, 32, 64, 128, 256, 256, 235]);
        assert_eq!((counts.windows, counts.split), (9, 8));
        assert_eq!(
            (counts.skipped, counts.wasted, computed),
            (0, 0, MATCHES as usize)
        );
        assert_eq!(heaps[0].n, 5);
    }

    /// Batches whose heaps never fill before one that fills: each query's cap on jobs at risk
    /// grows with its own work only, so the batch that fills drops at most its own cap's results.
    #[test]
    fn a_batch_that_fills_drops_at_most_its_own_cap_after_batches_that_never_fill() {
        const THREADS: usize = 4;
        // `-max_target_seqs 1`, each query alone in its batch: queries 0 and 1 keep a list above
        // the cutoff (E 1.0), so their heaps are never filled; query 2's first match fills its
        // heap, and the test ends the others.
        let weak = |count: u32| hit_list((0..count).map(|oid| list(oid, 1.0, false)).collect());
        let mut lists = vec![
            weak(200),
            weak(200),
            hit_list(
                (0..400)
                    .map(|oid| list(oid, if oid == 0 { 1e-20 } else { 1.0 }, false))
                    .collect(),
            ),
        ];
        let mut heaps = heaps_of(3, 1);
        let limits = XWindowLimits {
            jobs: X_WINDOW_JOBS_PER_THREAD * THREADS,
            at_risk_jobs: THREADS,
        };
        let mut at_risk = Vec::new();
        let (counts, computed) = drive(&mut lists, &mut heaps, &alone(3), limits, |_, window| {
            at_risk.push(window.at_risk.clone())
        });
        // The caps of queries 0 and 1 double over their own windows (50 jobs each; the first
        // match of each fits in the empty heap), and a window holds jobs at risk of two queries
        // where the first one's end within its cap. Query 2 starts at four: its first window
        // computes 16 matches, the test drops 15 of them, and the later 384 are ended before
        // their window.
        assert_eq!(
            at_risk,
            [
                vec![(0, 4)],
                vec![(0, 8)],
                vec![(0, 16)],
                vec![(0, 22), (1, 4)],
                vec![(1, 8)],
                vec![(1, 16)],
                vec![(1, 22), (2, 4)],
                vec![],
            ]
        );
        assert_eq!((counts.wasted, counts.skipped, computed), (15, 384, 416));
        assert_eq!((counts.windows, counts.split), (8, 7));
    }

    #[test]
    fn a_replayed_query_with_an_unfilled_heap_and_room_in_the_query_heap_leave_jobs_not_at_risk() {
        // A replayed query of the batch: query 0's one match is replayed when the window starts at
        // query 1 (eight matches, `-max_target_seqs 1`).
        let lists = vec![
            hit_list(vec![list(1, 1.0, false)]),
            hit_list((0..8).map(|oid| list(oid, 1.0, false)).collect()),
        ];
        let batches = vec![0..2, 0..2];
        let end_above = x_end_above(&x_query_matches(&lists, &[]), &heaps_of(2, 1), &batches);
        assert_eq!(end_above, [None, Some(COMPO_HEAP_EVALUE_STRETCH * ECUTOFF)]);
        let at_risk = |heaps: &[BlastCompoHeap], start| {
            window_from(&lists, heaps, &batches, &end_above, start, usize::MAX).at_risk_jobs()
        };
        // Query 0 kept a list above the cutoff: its heap is not filled and stays so, so no
        // match of query 1 can be ended.
        let mut heaps = heaps_of(2, 1);
        heaps[0].insert(list(1, 1.0, false), 1.0, 40, 1);
        assert_eq!(at_risk(&heaps, at(1, 0)), 0);
        // Within the cutoff: both jobs of query 1 are at risk (all its matches but the first).
        let mut heaps = heaps_of(2, 1);
        heaps[0].insert(list(1, 1e-9, false), 1e-9, 40, 1);
        assert_eq!(at_risk(&heaps, at(1, 0)), 2);
        // A window that starts at query 0 cannot count on its heap yet.
        assert_eq!(at_risk(&heaps_of(2, 1), at(0, 0)), 2);

        // Room in the query's heap: a query alone with 12 matches, `-max_target_seqs 8`, and a
        // window from its fifth match. The heap gains at most one list per match, so with `kept`
        // lists the window's first 8 - `kept` matches cannot be ended.
        let lists = vec![hit_list((0..12).map(|oid| list(oid, 1.0, false)).collect())];
        let end_above = x_end_above(&x_query_matches(&lists, &[]), &heaps_of(1, 8), &alone(1));
        for (kept, jobs_at_risk) in [(0, 0), (4, 1)] {
            let mut heaps = heaps_of(1, 8);
            for oid in 0..kept {
                heaps[0].insert(list(oid, 1.0, false), 1.0, 40, oid);
            }
            let window = window_from(&lists, &heaps, &alone(1), &end_above, at(0, 4), usize::MAX);
            assert_eq!(window.at_risk_jobs(), jobs_at_risk, "{kept} kept");
        }
    }

    #[test]
    fn a_query_that_spans_windows_keeps_its_list_until_its_last_match_is_replayed() {
        let matches = |count: u32| hit_list((0..count).map(|oid| list(oid, 1e-9, false)).collect());
        let mut lists = vec![
            matches(2),
            None,
            matches(9),
            hit_list(vec![list(0, 1e-9, true)]),
            matches(1),
        ];
        let mut heaps = heaps_of(5, 1);
        let mut present = Vec::new();
        let (counts, computed) = drive(
            &mut lists,
            &mut heaps,
            &alone(5),
            XWindowLimits {
                jobs: 1,
                at_risk_jobs: 1,
            },
            |lists, window| {
                present.push((
                    window.matches[0],
                    lists.iter().map(Option::is_some).collect::<Vec<_>>(),
                ))
            },
        );
        // One job per window: query 2 spans three. A query's list is there until the window
        // that holds its last match is replayed.
        let (yes, no) = (true, false);
        assert_eq!(
            present,
            [
                (at(0, 0), vec![yes, no, yes, yes, yes]),
                (at(2, 0), vec![no, no, yes, yes, yes]),
                (at(2, 4), vec![no, no, yes, yes, yes]),
                (at(2, 8), vec![no, no, yes, yes, yes]),
                (at(4, 0), vec![no, no, no, no, yes]),
            ]
        );
        assert_eq!((counts.windows, counts.split, computed), (5, 2, 12));
        assert_eq!(counts.released, 4);
        assert!(lists.iter().all(Option::is_none));
        // Every match is within the cutoff: each heap keeps all of its query's.
        assert_eq!(
            heaps.iter().map(|heap| heap.n).collect::<Vec<_>>(),
            [2, 0, 9, 0, 1]
        );
    }

    /// The compact list of a preliminary hit list: the same fields,
    /// subjects, best E-values and HSPs per list.
    fn compact_of(list: &Option<BlastpHitList>) -> Option<BlastpHitList<XPrelimHsp>> {
        list.as_ref().map(|list| BlastpHitList {
            hsplist_count: list.hsplist_count,
            hsplist_max: list.hsplist_max,
            worst_evalue: list.worst_evalue,
            low_score: list.low_score,
            heapified: list.heapified,
            hsplist_array: list
                .hsplist_array
                .iter()
                .map(|hsp_list| BlastpHspList {
                    oid: hsp_list.oid,
                    query_index: hsp_list.query_index,
                    hsps: hsp_list
                        .hsps
                        .iter()
                        .map(|hsp| XPrelimHsp {
                            hsp: BlastpPreliminaryHsp {
                                query_start: hsp.q_start as i32 - 1,
                                query_end: hsp.q_end as i32,
                                subject_start: hsp.s_start as i32 - 1,
                                subject_end: hsp.s_end as i32,
                                gapped_query_start: 0,
                                gapped_subject_start: 0,
                                raw_score: hsp.raw_score,
                                query_context: 0,
                                query_frame: 0,
                                subject_frame: 0,
                                query_length: hsp.query_length,
                                q_idx: 0,
                                s_idx: hsp_list.oid,
                            },
                            e_value: hsp.e_value,
                        })
                        .collect(),
                    best_evalue: hsp_list.best_evalue,
                })
                .collect(),
            hsplist_current: list.hsplist_current,
            num_hits: list.num_hits,
        })
    }

    /// A window as the tests record it.
    type Seen = (
        Vec<XMatch>,
        Vec<bool>,
        Vec<Range<usize>>,
        Option<XMatch>,
        Vec<(usize, usize)>,
    );

    fn seen(window: &XWindow) -> Seen {
        (
            window.matches.clone(),
            window.ended.clone(),
            window.jobs.clone(),
            window.next,
            window.at_risk.clone(),
        )
    }

    /// The driver on `lists` expanded before the windows and on their compact lists expanded ahead
    /// of the windows: the same windows, counts, matches computed and heaps. Returns, per window of
    /// the second run, which lists were there, and its counts.
    fn same_windows_expanding_ahead(
        lists: &[Option<BlastpHitList>],
        heap_threshold: usize,
        batches: &[Range<usize>],
        limits: XWindowLimits,
    ) -> (Vec<Vec<bool>>, XCounts) {
        let count = lists.len();
        let mut full_lists = lists.to_vec();
        let mut full_heaps = heaps_of(count, heap_threshold);
        let mut full_windows = Vec::new();
        let (full_counts, full_computed) = drive(
            &mut full_lists,
            &mut full_heaps,
            batches,
            limits,
            |_, window| full_windows.push(seen(window)),
        );
        let mut compact: Vec<_> = lists.iter().map(compact_of).collect();
        let mut ahead_lists: Vec<Option<BlastpHitList>> = vec![None; count];
        let mut expand = XExpandAhead::new(&ahead_lists, &mut compact);
        let mut heaps = heaps_of(count, heap_threshold);
        let mut windows = Vec::new();
        let mut present = Vec::new();
        let (counts, computed) = drive_expanding(
            &mut ahead_lists,
            &mut heaps,
            batches,
            limits,
            Some(&mut expand),
            |lists, window| {
                windows.push(seen(window));
                present.push(lists.iter().map(Option::is_some).collect());
            },
        );
        assert_eq!(windows, full_windows);
        assert_eq!(
            (
                counts.windows,
                counts.split,
                counts.matches,
                counts.skipped,
                counts.wasted,
                counts.released,
                computed
            ),
            (
                full_counts.windows,
                full_counts.split,
                full_counts.matches,
                full_counts.skipped,
                full_counts.wasted,
                full_counts.released,
                full_computed
            )
        );
        let ns = |heaps: &[BlastCompoHeap]| heaps.iter().map(|heap| heap.n).collect::<Vec<_>>();
        assert_eq!(ns(&heaps), ns(&full_heaps));
        drop(expand);
        assert!(compact.iter().all(Option::is_none), "every list expanded");
        assert!(
            ahead_lists.iter().all(Option::is_none),
            "every list released"
        );
        (present, counts)
    }

    /// Before each window the compact lists are expanded in query order until the queries from
    /// the window's start hold more jobs than a window, and a query's list is released after the
    /// window that holds its last match, so the lists there are those of about one window.
    #[test]
    fn windows_expand_the_compact_lists_as_they_reach_them_and_choose_the_same_windows() {
        let matches = |count: u32| hit_list((0..count).map(|oid| list(oid, 1e-9, false)).collect());
        // Jobs per query: 1, 0, 3, 0 (a list without HSPs), 1, 2, 2.
        let lists = vec![
            matches(2),
            None,
            matches(9),
            hit_list(vec![list(0, 1e-9, true)]),
            matches(1),
            matches(6),
            matches(5),
        ];
        // Heaps that keep every list: no job is at risk, and each window holds two jobs.
        let limits = XWindowLimits {
            jobs: 2,
            at_risk_jobs: 1,
        };
        let (present, _) = same_windows_expanding_ahead(&lists, 100, &alone(7), limits);
        let (yes, no) = (true, false);
        assert_eq!(
            present,
            [
                // Queries 0-2 (four jobs) for the window of query 0 and query 2's first job.
                vec![yes, no, yes, no, no, no, no],
                // Query 2's last two jobs, with queries 3 and 4 (one job) ahead.
                vec![no, no, yes, yes, yes, no, no],
                vec![no, no, no, no, yes, yes, no],
                vec![no, no, no, no, no, yes, yes],
                vec![no, no, no, no, no, no, yes],
            ]
        );
    }

    /// The same windows where the early termination test ends matches and the caps on jobs at
    /// risk end windows: queries 0 and 1 never fill their heaps, query 2's first match fills its
    /// heap (`-max_target_seqs 1`), and with one batch of three the test also reads the heaps of
    /// the queries before.
    #[test]
    fn windows_expanding_ahead_end_and_drop_the_matches_of_the_full_expansion() {
        const THREADS: usize = 4;
        let weak = |count: u32| hit_list((0..count).map(|oid| list(oid, 1.0, false)).collect());
        let strong = hit_list(
            (0..400)
                .map(|oid| list(oid, if oid == 0 { 1e-20 } else { 1.0 }, false))
                .collect(),
        );
        let lists = vec![weak(200), weak(200), strong];
        // Windows of 256 jobs (all 200 expanded before the first) and of 8 jobs, where each query's
        // list is expanded when the windows reach it and the caps on jobs at risk end windows.
        for (batches, limits) in [
            (
                alone(3),
                XWindowLimits {
                    jobs: X_WINDOW_JOBS_PER_THREAD * THREADS,
                    at_risk_jobs: THREADS,
                },
            ),
            (
                alone(3),
                XWindowLimits {
                    jobs: 8,
                    at_risk_jobs: 1,
                },
            ),
            (
                vec![0..3; 3],
                XWindowLimits {
                    jobs: 8,
                    at_risk_jobs: 1,
                },
            ),
        ] {
            let (present, counts) = same_windows_expanding_ahead(&lists, 1, &batches, limits);
            if limits.jobs == 8 && batches[2] == (2..3) {
                // Query 0's list alone at the first window; query 2's matches ended before
                // their window and at their replay.
                assert_eq!(present[0], [true, false, false]);
                assert!(counts.skipped > 0 && counts.wasted > 0, "{counts:?}");
            }
        }
    }

    #[test]
    fn results_differ_by_any_field_or_by_error() {
        let a = redone(1, 1e-10);
        assert_eq!(
            x_result_difference(&a, &redone(1, 1e-10).map_err(|e| e.to_string())),
            None
        );
        let mut b = redone(1, 1e-10).expect("a result");
        b.hsp_list.as_mut().expect("a list").hsps[0].gapped_s_start = 1;
        assert!(x_result_difference(&a, &Ok(b)).is_some());
        assert!(x_result_difference(&a, &Err("failed".to_string())).is_some());
        assert_eq!(
            x_result_difference(&Err(anyhow::anyhow!("failed")), &Err("failed".to_string())),
            None
        );
    }

    /// A fixed xorshift sequence.
    struct XorShift(u64);

    impl XorShift {
        fn below(&mut self, bound: usize) -> usize {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 % bound as u64) as usize
        }
    }

    #[test]
    fn the_lpt_order_is_a_permutation_largest_first_with_ties_in_index_order() {
        assert!(x_lpt_order(&[]).is_empty());
        assert_eq!(x_lpt_order(&[3, 9, 3, 0, 9, 4]), [1, 4, 5, 0, 2, 3]);
        let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
        for len in [1, 2, 7, 64, 300] {
            // Few distinct costs (many ties), and many.
            for distinct in [1, 3, 1000] {
                let costs: Vec<u64> = (0..len).map(|_| rng.below(distinct) as u64).collect();
                let order = x_lpt_order(&costs);
                let mut jobs = order.clone();
                jobs.sort_unstable();
                assert_eq!(jobs, (0..len).collect::<Vec<_>>(), "a permutation");
                for pair in order.windows(2) {
                    let (a, b) = (pair[0], pair[1]);
                    assert!(
                        costs[a] > costs[b] || (costs[a] == costs[b] && a < b),
                        "job {a} (cost {}) before job {b} (cost {})",
                        costs[a],
                        costs[b]
                    );
                }
            }
        }
    }

    /// The scope runs every job once, with at most one scratch per task, and returns the results
    /// by job index; also when the jobs run nested parallel work, where a thread that waits in it
    /// may run another task's loop (no task waits for another one).
    #[test]
    fn the_lpt_scope_runs_every_job_once_and_returns_the_results_by_job_index() {
        for threads in [2, 3, 7] {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                for jobs in [1_usize, 2, 5, 40] {
                    let costs: Vec<u64> = (0..jobs as u64).map(|job| (job * 7919) % 13).collect();
                    let order = x_lpt_order(&costs);
                    let tasks = threads.min(jobs);
                    let (scratches, runs) = (AtomicUsize::new(0), AtomicUsize::new(0));
                    let results = x_compute_window_lpt(
                        pool,
                        tasks,
                        &order,
                        || {
                            scratches.fetch_add(1, Ordering::Relaxed);
                            Vec::new()
                        },
                        |scratch: &mut Vec<usize>, job| {
                            runs.fetch_add(1, Ordering::Relaxed);
                            scratch.push(job);
                            let nested: u64 =
                                (0..64_u64).into_par_iter().map(|x| x * job as u64).sum();
                            (job, nested, rayon::current_thread_index().is_some())
                        },
                    );
                    let expected: Vec<(usize, u64, bool)> = (0..jobs)
                        .map(|job| (job, 2016 * job as u64, true))
                        .collect();
                    assert_eq!(results, expected, "threads {threads} jobs {jobs}");
                    assert_eq!(runs.load(Ordering::Relaxed), jobs);
                    let made = scratches.load(Ordering::Relaxed);
                    assert!(
                        (1..=tasks).contains(&made),
                        "{made} scratches for {tasks} tasks"
                    );
                }
                Ok(())
            })
            .expect("a pool");
        }
    }

    /// A job that panics: the panic comes out of the scope (the other tasks finish first, so
    /// nothing hangs) and the same pool computes the next call correctly.
    #[test]
    fn a_panicking_job_propagates_and_the_pool_stays_usable() {
        for threads in [2, 3, 7] {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                let jobs = 40_usize;
                let order: Vec<usize> = (0..jobs).collect();
                let run = |panicking: bool| {
                    x_compute_window_lpt(
                        pool,
                        threads.min(jobs),
                        &order,
                        || (),
                        |_: &mut (), job| {
                            assert!(!(panicking && (job == 5 || job == 30)), "job {job} panics");
                            job * 3
                        },
                    )
                };
                let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| run(true)));
                assert!(caught.is_err(), "threads {threads}: the panic was lost");
                let expected: Vec<usize> = (0..jobs).map(|job| job * 3).collect();
                assert_eq!(run(false), expected, "threads {threads}");
                Ok(())
            })
            .expect("a pool");
        }
    }

    const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

    /// A redo's inputs made up from related proteins: random queries; subjects that hold a
    /// mutated piece of one query between random flanks; per query a preliminary hit list with a
    /// list for each related subject (one HSP over the piece, its gapped start in the middle), a
    /// few lists of unrelated subjects, and one list without HSPs; one query without a list.
    struct XRedoFixture {
        contexts: Vec<QueryContext>,
        search_spaces: Vec<SearchSpace>,
        subjects: Vec<EncodedProtein>,
        lists: Vec<Option<BlastpHitList>>,
        batches: Vec<Range<usize>>,
        evalue_adjust_contexts: Vec<crate::algorithm::blastp::kappa::EvalueAdjustContext>,
        gapped_params: KarlinParams,
        gapped_gumbel: BlastGumbelBlk,
        matrix_info: crate::core::composition_adjustment::adjust_scores::BlastMatrixInfo,
    }

    /// A list of subject `oid` for query `q_idx` (of `query_length` letters) with one HSP of
    /// `length` letters from query offset `q0` and subject offset `s0` (0-based).
    fn aligned_list(
        q_idx: usize,
        oid: usize,
        (q0, s0, length): (usize, usize, usize),
        query_length: usize,
        e_value: f64,
    ) -> BlastpHspList {
        let mut hsp = hsp(oid as u32, e_value);
        hsp.q_start = q0 + 1;
        hsp.q_end = q0 + length;
        hsp.s_start = s0 + 1;
        hsp.s_end = s0 + length;
        hsp.gapped_q_start = (q0 + length / 2) as i32;
        hsp.gapped_s_start = (s0 + length / 2) as i32;
        hsp.length = length;
        hsp.raw_score = 2 * length as i32;
        hsp.query_length = query_length;
        hsp.q_idx = q_idx as u32;
        BlastpHspList {
            oid: oid as u32,
            query_index: q_idx as u32,
            hsps: vec![hsp],
            best_evalue: e_value,
        }
    }

    fn redo_fixture() -> XRedoFixture {
        let spec = crate::config::ProteinScoringSpec::default();
        let gapped_params = crate::stats::tables::lookup_protein_params(&spec);
        let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
        let random = |rng: &mut XorShift, len: usize| -> Vec<u8> {
            (0..len).map(|_| AMINO_ACIDS[rng.below(20)]).collect()
        };
        let queries: Vec<Vec<u8>> = (0..6)
            .map(|_| {
                let len = 120 + rng.below(240);
                random(&mut rng, len)
            })
            .collect();
        let mut texts = Vec::new();
        let mut related = Vec::new();
        for oid in 0..48 {
            let q_idx = rng.below(queries.len());
            let query = &queries[q_idx];
            let length = 50 + rng.below(query.len() - 60);
            let q0 = rng.below(query.len() - length + 1);
            let left = rng.below(80);
            let mut subject = random(&mut rng, left);
            for &residue in &query[q0..q0 + length] {
                // About a quarter of the letters changed.
                subject.push(if rng.below(4) == 0 {
                    AMINO_ACIDS[rng.below(20)]
                } else {
                    residue
                });
            }
            let right = rng.below(80);
            subject.extend(random(&mut rng, right));
            related.push((q_idx, oid, (q0, left, length)));
            texts.push(subject);
        }
        let mut per_query: Vec<Vec<BlastpHspList>> = vec![Vec::new(); queries.len()];
        for &(q_idx, oid, piece) in &related {
            let e_value = 10_f64.powi(-(3 + rng.below(40) as i32));
            per_query[q_idx].push(aligned_list(
                q_idx,
                oid,
                piece,
                queries[q_idx].len(),
                e_value,
            ));
        }
        for (q_idx, lists) in per_query.iter_mut().enumerate() {
            let unrelated: Vec<usize> = related
                .iter()
                .filter(|&&(other, ..)| other != q_idx)
                .map(|&(_, oid, _)| oid)
                .collect();
            for &oid in unrelated.iter().take(3) {
                let length = 30;
                let q0 = rng.below(queries[q_idx].len() - length);
                let s0 = rng.below(texts[oid].len() - length);
                lists.push(aligned_list(
                    q_idx,
                    oid,
                    (q0, s0, length),
                    queries[q_idx].len(),
                    0.5,
                ));
            }
            lists.push(BlastpHspList {
                oid: unrelated[3] as u32,
                query_index: q_idx as u32,
                hsps: Vec::new(),
                best_evalue: 1.0,
            });
            lists.sort_by(|a, b| a.best_evalue.total_cmp(&b.best_evalue));
        }
        let mut lists: Vec<Option<BlastpHitList>> = per_query.into_iter().map(hit_list).collect();
        lists[4] = None;
        let contexts: Vec<QueryContext> = queries
            .iter()
            .enumerate()
            .map(|(q_idx, query)| {
                let encoded = super::super::encoding::encode_protein_sequence(query);
                QueryContext {
                    q_idx: q_idx as u32,
                    frame: 0,
                    aa_seq: encoded.aa_seq,
                    aa_seq_nomask: None,
                    aa_len: encoded.aa_len,
                    orig_len: encoded.aa_len,
                    frame_base: 0,
                    is_valid: true,
                    karlin_params: gapped_params,
                }
            })
            .collect();
        let subjects: Vec<EncodedProtein> = texts
            .iter()
            .map(|text| super::super::encoding::encode_protein_sequence(text))
            .collect();
        let db_letters: usize = subjects.iter().map(|subject| subject.aa_len).sum();
        let search_spaces: Vec<SearchSpace> = contexts
            .iter()
            .map(|ctx| {
                let length_adjustment = 20_i64;
                let effective_query_len = (ctx.aa_len as i64 - length_adjustment) as f64;
                let effective_db_len =
                    (db_letters as i64 - subjects.len() as i64 * length_adjustment) as f64;
                SearchSpace {
                    effective_space: effective_query_len * effective_db_len,
                    length_adjustment,
                }
            })
            .collect();
        let ideal = crate::stats::compute_blosum62_ideal_karlin_params()
            .expect("BLOSUM62 ideal parameters");
        let batches = vec![0..3, 0..3, 0..3, 3..6, 3..6, 3..6];
        let evalue_adjust_contexts = crate::algorithm::blastp::kappa::evalue_adjust_contexts(
            &contexts,
            &search_spaces,
            &batches,
            0,
            None,
        )
        .expect("no refusal without unified P-values");
        XRedoFixture {
            contexts,
            search_spaces,
            subjects,
            lists,
            batches,
            evalue_adjust_contexts,
            gapped_params,
            gapped_gumbel: crate::stats::spouge::lookup_protein_gumbel_params(
                &spec,
                db_letters as i64,
            )
            .expect("BLOSUM62 11/1 Gumbel parameters"),
            matrix_info: crate::core::composition_adjustment::adjust_scores::build_matrix_info(
                ScoringMatrix::Blosum62,
                ideal.lambda / 32.0,
            )
            .expect("BLOSUM62 matrix information"),
        }
    }

    /// The redo parameters of blastp's defaults (composition matrix adjustment, scaled by 32), as
    /// `run_resolved_in_pool` builds them.
    fn redo_params(fixture: &XRedoFixture, q_idx: usize) -> BlastRedoAlignParams {
        use crate::core::composition_adjustment::redo_alignment::{
            BlastCompoAdjustMode, BlastCompoGappingParams,
        };
        const SCALE: f64 = 32.0;
        let scaled_lambda = fixture.gapped_params.lambda / SCALE;
        BlastRedoAlignParams {
            matrix_info: fixture.matrix_info.clone(),
            gapping_params: BlastCompoGappingParams {
                gap_open: 11 * SCALE as i32,
                gap_extend: SCALE as i32,
                x_dropoff: (25.0 * std::f64::consts::LN_2 / scaled_lambda) as i32,
                context: std::cell::Cell::new(None),
            },
            compo_adjust_mode: BlastCompoAdjustMode::CompositionMatrixAdjust,
            alphsize: crate::utils::matrix::BLASTAA_SIZE as i32,
            composition_test_index: 0,
            log_k: fixture.gapped_params.k.glibc_ln(),
            score_divisor: SCALE,
            smith_waterman: false,
            near_identical_cutoff: 1.74 * std::f64::consts::LN_2 / scaled_lambda,
            position_based: false,
            ccat_query_length: fixture.contexts[q_idx].aa_len as i32,
            query_is_translated: false,
            subject_is_translated: false,
            cutoff_score: 1,
            do_link_hsps: false,
        }
    }

    fn query_raw(ctx: &QueryContext) -> &[u8] {
        &ctx.aa_seq[1..ctx.aa_seq.len() - 1]
    }

    /// `x_compute_window` with the switch gives the results of the `par_iter` path, match for
    /// match (error or not, scores, E-value bits, every HSP field), on windows from random starts
    /// with random job budgets and matches ended at random, at 1, 2, 3 and 7 threads, with the
    /// search's SEG store, and when the parameters of one query fail.
    #[test]
    fn windows_computed_with_jobs_largest_first_equal_the_par_iter_windows() {
        let fixture = redo_fixture();
        let count = fixture.lists.len();
        let heaps = heaps_of(count, 500);
        let all = stream(&fixture.lists);
        let mut rng = XorShift(0x5851_f42d_4c95_7f2d);
        let windows: Vec<XWindow> = (0..10)
            .map(|index| {
                let start = all[rng.below(all.len())];
                let mut window = window_from(
                    &fixture.lists,
                    &heaps,
                    &fixture.batches,
                    &vec![None; count],
                    start,
                    2 + rng.below(10),
                );
                // None ended, about one in four, about one in two.
                let share = index % 3;
                for ended in &mut window.ended {
                    *ended = share > 0 && rng.below(4) < share;
                }
                window
            })
            .collect();
        let failing = all[all.len() / 2].q_idx;
        let (mut compared, mut with_list, mut errors, mut lpt_windows) = (0, 0, 0, 0);
        for threads in [1, 2, 3, 7] {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                for fail in [None, Some(failing)] {
                    let build = |q_idx: usize| -> Result<BlastRedoAlignParams> {
                        if Some(q_idx) == fail {
                            anyhow::bail!("no parameters for query {q_idx}");
                        }
                        Ok(redo_params(&fixture, q_idx))
                    };
                    let search = XRedoParSearch {
                        pool,
                        threads,
                        contexts: &fixture.contexts,
                        query_nomask: query_raw,
                        search_spaces: &fixture.search_spaces,
                        subjects: crate::algorithm::blastp::x_redo_subjects::XRedoSubjects::Encoded(
                            &fixture.subjects,
                        ),
                        query_batches: &fixture.batches,
                        evalue_adjust_contexts: &fixture.evalue_adjust_contexts,
                        matrix: ScoringMatrix::Blosum62,
                        gapped_params: &fixture.gapped_params,
                        gapped_gumbel: Some(&fixture.gapped_gumbel),
                        evalue: 10.0,
                        build_redo_align_params: &build,
                        timing: false,
                        engine_positives: false,
                        sw_errors: TEST_SW_ERRORS,
                    };
                    let shared = || Arc::new(XSharedSubjectRanges::new());
                    for window in &windows {
                        let (reference, ran) = x_compute_window::<_, false>(
                            &search,
                            &fixture.lists,
                            &heaps,
                            window,
                            &shared(),
                            false,
                        );
                        assert!(!ran, "the par_iter path");
                        let (lpt, ran) = x_compute_window::<_, false>(
                            &search,
                            &fixture.lists,
                            &heaps,
                            window,
                            &shared(),
                            true,
                        );
                        let acts = threads > 1
                            && window.jobs.len() > 1
                            && !window.ended.iter().all(|&ended| ended);
                        assert_eq!(ran, acts, "threads {threads}: {window:?}");
                        lpt_windows += usize::from(ran);
                        assert_eq!(lpt.len(), reference.len());
                        for (position, (computed, reference)) in
                            lpt.iter().zip(reference).enumerate()
                        {
                            let at = (threads, fail, window.matches[position]);
                            match (computed, reference) {
                                (None, None) => {}
                                (Some(computed), Some(reference)) => {
                                    let reference = reference.result;
                                    with_list += usize::from(
                                        reference
                                            .as_ref()
                                            .is_ok_and(|result| result.hsp_list.is_some()),
                                    );
                                    errors += usize::from(reference.is_err());
                                    let reference = reference.map_err(|error| format!("{error:#}"));
                                    assert_eq!(
                                        x_result_difference(&computed.result, &reference),
                                        None,
                                        "{at:?}"
                                    );
                                    assert!(computed.deferred.is_none(), "{at:?}");
                                    compared += 1;
                                }
                                (computed, reference) => panic!(
                                    "{at:?}: computed {} / {}",
                                    computed.is_some(),
                                    reference.is_some()
                                ),
                            }
                        }
                    }
                }
                Ok(())
            })
            .expect("a pool");
        }
        // Real work on every side: matches computed, alignments kept, the failing query's errors,
        // and windows the switch computed.
        assert!(
            compared > 200 && with_list > 50 && errors > 0 && lpt_windows > 30,
            "compared {compared} with a list {with_list} errors {errors} lpt windows {lpt_windows}"
        );
    }

    /// EXPERIMENT (LOSAT_X_BPSEGWIN): windows computed with the search's SEG store emptied after
    /// each one (`XSharedSubjectRanges::x_clear`) give the results of the same windows computed
    /// with one store kept for the whole search, match for match (error or not, scores, E-value
    /// bits, every HSP field), on windows from random starts with random job budgets and matches
    /// ended at random, at 1, 2 and 7 threads, and when the parameters of one query fail; the
    /// later windows mask ranges again that an earlier window had masked.
    #[test]
    fn windows_with_the_seg_store_emptied_after_each_equal_windows_with_one_store() {
        let fixture = redo_fixture();
        let count = fixture.lists.len();
        let heaps = heaps_of(count, 500);
        let all = stream(&fixture.lists);
        let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
        let windows: Vec<XWindow> = (0..10)
            .map(|index| {
                let start = all[rng.below(all.len())];
                let mut window = window_from(
                    &fixture.lists,
                    &heaps,
                    &fixture.batches,
                    &vec![None; count],
                    start,
                    2 + rng.below(10),
                );
                // None ended, about one in four, about one in two.
                let share = index % 3;
                for ended in &mut window.ended {
                    *ended = share > 0 && rng.below(4) < share;
                }
                window
            })
            .collect();
        let failing = all[all.len() / 2].q_idx;
        let (mut compared, mut with_list, mut errors, mut masked_again) = (0, 0, 0, 0);
        for threads in [1, 2, 7] {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                for fail in [None, Some(failing)] {
                    let build = |q_idx: usize| -> Result<BlastRedoAlignParams> {
                        if Some(q_idx) == fail {
                            anyhow::bail!("no parameters for query {q_idx}");
                        }
                        Ok(redo_params(&fixture, q_idx))
                    };
                    let search = XRedoParSearch {
                        pool,
                        threads,
                        contexts: &fixture.contexts,
                        query_nomask: query_raw,
                        search_spaces: &fixture.search_spaces,
                        subjects: crate::algorithm::blastp::x_redo_subjects::XRedoSubjects::Encoded(
                            &fixture.subjects,
                        ),
                        query_batches: &fixture.batches,
                        evalue_adjust_contexts: &fixture.evalue_adjust_contexts,
                        matrix: ScoringMatrix::Blosum62,
                        gapped_params: &fixture.gapped_params,
                        gapped_gumbel: Some(&fixture.gapped_gumbel),
                        evalue: 10.0,
                        build_redo_align_params: &build,
                        timing: false,
                        engine_positives: false,
                        sw_errors: TEST_SW_ERRORS,
                    };
                    let kept = Arc::new(XSharedSubjectRanges::new());
                    let emptied = Arc::new(XSharedSubjectRanges::new());
                    let mut window_entries = 0;
                    for window in &windows {
                        let (reference, _) = x_compute_window::<_, false>(
                            &search,
                            &fixture.lists,
                            &heaps,
                            window,
                            &kept,
                            false,
                        );
                        let (computed, _) = x_compute_window::<_, false>(
                            &search,
                            &fixture.lists,
                            &heaps,
                            window,
                            &emptied,
                            false,
                        );
                        window_entries += emptied.x_len();
                        emptied.x_clear();
                        assert_eq!(emptied.x_len(), 0);
                        assert_eq!(computed.len(), reference.len());
                        for (position, (computed, reference)) in
                            computed.iter().zip(reference).enumerate()
                        {
                            let at = (threads, fail, window.matches[position]);
                            match (computed, reference) {
                                (None, None) => {}
                                (Some(computed), Some(reference)) => {
                                    let reference = reference.result;
                                    with_list += usize::from(
                                        reference
                                            .as_ref()
                                            .is_ok_and(|result| result.hsp_list.is_some()),
                                    );
                                    errors += usize::from(reference.is_err());
                                    let reference = reference.map_err(|error| format!("{error:#}"));
                                    assert_eq!(
                                        x_result_difference(&computed.result, &reference),
                                        None,
                                        "{at:?}"
                                    );
                                    assert!(computed.deferred.is_none(), "{at:?}");
                                    compared += 1;
                                }
                                (computed, reference) => panic!(
                                    "{at:?}: computed {} / {}",
                                    computed.is_some(),
                                    reference.is_some()
                                ),
                            }
                        }
                    }
                    // The entries the emptied store made beyond the kept store's: ranges masked
                    // again by a later window.
                    masked_again += window_entries - kept.x_len();
                }
                Ok(())
            })
            .expect("a pool");
        }
        assert!(
            compared > 150 && with_list > 40 && errors > 0 && masked_again > 0,
            "compared {compared} with a list {with_list} errors {errors} masked again {masked_again}"
        );
    }

    // ---- The Smith-Waterman redo in windows (`-use_sw_tback`, item X4-mt). ----

    /// What the Smith-Waterman redo of a synthetic match gives: its first round's E-value and
    /// score (its heap test, when the test is asked) and the outcome after a TRUE test.
    #[derive(Clone, Copy, Debug, PartialEq)]
    enum SwKind {
        /// The first round's score is not significant: no test, no alignment.
        NotAsked,
        /// Alignments whose best E-value and score are the first round's.
        Kept,
        /// No alignment left after the redo (all above `-evalue`).
        Purged,
        /// A failed redo (`SwRedoFailed`).
        Failed,
        /// A redo that does not end (`SwRedoEndless`).
        Endless,
    }

    #[derive(Clone, Copy, Debug)]
    struct SwCase {
        evalue: f64,
        score: i32,
        kind: SwKind,
    }

    fn no_alignment() -> Result<BlastpPostprocessResult> {
        Ok(BlastpPostprocessResult {
            hsp_list: None,
            best_score: 0,
            best_evalue: f64::MAX,
        })
    }

    /// A job's computation of a synthetic match: the heap test answered as `test` says (as
    /// `x_redo_match_sw` answers it), the outcome of its case after a TRUE answer.
    fn sw_job_match(test: XSwTest<'_>, oid: u32, case: SwCase) -> XComputed {
        let mut deferred = None;
        let answer = |deferred: &mut Option<(f64, i32)>| match test {
            XSwTest::Heap(heap, before) => {
                x_sw_heap_test(heap, before, (case.evalue, case.score, oid), deferred)
            }
            XSwTest::Continued(asked) => {
                assert_eq!(asked, (case.evalue, case.score), "the continued test");
                true
            }
        };
        let result = if case.kind == SwKind::NotAsked || !answer(&mut deferred) {
            no_alignment()
        } else {
            match case.kind {
                SwKind::Kept => Ok(BlastpPostprocessResult {
                    hsp_list: Some(list(oid, case.evalue, false)),
                    best_score: case.score,
                    best_evalue: case.evalue,
                }),
                SwKind::Purged => no_alignment(),
                SwKind::Failed => Err(SwRedoFailed.into()),
                SwKind::Endless => Err(SwRedoEndless.into()),
                SwKind::NotAsked => unreachable!(),
            }
        };
        XComputed {
            result,
            deferred,
            first_round_s: 0.0,
        }
    }

    /// The windows of the Smith-Waterman redo on synthetic matches (`cases[q][position]`), each
    /// job's matches computed by `sw_job_match`; returns the redo's result and the matches
    /// computed.
    fn sw_drive(
        lists: &mut [Option<BlastpHitList>],
        cases: &[Vec<SwCase>],
        heaps: &mut [BlastCompoHeap],
        batches: &[Range<usize>],
        limits: XWindowLimits,
    ) -> (Result<XCounts>, usize) {
        let mut computed = 0;
        let counts = x_redo_windows::<true>(
            lists,
            heaps,
            batches,
            limits,
            false,
            None,
            TEST_SW_ERRORS,
            &mut |lists, heaps, window| {
                window
                    .jobs
                    .iter()
                    .flat_map(|job| job.clone())
                    .map(|index| {
                        let matched = window.matches[index];
                        (!window.ended[index]).then(|| {
                            computed += 1;
                            let first = window
                                .matches
                                .partition_point(|other| other.q_idx < matched.q_idx);
                            let test = match window
                                .continued
                                .binary_search_by_key(&index, |&(at, _)| at)
                            {
                                Ok(found) => XSwTest::Continued(window.continued[found].1),
                                Err(_) => XSwTest::Heap(&heaps[matched.q_idx], index - first),
                            };
                            sw_job_match(
                                test,
                                x_local_match(lists, matched).oid,
                                cases[matched.q_idx][matched.position],
                            )
                        })
                    })
                    .collect()
            },
        );
        (counts, computed)
    }

    /// NCBI's one-thread Smith-Waterman redo of the synthetic matches (the steps of
    /// `redo_run::one_thread_redo`): the early termination test, the heap test (a FALSE ends the
    /// match with status 0), the outcome, the insert, the batch's status and its check at the
    /// batch's last query.
    fn sw_one_thread(
        lists: &[Option<BlastpHitList>],
        cases: &[Vec<SwCase>],
        heaps: &mut [BlastCompoHeap],
        batches: &[Range<usize>],
    ) -> Result<()> {
        let mut failed = false;
        for q_idx in 0..lists.len() {
            let batch = batches[q_idx].clone();
            if batch.start == q_idx {
                failed = false;
            }
            if let Some(hit_list) = lists[q_idx].as_ref() {
                for (position, local_match) in hit_list.hsplist_array[..hit_list.hsplist_count]
                    .iter()
                    .enumerate()
                {
                    if local_match.hsps.is_empty()
                        || blast_compo_early_termination(
                            local_match.best_evalue,
                            &heaps[batch.clone()],
                        )
                    {
                        continue;
                    }
                    let (case, oid) = (cases[q_idx][position], local_match.oid);
                    if case.kind == SwKind::NotAsked
                        || !heaps[q_idx].would_insert(case.evalue, case.score, oid)
                    {
                        failed = false;
                        continue;
                    }
                    match case.kind {
                        SwKind::Kept => {
                            failed = false;
                            if heaps[q_idx].would_insert(case.evalue, case.score, oid) {
                                heaps[q_idx].insert(
                                    list(oid, case.evalue, false),
                                    case.evalue,
                                    case.score,
                                    oid,
                                );
                            }
                        }
                        SwKind::Purged => failed = false,
                        SwKind::Failed => failed = true,
                        SwKind::Endless => return Err(test_endless(batch.start)),
                        SwKind::NotAsked => unreachable!(),
                    }
                }
            }
            if failed && q_idx + 1 == batch.end {
                return Err(test_batch_failed(batch.start, &heaps[batch.start]));
            }
        }
        Ok(())
    }

    /// Each heap as its pops give it (the subjects, worst first), with its count and worst E.
    fn heap_contents(heaps: &[BlastCompoHeap]) -> Vec<(usize, u64, Vec<u32>)> {
        heaps
            .iter()
            .map(|heap| {
                let mut heap = heap.clone();
                let (n, worst) = (heap.n, heap.worst_evalue.to_bits());
                let mut oids = Vec::new();
                while let Some(list) = heap.pop() {
                    oids.push(list.oid);
                }
                (n, worst, oids)
            })
            .collect()
    }

    fn error_text(result: &Result<()>) -> Option<String> {
        result.as_ref().err().map(|error| format!("{error:#}"))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:251-276
    #[test]
    fn the_job_answers_the_heap_test_where_it_is_certain_and_defers_the_rest() {
        let mut heap = BlastCompoHeap::new(2, ECUTOFF);
        let mut deferred = None;
        let mut test = |heap: &BlastCompoHeap, before, evalue, score, oid| {
            deferred = None;
            let answer = x_sw_heap_test(heap, before, (evalue, score, oid), &mut deferred);
            (answer, deferred)
        };
        // Room for the matches before: TRUE for certain, whatever the E-value.
        assert_eq!(test(&heap, 1, 5.0, 30, 9), (true, None));
        // Within the cutoff: TRUE for certain, however full the heap may be.
        assert_eq!(test(&heap, 7, ECUTOFF, 30, 9), (true, None));
        // The heap may be full by then and the empty heap's test says TRUE: deferred (the job's
        // redo stops at the test, as on FALSE).
        assert_eq!(test(&heap, 2, 5.0, 30, 9), (false, Some((5.0, 30))));
        heap.insert(list(1, 0.5, false), 0.5, 40, 1);
        heap.insert(list(2, 0.01, false), 0.01, 50, 2);
        // Full above the cutoff: better than the root (E 0.5) is deferred, worse is FALSE.
        assert_eq!(test(&heap, 0, 0.1, 30, 9), (false, Some((0.1, 30))));
        assert_eq!(test(&heap, 0, 0.6, 30, 9), (false, None));
        // A tie of E-values: the higher score is better, the lower is not.
        assert_eq!(test(&heap, 0, 0.5, 41, 9), (false, Some((0.5, 41))));
        assert_eq!(test(&heap, 0, 0.5, 39, 9), (false, None));
        // Without side effect: the array was not turned into a heap.
        assert!(heap.heap_array.is_none());
    }

    #[test]
    fn smith_waterman_jobs_are_at_risk_only_as_x_drop_jobs_are() {
        // One query, six matches, two lists to keep, no match the early termination test may
        // end. The Smith-Waterman jobs' tests of the third to sixth matches are deferred (the heap
        // may hold two by then), which wastes nothing: neither redo has a job at risk, and a cap
        // of one job at risk does not end the window.
        let lists = vec![hit_list((0..6).map(|oid| list(oid, 1e-9, false)).collect())];
        let heaps = heaps_of(1, 2);
        let none = [None];
        for caps in [vec![usize::MAX], vec![1]] {
            let xdrop =
                x_window::<false>(&lists, &heaps, &alone(1), &none, &caps, None, at(0, 0), 8);
            let sw = x_window::<true>(&lists, &heaps, &alone(1), &none, &caps, None, at(0, 0), 8);
            assert_eq!((sw.jobs.len(), sw.at_risk_jobs()), (2, 0));
            assert_eq!(sw, xdrop);
        }
    }

    fn sw_case(evalue: f64, kind: SwKind) -> SwCase {
        SwCase {
            evalue,
            score: 40,
            kind,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1463-1476
    /// A heap test that the job defers and the replay answers FALSE ends the match as NCBI ends it
    /// (no alignment, status 0): the job's redo stopped at the test, so the failed or endless redo
    /// after it is never computed; where the test is TRUE for certain, the same redos fail the
    /// batch and end the search.
    #[test]
    fn a_deferred_heap_test_that_turns_false_masks_a_failed_and_an_endless_redo() {
        let lists = || vec![hit_list((0..3).map(|oid| list(oid, 1e-9, false)).collect())];
        let cases = vec![vec![
            sw_case(1e-20, SwKind::Kept),
            sw_case(0.5, SwKind::Failed),
            sw_case(0.5, SwKind::Endless),
        ]];
        let one_window = XWindowLimits {
            jobs: 8,
            at_risk_jobs: 8,
        };
        for (heap_threshold, expected) in [(1, None), (3, Some("endless in batch 0".to_string()))] {
            let mut heaps = heaps_of(1, heap_threshold);
            let (counts, computed) =
                sw_drive(&mut lists(), &cases, &mut heaps, &alone(1), one_window);
            assert_eq!(computed, 3);
            let mut model = heaps_of(1, heap_threshold);
            let reference = sw_one_thread(&lists(), &cases, &mut model, &alone(1));
            assert_eq!(
                error_text(&reference),
                expected,
                "threshold {heap_threshold}"
            );
            match counts {
                Ok(counts) => {
                    assert_eq!(expected, None);
                    // Both later tests were deferred and turned FALSE.
                    assert_eq!((counts.turned, counts.continued), (Some(2), 0));
                }
                Err(error) => assert_eq!(Some(format!("{error:#}")), expected),
            }
            assert_eq!(heap_contents(&heaps), heap_contents(&model));
        }
    }

    /// A heap test that the job defers and the replay answers TRUE stops the replay at its match:
    /// the next window's first job finishes the match (a continuation) and the carried matches
    /// after it are replayed with their results, before the jobs that follow them (which count
    /// them as the query's matches before theirs). Every match is computed once, and each
    /// continued one once more.
    #[test]
    fn a_deferred_heap_test_that_turns_true_is_finished_by_a_continuation_in_the_next_window() {
        // One query, one list to keep. The first match's redo asks no test (no insert), so the
        // second's deferred test is TRUE; the third is better than the second (TRUE); the others
        // are worse (FALSE), the last a failed redo.
        let lists = || {
            vec![hit_list(
                (0..10).map(|oid| list(oid, 1e-9, false)).collect(),
            )]
        };
        let mut query = vec![
            sw_case(0.5, SwKind::NotAsked),
            sw_case(0.5, SwKind::Kept),
            sw_case(0.4, SwKind::Kept),
        ];
        query.extend((3..9).map(|_| sw_case(0.6, SwKind::Kept)));
        query.push(sw_case(0.7, SwKind::Failed));
        let cases = vec![query];
        // One job per window: the fourth match is carried once, the later ones are FALSE for
        // certain at their windows. Two: the second window holds the continuation, the carried
        // matches and the last job. Eight: the first window holds all.
        for (jobs, turned) in [(1, 1), (2, 7), (8, 7)] {
            let limits = XWindowLimits {
                jobs,
                at_risk_jobs: jobs,
            };
            let mut heaps = heaps_of(1, 1);
            let (counts, computed) = sw_drive(&mut lists(), &cases, &mut heaps, &alone(1), limits);
            let counts = counts.expect("no error");
            let mut model = heaps_of(1, 1);
            sw_one_thread(&lists(), &cases, &mut model, &alone(1)).expect("no error");
            assert_eq!(heap_contents(&heaps), heap_contents(&model), "jobs {jobs}");
            assert_eq!(
                (counts.turned, counts.continued, computed),
                (Some(turned), 2, 12),
                "jobs {jobs}"
            );
        }
        // The deferred TRUE at a window's last match (nothing carried after it) and at the
        // stream's last match (the last window holds the continuation only).
        let lists = || {
            vec![
                hit_list(vec![list(0, 1e-9, false), list(1, 1e-9, false)]),
                hit_list(vec![list(2, 1e-9, false), list(3, 1e-9, false)]),
            ]
        };
        let pair = vec![sw_case(0.5, SwKind::NotAsked), sw_case(0.5, SwKind::Kept)];
        let cases = vec![pair.clone(), pair];
        for jobs in [1, 8] {
            let limits = XWindowLimits {
                jobs,
                at_risk_jobs: jobs,
            };
            let mut heaps = heaps_of(2, 1);
            let (counts, computed) = sw_drive(&mut lists(), &cases, &mut heaps, &alone(2), limits);
            let counts = counts.expect("no error");
            let mut model = heaps_of(2, 1);
            sw_one_thread(&lists(), &cases, &mut model, &alone(2)).expect("no error");
            assert_eq!(heap_contents(&heaps), heap_contents(&model), "jobs {jobs}");
            assert_eq!(
                (counts.turned, counts.continued, computed),
                (Some(0), 2, 6),
                "jobs {jobs}"
            );
        }
        // Two queries, each a match that inserts its list and a better one whose deferred test
        // is TRUE. At the stop in the first query, the second query's carried matches are
        // replayed ahead on a copy of its heap: its test is TRUE there, so both continuations run
        // in the second window.
        let cases = vec![vec![sw_case(0.5, SwKind::Kept), sw_case(0.4, SwKind::Kept)]; 2];
        let limits = XWindowLimits {
            jobs: 8,
            at_risk_jobs: 8,
        };
        let mut heaps = heaps_of(2, 1);
        let (counts, computed) = sw_drive(&mut lists(), &cases, &mut heaps, &alone(2), limits);
        let counts = counts.expect("no error");
        let mut model = heaps_of(2, 1);
        sw_one_thread(&lists(), &cases, &mut model, &alone(2)).expect("no error");
        assert_eq!(heap_contents(&heaps), heap_contents(&model));
        assert_eq!(
            (counts.turned, counts.continued, counts.windows, computed),
            (Some(0), 2, 2, 6)
        );
    }

    /// A continuation is replayed before the end of its batch: its redo sets the batch's status
    /// (clearing a failure before it, or failing the batch) or ends the search (an endless redo),
    /// before any match of a later batch.
    #[test]
    fn a_continuation_sets_its_batch_status_before_the_batch_ends() {
        let one_job = XWindowLimits {
            jobs: 1,
            at_risk_jobs: 1,
        };
        let one_window = XWindowLimits {
            jobs: 64,
            at_risk_jobs: 64,
        };
        let lists = || {
            vec![
                hit_list(vec![list(0, 1e-9, false), list(1, 1e-9, false)]),
                hit_list(vec![list(2, 1e-9, false)]),
            ]
        };
        let later = vec![sw_case(1e-20, SwKind::Endless)];
        for (first, second, expected) in [
            // A failed redo, then the continuation's redo: status 0, the later batch ends it.
            (
                sw_case(1e-20, SwKind::Failed),
                sw_case(0.5, SwKind::Kept),
                "endless in batch 1",
            ),
            // The continuation's redo fails: the batch fails at its end.
            (
                sw_case(0.5, SwKind::NotAsked),
                sw_case(0.5, SwKind::Failed),
                "batch 0 failed, first query holds 0",
            ),
            // The continuation's redo never ends.
            (
                sw_case(0.5, SwKind::NotAsked),
                sw_case(0.5, SwKind::Endless),
                "endless in batch 0",
            ),
        ] {
            let cases = vec![vec![first, second], later.clone()];
            for limits in [one_job, one_window] {
                let mut heaps = heaps_of(2, 1);
                let (result, computed) =
                    sw_drive(&mut lists(), &cases, &mut heaps, &alone(2), limits);
                let mut model = heaps_of(2, 1);
                let reference = sw_one_thread(&lists(), &cases, &mut model, &alone(2));
                let result = result.map(|_| ());
                assert_eq!(error_text(&reference), Some(expected.to_string()));
                assert_eq!(error_text(&result), error_text(&reference), "{limits:?}");
                assert_eq!(heap_contents(&heaps), heap_contents(&model), "{limits:?}");
                // The second match's job and its continuation.
                assert!(computed >= 3, "{limits:?}");
            }
        }
    }

    /// The batch's status crosses windows: a failed redo at the end of a batch fails it when the
    /// replay passes the batch's end (in a later window, or after the last one); a later match of
    /// the batch that the redo reaches clears it, one that the early termination test ends does
    /// not; and a failed batch comes before an endless redo of a later batch computed ahead.
    #[test]
    fn the_batch_status_crosses_windows_and_fails_the_batch_at_its_end() {
        let one_job = XWindowLimits {
            jobs: 1,
            at_risk_jobs: 1,
        };
        let one_window = XWindowLimits {
            jobs: 64,
            at_risk_jobs: 64,
        };
        // Batch 0..2: query 0 ends with a failed redo; query 1's matches are ended by the test
        // (its heap and query 0's filled to the cutoff, E 1.0) or reached (E 1e-9).
        let run = |second: f64, limits: XWindowLimits, batches: Vec<Range<usize>>| {
            let lists = vec![
                hit_list(vec![list(1, 1e-9, false), list(2, 1e-9, false)]),
                hit_list(vec![list(3, 1e-9, false), list(4, second, false)]),
                hit_list(vec![list(5, 1e-9, false)]),
            ];
            let cases = vec![
                vec![sw_case(1e-20, SwKind::Kept), sw_case(1e-20, SwKind::Failed)],
                vec![sw_case(1e-20, SwKind::Kept), sw_case(1e-20, SwKind::Kept)],
                vec![sw_case(1e-20, SwKind::Endless)],
            ];
            let mut heaps = heaps_of(3, 1);
            let (result, _) = sw_drive(&mut lists.clone(), &cases, &mut heaps, &batches, limits);
            let mut model = heaps_of(3, 1);
            let reference = sw_one_thread(&lists, &cases, &mut model, &batches);
            let result = result.map(|_| ());
            assert_eq!(error_text(&result), error_text(&reference));
            assert_eq!(heap_contents(&heaps), heap_contents(&model));
            error_text(&result)
        };
        let failed = Some("batch 0 failed, first query holds 1".to_string());
        let endless = Some("endless in batch 2".to_string());
        for limits in [one_job, one_window] {
            // Query 1's first match is reached and clears the status; batch 2's redo never ends.
            assert_eq!(run(1e-9, limits, vec![0..2, 0..2, 2..3]), endless);
            // Query 1 alone in a batch: query 0's failed status fails batch 0 before batch 1.
            assert_eq!(run(1e-9, limits, vec![0..1, 1..2, 2..3]), failed);
            // All three in one batch: the endless redo of query 2 ends it.
            assert_eq!(
                run(1e-9, limits, vec![0..3; 3]),
                Some("endless in batch 0".to_string())
            );
        }
        // A match ended by the test leaves the status that the match before it set: query 1's
        // third match (E 1.0) comes after both heaps are filled to the cutoff.
        let lists = vec![
            hit_list(vec![list(1, 1e-9, false), list(2, 1e-9, false)]),
            hit_list(vec![
                list(3, 1e-9, false),
                list(4, 1e-9, false),
                list(5, 1.0, false),
            ]),
        ];
        let cases = vec![
            vec![sw_case(1e-20, SwKind::Kept), sw_case(1e-20, SwKind::Kept)],
            vec![
                sw_case(1e-20, SwKind::Kept),
                sw_case(1e-20, SwKind::Failed),
                sw_case(1e-20, SwKind::Kept),
            ],
        ];
        for limits in [one_job, one_window] {
            let mut heaps = heaps_of(2, 1);
            let (result, _) = sw_drive(
                &mut lists.clone(),
                &cases,
                &mut heaps,
                &[0..2, 0..2],
                limits,
            );
            let error = result.expect_err("query 1's failed redo fails the batch");
            assert_eq!(format!("{error:#}"), "batch 0 failed, first query holds 2");
            let mut model = heaps_of(2, 1);
            let reference = sw_one_thread(&lists, &cases, &mut model, &[0..2, 0..2]);
            assert_eq!(error_text(&reference), Some(format!("{error:#}")));
        }
    }

    /// The windows of the Smith-Waterman redo give the one-thread loop's heaps and its error (the
    /// failed batch, the endless redo, or none) on random searches: queries in batches, matches
    /// whose first-round tests are certain, deferred and turned, failed and endless redos, early
    /// termination, `-max_target_seqs` 1 to 500, windows of one job up to a few, small caps on
    /// jobs at risk.
    #[test]
    fn smith_waterman_windows_give_the_one_thread_heaps_and_error_on_random_searches() {
        let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
        let prelim = [1e-30, 1e-5, 0.001, 0.005, 0.02, 1.0];
        let first_round = [1e-30, 1e-5, ECUTOFF, 0.004, 0.01, 0.5, 2.0];
        let (mut turned, mut continued, mut ended) = (0, 0, 0);
        let (mut batch_failures, mut endless, mut clean) = (0, 0, 0);
        for run in 0..1500 {
            let queries = 1 + rng.below(8);
            let mut batches = Vec::new();
            while batches.len() < queries {
                let start = batches.len();
                let end = (start + 1 + rng.below(3)).min(queries);
                batches.extend((start..end).map(|_| start..end));
            }
            let heap_threshold = [1, 2, 3, 5, 500][rng.below(5)];
            let mut lists = Vec::new();
            let mut cases = Vec::new();
            for _ in 0..queries {
                if rng.below(6) == 0 {
                    lists.push(None);
                    cases.push(Vec::new());
                    continue;
                }
                let count = rng.below(25);
                let mut evalues: Vec<f64> = (0..count)
                    .map(|_| prelim[rng.below(prelim.len())])
                    .collect();
                evalues.sort_by(f64::total_cmp);
                lists.push(hit_list(
                    evalues
                        .iter()
                        .enumerate()
                        .map(|(oid, &evalue)| list(oid as u32, evalue, rng.below(10) == 0))
                        .collect(),
                ));
                cases.push(
                    (0..count)
                        .map(|_| SwCase {
                            evalue: first_round[rng.below(first_round.len())],
                            score: 30 + rng.below(4) as i32,
                            kind: match rng.below(100) {
                                0..=9 => SwKind::NotAsked,
                                10..=19 => SwKind::Purged,
                                20..=29 => SwKind::Failed,
                                30 if run % 4 == 0 => SwKind::Endless,
                                _ => SwKind::Kept,
                            },
                        })
                        .collect(),
                );
            }
            let limits = XWindowLimits {
                jobs: 1 + rng.below(6),
                at_risk_jobs: 1 + rng.below(3),
            };
            let mut heaps = heaps_of(queries, heap_threshold);
            let (result, _) = sw_drive(&mut lists.clone(), &cases, &mut heaps, &batches, limits);
            let mut model = heaps_of(queries, heap_threshold);
            let reference = sw_one_thread(&lists, &cases, &mut model, &batches);
            let at = format!("run {run}: {limits:?} threshold {heap_threshold} {batches:?}");
            if let Ok(counts) = &result {
                turned += counts.turned.expect("the Smith-Waterman counts");
                continued += counts.continued;
                ended += counts.skipped + counts.wasted;
            }
            let result = result.map(|_| ());
            assert_eq!(error_text(&result), error_text(&reference), "{at}");
            assert_eq!(heap_contents(&heaps), heap_contents(&model), "{at}");
            match error_text(&reference) {
                None => clean += 1,
                Some(text) if text.starts_with("endless") => endless += 1,
                Some(_) => batch_failures += 1,
            }
        }
        assert!(
            turned > 100
                && continued > 100
                && ended > 100
                && batch_failures > 50
                && endless > 20
                && clean > 300,
            "turned {turned} continued {continued} ended {ended} failed batches {batch_failures} \
             endless {endless} clean {clean}"
        );
    }

    /// The Smith-Waterman jobs on real matches (`x_redo_match_sw`): at 1, 2, 3 and 7 threads, with
    /// and without the jobs largest first, each computed match equals the one-thread loop's call
    /// (`postprocess_preliminary_hits_smith_waterman` with fresh scratch) whose heap test the
    /// query's heap at the window's start answers, and the job defers exactly the tests that this
    /// heap cannot answer for certain.
    #[test]
    fn smith_waterman_jobs_equal_the_one_thread_call_on_the_window_start_heaps() {
        let fixture = redo_fixture();
        let count = fixture.lists.len();
        // Two sets of heaps that keep two lists. The first: by query, empty, one list, and full
        // above the cutoff (worst E 0.5). The second: full within the cutoff (worst E 1e-25), but
        // for each batch's last query, which stays empty so that no batch ends early: there every
        // test of an E-value above the cutoff is FALSE for certain (the fixture asks such tests of
        // the first query of each batch).
        let heap_sets: Vec<Vec<BlastCompoHeap>> = [false, true]
            .into_iter()
            .map(|within| {
                let mut heaps = heaps_of(count, 2);
                for (q_idx, heap) in heaps.iter_mut().enumerate() {
                    let kept: &[(f64, i32)] = if !within {
                        &[(1e-30, 300), (0.5, 60)][..q_idx % 3]
                    } else if q_idx + 1 == fixture.batches[q_idx].end {
                        &[]
                    } else {
                        &[(1e-30, 300), (1e-25, 250)]
                    };
                    for (index, &(evalue, score)) in kept.iter().enumerate() {
                        let oid = 1000 + index as u32;
                        heap.insert(list(oid, evalue, false), evalue, score, oid);
                    }
                }
                heaps
            })
            .collect();
        let build = |q_idx: usize| -> Result<BlastRedoAlignParams> {
            Ok(BlastRedoAlignParams {
                smith_waterman: true,
                ..redo_params(&fixture, q_idx)
            })
        };
        let (mut compared, mut with_list, mut deferred_tests, mut certain_false) = (0, 0, 0, 0);
        let mut continuations = 0;
        let stream_end = XMatch {
            q_idx: count,
            position: 0,
        };
        for heaps in &heap_sets {
            let all = stream(&fixture.lists);
            let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
            let caps = vec![usize::MAX; count];
            let windows: Vec<XWindow> = (0..8)
                .map(|index| {
                    let start = all[rng.below(all.len())];
                    let mut window = x_window::<true>(
                        &fixture.lists,
                        heaps,
                        &fixture.batches,
                        &vec![None; count],
                        &caps,
                        None,
                        start,
                        2 + rng.below(10),
                    );
                    for ended in &mut window.ended {
                        *ended = index % 2 == 1 && rng.below(4) == 0;
                    }
                    window
                })
                .collect();
            // Windows after a replay that stopped at a match of each window above whose test is
            // asked (built at the first thread count, with the test's E-value and score).
            let mut carried_windows: Vec<XWindow> = Vec::new();
            for threads in [1, 2, 3, 7] {
                crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                    let search = XRedoParSearch {
                        pool,
                        threads,
                        contexts: &fixture.contexts,
                        query_nomask: query_raw,
                        search_spaces: &fixture.search_spaces,
                        subjects: crate::algorithm::blastp::x_redo_subjects::XRedoSubjects::Encoded(
                            &fixture.subjects,
                        ),
                        query_batches: &fixture.batches,
                        evalue_adjust_contexts: &fixture.evalue_adjust_contexts,
                        matrix: ScoringMatrix::Blosum62,
                        gapped_params: &fixture.gapped_params,
                        gapped_gumbel: Some(&fixture.gapped_gumbel),
                        evalue: 10.0,
                        build_redo_align_params: &build,
                        timing: false,
                        engine_positives: false,
                        sw_errors: TEST_SW_ERRORS,
                    };
                    // The one-thread loop's call with fresh scratch, its heap test answered by
                    // `answer`; and the test's E-value and score, if it was asked.
                    let one_thread = |matched: XMatch, answer: &dyn Fn(f64, i32) -> bool| {
                        let q_idx = matched.q_idx;
                        let local_match = x_local_match(&fixture.lists, matched);
                        let hits: Vec<BlastpPreliminaryHsp> = local_match
                            .hsps
                            .iter()
                            .map(preliminary_hit_from_local_match_hsp)
                            .collect();
                        let mut scratch = XScratch::new(
                            BlastpKappaSubjectRangeCache::new(None),
                            ScoringMatrix::Blosum62,
                        );
                        let mut sw_scratch = SwMatchScratch::new();
                        let subject = search.subjects.read(local_match.oid).expect("a subject");
                        let mut asked = None;
                        let result = postprocess_preliminary_hits_smith_waterman(
                            &hits,
                            &x_query_workspace(&search, q_idx),
                            &mut scratch.composition_workspace,
                            query_raw(&fixture.contexts[q_idx]),
                            subject.residues(),
                            ScoringMatrix::Blosum62,
                            &fixture.gapped_params,
                            Some(&fixture.gapped_gumbel),
                            &build(q_idx).expect("parameters"),
                            10.0,
                            &mut scratch.gap_scratch,
                            &mut scratch.range_cache,
                            false,
                            &mut sw_scratch,
                            &mut |evalue, score| {
                                asked = Some((evalue, score));
                                answer(evalue, score)
                            },
                        );
                        (result.map_err(|error| format!("{error:#}")), asked)
                    };
                    if carried_windows.is_empty() {
                        for window in &windows {
                            let Some(index) = (0..window.matches.len()).find(|&index| {
                                !window.ended[index]
                                    && one_thread(window.matches[index], &|_: f64, _: i32| false)
                                        .1
                                        .is_some()
                            }) else {
                                continue;
                            };
                            let matched = window.matches[index];
                            let asked = one_thread(matched, &|_: f64, _: i32| false)
                                .1
                                .expect("asked");
                            let stopped = XComputed {
                                result: no_alignment(),
                                deferred: Some(asked),
                                first_round_s: 0.0,
                            };
                            let carry = XCarry {
                                matches: window.matches[index..].to_vec(),
                                results: std::iter::once(Some(stopped))
                                    .chain(window.ended[index + 1..].iter().map(|&ended| {
                                        (!ended).then(|| XComputed::from(no_alignment()))
                                    }))
                                    .collect(),
                            };
                            carried_windows.push(x_window::<true>(
                                &fixture.lists,
                                heaps,
                                &fixture.batches,
                                &vec![None; count],
                                &caps,
                                Some(&carry),
                                window.next.unwrap_or(stream_end),
                                window.jobs.len() + 2,
                            ));
                        }
                    }
                    for window in windows.iter().chain(&carried_windows) {
                        for lpt in [false, true] {
                            let shared = Arc::new(XSharedSubjectRanges::new());
                            let (computed, _) = x_compute_window::<_, true>(
                                &search,
                                &fixture.lists,
                                heaps,
                                window,
                                &shared,
                                lpt,
                            );
                            let in_jobs: Vec<usize> =
                                window.jobs.iter().flat_map(|job| job.clone()).collect();
                            assert_eq!(computed.len(), in_jobs.len());
                            assert_eq!(
                                in_jobs.len() + window.carried,
                                window.matches.len() + window.continued.len(),
                                "the carried matches but the continued ones are in no job"
                            );
                            for (index, computed) in in_jobs.into_iter().zip(&computed) {
                                let matched = window.matches[index];
                                let at = (threads, lpt, matched);
                                let Some(computed) = computed else {
                                    assert!(window.ended[index], "{at:?}");
                                    continue;
                                };
                                let q_idx = matched.q_idx;
                                let oid = x_local_match(&fixture.lists, matched).oid;
                                let heap = &heaps[q_idx];
                                let before = index
                                    - window.matches.partition_point(|other| other.q_idx < q_idx);
                                // A test the job cannot answer for certain on the window's start.
                                let deferrable = |evalue: f64| {
                                    evalue > heap.ecutoff && heap.n + before >= heap.heap_threshold
                                };
                                let continued = window
                                    .continued
                                    .iter()
                                    .find(|&&(at, _)| at == index)
                                    .map(|&(_, test)| test);
                                // The one-thread call answered as the job answers: TRUE for a
                                // continuation, FALSE for a deferred test, else the heap at the
                                // window's start.
                                let (reference, asked) =
                                    one_thread(matched, &|evalue: f64, score: i32| {
                                        continued.is_some()
                                            || (!deferrable(evalue)
                                                && heap.clone().would_insert(evalue, score, oid))
                                    });
                                with_list += usize::from(
                                    reference
                                        .as_ref()
                                        .is_ok_and(|result| result.hsp_list.is_some()),
                                );
                                assert_eq!(
                                    x_result_difference(&computed.result, &reference),
                                    None,
                                    "{at:?}"
                                );
                                if let Some(test) = continued {
                                    assert_eq!(asked, Some(test), "{at:?}");
                                    assert_eq!(computed.deferred, None, "{at:?}");
                                    continuations += 1;
                                    continue;
                                }
                                let expected = asked.filter(|&(evalue, score)| {
                                    deferrable(evalue) && heap.would_insert_view(evalue, score, oid)
                                });
                                assert_eq!(computed.deferred, expected, "{at:?}");
                                deferred_tests += usize::from(expected.is_some());
                                certain_false +=
                                    usize::from(asked.is_some_and(|(evalue, score)| {
                                        !heap.would_insert_view(evalue, score, oid)
                                    }));
                                compared += 1;
                            }
                        }
                    }
                    Ok(())
                })
                .expect("a pool");
            }
        }
        assert!(
            compared > 150
                && with_list > 30
                && deferred_tests > 10
                && certain_false > 10
                && continuations > 50,
            "compared {compared} with a list {with_list} deferred {deferred_tests} certain FALSE \
             {certain_false} continuations {continuations}"
        );
    }
}
