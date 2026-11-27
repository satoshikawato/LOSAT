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
//!    a query, so it holds at most 4 x 64 x threads results before their replay (two budgets with a
//!    `-use_sw_tback` carry, below). A window boundary
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
//!   logarithmic in q's work at risk. With `-use_sw_tback` a replay may stop inside its window
//!   and carry the rest, with their results, to the next one (below). q's results are still
//!   dropped in one window, the one whose replay reaches them after the heaps fill: a replay never
//!   stops within q after dropping a result of q (a stop is at a match the test did not end, and
//!   once the test ends a match of q it ends every later one: their E-values are not smaller); a
//!   carried match was at risk or not by the heaps at its own window's start, which stays valid.
//!   That window's results of q at risk are those of its new jobs (at most cap(q)) and those
//!   carried from the windows before it whose replay stopped before passing all of q's matches,
//!   each of which took at most cap(q) new jobs at risk of q with the same cap(q): cap(q) grows
//!   only by q's jobs at risk of a window whose replay passed all of q's matches in it. So it
//!   drops at most 4 x cap(q) x (1 + those stopped windows) results of q, and at most the results
//!   it holds with the carried ones (2 x 4 x 64 x threads, below).
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
//! window, which finishes the match (a continuation: the redo resumed at the test answered TRUE,
//! as the replay answered it, from the state that the job kept there, `SwMatchResume`, so it
//! repeats neither the composition adjustment nor the first round's score-only scan) and, in the
//! same way, each carried match whose deferred test is TRUE for certain on the
//! heaps where the replay stopped (`x_sw_continued`: the carried matches before it that keep no
//! list do not count), while its other jobs compute the matches after the carried ones. The carried
//! results stay valid: a certain answer stays certain, and a deferred one is asked at its replay.
//! A window with a carry takes at most a budget of new jobs, as one without, and the carried
//! matches, in jobs of four, and the new jobs together are at most two budgets, so a window holds
//! at most 2 x 4 x 64 x threads results before their replay with the carried ones, and a carry
//! fewer than that; the serial work of a stop (`x_sw_continued`, the carry's copy) is linear in
//! the carry (it reads no heap's root), at most that of a window's replay. The replay also keeps
//! the batch's status of NCBI's one-thread redo (`XSwBatch`, with core's `SwBatchStatus`, the
//! one-thread loop's too): a match that the redo reaches sets it, failed or not; a match the
//! early termination test ends leaves it; a list whose subject the redo's fetch excludes (line RG
//! X4-excl) sets it to 0 unless the early termination test ends that list, read when the replay
//! passes the list's query, with the one-thread loop's reader; a batch whose status says failed at
//! its end fails the search there, and an endless redo fails it at that match. A continuation is
//! replayed before every match after it, so a batch ends after its continuations. A Smith-Waterman
//! job is at risk only as an X-drop job is (the early termination test may end its match). With
//! `LOSAT_TIMING=1` the `[TIMING] x_bpredopar:` line of such a search ends with `turned=<t>
//! continued=<c> cont_stop=<a> cont_ahead=<b> cont_bound=<d> carried_max=<m>`: the deferred
//! answers that the replay answered FALSE, and the carried matches continued, by the reason why
//! their test is TRUE (the stopped match itself, the exact look-ahead in a later query, the bound
//! on the lists the heap holds); the most matches a replay carried. A job keeps the state of each
//! test it defers (about 3.2 KB, its matrix) with the match's result until the replay answers the
//! test: one state per match that a window holds with its carry at most (2 x 4 x 64 x threads),
//! freed by a FALSE at the replay or by the continuation's job.
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
//! queries of about one window instead of all queries. With `-use_sw_tback` they are also those of
//! the queries of the matches carried behind the window's start, fewer than two budgets' results;
//! the lists read ahead stay one budget's (a window whose carry fills its second budget takes no
//! new job), so the lists held span about three budgets: two with the carry and one read ahead,
//! with the lists of the first carried query and of the last query read ahead whole.
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
use crate::core::composition_adjustment::smith_waterman_match::{SwBatchStatus, SwRedoFailed};
use crate::stats::search_space::SearchSpace;
use crate::stats::{BlastGumbelBlk, KarlinParams};
use crate::utils::threading::SearchPool;

use super::blast_engine::{
    preliminary_hit_from_local_match_hsp, query_fetch_excluded_reached, sw_batch_failed_after,
    sw_endless_after, FetchExcludedMatch, FetchNotes,
};
use super::gapalign::{BlastpPreliminaryHsp, GapAlignScratch};
use super::hsp::{
    blast_compo_early_termination, BlastCompoHeap, BlastpHitList, BlastpHspList,
    COMPO_HEAP_EVALUE_STRETCH,
};
#[cfg(test)]
use super::hsp::{fill_results_from_compo_heaps, BlastpHsp};
use super::kappa::smith_waterman::{SwHeapAnswer, SwMatchResume, SwMatchScratch};
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
pub(super) fn x_window_jobs_per_thread() -> usize {
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
pub(super) const X_JOB_MATCHES: usize = 4;

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
    /// The lists of the Smith-Waterman redo's stream whose subjects the redo's fetch excludes (line
    /// RG X4-excl), by query and place: they left the preliminary hit lists before the redo, each
    /// leaving a mark of its place (`blast_engine::database::exclude_fetched`), and change the
    /// status of their query batch (`XSwBatch`); empty without a database or a negative seqid list
    /// with ids.
    pub sw_excluded: &'a [FetchExcludedMatch],
    /// NCBI's application name (`"blastp"` or `"psiblast"`) for Smith-Waterman batch errors.
    pub program: &'static str,
    /// Whether `-sum_stats true` (`do_link_hsps`) is enabled for Smith-Waterman endless errors.
    pub do_link_hsps: bool,
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
    /// The Smith-Waterman redo: the carried matches finished by a continuation (`XCarry`): the
    /// heap tests that a job deferred and that are TRUE at their replay.
    continued: usize,
    /// `continued` by the reason why the test is TRUE (`x_sw_continued`).
    continued_by: XContinuedBy,
    /// The most matches that a replay carried to the next window.
    carried_max: usize,
}

/// The carried matches that a window continues (`x_sw_continued`), by the reason why their heap
/// test is TRUE: the stopped match, whose test the replay answered; the exact look-ahead in a
/// query after the stopped match's; the bound on the lists that the heap holds.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct XContinuedBy {
    stop: usize,
    ahead: usize,
    bound: usize,
}

impl XContinuedBy {
    fn add(&mut self, other: Self) {
        self.stop += other.stop;
        self.ahead += other.ahead;
        self.bound += other.bound;
    }
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
    /// With `SW`, after a replay that stopped at a deferred heap test answered TRUE (`XCarry`): the
    /// carried matches whose test is TRUE (`x_sw_continued`), each the match of a job of its own (a
    /// continuation), with the test's E-value and score; by match index.
    continued: Vec<(usize, (f64, i32))>,
    /// `continued` by the reason why the test is TRUE.
    continued_by: XContinuedBy,
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
///   matches since that may insert one: those with alignments, those continued, and those with a
///   deferred test. A match ended before its window inserts none.
///
/// A deferred test that the rule does not continue counts as one that may insert a list, without
/// asking the heap: below the threshold the test without its side effect is TRUE whatever the key,
/// and at or past it the rule continues no later match of the query, whatever the count. So the
/// rule reads only the heap's count and threshold, never its root (in the array form a scan of its
/// lists, which a stop would otherwise make once per carried match; review 2 of X4-mt, N-5), and
/// the serial work is linear in the carry.
///
/// A match that the early termination test may end is not continued (its continuation could be
/// computed in vain); nor is any match after a redo that would end the search.
fn x_sw_continued(
    carry: &XCarry,
    lists: &[Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    end_above: &[Option<f64>],
) -> (Vec<(usize, (f64, i32))>, XContinuedBy) {
    let mut continued = Vec::new();
    let mut by = XContinuedBy::default();
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
            by.stop += 1;
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
                        by.ahead += 1;
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
                by.bound += 1;
                may_insert += 1;
            }
            // It may insert a list (the doc above: no heap test needed).
            Some(_) => may_insert += 1,
            None => may_insert += usize::from(inserts),
        }
    }
    (continued, by)
}

/// A match computed by its job: the result of its redo and, for a Smith-Waterman match whose heap
/// test the job deferred (`x_sw_heap_test`: the redo stopped at the test as on FALSE), that test's
/// E-value and score, which the replay asks the query's heap with, and the state of the redo at
/// the test (`resume`), from which a continuation finishes the match where the replay answers
/// TRUE. A FALSE at the replay frees the state with the result; a TRUE carries it to the
/// continuation's job.
#[derive(Debug)]
struct XComputed {
    result: Result<BlastpPostprocessResult>,
    deferred: Option<(f64, i32)>,
    resume: Option<Box<SwMatchResume>>,
}

impl From<Result<BlastpPostprocessResult>> for XComputed {
    fn from(result: Result<BlastpPostprocessResult>) -> Self {
        Self {
            result,
            deferred: None,
            resume: None,
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
/// the jobs from there while it holds fewer than `jobs` jobs (new jobs, with the carried matches
/// within two budgets, where a `SW` window has a carry: below) and, for a job at risk of query q,
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
///
/// With `SW` (the Smith-Waterman redo) and `carry` (the replay before stopped at a deferred heap
/// test answered TRUE), the window starts with the carried matches: those whose test is TRUE for
/// certain (`x_sw_continued`) are each a job of its own (a continuation), the others are in no job
/// (their results come with the carry); its other jobs start at `start`. As without a carry, the
/// window takes at most `jobs` new jobs; the carried matches, in jobs of `X_JOB_MATCHES`, and the
/// new jobs together are at most `2 x jobs` (so a window holds at most `2 x X_JOB_MATCHES x jobs`
/// matches with the carried ones, and a carry fewer: one budget of new jobs, and the carried
/// matches in place of new jobs past one budget). The carried matches count as the window's
/// matches of their query (each inserts at most one list); the last carried query's new jobs at
/// risk start from none, so a query's results at risk held before their replay grow by at most
/// `caps[q]` jobs per stopped window (`x_redo_windows` does not grow `caps[q]` there). A deferred
/// heap test that turns FALSE wastes nothing (`x_sw_heap_test`), so a Smith-Waterman job is at
/// risk only as an X-drop job is.
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
        continued_by: XContinuedBy::default(),
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
    // Whether the early termination test may end a match of query `q_idx` with `best_evalue`
    // at its replay, the window holding `before` of the query's matches before it.
    let may_be_ended = |q_idx: usize, batch_open: bool, best_evalue: f64, before: usize| {
        batch_open
            && end_above[q_idx].is_some_and(|above| best_evalue > above)
            && heaps[q_idx].may_fill(before)
    };
    // With `SW`: the carried matches, in jobs of `X_JOB_MATCHES`, which count toward the second
    // budget.
    let mut carried_jobs = 0;
    if let (true, Some(carry)) = (SW, carry) {
        first_query = carry.matches[0].q_idx;
        window.matches.extend_from_slice(&carry.matches);
        window
            .ended
            .extend(carry.results.iter().map(|computed| computed.is_none()));
        window.carried = carry.matches.len();
        carried_jobs = window.carried.div_ceil(X_JOB_MATCHES);
        (window.continued, window.continued_by) = x_sw_continued(carry, lists, heaps, end_above);
        window
            .jobs
            .extend(window.continued.iter().map(|&(index, _)| index..index + 1));
        let last = window.matches[window.matches.len() - 1].q_idx;
        let last_from = window
            .matches
            .partition_point(|matched| matched.q_idx < last);
        query = Some((
            last,
            batch_open(last, first_query),
            window.matches.len() - last_from,
            0,
        ));
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
        let mark = window.matches.len();
        let mut at_risk = false;
        for (index, &matched) in job.iter().enumerate() {
            let best_evalue = x_local_match(lists, matched).best_evalue;
            let may_end = may_be_ended(q_idx, batch_open, best_evalue, in_window + index);
            let ended = may_end
                && blast_compo_early_termination(best_evalue, &heaps[query_batches[q_idx].clone()]);
            at_risk |= may_end && !ended;
            window.ended.push(ended);
        }
        // The budget: `jobs` new jobs; with `SW` also at most `2 x jobs` with the carried matches'
        // (the continuations are carried matches, not new jobs).
        let full = if SW {
            let stream = window.stream_jobs();
            stream >= jobs || carried_jobs + stream >= 2 * jobs
        } else {
            window.jobs.len() >= jobs
        };
        // A window holds at least one job.
        if !window.jobs.is_empty() && (full || (at_risk && at_risk_jobs >= caps[q_idx])) {
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
/// FALSE; on TRUE a continuation in the next window finishes it (`XCarry`), resumed from the
/// state that the job keeps at the test (`x_redo_match_sw`).
///
/// Certain TRUE: the E-value is within the heap's cutoff (the test's second condition), or `n +
/// before` is below the threshold (the first condition at the replay: only this query's matches
/// insert into its heap, each at most one list, and an insert at or past the threshold leaves `n`).
///
/// Certain FALSE: the test without its side effect says FALSE on the window's start, so `n` is at
/// least `heapThreshold`, `e` is above `ecutoff` and at least `worstEvalue`, and the candidate is
/// not better than the root (`compo_heap_record_compare`, E-value, then score, then subject index).
/// A later insert keeps `n` at or past the threshold. It either adds a list while every kept
/// E-value and the new one are within the cutoff (so the root's E-value, the new worst, is within
/// it and below `e`), or replaces the root by a better record and restores the heap (the new root
/// is not worse than the old one, so the worst E-value does not grow and the candidate is not
/// better than the new root), or changes nothing. Each case keeps all four conditions, so the test
/// stays FALSE up to the replay (`hsp_tests.rs`,
/// `test_compo_heap_would_insert_view_equals_would_insert_and_a_false_stays_false`).
pub(super) fn x_sw_heap_test(
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
    /// A continuation (`XCarry`): TRUE, as the replay answered the test that the match's job left
    /// open; the redo resumes past the test (`x_redo_match_sw`), so it is not asked again.
    Continued,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3602-3621
/// The Smith-Waterman redo of one match (`-use_sw_tback`): the call of the one-thread loop, with its
/// heap test answered as `test` says. A test that the job defers keeps the redo's state there
/// (`XComputed::resume`); a continuation resumes from that state (`resume`, the state that the
/// match's job kept), past its test, without the composition adjustment and the first round's
/// score-only scan that the job made.
#[allow(clippy::too_many_arguments)]
fn x_redo_match_sw<P>(
    search: &XRedoParSearch<'_, P>,
    q_idx: usize,
    local_match: &BlastpHspList,
    query_workspace: &BlastpKappaQueryWorkspace,
    redo_align_params: &BlastRedoAlignParams,
    scratch: &mut XScratch,
    test: XSwTest<'_>,
    resume: Option<&SwMatchResume>,
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
                resume: None,
            }
        }
    };
    let subject_raw = subject.residues();
    let mut deferred = None;
    let mut kept = None;
    let resume = match test {
        XSwTest::Heap(..) => None,
        XSwTest::Continued => {
            Some(resume.expect("a continuation resumes from the state that its job kept"))
        }
    };
    let result = postprocess_preliminary_hits_smith_waterman(
        &scratch.hits,
        query_workspace,
        &mut scratch.composition_workspace,
        (search.query_nomask)(&search.contexts[q_idx]),
        subject_raw,
        search.matrix,
        None,
        search.gapped_params,
        search.gapped_gumbel,
        redo_align_params,
        search.evalue,
        &mut scratch.gap_scratch,
        &mut scratch.range_cache,
        search.engine_positives,
        scratch.sw.get_or_insert_with(SwMatchScratch::new),
        resume,
        &mut |evalue, score| match test {
            XSwTest::Heap(heap, before) => {
                let inserts = x_sw_heap_test(
                    heap,
                    before,
                    (evalue, score, local_match.oid),
                    &mut deferred,
                );
                if deferred.is_some() {
                    SwHeapAnswer::Deferred
                } else {
                    SwHeapAnswer::Known(inserts)
                }
            }
            // BLASTP's match has one window, whose one test the continuation resumes past.
            XSwTest::Continued => unreachable!("a continuation's heap test is asked again"),
        },
        &mut kept,
    );
    debug_assert_eq!(
        deferred.is_some(),
        kept.is_some(),
        "a deferred test keeps the redo's state"
    );
    XComputed {
        result,
        deferred,
        resume: kept,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3447-3449,3493-3511
/// The redo of a window's matches on the pool, in jobs of up to `X_JOB_MATCHES` consecutive
/// matches of one query; the results of the jobs' matches in order (all of the window's but the
/// carried ones), `None` for a match ended before the window (not computed). `heaps` are the heaps as the window's start left them; nothing changes them
/// while the jobs run. NCBI's loop also runs the early termination test and the heap insertion of
/// each match on its thread; here they are replayed afterwards (`x_replay`). With `lpt`
/// (`LOSAT_X_BPREDOLPT`), more than one thread on the search's pool and two or more jobs, the
/// jobs run as `x_compute_window_lpt` runs them; the second value says whether they did. With `SW`
/// each match runs the Smith-Waterman redo (`x_redo_match_sw`), its heap test answered from
/// `heaps`, or TRUE for a continuation (a match in `window.continued`), which resumes from the
/// state that its job kept (in `carried`, the carry's results: the window's carried matches').
fn x_compute_window<P, const SW: bool>(
    search: &XRedoParSearch<'_, P>,
    lists: &[Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    window: &XWindow,
    carried: &[Option<XComputed>],
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
                    .map(|&ended| {
                        (!ended).then(|| XComputed {
                            result: Err(anyhow::anyhow!("{error:#}")),
                            deferred: None,
                            resume: None,
                        })
                    })
                    .collect();
            }
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529
        // No NCBI counterpart: the Newton lanes solve the first Newton problems of the job's
        // matches, four per vector (kappa::x_newton_lanes_prefetch_matches); a result is used
        // only by a call on this thread with the same input bits. It does not change any value
        // NCBI computes. Given the heaps of the batch as the window's start left them, as the
        // one-thread loop gives its heaps, it leaves out the matches the test ends on them: those
        // not computed here. A continuation's job (a carried match) resumes past the composition
        // adjustment, which reads none.
        if !SW || job.start >= window.carried {
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
                        let at = job.start + index;
                        let (test, resume) =
                            match window.continued.binary_search_by_key(&at, |&(at, _)| at) {
                                Ok(found) => {
                                    let stopped = carried[at].as_ref();
                                    debug_assert_eq!(
                                        stopped.and_then(|stopped| stopped.deferred),
                                        Some(window.continued[found].1),
                                        "a continuation resumes from the test that was continued"
                                    );
                                    (
                                        XSwTest::Continued,
                                        stopped.and_then(|stopped| stopped.resume.as_deref()),
                                    )
                                }
                                Err(_) => (XSwTest::Heap(&heaps[q_idx], before_job + index), None),
                            };
                        x_redo_match_sw(
                            search,
                            q_idx,
                            x_local_match(lists, matched),
                            query_workspace,
                            &redo_align_params,
                            scratch,
                            test,
                            resume,
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
pub(super) fn x_lpt_order(costs: &[u64]) -> Vec<usize> {
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
pub(super) fn x_compute_window_lpt<S, R>(
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3644-3646,3735-3737,3774-3778
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3540
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3560-3575
// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:250-252
/// The status of a query batch in NCBI's one-thread Smith-Waterman redo (`status_code_tld[0]`),
/// replayed in stream order with the status of every program's redo loop, the one-thread loop's
/// too (core's `smith_waterman_match::SwBatchStatus`, the one owner of its rule): a match that the redo reaches sets it (`redone`: failed, or 0 for a redo
/// that ends without failing, with alignments or without, as when its heap test is FALSE); a match
/// the early termination test ends leaves it; a list whose subject the redo's fetch excludes sets it
/// to 0 unless the early termination test ends that list (line RG X4-excl,
/// `XRedoParSearch::sw_excluded`); at the batch's end a failed status fails the search there
/// (`redo_run::sw_batch_failed`), and the next batch starts at 0. It is kept across windows: a
/// window may end inside a batch. Both errors carry the warning of a header that the batch's
/// fetches read before them, or give way to LOSAT's refusal where such a fetch read a header that
/// cannot be read, with the one-thread loop's own functions (`redo_run::sw_batch_failed_after`,
/// `redo_run::sw_endless_after`, on the heaps as the replay left them, which are the loop's
/// there).
struct XSwBatch<'a> {
    /// The marks of the excluded lists of the queries that the replay has not passed yet.
    excluded: &'a [FetchExcludedMatch],
    /// The notes of the subjects that the search's fetches read.
    notes: &'a FetchNotes<'a>,
    /// NCBI's application name (`"blastp"` or `"psiblast"`).
    program: &'static str,
    /// Whether `-sum_stats true` (`do_link_hsps`) is enabled.
    do_link_hsps: bool,
    /// The query of the last match replayed, if any.
    query: Option<usize>,
    /// The place of that query's last match whose redo set the status, if any.
    status_position: Option<usize>,
    /// The status of that query's batch.
    status: SwBatchStatus,
}

impl<'a> XSwBatch<'a> {
    fn new(
        excluded: &'a [FetchExcludedMatch],
        notes: &'a FetchNotes<'a>,
        program: &'static str,
        do_link_hsps: bool,
    ) -> Self {
        Self {
            excluded,
            notes,
            program,
            do_link_hsps,
            query: None,
            status_position: None,
            status: SwBatchStatus::default(),
        }
    }

    /// Before the replay of a match of query `next`, or after the stream's last match (`None`),
    /// where `next` is in a later query than the last match replayed: the steps of the one-thread
    /// loop between them. First the excluded lists of the queries passed, from the last match's
    /// query up to `next` or the end of that query's batch, each query's read with the one-thread
    /// loop's reader after its loop (`redo_run::query_fetch_excluded_reached`) on the heaps as they
    /// are now, which are those after the query's loop (the queries after the last match's have no
    /// match). Then, if `next` is in a later batch or there is none, the end of the last match's
    /// batch: a failed status fails the search there, and otherwise the next batch starts at 0 (the
    /// batches between them have no match the redo reaches, so their status stays 0). The excluded
    /// lists of a batch's queries before its first match replayed change nothing (the status is 0
    /// there, and such a list only sets it to 0), nor do those of a batch with no match replayed.
    fn pass(
        &mut self,
        next: Option<usize>,
        query_batches: &[Range<usize>],
        heaps: &[BlastCompoHeap],
    ) -> Result<()> {
        if next.is_some() && next == self.query {
            return Ok(());
        }
        let to = next.unwrap_or(query_batches.len());
        let passed = self.excluded.partition_point(|mark| mark.query < to);
        let (passed, later) = self.excluded.split_at(passed);
        self.excluded = later;
        if let Some(query) = self.query {
            let batch = query_batches[query].clone();
            for marks in passed.chunk_by(|a, b| a.query == b.query) {
                let q_idx = marks[0].query;
                if (query..batch.end).contains(&q_idx) {
                    let status_position = if q_idx == query {
                        self.status_position
                    } else {
                        None
                    };
                    query_fetch_excluded_reached(
                        marks,
                        status_position,
                        &heaps[batch.clone()],
                        &mut self.status,
                    );
                }
            }
            if next.is_none_or(|next| !batch.contains(&next)) {
                if self.status.failed() {
                    return Err(sw_batch_failed_after(
                        self.program,
                        self.notes,
                        batch.clone(),
                        heaps,
                        batch.end - 1,
                    ));
                }
                self.status = SwBatchStatus::default();
            }
        }
        self.query = next;
        self.status_position = None;
        Ok(())
    }

    /// A match `matched` of the query of the last `pass` that the redo reaches, with its redo's
    /// `result`: it sets the batch's status (`SwBatchStatus::redone`). The match's result where
    /// its redo ended, none where it failed, the error of a redo that does not end in the
    /// match's batch (`redo_run::sw_endless`, on `heaps` as the replay left them before the
    /// match), any other error as it is.
    fn redone(
        &mut self,
        matched: XMatch,
        query_batches: &[Range<usize>],
        heaps: &[BlastCompoHeap],
        result: Result<BlastpPostprocessResult>,
    ) -> Result<Option<BlastpPostprocessResult>> {
        self.status_position = Some(matched.position);
        let notes = self.notes;
        let program = self.program;
        let do_link_hsps = self.do_link_hsps;
        let batch = query_batches[matched.q_idx].clone();
        self.status.redone(result, || {
            sw_endless_after(
                program,
                do_link_hsps,
                notes,
                batch.clone(),
                heaps,
                matched.q_idx,
                matched.position,
            )
        })
    }
}

/// The steps of the one-thread loop that depend on order, for a window's matches in stream order,
/// with the results computed for them: the early termination test on the heaps of the query's
/// batch, the match's error, and the heap insertion. With `SW` (the Smith-Waterman redo) also each
/// deferred heap test on the query's heap (a FALSE replays the match with the job's result, a redo
/// that stopped at the test; a TRUE stops the replay there and leaves the match and the matches
/// after it to the next window in `carry`), and the batch's status (`batch`): a failed redo sets
/// it, an endless one fails the search at that match, and the excluded lists of the queries that
/// the replay passes are read (`XSwBatch::pass`).
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
        if SW {
            batch.pass(Some(matched.q_idx), query_batches, heaps)?;
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
            resume,
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
        // match's redo goes on past the test: a continuation in the next window, resumed from the
        // job's state at the test, and the replay of the matches after it waits for it. A FALSE
        // frees that state.
        if let (true, Some(test)) = (SW, deferred) {
            if heaps[matched.q_idx].would_insert(test.0, test.1, local_match.oid) {
                let stopped = XComputed {
                    result,
                    deferred,
                    resume,
                };
                *carry = Some(XCarry {
                    matches: matches[index..].to_vec(),
                    results: std::iter::once(Some(stopped)).chain(results).collect(),
                });
                return Ok(replayed);
            }
            replayed.turned += 1;
        }
        let postprocessed = if SW {
            // A match the redo reaches sets the batch's status; a failed redo leaves it without
            // alignments.
            match batch.redone(matched, query_batches, heaps, result)? {
                Some(postprocessed) => postprocessed,
                None => continue,
            }
        } else {
            result?
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
/// The windows of a search, each computed (`compute`, with the lists and the heaps as the replay
/// of the windows before left them) and then replayed in stream order into `heaps`. A query's
/// preliminary list is released after the replay of the window that holds its last match. With
/// `expand`, the compact lists are expanded ahead of each window. With `SW` (the Smith-Waterman
/// redo) the replay keeps the query batches' status across the windows (`XSwBatch`, with the
/// excluded lists `sw_excluded` and the notes `notes` of the search's fetches), and the batch of
/// the last match ends after the last window.
#[allow(clippy::too_many_arguments)]
fn x_redo_windows<const SW: bool>(
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    query_batches: &[Range<usize>],
    limits: XWindowLimits,
    timing: bool,
    mut expand: Option<&mut XExpandAhead<'_>>,
    sw_excluded: &[FetchExcludedMatch],
    notes: &FetchNotes<'_>,
    program: &'static str,
    do_link_hsps: bool,
    compute: &mut dyn FnMut(
        &[Option<BlastpHitList>],
        &[BlastCompoHeap],
        &XWindow,
        &[Option<XComputed>],
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
    // at risk of a window whose replay passed all of q's matches in it (all used; it doubles when
    // it closed the window), up to the job budget. So it is at most `limits.at_risk_jobs` plus q's
    // jobs at risk computed and used before; and with `SW` a window whose replay stopped before
    // passing q's matches leaves it as it is, so each window that carries q's results at risk on
    // to the one that drops them took at most `caps[q]` new jobs at risk of q.
    let budget = limits.jobs.max(limits.at_risk_jobs);
    let mut caps = vec![limits.at_risk_jobs; lists.len()];
    let mut sw_batch = XSwBatch::new(sw_excluded, notes, program, do_link_hsps);
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
        // No NCBI counterpart: a check of the window's bounds (module header, `x_window`); it
        // changes no value NCBI computes. At most `jobs` new jobs; at most one budget of matches,
        // and with `SW` two with the carried ones; a carry fewer than two.
        debug_assert!(
            window.stream_jobs() <= limits.jobs
                && window.matches.len()
                    <= (X_JOB_MATCHES << usize::from(SW)).saturating_mul(limits.jobs)
                && window.carried < (2 * X_JOB_MATCHES).saturating_mul(limits.jobs),
            "a window over its bounds: {} new jobs, {} matches, {} carried (budget {} jobs)",
            window.stream_jobs(),
            window.matches.len(),
            window.carried,
            limits.jobs
        );
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
            // The carried results hold the states that the continuations resume from.
            let carried = carry.as_ref().map_or(&[][..], |carry| &carry.results[..]);
            let mut results = compute(lists, heaps, &window, carried);
            if let Some(par_start) = par_start {
                counts.par_s += par_start.elapsed().as_secs_f64();
            }
            // The carried results, the continued ones replaced by their continuations' (which frees
            // the states they resumed from), then the results of the jobs after them.
            if let Some(carry) = carry.take() {
                counts.carried_max = counts.carried_max.max(window.carried);
                counts.continued_by.add(window.continued_by);
                let mut computed = results.into_iter();
                results = carry.results;
                for index in window.jobs.iter().flat_map(|job| job.clone()) {
                    let result = computed.next().expect("one result per match of a job");
                    if index < window.carried {
                        results[index] = result;
                        counts.continued += 1;
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
            // With `SW`, a replay that stopped passed none of the matches from the stopped one's
            // query on (they go to the next window): the caps of those queries do not grow.
            let passed_to = match (SW, &carry) {
                (true, Some(carry)) => carry.matches[0].q_idx,
                _ => usize::MAX,
            };
            for &(q_idx, at_risk_jobs) in &window.at_risk {
                if q_idx < passed_to {
                    caps[q_idx] = caps[q_idx].saturating_add(at_risk_jobs).min(budget);
                }
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
        if let Err(error) = sw_batch.pass(None, query_batches, heaps) {
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
    notes: &FetchNotes<'_>,
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
        x_redo_with::<P, true>(search, lists, heaps, expand, notes)
    } else {
        x_redo_with::<P, false>(search, lists, heaps, expand, notes)
    }
}

/// `x_redo` for the redo `SW` (the Smith-Waterman redo or the X-drop redo).
fn x_redo_with<P, const SW: bool>(
    search: &XRedoParSearch<'_, P>,
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    expand: Option<&mut XExpandAhead<'_>>,
    notes: &FetchNotes<'_>,
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
        search.sw_excluded,
        notes,
        search.program,
        search.do_link_hsps,
        &mut |lists, heaps, window, carried| {
            let (results, lpt_ran) =
                x_compute_window::<P, SW>(search, lists, heaps, window, carried, &shared, lpt);
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
        counts.turned.map_or(String::new(), |turned| format!(
            " turned={turned} continued={} cont_stop={} cont_ahead={} cont_bound={} carried_max={}",
            counts.continued,
            counts.continued_by.stop,
            counts.continued_by.ahead,
            counts.continued_by.bound,
            counts.carried_max
        ))
    );
}

/// EXPERIMENT (LOSAT_X_BPREDOPAR): the redo of a search in jobs of a few matches, into `heaps`
/// (the redone heaps of the search). The preliminary lists are released after their redo, as the
/// one-thread loop releases them. With `compact` (one
/// entry per query), a query's compact list stands for its list until the windows expand it; none
/// is left afterwards. The Smith-Waterman redo's errors read the notes `notes` of the search's
/// fetches (`XSwBatch`). Returns the number of lists released.
pub(crate) fn x_redo_par_run<P>(
    search: &XRedoParSearch<'_, P>,
    lists: &mut [Option<BlastpHitList>],
    compact: Option<&mut [Option<BlastpHitList<XPrelimHsp>>]>,
    heaps: &mut [BlastCompoHeap],
    notes: &FetchNotes<'_>,
) -> Result<usize>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let mut expand = compact.map(|compact| XExpandAhead::new(lists, compact));
    let counts = x_redo(search, lists, heaps, expand.as_mut(), notes)?;
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
