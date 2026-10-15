//! EXPERIMENT (LOSAT_X_BPREDOPAR / LOSAT_X_BPREDOPARSHADOW): blastp's composition-based redo with
//! jobs of a few matches on the search's pool.
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
//!    consecutive of them of one query (one vector of `LOSAT_X_NEWTONLANES`);
//! 2. windows of consecutive jobs, each computed on the pool (one `pool.install`) and then
//!    replayed. A window holds at most 64 jobs per thread (`LOSAT_X_BPREDOPARWINDOW` sets another
//!    number, for tests) and at most one job per thread at risk (below), and it may end inside a
//!    query, so it holds at most 4 x 64 x threads results before their replay. A window boundary
//!    changes nothing: the replay is in stream order and the heaps carry over;
//! 3. a job builds its query's workspace (`build_query_workspace`) if no other job of the window has
//!    (one `OnceLock` per query and window, read only afterwards; a query that spans two windows
//!    has it built again from the same inputs, about 1 us), and its own parameters
//!    (`build_redo_align_params`), which it uses for its matches one after another as the one-thread
//!    loop uses a query's parameters for all its matches: the `Cell` they hold is set during a call
//!    and set back to its previous value at the end of it;
//! 4. with `LOSAT_X_NEWTONLANES`, the job prefetches the first Newton problems of its own matches
//!    against the heaps of its query's batch as the window's start left them
//!    (`x_newton_lanes_prefetch_matches`, as the one-thread loop calls it for the next matches);
//!    then each match runs `postprocess_preliminary_hits` with the arguments of the one-thread loop,
//!    with per-worker scratch (composition workspace, gapped-alignment scratch, and the
//!    subject-range SEG cache backed by one `LOSAT_X_SEGSHARE` store for the search);
//! 5. the replay, per match in order: the early termination test (`BlastCompo_EarlyTermination`)
//!    on the heaps of the query's batch; then the match's error, if it has one (the first error of
//!    a match that the test does not end is the one the one-thread loop reports); then
//!    `BlastCompo_HeapWouldInsert` / `BlastCompo_HeapInsert` on the query's heap. A query's
//!    preliminary list is released as the one-thread loop releases it (`LOSAT_X_BPPFREE`,
//!    `LOSAT_X_BPPCOMPACT`), after the replay of the window that holds its last match.
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
//!   queries all have at least `max_target_seqs` matches. A job with such a match, not ended at
//!   its window's start, is at risk; the other jobs are never ended and are computed ahead freely;
//! - a computed result is dropped only in the window during which its batch's heaps become filled
//!   (at the start of any later window the test is already true for it): at most 4 x threads
//!   results per batch.
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
//! Off unless `LOSAT_X_BPREDOPAR` is set (to any value). Not in the strict set. With
//! `LOSAT_TIMING=1` each search the path runs prints
//! `[TIMING] x_bpredopar: matches=<n> jobs=<j> windows=<w> split=<s> ended=<e> skipped=<k> wasted=<x> par_s=<s> replay_s=<s>`
//! (`split`: windows that end inside a query; `ended`: the matches the early termination test
//! ended, as many as the one-thread loop ends; `skipped`: those of them not computed; `wasted`:
//! those computed and dropped).
//!
//! `LOSAT_X_BPREDOPARSHADOW=1` (shadow mode, implies the switch) runs this path into separate
//! heaps without releasing the preliminary lists, and for every match the replay redoes computes
//! the result again on the calling thread with one set of scratch, as the one-thread loop does,
//! and compares the two results (error or not, scores, E-value bits, every HSP field). Then the
//! one-thread loop runs for the output, and the hit lists of both sets of heaps
//! (`fill_results_from_compo_heaps`) are compared field by field; when the one-thread loop ends
//! with an error instead, the errors are compared (`x_redo_par_shadow_error`). The first
//! difference stops the run (panic) after printing the line that is otherwise printed at exit:
//! `[X_SHADOW] BPREDOPAR searches= queries= matches_compared= lists_compared= diffs=`. With
//! `LOSAT_TIMING=1` the new path's `[TIMING] x_bpredopar:` line is printed as well.
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
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock, PoisonError};
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
use super::encoding::EncodedProtein;
use super::gapalign::{BlastpPreliminaryHsp, GapAlignScratch};
use super::hsp::{
    blast_compo_early_termination, fill_results_from_compo_heaps, BlastCompoHeap, BlastpHitList,
    BlastpHsp, BlastpHspList, COMPO_HEAP_EVALUE_STRETCH,
};
use super::kappa::{
    build_query_workspace, postprocess_preliminary_hits, x_newton_lanes_prefetch_matches,
    x_seg_share, BlastRedoAlignParams, BlastpKappaQueryWorkspace, BlastpKappaSubjectRangeCache,
    BlastpPostprocessResult, XSharedSubjectRanges,
};

/// EXPERIMENT (LOSAT_X_BPREDOPAR): the mode the variables select.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XRedoPar {
    Off,
    On,
    Shadow,
}

// No NCBI counterpart: reads the LOSAT_X_BPREDOPAR / LOSAT_X_BPREDOPARSHADOW switches once; it does
// not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDOPAR): off unless set to any value; `LOSAT_X_BPREDOPARSHADOW` (any
/// value) selects the shadow mode.
pub(crate) fn x_bp_redo_par() -> XRedoPar {
    static MODE: OnceLock<XRedoPar> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BPREDOPARSHADOW").is_some() {
            XRedoPar::Shadow
        } else if std::env::var_os("LOSAT_X_BPREDOPAR").is_some() {
            XRedoPar::On
        } else {
            XRedoPar::Off
        }
    })
}

/// Jobs per thread that a window holds at most, unless `LOSAT_X_BPREDOPARWINDOW` gives another
/// number.
const X_WINDOW_JOBS_PER_THREAD: usize = 64;

// No NCBI counterpart: reads the LOSAT_X_BPREDOPARWINDOW tuning variable once, when the path
// first runs; it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDOPARWINDOW): the jobs per thread of a window when the variable is a
/// positive integer (tests make small windows with it), `X_WINDOW_JOBS_PER_THREAD` otherwise.
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

// No NCBI counterpart: matches per job, one vector of LOSAT_X_NEWTONLANES (four problems); it
// changes no value NCBI computes.
/// Matches of one query per job at most.
const X_JOB_MATCHES: usize = 4;

/// What the redo of a search reads besides the preliminary hit lists and the heaps: the values
/// the one-thread loop of `run_resolved_in_pool` reads.
pub(crate) struct XRedoParSearch<'a, P> {
    pub pool: &'a SearchPool<'a>,
    pub threads: usize,
    pub contexts: &'a [QueryContext],
    /// `blast_engine::query_nomask_sequence`.
    pub query_nomask: fn(&QueryContext) -> &[u8],
    pub search_spaces: &'a [SearchSpace],
    pub subjects: &'a [EncodedProtein],
    /// Each query's batch (`kappa_query_batches`): the heaps its early termination test reads.
    pub query_batches: &'a [Range<usize>],
    pub matrix: ScoringMatrix,
    pub gapped_params: &'a KarlinParams,
    pub gapped_gumbel: &'a BlastGumbelBlk,
    pub evalue: f64,
    /// `build_redo_align_params` of `run_resolved_in_pool`.
    pub build_redo_align_params: &'a P,
    pub timing: bool,
}

/// One match: its query and the position of its list in the query's preliminary hit list.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct XMatch {
    q_idx: usize,
    position: usize,
}

/// The scratch of one worker (or of the shadow's one-thread recomputation).
struct XScratch {
    composition_workspace: BlastCompositionWorkspace,
    gap_scratch: GapAlignScratch,
    hits: Vec<BlastpPreliminaryHsp>,
    range_cache: BlastpKappaSubjectRangeCache,
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
    fn new(range_cache: BlastpKappaSubjectRangeCache) -> Self {
        Self {
            composition_workspace: BlastCompositionWorkspace::new_blosum62(),
            gap_scratch: GapAlignScratch::new(),
            hits: Vec::new(),
            range_cache,
        }
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
}

/// The matches of one query's preliminary hit list: its lists with HSPs.
fn x_list_matches(list: &BlastpHitList) -> usize {
    list.hsplist_array
        .iter()
        .take(list.hsplist_count)
        .filter(|local_match| !local_match.hsps.is_empty())
        .count()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3398-3399,3511
// ```c
//         while (BlastHSPStreamRead(hsp_stream, &localMatch)
//                 != kBlastHSPStream_Eof) {
// ...
//                 BlastHSPList* localMatch = theseMatches[b];
// ```
/// The matches that are redone: those whose list has HSPs.
pub(crate) fn x_matches(lists: &[Option<BlastpHitList>]) -> usize {
    lists.iter().flatten().map(x_list_matches).sum()
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
/// them. The test is true only when every heap of the batch is filled to its cutoff and the
/// E-value exceeds `EVALUE_STRETCH` times each heap's cutoff. In stream order a query's matches
/// are tested after those of the queries before it in its batch and before those of the queries
/// after it, so when they are tested the heaps up to the query hold at most `n` plus their
/// queries' matches, and those after it are as they are now.
fn x_end_above(
    lists: &[Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    query_batches: &[Range<usize>],
) -> Vec<Option<f64>> {
    let matches: Vec<usize> = lists
        .iter()
        .map(|list| list.as_ref().map_or(0, x_list_matches))
        .collect();
    (0..lists.len())
        .map(|q_idx| {
            let batch = query_batches[q_idx].clone();
            // The heaps after the query first: they are not filled unless the query is the
            // last one of its batch, so this ends the test at once for the other queries.
            let may_end = (q_idx + 1..batch.end).all(|other| heaps[other].filled_to_cutoff())
                && (batch.start..=q_idx).all(|other| {
                    heaps[other].n.saturating_add(matches[other]) >= heaps[other].heap_threshold
                });
            may_end.then(|| {
                batch
                    .map(|other| COMPO_HEAP_EVALUE_STRETCH * heaps[other].ecutoff)
                    .fold(f64::NEG_INFINITY, f64::max)
            })
        })
        .collect()
}

/// The bounds of a window: its jobs, and its jobs at risk (those with a match that the early
/// termination test may end and had not ended at the window's start).
#[derive(Clone, Copy, Debug)]
struct XWindowLimits {
    jobs: usize,
    at_risk_jobs: usize,
}

/// A window: consecutive jobs of the stream, computed together and then replayed.
#[derive(Debug)]
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
}

/// The results of a window's matches in order; `None` for a match not computed.
type XWindowResults = Vec<Option<Result<BlastpPostprocessResult>>>;

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
/// the jobs from there while it holds fewer than `limits.jobs` jobs and, for a job at risk, fewer
/// than `limits.at_risk_jobs` jobs at risk. A match whose test is already true is marked ended
/// (`end_above`: `x_end_above`); the test stays true until its replay.
fn x_window(
    lists: &[Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    query_batches: &[Range<usize>],
    end_above: &[Option<f64>],
    start: XMatch,
    limits: XWindowLimits,
) -> XWindow {
    let mut window = XWindow {
        matches: Vec::new(),
        ended: Vec::new(),
        jobs: Vec::new(),
        next: None,
    };
    let mut at_risk_jobs = 0;
    let mut cursor = start;
    let mut job = Vec::with_capacity(X_JOB_MATCHES);
    loop {
        x_next_job(lists, &mut cursor, &mut job);
        let Some(&first) = job.first() else {
            return window;
        };
        let mark = window.matches.len();
        let mut at_risk = false;
        for &matched in &job {
            let best_evalue = x_local_match(lists, matched).best_evalue;
            let may_end = end_above[matched.q_idx].is_some_and(|above| best_evalue > above);
            let ended = may_end
                && blast_compo_early_termination(
                    best_evalue,
                    &heaps[query_batches[matched.q_idx].clone()],
                );
            at_risk |= may_end && !ended;
            window.ended.push(ended);
        }
        // A window holds at least one job.
        if !window.jobs.is_empty()
            && (window.jobs.len() >= limits.jobs
                || (at_risk && at_risk_jobs >= limits.at_risk_jobs))
        {
            window.ended.truncate(mark);
            window.next = Some(first);
            return window;
        }
        at_risk_jobs += usize::from(at_risk);
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
    let subject = &search.subjects[local_match.oid as usize];
    let subject_raw = &subject.aa_seq[1..subject.aa_seq.len() - 1];
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
    )
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
/// matches of one query; the results in match order, `None` for a match ended before the window
/// (not computed). `heaps` are the heaps as the window's start left them; nothing changes them
/// while the jobs run. NCBI's loop also runs the early termination test and the heap insertion of
/// each match on its thread; here they are replayed afterwards (`x_replay`).
fn x_compute_window<P>(
    search: &XRedoParSearch<'_, P>,
    lists: &[Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    window: &XWindow,
    shared: &Option<Arc<XSharedSubjectRanges>>,
) -> XWindowResults
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let (Some(first), Some(last)) = (window.matches.first(), window.matches.last()) else {
        return Vec::new();
    };
    if window.ended.iter().all(|&ended| ended) {
        return window.ended.iter().map(|_| None).collect();
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
    // EXPERIMENT (LOSAT_X_NEWTONLANES): read as the one-thread loop reads it.
    let newton_lanes = crate::core::composition_adjustment::x_newton_lanes::mode() != 0;
    let per_job: Vec<XWindowResults> = search.pool.install(|| {
        window
            .jobs
            .par_iter()
            .map_init(
                // Each rayon work unit owns its scratch; how many are made depends on the
                // scheduler, and no value depends on it.
                || XScratch::new(BlastpKappaSubjectRangeCache::x_with_shared(shared.clone())),
                |scratch, job| {
                    let job_matches = &window.matches[job.clone()];
                    let job_ended = &window.ended[job.clone()];
                    if job_ended.iter().all(|&ended| ended) {
                        return job_ended.iter().map(|_| None).collect();
                    }
                    let q_idx = job_matches[0].q_idx;
                    let query_workspace = workspaces[q_idx - queries.start]
                        .get_or_init(|| x_query_workspace(search, q_idx));
                    // The parameters of the one-thread loop for this query, used for the job's
                    // matches one after another as that loop uses them for all of the query's.
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
                    // ```c
                    //                 if (BlastCompo_EarlyTermination(
                    //                         localMatch->best_evalue,
                    //                         redoneMatches,
                    //                         numQueries
                    //                 )) {
                    // ```
                    // No NCBI counterpart: LOSAT_X_NEWTONLANES solves the first Newton problems of
                    // the job's matches, four per vector (kappa::x_newton_lanes_prefetch_matches);
                    // a result is used only by a call on this thread with the same input bits. It
                    // does not change any value NCBI computes. Given the heaps of the batch as the
                    // window's start left them, as the one-thread loop gives its heaps, it leaves
                    // out the matches the test ends on them: those not computed here.
                    if newton_lanes {
                        let list = lists[q_idx]
                            .as_ref()
                            .expect("a match's query has its preliminary hit list");
                        x_newton_lanes_prefetch_matches(
                            &list.hsplist_array[job_matches[0].position
                                ..=job_matches[job_matches.len() - 1].position],
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
                },
            )
            .collect()
    });
    per_job.into_iter().flatten().collect()
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
/// The steps of the one-thread loop that depend on order, for a window's matches in stream order,
/// with the results computed for them: the early termination test on the heaps of the query's
/// batch, the match's error, and the heap insertion. `redone` sees each match that the test does
/// not end, before its result is used. Returns the numbers of matches that the test ended without
/// a result (not computed) and with one (computed and dropped).
fn x_replay(
    matches: &[XMatch],
    results: XWindowResults,
    lists: &[Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    query_batches: &[Range<usize>],
    redone: &mut dyn FnMut(XMatch, &BlastpHspList, &Result<BlastpPostprocessResult>),
) -> Result<(usize, usize)> {
    assert_eq!(matches.len(), results.len(), "one result per match");
    let (mut skipped, mut wasted) = (0, 0);
    for (&matched, postprocessed) in matches.iter().zip(results) {
        let local_match = x_local_match(lists, matched);
        if blast_compo_early_termination(
            local_match.best_evalue,
            &heaps[query_batches[matched.q_idx].clone()],
        ) {
            if postprocessed.is_some() {
                wasted += 1;
            } else {
                skipped += 1;
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
        redone(matched, local_match, &postprocessed);
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
    Ok((skipped, wasted))
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
/// of the windows before left them) and then replayed in stream order into `heaps`. With
/// `release`, a query's preliminary list is released after the replay of the window that holds
/// its last match. `redone` sees each match that the replay redoes.
#[allow(clippy::too_many_arguments)]
fn x_redo_windows(
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    query_batches: &[Range<usize>],
    limits: XWindowLimits,
    release: bool,
    timing: bool,
    compute: &mut dyn FnMut(
        &[Option<BlastpHitList>],
        &[BlastCompoHeap],
        &XWindow,
    ) -> XWindowResults,
    redone: &mut dyn FnMut(XMatch, &BlastpHspList, &Result<BlastpPostprocessResult>),
) -> Result<XCounts> {
    let end_above = x_end_above(lists, heaps, query_batches);
    let mut counts = XCounts::default();
    let mut start = Some(XMatch {
        q_idx: 0,
        position: 0,
    });
    let mut released_to = 0;
    while let Some(cursor) = start {
        let window = x_window(lists, heaps, query_batches, &end_above, cursor, limits);
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
            let (skipped, wasted) = x_replay(
                &window.matches,
                results,
                lists,
                heaps,
                query_batches,
                redone,
            )?;
            counts.skipped += skipped;
            counts.wasted += wasted;
            if let Some(replay_start) = replay_start {
                counts.replay_s += replay_start.elapsed().as_secs_f64();
            }
        }
        // EXPERIMENT (LOSAT_X_BPPFREE / LOSAT_X_BPPCOMPACT): nothing reads a query's preliminary
        // list after its redo; the one-thread loop releases each one after its query
        // (`x_prelim_free.rs`). The queries before the next window's first match are done.
        let done = window.next.map_or(lists.len(), |next| next.q_idx);
        if release {
            for list in &mut lists[released_to..done] {
                if list.take().is_some() {
                    counts.released += 1;
                }
            }
        }
        released_to = done;
        start = window.next;
    }
    Ok(counts)
}

/// The redo of a search: windows of jobs computed on the pool and replayed in stream order into
/// `heaps`; with `release`, each query's preliminary list is released after its redo.
fn x_redo<P>(
    search: &XRedoParSearch<'_, P>,
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    release: bool,
    mut reference: Option<&mut XSerialReference>,
) -> Result<XCounts>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    // The one-thread loop builds the parameters of query 0 before any match. Their only error
    // (the initial matrix, kept once for the search) is the same for every query, so it is
    // reported here, before any match, as there.
    (search.build_redo_align_params)(0)?;
    // EXPERIMENT (LOSAT_X_SEGSHARE): one store of masked subject ranges for the search, as in
    // the query-parallel redo (`XSharedSubjectRanges`).
    let shared = (x_seg_share() != 0).then(|| Arc::new(XSharedSubjectRanges::new()));
    let threads = search.threads.max(1);
    let limits = XWindowLimits {
        jobs: x_window_jobs_per_thread().saturating_mul(threads),
        at_risk_jobs: threads,
    };
    x_redo_windows(
        lists,
        heaps,
        search.query_batches,
        limits,
        release,
        search.timing,
        &mut |lists, heaps, window| x_compute_window(search, lists, heaps, window, &shared),
        &mut |matched, local_match, computed| {
            if let Some(reference) = reference.as_deref_mut() {
                reference.compare(search, matched, local_match, computed);
            }
        },
    )
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
fn x_print_timing(counts: &XCounts) {
    eprintln!(
        "[TIMING] x_bpredopar: matches={} jobs={} windows={} split={} ended={} skipped={} wasted={} par_s={:.6} replay_s={:.6}",
        counts.matches,
        counts.jobs,
        counts.windows,
        counts.split,
        counts.skipped + counts.wasted,
        counts.skipped,
        counts.wasted,
        counts.par_s,
        counts.replay_s
    );
}

/// EXPERIMENT (LOSAT_X_BPREDOPAR): the redo of a search in jobs of a few matches, into `heaps`
/// (the redone heaps of the search). With `release`, the preliminary lists are released after
/// their redo, as the one-thread loop releases them. Returns the number of lists released.
pub(crate) fn x_redo_par_run<P>(
    search: &XRedoParSearch<'_, P>,
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    release: bool,
) -> Result<usize>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let counts = x_redo(search, lists, heaps, release, None)?;
    if search.timing {
        x_print_timing(&counts);
    }
    Ok(counts.released)
}

static X_SHADOW_SEARCHES: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_QUERIES: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_MATCHES: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_LISTS: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_DIFFS: AtomicU64 = AtomicU64::new(0);
/// The new path's outcome for a search whose one-thread loop returned its error before the hit
/// lists were compared (`XRedoParShadow`'s drop), until `x_redo_par_shadow_error` compares it.
static X_SHADOW_UNCOMPARED: Mutex<Option<std::result::Result<(), String>>> = Mutex::new(None);

// No NCBI counterpart: the shadow mode's check stops the run on a difference.
fn x_shadow_difference(message: String) -> ! {
    X_SHADOW_DIFFS.fetch_add(1, Ordering::Relaxed);
    // The panic ends the run before `main` prints the line at exit: it is printed here, with the
    // difference counted.
    x_shadow_print_line();
    panic!("[X_SHADOW] BPREDOPAR {message}");
}

/// EXPERIMENT (LOSAT_X_BPREDOPARSHADOW): the one-thread loop's computation of the matches that
/// the replay redoes, in stream order, with one set of scratch for the search and the parameters
/// and workspace built once per query, as the one-thread loop holds them.
struct XSerialReference {
    scratch: XScratch,
    query: Option<(
        usize,
        std::result::Result<BlastRedoAlignParams, String>,
        BlastpKappaQueryWorkspace,
    )>,
    matches_compared: u64,
}

impl XSerialReference {
    fn new() -> Self {
        Self {
            // The one-thread loop's cache: `BlastpKappaSubjectRangeCache::new()`.
            scratch: XScratch::new(BlastpKappaSubjectRangeCache::new()),
            query: None,
            matches_compared: 0,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3624-3628
    // ```c
    //                         *pStatusCode =
    //                                 Blast_RedoOneMatch(
    //                                         alignments,             // thread-local
    //                                         redo_align_params,      // thread-local
    //                                         incoming_aligns,        // thread-local
    // ```
    // No NCBI counterpart for the comparison: it computes the match again and stops the run on a
    // difference; it changes no value NCBI computes.
    fn compare<P>(
        &mut self,
        search: &XRedoParSearch<'_, P>,
        matched: XMatch,
        local_match: &BlastpHspList,
        computed: &Result<BlastpPostprocessResult>,
    ) where
        P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
    {
        let Self { scratch, query, .. } = self;
        if query.as_ref().map(|(q_idx, ..)| *q_idx) != Some(matched.q_idx) {
            *query = Some((
                matched.q_idx,
                (search.build_redo_align_params)(matched.q_idx)
                    .map_err(|error| format!("{error:#}")),
                x_query_workspace(search, matched.q_idx),
            ));
        }
        let (_, params, query_workspace) = query.as_ref().expect("set above");
        let reference = match params {
            Ok(params) => x_redo_match(
                search,
                matched.q_idx,
                local_match,
                query_workspace,
                params,
                scratch,
            )
            .map_err(|error| format!("{error:#}")),
            Err(error) => Err(error.clone()),
        };
        if let Some(difference) = x_result_difference(computed, &reference) {
            x_shadow_difference(format!(
                "query {} match {} (subject {}): {difference}",
                matched.q_idx, matched.position, local_match.oid
            ));
        }
        self.matches_compared += 1;
    }
}

// No NCBI counterpart: the shadow mode's comparison (floating-point values by their bits).
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

// No NCBI counterpart: the shadow mode's comparison.
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

// No NCBI counterpart: the shadow mode's comparison of a computed result (first) with the
// one-thread computation's (second).
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

/// EXPERIMENT (LOSAT_X_BPREDOPARSHADOW): the new path's heaps of a search, kept until the
/// one-thread loop has filled the search's heaps.
pub(crate) struct XRedoParShadow {
    heaps: Vec<BlastCompoHeap>,
    outcome: std::result::Result<(), String>,
    hitlist_size: usize,
    compared: bool,
}

impl Drop for XRedoParShadow {
    // No NCBI counterpart: a check of the shadow mode; it changes no value NCBI computes.
    /// Dropped without `compare` (and not by a panic): the one-thread loop returned its error
    /// (`?`) before its hit lists were compared. Its caller compares the errors
    /// (`x_redo_par_shadow_error`).
    fn drop(&mut self) {
        if !self.compared && !std::thread::panicking() {
            *X_SHADOW_UNCOMPARED
                .lock()
                .unwrap_or_else(PoisonError::into_inner) =
                Some(std::mem::replace(&mut self.outcome, Ok(())));
        }
    }
}

// No NCBI counterpart: a check of the shadow mode; it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDOPARSHADOW): `result`, a search's result (`run_resolved_in_pool`),
/// against the new path's outcome when the one-thread loop of the search ended with an error
/// before the hit lists were compared: the new path must end with the same error.
pub(crate) fn x_redo_par_shadow_error(result: &Result<()>) {
    let Some(outcome) = X_SHADOW_UNCOMPARED
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .take()
    else {
        return;
    };
    // Only the one-thread loop's error leaves the search between the shadow and the comparison.
    let Err(error) = result else {
        return;
    };
    let serial = format!("{error:#}");
    match outcome {
        Ok(()) => x_shadow_difference(format!(
            "the one-thread loop ended with an error where the new path did not: {serial}"
        )),
        Err(computed) if computed != serial => {
            x_shadow_difference(format!("errors {computed} / {serial}"))
        }
        Err(_) => {
            X_SHADOW_SEARCHES.fetch_add(1, Ordering::Relaxed);
        }
    }
}

/// EXPERIMENT (LOSAT_X_BPREDOPARSHADOW): the new path into a copy of `heaps` (the search's heaps
/// before the redo), the preliminary lists kept, every match it redoes compared with the
/// one-thread computation. An error ends the path as it ends the one-thread loop and is kept for
/// the comparison.
pub(crate) fn x_redo_par_shadow_run<P>(
    search: &XRedoParSearch<'_, P>,
    lists: &mut [Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    hitlist_size: usize,
) -> XRedoParShadow
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
    let mut shadow_heaps = heaps.to_vec();
    let mut reference = XSerialReference::new();
    let outcome = x_redo(
        search,
        lists,
        &mut shadow_heaps,
        false,
        Some(&mut reference),
    )
    .map(|counts| {
        if search.timing {
            x_print_timing(&counts);
        }
    })
    .map_err(|error| format!("{error:#}"));
    X_SHADOW_MATCHES.fetch_add(reference.matches_compared, Ordering::Relaxed);
    XRedoParShadow {
        heaps: shadow_heaps,
        outcome,
        hitlist_size,
        compared: false,
    }
}

impl XRedoParShadow {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2493-2514
    // ```c
    // s_FillResultsFromCompoHeaps(BlastHSPResults * results,
    //                             BlastCompo_Heap heaps[],
    //                             Int4 hitlist_size)
    // ```
    // No NCBI counterpart for the comparison: it changes no value NCBI computes.
    /// The hit lists of the new path's heaps against those of `serial_heaps` (the one-thread
    /// loop's, which stay as they are), every field; floating-point values by their bits.
    pub(crate) fn compare(mut self, serial_heaps: &[BlastCompoHeap]) {
        self.compared = true;
        if let Err(error) = &self.outcome {
            x_shadow_difference(format!(
                "the new path ended with an error where the one-thread loop did not: {error}"
            ));
        }
        let mut serial_heaps = serial_heaps.to_vec();
        let reference = fill_results_from_compo_heaps(self.hitlist_size, &mut serial_heaps);
        let computed = fill_results_from_compo_heaps(self.hitlist_size, &mut self.heaps);
        if reference.len() != computed.len() {
            x_shadow_difference(format!("queries {} / {}", computed.len(), reference.len()));
        }
        let mut lists = 0u64;
        for (q_idx, (computed, reference)) in computed.iter().zip(&reference).enumerate() {
            let (computed, reference) = match (computed, reference) {
                (None, None) => continue,
                (Some(computed), Some(reference)) => (computed, reference),
                (computed, reference) => x_shadow_difference(format!(
                    "query {q_idx}: hit list {} / {}",
                    if computed.is_some() {
                        "present"
                    } else {
                        "absent"
                    },
                    if reference.is_some() {
                        "present"
                    } else {
                        "absent"
                    }
                )),
            };
            let header = |list: &BlastpHitList| {
                (
                    list.hsplist_count,
                    list.hsplist_max,
                    list.worst_evalue.to_bits(),
                    list.low_score,
                    list.heapified,
                    list.hsplist_array.len(),
                    list.hsplist_current,
                    list.num_hits,
                )
            };
            if header(computed) != header(reference) {
                x_shadow_difference(format!(
                    "query {q_idx}: hit list fields {:?} / {:?}",
                    header(computed),
                    header(reference)
                ));
            }
            for (position, (a, b)) in computed
                .hsplist_array
                .iter()
                .zip(&reference.hsplist_array)
                .enumerate()
            {
                if let Some(difference) = x_list_difference(a, b) {
                    x_shadow_difference(format!("query {q_idx} list {position}: {difference}"));
                }
            }
            lists += computed.hsplist_array.len() as u64;
        }
        X_SHADOW_SEARCHES.fetch_add(1, Ordering::Relaxed);
        X_SHADOW_QUERIES.fetch_add(reference.len() as u64, Ordering::Relaxed);
        X_SHADOW_LISTS.fetch_add(lists, Ordering::Relaxed);
    }
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
fn x_shadow_print_line() {
    eprintln!(
        "[X_SHADOW] BPREDOPAR searches={} queries={} matches_compared={} lists_compared={} diffs={}",
        X_SHADOW_SEARCHES.load(Ordering::Relaxed),
        X_SHADOW_QUERIES.load(Ordering::Relaxed),
        X_SHADOW_MATCHES.load(Ordering::Relaxed),
        X_SHADOW_LISTS.load(Ordering::Relaxed),
        X_SHADOW_DIFFS.load(Ordering::Relaxed),
    );
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDOPARSHADOW): the line `main` prints at exit.
pub(crate) fn x_bp_redo_par_print_summary() {
    if x_bp_redo_par() == XRedoPar::Shadow {
        x_shadow_print_line();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ECUTOFF: f64 = 0.002;

    fn hsp(oid: u32, e_value: f64) -> BlastpHsp {
        BlastpHsp {
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

    const UNLIMITED: XWindowLimits = XWindowLimits {
        jobs: usize::MAX,
        at_risk_jobs: usize::MAX,
    };

    /// All matches of `lists` in stream order.
    fn stream(lists: &[Option<BlastpHitList>]) -> Vec<XMatch> {
        let window = x_window(
            lists,
            &heaps_of(lists.len(), 1),
            &alone(lists.len()),
            &vec![None; lists.len()],
            at(0, 0),
            UNLIMITED,
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
        assert_eq!(x_matches(&lists), 9);
        let (heaps, batches, none) = (heaps_of(5, 1), alone(5), vec![None; 5]);
        let window = |start, jobs| {
            x_window(
                &lists,
                &heaps,
                &batches,
                &none,
                start,
                XWindowLimits {
                    jobs,
                    at_risk_jobs: 1,
                },
            )
        };
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
            Some(redone(10, 1e-20)),
            Some(redone(20, 1e-20)),
            Some(Err(anyhow::anyhow!("ended"))),
            Some(Err(anyhow::anyhow!("reported"))),
        ];
        let mut seen = Vec::new();
        let error = x_replay(
            &matches,
            results,
            &lists,
            &mut heaps,
            &batches,
            &mut |matched, _, _| seen.push(matched.position),
        )
        .expect_err("the redone match's error");
        assert_eq!(error.to_string(), "reported");
        assert_eq!(seen, [0, 0, 2]);
        assert_eq!((heaps[0].n, heaps[1].n), (1, 1));
    }

    #[test]
    fn the_replay_inserts_in_stream_order_and_counts_ended_results() {
        let (lists, heaps, batches) = two_query_batch();
        let matches = stream(&lists);
        // The third match is ended with a result (computed) or without one (not computed).
        for third in [Some(Err(anyhow::anyhow!("ended"))), None] {
            let computed = third.is_some();
            let mut heaps = heaps.clone();
            let results = vec![
                Some(redone(10, 1e-20)),
                Some(redone(20, 1e-20)),
                third,
                Some(redone(22, 0.005)),
            ];
            let ended = x_replay(
                &matches,
                results,
                &lists,
                &mut heaps,
                &batches,
                &mut |_, _, _| {},
            )
            .expect("no error is redone");
            assert_eq!(
                ended,
                if computed { (0, 1) } else { (1, 0) },
                "skipped, wasted"
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
        let _ = x_replay(
            &matches,
            results,
            &lists,
            &mut heaps,
            &batches,
            &mut |_, _, _| {},
        );
    }

    #[test]
    fn a_query_alone_in_its_batch_reads_only_its_own_heap() {
        // Query 0 has no match, so in a batch with query 1 the test is never true for query 1;
        // alone in its batch, query 1's second match is ended once its own heap is filled.
        let lists = vec![
            None,
            hit_list(vec![list(20, 1e-20, false), list(21, 1.0, false)]),
        ];
        for (batches, wasted) in [(vec![0..2, 0..2], 0), (vec![0..1, 1..2], 1)] {
            let mut heaps = heaps_of(2, 1);
            let matches = stream(&lists);
            let results = vec![Some(redone(20, 1e-20)), Some(redone(21, 1.0))];
            assert_eq!(
                x_replay(
                    &matches,
                    results,
                    &lists,
                    &mut heaps,
                    &batches,
                    &mut |_, _, _| {}
                )
                .expect("no error"),
                (0, wasted),
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
            x_end_above(&lists, &heaps_of(5, 1), &batches),
            [None, None, stretch, None, None]
        );
        // Two: query 1 cannot fill its heap.
        assert_eq!(x_end_above(&lists, &heaps_of(5, 2), &batches), [None; 5]);
        // A heap after a query that is filled already counts as such.
        let mut heaps = heaps_of(5, 1);
        heaps[2].insert(list(1, 1e-9, false), 1e-9, 40, 1);
        assert_eq!(
            x_end_above(&lists, &heaps, &batches),
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
        mut window_seen: impl FnMut(&[Option<BlastpHitList>], &XWindow),
    ) -> (XCounts, usize) {
        let mut computed = 0;
        let counts = x_redo_windows(
            lists,
            heaps,
            batches,
            limits,
            true,
            false,
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
                            redone(local_match.oid, local_match.best_evalue)
                        })
                    })
                    .collect()
            },
            &mut |_, _, _| {},
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
        // The first window holds the jobs at risk that the threads take (4 x 4 matches); the
        // test ends 15 of them at the replay, and every later match before its window.
        assert_eq!(computed, X_JOB_MATCHES * THREADS);
        assert_eq!(counts.wasted, X_JOB_MATCHES * THREADS - 1);
        assert_eq!(counts.skipped, MATCHES as usize - X_JOB_MATCHES * THREADS);
        assert_eq!(counts.skipped + counts.wasted, MATCHES as usize - 1);
        // The later windows hold their whole job budget (1,024 matches).
        assert_eq!((counts.windows, counts.split), (3, 2));
        assert_eq!((counts.matches, counts.released), (MATCHES as usize, 1));
        assert_eq!(heaps[0].n, 1);
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

    #[test]
    fn the_shadow_compares_the_errors_when_the_one_thread_loop_ends_with_one() {
        // A search's shadow dropped before its comparison: the one-thread loop returned an error.
        let uncompared = |outcome| {
            drop(XRedoParShadow {
                heaps: Vec::new(),
                outcome,
                hitlist_size: 1,
                compared: false,
            })
        };
        let stops = |result: Result<()>| {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                x_redo_par_shadow_error(&result)
            }))
            .is_err()
        };
        uncompared(Err("failed".to_string()));
        assert!(!stops(Err(anyhow::anyhow!("failed"))), "the same error");
        uncompared(Err("failed".to_string()));
        assert!(stops(Err(anyhow::anyhow!("other"))), "another error");
        uncompared(Ok(()));
        assert!(
            stops(Err(anyhow::anyhow!("failed"))),
            "an error of one path"
        );
        // Nothing is left to compare.
        assert!(!stops(Err(anyhow::anyhow!("failed"))));
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
}
