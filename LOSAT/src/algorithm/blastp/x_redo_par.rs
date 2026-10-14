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
//! 2. a job builds its query's workspace (`build_query_workspace`) if no other job of the window has
//!    (one `OnceLock` per query, read only afterwards), and its own parameters
//!    (`build_redo_align_params`), which it uses for its matches one after another as the one-thread
//!    loop uses a query's parameters for all its matches: the `Cell` they hold is set during a call
//!    and set back to its previous value at the end of it;
//! 3. with `LOSAT_X_NEWTONLANES`, the job prefetches the first Newton problems of its own matches
//!    (`x_newton_lanes_prefetch_matches`, as the one-thread loop calls it for the next matches);
//!    then each match runs `postprocess_preliminary_hits` with the arguments of the one-thread loop,
//!    with per-worker scratch (composition workspace, gapped-alignment scratch, and the
//!    subject-range SEG cache backed by one `LOSAT_X_SEGSHARE` store for the search);
//! 4. the replay, per match in order: the early termination test (`BlastCompo_EarlyTermination`)
//!    on the heaps of the query's batch; then the match's error, if it has one (the first error of
//!    a match that the test does not end is the one the one-thread loop reports); then
//!    `BlastCompo_HeapWouldInsert` / `BlastCompo_HeapInsert` on the query's heap. After a window's
//!    replay, its queries' preliminary lists are released as the one-thread loop releases them
//!    (`LOSAT_X_BPPFREE`, `LOSAT_X_BPPCOMPACT`);
//! 5. windows of whole queries (at least 64 jobs per thread, closed at a query boundary; one
//!    `pool.install` each) bound the results held before their replay. A window boundary changes
//!    nothing: the replay is in stream order.
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
//! error is reported) is replayed on the same heaps in the same sequence. A result that the test
//! ends is computed and dropped: the test can be true only for the last query of a batch, once
//! every heap of the batch is filled to its cutoff.
//!
//! Off unless `LOSAT_X_BPREDOPAR` is set (to any value). Not in the strict set. With
//! `LOSAT_TIMING=1` each search the path runs prints
//! `[TIMING] x_bpredopar: matches=<n> jobs=<j> windows=<w> ended=<e> par_s=<s> replay_s=<s>`
//! (`ended`: the computed results that the early termination test dropped).
//!
//! `LOSAT_X_BPREDOPARSHADOW=1` (shadow mode, implies the switch) runs this path into separate
//! heaps without releasing the preliminary lists, and for every match the replay redoes computes
//! the result again on the calling thread with one set of scratch, as the one-thread loop does,
//! and compares the two results (error or not, scores, E-value bits, every HSP field). Then the
//! one-thread loop runs for the output, and the hit lists of both sets of heaps
//! (`fill_results_from_compo_heaps`) are compared field by field. The first difference stops the
//! run (panic). One line at exit:
//! `[X_SHADOW] BPREDOPAR searches= queries= matches_compared= lists_compared= diffs=`.
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
use super::encoding::EncodedProtein;
use super::gapalign::{BlastpPreliminaryHsp, GapAlignScratch};
use super::hsp::{
    blast_compo_early_termination, fill_results_from_compo_heaps, BlastCompoHeap, BlastpHitList,
    BlastpHsp, BlastpHspList,
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

/// Jobs per thread that a window holds at least; it closes at the next query boundary.
const X_WINDOW_JOBS_PER_THREAD: usize = 64;

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

#[derive(Default)]
struct XCounts {
    matches: usize,
    jobs: usize,
    windows: usize,
    ended: usize,
    released: usize,
    par_s: f64,
    replay_s: f64,
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
    lists
        .iter()
        .flatten()
        .map(|list| {
            list.hsplist_array
                .iter()
                .take(list.hsplist_count)
                .filter(|local_match| !local_match.hsps.is_empty())
                .count()
        })
        .sum()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3398-3399,3421
// ```c
//         while (BlastHSPStreamRead(hsp_stream, &localMatch)
//                 != kBlastHSPStream_Eof) {
// ...
//             theseMatches[i] = head->match;
// ```
/// The next window: the queries from `start` until it holds at least `min_jobs` jobs (or the
/// last query). Returns the end of its query range, its matches in stream order, and its jobs as
/// ranges of those matches: up to `X_JOB_MATCHES` consecutive matches of one query each.
fn x_window(
    lists: &[Option<BlastpHitList>],
    start: usize,
    min_jobs: usize,
) -> (usize, Vec<XMatch>, Vec<Range<usize>>) {
    let mut matches = Vec::new();
    let mut jobs = Vec::new();
    let mut end = start;
    while end < lists.len() && (end == start || jobs.len() < min_jobs) {
        if let Some(list) = lists[end].as_ref() {
            let first = matches.len();
            for (position, local_match) in list
                .hsplist_array
                .iter()
                .take(list.hsplist_count)
                .enumerate()
            {
                if !local_match.hsps.is_empty() {
                    matches.push(XMatch {
                        q_idx: end,
                        position,
                    });
                }
            }
            let mut job_start = first;
            while job_start < matches.len() {
                let job_end = (job_start + X_JOB_MATCHES).min(matches.len());
                jobs.push(job_start..job_end);
                job_start = job_end;
            }
        }
        end += 1;
    }
    (end, matches, jobs)
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
/// matches of one query; the results in match order. NCBI's loop also runs the early termination
/// test and the heap insertion of each match on its thread; here they are replayed afterwards
/// (`x_replay`).
fn x_compute_window<P>(
    search: &XRedoParSearch<'_, P>,
    lists: &[Option<BlastpHitList>],
    window: Range<usize>,
    matches: &[XMatch],
    jobs: &[Range<usize>],
    shared: &Option<Arc<XSharedSubjectRanges>>,
) -> Vec<Result<BlastpPostprocessResult>>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
{
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
        window.clone().map(|_| OnceLock::new()).collect();
    // EXPERIMENT (LOSAT_X_NEWTONLANES): read as the one-thread loop reads it.
    let newton_lanes = crate::core::composition_adjustment::x_newton_lanes::mode() != 0;
    // A heap that no match has reached (not filled to its cutoff): the prefetch leaves out only
    // the matches that the test ends for the heaps it is given, and a job cannot read the heaps of
    // the replay; so it prefetches all its matches.
    let untouched_heap = [BlastCompoHeap::new(1, 0.0)];
    let per_job: Vec<Vec<Result<BlastpPostprocessResult>>> = search.pool.install(|| {
        jobs.par_iter()
            .map_init(
                // Each rayon work unit owns its scratch; how many are made depends on the
                // scheduler, and no value depends on it.
                || XScratch::new(BlastpKappaSubjectRangeCache::x_with_shared(shared.clone())),
                |scratch, job| {
                    let job_matches = &matches[job.clone()];
                    let q_idx = job_matches[0].q_idx;
                    let query_workspace = workspaces[q_idx - window.start]
                        .get_or_init(|| x_query_workspace(search, q_idx));
                    // The parameters of the one-thread loop for this query, used for the job's
                    // matches one after another as that loop uses them for all of the query's.
                    let redo_align_params = match (search.build_redo_align_params)(q_idx) {
                        Ok(params) => params,
                        Err(error) => {
                            return job_matches
                                .iter()
                                .map(|_| Err(anyhow::anyhow!("{error:#}")))
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
                    // does not change any value NCBI computes.
                    if newton_lanes {
                        let list = lists[q_idx]
                            .as_ref()
                            .expect("a match's query has its preliminary hit list");
                        x_newton_lanes_prefetch_matches(
                            &list.hsplist_array[job_matches[0].position
                                ..=job_matches[job_matches.len() - 1].position],
                            search.subjects,
                            &untouched_heap,
                            query_workspace,
                            &scratch.composition_workspace,
                            &redo_align_params,
                            &mut scratch.range_cache,
                        );
                    }
                    job_matches
                        .iter()
                        .map(|&matched| {
                            x_redo_match(
                                search,
                                q_idx,
                                x_local_match(lists, matched),
                                query_workspace,
                                &redo_align_params,
                                scratch,
                            )
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
/// not end, before its result is used. Returns the number of results that the test dropped.
fn x_replay(
    matches: &[XMatch],
    results: Vec<Result<BlastpPostprocessResult>>,
    lists: &[Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    query_batches: &[Range<usize>],
    redone: &mut dyn FnMut(XMatch, &BlastpHspList, &Result<BlastpPostprocessResult>),
) -> Result<usize> {
    assert_eq!(matches.len(), results.len(), "one result per match");
    let mut ended = 0;
    for (&matched, postprocessed) in matches.iter().zip(results) {
        let local_match = x_local_match(lists, matched);
        if blast_compo_early_termination(
            local_match.best_evalue,
            &heaps[query_batches[matched.q_idx].clone()],
        ) {
            ended += 1;
            continue;
        }
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
    Ok(ended)
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
/// The redo of a search: windows of jobs computed on the pool and replayed in stream order into
/// `heaps`; with `release`, each window's preliminary lists are released after its replay.
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
    let min_jobs = X_WINDOW_JOBS_PER_THREAD.saturating_mul(search.threads.max(1));
    let mut counts = XCounts::default();
    let mut start = 0;
    while start < lists.len() {
        let (end, matches, jobs) = x_window(lists, start, min_jobs);
        counts.windows += 1;
        counts.matches += matches.len();
        counts.jobs += jobs.len();
        if !jobs.is_empty() {
            let par_start = Instant::now();
            let results = x_compute_window(search, lists, start..end, &matches, &jobs, &shared);
            counts.par_s += par_start.elapsed().as_secs_f64();
            let replay_start = Instant::now();
            counts.ended += x_replay(
                &matches,
                results,
                lists,
                heaps,
                search.query_batches,
                &mut |matched, local_match, computed| {
                    if let Some(reference) = reference.as_deref_mut() {
                        reference.compare(search, matched, local_match, computed);
                    }
                },
            )?;
            counts.replay_s += replay_start.elapsed().as_secs_f64();
        }
        // EXPERIMENT (LOSAT_X_BPPFREE / LOSAT_X_BPPCOMPACT): nothing reads these queries'
        // preliminary lists after their redo; the one-thread loop releases each one after its
        // query (`x_prelim_free.rs`).
        if release {
            for list in &mut lists[start..end] {
                if list.take().is_some() {
                    counts.released += 1;
                }
            }
        }
        start = end;
    }
    Ok(counts)
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
        eprintln!(
            "[TIMING] x_bpredopar: matches={} jobs={} windows={} ended={} par_s={:.6} replay_s={:.6}",
            counts.matches,
            counts.jobs,
            counts.windows,
            counts.ended,
            counts.par_s,
            counts.replay_s
        );
    }
    Ok(counts.released)
}

static X_SHADOW_SEARCHES: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_QUERIES: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_MATCHES: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_LISTS: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_DIFFS: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: the shadow mode's check stops the run on a difference.
fn x_shadow_difference(message: String) -> ! {
    X_SHADOW_DIFFS.fetch_add(1, Ordering::Relaxed);
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
    .map(|_| ())
    .map_err(|error| format!("{error:#}"));
    X_SHADOW_MATCHES.fetch_add(reference.matches_compared, Ordering::Relaxed);
    XRedoParShadow {
        heaps: shadow_heaps,
        outcome,
        hitlist_size,
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
/// EXPERIMENT (LOSAT_X_BPREDOPARSHADOW): the line `main` prints at exit.
pub(crate) fn x_bp_redo_par_print_summary() {
    if x_bp_redo_par() == XRedoPar::Shadow {
        eprintln!(
            "[X_SHADOW] BPREDOPAR searches={} queries={} matches_compared={} lists_compared={} diffs={}",
            X_SHADOW_SEARCHES.load(Ordering::Relaxed),
            X_SHADOW_QUERIES.load(Ordering::Relaxed),
            X_SHADOW_MATCHES.load(Ordering::Relaxed),
            X_SHADOW_LISTS.load(Ordering::Relaxed),
            X_SHADOW_DIFFS.load(Ordering::Relaxed),
        );
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

    #[test]
    fn windows_close_at_query_boundaries_and_jobs_hold_up_to_four_matches_of_one_query() {
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
        let (end, matches, jobs) = x_window(&lists, 0, 3);
        assert_eq!(end, 3);
        let at = |q_idx, position| XMatch { q_idx, position };
        // The list without HSPs (query 0, position 1) is no match.
        assert_eq!(
            matches,
            [
                at(0, 0),
                at(0, 2),
                at(0, 3),
                at(0, 4),
                at(0, 5),
                at(0, 6),
                at(2, 0)
            ]
        );
        assert_eq!(jobs, [0..4, 4..6, 6..7]);
        // The last window takes the queries without matches after the last match.
        let (end, matches, jobs) = x_window(&lists, 3, 3);
        assert_eq!((end, matches.len()), (5, 2));
        assert_eq!(jobs, [0..2]);
        // A window holds at least one query.
        let (end, matches, jobs) = x_window(&lists, 1, 0);
        assert_eq!((end, matches.len(), jobs.len()), (2, 0, 0));
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
        let heaps = (0..2).map(|_| BlastCompoHeap::new(1, ECUTOFF)).collect();
        (lists, heaps, vec![0..2, 0..2])
    }

    #[test]
    fn the_replay_drops_what_the_early_termination_test_ends_and_reports_the_next_error() {
        let (lists, mut heaps, batches) = two_query_batch();
        let (_, matches, _) = x_window(&lists, 0, usize::MAX);
        assert_eq!(matches.len(), 4);
        // The match of best E 1.0 comes after both heaps are filled to the cutoff: the test
        // ends it, so its error is dropped. The next one (0.005 <= 5 * 0.002) is redone, and its
        // error is the one reported.
        let results = vec![
            redone(10, 1e-20),
            redone(20, 1e-20),
            Err(anyhow::anyhow!("ended")),
            Err(anyhow::anyhow!("reported")),
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
        let (lists, mut heaps, batches) = two_query_batch();
        let (_, matches, _) = x_window(&lists, 0, usize::MAX);
        let results = vec![
            redone(10, 1e-20),
            redone(20, 1e-20),
            Err(anyhow::anyhow!("ended")),
            redone(22, 0.005),
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
        assert_eq!(ended, 1);
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

    #[test]
    fn a_query_alone_in_its_batch_reads_only_its_own_heap() {
        // Query 0 has no match, so in a batch with query 1 the test is never true for query 1;
        // alone in its batch, query 1's second match is ended once its own heap is filled.
        let lists = vec![
            None,
            hit_list(vec![list(20, 1e-20, false), list(21, 1.0, false)]),
        ];
        for (batches, ended) in [(vec![0..2, 0..2], 0), (vec![0..1, 1..2], 1)] {
            let mut heaps: Vec<BlastCompoHeap> =
                (0..2).map(|_| BlastCompoHeap::new(1, ECUTOFF)).collect();
            let (_, matches, _) = x_window(&lists, 0, usize::MAX);
            let results = vec![redone(20, 1e-20), redone(21, 1.0)];
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
                ended,
                "{batches:?}"
            );
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
}
