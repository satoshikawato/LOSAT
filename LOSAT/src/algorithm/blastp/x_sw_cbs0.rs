//! BLASTP `-comp_based_stats 0 -use_sw_tback` (series RG, item X4-c0): NCBI's Smith-Waterman
//! redo without composition-based statistics, reading each query batch's preliminary HSP lists in
//! `(ordinal ASC, query_index DESC)` order.
//!
//! When composition-based statistics is disabled (`compositionBasedStats` is zero), NCBI does not
//! install `BlastHSPCBSStreamClose` (the score sort); `BlastHSPStreamClose` concatenates the
//! batch's lists in query order, sorts them by decreasing OID (preserving query order on equal
//! OIDs), and pops from the end. The resulting stream order within each query batch is therefore
//! ascending subject ordinal, breaking ties by descending query index, and
//! `Blast_RedoAlignmentCore_MT` still runs its early termination test, excluded-fetch status
//! reset, heap test before the first alignment, and batch status check on that stream with the
//! unadjusted scoring matrix (`localScalingFactor` of 1.0, no SEG masking).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:91-96,158-204,678-686
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1498
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3100-3109,3525-3575,3602-3646
//!
//! At more than one thread with `LOSAT_X_BPREDOPAR` on, each batch's OID-ordered event stream runs
//! in bounded windows on the search pool using LPT scheduling (`x_redo_par::x_compute_window_lpt`)
//! and replays early termination, deferred heap tests (`x_redo_par::x_sw_heap_test`), batch status,
//! and heap insertions in exact one-thread stream order. Default CBS 1–3 and ordinary traceback
//! paths do not enter this module.

use std::ops::Range;

use anyhow::Result;

use crate::algorithm::tblastx::lookup::QueryContext;
use crate::blastinput::seq_range::RangedBatch;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::{
    AdjustedProteinMatrix, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::smith_waterman_match::{
    unadjusted_square_matrix, SwBatchStatus,
};
use crate::stats::search_space::SearchSpace;
use crate::stats::{BlastGumbelBlk, KarlinParams};

use super::blast_engine::{
    blastp_use_matrix_workspace, preliminary_hit_from_local_match_hsp, sw_oid_stop_error,
    FetchExcludedMatch, FetchNotes,
};
use super::gapalign::GapAlignScratch;
use super::hsp::{blast_compo_early_termination, BlastCompoHeap, BlastpHitList, BlastpHspList};
use super::kappa::smith_waterman::{SwHeapAnswer, SwMatchScratch};
use super::kappa::{
    build_query_workspace, postprocess_preliminary_hits_smith_waterman, BlastRedoAlignParams,
    BlastpKappaQueryWorkspace, BlastpKappaSubjectRangeCache, EvalueAdjustContext,
};
use super::x_prelim_compact::{x_expand_hit_list, XPrelimHsp};
use super::x_redo_subjects::XRedoSubjects;
use crate::utils::threading::SearchPool;

/// Closure building per-query redo alignment parameters (`Fn(usize) -> Result<BlastRedoAlignParams> + Sync`).
pub(crate) trait RedoParamsFn: Fn(usize) -> Result<BlastRedoAlignParams> + Sync {}
impl<T: Fn(usize) -> Result<BlastRedoAlignParams> + Sync> RedoParamsFn for T {}

/// One entry in a query batch's OID-ordered preliminary stream (`-comp_based_stats 0 -use_sw_tback`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:91-96,158-204
#[derive(Clone, Copy, Debug, PartialEq)]
enum Cbs0Event {
    Kept {
        ordinal: u32,
        q_idx: usize,
        position: usize,
    },
    Excluded {
        ordinal: u32,
        q_idx: usize,
        best_evalue: f64,
    },
}

impl Cbs0Event {
    #[inline]
    fn ordinal(self) -> u32 {
        match self {
            Self::Kept { ordinal, .. } | Self::Excluded { ordinal, .. } => ordinal,
        }
    }

    #[inline]
    fn q_idx(self) -> usize {
        match self {
            Self::Kept { q_idx, .. } | Self::Excluded { q_idx, .. } => q_idx,
        }
    }
}

/// Parallel execution knobs for [`run`] (unused on non-threaded targets).
pub(crate) struct SwCbs0Par<'a> {
    pub pool: &'a SearchPool<'a>,
    pub threads: usize,
    pub allow_parallel: bool,
    pub timing: bool,
}

/// Inputs read by [`run`] for BLASTP `-comp_based_stats 0 -use_sw_tback`.
pub(crate) struct SwCbs0Search<'a, P> {
    pub par: SwCbs0Par<'a>,
    pub contexts: &'a [QueryContext],
    pub query_nomask: fn(&QueryContext) -> &[u8],
    pub search_spaces: &'a [SearchSpace],
    pub subjects: XRedoSubjects<'a>,
    /// Database ordinal ids by rank (`ListedSubjects::ordinals`), or `None` for `-subject`.
    pub ordinals: Option<&'a [u32]>,
    pub query_batches: &'a [RangedBatch],
    pub evalue_adjust_contexts: &'a [EvalueAdjustContext],
    pub matrix: ScoringMatrix,
    pub gapped_params: &'a KarlinParams,
    pub gapped_gumbel: Option<&'a BlastGumbelBlk>,
    pub evalue: f64,
    pub build_redo_align_params: &'a P,
    pub engine_positives: bool,
    pub fetch_excluded: &'a [FetchExcludedMatch],
}

/// Collects and sorts the kept matches and fetch-excluded marks of `batch` in NCBI's OID stream
/// order: `(ordinal ASC, query_index DESC)`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:91-96,158-204,321-329
fn collect_batch_events(
    batch: Range<usize>,
    lists: &[Option<BlastpHitList>],
    ordinals: Option<&[u32]>,
    fetch_excluded: &[FetchExcludedMatch],
) -> Vec<Cbs0Event> {
    let mut events = Vec::new();
    for q_idx in batch.clone() {
        if let Some(hit_list) = lists[q_idx].as_ref() {
            for (position, local_match) in hit_list
                .hsplist_array
                .iter()
                .take(hit_list.hsplist_count)
                .enumerate()
            {
                if local_match.hsps.is_empty() {
                    continue;
                }
                let ordinal =
                    ordinals.map_or(local_match.oid, |ords| ords[local_match.oid as usize]);
                events.push(Cbs0Event::Kept {
                    ordinal,
                    q_idx,
                    position,
                });
            }
        }
    }
    for mark in fetch_excluded {
        if batch.contains(&mark.query) {
            events.push(Cbs0Event::Excluded {
                ordinal: mark.ordinal,
                q_idx: mark.query,
                best_evalue: mark.best_evalue,
            });
        }
    }
    events.sort_unstable_by(|a, b| {
        a.ordinal()
            .cmp(&b.ordinal())
            .then_with(|| b.q_idx().cmp(&a.q_idx()))
    });
    events
}

fn query_workspace_for<P>(search: &SwCbs0Search<'_, P>, q_idx: usize) -> BlastpKappaQueryWorkspace {
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

#[inline]
fn local_match_at(
    lists: &[Option<BlastpHitList>],
    q_idx: usize,
    position: usize,
) -> &BlastpHspList {
    &lists[q_idx]
        .as_ref()
        .expect("a match's query has its preliminary hit list")
        .hsplist_array[position]
}

fn clear_batch_lists(
    lists: &mut [Option<BlastpHitList>],
    range: Range<usize>,
    x_prelim_freed: &mut usize,
) {
    for q_idx in range {
        if lists[q_idx].take().is_some() {
            *x_prelim_freed += 1;
        }
    }
}

#[inline]
fn record_excluded(
    sw_status: &mut SwBatchStatus,
    best_evalue: f64,
    batch_heaps: &[BlastCompoHeap],
) {
    let terminates_early = blast_compo_early_termination(best_evalue, batch_heaps);
    let redone = sw_status.reaches_redo(terminates_early, true);
    debug_assert!(!redone);
}

#[inline]
fn insert_postprocessed(
    heap: &mut BlastCompoHeap,
    oid: u32,
    postprocessed: super::kappa::BlastpPostprocessResult,
) {
    let Some(hsp_list) = postprocessed.hsp_list else {
        return;
    };
    if heap.would_insert(postprocessed.best_evalue, postprocessed.best_score, oid) {
        heap.insert(
            hsp_list,
            postprocessed.best_evalue,
            postprocessed.best_score,
            oid,
        );
    }
}

#[inline]
fn check_batch_end_status(
    sw_status: SwBatchStatus,
    fetch_notes: &FetchNotes<'_>,
    batch: Range<usize>,
    heaps: &[BlastCompoHeap],
) -> Result<()> {
    if sw_status.failed() {
        return Err(sw_oid_stop_error(fetch_notes, batch, heaps, None));
    }
    Ok(())
}

#[inline]
fn endless_oid_error(
    fetch_notes: &FetchNotes<'_>,
    batch: &Range<usize>,
    heaps: &[BlastCompoHeap],
    ordinal: u32,
    q_idx: usize,
) -> anyhow::Error {
    sw_oid_stop_error(fetch_notes, batch.clone(), heaps, Some((ordinal, q_idx)))
}

#[allow(clippy::too_many_arguments)]
fn run_sw_on_match<P>(
    search: &SwCbs0Search<'_, P>,
    q_idx: usize,
    local_match: &BlastpHspList,
    query_workspace: &BlastpKappaQueryWorkspace,
    redo_align_params: &BlastRedoAlignParams,
    unadjusted: &AdjustedProteinMatrix,
    scratch: &mut SwScratch,
    resume: Option<&super::kappa::smith_waterman::SwMatchResume>,
    would_insert: &mut dyn FnMut(f64, i32) -> SwHeapAnswer,
    save_resume: &mut Option<Box<super::kappa::smith_waterman::SwMatchResume>>,
) -> Result<super::kappa::BlastpPostprocessResult> {
    scratch.hits.clear();
    scratch.hits.extend(
        local_match
            .hsps
            .iter()
            .map(preliminary_hit_from_local_match_hsp),
    );
    let subject = search.subjects.read(local_match.oid)?;
    let subject_raw = subject.residues();
    postprocess_preliminary_hits_smith_waterman(
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
        Some(unadjusted),
        &mut scratch.sw_scratch,
        resume,
        would_insert,
        save_resume,
    )
}

/// Runs BLASTP's `-comp_based_stats 0 -use_sw_tback` redo across all query batches in NCBI's
/// OID stream order, expanding and releasing each batch's preliminary lists one batch at a time.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1498
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3100-3109,3525-3781
pub(crate) fn run<P: RedoParamsFn>(
    search: &SwCbs0Search<'_, P>,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    redone_matches: &mut [BlastCompoHeap],
    fetch_notes: &FetchNotes<'_>,
    x_prelim_freed: &mut usize,
) -> Result<()> {
    let unadjusted = unadjusted_square_matrix(search.matrix);
    let mut serial_scratch = SwScratch::new(search.matrix);

    for batch in search.query_batches {
        if batch.searched.is_empty() {
            continue;
        }
        let range = batch.searched.clone();
        // Expand only this batch's compact lists now, and free them as soon as the batch finishes.
        for q_idx in range.clone() {
            if let Some(compact) = x_compact_hit_lists[q_idx].take() {
                preliminary_hit_lists[q_idx] = Some(x_expand_hit_list(compact));
            }
        }
        let events = collect_batch_events(
            range.clone(),
            preliminary_hit_lists,
            search.ordinals,
            search.fetch_excluded,
        );
        let batch_result = if !events.is_empty() {
            // Validate initial matrix / parameters once at batch start, as `one_thread_redo` does.
            let _ = (search.build_redo_align_params)(range.start)?;
            run_batch(
                search,
                range.clone(),
                &events,
                preliminary_hit_lists,
                redone_matches,
                fetch_notes,
                &unadjusted,
                &mut serial_scratch,
            )
        } else {
            Ok(())
        };
        // Release this batch's preliminary lists even before propagating an error.
        clear_batch_lists(preliminary_hit_lists, range, x_prelim_freed);
        batch_result?;
    }
    Ok(())
}

struct SwScratch {
    composition_workspace: BlastCompositionWorkspace,
    gap_scratch: GapAlignScratch,
    range_cache: BlastpKappaSubjectRangeCache,
    sw_scratch: SwMatchScratch,
    hits: Vec<super::gapalign::BlastpPreliminaryHsp>,
}

impl SwScratch {
    fn new(matrix: ScoringMatrix) -> Self {
        let mut composition_workspace = BlastCompositionWorkspace::new_blosum62();
        if matrix != ScoringMatrix::Blosum62 {
            blastp_use_matrix_workspace(&mut composition_workspace, matrix);
        }
        Self {
            composition_workspace,
            gap_scratch: GapAlignScratch::new(),
            range_cache: BlastpKappaSubjectRangeCache::new(None),
            sw_scratch: SwMatchScratch::new(),
            hits: Vec::new(),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn run_batch<P: RedoParamsFn>(
    search: &SwCbs0Search<'_, P>,
    batch: Range<usize>,
    events: &[Cbs0Event],
    lists: &[Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    fetch_notes: &FetchNotes<'_>,
    unadjusted: &AdjustedProteinMatrix,
    serial_scratch: &mut SwScratch,
) -> Result<()> {
    if par::try_par(
        search,
        batch.clone(),
        events,
        lists,
        heaps,
        fetch_notes,
        unadjusted,
    )? {
        return Ok(());
    }
    run_batch_serial(
        search,
        batch,
        events,
        lists,
        heaps,
        fetch_notes,
        unadjusted,
        serial_scratch,
    )
}

/// One-thread OID-ordered Smith-Waterman redo of one query batch.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3781
#[allow(clippy::too_many_arguments)]
fn run_batch_serial<P>(
    search: &SwCbs0Search<'_, P>,
    batch: Range<usize>,
    events: &[Cbs0Event],
    lists: &[Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    fetch_notes: &FetchNotes<'_>,
    unadjusted: &AdjustedProteinMatrix,
    scratch: &mut SwScratch,
) -> Result<()>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams>,
{
    let mut workspaces: Vec<Option<BlastpKappaQueryWorkspace>> = vec![None; batch.len()];
    let mut params_by_q: Vec<Option<BlastRedoAlignParams>> =
        (0..batch.len()).map(|_| None).collect();
    let mut sw_status = SwBatchStatus::default();

    for &event in events {
        match event {
            Cbs0Event::Excluded { best_evalue, .. } => {
                record_excluded(&mut sw_status, best_evalue, &heaps[batch.clone()]);
            }
            Cbs0Event::Kept {
                ordinal,
                q_idx,
                position,
            } => {
                let m = local_match_at(lists, q_idx, position);
                if blast_compo_early_termination(m.best_evalue, &heaps[batch.clone()]) {
                    continue;
                }
                let slot = q_idx - batch.start;
                let query_workspace =
                    workspaces[slot].get_or_insert_with(|| query_workspace_for(search, q_idx));
                let redo_align_params = match params_by_q[slot].as_ref() {
                    Some(params) => params,
                    None => {
                        let built = (search.build_redo_align_params)(q_idx)?;
                        params_by_q[slot].insert(built)
                    }
                };
                let redone_match = &mut heaps[q_idx];
                let redone = run_sw_on_match(
                    search,
                    q_idx,
                    m,
                    query_workspace,
                    redo_align_params,
                    unadjusted,
                    scratch,
                    None,
                    &mut |evalue, score| {
                        SwHeapAnswer::Known(redone_match.would_insert(evalue, score, m.oid))
                    },
                    &mut None,
                );
                let Some(postprocessed) = sw_status.redone(redone, || {
                    endless_oid_error(fetch_notes, &batch, heaps, ordinal, q_idx)
                })?
                else {
                    continue;
                };
                insert_postprocessed(&mut heaps[q_idx], m.oid, postprocessed);
            }
        }
    }

    check_batch_end_status(sw_status, fetch_notes, batch, heaps)
}

// ---- Parallel windows for OID-ordered CBS=0 Smith-Waterman redo ----

#[cfg(any(
    not(feature = "parallel"),
    all(target_arch = "wasm32", not(feature = "wasm-threads"))
))]
mod par {
    use super::*;

    #[allow(clippy::too_many_arguments)]
    pub(super) fn try_par<P>(
        search: &SwCbs0Search<'_, P>,
        _batch: Range<usize>,
        _events: &[Cbs0Event],
        _lists: &[Option<BlastpHitList>],
        _heaps: &mut [BlastCompoHeap],
        _fetch_notes: &FetchNotes<'_>,
        _unadjusted: &AdjustedProteinMatrix,
    ) -> Result<bool> {
        let _ = (
            search.par.pool,
            search.par.threads,
            search.par.allow_parallel,
            search.par.timing,
        );
        Ok(false)
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod par {
    use super::super::hsp::COMPO_HEAP_EVALUE_STRETCH;
    use super::super::kappa::smith_waterman::SwMatchResume;
    use super::super::kappa::BlastpPostprocessResult;
    use super::super::x_redo_par::{
        x_bp_redo_lpt, x_compute_window_lpt, x_lpt_order, x_sw_heap_test, x_window_jobs_per_thread,
    };
    use super::*;
    use crate::core::composition_adjustment::smith_waterman_match::SwRedoFailed;
    use rayon::prelude::*;
    use std::sync::OnceLock;
    use std::time::Instant;

    /// Maximum matches grouped into one new parallel job.
    const OID_JOB_MATCHES: usize = 4;

    #[derive(Debug)]
    pub(super) struct OidComputed {
        pub result: Result<BlastpPostprocessResult>,
        pub deferred: Option<(f64, i32)>,
        pub resume: Option<Box<SwMatchResume>>,
    }

    #[derive(Debug)]
    pub(super) struct OidCarry {
        pub events: Vec<Cbs0Event>,
        pub results: Vec<Option<OidComputed>>,
    }

    #[derive(Debug, Default)]
    pub(super) struct OidCounts {
        pub matches: usize,
        pub windows: usize,
        pub skipped: usize,
        pub wasted: usize,
        pub turned: usize,
        pub continued: usize,
        pub carried_max: usize,
        pub par_s: f64,
        pub replay_s: f64,
    }

    #[derive(Clone, Copy, Debug)]
    pub(super) struct OidLimits {
        pub jobs: usize,
        pub at_risk_matches: usize,
    }

    #[derive(Debug)]
    pub(super) struct OidWindow {
        /// All events in this window (carried events followed by new events).
        pub events: Vec<Cbs0Event>,
        /// For each new kept event in `events[carried..]`, count of earlier non-ended kept matches
        /// of the same query in this window.
        pub before_same_query: Vec<usize>,
        /// Carried indices whose deferred test is now known TRUE (sorted by index).
        pub continued: Vec<(usize, (f64, i32))>,
        /// Number of carried events at the start of `events`.
        pub carried: usize,
        /// Jobs to execute in parallel: each entry is a list of event indices into `events`.
        /// Continued matches are 1-element jobs; new kept matches are grouped up to `OID_JOB_MATCHES`.
        pub jobs: Vec<Vec<usize>>,
        /// Next cursor in the batch's `events` slice.
        pub next_cursor: usize,
    }

    /// Computes the E-value above which early termination may trigger in `batch` given the current
    /// `heaps` and remaining kept match counts per query (`remaining_per_q`), or `None` if at least
    /// one query in `batch` cannot fill its heap anymore.
    fn batch_end_above(
        batch: Range<usize>,
        heaps: &[BlastCompoHeap],
        remaining_per_q: &[usize],
    ) -> Option<f64> {
        let may_fill = batch
            .clone()
            .all(|q| heaps[q].may_fill(remaining_per_q[q - batch.start]));
        may_fill.then(|| {
            batch
                .map(|q| COMPO_HEAP_EVALUE_STRETCH * heaps[q].ecutoff)
                .fold(f64::NEG_INFINITY, f64::max)
        })
    }

    /// Determines which carried deferred matches are known TRUE at their replay in OID order.
    ///
    /// Event 0 (the stopped match) is always continued. For any query other than the stopped
    /// match's query, a copy of its heap replays its carried matches ahead exactly until its first
    /// deferred TRUE or a match that early termination may end; after that (and for the stopped
    /// match's query) the heap-size bound (`heap.n + may_insert < heap.heap_threshold`) applies.
    pub(super) fn oid_sw_continued(
        carry: &OidCarry,
        batch: Range<usize>,
        lists: &[Option<BlastpHitList>],
        heaps: &[BlastCompoHeap],
        end_above: Option<f64>,
    ) -> Vec<(usize, (f64, i32))> {
        let mut continued = Vec::new();
        let nq = batch.len();
        let mut ahead: Vec<Option<BlastCompoHeap>> = vec![None; nq];
        let mut exact: Vec<bool> = vec![true; nq];
        let mut may_insert: Vec<usize> = vec![0; nq];

        for (index, (&event, computed)) in carry.events.iter().zip(&carry.results).enumerate() {
            let Cbs0Event::Kept {
                q_idx, position, ..
            } = event
            else {
                continue;
            };
            let Some(computed) = computed.as_ref() else {
                continue;
            };
            let slot = q_idx - batch.start;
            if index == 0 {
                let test = computed
                    .deferred
                    .expect("the stopped match's test was deferred");
                continued.push((index, test));
                exact[slot] = false;
                may_insert[slot] = 1;
                continue;
            }
            let m = local_match_at(lists, q_idx, position);
            let inserts = match &computed.result {
                Ok(postprocessed) => postprocessed.hsp_list.is_some(),
                Err(error) if error.is::<SwRedoFailed>() => false,
                Err(_) => break,
            };
            let may_end = end_above.is_some_and(|above| m.best_evalue > above);
            if exact[slot] && !may_end {
                let heap = ahead[slot].get_or_insert_with(|| heaps[q_idx].keys());
                match (computed.deferred, &computed.result) {
                    (Some((evalue, score)), _) => {
                        if heap.would_insert(evalue, score, m.oid) {
                            continued.push((index, (evalue, score)));
                            exact[slot] = false;
                            may_insert[slot] = 1;
                        }
                    }
                    (None, Ok(postprocessed)) if inserts => {
                        let (ev, sc) = (postprocessed.best_evalue, postprocessed.best_score);
                        if heap.would_insert(ev, sc, m.oid) {
                            let keys = BlastpHspList {
                                oid: m.oid,
                                query_index: m.query_index,
                                hsps: Vec::new(),
                                best_evalue: ev,
                            };
                            heap.insert(keys, ev, sc, m.oid);
                        }
                    }
                    (None, _) => {}
                }
                continue;
            }
            exact[slot] = false;
            let heap = ahead[slot].as_ref().unwrap_or(&heaps[q_idx]);
            match computed.deferred {
                Some(test)
                    if !may_end
                        && heap.n.saturating_add(may_insert[slot]) < heap.heap_threshold =>
                {
                    continued.push((index, test));
                    may_insert[slot] += 1;
                }
                Some(_) => may_insert[slot] += 1,
                None => may_insert[slot] += usize::from(inserts),
            }
        }
        continued
    }

    /// Builds the next window of carried + new OID-ordered events.
    #[allow(clippy::too_many_arguments)]
    fn build_oid_window(
        batch: Range<usize>,
        all_events: &[Cbs0Event],
        lists: &[Option<BlastpHitList>],
        heaps: &[BlastCompoHeap],
        end_above: Option<f64>,
        carry: Option<&OidCarry>,
        cursor: usize,
        limits: OidLimits,
        risk_cap: usize,
    ) -> OidWindow {
        let mut events = Vec::new();
        let mut before_same_query = Vec::new();
        let mut continued = Vec::new();
        let mut jobs: Vec<Vec<usize>> = Vec::new();
        let mut in_window_per_q = vec![0usize; batch.len()];

        let mut carried = 0usize;
        let mut carried_jobs = 0usize;
        if let Some(carry) = carry {
            carried = carry.events.len();
            events.extend_from_slice(&carry.events);
            before_same_query.resize(carried, 0);
            let mut carried_kept = 0usize;
            for (&ev, res) in carry.events.iter().zip(&carry.results) {
                if let (Cbs0Event::Kept { q_idx, .. }, Some(_)) = (ev, res) {
                    in_window_per_q[q_idx - batch.start] += 1;
                    carried_kept += 1;
                }
            }
            carried_jobs = carried_kept.div_ceil(OID_JOB_MATCHES);
            continued = oid_sw_continued(carry, batch.clone(), lists, heaps, end_above);
            for &(idx, _) in &continued {
                jobs.push(vec![idx]);
            }
        }

        let max_new_matches = limits.jobs.saturating_mul(OID_JOB_MATCHES);
        let max_total_new_jobs = (2 * limits.jobs)
            .saturating_sub(carried_jobs)
            .min(limits.jobs);
        let max_new_kept = max_total_new_jobs
            .saturating_mul(OID_JOB_MATCHES)
            .min(max_new_matches);

        let mut new_kept = 0usize;
        let mut new_at_risk = 0usize;
        let mut new_indices = Vec::new();
        let mut next_cursor = cursor;

        while next_cursor < all_events.len() {
            let event = all_events[next_cursor];
            match event {
                Cbs0Event::Excluded { .. } => {
                    events.push(event);
                    before_same_query.push(0);
                    next_cursor += 1;
                }
                Cbs0Event::Kept {
                    q_idx, position, ..
                } => {
                    let m = local_match_at(lists, q_idx, position);
                    let ended_now =
                        blast_compo_early_termination(m.best_evalue, &heaps[batch.clone()]);
                    let at_risk =
                        !ended_now && end_above.is_some_and(|above| m.best_evalue > above);
                    if (!jobs.is_empty() || !new_indices.is_empty())
                        && (new_kept >= max_new_kept || (at_risk && new_at_risk >= risk_cap))
                    {
                        break;
                    }
                    let idx = events.len();
                    events.push(event);
                    let slot = q_idx - batch.start;
                    before_same_query.push(in_window_per_q[slot]);
                    if !ended_now {
                        in_window_per_q[slot] += 1;
                        new_indices.push(idx);
                        new_kept += 1;
                        if at_risk {
                            new_at_risk += 1;
                        }
                    }
                    next_cursor += 1;
                }
            }
        }

        // Group new kept matches into jobs of up to OID_JOB_MATCHES (1 per job when small so all
        // threads stay busy on small batches).
        let chunk_size = if new_indices.len() <= limits.at_risk_matches {
            1
        } else {
            OID_JOB_MATCHES
        };
        for chunk in new_indices.chunks(chunk_size) {
            jobs.push(chunk.to_vec());
        }

        OidWindow {
            events,
            before_same_query,
            continued,
            carried,
            jobs,
            next_cursor,
        }
    }

    /// Replays a computed window in exact OID event order.
    #[allow(clippy::too_many_arguments)]
    fn replay_oid_window(
        batch: Range<usize>,
        events: &[Cbs0Event],
        results: Vec<Option<OidComputed>>,
        lists: &[Option<BlastpHitList>],
        heaps: &mut [BlastCompoHeap],
        remaining_per_q: &mut [usize],
        sw_status: &mut SwBatchStatus,
        fetch_notes: &FetchNotes<'_>,
        carry: &mut Option<OidCarry>,
    ) -> Result<(usize, usize, usize)> {
        assert_eq!(events.len(), results.len(), "one result per window event");
        let mut skipped = 0usize;
        let mut wasted = 0usize;
        let mut turned = 0usize;
        let mut results = results.into_iter();

        for (index, &event) in events.iter().enumerate() {
            let computed = results.next().expect("one result per event");
            match event {
                Cbs0Event::Excluded { best_evalue, .. } => {
                    record_excluded(sw_status, best_evalue, &heaps[batch.clone()]);
                }
                Cbs0Event::Kept {
                    ordinal,
                    q_idx,
                    position,
                } => {
                    let m = local_match_at(lists, q_idx, position);
                    let slot = q_idx - batch.start;
                    if blast_compo_early_termination(m.best_evalue, &heaps[batch.clone()]) {
                        if computed.is_some() {
                            wasted += 1;
                        } else {
                            skipped += 1;
                        }
                        remaining_per_q[slot] -= 1;
                        continue;
                    }
                    let Some(OidComputed {
                        result,
                        deferred,
                        resume,
                    }) = computed
                    else {
                        unreachable!(
                            "query {q_idx} ordinal {ordinal}: ended before its window, not ended at its replay"
                        );
                    };
                    if let Some((evalue, score)) = deferred {
                        if heaps[q_idx].would_insert(evalue, score, m.oid) {
                            let stopped = OidComputed {
                                result,
                                deferred,
                                resume,
                            };
                            *carry = Some(OidCarry {
                                events: events[index..].to_vec(),
                                results: std::iter::once(Some(stopped)).chain(results).collect(),
                            });
                            return Ok((skipped, wasted, turned));
                        }
                        turned += 1;
                    }
                    remaining_per_q[slot] -= 1;
                    let Some(postprocessed) = sw_status.redone(result, || {
                        endless_oid_error(fetch_notes, &batch, heaps, ordinal, q_idx)
                    })?
                    else {
                        continue;
                    };
                    insert_postprocessed(&mut heaps[q_idx], m.oid, postprocessed);
                }
            }
        }
        Ok((skipped, wasted, turned))
    }

    /// Generic OID-window driver shared between production and unit tests.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn drive_oid_windows(
        batch: Range<usize>,
        all_events: &[Cbs0Event],
        lists: &[Option<BlastpHitList>],
        heaps: &mut [BlastCompoHeap],
        limits: OidLimits,
        timing: bool,
        fetch_notes: &FetchNotes<'_>,
        compute: &mut dyn FnMut(
            &[BlastCompoHeap],
            &OidWindow,
            &[Option<OidComputed>],
        ) -> Vec<(usize, OidComputed)>,
    ) -> Result<OidCounts> {
        let mut remaining_per_q = vec![0usize; batch.len()];
        for &event in all_events {
            if let Cbs0Event::Kept { q_idx, .. } = event {
                remaining_per_q[q_idx - batch.start] += 1;
            }
        }
        let mut counts = OidCounts::default();
        let mut status = SwBatchStatus::default();
        let mut carry: Option<OidCarry> = None;
        let mut cursor = 0usize;
        let max_risk_cap = limits
            .jobs
            .saturating_mul(OID_JOB_MATCHES)
            .max(limits.at_risk_matches);
        let mut risk_cap = limits.at_risk_matches.max(1);

        while cursor < all_events.len() || carry.is_some() {
            let end_above = batch_end_above(batch.clone(), heaps, &remaining_per_q);
            let window = build_oid_window(
                batch.clone(),
                all_events,
                lists,
                heaps,
                end_above,
                carry.as_ref(),
                cursor,
                limits,
                risk_cap,
            );
            counts.windows += 1;
            counts.matches += window.events.len() - window.carried;

            let par_start = timing.then(Instant::now);
            let carried_slice = carry.as_ref().map_or(&[][..], |c| &c.results[..]);
            let computed_pairs = if window.jobs.is_empty() {
                Vec::new()
            } else {
                compute(heaps, &window, carried_slice)
            };
            if let Some(par_start) = par_start {
                counts.par_s += par_start.elapsed().as_secs_f64();
            }

            let mut results: Vec<Option<OidComputed>> = match carry.take() {
                Some(c) => {
                    counts.carried_max = counts.carried_max.max(window.carried);
                    let mut v = c.results;
                    v.resize_with(window.events.len(), || None);
                    v
                }
                None => (0..window.events.len()).map(|_| None).collect(),
            };
            for (idx, comp) in computed_pairs {
                if idx < window.carried {
                    counts.continued += 1;
                }
                results[idx] = Some(comp);
            }

            let replay_start = timing.then(Instant::now);
            let (skipped, wasted, turned) = replay_oid_window(
                batch.clone(),
                &window.events,
                results,
                lists,
                heaps,
                &mut remaining_per_q,
                &mut status,
                fetch_notes,
                &mut carry,
            )?;
            if let Some(replay_start) = replay_start {
                counts.replay_s += replay_start.elapsed().as_secs_f64();
            }
            counts.skipped += skipped;
            counts.wasted += wasted;
            counts.turned += turned;
            if wasted > 0 {
                risk_cap = limits.at_risk_matches.max(1);
            } else if carry.is_none() {
                risk_cap = risk_cap.saturating_mul(2).min(max_risk_cap);
            }
            cursor = window.next_cursor;
        }

        check_batch_end_status(status, fetch_notes, batch, heaps)?;
        Ok(counts)
    }

    #[allow(clippy::too_many_arguments)]
    fn redo_one_oid_match<P>(
        search: &SwCbs0Search<'_, P>,
        q_idx: usize,
        m: &BlastpHspList,
        query_workspace: &BlastpKappaQueryWorkspace,
        redo_align_params: &BlastRedoAlignParams,
        unadjusted: &AdjustedProteinMatrix,
        scratch: &mut SwScratch,
        heap_before: Option<(&BlastCompoHeap, usize)>,
        resume: Option<&SwMatchResume>,
    ) -> OidComputed {
        let mut deferred = None;
        let mut kept = None;
        let result = run_sw_on_match(
            search,
            q_idx,
            m,
            query_workspace,
            redo_align_params,
            unadjusted,
            scratch,
            resume,
            &mut |evalue, score| match heap_before {
                Some((heap, before)) => {
                    let inserts =
                        x_sw_heap_test(heap, before, (evalue, score, m.oid), &mut deferred);
                    if deferred.is_some() {
                        SwHeapAnswer::Deferred
                    } else {
                        SwHeapAnswer::Known(inserts)
                    }
                }
                None => unreachable!("a continuation's heap test is asked again"),
            },
            &mut kept,
        );
        debug_assert_eq!(
            deferred.is_some(),
            kept.is_some(),
            "a deferred test keeps the redo's state"
        );
        OidComputed {
            result,
            deferred,
            resume: kept,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn try_par<P: RedoParamsFn>(
        search: &SwCbs0Search<'_, P>,
        batch: Range<usize>,
        events: &[Cbs0Event],
        lists: &[Option<BlastpHitList>],
        heaps: &mut [BlastCompoHeap],
        fetch_notes: &FetchNotes<'_>,
        unadjusted: &AdjustedProteinMatrix,
    ) -> Result<bool> {
        if !search.par.allow_parallel || !search.par.pool.enabled() || search.par.threads <= 1 {
            return Ok(false);
        }
        let kept = events
            .iter()
            .filter(|e| matches!(e, Cbs0Event::Kept { .. }))
            .count();
        if kept < 2 {
            return Ok(false);
        }
        let workspaces: Vec<OnceLock<BlastpKappaQueryWorkspace>> =
            (0..batch.len()).map(|_| OnceLock::new()).collect();
        let threads = search.par.threads.max(1);
        let limits = OidLimits {
            jobs: x_window_jobs_per_thread().saturating_mul(threads),
            at_risk_matches: threads.saturating_mul(OID_JOB_MATCHES),
        };
        let lpt = x_bp_redo_lpt();
        let counts = drive_oid_windows(
            batch.clone(),
            events,
            lists,
            heaps,
            limits,
            search.par.timing,
            fetch_notes,
            &mut |heaps, window, carried| {
                let run_job = |scratch: &mut SwScratch, job_indices: &[usize]| {
                    let mut out = Vec::with_capacity(job_indices.len());
                    let mut last_q: Option<(usize, BlastRedoAlignParams)> = None;
                    for &idx in job_indices {
                        let Cbs0Event::Kept {
                            q_idx, position, ..
                        } = window.events[idx]
                        else {
                            continue;
                        };
                        let m = local_match_at(lists, q_idx, position);
                        let ws = workspaces[q_idx - batch.start]
                            .get_or_init(|| query_workspace_for(search, q_idx));
                        if last_q.as_ref().is_none_or(|(q, _)| *q != q_idx) {
                            match (search.build_redo_align_params)(q_idx) {
                                Ok(params) => last_q = Some((q_idx, params)),
                                Err(error) => {
                                    out.push((
                                        idx,
                                        OidComputed {
                                            result: Err(anyhow::anyhow!("{error:#}")),
                                            deferred: None,
                                            resume: None,
                                        },
                                    ));
                                    continue;
                                }
                            }
                        }
                        let redo_align_params = &last_q.as_ref().expect("built above").1;
                        let (heap_before, resume) = match window
                            .continued
                            .binary_search_by_key(&idx, |&(at, _)| at)
                        {
                            Ok(found) => {
                                let stopped = carried[idx].as_ref();
                                debug_assert_eq!(
                                    stopped.and_then(|s| s.deferred),
                                    Some(window.continued[found].1),
                                    "a continuation resumes from the test that was continued"
                                );
                                (None, stopped.and_then(|s| s.resume.as_deref()))
                            }
                            Err(_) => (Some((&heaps[q_idx], window.before_same_query[idx])), None),
                        };
                        let comp = redo_one_oid_match(
                            search,
                            q_idx,
                            m,
                            ws,
                            redo_align_params,
                            unadjusted,
                            scratch,
                            heap_before,
                            resume,
                        );
                        out.push((idx, comp));
                    }
                    out
                };
                let new_scratch = || SwScratch::new(search.matrix);
                let per_job: Vec<Vec<(usize, OidComputed)>> =
                    if lpt && search.par.pool.enabled() && threads > 1 && window.jobs.len() > 1 {
                        let costs: Vec<u64> = window
                            .jobs
                            .iter()
                            .map(|job| {
                                job.iter()
                                    .map(|&idx| {
                                        if let Cbs0Event::Kept {
                                            q_idx, position, ..
                                        } = window.events[idx]
                                        {
                                            let q_len = search.contexts[q_idx].aa_len as u64;
                                            let s_len = search
                                                .subjects
                                                .aa_len(
                                                    local_match_at(lists, q_idx, position).oid
                                                        as usize,
                                                )
                                                .map_or(0, |l| l as u64);
                                            q_len.saturating_mul(s_len)
                                        } else {
                                            0
                                        }
                                    })
                                    .fold(0, u64::saturating_add)
                            })
                            .collect();
                        let order = x_lpt_order(&costs);
                        x_compute_window_lpt(
                            search.par.pool,
                            threads.min(window.jobs.len()),
                            &order,
                            new_scratch,
                            |scratch, job_idx| run_job(scratch, &window.jobs[job_idx]),
                        )
                    } else {
                        search.par.pool.install(|| {
                            window
                                .jobs
                                .par_iter()
                                .map_init(new_scratch, |scratch, job| run_job(scratch, job))
                                .collect()
                        })
                    };
                per_job.into_iter().flatten().collect()
            },
        )?;
        if search.par.timing {
            eprintln!(
                "[TIMING] x_sw_cbs0: batch={}..{} matches={} windows={} skipped={} wasted={} turned={} continued={} carried_max={} par_s={:.6} replay_s={:.6}",
                batch.start,
                batch.end,
                counts.matches,
                counts.windows,
                counts.skipped,
                counts.wasted,
                counts.turned,
                counts.continued,
                counts.carried_max,
                counts.par_s,
                counts.replay_s,
            );
        }
        Ok(true)
    }
}

#[cfg(test)]
#[path = "x_sw_cbs0_tests.rs"]
mod tests;
