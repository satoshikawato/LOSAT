//! BLASTP's Smith-Waterman redo without composition-based statistics (`-comp_based_stats 0
//! -use_sw_tback`, line RG item X4-c0).
//!
//! Without composition-based statistics, NCBI's `Blast_RedoAlignmentCore_MT` still runs for
//! `-use_sw_tback` (`blast_traceback.c:1486-1487`), with two differences from modes 1 to 3:
//!
//! 1. Each query batch's preliminary stream is not sorted by score (`blast_hspstream.c:678-686`):
//!    `BlastHSPStreamClose` collects all queries' lists in the batch, sorts by decreasing OID
//!    (stable over ascending query index on glibc), and `BlastHSPStreamRead` pops from the end,
//!    yielding `(ordinal ASC, query DESC)` across all queries of the batch.
//! 2. No composition adjustment or subject SEG runs (a scaling factor of 1.0 and rule
//!    `eDontAdjustMatrix`), and every window scores against the unadjusted search matrix
//!    (`smith_waterman_match::unadjusted_square_matrix`).
//!
//! Both the one-thread loop (`sw_c0_one_thread`) and the windowed parallel path
//! (`sw_c0_windows`, under `LOSAT_X_BPREDOPAR` at 2+ threads with 2+ matches) reuse core's
//! `smith_waterman_match` (via `kappa::postprocess_preliminary_hits_smith_waterman`),
//! `SwBatchStatus`, `redo_run::sw_batch_failed_after_oid` / `sw_endless_after_oid`, and
//! `x_redo_par`'s `x_sw_heap_test` / `x_compute_window_lpt`. Because queries are interleaved in
//! OID order and preliminary E-values are not monotone, `sw_c0_windows` tracks per-query `before`
//! counts across each window, interleaves `FetchExcludedMatch` marks at their exact database
//! ordinals, speculates `Known(true)` when `x_sw_heap_test` defers (still short-circuiting
//! `Known(false)` once heaps fill via a compact fill window for small hit lists), and replays all
//! order-dependent steps (early termination, deferred heap tests, batch status, heap insertions)
//! in exact `(batch, ordinal ASC, query DESC)` order.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1487
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:91-96,158-207,678-686
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3100-3109,3525-3781

use std::cmp::Reverse;
use std::ops::Range;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use {std::collections::VecDeque, std::sync::OnceLock, std::time::Instant};

use anyhow::Result;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

use super::database::FetchExcludedMatch;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use super::database::FetchNotes;
use super::redo_run::{sw_batch_failed_after_oid, sw_endless_after_oid, RedoInputs};
use super::{preliminary_hit_from_local_match_hsp, query_nomask_sequence};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::gapalign::GapAlignScratch;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastp::hsp::COMPO_HEAP_EVALUE_STRETCH;
use crate::algorithm::blastp::hsp::{
    blast_compo_early_termination, BlastCompoHeap, BlastpHitList, BlastpHspList,
};
use crate::algorithm::blastp::kappa::smith_waterman::{SwHeapAnswer, SwMatchScratch};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastp::kappa::BlastpPostprocessResult;
use crate::algorithm::blastp::kappa::{
    self, build_query_workspace, BlastRedoAlignParams, BlastpKappaQueryWorkspace,
    BlastpKappaSubjectRangeCache,
};
use crate::algorithm::blastp::x_prelim_compact::{x_expand_hit_list, XPrelimHsp};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastp::x_redo_par::{
    self, x_compute_window_lpt, x_lpt_order, x_sw_heap_test, x_window_jobs_per_thread,
    X_JOB_MATCHES,
};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::core::composition_adjustment::adjust_scores::AdjustedProteinMatrix;
use crate::core::composition_adjustment::adjust_scores::BlastCompositionWorkspace;
use crate::core::composition_adjustment::smith_waterman_match::{
    unadjusted_square_matrix, SwBatchStatus,
};

/// One item in a batch's OID-ordered stream for `-comp_based_stats 0 -use_sw_tback`.
#[derive(Clone, Copy, Debug, PartialEq)]
enum SwC0Kind {
    /// A kept preliminary HSP list at `preliminary_hit_lists[q_idx].hsplist_array[position]`.
    Kept { position: usize },
    /// A preliminary HSP list whose subject NCBI's fetch excluded (`database::exclude_fetched`).
    Excluded { best_evalue: f64 },
}

/// An entry of a query batch's OID-ordered stream (`(batch_idx, ordinal ASC, Reverse(q_idx))`).
#[derive(Clone, Copy, Debug, PartialEq)]
struct SwC0Entry {
    batch_idx: usize,
    ordinal: u32,
    q_idx: usize,
    kind: SwC0Kind,
}

/// Disjoint contiguous query batches in order from `kappa_query_batches`.
fn batch_ranges(kappa_query_batches: &[Range<usize>]) -> Vec<Range<usize>> {
    let mut batches = Vec::new();
    let mut q = 0;
    while q < kappa_query_batches.len() {
        let batch = kappa_query_batches[q].clone();
        debug_assert!(batch.start <= q && batch.end > q);
        let end = batch.end;
        batches.push(batch);
        q = end;
    }
    batches
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:91-96,158-207,678-686
/// Builds one query batch's OID-ordered stream, merging kept preliminary lists and excluded fetch
/// marks in NCBI's `(ordinal ASC, query DESC)` order.
fn build_batch_stream(
    batch_idx: usize,
    batch: &Range<usize>,
    lists: &[Option<BlastpHitList>],
    excluded: &[FetchExcludedMatch],
    listed_ordinals: Option<&[u32]>,
) -> Vec<SwC0Entry> {
    let mut entries = Vec::new();
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
                    listed_ordinals.map_or(local_match.oid, |ords| ords[local_match.oid as usize]);
                entries.push(SwC0Entry {
                    batch_idx,
                    ordinal,
                    q_idx,
                    kind: SwC0Kind::Kept { position },
                });
            }
        }
    }
    let start = excluded.partition_point(|mark| mark.query < batch.start);
    for mark in &excluded[start..] {
        if mark.query >= batch.end {
            break;
        }
        entries.push(SwC0Entry {
            batch_idx,
            ordinal: mark.ordinal,
            q_idx: mark.query,
            kind: SwC0Kind::Excluded {
                best_evalue: mark.best_evalue,
            },
        });
    }
    entries.sort_by_key(|entry| (entry.ordinal, Reverse(entry.q_idx)));
    entries
}

/// Expands every compact preliminary hit list of `batch` into `preliminary_hit_lists`.
fn expand_batch(
    batch: &Range<usize>,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) {
    for q_idx in batch.clone() {
        if let Some(compact) = x_compact_hit_lists[q_idx].take() {
            preliminary_hit_lists[q_idx] = Some(x_expand_hit_list(compact));
        }
    }
}

/// Releases every expanded preliminary hit list of `batch`, returning how many were freed.
fn release_batch(
    batch: &Range<usize>,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
) -> usize {
    let mut released = 0;
    for q_idx in batch.clone() {
        if preliminary_hit_lists[q_idx].take().is_some() {
            released += 1;
        }
    }
    released
}

fn local_match_of(
    lists: &[Option<BlastpHitList>],
    q_idx: usize,
    position: usize,
) -> &BlastpHspList {
    &lists[q_idx]
        .as_ref()
        .expect("a kept stream entry's query has its preliminary hit list")
        .hsplist_array[position]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1487
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:678-686
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3100-3109,3525-3781
/// One-thread Smith-Waterman redo without composition-based statistics (`-comp_based_stats 0
/// -use_sw_tback`): processes each query batch in `(ordinal ASC, query DESC)` order with the
/// unadjusted search matrix, expanding and releasing one batch at a time.
pub(super) fn sw_c0_one_thread<P>(
    redo: &RedoInputs<'_, P>,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    redone_matches: &mut [BlastCompoHeap],
    x_prelim_freed: &mut usize,
) -> Result<()>
where
    P: Fn(usize) -> Result<BlastRedoAlignParams>,
{
    let RedoInputs {
        args,
        needs,
        subjects,
        contexts,
        search_spaces,
        gapped_params,
        gapped_gumbel,
        kappa_query_batches,
        evalue_adjust_contexts,
        build_redo_align_params,
        fetch_excluded,
        fetch_notes,
        listed_ordinals,
        ..
    } = *redo;
    if kappa_query_batches.is_empty() {
        return Ok(());
    }
    let mut redo_align_params = build_redo_align_params(0)?;
    let unadjusted = unadjusted_square_matrix(args.scoring.matrix);
    let mut composition_workspace = BlastCompositionWorkspace::new_blosum62();
    let mut kappa_preliminary_hits = Vec::new();
    let mut kappa_gap_scratch = GapAlignScratch::new();
    let mut kappa_subject_range_cache = BlastpKappaSubjectRangeCache::new(None);
    let mut sw_scratch = SwMatchScratch::new();
    let batches = batch_ranges(kappa_query_batches);

    for (batch_idx, batch) in batches.iter().enumerate() {
        expand_batch(batch, preliminary_hit_lists, x_compact_hit_lists);
        let stream = build_batch_stream(
            batch_idx,
            batch,
            preliminary_hit_lists,
            fetch_excluded,
            listed_ordinals,
        );
        let mut workspaces: Vec<Option<BlastpKappaQueryWorkspace>> =
            std::iter::repeat_with(|| None).take(batch.len()).collect();
        let mut sw_status = SwBatchStatus::default();
        for entry in stream {
            match entry.kind {
                SwC0Kind::Excluded { best_evalue } => {
                    let terminates_early =
                        blast_compo_early_termination(best_evalue, &redone_matches[batch.clone()]);
                    let redone = sw_status.reaches_redo(terminates_early, true);
                    debug_assert!(!redone);
                }
                SwC0Kind::Kept { position } => {
                    let local_match = local_match_of(preliminary_hit_lists, entry.q_idx, position);
                    if blast_compo_early_termination(
                        local_match.best_evalue,
                        &redone_matches[batch.clone()],
                    ) {
                        continue;
                    }
                    kappa_preliminary_hits.clear();
                    kappa_preliminary_hits.extend(
                        local_match
                            .hsps
                            .iter()
                            .map(preliminary_hit_from_local_match_hsp),
                    );
                    let s_idx = local_match.oid;
                    let subject = subjects.read(s_idx)?;
                    let subject_raw = subject.residues();
                    let ctx = &contexts[entry.q_idx];
                    let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
                    let query_nomask = query_nomask_sequence(ctx);
                    let query_workspace =
                        workspaces[entry.q_idx - batch.start].get_or_insert_with(|| {
                            build_query_workspace(
                                query_raw,
                                ctx.aa_len as i32,
                                search_spaces[entry.q_idx].length_adjustment as i32,
                                search_spaces[entry.q_idx].effective_space,
                                evalue_adjust_contexts[entry.q_idx],
                            )
                        });
                    redo_align_params.ccat_query_length = ctx.aa_len as i32;
                    let redone_match = &mut redone_matches[entry.q_idx];
                    let redone = kappa::postprocess_preliminary_hits_smith_waterman(
                        &kappa_preliminary_hits,
                        query_workspace,
                        &mut composition_workspace,
                        query_nomask,
                        subject_raw,
                        args.scoring.matrix,
                        Some(&unadjusted),
                        gapped_params,
                        gapped_gumbel.as_ref(),
                        &redo_align_params,
                        args.evalue,
                        &mut kappa_gap_scratch,
                        &mut kappa_subject_range_cache,
                        needs.seqalign,
                        &mut sw_scratch,
                        None,
                        &mut |evalue, score| {
                            SwHeapAnswer::Known(redone_match.would_insert(evalue, score, s_idx))
                        },
                        &mut None,
                    );
                    let Some(postprocessed) = sw_status.redone(redone, || {
                        sw_endless_after_oid(
                            fetch_notes,
                            batch.clone(),
                            redone_matches,
                            entry.ordinal,
                            entry.q_idx,
                        )
                    })?
                    else {
                        continue;
                    };
                    let best_score = postprocessed.best_score;
                    let best_evalue = postprocessed.best_evalue;
                    let Some(hsp_list) = postprocessed.hsp_list else {
                        continue;
                    };
                    let redone_match = &mut redone_matches[entry.q_idx];
                    if redone_match.would_insert(best_evalue, best_score, s_idx) {
                        redone_match.insert(hsp_list, best_evalue, best_score, s_idx);
                    }
                }
            }
        }
        *x_prelim_freed += release_batch(batch, preliminary_hit_lists);
        if sw_status.failed() {
            return Err(sw_batch_failed_after_oid(
                fetch_notes,
                batch.clone(),
                redone_matches,
            ));
        }
    }
    Ok(())
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod par {
    use super::*;

    /// Bounds of a CBS 0 Smith-Waterman window: normal window size, compact fill window for small
    /// unfilled heaps, compact steady-state window for small heaps, and initial at-risk cap per
    /// query batch.
    #[derive(Clone, Copy, Debug)]
    pub(super) struct SwC0WindowLimits {
        pub(super) matches: usize,
        pub(super) fill_matches: usize,
        pub(super) small_heap_matches: usize,
        pub(super) at_risk: usize,
    }

    /// One computed kept match: its redo result and, when its heap test could not be decided from
    /// the window-start heap (`x_sw_heap_test`), the first-round `(evalue, score)` test inputs.
    /// Workers speculate `Known(true)` on deferred tests so OID-order replay never stalls on
    /// improving matches; replay checks `heaps[q].would_insert` in exact stream order and discards
    /// speculative results (including failures or endless-redo errors) when FALSE.
    #[derive(Debug)]
    pub(super) struct SwC0Computed {
        pub(super) result: Result<BlastpPostprocessResult>,
        pub(super) deferred: Option<(f64, i32)>,
    }

    /// One window of OID-ordered stream entries, computed on the pool and replayed in stream order.
    #[derive(Debug)]
    pub(super) struct SwC0Window {
        pub(super) entries: Vec<SwC0Entry>,
        pub(super) ended: Vec<bool>,
        /// Entry indices in `entries` to compute on the pool.
        pub(super) jobs: Vec<usize>,
        /// Per entry `i`: non-ended kept matches of `entries[i].q_idx` before `i` in this window.
        pub(super) before: Vec<usize>,
        /// Per-batch count of at-risk kept matches in this window.
        at_risk: Vec<(usize, usize)>,
    }

    #[derive(Debug, Default)]
    pub(super) struct SwC0Counts {
        pub(super) matches: usize,
        pub(super) jobs: usize,
        pub(super) windows: usize,
        pub(super) skipped: usize,
        pub(super) wasted: usize,
        pub(super) released: usize,
        pub(super) par_s: f64,
        pub(super) replay_s: f64,
        pub(super) lpt_windows: Option<usize>,
        pub(super) turned: usize,
        pub(super) continued: usize,
    }

    #[derive(Debug, Default)]
    struct SwC0Replayed {
        skipped: usize,
        wasted: usize,
        turned: usize,
        continued: usize,
        dropped_batches: Vec<usize>,
    }

    impl SwC0Replayed {
        fn drop_batch(&mut self, batch_idx: usize) {
            if self.dropped_batches.last() != Some(&batch_idx) {
                self.dropped_batches.push(batch_idx);
            }
        }
    }

    /// Batch status carried across windows during replay.
    #[derive(Default)]
    struct SwC0BatchState {
        current: Option<usize>,
        status: SwBatchStatus,
    }

    impl SwC0BatchState {
        fn finish_if_needed(
            &mut self,
            next_batch: Option<usize>,
            batches: &[Range<usize>],
            heaps: &[BlastCompoHeap],
            notes: &FetchNotes<'_>,
        ) -> Result<()> {
            if self.current == next_batch {
                return Ok(());
            }
            if let Some(prev) = self.current {
                if self.status.failed() {
                    return Err(sw_batch_failed_after_oid(
                        notes,
                        batches[prev].clone(),
                        heaps,
                    ));
                }
            }
            self.current = next_batch;
            self.status = SwBatchStatus::default();
            Ok(())
        }
    }

    /// Batch-ahead stream queue: expands one batch's compact lists at a time into `pending`.
    struct SwC0Queue<'a> {
        batches: &'a [Range<usize>],
        excluded: &'a [FetchExcludedMatch],
        listed_ordinals: Option<&'a [u32]>,
        expanded_batch: usize,
        pending: VecDeque<SwC0Entry>,
        pending_kept: usize,
    }

    impl<'a> SwC0Queue<'a> {
        fn new(
            batches: &'a [Range<usize>],
            excluded: &'a [FetchExcludedMatch],
            listed_ordinals: Option<&'a [u32]>,
        ) -> Self {
            Self {
                batches,
                excluded,
                listed_ordinals,
                expanded_batch: 0,
                pending: VecDeque::new(),
                pending_kept: 0,
            }
        }

        fn ensure_kept(
            &mut self,
            needed: usize,
            lists: &mut [Option<BlastpHitList>],
            compact: &mut [Option<BlastpHitList<XPrelimHsp>>],
        ) {
            while self.pending_kept < needed && self.expanded_batch < self.batches.len() {
                let batch_idx = self.expanded_batch;
                let batch = &self.batches[batch_idx];
                expand_batch(batch, lists, compact);
                let entries = build_batch_stream(
                    batch_idx,
                    batch,
                    lists,
                    self.excluded,
                    self.listed_ordinals,
                );
                self.pending_kept += entries
                    .iter()
                    .filter(|entry| matches!(entry.kind, SwC0Kind::Kept { .. }))
                    .count();
                self.pending.extend(entries);
                self.expanded_batch += 1;
            }
        }

        fn is_done(&self) -> bool {
            self.pending.is_empty() && self.expanded_batch == self.batches.len()
        }
    }

    fn count_list_matches<H>(list: &BlastpHitList<H>) -> usize {
        list.hsplist_array
            .iter()
            .take(list.hsplist_count)
            .filter(|local_match| !local_match.hsps.is_empty())
            .count()
    }

    fn query_remaining_matches(
        lists: &[Option<BlastpHitList>],
        compact: &[Option<BlastpHitList<XPrelimHsp>>],
    ) -> Vec<usize> {
        (0..lists.len())
            .map(|q_idx| match compact.get(q_idx).and_then(Option::as_ref) {
                Some(list) => count_list_matches(list),
                None => lists[q_idx].as_ref().map_or(0, count_list_matches),
            })
            .collect()
    }

    /// For each batch, the E-value threshold above which `BlastCompo_EarlyTermination` can still
    /// become true during replay given `heaps` and each query's `remaining` un-replayed matches.
    fn batch_end_above(
        batches: &[Range<usize>],
        heaps: &[BlastCompoHeap],
        remaining: &[usize],
    ) -> Vec<Option<f64>> {
        batches
            .iter()
            .map(|batch| {
                let may_end = batch
                    .clone()
                    .all(|q| heaps[q].x_ends_early && heaps[q].may_fill(remaining[q]));
                may_end.then(|| {
                    batch
                        .clone()
                        .map(|q| COMPO_HEAP_EVALUE_STRETCH * heaps[q].ecutoff)
                        .fold(f64::NEG_INFINITY, f64::max)
                })
            })
            .collect()
    }

    /// Chooses the target match count for the next window based on the active batch's hit list
    /// size: small hit lists (`< 64`) fill empty heaps with a compact window first so subsequent
    /// windows immediately prune weak matches via `x_sw_heap_test` (`Known(false)`).
    fn window_match_budget(
        batches: &[Range<usize>],
        front_batch: Option<usize>,
        heaps: &[BlastCompoHeap],
        remaining: &[usize],
        limits: SwC0WindowLimits,
    ) -> usize {
        let Some(batch_idx) = front_batch else {
            return limits.matches;
        };
        let batch = &batches[batch_idx];
        let threshold = batch
            .clone()
            .map(|q| heaps[q].heap_threshold)
            .min()
            .unwrap_or(usize::MAX);
        if threshold < 64 {
            let needs_fill = batch
                .clone()
                .any(|q| heaps[q].n < heaps[q].heap_threshold && remaining[q] > 0);
            if needs_fill {
                limits.fill_matches
            } else {
                limits.small_heap_matches
            }
        } else {
            limits.matches
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn build_window(
        queue: &mut SwC0Queue<'_>,
        lists: &mut [Option<BlastpHitList>],
        compact: &mut [Option<BlastpHitList<XPrelimHsp>>],
        heaps: &[BlastCompoHeap],
        remaining: &[usize],
        end_above: &[Option<f64>],
        caps: &[usize],
        limits: SwC0WindowLimits,
    ) -> SwC0Window {
        let mut window = SwC0Window {
            entries: Vec::new(),
            ended: Vec::new(),
            jobs: Vec::new(),
            before: Vec::new(),
            at_risk: Vec::new(),
        };
        if queue.pending.is_empty() {
            queue.ensure_kept(1, lists, compact);
        }
        let max_new = window_match_budget(
            queue.batches,
            queue.pending.front().map(|entry| entry.batch_idx),
            heaps,
            remaining,
            limits,
        );
        if max_new > 0 {
            queue.ensure_kept(max_new, lists, compact);
        }
        let mut in_window_by_q = vec![0usize; heaps.len()];
        let mut new_kept = 0usize;
        while let Some(&front) = queue.pending.front() {
            match front.kind {
                SwC0Kind::Excluded { .. } => {
                    if !window.entries.is_empty() && new_kept >= max_new {
                        break;
                    }
                    queue.pending.pop_front();
                    window.entries.push(front);
                    window.ended.push(false);
                }
                SwC0Kind::Kept { position } => {
                    if !window.entries.is_empty() && new_kept >= max_new {
                        break;
                    }
                    let batch = &queue.batches[front.batch_idx];
                    let best_evalue = local_match_of(lists, front.q_idx, position).best_evalue;
                    let may_end = end_above[front.batch_idx].is_some_and(|above| {
                        best_evalue > above
                            && batch.clone().all(|q| heaps[q].may_fill(in_window_by_q[q]))
                    });
                    let ended = may_end
                        && blast_compo_early_termination(best_evalue, &heaps[batch.clone()]);
                    let at_risk = may_end && !ended;
                    let current_at_risk = window
                        .at_risk
                        .iter()
                        .find_map(|&(b, count)| (b == front.batch_idx).then_some(count))
                        .unwrap_or(0);
                    if !window.jobs.is_empty()
                        && at_risk
                        && current_at_risk >= caps[front.batch_idx]
                    {
                        break;
                    }
                    queue.pending.pop_front();
                    queue.pending_kept -= 1;
                    let at = window.entries.len();
                    window.entries.push(front);
                    window.ended.push(ended);
                    new_kept += 1;
                    if !ended {
                        window.jobs.push(at);
                        in_window_by_q[front.q_idx] += 1;
                        if at_risk {
                            match window
                                .at_risk
                                .iter_mut()
                                .find(|(b, _)| *b == front.batch_idx)
                            {
                                Some((_, count)) => *count += 1,
                                None => window.at_risk.push((front.batch_idx, 1)),
                            }
                        }
                    }
                }
            }
        }
        // Compute `before[i]`: earlier non-ended kept matches of `q_idx` in `window.entries[..i]`.
        window.before.resize(window.entries.len(), 0);
        let mut seen_by_q = vec![0usize; heaps.len()];
        for (i, (&entry, &ended)) in window.entries.iter().zip(&window.ended).enumerate() {
            if matches!(entry.kind, SwC0Kind::Kept { .. }) && !ended {
                window.before[i] = seen_by_q[entry.q_idx];
                seen_by_q[entry.q_idx] += 1;
            }
        }
        window
    }

    #[allow(clippy::too_many_arguments)]
    fn sw_c0_replay(
        entries: &[SwC0Entry],
        results: Vec<Option<SwC0Computed>>,
        lists: &[Option<BlastpHitList>],
        heaps: &mut [BlastCompoHeap],
        batches: &[Range<usize>],
        remaining_by_q: &mut [usize],
        batch_state: &mut SwC0BatchState,
        notes: &FetchNotes<'_>,
    ) -> Result<SwC0Replayed> {
        assert_eq!(entries.len(), results.len(), "one result slot per entry");
        let mut replayed = SwC0Replayed::default();
        for (&entry, computed) in entries.iter().zip(results) {
            batch_state.finish_if_needed(Some(entry.batch_idx), batches, heaps, notes)?;
            let batch = &batches[entry.batch_idx];
            match entry.kind {
                SwC0Kind::Excluded { best_evalue } => {
                    let terminates_early =
                        blast_compo_early_termination(best_evalue, &heaps[batch.clone()]);
                    let redone = batch_state.status.reaches_redo(terminates_early, true);
                    debug_assert!(!redone);
                }
                SwC0Kind::Kept { position } => {
                    let local_match = local_match_of(lists, entry.q_idx, position);
                    if blast_compo_early_termination(local_match.best_evalue, &heaps[batch.clone()])
                    {
                        if computed.is_some() {
                            replayed.wasted += 1;
                            replayed.drop_batch(entry.batch_idx);
                        } else {
                            replayed.skipped += 1;
                        }
                        remaining_by_q[entry.q_idx] -= 1;
                        continue;
                    }
                    let Some(SwC0Computed { result, deferred }) = computed else {
                        unreachable!(
                            "CBS 0 SW query {} ordinal {}: ended before window, not ended at replay",
                            entry.q_idx, entry.ordinal
                        );
                    };
                    remaining_by_q[entry.q_idx] -= 1;
                    if let Some((evalue, score)) = deferred {
                        if !heaps[entry.q_idx].would_insert(evalue, score, local_match.oid) {
                            // NCBI's Smith-Waterman redo turns this match away at the heap test,
                            // returning status 0 with no alignments and discarding any speculative
                            // alignments, failure, or endless-redo error.
                            let _ = batch_state
                                .status
                                .redone(Ok::<(), anyhow::Error>(()), || unreachable!())?;
                            replayed.turned += 1;
                            continue;
                        }
                        replayed.continued += 1;
                    }
                    let Some(postprocessed) = batch_state.status.redone(result, || {
                        sw_endless_after_oid(
                            notes,
                            batch.clone(),
                            heaps,
                            entry.ordinal,
                            entry.q_idx,
                        )
                    })?
                    else {
                        continue;
                    };
                    let best_score = postprocessed.best_score;
                    let best_evalue = postprocessed.best_evalue;
                    let Some(hsp_list) = postprocessed.hsp_list else {
                        continue;
                    };
                    let redone_match = &mut heaps[entry.q_idx];
                    if redone_match.would_insert(best_evalue, best_score, local_match.oid) {
                        redone_match.insert(hsp_list, best_evalue, best_score, local_match.oid);
                    }
                }
            }
        }
        Ok(replayed)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn sw_c0_run_windows(
        lists: &mut [Option<BlastpHitList>],
        compact: &mut [Option<BlastpHitList<XPrelimHsp>>],
        heaps: &mut [BlastCompoHeap],
        batches: &[Range<usize>],
        excluded: &[FetchExcludedMatch],
        listed_ordinals: Option<&[u32]>,
        notes: &FetchNotes<'_>,
        limits: SwC0WindowLimits,
        timing: bool,
        compute: &mut dyn FnMut(
            &[Option<BlastpHitList>],
            &[BlastCompoHeap],
            &SwC0Window,
        ) -> Vec<SwC0Computed>,
    ) -> Result<SwC0Counts> {
        let mut remaining_by_q = query_remaining_matches(lists, compact);
        let mut queue = SwC0Queue::new(batches, excluded, listed_ordinals);
        let mut counts = SwC0Counts::default();
        let budget = limits.matches.max(limits.at_risk);
        let mut caps = vec![limits.at_risk; batches.len()];
        let mut batch_state = SwC0BatchState::default();
        let mut released_batch = 0usize;

        while !queue.is_done() {
            let end_above = batch_end_above(batches, heaps, &remaining_by_q);
            let window = build_window(
                &mut queue,
                lists,
                compact,
                heaps,
                &remaining_by_q,
                &end_above,
                &caps,
                limits,
            );
            if window.entries.is_empty() {
                break;
            }
            let new_kept = window
                .entries
                .iter()
                .filter(|entry| matches!(entry.kind, SwC0Kind::Kept { .. }))
                .count();
            counts.matches += new_kept;
            if !window.jobs.is_empty() {
                counts.windows += 1;
                counts.jobs += window.jobs.len();
            }
            let par_start = timing.then(Instant::now);
            let job_results = if window.jobs.is_empty() {
                Vec::new()
            } else {
                compute(lists, heaps, &window)
            };
            if let Some(par_start) = par_start {
                counts.par_s += par_start.elapsed().as_secs_f64();
            }
            let mut results: Vec<Option<SwC0Computed>> = std::iter::repeat_with(|| None)
                .take(window.entries.len())
                .collect();
            for (&at, computed) in window.jobs.iter().zip(job_results) {
                results[at] = Some(computed);
            }
            let replay_start = timing.then(Instant::now);
            let replayed = match sw_c0_replay(
                &window.entries,
                results,
                lists,
                heaps,
                batches,
                &mut remaining_by_q,
                &mut batch_state,
                notes,
            ) {
                Ok(replayed) => replayed,
                Err(error) => {
                    if timing {
                        print_timing(&counts);
                    }
                    return Err(error);
                }
            };
            counts.skipped += replayed.skipped;
            counts.wasted += replayed.wasted;
            counts.turned += replayed.turned;
            counts.continued += replayed.continued;
            if let Some(replay_start) = replay_start {
                counts.replay_s += replay_start.elapsed().as_secs_f64();
            }
            for &(batch_idx, at_risk_count) in &window.at_risk {
                caps[batch_idx] = caps[batch_idx].saturating_add(at_risk_count).min(budget);
            }
            for &batch_idx in &replayed.dropped_batches {
                caps[batch_idx] = limits.at_risk;
            }
            let done_batch = queue
                .pending
                .front()
                .map_or(queue.expanded_batch, |entry| entry.batch_idx);
            for batch in &batches[released_batch..done_batch] {
                counts.released += release_batch(batch, lists);
            }
            released_batch = done_batch;
        }
        if let Err(error) = batch_state.finish_if_needed(None, batches, heaps, notes) {
            if timing {
                print_timing(&counts);
            }
            return Err(error);
        }
        for batch in &batches[released_batch..queue.expanded_batch] {
            counts.released += release_batch(batch, lists);
        }
        Ok(counts)
    }

    struct SwC0Scratch {
        composition_workspace: BlastCompositionWorkspace,
        gap_scratch: GapAlignScratch,
        hits: Vec<BlastpPreliminaryHsp>,
        range_cache: BlastpKappaSubjectRangeCache,
        sw: SwMatchScratch,
        redo_params: BlastRedoAlignParams,
    }

    /// The `Sync` slice of `RedoInputs` read by parallel worker threads (excludes `!Sync`
    /// `FetchNotes`).
    struct SwC0WorkerInputs<'a> {
        subjects: crate::algorithm::blastp::x_redo_subjects::XRedoSubjects<'a>,
        contexts: &'a [crate::algorithm::tblastx::lookup::QueryContext],
        search_spaces: &'a [crate::stats::search_space::SearchSpace],
        evalue_adjust_contexts: &'a [kappa::EvalueAdjustContext],
        matrix: crate::config::ScoringMatrix,
        gapped_params: &'a crate::stats::KarlinParams,
        gapped_gumbel: Option<&'a crate::stats::BlastGumbelBlk>,
        evalue: f64,
        engine_positives: bool,
    }

    #[allow(clippy::too_many_arguments)]
    fn run_match(
        inputs: &SwC0WorkerInputs<'_>,
        unadjusted: &AdjustedProteinMatrix,
        lists: &[Option<BlastpHitList>],
        heaps: &[BlastCompoHeap],
        window: &SwC0Window,
        workspaces: &[OnceLock<BlastpKappaQueryWorkspace>],
        scratch: &mut SwC0Scratch,
        at: usize,
    ) -> SwC0Computed {
        let entry = window.entries[at];
        let SwC0Kind::Kept { position } = entry.kind else {
            unreachable!("only kept entries become compute jobs");
        };
        let q_idx = entry.q_idx;
        let local_match = local_match_of(lists, q_idx, position);
        scratch.hits.clear();
        scratch.hits.extend(
            local_match
                .hsps
                .iter()
                .map(preliminary_hit_from_local_match_hsp),
        );
        let subject = match inputs.subjects.read(local_match.oid) {
            Ok(subject) => subject,
            Err(error) => {
                return SwC0Computed {
                    result: Err(error),
                    deferred: None,
                };
            }
        };
        let subject_raw = subject.residues();
        let ctx = &inputs.contexts[q_idx];
        let query_workspace = workspaces[q_idx].get_or_init(|| {
            let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
            build_query_workspace(
                query_raw,
                ctx.aa_len as i32,
                inputs.search_spaces[q_idx].length_adjustment as i32,
                inputs.search_spaces[q_idx].effective_space,
                inputs.evalue_adjust_contexts[q_idx],
            )
        });
        scratch.redo_params.ccat_query_length = ctx.aa_len as i32;
        let mut deferred = None;
        let before = window.before[at];
        let heap = &heaps[q_idx];
        let result = kappa::postprocess_preliminary_hits_smith_waterman(
            &scratch.hits,
            query_workspace,
            &mut scratch.composition_workspace,
            query_nomask_sequence(ctx),
            subject_raw,
            inputs.matrix,
            Some(unadjusted),
            inputs.gapped_params,
            inputs.gapped_gumbel,
            &scratch.redo_params,
            inputs.evalue,
            &mut scratch.gap_scratch,
            &mut scratch.range_cache,
            inputs.engine_positives,
            &mut scratch.sw,
            None,
            &mut |evalue, score| {
                let inserts = x_sw_heap_test(
                    heap,
                    before,
                    (evalue, score, local_match.oid),
                    &mut deferred,
                );
                // Speculate `Known(true)` when deferred: replay checks `heaps[q].would_insert` in
                // exact OID order and discards the speculative result if FALSE.
                SwHeapAnswer::Known(if deferred.is_some() { true } else { inserts })
            },
            &mut None,
        );
        SwC0Computed { result, deferred }
    }

    fn print_timing(counts: &SwC0Counts) {
        eprintln!(
            "[TIMING] x_bpredopar: matches={} jobs={} windows={} split=0 ended={} skipped={} wasted={} par_s={:.6} replay_s={:.6}{} turned={} continued={} cont_stop=0 cont_ahead=0 cont_bound=0 carried_max=0",
            counts.matches,
            counts.jobs,
            counts.windows,
            counts.skipped + counts.wasted,
            counts.skipped,
            counts.wasted,
            counts.par_s,
            counts.replay_s,
            counts
                .lpt_windows
                .map_or(String::new(), |w| format!(" lpt_windows={w}")),
            counts.turned,
            counts.continued,
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1487
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:678-686
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3447-3449,3525-3781
    /// Windowed parallel Smith-Waterman redo without composition-based statistics
    /// (`-comp_based_stats 0 -use_sw_tback`, item X4-c0): computes matches on the search pool
    /// (with LPT scheduling when `LOSAT_X_BPREDOLPT` is on) and replays early termination,
    /// deferred heap tests, excluded marks, batch status, and heap insertions in exact
    /// `(batch, ordinal ASC, query DESC)` order.
    pub(in super::super) fn sw_c0_windows<P>(
        redo: &RedoInputs<'_, P>,
        preliminary_hit_lists: &mut [Option<BlastpHitList>],
        x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
        redone_matches: &mut [BlastCompoHeap],
        x_prelim_freed: &mut usize,
    ) -> Result<()>
    where
        P: Fn(usize) -> Result<BlastRedoAlignParams> + Sync,
    {
        if redo.kappa_query_batches.is_empty() {
            return Ok(());
        }
        // Build query 0's parameters once up front so `OnceLock` caches `initial_matrix_info`
        // before worker threads clone it, reporting any error before any match.
        let _ = (redo.build_redo_align_params)(0)?;
        let build_redo_align_params = redo.build_redo_align_params;
        let pool = redo.blastp_parallel_pool;
        let worker_inputs = SwC0WorkerInputs {
            subjects: redo.subjects,
            contexts: redo.contexts,
            search_spaces: redo.search_spaces,
            evalue_adjust_contexts: redo.evalue_adjust_contexts,
            matrix: redo.args.scoring.matrix,
            gapped_params: redo.gapped_params,
            gapped_gumbel: redo.gapped_gumbel.as_ref(),
            evalue: redo.args.evalue,
            engine_positives: redo.needs.seqalign,
        };
        let unadjusted = unadjusted_square_matrix(redo.args.scoring.matrix);
        let batches = batch_ranges(redo.kappa_query_batches);
        let threads = redo.blastp_num_threads.max(1);
        let jobs_per_thread = x_window_jobs_per_thread();
        let full_matches = jobs_per_thread
            .saturating_mul(threads)
            .saturating_mul(X_JOB_MATCHES);
        let limits = SwC0WindowLimits {
            matches: full_matches,
            fill_matches: full_matches.min(threads.saturating_mul(32)),
            small_heap_matches: full_matches.min(threads.saturating_mul(64)),
            at_risk: threads.saturating_mul(X_JOB_MATCHES),
        };
        let workspaces: Vec<OnceLock<BlastpKappaQueryWorkspace>> =
            (0..redo.contexts.len()).map(|_| OnceLock::new()).collect();
        let lpt_on = x_redo_par::x_bp_redo_lpt();
        let mut lpt_windows = 0usize;
        let mut counts = sw_c0_run_windows(
            preliminary_hit_lists,
            x_compact_hit_lists,
            redone_matches,
            &batches,
            redo.fetch_excluded,
            redo.listed_ordinals,
            redo.fetch_notes,
            limits,
            redo.timing_enabled,
            &mut |lists, heaps, window| {
                let new_scratch = || SwC0Scratch {
                    composition_workspace: BlastCompositionWorkspace::new_blosum62(),
                    gap_scratch: GapAlignScratch::new(),
                    hits: Vec::new(),
                    range_cache: BlastpKappaSubjectRangeCache::new(None),
                    sw: SwMatchScratch::new(),
                    redo_params: build_redo_align_params(0)
                        .expect("initial matrix info was initialized before the window loop"),
                };
                let use_lpt = lpt_on && {
                    let acts = pool.enabled() && threads > 1 && window.jobs.len() > 1;
                    crate::utils::xreach::x_reach!("LOSAT_X_BPREDOLPT", acts);
                    acts
                };
                if use_lpt {
                    lpt_windows += 1;
                    let costs: Vec<u64> = window
                        .jobs
                        .iter()
                        .map(|&at| {
                            let entry = window.entries[at];
                            let SwC0Kind::Kept { position } = entry.kind else {
                                return 0;
                            };
                            let query_len = worker_inputs.contexts[entry.q_idx].aa_len as u64;
                            let subject_len = worker_inputs
                                .subjects
                                .aa_len(local_match_of(lists, entry.q_idx, position).oid as usize)
                                .map_or(0, |length| length as u64);
                            query_len.saturating_mul(subject_len)
                        })
                        .collect();
                    let order = x_lpt_order(&costs);
                    x_compute_window_lpt(
                        pool,
                        threads.min(window.jobs.len()),
                        &order,
                        new_scratch,
                        |scratch, job_idx| {
                            run_match(
                                &worker_inputs,
                                &unadjusted,
                                lists,
                                heaps,
                                window,
                                &workspaces,
                                scratch,
                                window.jobs[job_idx],
                            )
                        },
                    )
                } else {
                    pool.install(|| {
                        window
                            .jobs
                            .par_iter()
                            .map_init(new_scratch, |scratch, &at| {
                                run_match(
                                    &worker_inputs,
                                    &unadjusted,
                                    lists,
                                    heaps,
                                    window,
                                    &workspaces,
                                    scratch,
                                    at,
                                )
                            })
                            .collect()
                    })
                }
            },
        )?;
        counts.lpt_windows = lpt_on.then_some(lpt_windows);
        if redo.timing_enabled {
            print_timing(&counts);
        }
        *x_prelim_freed += counts.released;
        Ok(())
    }

    #[cfg(test)]
    pub(super) use {
        sw_c0_run_windows as test_run_windows, SwC0Computed as TestComputed,
        SwC0Window as TestWindow, SwC0WindowLimits as TestWindowLimits,
    };
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) use par::sw_c0_windows;

#[cfg(test)]
#[path = "sw_c0_tests.rs"]
mod tests;
