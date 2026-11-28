//! The composition-based-statistics redo of one query batch's matches (a gapped search with
//! `-comp_based_stats` other than 0): the kappa branch of `Runtime::batch_end`. It reads the batch's
//! matches from the stream, redoes them (a batch of them ahead of time on the search pool, then each
//! in stream order), keeps the redone lists in the per-query heaps and pours the heaps into the
//! queries' hit lists.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2980-2981
//! (`Blast_RedoAlignmentCore_MT`, which the loop over the matches of the stream and the heaps'
//! pour-out belong to)
//!
//! This module is a child of `runtime` for the types it keeps private (`Runtime`, `XKappaWorker`,
//! `XMatrixLeft`, `finalize_redone`, `x_speculate_match`, the parallel rules). `kappa_redo_batch` is
//! the statements of `batch_end`'s kappa branch, in their order; what it reads of the batch is
//! `KappaBatch`.

use super::{
    finalize_redone, CompoHeap, ContextParameters, HeapEntry, HitList, HspList, HspStream,
    KappaState, PreparedQueryBatch, Result, Runtime, XMatrixLeft,
};
#[cfg(feature = "parallel")]
use super::{x_bx_par_rule, x_inner_serial, x_kappa_workers, x_speculate_match, XKappaWorker};
use crate::algorithm::blastx::kappa::{sw_batch_failed, sw_endless};
use crate::algorithm::blastx::x_timing::{self, End as XEnd};
use crate::algorithm::common::database::{FetchMarks, FetchedUtf8, Utf8Read};
use crate::blastinput::app::ErrorAfterBatchReports;
use crate::core::composition_adjustment::smith_waterman_match::SwBatchStatus;

/// No NCBI counterpart (scheduling only): EXPERIMENT (LOSAT_X_BXREDOFLOW, in the strict set: on
/// unless `=0`), the matches redone ahead as a flow instead of in fork-join batches (see
/// `kappa_redo_batch`). It acts wherever the batches would run (LOSAT_X_BXPAR's rule: a pool of
/// two or more threads, two or more matches, no recording observer, not e-value mode of
/// `-use_sw_tback`), and has no rule of its own; its workers are the batches': one per thread of
/// the pool.
#[cfg(feature = "parallel")]
fn x_redo_flow() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXREDOFLOW"))
}

/// No NCBI counterpart (scheduling only): the guess of which matches the redo makes ahead (the flow's
/// offers and `=0`'s batches, `kappa_redo_batch`), with the matches offered ahead of the loop that
/// it counts. A match is offered when the loop's early termination test (`CompoHeap::early`) will
/// let it reach its redo even after the matches offered before it have filled their heaps: each
/// offered match whose best E-value is at or below its heap's cutoff (`CompoHeap::x_counts_ahead`)
/// counts as one more entry there, from the loop's match on (`CompoHeap::x_early_ahead`). A match
/// that only those counted matches would end early is "not yet": nothing at or past it is offered
/// until the loop has passed them and the heaps themselves decide. Why: with a short hit list a few
/// inserts fill the heaps, while the heaps were empty when the whole window was offered (1 x 9,000
/// x 595 `-use_sw_tback -max_target_seqs 3`: the loop redoes 32 of 205 matches; the flow redid
/// 130 ahead at 4 threads and 204 at 8, its loop up to 6 times one thread's; with the matches
/// counted, the 31 or 32 that the loop takes, in 0.6 times one thread's time). With a long hit
/// list (the default 500) the counted matches do not reach a heap's threshold but at its very end,
/// so the matches offered are those of the heaps as they are. The guess changes the time only: the
/// loop tests and redoes every match in stream order.
#[cfg(feature = "parallel")]
struct XPending {
    /// Per query: the counted matches offered from the loop's match on.
    queries: Vec<usize>,
    /// Per match: whether it is counted in `queries`.
    counted: Vec<bool>,
}

#[cfg(feature = "parallel")]
impl XPending {
    fn new(queries: usize, matches: usize) -> Self {
        Self {
            queries: vec![0; queries],
            counted: vec![false; matches],
        }
    }

    /// The loop is at match `index` of `matches`: the match before it, inserted or skipped, is no
    /// longer ahead of the loop.
    fn at(&mut self, index: usize, matches: &[HspList]) {
        if let Some(before) = index.checked_sub(1) {
            if std::mem::take(&mut self.counted[before]) {
                self.queries[matches[before].query_index] -= 1;
            }
        }
    }

    /// The guess for match `index`, `list` (`excluded`: `-use_sw_tback`'s fetch excludes its
    /// subject, so it is never redone), on `heaps` as they are now: `Some(true)` offer it (and
    /// count it), `Some(false)` the loop will skip it, `None` not yet (see `XPending`).
    fn offer(
        &mut self,
        index: usize,
        list: &HspList,
        excluded: bool,
        heaps: &[CompoHeap],
    ) -> Option<bool> {
        let e = list.best_evalue;
        if excluded {
            return Some(false);
        }
        if !CompoHeap::x_early_ahead(e, heaps, &self.queries) {
            if heaps[list.query_index].x_counts_ahead(e) {
                self.queries[list.query_index] += 1;
                self.counted[index] = true;
            }
            return Some(true);
        }
        // Early on the heaps as they are: it stays so (no insert takes an entry away, `early`).
        // Otherwise only the counted matches would end it early: not yet.
        if CompoHeap::early(e, heaps) {
            Some(false)
        } else {
            None
        }
    }
}

/// What the kappa branch of `batch_end` reads and writes of the batch: the batch with the
/// parameters of its contexts, the Kappa state, the stream the matches are read from, and the
/// queries' hit lists with the heaps that fill them (`heaps`) and the inactive ones (`inactive_heaps`,
/// only for their trace numbers).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2990,2995
/// (`Blast_RedoAlignmentCore_MT`'s `hsp_stream` and `results`)
pub(super) struct KappaBatch<'a> {
    pub(super) batch: &'a PreparedQueryBatch,
    pub(super) full: &'a [ContextParameters],
    pub(super) state: &'a mut KappaState,
    pub(super) stream: &'a mut HspStream,
    pub(super) queries: &'a mut [HitList],
    pub(super) heaps: &'a mut [CompoHeap],
    pub(super) inactive_heaps: &'a mut [CompoHeap],
    /// The subjects of the stream whose header holds a UTF8String where a VisibleString is
    /// expected and those whose header cannot be read (`FetchExcluded::retain_noting`: with the
    /// query of the list, its best E-value and whether it has HSPs), the lists left in the stream
    /// (query and subject, in the stream's order) where there is one, and the first of the former
    /// that the redo fetches (`redo_fetched`; the Smith-Waterman redo's, `sw_fetched`, which reads
    /// its own stream's). A fetch of one of the latter is LOSAT's refusal (`fetch_refused`).
    pub(super) fetched: (
        &'a [FetchedRead],
        &'a [FetchedRead],
        &'a [(usize, usize)],
        &'a mut Option<u32>,
    ),
}

/// A subject of the batch's lists whose header holds a UTF8String where a VisibleString is
/// expected (`FetchMarks::of_noting`), with its query, its list's best E-value and whether the
/// list has HSPs.
pub(super) type FetchedRead = (FetchedUtf8, (usize, f64, bool));

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3575 (cited at
// `blastp::blast_engine::database::exclude_fetched`)
/// The first of `reads` (the Smith-Waterman redo's, `FetchMarks::of_noting` of its stream, whose
/// lists the fetch excludes stay in it as markers) that NCBI's redo fetches, of the matches up to
/// the stream's place `until` (the match whose redo fails the batch, or every match):
/// `common::database::first_redo_fetch`, with each query's heap as at the turn (the records of
/// `heaps`, as the redo left them, from the lists of `stream` before it).
#[cold]
#[inline(never)]
fn sw_fetched(
    reads: &[FetchedRead],
    stream: &[(usize, usize)],
    heaps: &[CompoHeap],
    until: usize,
) -> Option<u32> {
    crate::algorithm::common::database::first_redo_fetch(reads, stream, until, |before| {
        heaps
            .iter()
            .enumerate()
            .map(|(query, heap)| heap.redo_heap(|oid| before.contains(&(query, oid))))
            .collect()
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3575 (cited at
// `blastp::blast_engine::database::exclude_fetched`)
/// The first of `reads` that NCBI's composition-based redo fetches, in its order: the stream's
/// read from its end (`order`: the lists that the fetches left, in the stream's order), each
/// match's turn where its early termination test leaves it (`common::database::redo_fetches`, the
/// heaps `heaps` as the redo left them: a record counts at a match's turn where its list came
/// before the match's).
#[cold]
#[inline(never)]
fn redo_fetched(
    reads: &[FetchedRead],
    order: &[(usize, usize)],
    heaps: &[CompoHeap],
) -> Option<u32> {
    // A list's turn: the kept lists the redo takes before it (it reads the stream from its end).
    let turn_of = |query: usize, oid: usize| {
        order
            .iter()
            .position(|&list| list == (query, oid))
            .map(|at| order.len() - 1 - at)
    };
    reads
        .iter()
        .filter(|(_, (_, _, has_hsps))| *has_hsps)
        .filter(|(read, (_, evalue, _))| {
            let turn = order.len() - read.before - usize::from(read.kept);
            let redo_heaps: Vec<_> = heaps
                .iter()
                .enumerate()
                .map(|(query, heap)| {
                    heap.redo_heap(|oid| turn_of(query, oid).is_some_and(|at| at < turn))
                })
                .collect();
            crate::algorithm::common::database::redo_fetches(*evalue, &redo_heaps)
        })
        .max_by_key(|(read, _)| read.index)
        .map(|(read, _)| read.ordinal)
}

impl Runtime<'_> {
    /// LOSAT's refusal of the search whose redo fetches the sequence `ordinal`, whose header
    /// cannot be read: NCBI's search ends there (`FetchExcluded::refusal_at`).
    #[cold]
    #[inline(never)]
    fn fetch_refused(&self, ordinal: u32) -> anyhow::Error {
        match self.subjects.seqdb() {
            Some(db) => {
                crate::algorithm::common::database::FetchExcluded::refusal_at(db, ordinal, "BLASTX")
            }
            None => anyhow::anyhow!("a fetch of sequence {ordinal} without a database"),
        }
    }

    /// `error`, which fails the batch after its Smith-Waterman redo fetched `fetched` (the first
    /// of the batch's fetches that read a header with a UTF8String where a VisibleString is
    /// expected, `sw_fetched`): NCBI posted its warning before the error (rows
    /// `db.t5d.blastx.nsw.fail.*`), which the batch's error then carries
    /// (`ErrorAfterBatchReports::utf8`, written before it by `native::failed_batch_reads`).
    #[cold]
    #[inline(never)]
    fn sw_failed_after(&self, error: anyhow::Error, fetched: Option<u32>) -> anyhow::Error {
        let Some(db) = self.subjects.seqdb() else {
            return error;
        };
        let Some(path) = fetched
            .and_then(|ordinal| db.location_of(ordinal))
            .and_then(|at| db.utf8_string_path_at(at))
        else {
            return error;
        };
        // The batch's error (`kappa::sw_batch_failed`, `kappa::sw_endless`, or the redo's own)
        // with the read before it; `native::failed_batch_reads` writes the warning, then that
        // error.
        ErrorAfterBatchReports {
            query: Some(self.ordinal),
            error,
            utf8: Some(Utf8Read {
                query: self.ordinal,
                path,
            }),
            title_reads: Vec::new(),
        }
        .into()
    }

    /// Redo the matches of the batch's stream with composition-based statistics and leave the
    /// result in the queries' hit lists; returns the clock of the finish interval of
    /// `LOSAT_TIMING` (`None` when it is off).
    ///
    /// `SW`: `-use_sw_tback` (item X4b-2), the Smith-Waterman redo of each match, with or without
    /// composition-based statistics (`KappaState`'s redo), chosen once per batch. Its batch
    /// keeps NCBI's status of the last match the redo reached (core's `SwBatchStatus`, the rule's
    /// one owner): a list whose subject the fetch excludes stays in the stream as a marker that
    /// clears it (`FetchMarks`); a failed last match fails the batch (`sw_batch_failed`). In
    /// e-value mode the redo of a window tests the heap that its query context indexes (NCBI's
    /// heap of that query of the batch, or an empty heap past the batch's heaps,
    /// `CompoHeap::would_insert_at`), so the matches are redone one by one in stream order; in
    /// cutoff-score mode (sum statistics) a match's redo reads no heap and no state of the match
    /// before it (every window adjusts its own matrix), so the pool redoes them ahead as for the
    /// X-drop redo. The result is NCBI's `-num_threads 1` at every thread count.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2980-2981
    pub(super) fn kappa_redo_batch<const SW: bool>(
        &mut self,
        args: KappaBatch<'_>,
    ) -> Result<Option<std::time::Instant>> {
        let KappaBatch {
            batch,
            full,
            state,
            stream,
            queries,
            heaps,
            inactive_heaps,
            fetched: (fetched_reads, fetched_unreadable, fetched_order, fetched_first),
        } = args;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3392-3424
        let x_t = x_timing::lap();
        let mut matches = Vec::new();
        while let Some(list) = stream.read_observed(self.observer) {
            matches.push(list);
        }
        x_timing::end(XEnd::Read, x_t);
        // `-use_sw_tback`: the lists that NCBI's fetch excludes (a `-negative_seqidlist` that keeps
        // a subject without a line), by their place in the stream, and the batch's status.
        // The subjects whose header holds a UTF8String where a VisibleString is expected come
        // from the same read of the headers (`sw_reads`), with the stream's lists, where there is
        // one (`sw_stream`: query and subject, markers too).
        // Those whose header cannot be read stay in the stream (`sw_unreadable`): a fetch of one
        // is LOSAT's refusal (`fetch_refused`).
        let (marks, sw_reads, sw_unreadable, do_link_hsps) = if SW {
            let excluded = self
                .subjects
                .seqdb()
                .and_then(crate::algorithm::common::database::FetchExcluded::of);
            let (marks, reads) = FetchMarks::of_noting(
                excluded.as_ref(),
                &matches,
                |list| u32::try_from(list.oid).ok(),
                |list| (list.query_index, list.best_evalue, !list.hsps.is_empty()),
            );
            (marks, reads.utf8, reads.unreadable, state.do_link_hsps())
        } else {
            (FetchMarks::default(), Vec::new(), Vec::new(), false)
        };
        let sw_stream: Vec<(usize, usize)> = if sw_reads.is_empty() && sw_unreadable.is_empty() {
            Vec::new()
        } else {
            matches
                .iter()
                .map(|list| (list.query_index, list.oid))
                .collect()
        };
        // Whether the Smith-Waterman redo reads the heaps of the batch (e-value mode), so that
        // a match's outcome depends on the matches before it: no speculation then.
        #[cfg(feature = "parallel")]
        let reads_heaps = SW && !do_link_hsps;
        let mut sw_status = SwBatchStatus::default();
        // Redo the matches of one batch on the search pool before walking them in stream order.
        // A match is a function of its own HSP list and of the matrix left by the match before
        // it; `x_speculate_match` returns a result only when the latter was never read, and
        // reports the matrix the match leaves.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3448-3451,3499-3502
        // NCBI redoes the matches in this omp loop, each thread with its own state. The batch
        // of speculative redoes below does the same per-match work on the pool. A result is
        // used only if the match did not read the matrix of its predecessor; every other match
        // is redone in stream order (see x_speculate_match).
        #[cfg(feature = "parallel")]
        let mut speculated: Vec<Option<(HspList, i32, Option<XMatrixLeft>)>> = Vec::new();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334,3499-3502
        // One Kappa state per pool thread, as NCBI keeps one NRrecord, redo_align_params and
        // matrix per thread.
        // EXPERIMENT (LOSAT_X_BXREDOFLOW): the form of the redo made ahead (below); a test's
        // observer may choose it (`RuntimeObserver::x_redo_flow`).
        #[cfg(feature = "parallel")]
        let x_flow = self.observer.x_redo_flow().unwrap_or_else(x_redo_flow);
        #[cfg(feature = "parallel")]
        let x_workers = if x_bx_par_rule(|| {
            self.observer.x_passive()
                && !x_inner_serial()
                && rayon::current_thread_index().is_some()
                && rayon::current_num_threads() > 1
                && matches.len() > 1
                && !reads_heaps
        }) {
            // The workers of the flow and of the batches: one per thread of the pool (at most one
            // per match). Not the share the cap leaves a batch of one task for its subject flow
            // (`x_bx_mem::flow_threads`): a redo that ends a search alone would keep that share
            // (4 of 32 threads at the end of q200k `-use_sw_tback`, 2.5 times the batches' time).
            (0..x_kappa_workers(rayon::current_num_threads(), matches.len()))
                .map(|_| {
                    Ok(std::sync::Mutex::new(XKappaWorker(state.x_for_thread(
                        batch,
                        full,
                        self.options,
                    )?)))
                })
                .collect::<Result<Vec<_>>>()?
        } else {
            Vec::new()
        };
        // `x_workers` has one entry per pool thread, or one per match where there are fewer
        // matches (`x_kappa_workers`); none outside a pool (asking Rayon for the thread count
        // there would start its global pool).
        // No NCBI counterpart: size of one speculative batch (scheduling only); it does not
        // change any value NCBI computes.
        #[cfg(feature = "parallel")]
        let x_batch = x_workers.len().saturating_mul(32).max(64);
        // No NCBI counterpart: stage report for the thread-use statistics; it does not change
        // any value NCBI computes.
        #[cfg(feature = "parallel")]
        crate::utils::threading::report_stage(
            "blastx",
            "kappa_redo",
            matches.len(),
            !x_workers.is_empty(),
        );
        // EXPERIMENT (LOSAT_X_BXREDOFLOW, in the strict set: on unless `=0`): the matches ahead
        // are redone as a flow (`utils::xahead`) instead of in the fork-join batches below: the
        // loop's thread offers the next `x_batch` matches as it walks the stream (those the guess
        // `XPending` lets through, up to the first it cannot tell yet; the batches take the same
        // guess), each worker takes the next one offered, and the loop takes a match's result when
        // it reaches it (redoing other offered matches while one is still being redone). No
        // thread waits for the longest match of a batch. The workers (`x_kappa_workers`: one per
        // thread of the pool, at most one per match) are the loop's thread's and one per helper.
        // A thread with nothing to redo spins, yields, then sleeps 50 us at a time
        // (`Ahead::x_new_napping`). A match's result does not depend on who redoes it or when
        // (`x_speculate_match`), so the schedule changes the time only.
        #[cfg(feature = "parallel")]
        let (x_workers, mut x_own, x_helpers) = if !x_workers.is_empty() && x_flow {
            let mut all = x_workers
                .into_iter()
                .map(|worker| worker.into_inner().unwrap_or_else(|p| p.into_inner()));
            let own = all.next();
            (Vec::new(), own, all.collect::<Vec<_>>())
        } else {
            (x_workers, None, Vec::new())
        };
        #[cfg(feature = "parallel")]
        let x_ahead = x_own
            .is_some()
            .then(|| crate::utils::xahead::Ahead::x_new_napping(x_batch));
        #[cfg(feature = "parallel")]
        let (x_options, x_subjects, x_db_length) = (self.options, self.subjects, self.db_length);
        // No NCBI counterpart: the flow's redo of one match ahead, by the statements of the
        // batches' (`x_speculate_match`).
        #[cfg(feature = "parallel")]
        let x_compute = |index: usize, worker: &mut XKappaWorker| {
            x_speculate_match(
                &mut worker.0,
                &matches[index],
                batch,
                full,
                x_options,
                x_subjects,
                x_db_length,
            )
        };
        // No NCBI counterpart: counters of speculative results used and matches redone in
        // order; it does not change any value NCBI computes.
        let mut x_reused = 0usize;
        let mut x_serial = 0usize;
        // No NCBI counterpart (scheduling only): the guess of the matches made ahead, with the
        // matches offered ahead of the loop that it counts (`XPending`), and (`=0`) the match where
        // the next batch starts and the first match of the batch in `speculated`.
        #[cfg(feature = "parallel")]
        let mut x_pending = XPending::new(heaps.len(), matches.len());
        #[cfg(feature = "parallel")]
        let (mut x_batch_at, mut x_batch_from) = (0usize, 0usize);
        let x_t_loop = x_timing::lap();
        let mut x_walk = || -> Result<()> {
            for (match_index, list) in matches.iter().enumerate() {
                #[cfg(feature = "parallel")]
                x_pending.at(match_index, &matches);
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3448-3451
                // The flow offers the matches up to `x_batch` past this one that it has not offered
                // yet, with the guess (`XPending`) on the heaps as they are now; it stops at the
                // first match the guess cannot tell yet and goes on from there at a later match.
                #[cfg(feature = "parallel")]
                if let Some(ahead) = &x_ahead {
                    ahead.x_announce_while(match_index, matches.len(), |next| {
                        x_pending.offer(next, &matches[next], SW && marks.excluded(next), heaps)
                    });
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3448-3451
                // A batch of the following matches is redone ahead of time on the pool. The early-
                // termination test (BlastCompo_EarlyTermination) is guessed here on the heaps as
                // they are now (`XPending`), and evaluated, in stream order, below; a result is
                // used only when the second test lets the match reach its redo.
                #[cfg(feature = "parallel")]
                if !x_workers.is_empty() && match_index == x_batch_at {
                    use rayon::prelude::*;
                    let window = &matches[match_index..(match_index + x_batch).min(matches.len())];
                    // The batch ends before the first match the guess cannot tell yet; the next
                    // batch starts there when the loop reaches it (every match before it then
                    // done). The first match is always told: nothing is counted ahead of it.
                    let mut wanted: Vec<bool> = Vec::with_capacity(window.len());
                    for (position, list) in window.iter().enumerate() {
                        let index = match_index + position;
                        match x_pending.offer(index, list, SW && marks.excluded(index), heaps) {
                            Some(offered) => wanted.push(offered),
                            None => break,
                        }
                    }
                    let batch_lists = &window[..wanted.len()];
                    (x_batch_at, x_batch_from) = (match_index + wanted.len().max(1), match_index);
                    let (options, subjects, db_length) =
                        (self.options, self.subjects, self.db_length);
                    let full = &full;
                    let workers = &x_workers;
                    let x_t = x_timing::lap();
                    // One worker per thread: the thread's own. One per match (fewer matches than
                    // threads, so this batch holds them all): the match's own.
                    let per_thread = workers.len() == rayon::current_num_threads();
                    speculated = batch_lists
                        .par_iter()
                        .zip(wanted.par_iter())
                        .enumerate()
                        .map(|(position, (list, &wanted))| {
                            if !wanted {
                                return None;
                            }
                            let slot = if per_thread {
                                rayon::current_thread_index().unwrap_or(0) % workers.len()
                            } else {
                                position % workers.len()
                            };
                            let mut guard = workers[slot].lock().expect("BLASTX redo worker state");
                            x_speculate_match(
                                &mut guard.0,
                                list,
                                batch,
                                full,
                                options,
                                subjects,
                                db_length,
                            )
                        })
                        .collect();
                    x_timing::end(XEnd::KappaSpec, x_t);
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1580
                let early = CompoHeap::early(list.best_evalue, heaps);
                self.observer.early(list.best_evalue, heaps.len(), early);
                // The early termination test, then (`-use_sw_tback`) the fetch, which may exclude
                // the subject: a list that reaches the fetch sets the batch's status
                // (`SwBatchStatus`).
                let skipped = if SW {
                    !sw_status.reaches_redo(early, marks.excluded(match_index))
                } else {
                    early
                };
                if skipped {
                    // The flow's workers need not redo it.
                    #[cfg(feature = "parallel")]
                    if let Some(ahead) = &x_ahead {
                        ahead.skip(match_index);
                    }
                    continue;
                }
                #[cfg(feature = "parallel")]
                let x_ready = if let Some(ahead) = &x_ahead {
                    let x_t = x_timing::lap();
                    let own = x_own.as_mut().expect("the flow's loop keeps a redo worker");
                    let ready = ahead.take(match_index, own, &x_compute).flatten();
                    x_timing::end(XEnd::KappaSpec, x_t);
                    ready
                } else if x_workers.is_empty() {
                    None
                } else {
                    speculated
                        .get_mut(match_index - x_batch_from)
                        .and_then(Option::take)
                };
                #[cfg(not(feature = "parallel"))]
                let x_ready: Option<(HspList, i32, Option<XMatrixLeft>)> = None;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3626
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1232-1234
                // Dispatch point: the reference path is the else branch, which redoes the match in
                // stream order (Blast_RedoOneMatch). A speculative result is taken instead only
                // when it was computed without reading its predecessor's matrix; its matrix, if it
                // leaves one, is handed on.
                let (output, best_score) = if let Some((output, best_score, matrix)) = x_ready {
                    x_reused += 1;
                    if let Some(matrix) = matrix {
                        state.x_set_matrix(Some(*matrix));
                    }
                    if SW {
                        // A match redone ahead ended (a failed or endless redo is redone here).
                        sw_status.redone(Ok(()), || sw_endless(do_link_hsps))?;
                    }
                    (output, best_score)
                } else {
                    x_serial += 1;
                    let x_t = x_timing::lap();
                    let encoded = self.subjects.read_searched(list.oid, self.options)?;
                    let input: Vec<_> = list.hsps.iter().map(|h| h.hsp.clone()).collect();
                    let redone = if SW {
                        // NCBI's heap test in e-value mode, with the subject's OID
                        // (`matchingSeq->index`), on the heap that the window's query context
                        // indexes (`CompoHeap::would_insert_at`).
                        let oid = list.oid;
                        let result = state.redo_list_tested(
                            &input,
                            batch,
                            &encoded,
                            &mut |context, evalue, score| {
                                Ok(CompoHeap::would_insert_at(
                                    heaps, context, evalue, score, oid,
                                ))
                            },
                        );
                        match sw_status
                            .redone(result, || sw_endless(do_link_hsps))
                            .map_err(|error| {
                                // A fetch before this match that read a header that cannot be
                                // read ended NCBI's search there.
                                match sw_fetched(&sw_unreadable, &sw_stream, heaps, match_index) {
                                    Some(ordinal) => self.fetch_refused(ordinal),
                                    None => self.sw_failed_after(
                                        error,
                                        sw_fetched(&sw_reads, &sw_stream, heaps, match_index),
                                    ),
                                }
                            })? {
                            Some(redone) => redone,
                            None => continue,
                        }
                    } else {
                        state.redo_list_observed(&input, batch, &encoded)?
                    };
                    let finalized = finalize_redone(
                        &redone,
                        batch,
                        full,
                        self.options,
                        &encoded,
                        self.db_length,
                        list.query_index,
                        list.oid,
                        self.observer,
                    );
                    x_timing::end(XEnd::KappaSerial, x_t);
                    finalized?
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3687
                // The match is kept when its best e-value is at most `-evalue`; an `-evalue` of NaN
                // (NCBI's option check lets it pass) keeps none, though the reaping keeps every
                // HSP. (An empty list is not kept either: its best e-value is DBL_MAX, above every
                // `-evalue` that LOSAT runs with the redo.)
                let significant = output.best_evalue <= self.options.evalue;
                if output.hsps.is_empty() || !significant {
                    continue;
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
                let heap = &mut heaps[list.query_index];
                let e = output.best_evalue;
                let insert = heap.would_insert(e, best_score, list.oid);
                self.observer.heap(
                    list.query_index + inactive_heaps.len(),
                    "WOULD",
                    e,
                    best_score,
                    list.oid as i64,
                    i64::from(insert),
                    heap,
                );
                if insert {
                    let discarded = heap.insert(HeapEntry {
                        list: output,
                        best_score,
                    });
                    self.observer.heap(
                        list.query_index + inactive_heaps.len(),
                        "INSERT",
                        e,
                        best_score,
                        list.oid as i64,
                        discarded.as_ref().map_or(-1, |h| h.list.oid as i64),
                        heap,
                    );
                }
            }
            Ok(())
        };
        // The loop on this thread; with the flow, its helpers on the pool's other threads.
        #[cfg(feature = "parallel")]
        crate::utils::xahead::with_helpers(x_ahead.as_ref(), x_helpers, &x_compute, &mut x_walk)?;
        #[cfg(not(feature = "parallel"))]
        x_walk()?;
        #[cfg(feature = "parallel")]
        if let Some(ahead) = &x_ahead {
            x_timing::end_add(XEnd::Ahead, ahead.counts().0);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2493-2513
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:444-464
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3775-3779
        // The status of the last match the redo reached fails the batch.
        // A fetch of the redo that read a header that cannot be read ended NCBI's search before.
        let refused = if SW {
            sw_fetched(&sw_unreadable, &sw_stream, heaps, usize::MAX)
        } else if !fetched_unreadable.is_empty() {
            redo_fetched(fetched_unreadable, fetched_order, heaps)
        } else {
            None
        };
        if let Some(ordinal) = refused {
            return Err(self.fetch_refused(ordinal));
        }
        if SW && sw_status.failed() {
            let error = sw_batch_failed(heaps.first().is_some_and(|heap| !heap.entries.is_empty()));
            return Err(
                self.sw_failed_after(error, sw_fetched(&sw_reads, &sw_stream, heaps, usize::MAX))
            );
        }
        x_timing::end(XEnd::Loop, x_t_loop);
        x_timing::end_add(XEnd::Matches, matches.len() as u64);
        x_timing::end_add(XEnd::Reused, x_reused as u64);
        x_timing::end_add(XEnd::Serial, x_serial as u64);
        let x_t_finish = x_timing::lap();
        let _ = (x_reused, x_serial);
        if SW {
            *fetched_first = sw_fetched(&sw_reads, &sw_stream, heaps, usize::MAX);
        } else if !fetched_reads.is_empty() {
            *fetched_first = redo_fetched(fetched_reads, fetched_order, heaps);
        }
        for (q, (heap, query)) in heaps.iter_mut().zip(&mut *queries).enumerate() {
            loop {
                let entry = heap.pop();
                self.observer.heap(
                    q + inactive_heaps.len(),
                    "POP",
                    0.0,
                    0,
                    entry.as_ref().map_or(-1, |e| e.list.oid as i64),
                    0,
                    heap,
                );
                let Some(entry) = entry else {
                    break;
                };
                query.update_observed(entry.list, self.observer);
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2503-2513
        for query in &mut *queries {
            query.lists.reverse();
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3787-3814
        if let Some(heap) = heaps.first_mut() {
            let entry = heap.pop();
            self.observer.heap(
                inactive_heaps.len(),
                "POP",
                0.0,
                0,
                entry.as_ref().map_or(-1, |e| e.list.oid as i64),
                0,
                heap,
            );
        }
        if let Some(heap) = inactive_heaps.first_mut() {
            let entry = heap.pop();
            self.observer.heap(
                0,
                "POP",
                0.0,
                0,
                entry.as_ref().map_or(-1, |e| e.list.oid as i64),
                0,
                heap,
            );
        }
        Ok(x_t_finish)
    }
}

#[cfg(all(test, feature = "parallel"))]
#[path = "kappa_batch_tests.rs"]
mod tests;
