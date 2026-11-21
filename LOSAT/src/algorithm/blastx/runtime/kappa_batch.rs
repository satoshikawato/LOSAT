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
use crate::algorithm::blastx::x_timing::{self, End as XEnd};

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
    /// expected (`FetchExcluded::retain_noting`: with the query of the list, its best E-value and
    /// whether it has HSPs), the lists left in the stream (query and subject, in the stream's
    /// order) where there is one, and the first that the redo fetches (`redo_fetched`).
    pub(super) fetched: (
        &'a [(
            crate::algorithm::common::database::FetchedUtf8,
            (usize, f64, bool),
        )],
        &'a [(usize, usize)],
        &'a mut Option<u32>,
    ),
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
    reads: &[(
        crate::algorithm::common::database::FetchedUtf8,
        (usize, f64, bool),
    )],
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
    /// Redo the matches of the batch's stream with composition-based statistics and leave the
    /// result in the queries' hit lists; returns the clock of the finish interval of
    /// `LOSAT_TIMING` (`None` when it is off).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2980-2981
    pub(super) fn kappa_redo_batch(
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
            fetched: (fetched_reads, fetched_order, fetched_first),
        } = args;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3392-3424
        let x_t = x_timing::lap();
        let mut matches = Vec::new();
        while let Some(list) = stream.read_observed(self.observer) {
            matches.push(list);
        }
        x_timing::end(XEnd::Read, x_t);
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
        #[cfg(feature = "parallel")]
        let x_workers = if x_bx_par_rule(|| {
            self.observer.x_passive()
                && !x_inner_serial()
                && rayon::current_thread_index().is_some()
                && rayon::current_num_threads() > 1
                && matches.len() > 1
        }) {
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
        // No NCBI counterpart: counters of speculative results used and matches redone in
        // order; it does not change any value NCBI computes.
        let mut x_reused = 0usize;
        let mut x_serial = 0usize;
        let x_t_loop = x_timing::lap();
        #[cfg_attr(
            not(feature = "parallel"),
            expect(unused_variables, reason = "used only with parallel")
        )]
        for (match_index, list) in matches.iter().enumerate() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3448-3451
            // A batch of the following matches is redone ahead of time on the pool. The early-
            // termination test (BlastCompo_EarlyTermination) is evaluated here on the heaps as
            // they are now, and again, in stream order, below; a result is used only when the
            // second test agrees.
            #[cfg(feature = "parallel")]
            if !x_workers.is_empty() && match_index % x_batch == 0 {
                use rayon::prelude::*;
                let batch_lists = &matches[match_index..(match_index + x_batch).min(matches.len())];
                // Early termination only becomes true as the heaps fill.
                let wanted: Vec<bool> = batch_lists
                    .iter()
                    .map(|list| !CompoHeap::early(list.best_evalue, heaps))
                    .collect();
                let (options, subjects, db_length) = (self.options, self.subjects, self.db_length);
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
            if early {
                continue;
            }
            #[cfg(feature = "parallel")]
            let x_ready = if x_workers.is_empty() {
                None
            } else {
                speculated[match_index % x_batch].take()
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
                (output, best_score)
            } else {
                x_serial += 1;
                let x_t = x_timing::lap();
                let encoded = self.subjects.read_searched(list.oid, self.options)?;
                let input: Vec<_> = list.hsps.iter().map(|h| h.hsp.clone()).collect();
                let redone = state.redo_list_observed(&input, batch, &encoded)?;
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
            if output.hsps.is_empty() {
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2493-2513
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:444-464
        x_timing::end(XEnd::Loop, x_t_loop);
        x_timing::end_add(XEnd::Matches, matches.len() as u64);
        x_timing::end_add(XEnd::Reused, x_reused as u64);
        x_timing::end_add(XEnd::Serial, x_serial as u64);
        let x_t_finish = x_timing::lap();
        let _ = (x_reused, x_serial);
        if !fetched_reads.is_empty() {
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
