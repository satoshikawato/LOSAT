//! The subject passes of one blastn query batch (`SubjectPasses`): the search of one subject by stage
//! (`search_subject`, which was the closure `process_subject` of `search_query_batch`), the pass over a
//! subject set's subjects in subject order (`search_subjects`, serial, on the pool, or with
//! `LOSAT_X_BNSTREAM`), the preliminary pass (`prelim_pass`: a subject set or a database, into the HSP
//! collector or, with `-culling_limit`, NCBI's culling writer, with the `-ungapped` finish of each
//! subject's list) and the traceback pass (`traceback_pass`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1427 (the loop of
//! `BLAST_PreliminarySearchEngine` on the subjects)
//! ```c
//!     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//!        if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
//!            continue;
//!        }
//! ```
//! and c++/src/algo/blast/core/blast_traceback.c:1521-1523 (the loop of `BLAST_ComputeTraceback_MT`
//! on the subjects' batches)
//! ```c
//!         for (i = 0; i < batches->num_batches; i++) {
//! ```
//!
//! This module is a child of `run` for the types it keeps private (`SubjectSource`,
//! `PreparedSubjects`, `PrelimHit`, `SubjectHits`, `SubjectScratch`, `SubjectCtx`, the hit-list
//! helpers, `database`). Its methods are the statements that `search_query_batch` held between the
//! batch's set-up and its result, in their order; what they read of the batch is `SubjectPasses`, built
//! once per batch beside `prelim_subject::PrelimBatch` and `traceback_subject::TraceBatch`.

use super::batch_setup::BatchCutoffs;
use super::calculate_blastn_context_statistics;
use super::prelim_subject::{self, PrelimBatch};
use super::traceback_subject::{self, TraceBatch};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use super::X_BNSTREAM_SUBJECTS_PER_THREAD;
use super::{
    culling, database, listed_prelim_hits, prelim_hits_by_subject, ungapped_reap_query_coverage,
    ungapped_subject_stage, update_hitlists_with_subject_hits, x_bnstream, BlastnArgs,
    BlastnHitList, BlastnSubjectRef, BlastnTiming, FastaRecord, GapAlignScratch, LookupTables,
    PrelimHit, PrelimHitCollector, PrelimHitLists, PreparedSubjects, QueryContext, SubjectCtx,
    SubjectHits, SubjectScratch, SubjectSource,
};
use crate::algorithm::blastn::coordination::subject_trace_id;
use anyhow::Result;

// NCBI reference: ncbi-blast/c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
// ```c
// CMTArgs(...)
// {
// #ifdef NCBI_NO_THREADS
//     m_NumThreads = CThreadable::kMinNumThreads;
//     m_MTMode = eNotSupported;
// #endif
// }
// ```
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

/// The scratch of one worker: the gapped aligner's and the subject search's.
type SubjectState = (GapAlignScratch, SubjectScratch);

/// What the passes over a batch's subjects read of the batch, built once by `search_query_batch`
/// (shared references and `Copy` values, so the struct is `Sync`: the passes run on the pool).
///
/// NCBI's preliminary stage keeps, for each query, the preliminary HSPs of at most
/// `prelim_hitlist_size` subjects (a split batch merges those of its query chunks), and
/// its traceback reads only those.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1674-1783
/// ```c
/// BLAST_PreliminarySearchEngine(program_number, query, query_info,
///    seq_src, gap_align, score_params, lookup_wrap, word_options,
///    ext_params, hit_params, eff_len_params, psi_options,
///    db_options, hsp_stream, diagnostics, interrupt_search,
///    progress_info);
/// ...
/// BLAST_ComputeTraceback(program_number, hsp_stream, query, query_info,
///    seq_src, gap_align, score_params, ext_params, hit_params,
///    eff_len_params, db_options, psi_options, rps_info, pattern_blk,
///    results, interrupt_search, progress_info);
/// ```
/// A database has its own drivers of the two passes (`database::search_blocks`,
/// `database::traceback_listed`), with a scratch per worker made as the subject set's
/// parallel pass makes them.
pub(super) struct SubjectPasses<'a> {
    pub(super) args: &'a BlastnArgs,
    /// The subject set or the database.
    pub(super) subjects: SubjectSource<'a>,
    /// The records of a subject set and their letters (none for a database, which has drivers of its own).
    pub(super) prepared: Option<&'a PreparedSubjects<'a>>,
    pub(super) subject_records: &'a [FastaRecord],
    /// The values the preliminary and the traceback stage read (`search_subject` dispatches).
    pub(super) prelim_batch: &'a PrelimBatch<'a>,
    pub(super) trace_batch: &'a TraceBatch<'a>,
    pub(super) any_valid_context: bool,
    pub(super) queries: &'a [FastaRecord],
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) query_eff_searchsp: &'a Vec<i64>,
    pub(super) batch_cutoffs: &'a BatchCutoffs,
    pub(super) timing: Option<&'a BlastnTiming>,
    /// The sizes of the preliminary and of the final hit lists.
    pub(super) prelim_hitlist_size: usize,
    pub(super) hitlist_size: usize,
}

impl SubjectPasses<'_> {
    /// One subject for one stage: the preliminary search (`prelim_source` is `None`: the subject's
    /// preliminary HSPs, `lookups` its tables) or the traceback of the HSPs the preliminary stage kept
    /// (`Some`: the subject's final HSPs). `state` is the worker's scratch. Nothing is searched without a
    /// valid context.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1469-1470 (the subject's
    /// preliminary search) and c++/src/algo/blast/core/blast_traceback.c:1521,1632-1633 (its traceback)
    /// ```c
    ///       status =
    ///           s_BlastSearchEngineCore(program_number, query, query_info,
    /// ```
    pub(super) fn search_subject<const UNGAPPED: bool>(
        &self,
        s_idx: usize,
        subject: BlastnSubjectRef<'_>,
        state: &mut SubjectState,
        lookups: Option<&LookupTables>,
        prelim_source: Option<&Vec<PrelimHit>>,
    ) -> (Option<SubjectHits>, Vec<PrelimHit>) {
        let Self {
            any_valid_context,
            queries: queries_ref,
            query_contexts: query_contexts_ref,
            query_eff_searchsp: query_eff_searchsp_ref,
            batch_cutoffs: x_batch_cutoffs,
            timing,
            subject_records: trace_subject_records,
            prelim_batch,
            trace_batch,
            ..
        } = *self;
        let (blastn_debug, blastn_trace_enabled) =
            (prelim_batch.blastn_debug, prelim_batch.blastn_trace_enabled);
        let (gap_scratch, subject_scratch) = state;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:177-180
        // ```c
        //     int status = m_PrelimSearch->CheckInternalData();
        //     if (status != 0)
        //     {
        //          // Search was not run, but we send back an empty CSearchResultSet.
        // ```
        // With no valid context, NCBI does not search.
        if !any_valid_context {
            return (None, Vec::new());
        }
        let queries = queries_ref;
        let query_contexts = query_contexts_ref;
        let query_eff_searchsp = query_eff_searchsp_ref;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:493-526
        // ```c
        // if (aux_struct->WordFinder) {
        //     aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
        //                            word_params, aux_struct->ewp,
        //                            aux_struct->offset_pairs,
        //                            kScanSubjectOffsetArraySize,
        //                            init_hitlist, ungapped_stats);
        //     if (init_hitlist->total == 0) continue;
        // }
        // if (aux_struct->GetGappedScore) {
        //     status = aux_struct->GetGappedScore(program_number, query,
        //             query_info, subject, gap_align, score_params, ext_params,
        //             hit_params, word_params, init_hitlist, &hsp_list,
        //             gapped_stats, NULL);
        // }
        // ```
        let timing_ref = timing;
        let timing_enabled = timing_ref.is_some();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
        // ```c
        // hsp_list = Blast_HSPListFree(hsp_list);
        // BlastInitHitListReset(init_hitlist);
        // ```
        let s_len_full = subject.length;
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to. */
        //    BlastHSP** hsp_array;
        //    Int4 hspcnt;
        //    ...
        // } BlastHSPList;
        // ```
        // Made only for the debug and trace lines, which alone read it.
        let s_id_trace = if blastn_debug || blastn_trace_enabled {
            trace_subject_records.get(s_idx).map_or_else(
                || "unknown".to_string(),
                |record| subject_trace_id(record.title()),
            )
        } else {
            String::new()
        };
        let s_id = s_id_trace.as_str();

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:336-377
        // ```c
        // if (!gapped_calculation || sbp->matrix_only_scoring) {
        //     double cutoff_e = s_GetCutoffEvalue(program_number);
        //     Int4 query_length = query_info->contexts[context].query_length;
        //     if (program_number == eBlastTypeBlastn ||
        //         program_number == eBlastTypeMapping)
        //         query_length *= 2;
        //     BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp,
        //                   MIN((Uint8)subj_length,
        //                       (Uint8)query_length) * ((Uint8)subj_length),
        //                   TRUE, gap_decay_rate);
        // } else {
        //     new_cutoff = gap_trigger;
        // }
        // new_cutoff = MIN(new_cutoff,
        //                  hit_params->cutoffs[context].cutoff_score_max);
        // curr_cutoffs->cutoff_score = new_cutoff;
        // ```
        // NCBI uses UNGAPPED params (kbp_std) for gap_trigger calculation
        // and GAPPED params (kbp_gap) for cutoff_score_max calculation
        // The values depend on the context only (the gapped calculation's branch above: gap_trigger
        // capped by cutoff_score_max), so the batch computes them once (`x_batch_cutoffs`).
        let (cutoff_scores, hit_saving_cutoff_scores, x_dropoff_scores, reduced_cutoff_scores) = (
            &x_batch_cutoffs.0,
            &x_batch_cutoffs.1,
            &x_batch_cutoffs.2,
            &x_batch_cutoffs.3,
        );

        // What both stages read of this call (`SubjectCtx`).
        let subject_ctx = SubjectCtx {
            s_idx,
            s_id,
            s_len_full,
            queries,
            query_contexts,
            query_eff_searchsp,
            cutoff_scores,
            hit_saving_cutoff_scores,
            x_dropoff_scores,
            reduced_cutoff_scores,
            timing_ref,
            timing_enabled,
        };
        if let Some(kept) = prelim_source {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1521,1632-1633
            // ```c
            //         for (i = 0; i < batches->num_batches; i++) {
            // ...
            //             /* process all the hits to this subject sequence, one list at a time */
            //             for (hsplist_itr = 0; hsplist_itr < batch->num_hsplists; hsplist_itr++) {
            // ```
            // The traceback of the HSPs that the preliminary stage kept
            // (`traceback_subject::traceback_subject`).
            let mut subject_hits: Option<SubjectHits> = None;
            traceback_subject::traceback_subject::<UNGAPPED>(
                trace_batch,
                &subject_ctx,
                subject,
                kept,
                gap_scratch,
                subject_scratch,
                &mut subject_hits,
            );
            (subject_hits, Vec::new())
        } else {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1469-1470
            // ```c
            //       status =
            //           s_BlastSearchEngineCore(program_number, query, query_info,
            // ```
            // The preliminary search of the subject (`prelim_subject::prelim_subject`).
            // The preliminary stage gives the subject's HSPs to the collector (`search_subjects`).
            (
                None,
                prelim_subject::prelim_subject::<UNGAPPED>(
                    prelim_batch,
                    &subject_ctx,
                    subject,
                    lookups,
                    gap_scratch,
                    subject_scratch,
                ),
            )
        }
    }

    /// One pass over the subjects in subject order: the preliminary stage gives each subject's
    /// preliminary HSPs (`kept` is `None`), the traceback saves the final HSPs of the kept
    /// preliminary HSPs in the final hit lists of `hitlist_size` subjects
    /// (`Blast_HSPResultsInsertHSPList`).
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1427
    /// ```c
    /// itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    /// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    ///        != BLAST_SEQSRC_EOF) {
    ///    if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    ///        continue;
    ///    }
    /// }
    /// ```
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:1700-1707
    /// ```c
    ///                 batch->hsplist_array[hsplist_itr] = NULL;
    ///                 if (hsp_list->hspcnt == 0) {
    ///                     hsp_list = Blast_HSPListFree(hsp_list);
    ///                 }
    ///                 else {
    ///                     Blast_HSPResultsInsertHSPList(thread_data->tld[tid]->results, hsp_list,
    ///                                   hit_params->options->hitlist_size);
    ///                 }
    /// ```
    /// Each subject's results go to `sink` in OID order, which fills the hit lists as NCBI's
    /// `Blast_HitListUpdate` does, one HSP list at a time. One thread hands each subject over as
    /// soon as it is searched, so only the HSPs the hit lists keep stay in memory. Two or more
    /// threads search all subjects of the batch and then hand them over, or, with
    /// `LOSAT_X_BNSTREAM`, hand each subject over in OID order while they search the next ones.
    pub(super) fn search_subjects<const UNGAPPED: bool>(
        &self,
        kept: Option<&[Vec<PrelimHit>]>,
        lookups: Option<&LookupTables>,
        sink: &mut (dyn FnMut(usize, Option<SubjectHits>, Vec<PrelimHit>) + Send),
    ) {
        let use_parallel = self.prelim_batch.use_parallel;
        let new_state = || self.new_scratch();
        let search_one = |state: &mut SubjectState, s_idx: usize| {
            self.search_subject::<UNGAPPED>(
                s_idx,
                self.prepared
                    .expect("the subjects of a subject set")
                    .subject(s_idx),
                state,
                lookups,
                kept.map(|kept| &kept[s_idx]),
            )
        };
        let num_subjects = self.subject_records.len();
        if use_parallel {
            #[cfg(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))]
            {
                if x_bnstream() {
                    crate::utils::xreach::x_reach!("LOSAT_X_BNSTREAM", true);
                    let threads = self.prelim_batch.parallel_pool.threads().max(1);
                    let window = threads * X_BNSTREAM_SUBJECTS_PER_THREAD;
                    self.prelim_batch.parallel_pool.install(|| {
                        x_bnstream::search_in_order(
                            num_subjects,
                            window,
                            false,
                            &new_state,
                            &search_one,
                            &mut |s_idx, (hits, prelim)| sink(s_idx, hits, prelim),
                        )
                    });
                } else {
                    let results: Vec<(Option<SubjectHits>, Vec<PrelimHit>)> =
                        self.prelim_batch.parallel_pool.install(|| {
                            (0..num_subjects)
                                .into_par_iter()
                                .map_init(new_state, search_one)
                                .collect()
                        });
                    for (s_idx, (hits, prelim)) in results.into_iter().enumerate() {
                        sink(s_idx, hits, prelim);
                    }
                }
            }
            #[cfg(any(
                not(feature = "parallel"),
                all(target_arch = "wasm32", not(feature = "wasm-threads"))
            ))]
            {
                unreachable!("use_parallel is false when parallel threads are disabled");
            }
        } else {
            if x_bnstream() {
                crate::utils::xreach::x_reach!("LOSAT_X_BNSTREAM", false);
            }
            let mut state = new_state();
            for s_idx in 0..num_subjects {
                let (hits, prelim) = search_one(&mut state, s_idx);
                sink(s_idx, hits, prelim);
            }
        }
    }

    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:991-1041
    /// ```c
    /// Int4 offset_array_size = GetOffsetArraySize(lookup_wrap);
    /// ...
    /// aux_struct->offset_pairs =
    ///   (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
    /// ```
    fn new_scratch(&self) -> SubjectState {
        (
            GapAlignScratch::new(),
            SubjectScratch::new(
                self.prelim_batch.offset_array_size,
                self.prelim_batch.window_size,
            ),
        )
    }

    /// An ungapped (`-ungapped`) search only (the other instance never calls them): the steps of
    /// NCBI's preliminary stage that finish a subject's HSP list and that `search_subject` leaves
    /// to its caller. The reap by the query coverage is made as the HSPs are given to the
    /// collector (it takes nothing of the subject). A database sequence is read again, unpacked
    /// (BLASTNA), when its HSP list is not empty only, and re-evaluated there.
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1482
    /// ```c
    ///           if (hsp_list && hsp_list->hspcnt > 0) {
    ///              ...  Blast_HSPListReevaluateUngapped(hsp_list, query, subject, ...
    /// ```
    fn x_ungapped_query_coverage(&self, hits: &mut Vec<PrelimHit>) {
        ungapped_reap_query_coverage(
            hits,
            self.args.qcov_hsp_perc.unwrap_or(0.0),
            |context_idx| self.query_contexts[context_idx].len,
        );
    }

    /// The finish of a database sequence's HSP list in an `-ungapped` search: the sequence is read again,
    /// unpacked (BLASTNA), when its HSP list is not empty, and re-evaluated there.
    fn x_ungapped_database_stage(
        &self,
        database: &database::BlastnDatabase,
        ordinal: u32,
        hits: &mut Vec<PrelimHit>,
    ) -> Result<()> {
        let &PrelimBatch {
            encoded_queries_blastna,
            encoded_queries_nomask,
            score_matrix,
            search_karlin_ref,
            round_down_evalue_score,
            evalue_threshold,
            percent_identity,
            ..
        } = self.prelim_batch;
        if hits.is_empty() {
            return Ok(());
        }
        database.with_blastna(ordinal, |s_blastna| {
            ungapped_subject_stage(
                hits,
                s_blastna,
                encoded_queries_blastna,
                encoded_queries_nomask,
                score_matrix,
                &self.batch_cutoffs.0,
                percent_identity,
                |prelim| {
                    let eff_searchsp = self
                        .query_eff_searchsp
                        .get(prelim.context_idx as usize)
                        .copied()
                        .unwrap_or(0);
                    calculate_blastn_context_statistics(
                        prelim.prelim_score,
                        &search_karlin_ref[prelim.context_idx as usize].gapped,
                        eff_searchsp,
                        round_down_evalue_score,
                    )
                    .1
                },
                evalue_threshold,
            );
            self.x_ungapped_query_coverage(hits);
        })
    }

    /// The preliminary pass of the batch: the lists of `split_prelim_lists` (a split batch's query chunks
    /// searched before), else the HSPs of every subject, in the order of their ordinal ids, into the
    /// preliminary hit lists of `prelim_hitlist_size` subjects per query. With `-culling_limit` NCBI's
    /// culling writer takes the place of the HSP collector. A database's read error is returned as it is
    /// (the caller reports it with `database_search_error`).
    pub(super) fn prelim_pass<const UNGAPPED: bool>(
        &self,
        split_prelim_lists: Option<PrelimHitLists>,
        lookup_tables: &LookupTables,
    ) -> Result<PrelimHitLists> {
        if let Some(lists) = split_prelim_lists {
            return Ok(lists);
        }
        // `-culling_limit`: NCBI's culling writer takes the place of the HSP collector
        // (`BlastHSPWriter` of `CreateHspWriter`), for a batch and for a query chunk alike
        // (a chunk's queries are its parts, so the lengths are the chunk's). Chosen once;
        // the default `else` arm below is unchanged.
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/setup_factory.cpp:318-342
        // ```c
        //     else if (filt_opts) {
        //         if (filt_opts->best_hit && (filt_opts->best_hit_stage & ePrelimSearch))
        //         ...
        //         else if (filt_opts->culling_opts &&
        //                  (filt_opts->culling_stage & ePrelimSearch))
        //         {
        // ```
        if self.args.culling_limit_value() > 0 {
            let mut writer = culling::PrelimCulling::new(
                self.prelim_batch.query_lengths,
                self.prelim_hitlist_size,
                self.args.culling_limit_value(),
            );
            self.prelim_into::<UNGAPPED, _>(lookup_tables, &mut |oid, prelim| {
                writer.add(oid, prelim)
            })?;
            Ok(writer.finish())
        } else {
            let mut collector =
                PrelimHitCollector::new(self.queries.len(), self.prelim_hitlist_size);
            self.prelim_into::<UNGAPPED, _>(lookup_tables, &mut |oid, prelim| {
                collector.add(oid, prelim)
            })?;
            Ok(collector.hit_lists)
        }
    }

    /// The preliminary search of every subject, each subject's HSPs (with `-ungapped`, finished) given to
    /// `sink` in the order of the subjects' ids: a subject set's by `search_subjects`, a database's by
    /// `database::search_blocks` (one instance of it for each sink and each `UNGAPPED`). The subject set's
    /// and the database's drivers are the same for the HSP collector and the culling writer.
    fn prelim_into<const UNGAPPED: bool, M: FnMut(usize, Vec<PrelimHit>) + Send>(
        &self,
        lookup_tables: &LookupTables,
        sink: &mut M,
    ) -> Result<()> {
        match self.subjects {
            SubjectSource::Prepared(_) => {
                self.search_subjects::<UNGAPPED>(
                    None,
                    Some(lookup_tables),
                    &mut |oid, _, mut prelim| {
                        if UNGAPPED {
                            self.x_ungapped_query_coverage(&mut prelim);
                        }
                        sink(oid, prelim)
                    },
                );
                Ok(())
            }
            SubjectSource::Database(database) => {
                // The subjects come one at a time in the order of their ordinal ids
                // (`database::search_blocks`), as the subject set's come in OID order.
                // NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:177-180
                // ```c
                //     int status = m_PrelimSearch->CheckInternalData();
                //     if (status != 0)
                //     {
                //          // Search was not run, but we send back an empty CSearchResultSet.
                // ```
                // With no valid context NCBI does not search, and reads no sequence.
                if !self.any_valid_context {
                    return Ok(());
                }
                let mut add =
                    |ordinal: u32, prelim_hits: Vec<PrelimHit>| sink(ordinal as usize, prelim_hits);
                let (parallel_pool, use_parallel) = (
                    self.prelim_batch.parallel_pool,
                    self.prelim_batch.use_parallel,
                );
                let new_scratch = || self.new_scratch();
                if UNGAPPED {
                    // The first error of reading a database sequence again for the
                    // ungapped re-evaluation (`BlastnDatabase::with_blastna`).
                    let x_ungapped_error: std::sync::Mutex<Option<(u32, anyhow::Error)>> =
                        std::sync::Mutex::new(None);
                    // The callback of the database's preliminary pass in an ungapped search: the search of the
                    // sequence as for a default search, then the finish of its HSP list when it has HSPs left
                    // (`x_ungapped_database_stage`). A read error goes to the cell `x_ungapped_error` (the one of the
                    // smallest ordinal stays) and the sequence has no HSPs. The default search has its own
                    // callback (a different closure, so that the default instance of this function holds no
                    // capture of this one).
                    database::search_blocks(
                        database,
                        parallel_pool,
                        use_parallel,
                        &new_scratch,
                        &|state: &mut SubjectState, ordinal: u32, subject: BlastnSubjectRef<'_>| {
                            let (_, mut prelim_hits) = self.search_subject::<UNGAPPED>(
                                ordinal as usize,
                                subject,
                                state,
                                Some(lookup_tables),
                                None,
                            );
                            if let Err(error) =
                                self.x_ungapped_database_stage(database, ordinal, &mut prelim_hits)
                            {
                                let mut first = x_ungapped_error
                                    .lock()
                                    .expect("the error cell is not poisoned");
                                if first
                                    .as_ref()
                                    .is_none_or(|(smallest, _)| ordinal < *smallest)
                                {
                                    *first = Some((ordinal, error));
                                }
                                return None;
                            }
                            (!prelim_hits.is_empty()).then_some(prelim_hits)
                        },
                        &mut add,
                    )
                    .and_then(|()| {
                        match x_ungapped_error
                            .lock()
                            .expect("the error cell is not poisoned")
                            .take()
                        {
                            Some((_, error)) => Err(error),
                            None => Ok(()),
                        }
                    })
                } else {
                    database::search_blocks(
                        database,
                        parallel_pool,
                        use_parallel,
                        &new_scratch,
                        &|state: &mut SubjectState, ordinal: u32, subject: BlastnSubjectRef<'_>| {
                            let (_, prelim_hits) = self.search_subject::<UNGAPPED>(
                                ordinal as usize,
                                subject,
                                state,
                                Some(lookup_tables),
                                None,
                            );
                            (!prelim_hits.is_empty()).then_some(prelim_hits)
                        },
                        &mut add,
                    )
                }
            }
        }
    }

    /// The traceback pass of the batch: the final HSPs of the HSPs the preliminary pass kept
    /// (`prelim_lists`), in the final hit lists of `hitlist_size` subjects per query. A database's read
    /// error is returned as it is (the caller reports it with `database_search_error`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1700-1707
    /// ```c
    ///                 batch->hsplist_array[hsplist_itr] = NULL;
    ///                 if (hsp_list->hspcnt == 0) {
    ///                     hsp_list = Blast_HSPListFree(hsp_list);
    ///                 }
    ///                 else {
    ///                     Blast_HSPResultsInsertHSPList(thread_data->tld[tid]->results, hsp_list,
    ///                                   hit_params->options->hitlist_size);
    ///                 }
    /// ```
    pub(super) fn traceback_pass<const UNGAPPED: bool>(
        &self,
        prelim_lists: PrelimHitLists,
    ) -> Result<Vec<Option<BlastnHitList>>> {
        let (queries_ref, hitlist_size) = (self.queries, self.hitlist_size);
        match self.subjects {
            SubjectSource::Prepared(_) => {
                let kept = prelim_hits_by_subject(prelim_lists, self.subject_records.len());
                let mut hit_lists: Vec<Option<BlastnHitList>> =
                    Vec::with_capacity(queries_ref.len());
                hit_lists.resize_with(queries_ref.len(), || None);
                self.search_subjects::<UNGAPPED>(Some(&kept), None, &mut |_, hits, _| {
                    if let Some(hits) = hits {
                        update_hitlists_with_subject_hits(&mut hit_lists, hits, hitlist_size);
                    }
                });
                Ok(hit_lists)
            }
            SubjectSource::Database(database) => {
                let (parallel_pool, use_parallel) = (
                    self.prelim_batch.parallel_pool,
                    self.prelim_batch.use_parallel,
                );
                let new_scratch = || self.new_scratch();
                let listed = listed_prelim_hits(prelim_lists);
                let found = database::traceback_listed(
                    database,
                    parallel_pool,
                    use_parallel,
                    &listed,
                    &new_scratch,
                    &|state: &mut SubjectState,
                      ordinal: u32,
                      subject: BlastnSubjectRef<'_>,
                      kept: &Vec<PrelimHit>| {
                        let (mut subject_hits, _) = self.search_subject::<UNGAPPED>(
                            ordinal as usize,
                            subject,
                            state,
                            None,
                            Some(kept),
                        );
                        // Every listed sequence's lists are kept until the last one is traced
                        // (`found` below), so each keeps only its length (`SubjectHits::fit`). A
                        // subject set's lists go to the hit lists at once (`search_subjects` above),
                        // whose copies and `update_hitlists_with_subject_hits` hold their counts.
                        if let Some(hits) = subject_hits.as_mut() {
                            hits.fit();
                        }
                        subject_hits
                    },
                    &|kept: &Vec<PrelimHit>| {
                        kept.iter()
                            .map(|hit| (hit.prelim_ss, hit.prelim_se))
                            .collect()
                    },
                )?;
                // The final HSPs of the listed sequences in the order of their ordinal ids, as
                // the subject set's traceback saves those of its subjects
                // (`Blast_HSPResultsInsertHSPList`, `search_subjects`).
                let mut hit_lists: Vec<Option<BlastnHitList>> =
                    Vec::with_capacity(queries_ref.len());
                hit_lists.resize_with(queries_ref.len(), || None);
                for hits in found {
                    update_hitlists_with_subject_hits(&mut hit_lists, hits, hitlist_size);
                }
                Ok(hit_lists)
            }
        }
    }
}
