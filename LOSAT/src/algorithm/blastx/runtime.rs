//! Native serial BLASTX runtime returning actual final HSPs and ancillary state.
use super::{
    args::ResolvedOptions,
    input::FastaRecord,
    kappa::KappaState,
    parameters::{subject_parameters, ContextParameters},
    preliminary::PreliminaryHsp,
    query_setup::PreparedQueryBatch,
    results::{
        prelim_hitlist_size, subject_besthit, CompoHeap, HeapEntry, HitList, Hsp, HspList,
        HspStream,
    },
    search::{search_core, PreliminarySink},
    split::QueryChunk,
    statistics::{
        evaluate_gapped, identity_and_positive, normalize_scores, reap_final, reap_preliminary,
    },
    traceback::ordinary_traceback,
};
use crate::utils::matrix::aa_char_to_ncbistdaa;
use anyhow::Result;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1777
// ```c++
//     if(hit_params->options->query_cov_hsp_perc > 0 || hit_params->options->max_hsps_per_subject > 0 ||
//        (hit_params->options->hsp_filt_opt != NULL && hit_params->options->hsp_filt_opt->subject_besthit_opts != NULL)) {
//     	s_FilterBlastResults(results, hit_params->options, query_info, program_number);
//     }
//
//     /* Re-sort the hit lists according to their best e-values, because they
//        could have changed. Only do this for a database search. */
//     if (BlastSeqSrcGetTotLen(seq_src) > 0) {
//         Blast_HSPResultsSortByEvalue(results);
//     }
//
//
//     /* Eliminate extra hits from results, if preliminary hit list size is
//        larger than the final hit list size */
//     s_BlastPruneExtraHits(results, hit_params->options->hitlist_size);
// ```
pub struct BatchResults {
    pub query_ordinal: usize,
    pub queries: Vec<Vec<HspList>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-101
    // ```c++
    //     // find the first valid context corresponding to this query
    //     for (i = 0; i < context_per_query; i++) {
    //         BlastContextInfo *ctx = query_info->contexts +
    //                                 query_number * context_per_query + i;
    //         if (ctx->is_valid) {
    //             m_SearchSpace = ctx->eff_searchsp;
    // 	    m_LengthAdjustment = ctx->length_adjustment;
    //             break;
    //         }
    //     }
    //     if (i >= context_per_query) {
    //         return; // we didn't find a valid context :(
    //     }
    //
    //     // fill in the Karlin blocks for that context, if they
    //     // are valid
    //     const int ctx_index = query_number * context_per_query + i;
    //     if (sbp->kbp_std) {
    //         s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
    //     }
    // ```
    pub prepared: PreparedQueryBatch,
    pub parameters: Vec<ContextParameters>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-181,202-210
    // ```c++
    //
    //     int status = m_PrelimSearch->CheckInternalData();
    //     if (status != 0)
    //     {
    //          // Search was not run, but we send back an empty CSearchResultSet.
    //          CRef<ILocalQueryData> local_query_data = m_QueryFactory->MakeLocalQueryData(m_Opts);
    //               msg_vec.push_back(q_msg);
    //               seqid_vec.push_back(query_id);
    //               CRef<objects::CSeq_align_set> tmp_align;
    //               sa_vec.push_back(tmp_align);
    //               pair<double, double> tmp_pair(-1.0, -1.0);
    //               CRef<CBlastAncillaryData>  tmp_ancillary_data(new CBlastAncillaryData(tmp_pair, tmp_pair, tmp_pair, 0));
    //               ancill_vec.push_back(tmp_ancillary_data);
    //
    //               for(unsigned int i =1; i < num_subjects; i++)
    // ```
    pub search_skipped: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1539-1554
// ```c++
//             /* Calculate and fill the bit scores, since there will be no
//                traceback stage where this can be done. */
//             Blast_HSPListGetBitScores(hsp_list, gapped_calculation, sbp);
//          }
//
//          // This should only happen for sra searches
//          if((seq_arg.seq->bases_offset > 0) && (gapped_calculation))
//          {
//         	 if (Blast_SubjectIsTranslated(program_number))
//         		 s_AdjustSubjectForTranslatedSraSearch(hsp_list, seq_arg.seq->bases_offset, seq_arg.seq->length);
//         	 else
//         		 s_AdjustSubjectForSraSearch(hsp_list, seq_arg.seq->bases_offset);
//          }
//
//          /* Save the results. */
//          status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3647
// ```c++
//                                 Blast_RedoOneMatch(
//                                         alignments,             // thread-local
//                                         redo_align_params,      // thread-local
//                                         incoming_aligns,        // thread-local
//                                         numAligns[frame_index], // local
//                                         kbp->Lambda,            // thread-local
//                                         &matchingSeq,           // thread-local
//                                         -1,                     // const
//                                         query_info,             // thread-local
//                                         numContexts,            // thread-local
//                                         matrix,                 // thread-local
//                                         BLASTAA_SIZE,           // const
//                                         NRrecord,               // thread-local
//                                         &pvalueForThisPair,     // local
//                                         compositionTestIndex,   // thread-local
//                                         &LambdaRatio            // local
//                                 );
//                     }
//
//                     if (*pStatusCode != 0) {
//                         goto match_loop_cleanup;
//                     }
//
// ```
pub trait RuntimeObserver: super::results::ResultsObserver + Send {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:292-301
    // ```c++
    //         // Restore the full query sequence for the traceback stage!
    //         if (m_InternalData->m_Queries == NULL) {
    //             CRef<ILocalQueryData> query_data
    //                 (m_QueryFactory->MakeLocalQueryData(&*m_Options));
    //             // Query masking info is calculated as a side-effect
    //             CBlastScoreBlk sbp
    //                 (CSetupFactory::CreateScoreBlock(opts_memento.get(), query_data,
    //                                         NULL, m_Messages, NULL, NULL));
    //             m_InternalData->m_Queries = query_data->GetSequenceBlk();
    //         }
    // ```
    fn restored(&mut self, _batch: &PreparedQueryBatch, _params: &[ContextParameters]) {}

    fn matrix(
        &mut self,
        _oid: usize,
        _event: crate::core::composition_adjustment::redo_alignment::BlastRedoTraceEvent<'_>,
    ) {
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
    // ```c++
    // BlastCompo_HeapWouldInsert(BlastCompo_Heap * self,
    //                            double eValue,
    //                            int score,
    //                            int subject_index)
    // {
    //     if (self->n < self->heapThreshold ||
    //         eValue <= self->ecutoff ||
    //         eValue <  self->worstEvalue) {
    //         return TRUE;
    //     } else {
    //         /* self is either currently a heap, or must be converted to
    //          * one; use s_CompoHeapRecordCompare to compare against
    //          * the worst element in the heap */
    //         BlastCompo_HeapRecord heapRecord; /* temporary record to
    //                                              compare against */
    //         if (self->heapArray == NULL) s_ConvertToHeap(self);
    //
    //         heapRecord.bestEvalue       = eValue;
    //         heapRecord.bestScore        = score;
    //         heapRecord.subject_index    = subject_index;
    //         heapRecord.theseAlignments  = NULL;
    //
    //         return s_CompoHeapRecordCompare(&self->heapArray[1], &heapRecord);
    //     }
    // ```
    #[allow(clippy::too_many_arguments)]
    fn heap(
        &mut self,
        _query: usize,
        _stage: &str,
        _evalue: f64,
        _score: i32,
        _oid: i64,
        _decision: i64,
        _heap: &CompoHeap,
    ) {
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:757-764
    // ```c++
    //       gapped_calculation = FALSE;
    //
    //    if (options->do_sum_stats && gapped_calculation && avg_subj_length <= 0)
    //        return 1;
    //
    //
    //    /* If parameters have not yet been created, allocate and fill all
    //       parameters that are constant throughout the search */
    // ```
    fn link_setup(
        &mut self,
        _subject_length: i32,
        _db_length: i64,
        _state: ([i32; 2], f64),
        _link: &super::linking::LinkParameters,
    ) {
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
    // ```c++
    // {
    //     Int4 index;
    //
    //     if (!hsp_list || hsp_list->hspcnt == 0)
    //         return 0;
    //
    //     ASSERT(link_hsp_params);
    //
    //     /* Remove any information on number of linked HSPs from previous
    //        linking. */
    //     for (index = 0; index < hsp_list->hspcnt; ++index)
    //         hsp_list->hsp_array[index]->num = 1;
    //
    //     /* Link up the HSP's for this hsp_list. */
    //     if (link_hsp_params->longest_intron <= 0) {
    //         s_BlastEvenGapLinkHSPs(program_number, hsp_list, query_info,
    //                               subject_length, sbp, link_hsp_params,
    //                               gapped_calculation);
    //         /* The HSP's may be in a different order than they were before,
    //            but hsp contains the first one. */
    //     } else {
    //         Blast_HSPListAdjustOddBlastnScores(hsp_list, gapped_calculation, sbp);
    //         /* Calculate individual HSP e-values first - they'll be needed to
    //            compare with sum e-values. Use decay rate to compensate for
    //            multiple tests. */
    // ```
    fn link_run(
        &mut self,
        _subject_length: i32,
        _gapped: bool,
        _link: &super::linking::LinkParameters,
    ) {
    }
    fn early(&mut self, _evalue: f64, _queries: usize, _decision: bool) {}
    // No NCBI counterpart: tells the search whether any observer callback records something; it
    // does not change any value NCBI computes.
    /// EXPERIMENT (LOSAT_X_BXPAR): true when no callback records anything, so
    /// that work may be done out of order (and some of it twice).
    fn x_passive(&self) -> bool {
        false
    }
}
struct Noop;
impl super::results::ResultsObserver for Noop {}
impl RuntimeObserver for Noop {
    fn x_passive(&self) -> bool {
        true
    }
}

// No NCBI counterpart: reads the LOSAT_X_BXLEAN switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXLEAN): leave out work whose result nothing reads: the
/// copy of every candidate HSP kept for an observer when there is none, and
/// the tree and DP scratch of a (chunk, subject) pair without an initial HSP.
pub(crate) fn x_bx_lean() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BXLEAN").is_some())
}

// No NCBI counterpart: reads the LOSAT_X_BXPAR switch once; it does not change any value NCBI
// computes.
// EXPERIMENT (LOSAT_X_BXPAR): parallel BLASTX stages.
pub(crate) fn x_bx_parallel() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BXPAR").is_some())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334,3499-3502
// ```c
//         if ((int) compo_adjust_mode > 1 && !positionBased) {
//             NRrecord_tld[i] = Blast_CompositionWorkspaceNew();
//             status_code = Blast_CompositionWorkspaceInit(
//                     NRrecord_tld[i],
//                     scoringParams->options->matrix
//             );
// ...
//                 NRrecord             = NRrecord_tld[tid];
//                 sbp                  = sbp_tld[tid];
//                 redo_align_params    = redo_align_params_tld[tid];
//                 matrix               = matrix_tld[tid];
// ```
// One XKappaWorker is the Rust counterpart of the per-thread state NCBI keeps in NRrecord_tld,
// redo_align_params_tld and matrix_tld. It is used by one thread at a time.
// One pool slot's Kappa state.  SAFETY of `Send`: the only non-Send field is the
// `context` cell of the gapping parameters, which `redo_context_observed` fills
// for one synchronous call and restores before returning; the state is only
// reached through its mutex.
#[cfg(feature = "parallel")]
struct XKappaWorker(KappaState);
#[cfg(feature = "parallel")]
unsafe impl Send for XKappaWorker {}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1232-1234
// ```c
//                 if (compo_adjust_mode != eNoCompositionBasedStats &&
//                         (subject_is_translated || hsp_index == 0
//                                 || (nearIdenticalStatus != oldNearIdenticalStatus))) {
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3626
// ```c
//                                 Blast_RedoOneMatch(
//                                         alignments,             // thread-local
// ```
// x_speculate_match redoes one match with the same steps. Blast_AdjustScores runs for a protein
// subject only when hsp_index == 0 (line 1233), so a match can read the matrix its predecessor left
// behind. X_REDO_PROBE records whether that happened; if it did the result is dropped and the match
// is redone in stream order.
/// EXPERIMENT (LOSAT_X_BXPAR): redo one match without knowing the matrix its
/// predecessor leaves behind.  `None` when the match would have used that
/// matrix (it is then redone in order); otherwise the finished list, its best
/// score, and the matrix this match leaves (`None` if it leaves the old one).
#[cfg(feature = "parallel")]
#[allow(clippy::type_complexity)]
fn x_speculate_match(
    worker: &mut KappaState,
    list: &HspList,
    batch: &PreparedQueryBatch,
    full: &[ContextParameters],
    options: &ResolvedOptions,
    subjects: &[FastaRecord],
    db_length: i64,
) -> Option<(
    HspList,
    i32,
    Option<crate::core::composition_adjustment::adjust_scores::AdjustedProteinMatrix>,
)> {
    use crate::core::composition_adjustment::redo_alignment::X_REDO_PROBE;
    worker.x_set_matrix(None);
    X_REDO_PROBE.with(|probe| probe.set(0));
    let encoded: Vec<_> = subjects[list.oid]
        .sequence
        .iter()
        .copied()
        .map(aa_char_to_ncbistdaa)
        .collect();
    let input: Vec<_> = list.hsps.iter().map(|h| h.hsp.clone()).collect();
    let redone = worker
        .redo_list_observed(&input, batch, &encoded, &mut |_| {})
        .ok()?;
    let probe = X_REDO_PROBE.with(|probe| probe.get());
    if probe & 2 != 0 {
        return None;
    }
    let (output, best_score) = finalize_redone(
        &redone,
        batch,
        full,
        options,
        &encoded,
        db_length,
        list.query_index,
        list.oid,
        &mut Noop,
    )
    .ok()?;
    let matrix = if probe & 1 != 0 {
        worker.x_take_matrix()
    } else {
        None
    };
    Some((output, best_score, matrix))
}
pub fn search_internal(
    records: &[FastaRecord],
    subjects: &[FastaRecord],
    options: &ResolvedOptions,
) -> Result<Vec<BatchResults>> {
    search_internal_observed(records, subjects, options, &mut Noop)
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// ```c
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
// ```
// NCBI builds a CLocalBlast, and so its worker threads, for every query batch. LOSAT_X_BXPOOL keeps
// one pool for all batches. Reuse of threads only.
/// EXPERIMENT (LOSAT_X_BXPOOL): `LOSAT_X_BXPOOL` set.
pub(crate) fn x_bx_pool() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BXPOOL").is_some())
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-262,289-291
// ```c
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
// ...
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
// ...
//             	ITERATE(CSearchResultSet, result, *results) {
//                	    formatter.PrintOneResultSet(**result, query_batch);
//             	}
// ```
// NCBI searches one query batch and prints its results before it takes the next. With
// LOSAT_X_BXBATCH several small batches are searched together and printed in input order
// afterwards.
/// EXPERIMENT (LOSAT_X_BXBATCH): `LOSAT_X_BXBATCH` set (needs LOSAT_X_BXPOOL).
pub(crate) fn x_bx_batch() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BXBATCH").is_some())
}

thread_local! {
    // No NCBI counterpart: flag for the batched search (a thread that runs one whole batch must not
    // start parallel work of its own); it does not change any value NCBI computes.
    // EXPERIMENT (LOSAT_X_BXBATCH): set while this thread runs a search that
    // must not start parallel work of its own.
    static X_INNER_SERIAL: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

// No NCBI counterpart: reads the flag above; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXBATCH): true inside `x_search_internal_serial`.
pub(crate) fn x_inner_serial() -> bool {
    X_INNER_SERIAL.with(std::cell::Cell::get)
}

// No NCBI counterpart: scope guard that sets and restores the flag above; it does not change any
// value NCBI computes.
struct XInnerSerial(bool);

impl XInnerSerial {
    fn enter() -> Self {
        Self(X_INNER_SERIAL.with(|flag| flag.replace(true)))
    }
}

impl Drop for XInnerSerial {
    fn drop(&mut self) {
        X_INNER_SERIAL.with(|flag| flag.set(self.0));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:256-265
// ```c
//                 if (IsMultiThreaded()) {
//                      x_LaunchMultiThreadedSearch(*chunk_data);
//                 } else {
//                     retval =
//                         CPrelimSearchRunner(*chunk_data, opts_memento.get())();
// ```
// NCBI runs the single-threaded search (CPrelimSearchRunner) when IsMultiThreaded() is false. This
// is search_internal run on one thread, whatever thread that is, so that several batches can run
// side by side.
/// EXPERIMENT (LOSAT_X_BXBATCH): `search_internal` as the one-thread search,
/// whatever thread it is called on. The caller runs several of these side by
/// side; each takes the serial path of every stage.
pub fn x_search_internal_serial(
    records: &[FastaRecord],
    subjects: &[FastaRecord],
    options: &ResolvedOptions,
) -> Result<Vec<BatchResults>> {
    let _serial = XInnerSerial::enter();
    x_search_internal_in_pool(
        records,
        subjects,
        options,
        &crate::utils::threading::SearchPool::x_serial(),
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// ```c
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
// ```
// The search of one query batch (lcl_blast.Run()) on a pool that the caller keeps. The search
// itself is search_internal.
/// EXPERIMENT (LOSAT_X_BXPOOL): `search_internal` on a pool that the caller
/// keeps for all its query batches.
pub fn x_search_internal_in_pool(
    records: &[FastaRecord],
    subjects: &[FastaRecord],
    options: &ResolvedOptions,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Result<Vec<BatchResults>> {
    let mut observer = Noop;
    let mut pipeline = Runtime {
        observer: &mut observer,
        ungapped_link_state: None,
        options,
        subjects,
        db_length: subjects.iter().map(|s| s.sequence.len() as i64).sum(),
        ordinal: 0,
        stream: None,
        chunk: None,
        output: Vec::new(),
    };
    search_core(records, subjects, options, false, &mut pipeline, pool)?;
    Ok(pipeline.output)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3647
// ```c++
//                                 Blast_RedoOneMatch(
//                                         alignments,             // thread-local
//                                         redo_align_params,      // thread-local
//                                         incoming_aligns,        // thread-local
//                                         numAligns[frame_index], // local
//                                         kbp->Lambda,            // thread-local
//                                         &matchingSeq,           // thread-local
//                                         -1,                     // const
//                                         query_info,             // thread-local
//                                         numContexts,            // thread-local
//                                         matrix,                 // thread-local
//                                         BLASTAA_SIZE,           // const
//                                         NRrecord,               // thread-local
//                                         &pvalueForThisPair,     // local
//                                         compositionTestIndex,   // thread-local
//                                         &LambdaRatio            // local
//                                 );
//                     }
//
//                     if (*pStatusCode != 0) {
//                         goto match_loop_cleanup;
//                     }
//
// ```
pub fn search_internal_observed(
    records: &[FastaRecord],
    subjects: &[FastaRecord],
    options: &ResolvedOptions,
    observer: &mut dyn RuntimeObserver,
) -> Result<Vec<BatchResults>> {
    let mut pipeline = Runtime {
        observer,
        ungapped_link_state: None,
        options,
        subjects,
        db_length: subjects.iter().map(|s| s.sequence.len() as i64).sum(),
        ordinal: 0,
        stream: None,
        chunk: None,
        output: Vec::new(),
    };
    crate::utils::threading::with_search_pool(options.num_threads as usize, "BLASTX", |pool| {
        search_core(records, subjects, options, false, &mut pipeline, pool)
    })?;
    Ok(pipeline.output)
}
struct Runtime<'a> {
    ungapped_link_state: Option<([i32; 2], f64)>,
    observer: &'a mut dyn RuntimeObserver,
    options: &'a ResolvedOptions,
    subjects: &'a [FastaRecord],
    db_length: i64,
    ordinal: usize,
    stream: Option<HspStream>,
    chunk: Option<HspStream>,
    output: Vec<BatchResults>,
}
impl PreliminarySink for Runtime<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-181,202-210
    // ```c++
    //
    //     int status = m_PrelimSearch->CheckInternalData();
    //     if (status != 0)
    //     {
    //          // Search was not run, but we send back an empty CSearchResultSet.
    //          CRef<ILocalQueryData> local_query_data = m_QueryFactory->MakeLocalQueryData(m_Opts);
    //               msg_vec.push_back(q_msg);
    //               seqid_vec.push_back(query_id);
    //               CRef<objects::CSeq_align_set> tmp_align;
    //               sa_vec.push_back(tmp_align);
    //               pair<double, double> tmp_pair(-1.0, -1.0);
    //               CRef<CBlastAncillaryData>  tmp_ancillary_data(new CBlastAncillaryData(tmp_pair, tmp_pair, tmp_pair, 0));
    //               ancill_vec.push_back(tmp_ancillary_data);
    //
    //               for(unsigned int i =1; i < num_subjects; i++)
    // ```
    fn unsearched_batch(
        &mut self,
        ordinal: usize,
        batch: &PreparedQueryBatch,
        params: &[ContextParameters],
    ) -> Result<()> {
        self.output.push(BatchResults {
            query_ordinal: ordinal,
            queries: vec![Vec::new(); batch.original_lengths.len()],
            prepared: batch.clone(),
            parameters: params.to_vec(),
            search_skipped: true,
        });
        Ok(())
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-235
    // ```c++
    //     CEffectiveSearchSpacesMemento eff_memento(m_Options);
    //     SplitQuery_SetEffectiveSearchSpace(m_Options, m_QueryFactory,
    //                                        m_InternalData);
    //     int retval = 0;
    //
    //     unique_ptr<const CBlastOptionsMemento> opts_memento
    //         (m_Options->CreateSnapshot());
    //     BLAST_SequenceBlk* queries = m_InternalData->m_Queries;
    //     LookupTableOptions * lut_options = opts_memento->m_LutOpts;
    //     BlastInitialWordOptions * word_options = opts_memento->m_InitWordOpts;
    //
    //     // Query splitting data structure (used only if applicable)
    //     CRef<SBlastSetupData> setup_data(new SBlastSetupData(m_QueryFactory, m_Options));
    //     CRef<CQuerySplitter> query_splitter = setup_data->m_QuerySplitter;
    //     if (query_splitter->IsQuerySplit()) {
    //
    //         CRef<CSplitQueryBlk> split_query_blk = query_splitter->Split();
    // ```
    fn batch_start(
        &mut self,
        ordinal: usize,
        batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
    ) -> Result<()> {
        self.ordinal = ordinal;
        self.stream = Some(if self.options.culling_limit > 0 {
            HspStream::with_culling(
                batch.contexts.iter().map(|c| c.length as i32).collect(),
                self.options.culling_limit,
                prelim_hitlist_size(self.options),
            )
        } else {
            HspStream::new(
                batch.original_lengths.len(),
                prelim_hitlist_size(self.options),
            )
        });
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:237-244
    // ```c++
    //         for (Uint4 i = 0; i < query_splitter->GetNumberOfChunks(); i++) {
    //             try {
    //                 CRef<IQueryFactory> chunk_qf =
    //                     query_splitter->GetQueryFactoryForChunk(i);
    //                 _TRACE("Query chunk " << i << "/" <<
    //                        query_splitter->GetNumberOfChunks());
    //                 CRef<SInternalData> chunk_data =
    //                     SplitQuery_CreateChunkData(chunk_qf, m_Options,
    // ```
    fn chunk_start(&mut self, chunk: &QueryChunk) -> Result<()> {
        self.chunk = Some(if self.options.culling_limit > 0 {
            HspStream::with_culling(
                chunk
                    .prepared
                    .contexts
                    .iter()
                    .map(|c| c.length as i32)
                    .collect(),
                self.options.culling_limit,
                prelim_hitlist_size(self.options),
            )
        } else {
            HspStream::new(
                chunk.prepared.original_lengths.len(),
                prelim_hitlist_size(self.options),
            )
        });
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:583-593
    // ```c++
    //                      overlap, score_options->gapped_calculation,
    //                      Blast_ProgramIsMapping(program_number));
    //
    //         if((hit_params->options->hsp_filt_opt != NULL) &&
    //            (hit_params->options->hsp_filt_opt->subject_besthit_opts != NULL)) {
    //            	Blast_HSPListSubjectBestHit(program_number,
    //            								hit_params->options->hsp_filt_opt->subject_besthit_opts,
    //            								query_info, combined_hsp_list);
    //         }
    //
    //     } /* End loop on chunks of subject sequence */
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1450-1458
    // ```c++
    //          searches. */
    //       if (hit_params->link_hsp_params && !kNucleotide &&
    //           !gapped_calculation) {
    //           CalculateLinkHSPCutoffs(program_number, query_info, sbp,
    //             hit_params->link_hsp_params, word_params, db_length,
    //             seq_arg.seq->length);
    //       }
    //
    //       if (Blast_SubjectIsTranslated(program_number)) {
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:550-552
    // ```c++
    //             Blast_HSPListAdjustOddBlastnScores(hsp_list, score_options->gapped_calculation, gap_align->sbp);
    // #endif
    //
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3916-3921
    // ```c++
    //
    //       /* use priate interval tree when recomputing alignments */
    //       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
    //                                         hit_options->min_diag_separation))
    //       {
    //          BlastHSP* new_hsp;
    // ```
    fn containment(&mut self, h: &PreliminaryHsp, contained: bool) {
        self.observer.containment(h, contained);
    }
    // No NCBI counterpart: tells the search whether the containment callback is read; it does not
    // change any value NCBI computes.
    fn x_wants_containment(&self) -> bool {
        !(x_bx_lean() && self.observer.x_passive())
    }
    fn preliminary_purge(
        &mut self,
        oid: usize,
        before: &[PreliminaryHsp],
        after: &[PreliminaryHsp],
    ) -> Result<()> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:550-568
        // ```c++
        //             Blast_HSPListAdjustOddBlastnScores(hsp_list, score_options->gapped_calculation, gap_align->sbp);
        // #endif
        //
        //         }
        //
        //         Blast_HSPListSortByScore(hsp_list);
        //
        //
        //         if (score_options->is_ooframe && kTranslatedSubject)
        //             subject->length = prot_length;
        //         } else {
        //             BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
        //                     hit_params->options, &hsp_list);
        //         }
        //
        //         if (hsp_list->hspcnt == 0) continue;
        //
        //         /* The subject ordinal id is not yet filled in this HSP list */
        //         hsp_list->oid = subject->oid;
        // ```
        self.observer.preliminary_call(oid);
        // The score-only list is newly allocated with oid=0; the engine sets
        // its real subject OID only after endpoint purge/score sort.
        self.observer.numeric("PRE_ENDPOINT_IN", 0, before);
        self.observer.numeric("PRE_ENDPOINT_OUT", 0, after);
        Ok(())
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1451-1455
    // ```c++
    //       if (hit_params->link_hsp_params && !kNucleotide &&
    //           !gapped_calculation) {
    //           CalculateLinkHSPCutoffs(program_number, query_info, sbp,
    //             hit_params->link_hsp_params, word_params, db_length,
    //             seq_arg.seq->length);
    // ```
    // search::compute_subject has performed this exact preparation before
    // WordFinder with worker-local state. Adopt it without recomputation,
    // retaining serial observer and linking state in input OID order.
    fn subject_start(
        &mut self,
        oid: usize,
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
        prepared: Option<([i32; 2], f64)>,
    ) -> Result<()> {
        self.ungapped_link_state = prepared;
        if let Some(prepared) = prepared {
            let length = self.subjects[oid].sequence.len() as i32;
            let link = super::linking::link_parameters(self.options)
                .expect("ungapped sum-stat parameters");
            self.observer
                .link_setup(length, self.db_length, prepared, &link);
        }
        Ok(())
    }
    fn subject(
        &mut self,
        oid: usize,
        input: &[PreliminaryHsp],
        batch: &PreparedQueryBatch,
        params: &[ContextParameters],
    ) -> Result<()> {
        let mut hsps: Vec<_> = input
            .iter()
            .map(|h| Hsp {
                hsp: h.clone(),
                num: 0,
                evalue: 0.0,
                bit_score: 0.0,
                identity: 0,
                positive: 0,
                edit_script: Vec::new(),
                composition_method: 0,
            })
            .collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2538-2588
        // ```c++
        // Blast_HSPListSubjectBestHit(EBlastProgramType program,
        // 		                    const BlastHSPSubjectBestHitOptions* subject_besthit_opts,
        // 		                    const BlastQueryInfo *query_info,
        //                             BlastHSPList* hsp_list)
        // {
        //    BlastHSP** hsp_array;  /* hsp_array to purge. */
        //    const int range_diff = subject_besthit_opts->max_range_diff;
        //    Boolean isBlastn = (program == eBlastTypeBlastn);
        //    unsigned int i, j;
        //    int o, e;
        //    int curr_context, target_context;
        //    int qlen;
        //
        //    /* If HSP list is empty, return immediately. */
        //    if (hsp_list == NULL || hsp_list->hspcnt == 0)
        //        return 0;
        //
        //    if (Blast_ProgramIsPhiBlast(program))
        //        return hsp_list->hspcnt;
        //
        //    hsp_array = hsp_list->hsp_array;
        //
        //    // The hsp list is sorted by score
        //    for(i=0; i < hsp_list->hspcnt -1; i++) {
        // 	  if(hsp_array[i] == NULL){
        // 		  continue;
        // 	  }
        //       j = 1;
        //       o = hsp_array[i]->query.offset - range_diff;
        //       e = hsp_array[i]->query.end + range_diff;
        //       if (o < 0) o = 0;
        //       if (e < 0) e = hsp_array[i]->query.end;
        //       while (i+j < hsp_list->hspcnt) {
        //           if (hsp_array[i+j] && hsp_array[i]->context == hsp_array[i+j]->context &&
        //               ((hsp_array[i+j]->query.offset >= o) &&
        //                (hsp_array[i+j]->query.end <= e))){
        //        	      hsp_array[i+j] = Blast_HSPFree(hsp_array[i+j]);
        //           }
        //           j++;
        //       }
        //    }
        //
        //    Blast_HSPListPurgeNullHSPs(hsp_list);
        //
        //    if(isBlastn) {
        // 	   for(i=0; i < hsp_list->hspcnt -1; i++) {
        // 	   	  if(hsp_array[i] == NULL){
        // 	   		  continue;
        // 	   	  }
        // 	   	  // Flip query offsets of current hsp to target context frame
        // 	   	  j = 1;
        // ```
        if self.options.subject_besthit && !hsps.is_empty() {
            self.observer.hsps("SUBJECT_BEST_IN", oid, &hsps);
            subject_besthit(&mut hsps);
            self.observer.hsps("SUBJECT_BEST_OUT", oid, &hsps);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-899
        // ```c++
        //     if (hit_params->link_hsp_params) {
        //         status = BLAST_LinkHsps(program_number, hsp_list_out, query_info,
        //                   subject->length, gap_align->sbp, hit_params->link_hsp_params,
        //                   score_options->gapped_calculation);
        //     } else if (!Blast_ProgramIsPhiBlast(program_number)
        //            && !(isRPS && !sbp->gbp)
        //            /* do not calculate E-values for mapping */
        //            && program_number != eBlastTypeMapping ) {
        //         /* Calculate e-values for all HSPs. Skip this step
        //            for PHI or RPS with old FSC, since calculating the E values
        //            requires precomputation that has not been done yet */
        //         double scale_factor = 1.0;
        //         if (isRPS) {
        //             scale_factor = score_params->scale_factor;
        //         }
        //         Blast_HSPListGetEvalues(program_number, query_info,
        //                                          stat_length, hsp_list_out,
        //                                          score_options->gapped_calculation,
        //                                          isRPS, gap_align->sbp, 0, scale_factor);
        //     }
        //
        //    /* Use score threshold rather than evalue if
        //     * matrix_only_scoring is used.  -RMH-
        //     */
        //     if ( sbp->matrix_only_scoring )
        //     {
        //         status = Blast_HSPListReapByRawScore(hsp_list_out, hit_options);
        //     }else {
        //        /* Discard HSPs that don't pass the e-value test. */
        //         status = s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
        // ```
        let input: Vec<_> = hsps.iter().map(|h| h.hsp.clone()).collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1490-1542
        // ```c++
        //             /* The following must be performed for any ungapped
        //                search with a nucleotide database. */
        //                status =
        //                   Blast_HSPListReevaluateUngapped(
        //                             program_number, hsp_list, query,
        //                             seq_arg.seq, word_params, hit_params,
        //                             query_info, sbp, score_params, seq_src,
        //                             seq_arg.seq->gen_code_string);
        //                if (status) {
        //                   /* Tell the indexing library that this thread is done with
        //                      preliminary search.
        //                   */
        //                   if( check_index_oid != 0 ) {
        //                     ((T_MB_IdxEndSearchIndication)(
        //                         lookup_wrap->end_search_indication))( last_vol_idx );
        //                   }
        //
        //                   BlastSeqSrcReleaseSequence(seq_src, &seq_arg);
        //                   return status;
        //                }
        //                /* Relink HSPs if sum statistics is used, because scores might
        //                 * have changed after reevaluation with ambiguities, and there
        //                 * will be no traceback stage where relinking is done normally.
        //                 * If sum statistics are not used, just recalculate e-values.
        //                 */
        //                if (hit_params->link_hsp_params) {
        //                    status =
        //                        BLAST_LinkHsps(program_number, hsp_list, query_info,
        //                                       seq_arg.seq->length, sbp,
        //                                       hit_params->link_hsp_params,
        //                                       gapped_calculation);
        //                } else {
        //                   Blast_HSPListGetEvalues(program_number, query_info,
        //                                           stat_length, hsp_list,
        //                                           gapped_calculation, FALSE,
        //                                           sbp, 0, 1.0);
        //                }
        //                /* Use score threshold rather than evalue if
        //                 * matrix_only_scoring is used.  -RMH-
        //                 */
        //                if ( sbp->matrix_only_scoring )
        //                {
        //                    status = Blast_HSPListReapByRawScore(hsp_list,
        //                                           hit_params->options);
        //                }else {
        //         	   status = s_Blast_HSPListReapByPrelimEvalue(hsp_list, hit_params);
        //                }
        //
        //                Blast_HSPListReapByQueryCoverage(hsp_list,hit_params->options, query_info, program_number);
        //             /* Calculate and fill the bit scores, since there will be no
        //                traceback stage where this can be done. */
        //             Blast_HSPListGetBitScores(hsp_list, gapped_calculation, sbp);
        //          }
        // ```
        if !self.options.gapped {
            let encoded: Vec<_> = self.subjects[oid]
                .sequence
                .iter()
                .copied()
                .map(aa_char_to_ncbistdaa)
                .collect();
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
            // ```c++
            // {
            //     Int4 index;
            //
            //     if (!hsp_list || hsp_list->hspcnt == 0)
            //         return 0;
            //
            //     ASSERT(link_hsp_params);
            //
            //     /* Remove any information on number of linked HSPs from previous
            //        linking. */
            //     for (index = 0; index < hsp_list->hspcnt; ++index)
            //         hsp_list->hsp_array[index]->num = 1;
            //
            //     /* Link up the HSP's for this hsp_list. */
            //     if (link_hsp_params->longest_intron <= 0) {
            //         s_BlastEvenGapLinkHSPs(program_number, hsp_list, query_info,
            //                               subject_length, sbp, link_hsp_params,
            //                               gapped_calculation);
            //         /* The HSP's may be in a different order than they were before,
            //            but hsp contains the first one. */
            //     } else {
            //         Blast_HSPListAdjustOddBlastnScores(hsp_list, gapped_calculation, sbp);
            //         /* Calculate individual HSP e-values first - they'll be needed to
            //            compare with sum e-values. Use decay rate to compensate for
            //            multiple tests. */
            // ```
            if !input.is_empty() {
                if let Some(link) = super::linking::link_parameters(self.options) {
                    self.observer
                        .link_run(self.subjects[oid].sequence.len() as i32, false, &link);
                }
            }
            let mut linked = super::statistics::evaluate_ungapped(
                &input,
                batch,
                params,
                self.options,
                encoded.len() as i32,
                self.ungapped_link_state,
            )?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-676
            // ```c++
            // s_Blast_HSPListReapByPrelimEvalue(BlastHSPList* hsp_list, const BlastHitSavingParameters* hit_params)
            // {
            //    BlastHSP* hsp;
            //    BlastHSP** hsp_array;
            //    Int4 hsp_cnt = 0;
            //    Int4 index;
            //    double cutoff;
            //
            //    if (hsp_list == NULL)
            //       return 0;
            //
            //    cutoff = hit_params->prelim_evalue;
            //
            //    hsp_array = hsp_list->hsp_array;
            //    for (index = 0; index < hsp_list->hspcnt; index++) {
            //       hsp = hsp_array[index];
            //
            //       ASSERT(hsp != NULL);
            //
            //       if (hsp->evalue > cutoff) {
            //          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
            //       } else {
            //          if (index > hsp_cnt)
            //             hsp_array[hsp_cnt] = hsp_array[index];
            //          hsp_cnt++;
            //       }
            //    }
            //
            //    hsp_list->hspcnt = hsp_cnt;
            //
            //    return 0;
            // }
            //
            // /** The core of the BLAST search: comparison between the (concatenated)
            // ```
            let reap_cutoff = self.options.evalue
                * if self.options.composition > 0 {
                    5.0
                } else {
                    1.0
                };
            self.observer.reap("PRE_REAP_IN", oid, &linked, reap_cutoff);
            reap_preliminary(&mut linked, self.options);
            self.observer
                .reap("PRE_REAP_OUT", oid, &linked, reap_cutoff);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2708-2736
            // ```c++
            //       if (!delete_hsp) {
            //           const Uint1* query_nomask = query_blk->sequence_nomask +
            //               query_info->contexts[context].query_offset;
            //           Int4 align_length = 0;
            //           Blast_HSPGetNumIdentitiesAndPositives(query_nomask,
            //                    							    subject_start,
            //                                					hsp,
            //                                					score_params->options,
            //                                					&align_length,
            //                                					sbp);
            //
            //            delete_hsp = Blast_HSPTest(hsp, hit_params->options, align_length);
            //       }
            //       if (delete_hsp) { /* This HSP is now below the cutoff */
            //          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
            //          purge = TRUE;
            //       }
            //    }
            //
            //    if (target_t)
            //       target_t = BlastTargetTranslationFree(target_t);
            //
            //    if (purge)
            //       Blast_HSPListPurgeNullHSPs(hsp_list);
            //
            //    /* Sort the HSP array by score (scores may have changed!) */
            //    Blast_HSPListSortByScore(hsp_list);
            //    Blast_HSPListAdjustOddBlastnScores(hsp_list, FALSE, sbp);
            //    return 0;
            // ```
            // BLASTX has protein subjects: ambiguity rescoring is inactive.
            // The accepted option profile has percent_identity/min_hit_length
            // at zero, so Blast_HSPTest cannot delete an HSP here.
            let mut payloads: Vec<_> = linked
                .hsps
                .into_iter()
                .map(|h| {
                    let redone = super::kappa::RedoneHsp {
                        hsp: h.hsp.clone(),
                        edit_script: vec![crate::common::GapEditOp::Sub(
                            (h.hsp.q_end - h.hsp.q_start) as u32,
                        )],
                        composition_method: 0,
                        identity_matrix: None,
                    };
                    let (identity, positive) = identity_and_positive(&redone, batch, &encoded);
                    Hsp {
                        hsp: h.hsp,
                        num: h.num,
                        evalue: h.evalue,
                        bit_score: 0.0,
                        identity,
                        positive,
                        edit_script: Vec::new(),
                        composition_method: 0,
                    }
                })
                .collect();
            payloads.sort_by(|a, b| super::preliminary::compare_score(&a.hsp, &b.hsp));
            let input: Vec<_> = payloads.iter().map(|h| h.hsp.clone()).collect();
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
            // ```c++
            // {
            //     Int4 index;
            //
            //     if (!hsp_list || hsp_list->hspcnt == 0)
            //         return 0;
            //
            //     ASSERT(link_hsp_params);
            //
            //     /* Remove any information on number of linked HSPs from previous
            //        linking. */
            //     for (index = 0; index < hsp_list->hspcnt; ++index)
            //         hsp_list->hsp_array[index]->num = 1;
            //
            //     /* Link up the HSP's for this hsp_list. */
            //     if (link_hsp_params->longest_intron <= 0) {
            //         s_BlastEvenGapLinkHSPs(program_number, hsp_list, query_info,
            //                               subject_length, sbp, link_hsp_params,
            //                               gapped_calculation);
            //         /* The HSP's may be in a different order than they were before,
            //            but hsp contains the first one. */
            //     } else {
            //         Blast_HSPListAdjustOddBlastnScores(hsp_list, gapped_calculation, sbp);
            //         /* Calculate individual HSP e-values first - they'll be needed to
            //            compare with sum e-values. Use decay rate to compensate for
            //            multiple tests. */
            // ```
            if !input.is_empty() {
                if let Some(link) = super::linking::link_parameters(self.options) {
                    self.observer
                        .link_run(self.subjects[oid].sequence.len() as i32, false, &link);
                }
            }
            let mut linked = super::statistics::evaluate_ungapped(
                &input,
                batch,
                params,
                self.options,
                encoded.len() as i32,
                self.ungapped_link_state,
            )?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-676
            // ```c++
            // s_Blast_HSPListReapByPrelimEvalue(BlastHSPList* hsp_list, const BlastHitSavingParameters* hit_params)
            // {
            //    BlastHSP* hsp;
            //    BlastHSP** hsp_array;
            //    Int4 hsp_cnt = 0;
            //    Int4 index;
            //    double cutoff;
            //
            //    if (hsp_list == NULL)
            //       return 0;
            //
            //    cutoff = hit_params->prelim_evalue;
            //
            //    hsp_array = hsp_list->hsp_array;
            //    for (index = 0; index < hsp_list->hspcnt; index++) {
            //       hsp = hsp_array[index];
            //
            //       ASSERT(hsp != NULL);
            //
            //       if (hsp->evalue > cutoff) {
            //          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
            //       } else {
            //          if (index > hsp_cnt)
            //             hsp_array[hsp_cnt] = hsp_array[index];
            //          hsp_cnt++;
            //       }
            //    }
            //
            //    hsp_list->hspcnt = hsp_cnt;
            //
            //    return 0;
            // }
            //
            // /** The core of the BLAST search: comparison between the (concatenated)
            // ```
            let reap_cutoff = self.options.evalue
                * if self.options.composition > 0 {
                    5.0
                } else {
                    1.0
                };
            self.observer.reap("PRE_REAP_IN", oid, &linked, reap_cutoff);
            reap_preliminary(&mut linked, self.options);
            self.observer
                .reap("PRE_REAP_OUT", oid, &linked, reap_cutoff);
            let hsps = linked
                .hsps
                .into_iter()
                .map(|h| {
                    let p = &params[h.hsp.context];
                    let mut payload = payloads[h.source_index].clone();
                    payload.bit_score = (h.hsp.score as f64 * p.ungapped.lambda
                        - p.ungapped.k.ln())
                        / std::f64::consts::LN_2;
                    payload.hsp = h.hsp;
                    payload.num = h.num;
                    payload.evalue = h.evalue;
                    payload
                })
                .collect();
            self.chunk
                .as_mut()
                .expect("chunk stream")
                .write_observed(oid, hsps, self.observer)?;
            return Ok(());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
        // ```c++
        // {
        //     Int4 index;
        //
        //     if (!hsp_list || hsp_list->hspcnt == 0)
        //         return 0;
        //
        //     ASSERT(link_hsp_params);
        //
        //     /* Remove any information on number of linked HSPs from previous
        //        linking. */
        //     for (index = 0; index < hsp_list->hspcnt; ++index)
        //         hsp_list->hsp_array[index]->num = 1;
        //
        //     /* Link up the HSP's for this hsp_list. */
        //     if (link_hsp_params->longest_intron <= 0) {
        //         s_BlastEvenGapLinkHSPs(program_number, hsp_list, query_info,
        //                               subject_length, sbp, link_hsp_params,
        //                               gapped_calculation);
        //         /* The HSP's may be in a different order than they were before,
        //            but hsp contains the first one. */
        //     } else {
        //         Blast_HSPListAdjustOddBlastnScores(hsp_list, gapped_calculation, sbp);
        //         /* Calculate individual HSP e-values first - they'll be needed to
        //            compare with sum e-values. Use decay rate to compensate for
        //            multiple tests. */
        // ```
        if !input.is_empty() {
            if let Some(link) = super::linking::link_parameters(self.options) {
                self.observer
                    .link_run(self.subjects[oid].sequence.len() as i32, true, &link);
            }
        }
        let mut linked = evaluate_gapped(
            &input,
            batch,
            params,
            self.options,
            self.subjects[oid].sequence.len() as i32,
            self.db_length,
            1.0,
        )?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-676
        // ```c++
        // s_Blast_HSPListReapByPrelimEvalue(BlastHSPList* hsp_list, const BlastHitSavingParameters* hit_params)
        // {
        //    BlastHSP* hsp;
        //    BlastHSP** hsp_array;
        //    Int4 hsp_cnt = 0;
        //    Int4 index;
        //    double cutoff;
        //
        //    if (hsp_list == NULL)
        //       return 0;
        //
        //    cutoff = hit_params->prelim_evalue;
        //
        //    hsp_array = hsp_list->hsp_array;
        //    for (index = 0; index < hsp_list->hspcnt; index++) {
        //       hsp = hsp_array[index];
        //
        //       ASSERT(hsp != NULL);
        //
        //       if (hsp->evalue > cutoff) {
        //          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
        //       } else {
        //          if (index > hsp_cnt)
        //             hsp_array[hsp_cnt] = hsp_array[index];
        //          hsp_cnt++;
        //       }
        //    }
        //
        //    hsp_list->hspcnt = hsp_cnt;
        //
        //    return 0;
        // }
        //
        // /** The core of the BLAST search: comparison between the (concatenated)
        // ```
        let reap_cutoff = self.options.evalue
            * if self.options.composition > 0 {
                5.0
            } else {
                1.0
            };
        self.observer.reap("PRE_REAP_IN", oid, &linked, reap_cutoff);
        reap_preliminary(&mut linked, self.options);
        self.observer
            .reap("PRE_REAP_OUT", oid, &linked, reap_cutoff);
        let hsps = linked
            .hsps
            .into_iter()
            .map(|h| {
                let mut payload = hsps[h.source_index].clone();
                payload.num = h.num;
                payload.evalue = h.evalue;
                payload
            })
            .collect();
        self.chunk
            .as_mut()
            .expect("chunk stream")
            .write_observed(oid, hsps, self.observer)?;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:268-271
    // ```c++
    //                 _ASSERT(chunk_data->m_HspStream->GetPointer());
    //                 BlastHSPStreamMerge(split_query_blk->GetCStruct(), i,
    //                                 chunk_data->m_HspStream->GetPointer(),
    //                                 m_InternalData->m_HspStream->GetPointer());
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:307-319
    // ```c++
    //         if (IsMultiThreaded()) {
    //              x_LaunchMultiThreadedSearch(*m_InternalData);
    //         } else {
    //             retval = CPrelimSearchRunner(*m_InternalData, opts_memento.get())();
    //             if (retval) {
    //                 NCBI_THROW(CBlastException, eCoreBlastError,
    //                            BlastErrorCode2String(retval));
    //             }
    //         }
    //     }
    //
    //     return m_InternalData;
    // }
    // ```
    fn chunk_end(&mut self, chunk: &QueryChunk, is_split: bool) -> Result<()> {
        let incoming = self.chunk.take().expect("chunk stream");
        if is_split {
            self.stream.as_mut().expect("batch stream").merge_observed(
                incoming,
                chunk,
                self.observer,
            )?;
        } else {
            self.stream = Some(incoming);
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:292-301
    // ```c++
    //         // Restore the full query sequence for the traceback stage!
    //         if (m_InternalData->m_Queries == NULL) {
    //             CRef<ILocalQueryData> query_data
    //                 (m_QueryFactory->MakeLocalQueryData(&*m_Options));
    //             // Query masking info is calculated as a side-effect
    //             CBlastScoreBlk sbp
    //                 (CSetupFactory::CreateScoreBlock(opts_memento.get(), query_data,
    //                                         NULL, m_Messages, NULL, NULL));
    //             m_InternalData->m_Queries = query_data->GetSequenceBlk();
    //         }
    // ```
    fn batch_end(
        &mut self,
        batch: &PreparedQueryBatch,
        params: &[ContextParameters],
    ) -> Result<()> {
        let mut stream = self.stream.take().expect("batch stream");

        let mut full = params.to_vec();
        subject_parameters(
            batch,
            &mut full,
            self.options,
            self.db_length,
            self.subjects.len() as i64,
            self.subjects
                .iter()
                .map(|s| s.sequence.len())
                .min()
                .expect("subjects"),
            None,
        )?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1842-1846
        // ```c++
        //
        //     }
        //     else {
        //     	BlastHSPStreamClose(hsp_stream);
        //     }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1842-1846
        // ```c++
        //
        //     }
        //     else {
        //     	BlastHSPStreamClose(hsp_stream);
        //     }
        // ```
        stream.close_observed(self.observer);
        self.observer.restored(batch, &full);

        let mut queries: Vec<_> = (0..batch.original_lengths.len())
            .map(|_| HitList::new(self.options.hitlist_size as usize))
            .collect();

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3281-3296
        // ```c++
        //
        //         redoneMatches_tld[i] =
        //                 (BlastCompo_Heap*) calloc(numQueries, sizeof(BlastCompo_Heap));
        //         if (redoneMatches_tld[i] == NULL) {
        //             status_code = -1;
        //             goto function_cleanup;
        //         }
        //         for (query_index = 0; query_index < numQueries; query_index++) {
        //             status_code =
        //                 BlastCompo_HeapInitialize(&redoneMatches_tld[i][query_index],
        //                                           hitParams->options->hitlist_size,
        //                                           inclusion_ethresh);
        //             if (status_code != 0) {
        //                 goto function_cleanup;
        //             }
        //         }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3120-3132
        // ```c++
        //     if (redoneMatches == NULL) {
        //         status_code = -1;
        //         goto function_cleanup;
        //     }
        //     for (query_index = 0;  query_index < numQueries;  query_index++) {
        //         status_code =
        //             BlastCompo_HeapInitialize(&redoneMatches[query_index],
        //                                       hitParams->options->hitlist_size,
        //                                       inclusion_ethresh);
        //         if (status_code != 0) {
        //             goto function_cleanup;
        //         }
        //     }
        // ```
        let mut inactive_heaps: Vec<_> = if self.options.composition == 2 {
            (0..queries.len())
                .map(|_| CompoHeap::new(self.options.hitlist_size as usize))
                .collect()
        } else {
            Vec::new()
        };
        for (q, heap) in inactive_heaps.iter().enumerate() {
            self.observer.heap(q, "INIT", 0.0, 0, -1, 0, heap);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3281-3296
        // ```c++
        //
        //         redoneMatches_tld[i] =
        //                 (BlastCompo_Heap*) calloc(numQueries, sizeof(BlastCompo_Heap));
        //         if (redoneMatches_tld[i] == NULL) {
        //             status_code = -1;
        //             goto function_cleanup;
        //         }
        //         for (query_index = 0; query_index < numQueries; query_index++) {
        //             status_code =
        //                 BlastCompo_HeapInitialize(&redoneMatches_tld[i][query_index],
        //                                           hitParams->options->hitlist_size,
        //                                           inclusion_ethresh);
        //             if (status_code != 0) {
        //                 goto function_cleanup;
        //             }
        //         }
        // ```
        let mut heaps: Vec<_> = if self.options.composition == 2 {
            (0..queries.len())
                .map(|_| CompoHeap::new(self.options.hitlist_size as usize))
                .collect()
        } else {
            Vec::new()
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:414-426
        // ```c++
        // BlastCompo_HeapInitialize(BlastCompo_Heap * self, int heapThreshold,
        //                           double ecutoff)
        // {
        //     self->n             = 0;
        //     self->heapThreshold = heapThreshold;
        //     self->ecutoff       = ecutoff;
        //     self->heapArray     = NULL;
        //     self->capacity      = MIN(HEAP_INITIAL_CAPACITY, heapThreshold);
        //     self->worstEvalue   = 0;
        //     /* Begin life as a list */
        //     self->array = calloc(self->capacity + 1, sizeof(BlastCompo_HeapRecord));
        //
        //     return self->array != NULL ? 0 : -1;
        // ```
        for (q, heap) in heaps.iter().enumerate() {
            self.observer
                .heap(q + inactive_heaps.len(), "INIT", 0.0, 0, -1, 0, heap);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3120-3132
        // ```c++
        //     if (redoneMatches == NULL) {
        //         status_code = -1;
        //         goto function_cleanup;
        //     }
        //     for (query_index = 0;  query_index < numQueries;  query_index++) {
        //         status_code =
        //             BlastCompo_HeapInitialize(&redoneMatches[query_index],
        //                                       hitParams->options->hitlist_size,
        //                                       inclusion_ethresh);
        //         if (status_code != 0) {
        //             goto function_cleanup;
        //         }
        //     }
        // ```

        let mut kappa = if self.options.composition == 2 {
            Some(KappaState::new(batch, &full, self.options)?)
        } else {
            None
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3543
        // ```c++
        //                 if (BlastCompo_EarlyTermination(
        //                         localMatch->best_evalue,
        //                         redoneMatches,
        //                         numQueries
        //                 )) {
        //                     Blast_HSPListFree(localMatch);
        //                     if (seqSrc) {
        //                         continue;
        //                     }
        //                     if(actual_num_threads > 1) {
        // #pragma omp critical(intrpt)
        //                     	interrupt = TRUE;
        // #pragma omp flush(interrupt)
        //                     	continue;
        //                     }
        //                 }
        //
        //                 query_index = localMatch->query_index;
        //                 context_index = query_index * numFrames;
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1688-1707
        // ```c++
        //                         } /* fence_hit */
        //                     }    /* !phi_blast */
        //
        //                 } else {
        //                     /* traceback skipped; compute bit scores for searches
        //                        where the traceback phase is seperated from the
        //                        preliminary search. */
        //                     Blast_HSPListGetBitScores(hsp_list, FALSE, sbp);
        //                 }
        //
        //                 /* Free HSP list if all HSPs have been deleted. */
        //
        //                 batch->hsplist_array[hsplist_itr] = NULL;
        //                 if (hsp_list->hspcnt == 0) {
        //                     hsp_list = Blast_HSPListFree(hsp_list);
        //                 }
        //                 else {
        //                     Blast_HSPResultsInsertHSPList(thread_data->tld[tid]->results, hsp_list,
        //                                   hit_params->options->hitlist_size);
        //                 }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream_mt_utils.c:145-170
        // ```c++
        // int BlastHSPStreamToHSPStreamResultsBatch(BlastHSPStream* hsp_stream,
        //                                           BlastHSPStreamResultsBatchArray** batches)
        // {
        //     BlastHSPStreamResultBatch *batch = NULL;
        //
        //     if (!batches || !hsp_stream) {
        //         return BLASTERR_INVALIDPARAM;
        //     }
        //
        //     *batches =
        //         s_BlastHSPStreamResultsBatchArrayNew(s_BlastHSPStreamCountNumOids(hsp_stream));
        //     if ( !*batches ) {
        //         return BLASTERR_MEMORY;
        //     }
        //
        //     for (batch = Blast_HSPStreamResultBatchInit(hsp_stream->results->num_queries);
        //          BlastHSPStreamBatchRead(hsp_stream, batch) != kBlastHSPStream_Eof;
        //          batch = Blast_HSPStreamResultBatchInit(hsp_stream->results->num_queries)) {
        //         if (s_BlastHSPStreamResultsBatchArrayAppend(*batches, batch) != 0) {
        //             s_BlastHSPStreamResultsBatchArrayReset(*batches);
        //             *batches = BlastHSPStreamResultsBatchArrayFree(*batches);
        //             return BLASTERR_MEMORY;
        //         }
        //     }
        //     batch = Blast_HSPStreamResultBatchFree(batch);
        //     return kBlastHSPStream_Success;
        // ```
        let mut result_batches = Vec::new();
        if !self.options.gapped || self.options.composition == 0 {
            loop {
                let lists = stream.batch_read_observed(self.observer);
                if lists.is_empty() {
                    break;
                }
                result_batches.push(lists);
            }
        }
        if !self.options.gapped {
            for lists in result_batches {
                for mut list in lists {
                    for h in &mut list.hsps {
                        let p = &full[h.hsp.context];
                        h.bit_score = (h.hsp.score as f64 * p.ungapped.lambda - p.ungapped.k.ln())
                            / std::f64::consts::LN_2;
                    }
                    queries[list.query_index].update_observed(list, self.observer);
                }
            }
        } else if let Some(state) = &mut kappa {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3392-3424
            // ```c++
            //         BlastHSPListLinkedList* head = NULL;
            //         BlastHSPListLinkedList* tail = NULL;
            //         /*
            //          * Collect matches from stream into linked list, counting them
            //          * along the way.
            //          */
            //         while (BlastHSPStreamRead(hsp_stream, &localMatch)
            //                 != kBlastHSPStream_Eof) {
            //             BlastHSPListLinkedList* entry =
            //                     (BlastHSPListLinkedList*) calloc(
            //                             1,
            //                             sizeof(BlastHSPListLinkedList)
            //                     );
            //             entry->match = localMatch;
            //             if (head == NULL) {
            //                 head = entry;
            //             } else {
            //                 tail->next = entry;
            //             }
            //             tail = entry;
            //             ++numMatches;
            //         }
            //         /*
            //          * Convert linked list of matches into array.
            //          */
            //         theseMatches =
            //                 (BlastHSPList**) calloc(numMatches, sizeof(BlastHSPList*));
            //         int i;
            //         for (i = 0; i < numMatches; ++i) {
            //             theseMatches[i] = head->match;
            //             BlastHSPListLinkedList* here = head;
            //             head = head->next;
            //             sfree(here);
            // ```
            let mut matches = Vec::new();
            while let Some(list) = stream.read_observed(self.observer) {
                matches.push(list);
            }
            // EXPERIMENT (LOSAT_X_BXPAR): redo the matches of one batch on the
            // search pool before walking them in stream order.  A match is a
            // function of its own HSP list and of the matrix left by the match
            // before it; `x_speculate_match` returns a result only when the
            // latter was never read, and reports the matrix the match leaves.
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3448-3451,3499-3502
            // ```c
            // #pragma omp for schedule(static)
            //         for (b = 0; b < numMatches; ++b) {
            // #pragma omp flush(interrupt)
            //             if (!interrupt) {
            // ...
            //                 NRrecord             = NRrecord_tld[tid];
            //                 sbp                  = sbp_tld[tid];
            //                 redo_align_params    = redo_align_params_tld[tid];
            //                 matrix               = matrix_tld[tid];
            // ```
            // NCBI redoes the matches in this omp loop, each thread with its own state. The batch
            // of speculative redoes below does the same per-match work on the pool. A result is
            // used only if the match did not read the matrix of its predecessor; every other match
            // is redone in stream order (see x_speculate_match).
            #[cfg(feature = "parallel")]
            let mut speculated: Vec<Option<(HspList, i32, Option<_>)>> = Vec::new();
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334,3499-3502
            // ```c
            //         if ((int) compo_adjust_mode > 1 && !positionBased) {
            //             NRrecord_tld[i] = Blast_CompositionWorkspaceNew();
            //             status_code = Blast_CompositionWorkspaceInit(
            //                     NRrecord_tld[i],
            //                     scoringParams->options->matrix
            //             );
            // ...
            //                 NRrecord             = NRrecord_tld[tid];
            //                 sbp                  = sbp_tld[tid];
            //                 redo_align_params    = redo_align_params_tld[tid];
            //                 matrix               = matrix_tld[tid];
            // ```
            // One Kappa state per pool thread, as NCBI keeps one NRrecord, redo_align_params and
            // matrix per thread.
            #[cfg(feature = "parallel")]
            let x_workers = if x_bx_parallel()
                && self.observer.x_passive()
                && !x_inner_serial()
                && rayon::current_thread_index().is_some()
                && rayon::current_num_threads() > 1
                && matches.len() > 1
            {
                (0..rayon::current_num_threads())
                    .map(|_| {
                        Ok(std::sync::Mutex::new(XKappaWorker(state.x_for_thread(
                            batch,
                            &full,
                            self.options,
                        )?)))
                    })
                    .collect::<Result<Vec<_>>>()?
            } else {
                Vec::new()
            };
            // `x_workers` has one entry per pool thread (none outside a pool; asking
            // Rayon for the thread count there would start its global pool).
            // No NCBI counterpart: size of one speculative batch (scheduling only); it does not
            // change any value NCBI computes.
            #[cfg(feature = "parallel")]
            let x_batch = x_workers.len().saturating_mul(32).max(64);
            // No NCBI counterpart: stage report for the thread-use statistics; it does not change
            // any value NCBI computes.
            #[cfg(feature = "parallel")]
            if x_bx_parallel() {
                crate::utils::threading::report_stage(
                    "blastx",
                    "kappa_redo",
                    matches.len(),
                    !x_workers.is_empty(),
                );
            }
            // No NCBI counterpart: counters of speculative results used and matches redone in
            // order; it does not change any value NCBI computes.
            let mut x_reused = 0usize;
            let mut x_serial = 0usize;
            for (match_index, list) in matches.iter().enumerate() {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3448-3451
                // ```c
                // #pragma omp for schedule(static)
                //         for (b = 0; b < numMatches; ++b) {
                // #pragma omp flush(interrupt)
                //             if (!interrupt) {
                // ```
                // A batch of the following matches is redone ahead of time on the pool. The early-
                // termination test (BlastCompo_EarlyTermination) is evaluated here on the heaps as
                // they are now, and again, in stream order, below; a result is used only when the
                // second test agrees.
                #[cfg(feature = "parallel")]
                if !x_workers.is_empty() && match_index % x_batch == 0 {
                    use rayon::prelude::*;
                    let batch_lists =
                        &matches[match_index..(match_index + x_batch).min(matches.len())];
                    // Early termination only becomes true as the heaps fill.
                    let wanted: Vec<bool> = batch_lists
                        .iter()
                        .map(|list| !CompoHeap::early(list.best_evalue, &heaps))
                        .collect();
                    let (options, subjects, db_length) =
                        (self.options, self.subjects, self.db_length);
                    let full = &full;
                    let workers = &x_workers;
                    speculated = batch_lists
                        .par_iter()
                        .zip(wanted.par_iter())
                        .map(|(list, &wanted)| {
                            if !wanted {
                                return None;
                            }
                            let slot = rayon::current_thread_index().unwrap_or(0) % workers.len();
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
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1582
                // ```c++
                // BlastCompo_EarlyTermination(double evalue,
                //                             BlastCompo_Heap significantMatches[],
                //                             int numQueries)
                // {
                //     int i;
                //     for (i = 0;  i < numQueries;  i++) {
                //         if (BlastCompo_HeapFilledToCutoff(&significantMatches[i])) {
                //             double ecutoff = significantMatches[i].ecutoff;
                //             /* Only matches with evalue <= ethresh will be saved. */
                //             if (evalue <= EVALUE_STRETCH * ecutoff) {
                //                 /* The evalue if this match is sufficiently small
                //                  * that we want to redo it to try to obtain an
                //                  * alignment with evalue smaller than ecutoff. */
                //                 return FALSE;
                //             }
                //         } else {
                //             return FALSE;
                //         }
                //     }
                //     return TRUE;
                // }
                // ```
                let early = CompoHeap::early(list.best_evalue, &heaps);
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
                let x_ready: Option<(HspList, i32, Option<_>)> = None;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3626
                // ```c
                //                                 Blast_RedoOneMatch(
                //                                         alignments,             // thread-local
                // ```
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1232-1234
                // ```c
                //                 if (compo_adjust_mode != eNoCompositionBasedStats &&
                //                         (subject_is_translated || hsp_index == 0
                //                                 || (nearIdenticalStatus != oldNearIdenticalStatus))) {
                // ```
                // Dispatch point: the reference path is the else branch, which redoes the match in
                // stream order (Blast_RedoOneMatch). A speculative result is taken instead only
                // when it was computed without reading its predecessor's matrix; its matrix, if it
                // leaves one, is handed on.
                let (output, best_score) = if let Some((output, best_score, matrix)) = x_ready {
                    x_reused += 1;
                    if matrix.is_some() {
                        state.x_set_matrix(matrix);
                    }
                    (output, best_score)
                } else {
                    x_serial += 1;
                    let encoded: Vec<_> = self.subjects[list.oid]
                        .sequence
                        .iter()
                        .copied()
                        .map(aa_char_to_ncbistdaa)
                        .collect();
                    let input: Vec<_> = list.hsps.iter().map(|h| h.hsp.clone()).collect();
                    let redone =
                        state.redo_list_observed(&input, batch, &encoded, &mut |event| {
                            self.observer.matrix(list.oid, event)
                        })?;
                    finalize_redone(
                        &redone,
                        batch,
                        &full,
                        self.options,
                        &encoded,
                        self.db_length,
                        list.query_index,
                        list.oid,
                        self.observer,
                    )?
                };
                if output.hsps.is_empty() {
                    continue;
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
                // ```c++
                // BlastCompo_HeapWouldInsert(BlastCompo_Heap * self,
                //                            double eValue,
                //                            int score,
                //                            int subject_index)
                // {
                //     if (self->n < self->heapThreshold ||
                //         eValue <= self->ecutoff ||
                //         eValue <  self->worstEvalue) {
                //         return TRUE;
                //     } else {
                //         /* self is either currently a heap, or must be converted to
                //          * one; use s_CompoHeapRecordCompare to compare against
                //          * the worst element in the heap */
                //         BlastCompo_HeapRecord heapRecord; /* temporary record to
                //                                              compare against */
                //         if (self->heapArray == NULL) s_ConvertToHeap(self);
                //
                //         heapRecord.bestEvalue       = eValue;
                //         heapRecord.bestScore        = score;
                //         heapRecord.subject_index    = subject_index;
                //         heapRecord.theseAlignments  = NULL;
                //
                //         return s_CompoHeapRecordCompare(&self->heapArray[1], &heapRecord);
                //     }
                // ```
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
            // ```c++
            // s_FillResultsFromCompoHeaps(BlastHSPResults * results,
            //                             BlastCompo_Heap heaps[],
            //                             Int4 hitlist_size)
            // {
            //     int query_index;   /* loop index */
            //     int num_queries;   /* Number of queries in this search */
            //
            //     num_queries = results->num_queries;
            //     for (query_index = 0;  query_index < num_queries;  query_index++) {
            //         BlastHSPList* hsp_list;
            //         BlastHitList* hitlist;
            //         BlastCompo_Heap * heap = &heaps[query_index];
            //
            //         results->hitlist_array[query_index] = Blast_HitListNew(hitlist_size);
            //         hitlist = results->hitlist_array[query_index];
            //
            //         while (NULL != (hsp_list = BlastCompo_HeapPop(heap))) {
            //             Blast_HitListUpdate(hitlist, hsp_list);
            //         }
            //     }
            //     Blast_HSPResultsReverseOrder(results);
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:444-464
            // ```c++
            // BlastCompo_HeapPop(BlastCompo_Heap * self)
            // {
            //     void * results = NULL;   /* the list of SeqAligns to be returned */
            //
            //     s_ConvertToHeap(self);
            //     if (self->n > 0) { /* The heap is not empty */
            //         BlastCompo_HeapRecord *first, *last; /* The first and last
            //                                                 elements of the array
            //                                                 that represents the
            //                                                 heap.  */
            //         first = &self->heapArray[1];
            //         last  = &self->heapArray[self->n];
            //
            //         results = first->theseAlignments;
            //         if (--self->n > 0) {
            //             /* The heap is still not empty */
            //             memcpy(first, last, sizeof(BlastCompo_HeapRecord));
            //             s_CompoHeapifyDown(self->heapArray, 1, self->n);
            //         }
            //     }
            //     if (COMPO_INTENSE_DEBUG) {
            // ```
            let _ = (x_reused, x_serial);
            for (q, (heap, query)) in heaps.iter_mut().zip(&mut queries).enumerate() {
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
            // ```c++
            //         BlastHitList* hitlist;
            //         BlastCompo_Heap * heap = &heaps[query_index];
            //
            //         results->hitlist_array[query_index] = Blast_HitListNew(hitlist_size);
            //         hitlist = results->hitlist_array[query_index];
            //
            //         while (NULL != (hsp_list = BlastCompo_HeapPop(heap))) {
            //             Blast_HitListUpdate(hitlist, hsp_list);
            //         }
            //     }
            //     Blast_HSPResultsReverseOrder(results);
            // ```
            for query in &mut queries {
                query.lists.reverse();
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3787-3814
            // ```c++
            //             if (redoneMatches_tld[i] != NULL) {
            //                 int qi;
            //                 for (qi = 0; qi < numQueries; ++qi) {
            //                     sfree(redoneMatches_tld[i][qi].array);
            //                     sfree(redoneMatches_tld[i][qi].heapArray);
            //                 }
            //                 s_ClearHeap(redoneMatches_tld[i]);
            //             }
            //         } else {
            //             if (redoneMatches_tld[i] != NULL) {
            //                 int qi;
            //                 for (qi = 0; qi < numQueries; ++qi) {
            //                     sfree(redoneMatches_tld[i][qi].array);
            //                     sfree(redoneMatches_tld[i][qi].heapArray);
            //                 }
            //                 s_ClearHeap(redoneMatches_tld[i]);
            //             }
            //         }
            //         sfree(redoneMatches_tld[i]);
            //     }
            //     if (redoneMatches != NULL) {
            //         int qi;
            //         for (qi = 0; qi < numQueries; ++qi) {
            //             sfree(redoneMatches[qi].array);
            //             sfree(redoneMatches[qi].heapArray);
            //         }
            //         s_ClearHeap(redoneMatches);
            //     }
            // ```
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
        } else {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1640-1707
            // ```c++
            //                                         score_params, hit_params,
            //                                         query_info, pattern_blk);
            //                     } else {
            //                         Boolean fence_hit = FALSE;
            //                         Blast_TracebackFromHSPList(program_number, hsp_list, query,
            //                                          seq_arg.seq, query_info,
            //                                          gap_align, sbp, score_params,
            //                                          ext_params->options, hit_params,
            //                                          seq_arg.seq->gen_code_string,
            //                                          &fence_hit);
            //
            //                         if (fence_hit) {
            //                             /* Disable range support and refetch the
            //                             (whole) subject sequence */
            //
            //                             seq_arg.reset_ranges = TRUE;
            //                             BlastSeqSrcReleaseSequence(seqsrc, &seq_arg);
            //                             BlastSeqSrcGetSequence(seqsrc, &seq_arg);
            //
            //                             /* The C toolkit will erase genetic_code, so do it again */
            //                             if (Blast_SubjectIsTranslated(program_number) &&
            //                                 seq_arg.seq->gen_code_string == NULL) {
            //                             	if (actual_num_threads > 1) {
            // #pragma omp critical(tback_gen_code)
            //                                     seq_arg.seq->gen_code_string =
            //                                         GenCodeSingletonFind(default_db_genetic_code);
            // #ifndef _OPENMP
            //                                     ASSERT(seq_arg.seq->gen_code_string);
            // #endif
            //                                 }
            //                             	else {
            //                                     seq_arg.seq->gen_code_string =
            //                                         GenCodeSingletonFind(default_db_genetic_code);
            //                             	}
            //                             }
            //
            //                             /* Retry the alignment with fence_hit set*/
            //                             Blast_TracebackFromHSPList(program_number, hsp_list,
            //                                                 query, seq_arg.seq,
            //                                                 query_info, gap_align,
            //                                                 sbp, score_params,
            //                                                 ext_params->options,
            //                                                 hit_params,
            //                                                 seq_arg.seq->gen_code_string,
            //                                                 &fence_hit);
            // #ifndef _OPENMP
            //                             ASSERT(fence_hit == FALSE);
            // #endif
            //                         } /* fence_hit */
            //                     }    /* !phi_blast */
            //
            //                 } else {
            //                     /* traceback skipped; compute bit scores for searches
            //                        where the traceback phase is seperated from the
            //                        preliminary search. */
            //                     Blast_HSPListGetBitScores(hsp_list, FALSE, sbp);
            //                 }
            //
            //                 /* Free HSP list if all HSPs have been deleted. */
            //
            //                 batch->hsplist_array[hsplist_itr] = NULL;
            //                 if (hsp_list->hspcnt == 0) {
            //                     hsp_list = Blast_HSPListFree(hsp_list);
            //                 }
            //                 else {
            //                     Blast_HSPResultsInsertHSPList(thread_data->tld[tid]->results, hsp_list,
            //                                   hit_params->options->hitlist_size);
            //                 }
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1500-1508
            // ```c++
            //         Int4 i;
            //         Uint4 actual_num_threads = 0;
            //         BlastHSPStreamResultsBatchArray* batches = NULL;
            //         Boolean has_been_interrupted = FALSE;
            //
            //         if ( (retval = BlastHSPStreamToHSPStreamResultsBatch(hsp_stream, &batches))) {
            //             return retval;
            //         }
            //         ASSERT(batches);
            // ```
            for lists in result_batches {
                let encoded: Vec<_> = self.subjects[lists[0].oid]
                    .sequence
                    .iter()
                    .copied()
                    .map(aa_char_to_ncbistdaa)
                    .collect();
                for list in lists {
                    let input: Vec<_> = list.hsps.iter().map(|h| h.hsp.clone()).collect();
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:402-407
                    // ```c++
                    //        */
                    //       if (program_number == eBlastTypeRpsBlast ||
                    //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                    //                                hit_options->min_diag_separation)) {
                    //
                    //          Int4 start_shift = 0;
                    // ```
                    let redone = super::traceback::ordinary_traceback_observed(
                        &input,
                        batch,
                        &full,
                        self.options,
                        &encoded,
                        &mut |h, contained| self.observer.containment(h, contained),
                    )?;
                    let (output, _) = finalize_redone(
                        &redone,
                        batch,
                        &full,
                        self.options,
                        &encoded,
                        self.db_length,
                        list.query_index,
                        list.oid,
                        self.observer,
                    )?;
                    if !output.hsps.is_empty() {
                        queries[list.query_index].update_observed(output, self.observer);
                    }
                }
            }
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:849-865
        // ```c++
        //     	  	  if(hit_options->max_hsps_per_subject) {
        //     	  		  Blast_TrimHSPListByMaxHsps(hsp_list, hit_options);
        //     	  	  }
        //     	  	  if(hit_options->query_cov_hsp_perc) {
        //     	  		  Blast_HSPListReapByQueryCoverage(hsp_list, hit_options, query_info, program_number);
        //     	  		  if(hsp_list->hspcnt == 0){
        //     	  			hit_list->hsplist_array[subject_index] =  Blast_HSPListFree(hsp_list);
        //     	  		  }
        //     	  	  }
        //     	  	  if((hit_options->hsp_filt_opt != NULL) && (hit_options->hsp_filt_opt->subject_besthit_opts != NULL)) {
        //     	  		  Blast_HSPListSubjectBestHit(program_number,
        //     	  		  				           hit_options->hsp_filt_opt->subject_besthit_opts,
        //     	  		  				           query_info, hsp_list);
        //     	  	  }
        //        }
        //        if(hit_options->query_cov_hsp_perc) {
        //             Blast_HitListPurgeNullHSPLists(hit_list);
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1717-1723
        // ```c++
        //         results = SThreadLocalDataArrayConsolidateResults(thread_data);
        //         ASSERT(results);
        //
        //         /* post-traceback pipes */
        //         BlastHSPStreamTBackClose(hsp_stream, results);
        //
        //     } /* end of else */
        // ```
        if self.options.culling_limit > 0 {
            super::culling::traceback_pipe(
                &mut queries,
                batch.contexts.iter().map(|c| c.length as i32).collect(),
                self.options.culling_limit,
                prelim_hitlist_size(self.options),
                self.observer,
            );
        }
        for query in &mut queries {
            for list in &mut query.lists {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2049-2070
                // ```c++
                // Int2 Blast_TrimHSPListByMaxHsps(BlastHSPList* hsp_list,
                //                                 const BlastHitSavingOptions* hit_options)
                // {
                //    BlastHSP** hsp_array;
                //    Int4 index;
                //    Int4 hsp_max;
                //
                //    if ((hsp_list == NULL) ||
                // 	   (hit_options->max_hsps_per_subject == 0) ||
                // 	   (hsp_list->hspcnt <= hit_options->max_hsps_per_subject))
                //       return 0;
                //
                //    hsp_max = hit_options->max_hsps_per_subject;
                //    hsp_array = hsp_list->hsp_array;
                //    for (index = hsp_max; index < hsp_list->hspcnt; index++) {
                //       hsp_array[index] = Blast_HSPFree(hsp_array[index]);
                //    }
                //
                //    hsp_list->hspcnt = hsp_max;
                //    return 0;
                // }
                //
                // ```
                if self.options.max_hsps > 0 {
                    self.observer.hsps("MAX_HSPS_IN", list.oid, &list.hsps);
                    list.hsps.truncate(self.options.max_hsps as usize);
                    self.observer.hsps("MAX_HSPS_OUT", list.oid, &list.hsps);
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2538-2588
                // ```c++
                // Blast_HSPListSubjectBestHit(EBlastProgramType program,
                // 		                    const BlastHSPSubjectBestHitOptions* subject_besthit_opts,
                // 		                    const BlastQueryInfo *query_info,
                //                             BlastHSPList* hsp_list)
                // {
                //    BlastHSP** hsp_array;  /* hsp_array to purge. */
                //    const int range_diff = subject_besthit_opts->max_range_diff;
                //    Boolean isBlastn = (program == eBlastTypeBlastn);
                //    unsigned int i, j;
                //    int o, e;
                //    int curr_context, target_context;
                //    int qlen;
                //
                //    /* If HSP list is empty, return immediately. */
                //    if (hsp_list == NULL || hsp_list->hspcnt == 0)
                //        return 0;
                //
                //    if (Blast_ProgramIsPhiBlast(program))
                //        return hsp_list->hspcnt;
                //
                //    hsp_array = hsp_list->hsp_array;
                //
                //    // The hsp list is sorted by score
                //    for(i=0; i < hsp_list->hspcnt -1; i++) {
                // 	  if(hsp_array[i] == NULL){
                // 		  continue;
                // 	  }
                //       j = 1;
                //       o = hsp_array[i]->query.offset - range_diff;
                //       e = hsp_array[i]->query.end + range_diff;
                //       if (o < 0) o = 0;
                //       if (e < 0) e = hsp_array[i]->query.end;
                //       while (i+j < hsp_list->hspcnt) {
                //           if (hsp_array[i+j] && hsp_array[i]->context == hsp_array[i+j]->context &&
                //               ((hsp_array[i+j]->query.offset >= o) &&
                //                (hsp_array[i+j]->query.end <= e))){
                //        	      hsp_array[i+j] = Blast_HSPFree(hsp_array[i+j]);
                //           }
                //           j++;
                //       }
                //    }
                //
                //    Blast_HSPListPurgeNullHSPs(hsp_list);
                //
                //    if(isBlastn) {
                // 	   for(i=0; i < hsp_list->hspcnt -1; i++) {
                // 	   	  if(hsp_array[i] == NULL){
                // 	   		  continue;
                // 	   	  }
                // 	   	  // Flip query offsets of current hsp to target context frame
                // 	   	  j = 1;
                // ```
                if self.options.subject_besthit && !list.hsps.is_empty() {
                    self.observer.hsps("SUBJECT_BEST_IN", list.oid, &list.hsps);
                    subject_besthit(&mut list.hsps);
                    self.observer.hsps("SUBJECT_BEST_OUT", list.oid, &list.hsps);
                }
            }
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1768-1777
        // ```c++
        //     /* Re-sort the hit lists according to their best e-values, because they
        //        could have changed. Only do this for a database search. */
        //     if (BlastSeqSrcGetTotLen(seq_src) > 0) {
        //         Blast_HSPResultsSortByEvalue(results);
        //     }
        //
        //
        //     /* Eliminate extra hits from results, if preliminary hit list size is
        //        larger than the final hit list size */
        //     s_BlastPruneExtraHits(results, hit_params->options->hitlist_size);
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1768-1777
        // ```c++
        //     /* Re-sort the hit lists according to their best e-values, because they
        //        could have changed. Only do this for a database search. */
        //     if (BlastSeqSrcGetTotLen(seq_src) > 0) {
        //         Blast_HSPResultsSortByEvalue(results);
        //     }
        //
        //
        //     /* Eliminate extra hits from results, if preliminary hit list size is
        //        larger than the final hit list size */
        //     s_BlastPruneExtraHits(results, hit_params->options->hitlist_size);
        // ```
        for query in &mut queries {
            query.lists.sort_by(super::results::compare_list);
        }
        for (q, query) in queries.iter().enumerate() {
            self.observer.targets("IN", q, &query.lists);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:878-894
        // ```c++
        // s_BlastPruneExtraHits(BlastHSPResults* results, Int4 hitlist_size)
        // {
        //    Int4 query_index, subject_index;
        //    BlastHitList* hit_list;
        //
        //    for (query_index = 0; query_index < results->num_queries; ++query_index) {
        //       if (!(hit_list = results->hitlist_array[query_index]))
        //          continue;
        //       for (subject_index = hitlist_size;
        //            subject_index < hit_list->hsplist_count; ++subject_index) {
        //          hit_list->hsplist_array[subject_index] =
        //          Blast_HSPListFree(hit_list->hsplist_array[subject_index]);
        //       }
        //       hit_list->hsplist_count = MIN(hit_list->hsplist_count, hitlist_size);
        //    }
        // }
        //
        // ```
        for query in &mut queries {
            query.lists.truncate(self.options.hitlist_size as usize);
        }
        for (q, query) in queries.iter().enumerate() {
            self.observer.targets("OUT", q, &query.lists);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-101
        // ```c++
        //     // find the first valid context corresponding to this query
        //     for (i = 0; i < context_per_query; i++) {
        //         BlastContextInfo *ctx = query_info->contexts +
        //                                 query_number * context_per_query + i;
        //         if (ctx->is_valid) {
        //             m_SearchSpace = ctx->eff_searchsp;
        // 	    m_LengthAdjustment = ctx->length_adjustment;
        //             break;
        //         }
        //     }
        //     if (i >= context_per_query) {
        //         return; // we didn't find a valid context :(
        //     }
        //
        //     // fill in the Karlin blocks for that context, if they
        //     // are valid
        //     const int ctx_index = query_number * context_per_query + i;
        //     if (sbp->kbp_std) {
        //         s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
        //     }
        // ```
        self.output.push(BatchResults {
            query_ordinal: self.ordinal,
            queries: queries.into_iter().map(|q| q.lists).collect(),
            prepared: batch.clone(),
            parameters: full,
            search_skipped: false,
        });
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3670-3703
// ```c++
//                 if (hsp_list->hspcnt > 1) {
//                     s_HitlistReapContained(hsp_list->hsp_array,
//                             &hsp_list->hspcnt);
//                 }
//                 *pStatusCode =
//                         s_HitlistEvaluateAndPurge(&best_score, &best_evalue,
//                                 hsp_list,
//                                 seqSrc,
//                                 matchingSeq.length,
//                                 program_number,
//                                 queryInfo, context_index,
//                                 sbp, hitParams,
//                                 pvalueForThisPair, LambdaRatio,
//                                 matchingSeq.index);
//                 if (*pStatusCode != 0) {
//                     goto query_loop_cleanup;
//                 }
//                 if (best_evalue <= hitParams->options->expect_value) {
//                     /* The best alignment is significant */
//                     s_HSPListNormalizeScores(hsp_list, kbp->Lambda, kbp->logK,
//                             localScalingFactor);
//                     s_ComputeNumIdentities(
//                             queryBlk,
//                             queryInfo,
//                             subjectBlk,
//                             seqSrc,
//                             hsp_list,
//                             scoringParams->options,
//                             genetic_code_string,
//                             sbp,
//                             ranges
//                     );
//                     if (!seqSrc) {
//                         goto query_loop_cleanup;
// ```
#[allow(clippy::too_many_arguments)]
fn finalize_redone(
    redone: &[super::kappa::RedoneHsp],
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    db_length: i64,
    query_index: usize,
    oid: usize,
    observer: &mut dyn RuntimeObserver,
) -> Result<(HspList, i32)> {
    let input: Vec<_> = redone.iter().map(|h| h.hsp.clone()).collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
    // ```c++
    // {
    //     Int4 index;
    //
    //     if (!hsp_list || hsp_list->hspcnt == 0)
    //         return 0;
    //
    //     ASSERT(link_hsp_params);
    //
    //     /* Remove any information on number of linked HSPs from previous
    //        linking. */
    //     for (index = 0; index < hsp_list->hspcnt; ++index)
    //         hsp_list->hsp_array[index]->num = 1;
    //
    //     /* Link up the HSP's for this hsp_list. */
    //     if (link_hsp_params->longest_intron <= 0) {
    //         s_BlastEvenGapLinkHSPs(program_number, hsp_list, query_info,
    //                               subject_length, sbp, link_hsp_params,
    //                               gapped_calculation);
    //         /* The HSP's may be in a different order than they were before,
    //            but hsp contains the first one. */
    //     } else {
    //         Blast_HSPListAdjustOddBlastnScores(hsp_list, gapped_calculation, sbp);
    //         /* Calculate individual HSP e-values first - they'll be needed to
    //            compare with sum e-values. Use decay rate to compensate for
    //            multiple tests. */
    // ```
    if !input.is_empty() {
        if let Some(link) = super::linking::link_parameters(options) {
            observer.link_run(subject.len() as i32, true, &link);
        }
    }
    let mut linked = evaluate_gapped(
        &input,
        batch,
        params,
        options,
        subject.len() as i32,
        db_length,
        if options.composition == 2 { 32.0 } else { 1.0 },
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1976-2008
    // ```c++
    // Int2 Blast_HSPListReapByEvalue(BlastHSPList* hsp_list,
    //         const BlastHitSavingOptions* hit_options)
    // {
    //    BlastHSP* hsp;
    //    BlastHSP** hsp_array;
    //    Int4 hsp_cnt = 0;
    //    Int4 index;
    //    double cutoff;
    //
    //    if (hsp_list == NULL)
    //       return 0;
    //
    //    cutoff = hit_options->expect_value;
    //
    //    hsp_array = hsp_list->hsp_array;
    //    for (index = 0; index < hsp_list->hspcnt; index++) {
    //       hsp = hsp_array[index];
    //
    //       ASSERT(hsp != NULL);
    //
    //       if (hsp->evalue > cutoff) {
    //          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
    //       } else {
    //          if (index > hsp_cnt)
    //             hsp_array[hsp_cnt] = hsp_array[index];
    //          hsp_cnt++;
    //       }
    //    }
    //
    //    hsp_list->hspcnt = hsp_cnt;
    //
    //    return 0;
    // }
    // ```
    let observe_reap = !linked.hsps.is_empty();
    if observe_reap {
        observer.reap("FINAL_REAP_IN", oid, &linked, options.evalue);
    }
    reap_final(&mut linked, options);
    if observe_reap {
        observer.reap("FINAL_REAP_OUT", oid, &linked, options.evalue);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:438-444
    // ```c++
    //     if (status == 0) {
    //         Blast_HSPListReapByEvalue(hsp_list, hitParams->options);
    //         if (hsp_list->hspcnt > 0) {
    //             *pbestEvalue = hsp_list->best_evalue;
    //             *pbestScore  = hsp_list->hsp_array[0]->score;
    //         }
    //     }
    // ```
    let best_score = linked.hsps.first().map(|h| h.hsp.score).unwrap_or(0);
    let bits = normalize_scores(&mut linked, options);
    let hsps = linked
        .hsps
        .iter()
        .zip(bits)
        .map(|(h, bit_score)| {
            let r = &redone[h.source_index];
            let (identity, positive) = identity_and_positive(r, batch, subject);
            Hsp {
                hsp: h.hsp.clone(),
                num: h.num,
                evalue: h.evalue,
                bit_score,
                identity,
                positive,
                edit_script: r.edit_script.clone(),
                composition_method: r.composition_method,
            }
        })
        .collect();
    let list = HspList {
        query_index,
        oid,
        hsps,
        best_evalue: linked.best_evalue,
    };
    Ok((list, best_score))
}
