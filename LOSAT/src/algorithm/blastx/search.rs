//! Internal native preliminary boundary. Linking/Kappa/traceback/report remain separate.
use super::{
    args::ResolvedOptions,
    input::FastaRecord,
    parameters::{effective_lengths, score_block, subject_parameters, ContextParameters},
    preliminary::{gapped, purge_endpoints, sort_by_score, ungapped, PreliminaryHsp},
    query_setup::{batch_ranges, prepare_queries, PreparedQueryBatch},
    seed::{word_finder, Diagonals, InitHsp, Lookup},
    split::{adjust_chunk_hsps, merge_chunk_hsps, split_queries, QueryChunk},
};
use crate::utils::matrix::aa_char_to_ncbistdaa;
use anyhow::{ensure, Result};

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:235-271
// ```c++
//         CRef<CSplitQueryBlk> split_query_blk = query_splitter->Split();
//
//         for (Uint4 i = 0; i < query_splitter->GetNumberOfChunks(); i++) {
//             try {
//                 CRef<IQueryFactory> chunk_qf =
//                     query_splitter->GetQueryFactoryForChunk(i);
//                 _TRACE("Query chunk " << i << "/" <<
//                        query_splitter->GetNumberOfChunks());
//                 CRef<SInternalData> chunk_data =
//                     SplitQuery_CreateChunkData(chunk_qf, m_Options,
//                                                m_InternalData,
//                                                GetNumberOfThreads());
//
//                 CRef<ILocalQueryData> query_data(
//                         chunk_qf->MakeLocalQueryData( &*m_Options ) );
//                 BLAST_SequenceBlk * chunk_queries =
//                     query_data->GetSequenceBlk();
//                 GetDbIndexSetUsingThreadsFn()( IsMultiThreaded() );
//                 GetDbIndexRunSearchFn()(
//                         chunk_queries, lut_options, word_options );
//
//                 if (IsMultiThreaded()) {
//                      x_LaunchMultiThreadedSearch(*chunk_data);
//                 } else {
//                     retval =
//                         CPrelimSearchRunner(*chunk_data, opts_memento.get())();
//                     if (retval) {
//                         NCBI_THROW(CBlastException, eCoreBlastError,
//                                    BlastErrorCode2String(retval));
//                     }
//                 }
//
//
//                 _ASSERT(chunk_data->m_HspStream->GetPointer());
//                 BlastHSPStreamMerge(split_query_blk->GetCStruct(), i,
//                                 chunk_data->m_HspStream->GetPointer(),
//                                 m_InternalData->m_HspStream->GetPointer());
// ```
pub(crate) trait PreliminarySink {
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
        _ordinal: usize,
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
    ) -> Result<()> {
        Ok(())
    }
    fn batch_start(
        &mut self,
        _ordinal: usize,
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
    ) -> Result<()> {
        Ok(())
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3916-3921
    // ```c++
    //
    //       /* use priate interval tree when recomputing alignments */
    //       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
    //                                         hit_options->min_diag_separation))
    //       {
    //          BlastHSP* new_hsp;
    // ```
    fn containment(&mut self, _h: &PreliminaryHsp, _contained: bool) {}
    fn chunk_start(&mut self, _chunk: &QueryChunk) -> Result<()> {
        Ok(())
    }
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
    fn preliminary_purge(
        &mut self,
        _oid: usize,
        _before: &[PreliminaryHsp],
        _after: &[PreliminaryHsp],
    ) -> Result<()> {
        Ok(())
    }
    fn subject_start(
        &mut self,
        _oid: usize,
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
        _prepared_link: Option<([i32; 2], f64)>,
    ) -> Result<()> {
        Ok(())
    }
    fn subject(
        &mut self,
        _oid: usize,
        _hsps: &[PreliminaryHsp],
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
    ) -> Result<()> {
        Ok(())
    }
    fn chunk_end(&mut self, _chunk: &QueryChunk, _is_split: bool) -> Result<()> {
        Ok(())
    }
    fn batch_end(
        &mut self,
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
    ) -> Result<()> {
        Ok(())
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491-505
// ```c++
//         BlastInitHitListReset(init_hitlist);
//
//         if (aux_struct->WordFinder) {
//             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
//                                    word_params, aux_struct->ewp,
//                                    aux_struct->offset_pairs,
//                                    kScanSubjectOffsetArraySize,
//                                    init_hitlist, ungapped_stats);
//
//             if (init_hitlist->total == 0) continue;
//         }
//
//         if (score_options->gapped_calculation) {
//             Int4 prot_length = 0;
//             if (score_options->is_ooframe) {
// ```
struct DiagnosticOnly;
impl PreliminarySink for DiagnosticOnly {}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491-526
// ```c++
//         BlastInitHitListReset(init_hitlist);
//
//         if (aux_struct->WordFinder) {
//             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
//                                    word_params, aux_struct->ewp,
//                                    aux_struct->offset_pairs,
//                                    kScanSubjectOffsetArraySize,
//                                    init_hitlist, ungapped_stats);
//
//             if (init_hitlist->total == 0) continue;
//         }
//
//         if (score_options->gapped_calculation) {
//             Int4 prot_length = 0;
//             if (score_options->is_ooframe) {
//                 /* Convert query offsets in all HSPs into the mixed-frame
//                    coordinates */
//                 s_TranslateHSPsToDNAPCoord(program_number, init_hitlist,
//                        query_info, subject->frame, orig_length, backup.offset);
//                 if (kTranslatedSubject) {
//                     prot_length = subject->length;
//                     subject->length = orig_length;
//                 }
//             }
//         /** NB: If queries are concatenated, HSP offsets must be adjusted
//           * inside the following function call, so coordinates are
//           * relative to the individual contexts (i.e. queries, strands or
//           * frames). Contexts should also be filled in HSPs when they
//           * are saved.
//           */
//         /* fence_hit is null, since this is only for prelim stage. */
//         if (aux_struct->GetGappedScore) {
//             status = aux_struct->GetGappedScore(program_number, query,
//                     query_info,
//                     subject, gap_align, score_params, ext_params, hit_params,
//                     word_params, init_hitlist, &hsp_list, gapped_stats, NULL);
// ```
pub struct SubjectStage {
    pub oid: usize,
    pub call: usize,
    pub parameters: Vec<ContextParameters>,
    pub gap_xdrop: i32,
    pub seeds: Vec<(i32, i32)>,
    pub initial: Vec<InitHsp>,
    pub raw: Vec<PreliminaryHsp>,
    pub purged: Vec<PreliminaryHsp>,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:235-271
// ```c++
//         CRef<CSplitQueryBlk> split_query_blk = query_splitter->Split();
//
//         for (Uint4 i = 0; i < query_splitter->GetNumberOfChunks(); i++) {
//             try {
//                 CRef<IQueryFactory> chunk_qf =
//                     query_splitter->GetQueryFactoryForChunk(i);
//                 _TRACE("Query chunk " << i << "/" <<
//                        query_splitter->GetNumberOfChunks());
//                 CRef<SInternalData> chunk_data =
//                     SplitQuery_CreateChunkData(chunk_qf, m_Options,
//                                                m_InternalData,
//                                                GetNumberOfThreads());
//
//                 CRef<ILocalQueryData> query_data(
//                         chunk_qf->MakeLocalQueryData( &*m_Options ) );
//                 BLAST_SequenceBlk * chunk_queries =
//                     query_data->GetSequenceBlk();
//                 GetDbIndexSetUsingThreadsFn()( IsMultiThreaded() );
//                 GetDbIndexRunSearchFn()(
//                         chunk_queries, lut_options, word_options );
//
//                 if (IsMultiThreaded()) {
//                      x_LaunchMultiThreadedSearch(*chunk_data);
//                 } else {
//                     retval =
//                         CPrelimSearchRunner(*chunk_data, opts_memento.get())();
//                     if (retval) {
//                         NCBI_THROW(CBlastException, eCoreBlastError,
//                                    BlastErrorCode2String(retval));
//                     }
//                 }
//
//
//                 _ASSERT(chunk_data->m_HspStream->GetPointer());
//                 BlastHSPStreamMerge(split_query_blk->GetCStruct(), i,
//                                 chunk_data->m_HspStream->GetPointer(),
//                                 m_InternalData->m_HspStream->GetPointer());
// ```
pub struct ChunkStage {
    pub statistically_valid: bool,
    pub chunk: QueryChunk,
    pub subjects: Vec<SubjectStage>,
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
pub struct BatchStage {
    pub query_ordinal: usize,
    pub restored: PreparedQueryBatch,
    pub full_parameters: Vec<ContextParameters>,
    pub chunks: Vec<ChunkStage>,
    pub merged_by_query: Vec<Vec<Vec<PreliminaryHsp>>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-246
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
//
//         for (Uint4 i = 0; i < query_splitter->GetNumberOfChunks(); i++) {
//             try {
//                 CRef<IQueryFactory> chunk_qf =
//                     query_splitter->GetQueryFactoryForChunk(i);
//                 _TRACE("Query chunk " << i << "/" <<
//                        query_splitter->GetNumberOfChunks());
//                 CRef<SInternalData> chunk_data =
//                     SplitQuery_CreateChunkData(chunk_qf, m_Options,
//                                                m_InternalData,
//                                                GetNumberOfThreads());
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1428-1443
// ```c++
//
//        if (seq_arg.seq->length < min_subj_seq_length) {
//            BlastSeqSrcReleaseSequence(seq_src, &seq_arg);
//            continue;
//        }
//
//        if (db_length == 0) {
//            /* This is not a database search, hence need to recalculate and save
//             the effective search spaces and length adjustments for all
//             queries based on the length of the current single subject
//             sequence. */
//            if ((status = BLAST_OneSubjectUpdateParameters(program_number,
//                           seq_arg.seq->length, score_options, query_info,
//                           sbp, hit_params, word_params,
//                           eff_len_params)) != 0)
//               return status;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491-555
// ```c++
//         BlastInitHitListReset(init_hitlist);
//
//         if (aux_struct->WordFinder) {
//             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
//                                    word_params, aux_struct->ewp,
//                                    aux_struct->offset_pairs,
//                                    kScanSubjectOffsetArraySize,
//                                    init_hitlist, ungapped_stats);
//
//             if (init_hitlist->total == 0) continue;
//         }
//
//         if (score_options->gapped_calculation) {
//             Int4 prot_length = 0;
//             if (score_options->is_ooframe) {
//                 /* Convert query offsets in all HSPs into the mixed-frame
//                    coordinates */
//                 s_TranslateHSPsToDNAPCoord(program_number, init_hitlist,
//                        query_info, subject->frame, orig_length, backup.offset);
//                 if (kTranslatedSubject) {
//                     prot_length = subject->length;
//                     subject->length = orig_length;
//                 }
//             }
//         /** NB: If queries are concatenated, HSP offsets must be adjusted
//           * inside the following function call, so coordinates are
//           * relative to the individual contexts (i.e. queries, strands or
//           * frames). Contexts should also be filled in HSPs when they
//           * are saved.
//           */
//         /* fence_hit is null, since this is only for prelim stage. */
//         if (aux_struct->GetGappedScore) {
//             status = aux_struct->GetGappedScore(program_number, query,
//                     query_info,
//                     subject, gap_align, score_params, ext_params, hit_params,
//                     word_params, init_hitlist, &hsp_list, gapped_stats, NULL);
//         }
//         else if (aux_struct->JumperGapped) {
//             status = aux_struct->JumperGapped(subject, query, query_info,
//                                               lookup, word_params,
//                                               score_params, hit_params,
//                                               aux_struct->offset_pairs,
//                                               aux_struct->mapper_wordhits,
//                                               kScanSubjectOffsetArraySize,
//                                               gap_align, init_hitlist,
//                                               &hsp_list, ungapped_stats,
//                                               gapped_stats);
//         }
//         if (status) break;
//
//         /* No need to do this for short reads */
//         if (aux_struct->GetGappedScore) {
//
//             /* Removes redundant HSPs. */
//             Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
//
//             /* For nucleotide search, if match score is = 2, the odd scores
//                are rounded down to the nearest even number. */
// #if 0
//             Blast_HSPListAdjustOddBlastnScores(hsp_list, score_options->gapped_calculation, gap_align->sbp);
// #endif
//
//         }
//
//         Blast_HSPListSortByScore(hsp_list);
// ```
pub fn search_preliminary(
    records: &[FastaRecord],
    subjects: &[FastaRecord],
    options: &ResolvedOptions,
) -> Result<Vec<BatchStage>> {
    crate::utils::threading::with_search_pool(options.num_threads as usize, "BLASTX", |pool| {
        search_core(records, subjects, options, true, &mut DiagnosticOnly, pool)
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491-526
// ```c++
//         BlastInitHitListReset(init_hitlist);
//
//         if (aux_struct->WordFinder) {
//             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
//                                    word_params, aux_struct->ewp,
//                                    aux_struct->offset_pairs,
//                                    kScanSubjectOffsetArraySize,
//                                    init_hitlist, ungapped_stats);
//
//             if (init_hitlist->total == 0) continue;
//         }
//
//         if (score_options->gapped_calculation) {
//             Int4 prot_length = 0;
//             if (score_options->is_ooframe) {
//                 /* Convert query offsets in all HSPs into the mixed-frame
//                    coordinates */
//                 s_TranslateHSPsToDNAPCoord(program_number, init_hitlist,
//                        query_info, subject->frame, orig_length, backup.offset);
//                 if (kTranslatedSubject) {
//                     prot_length = subject->length;
//                     subject->length = orig_length;
//                 }
//             }
//         /** NB: If queries are concatenated, HSP offsets must be adjusted
//           * inside the following function call, so coordinates are
//           * relative to the individual contexts (i.e. queries, strands or
//           * frames). Contexts should also be filled in HSPs when they
//           * are saved.
//           */
//         /* fence_hit is null, since this is only for prelim stage. */
//         if (aux_struct->GetGappedScore) {
//             status = aux_struct->GetGappedScore(program_number, query,
//                     query_info,
//                     subject, gap_align, score_params, ext_params, hit_params,
//                     word_params, init_hitlist, &hsp_list, gapped_stats, NULL);
// ```
pub(crate) fn search_core(
    records: &[FastaRecord],
    subjects: &[FastaRecord],
    options: &ResolvedOptions,
    diagnostics: bool,
    sink: &mut impl PreliminarySink,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Result<Vec<BatchStage>> {
    ensure!(
        !subjects.is_empty(),
        "BLASTX preliminary requires a protein subject"
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181
    // ```c++
    //     if(dbscan_mode)
    //     {
    // 	ITERATE(vector<BLAST_SequenceBlk*>, iter, m_ivSeqBlkVec)
    // 	{
    //             m_iTotalLength += (Int8) (*iter)->length;
    // 	}
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:53-54
    // ```c++
    // /** Split subject sequences if longer than this */
    // #define MAX_DBSEQ_LEN 5000000
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-250
    // ```c++
    //     if (backup->offset + MAX_DBSEQ_LEN <
    //         backup->hard_ranges[backup->hm_index].right) {
    //
    //         subject->length = MAX_DBSEQ_LEN;
    //         backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
    // ```
    ensure!(
        subjects.iter().all(|s| s.sequence.len() <= 5_000_000),
        "BLASTX preliminary protein-subject splitting above 5000000 residues is not implemented"
    );
    let db_length = subjects
        .iter()
        .try_fold(0i64, |sum, s| {
            sum.checked_add(i64::try_from(s.sequence.len()).ok()?)
        })
        .ok_or_else(|| anyhow::anyhow!("BLASTX subject-set length exceeds Int8"))?;
    let db_count = i64::try_from(subjects.len())?;
    let min_subject_length = subjects
        .iter()
        .map(|s| s.sequence.len())
        .min()
        .expect("protein subject");
    let mut result = Vec::new();
    let mut call = 0;
    for range in batch_ranges(records) {
        let queries = &records[range.clone()];
        let mut full = prepare_queries(queries, options)?;
        let mut full_params = score_block(&full);
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
        if !full_params.iter().any(|p| p.valid) {
            for context in &mut full.contexts {
                context.is_valid = false;
            }
            sink.unsearched_batch(range.start, &full, &full_params)?;
            if diagnostics {
                result.push(BatchStage {
                    query_ordinal: range.start,
                    restored: full,
                    full_parameters: full_params,
                    chunks: Vec::new(),
                    merged_by_query: vec![vec![Vec::new(); subjects.len()]; queries.len()],
                });
            }
            continue;
        }
        effective_lengths(&full, &mut full_params, options, db_length, db_count, None)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2781-2787
        // ```c++
        //       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
        //       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
        //       if (loop_status) {
        //           contexts[context].is_valid = FALSE;
        //           sbp->sfp[context] = Blast_ScoreFreqFree(sbp->sfp[context]);
        //           sbp->kbp_std[context] = Blast_KarlinBlkFree(sbp->kbp_std[context]);
        //           if (!Blast_QueryIsTranslated(program) ) {
        // ```
        for (context, p) in full.contexts.iter_mut().zip(&full_params) {
            context.is_valid = p.valid;
        }
        let spaces: Vec<_> = full_params.iter().map(|p| p.search_space).collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-222
        // ```c++
        //     CEffectiveSearchSpacesMemento eff_memento(m_Options);
        //     SplitQuery_SetEffectiveSearchSpace(m_Options, m_QueryFactory,
        //                                        m_InternalData);
        //     int retval = 0;
        // ```
        sink.batch_start(range.start, &full, &full_params)?;
        let chunks = split_queries(queries, &full, options)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:233-236
        // ```c++
        //     if (query_splitter->IsQuerySplit()) {
        //
        //         CRef<CSplitQueryBlk> split_query_blk = query_splitter->Split();
        //
        // ```
        let is_split = chunks.len() > 1;
        let mut stages = Vec::new();
        let mut merged = vec![vec![Vec::new(); subjects.len()]; queries.len()];
        for chunk in chunks {
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
            sink.chunk_start(&chunk)?;
            let batch = &chunk.prepared;
            let mut params = score_block(batch);
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:277-288
            // ```c++
            //             } catch (const CBlastException& e) {
            //                 // This error message is safe to ignore for a given chunk,
            //                 // because the chunks might end up producing a region of
            //                 // the query for which ungapped Karlin-Altschul blocks
            //                 // cannot be calculated
            //                 const string err_msg1("search cannot proceed due to errors "
            //                                      "in all contexts/frames of query "
            //                                      "sequences");
            //                 const string err_msg2(kBlastErrMsg_CantCalculateUngappedKAParams);
            //                 if (e.GetMsg().find(err_msg1) == NPOS && e.GetMsg().find(err_msg2) == NPOS) {
            //                     throw;
            //                 }
            // ```
            if !params.iter().any(|p| p.valid) {
                if diagnostics {
                    stages.push(ChunkStage {
                        chunk,
                        statistically_valid: false,
                        subjects: Vec::new(),
                    });
                }
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:685-696
            // ```c++
            //         if (context_index != 0) {
            //             Blast_MessageWrite(blast_message, eBlastSevWarning, context_index,
            //                     "One search space is being used for multiple sequences");
            //         }
            //         retval = eff_len_options->searchsp_eff[0];
            //     } else if (eff_len_options->num_searchspaces > 1) {
            //         ASSERT(context_index < eff_len_options->num_searchspaces);
            //         retval = eff_len_options->searchsp_eff[context_index];
            //     } else {
            //         abort();    /* should never happen */
            //     }
            //     return retval;
            // ```
            let fixed_spaces = &spaces;
            let lookup = Lookup::new(batch, &params, options)?;
            let last = batch.contexts.last().expect("BLASTX contexts");
            let mut diagonals = Diagonals::new(last.offset + last.length, options.window_size)?;
            let xdrop = subject_parameters(
                batch,
                &mut params,
                options,
                db_length,
                db_count,
                min_subject_length,
                Some(fixed_spaces),
            )?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:172-180
            // ```c++
            //     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
            //         (*thread)->Run();
            //     }
            //
            //     // ... and wait for the threads to finish
            //     Uint8 retv(0);
            //     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
            //         void * result(0);
            //         (*thread)->Join(&result);
            // ```
            // utils::threading owns the search pool; each slot retains its own
            // diagonal table across bounded O(worker-count) reduction windows.
            // Caller/collector/Kappa state remains serial in the original OID order.
            let wave_size = if pool.enabled() { pool.threads() } else { 1 };
            #[cfg(feature = "parallel")]
            let mut worker_diagonals = if pool.enabled() {
                (0..subjects.len().min(wave_size))
                    .map(|_| Diagonals::new(last.offset + last.length, options.window_size))
                    .collect::<Result<Vec<_>>>()?
            } else {
                Vec::new()
            };
            let mut subject_stages = Vec::new();
            for (wave_index, wave) in subjects.chunks(wave_size).enumerate() {
                let mut scheduled = {
                    #[cfg(feature = "parallel")]
                    {
                        use rayon::prelude::*;
                        if pool.enabled() {
                            pool.install(|| {
                                wave.par_iter()
                                    .zip(worker_diagonals.par_iter_mut())
                                    .enumerate()
                                    .map(|(within, (subject, scratch))| {
                                        compute_subject(
                                            subject,
                                            batch,
                                            &params,
                                            options,
                                            &lookup,
                                            scratch,
                                            xdrop,
                                            diagnostics,
                                            db_length,
                                            wave_index * wave_size + within,
                                        )
                                    })
                                    .collect::<Vec<_>>()
                            })
                        } else {
                            Vec::new()
                        }
                    }
                    #[cfg(not(feature = "parallel"))]
                    {
                        Vec::<Result<SubjectComputation>>::new()
                    }
                }
                .into_iter();
                for (within, subject) in wave.iter().enumerate() {
                    let oid = wave_index * wave_size + within;
                    let computed = if let Some(computed) = scheduled.next() {
                        computed?
                    } else {
                        compute_subject(
                            subject,
                            batch,
                            &params,
                            options,
                            &lookup,
                            &mut diagonals,
                            xdrop,
                            diagnostics,
                            db_length,
                            oid,
                        )?
                    };
                    sink.subject_start(oid, batch, &params, computed.prepared_link)?;
                    for (h, contained) in &computed.containment {
                        sink.containment(h, *contained);
                    }
                    let SubjectComputation {
                        seeds,
                        initial,
                        raw,
                        purged,
                        ..
                    } = computed;
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
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:550-552
                    // ```c++
                    //             Blast_HSPListAdjustOddBlastnScores(hsp_list, score_options->gapped_calculation, gap_align->sbp);
                    // #endif
                    //
                    // ```
                    if options.gapped && !raw.is_empty() {
                        sink.preliminary_purge(oid, &raw, &purged)?;
                    }
                    sink.subject(oid, &purged, batch, &params)?;
                    if diagnostics {
                        let mut adjusted = purged.clone();
                        let points = adjust_chunk_hsps(&mut adjusted, &chunk)?;
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:456-487
                        // ```c++
                        //    for (i = 0; i < results1->num_queries; i++) {
                        //        BlastHitList *hitlist = results1->hitlist_array[i];
                        //        Int4 global_query = query_list[i];
                        //        Int4 split_points[NUM_FRAMES];
                        // #ifdef _DEBUG
                        //        ASSERT(i < num_queries);
                        // #endif
                        //
                        //        if (hitlist == NULL) {
                        // #if defined(_DEBUG_VERBOSE)
                        // fprintf(stderr, "No hits to query %d\n", global_query);
                        // #endif
                        //            continue;
                        //        }
                        //
                        //        /* we will be mapping HSPs from the local context to
                        //           their place on the unsplit concatenated query. Once
                        //           that's done, overlapping HSPs need to get merged, and
                        //           to do that we must know the offset within each context
                        //           where the last chunk ended and the current chunk begins */
                        //        for (j = 0; j < contexts_per_query; j++) {
                        //            split_points[j] = -1;
                        //        }
                        //
                        //        for (j = 0; j < contexts_per_query; j++) {
                        //            Int4 local_context = i * contexts_per_query + j;
                        //            if (context_list[local_context] >= 0) {
                        //                split_points[context_list[local_context] % contexts_per_query] =
                        //                                 offset_list[local_context];
                        //            }
                        //        }
                        //
                        // ```
                        let mut by_query = vec![Vec::new(); queries.len()];
                        for h in adjusted {
                            by_query[h.context / 6].push(h);
                        }
                        for (query, hsps) in by_query.into_iter().enumerate() {
                            merge_chunk_hsps(
                                &mut merged[query][oid],
                                hsps,
                                &points,
                                i32::MAX as usize,
                            );
                        }
                        subject_stages.push(SubjectStage {
                            oid,
                            call,
                            parameters: params.clone(),
                            gap_xdrop: xdrop,
                            seeds,
                            initial,
                            raw,
                            purged,
                        });
                    }
                    call += 1;
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:268-271
            // ```c++
            //                 _ASSERT(chunk_data->m_HspStream->GetPointer());
            //                 BlastHSPStreamMerge(split_query_blk->GetCStruct(), i,
            //                                 chunk_data->m_HspStream->GetPointer(),
            //                                 m_InternalData->m_HspStream->GetPointer());
            // ```
            sink.chunk_end(&chunk, is_split)?;
            if diagnostics {
                stages.push(ChunkStage {
                    statistically_valid: true,
                    chunk,
                    subjects: subject_stages,
                });
            }
        }
        // full preparation retains full DNA lengths, masked/nomask bytes and ordinal.
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
        sink.batch_end(&full, &full_params)?;
        if diagnostics {
            result.push(BatchStage {
                query_ordinal: range.start,
                restored: full,
                full_parameters: full_params,
                chunks: stages,
                merged_by_query: merged,
            });
        }
    }
    Ok(result)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:172-179
// ```c++
//     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
//         (*thread)->Run();
//     }
//
//     // ... and wait for the threads to finish
//     Uint8 retv(0);
//     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
//         void * result(0);
// ```
// Only independent seed/extension state belongs to workers. Existing LOSAT
// utils::threading owns pool lifetime; all sink state is replayed by input OID.
struct SubjectComputation {
    seeds: Vec<(i32, i32)>,
    initial: Vec<InitHsp>,
    raw: Vec<PreliminaryHsp>,
    purged: Vec<PreliminaryHsp>,
    containment: Vec<(PreliminaryHsp, bool)>,
    prepared_link: Option<([i32; 2], f64)>,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:172-179
// ```c++
//     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
//         (*thread)->Run();
//     }
//
//     // ... and wait for the threads to finish
//     Uint8 retv(0);
//     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
//         void * result(0);
// ```
// Only independent seed/extension state belongs to workers. Existing LOSAT
// utils::threading owns pool lifetime; all sink state is replayed by input OID.
#[allow(clippy::too_many_arguments)]
fn compute_subject(
    subject: &FastaRecord,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    diagnostics: bool,
    db_length: i64,
    _oid: usize,
) -> Result<SubjectComputation> {
    #[cfg(feature = "blastx-worker-probe")]
    worker_probe("start", _oid);
    let mut containment = Vec::new();
    let mut encoded: Vec<_> = subject
        .sequence
        .iter()
        .copied()
        .map(aa_char_to_ncbistdaa)
        .collect();
    encoded.push(0);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1451-1455
    // ```c++
    //       if (hit_params->link_hsp_params && !kNucleotide &&
    //           !gapped_calculation) {
    //           CalculateLinkHSPCutoffs(program_number, query_info, sbp,
    //             hit_params->link_hsp_params, word_params, db_length,
    //             seq_arg.seq->length);
    // ```
    // Worker-local cutoff preparation precedes WordFinder as in NCBI; only
    // adoption/observer notification is replayed on the serial OID collector.
    let prepared_link = if !options.gapped && options.sum_stats {
        Some(super::even_gap::cutoffs(
            batch,
            params,
            subject.sequence.len() as i32,
            db_length,
        ))
    } else {
        None
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-515
    // ```c++
    //         scan_range[2] = scan_range[1];
    //
    //     while (scan_range[1] <= scan_range[2]) {
    //         /* scan the subject sequence for hits */
    //         hits = scansub(lookup_wrap, subject,
    //                                   offset_pairs, array_size, scan_range);
    //
    //         totalhits += hits;
    //         /* for each hit, */
    //         for (i = 0; i < hits; ++i) {
    //             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
    //             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
    //
    //             /* calculate the diagonal associated with this query-subject pair
    //              */
    //
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
    let mut seeds = Vec::new();
    let initial = word_finder(
        batch,
        params,
        options,
        &encoded,
        lookup,
        diagonals,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:494-500
        // ```c++
        //             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
        //                                    word_params, aux_struct->ewp,
        //                                    aux_struct->offset_pairs,
        //                                    kScanSubjectOffsetArraySize,
        //                                    init_hitlist, ungapped_stats);
        //
        //             if (init_hitlist->total == 0) continue;
        // ```
        |pairs| {
            if diagnostics {
                seeds.extend_from_slice(pairs);
            }
        },
    )?;
    let raw = if options.gapped {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3916-3921
        // ```c++
        //
        //       /* use priate interval tree when recomputing alignments */
        //       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
        //                                         hit_options->min_diag_separation))
        //       {
        //          BlastHSP* new_hsp;
        // ```
        super::preliminary::gapped_observed(
            batch,
            params,
            options,
            &encoded,
            &initial,
            xdrop,
            &mut |h, contained| containment.push((h.clone(), contained)),
        )?
    } else {
        ungapped(batch, &initial)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:552-563
    // ```c++
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
    // ```
    let mut purged = raw.clone();
    if options.gapped {
        purge_endpoints(&mut purged);
    } else {
        sort_by_score(&mut purged);
    }
    #[cfg(feature = "blastx-worker-probe")]
    worker_probe("finish", _oid);
    Ok(SubjectComputation {
        seeds,
        initial,
        raw,
        purged,
        containment,
        prepared_link,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1659-1664
// ```c++
//     /* Use a local diagnostics structure, because the one passed in an input
//       argument can be shared between multiple threads, so we don't want to pass
//       it to the engine and have a lot of mutex contention. */
//     BlastDiagnostics* local_diagnostics = Blast_DiagnosticsInit();
//
//     if ((status =
// ```
// Diagnostic-only host clock/worker identity; never read by search or reporting.
#[cfg(feature = "blastx-worker-probe")]
fn worker_probe(event: &str, oid: usize) {
    let ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("host clock after epoch")
        .as_nanos();
    crate::utils::threading::probe_write(format!(
        "[blastx-worker-probe] event={event} oid={oid} thread={:?} ns={ns}\n",
        std::thread::current().id()
    ));
}
