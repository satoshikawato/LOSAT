//! BLASTX native serial local-subject preparation, search and reporting.
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:49-62
// ```c++
//                                   "Translated Query-Protein Subject BLAST"));
//     const bool kQueryIsProtein = false;
//     m_Args.push_back(arg);
//     m_ClientId = kProgram + " " + CBlastVersion().Print();
//
//     static const char kDefaultTask[] = "blastx";
//     SetTask(kDefaultTask);
//     set<string> tasks;
//     tasks.insert(kDefaultTask);
//     tasks.insert("blastx-fast");
//     arg.Reset(new CTaskCmdLineArgs(tasks, kDefaultTask));
//     m_Args.push_back(arg);
//
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
// ```
pub mod args;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:319-326
// ```c++
//     CFastaReader::TFlags flags = m_Config.GetBelieveDeflines() ?
//                                     CFastaReader::fParseRawID:
//                                     (CFastaReader::fNoParseID |
//                                      CFastaReader::fDLOptional);
//
//     // Allow CFastaReader fSkipCheck flag to be set based
//     // on new CBlastInputSourceConfig property - GetSkipSeqCheck() -RMH-
//     flags += ( m_Config.GetSkipSeqCheck() ? CFastaReader::fSkipCheck : 0 );
// ```
pub mod input;
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:577-587
// ```c++
//                 for (unsigned int i = 0; i < kNumContexts; i++) {
//                     if (qinfo->contexts[ctx_index + i].query_length <= 0) {
//                         continue;
//                     }
//
//                     int offset = qinfo->contexts[ctx_index + i].query_offset;
//                     BLAST_GetTranslation(sequence.data.get() + 1,
//                                          seqbuf_rev,
//                                          na_length,
//                                          qinfo->contexts[ctx_index + i].frame,
//                                          & buf.get()[offset], gc);
// ```
pub mod query_setup;
pub use args::BlastxArgs;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:960-985
// ```c++
//    if((status=BlastExtensionParametersNew(program_number, ext_options, sbp,
//                                query_info, ext_params)) != 0)
//    {
//       *eff_len_params = BlastEffectiveLengthsParametersFree(*eff_len_params);
//       *score_params = BlastScoringParametersFree(*score_params);
//       *ext_params = BlastExtensionParametersFree(*ext_params);
//       return status;
//    }
//
//    if (sbp->gbp) {
//        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
//        if (Blast_SubjectIsTranslated(program_number)) {
//            min_subject_length/=3;
//        }
//    } else {
//        min_subject_length = (Int4) (total_length/num_seqs);
//    }
//
//    if(min_subject_length <=0) {
// 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
//    }
//
//    if ((status = BlastHitSavingParametersNew(program_number, hit_options, sbp,
// 		                                     query_info, min_subject_length,
// 		                                     (*ext_params)->options->compositionBasedStats,
// 		                                     hit_params)) != 0){
// ```
pub mod parameters;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:200-234
// ```c++
// Int2 BlastAaWordFinder(BLAST_SequenceBlk * subject,
//                        BLAST_SequenceBlk * query,
//                        BlastQueryInfo * query_info,
//                        LookupTableWrap * lut_wrap,
//                        Int4 ** matrix,
//                        const BlastInitialWordParameters * word_params,
//                        Blast_ExtendWord * ewp,
//                        BlastOffsetPair * NCBI_RESTRICT offset_pairs,
//                        Int4 offset_array_size,
//                        BlastInitHitList * init_hitlist,
//                        BlastUngappedStats * ungapped_stats)
// {
//     Int2 status = 0;
//
//     if (ewp->diag_table->multiple_hits) {
//         status = s_BlastAaWordFinder_TwoHit(subject, query,
//                                             lut_wrap, ewp,
//                                             matrix,
//                                             word_params,
//                                             query_info,
//                                             offset_pairs,
//                                             offset_array_size,
//                                             init_hitlist, ungapped_stats);
//     } else {
//         status = s_BlastAaWordFinder_OneHit(subject, query,
//                                             lut_wrap, ewp,
//                                             matrix,
//                                             word_params,
//                                             query_info,
//                                             offset_pairs,
//                                             offset_array_size,
//                                             init_hitlist, ungapped_stats);
//     }
//
//     Blast_InitHitListSortByScore(init_hitlist);
// ```
pub mod seed;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:522-545
// ```c++
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
// ```
pub mod preliminary;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:50-62
// ```c++
// CQuerySplitter::CQuerySplitter(CRef<IQueryFactory> query_factory,
//                                const CBlastOptions* options)
//     : m_QueryFactory(query_factory), m_Options(options), m_NumChunks(0),
//     m_LocalQueryData(0), m_TotalQueryLength(0), m_ChunkSize(0)
// {
//     m_ChunkSize = SplitQuery_GetChunkSize(m_Options->GetProgram());
//     m_LocalQueryData = m_QueryFactory->MakeLocalQueryData(m_Options);
//     m_TotalQueryLength = m_LocalQueryData->GetSumOfSequenceLengths();
//     m_NumChunks = SplitQuery_CalculateNumChunks(m_Options->GetProgramType(),
//         &m_ChunkSize, m_TotalQueryLength, m_LocalQueryData->GetNumQueries());
//     /* No split for ungapped mode JIRA SB-1082 */
//     if (!options->GetGappedMode()) m_NumChunks = 1;
//     x_ExtractCScopesAndMasks();
// ```
pub mod split;

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
pub mod search;

// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/split_query_unit_test.cpp:1901-1931
// ```c++
// BOOST_AUTO_TEST_CASE(CalculateNumberChunks)
// {
//     EBlastProgramType program = eBlastTypeBlastx;
//     size_t chunk_size = 10002;
//     Uint4 retval = SplitQuery_CalculateNumChunks(program,
//                        &chunk_size, 10240000, 1);
//     BOOST_REQUIRE_EQUAL(1055, retval);
//
//     retval = SplitQuery_CalculateNumChunks(eBlastTypeBlastx,
//                        &chunk_size, chunk_size/2, 1);
//
//     BOOST_REQUIRE_EQUAL(1, retval);
//
//     retval = SplitQuery_CalculateNumChunks(program,
//                        &chunk_size,
//                        3*chunk_size-2*SplitQuery_GetOverlapChunkSize(program), 1);
//
//     BOOST_REQUIRE_EQUAL(3, retval);
//
//     retval = SplitQuery_CalculateNumChunks(program,
//                        &chunk_size,
//                        1+2*chunk_size+SplitQuery_GetOverlapChunkSize(program), 1);
//
//     BOOST_REQUIRE_EQUAL(2, retval);
// }
//
// BOOST_AUTO_TEST_CASE(InvalidChunkSizeBlastx)
// {
//     CAutoEnvironmentVariable tmp_env("CHUNK_SIZE", "40000");
//     BOOST_REQUIRE_THROW(SplitQuery_GetChunkSize(blast::eBlastx), CBlastException);
// }
// ```
#[cfg(test)]
mod stage_c_tests;

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
pub mod linking;

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
pub mod statistics;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1498
// ```c++
//     } else if (ext_params->options->compositionBasedStats > 0 ||
//                ext_params->options->eTbackExt == eSmithWatermanTbck) {
//         Uint4 num_threads = MAX(1, thread_data->num_elems);
//         results = Blast_HSPResultsNew(query_info->num_queries);
//         /* FIXME partial sequence fetching/translation could lead to fence hit
//            and seg fault */
//         retval =
//                 Blast_RedoAlignmentCore_MT(program_number,
//                         num_threads,
//                         query, query_info, sbp,
//                         NULL, seq_src, default_db_genetic_code,
//                         NULL, hsp_stream, score_params, ext_params,
//                         hit_params, psi_options, results);
// ```
pub mod kappa;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636
// ```c
//                 status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
//                                              subject_maybe_biased);
// ```
// EXPERIMENT (LOSAT_X_BXSEGMEMO): the redo-stage subject SEG memoised per subject.
mod x_seg_memo;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1508
// ```c++
//     } else if (ext_params->options->compositionBasedStats > 0 ||
//                ext_params->options->eTbackExt == eSmithWatermanTbck) {
//         Uint4 num_threads = MAX(1, thread_data->num_elems);
//         results = Blast_HSPResultsNew(query_info->num_queries);
//         /* FIXME partial sequence fetching/translation could lead to fence hit
//            and seg fault */
//         retval =
//                 Blast_RedoAlignmentCore_MT(program_number,
//                         num_threads,
//                         query, query_info, sbp,
//                         NULL, seq_src, default_db_genetic_code,
//                         NULL, hsp_stream, score_params, ext_params,
//                         hit_params, psi_options, results);
//     } else {
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
pub mod traceback;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:105-114
// ```c++
//    /* Rearrange HSPs into multiple hit lists if more than one query */
//    if (results->num_queries > 1) {
//       BlastHSP* hsp;
//       BlastHSPList** hsp_list_array;
//       BlastHSPList* tmp_hsp_list;
//       Int4 index;
//
//       hsp_list_array = calloc(results->num_queries, sizeof(BlastHSPList*));
//       if (hsp_list_array == NULL)
//          return -1;
// ```
pub mod results;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:268-271
// ```c++
//                 _ASSERT(chunk_data->m_HspStream->GetPointer());
//                 BlastHSPStreamMerge(split_query_blk->GetCStruct(), i,
//                                 chunk_data->m_HspStream->GetPointer(),
//                                 m_InternalData->m_HspStream->GetPointer());
// ```
pub mod runtime;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:386-394
// ```c++
//         } else if (filt_opts->culling_opts &&
//                    (filt_opts->culling_stage & eTracebackSearch)) {
//             BlastHSPCullingParams* params =
//                 BlastHSPCullingParamsNew(opts_memento->m_HitSaveOpts,
//                      filt_opts->culling_opts,
//                      opts_memento->m_ExtnOpts->compositionBasedStats,
//                      opts_memento->m_ScoringOpts->gapped_calculation);
//             BlastHSPPipeInfo_Add(&pipe_info,
//                                  BlastHSPCullingPipeInfoNew(params));
// ```
mod culling;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1777-1784
// ```c++
//
//     /* Link up the HSP's for this hsp_list. */
//     if (link_hsp_params->longest_intron <= 0) {
//         s_BlastEvenGapLinkHSPs(program_number, hsp_list, query_info,
//                               subject_length, sbp, link_hsp_params,
//                               gapped_calculation);
//         /* The HSP's may be in a different order than they were before,
//            but hsp contains the first one. */
// ```
mod even_gap;

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:279-297
// ```c++
//                 results = lcl_blast.Run();
// 	        BLAST_PROF_STOP( APP.LOOP.BLAST );
//             }
// 	    BLAST_PROF_START( APP.LOOP.FMT );
//             if (fmt_args->ArchiveFormatRequested(args)) {
//                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
//                 m_Bah.ResetMessages();
//             } else {
//                 BlastFormatter_PreFetchSequenceData(*results, scope,
//                 		                            fmt_args->GetFormattedOutputChoice());
//             	ITERATE(CSearchResultSet, result, *results) {
//                	    formatter.PrintOneResultSet(**result, query_batch);
//             	}
//             }
// 	    BLAST_PROF_STOP( APP.LOOP.FMT );
// 	    batch_num++;
//         }
//         BLAST_PROF_START( APP.POST );
//         formatter.PrintEpilog(opt);
// ```
pub mod report;

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
// ```c++
//         formatter.PrintProlog();
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// ```c++
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:289-291
// ```c++
//             	ITERATE(CSearchResultSet, result, *results) {
//                	    formatter.PrintOneResultSet(**result, query_batch);
//             	}
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
// ```c++
//         formatter.PrintEpilog(opt);
// ```
pub mod native;

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// ```c++
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
// ```
pub(crate) mod web;
