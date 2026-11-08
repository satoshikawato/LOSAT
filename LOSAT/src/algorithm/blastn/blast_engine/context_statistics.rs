//! The statistics of a batch's query contexts, set up before its subjects are searched: the
//! Karlin-Altschul blocks of every context (`karlin_blocks`), which contexts are valid, the checks
//! of the set-up on the subjects and the gap costs, and the effective search space of every
//! context (`effective_search_spaces`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:455-532
//! (`BlastSetup_ScoreBlkInit`), which calls c++/src/algo/blast/core/blast_stat.c:2736-2829
//! (`Blast_ScoreBlkKbpUngappedCalc`) and blast_setup.c:40-128 (`Blast_ScoreBlkKbpGappedCalc`);
//! blast_setup.c:887-999 (`BLAST_GapAlignSetUp`), which calls blast_setup.c:699-850
//! (`BLAST_CalcEffLengths`).
//!
//! This module is a child of `run` for the types it keeps private (`QueryContext`,
//! `SubjectSource`, `BatchStage`, `ShownSubjects`, `BatchOutputs`). The functions are statements
//! of `search_query_batch` moved here in their order; what they read is a struct of references or
//! `Copy` values, and what they make is what `search_query_batch` then reads under the same names.
//! An error of the set-up writes the outfmt 0 prolog of a file sink of the first batch through the
//! batch's outputs before it is returned, as it did in `search_query_batch`.

use super::{
    context_lengths, database, database_search_error, prolog_with_report, reborrow_outputs,
    write_format_prolog, write_pairwise_prologs, BatchOutputs, BatchStage, BlastnArgs,
    ContextKarlin, ErrorAfterInvalidBatches, KarlinParams, NuclScoringSpec, QueryBatching,
    QueryContext, SequenceData, ShownSubjects, SubjectSource, FIRST_BATCH_OUTPUTS,
};
use crate::algorithm::blastn::scoring::{
    check_greedy_gap_costs, context_blocks, context_blocks_ungapped, karlin_error,
};
use anyhow::Result;

/// What `context_statistics` reads: the batch's contexts with the ungapped blocks of their
/// compositions, the subjects, the stage of the batch, and what the outfmt 0 prolog of an error
/// shows.
pub(super) struct ContextInputs<'a> {
    pub(super) args: &'a BlastnArgs,
    /// The ungapped block of each context's composition (`None`: none could be computed).
    pub(super) context_ungapped: &'a [Option<KarlinParams>],
    pub(super) scoring_spec: NuclScoringSpec,
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) query_lengths: &'a [usize],
    pub(super) batching: QueryBatching,
    pub(super) first_batch: bool,
    pub(super) stage: BatchStage<'a>,
    pub(super) subjects: SubjectSource<'a>,
    pub(super) seq_data: &'a SequenceData<'a>,
    pub(super) megablast: bool,
    pub(super) shown: &'a ShownSubjects,
}

/// The statistics of the batch's contexts, which `search_query_batch` reads under the same names.
pub(super) struct ContextStatistics {
    /// The blocks of each context, `None` for an invalid context.
    pub(super) context_karlin: Vec<Option<ContextKarlin>>,
    pub(super) round_down_evalue_score: bool,
    pub(super) context_valid: Vec<bool>,
    pub(super) any_valid_context: bool,
    /// The blocks the search reads for each context (an invalid context's slot holds those of a
    /// valid one).
    pub(super) search_karlin: Vec<ContextKarlin>,
    /// The effective search space of each context.
    pub(super) query_eff_searchsp: Vec<i64>,
}

/// The Karlin-Altschul blocks, the validity and the effective search space of every context of the
/// batch, after the checks of the set-up; an error writes what NCBI writes before it (the outfmt 0
/// prolog of a file sink of the first batch) through `outputs`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:455-532
/// (`BlastSetup_ScoreBlkInit`, which `BLAST_MainSetUp` calls at :652 for
/// `CSetupFactory::CreateScoreBlock`, c++/src/algo/blast/api/setup_factory.cpp:121-158) and
/// blast_setup.c:887-999 (`BLAST_GapAlignSetUp`, which c++/src/algo/blast/core/blast_engine.c:1665
/// calls first in `Blast_RunPreliminarySearchWithInterrupt`, before the subjects).
pub(super) fn context_statistics<const UNGAPPED: bool>(
    inputs: &ContextInputs<'_>,
    mut outputs: BatchOutputs<'_, '_>,
) -> Result<ContextStatistics> {
    let (context_karlin, round_down_evalue_score) =
        karlin_blocks::<UNGAPPED>(inputs, reborrow_outputs(&mut outputs))?;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:324-332
    // ```c
    //        if (!(query_info->contexts[context].is_valid)) {
    //           /* either this context was never valid, or it was
    //              valid at the beginning of the search but is not
    //              valid now. The latter means that ungapped
    //              alignments can still occur to this context,
    //              so we set the cutoff score to be infinite */
    //           curr_cutoffs->cutoff_score = INT4_MAX;
    //           continue;
    //        }
    // ```
    // An invalid context has an infinite cutoff (`batch_setup::batch_cutoffs`), so it gets no hit.
    // Its slot holds the blocks of a valid context (or placeholder blocks when none is valid),
    // which no hit reads.
    let context_valid: Vec<bool> = context_karlin.iter().map(Option::is_some).collect();
    let any_valid_context = context_valid.iter().any(|&valid| valid);
    let (search_karlin, query_eff_searchsp) = effective_search_spaces::<UNGAPPED>(
        inputs,
        outputs,
        &context_karlin,
        &context_valid,
        any_valid_context,
    )?;
    Ok(ContextStatistics {
        context_karlin,
        round_down_evalue_score,
        context_valid,
        any_valid_context,
        search_karlin,
        query_eff_searchsp,
    })
}

/// The ungapped and gapped Karlin-Altschul blocks of every context (`None` for an invalid
/// context) and whether the e-values round the score down; the error of scoring options without
/// Karlin-Altschul values.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:455-532
/// (`BlastSetup_ScoreBlkInit`), which calls c++/src/algo/blast/core/blast_stat.c:2736-2829
/// (`Blast_ScoreBlkKbpUngappedCalc`) and, for a gapped search, blast_setup.c:40-128
/// (`Blast_ScoreBlkKbpGappedCalc`).
fn karlin_blocks<const UNGAPPED: bool>(
    inputs: &ContextInputs<'_>,
    outputs: BatchOutputs<'_, '_>,
) -> Result<(Vec<Option<ContextKarlin>>, bool)> {
    let &ContextInputs {
        context_ungapped,
        scoring_spec,
        query_contexts,
        query_lengths,
        batching,
        first_batch,
        megablast,
        shown,
        ..
    } = inputs;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2758-2792
    // ```c
    //       Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
    //       sbp->sfp[context] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
    //       BlastScoreFreqCalc(sbp, sbp->sfp[context], rfp, stdrfp);
    //       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
    //       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    //       if (loop_status) {
    //           contexts[context].is_valid = FALSE;
    // ```
    // Every query context (strand) has the ungapped block of its composition (kbp_std:
    // the gap trigger and the ungapped X-drop) and a gapped block (kbp_gap: e-values,
    // bit scores, cutoffs and the length adjustment); `scoring.rs` has the NCBI
    // references. NCBI raises an unsupported scoring system when it sets up the first
    // batch with a valid query, after the outfmt 0 prolog.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:513-529
    // ```c
    //     if (scoring_options->gapped_calculation) { ... s_BlastScoreBlkGappedFill ... }
    // ```
    // An ungapped (`-ungapped`) search has the ungapped blocks only (`context_blocks_ungapped`).
    let context_blocks_result = if UNGAPPED {
        Ok(context_blocks_ungapped(context_ungapped, &scoring_spec))
    } else {
        context_blocks(context_ungapped, &scoring_spec)
    };
    let (context_karlin, round_down_evalue_score) = match context_blocks_result {
        Ok(blocks) => blocks,
        Err(message) => {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/setup_factory.cpp:170-172
            // ```c
            //     Blast_Message2TSearchMessages(blast_msg.Get(), query_info, search_messages);
            //     if (status != 0 &&
            //         (!(blast_msg.Get()) || (blast_msg.Get() && blast_msg.Get()->severity == eBlastSevError)))
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2784-2791
            // ```c
            //           contexts[context].is_valid = FALSE;
            //     ...
            //              Blast_MessageWrite(blast_message, eBlastSevWarning, context,
            //              kBlastErrMsg_CantCalculateUngappedKAParams);
            // ```
            // NCBI raises the error only when the first message of the batch is it: an
            // invalid context puts its warning first, and NCBI then goes on without the
            // gapped blocks (it crashes; with a first batch of invalid queries only, a later
            // batch raises it). LOSAT reproduces the error only for a first batch of valid
            // queries (the error does not depend on the queries, so the first batch with a valid
            // query meets it; NCBI writes the results and warnings of the batches before it).
            // The contexts of a query without letters are invalid from their set-up and put
            // no message here (`context_ungapped_blocks`).
            if query_contexts
                .iter()
                .zip(context_ungapped)
                .any(|(context, block)| context.len != 0 && block.is_none())
            {
                anyhow::bail!(
                    "these scoring options have no Karlin-Altschul values, and the query batch (up to {} residues) has an invalid query; NCBI BLAST+ goes on without the gapped values and crashes, which is not supported by LOSAT's BLASTN",
                    batching.batch_size
                );
            }
            // A query without letters (also a `-query_loc` interval that holds none) keeps its
            // warning in the error (`karlin_error`).
            let empty_queries = query_lengths.iter().map(|&length| length == 0);
            // A batch after the first: every batch before it had only invalid queries (one
            // with a valid query would have raised the error), and NCBI has written their
            // reports (`run_in_pool`).
            if !first_batch {
                return Err(ErrorAfterInvalidBatches(karlin_error(&message, empty_queries)).into());
            }
            // The prolog of a file sink (the others were written before the search).
            let (outputs, output_formats) = outputs.expect(FIRST_BATCH_OUTPUTS);
            for (format, &output_format) in outputs.formats.iter_mut().zip(output_formats) {
                if prolog_with_report(format, output_format) {
                    write_format_prolog(format, output_format, megablast, shown)?;
                }
            }
            return Err(karlin_error(&message, empty_queries));
        }
    };
    Ok((context_karlin, round_down_evalue_score))
}

/// The checks of the set-up on the subjects (a database without sequences, the average subject
/// length) and on the greedy gap costs, then the blocks the search reads for each context and the
/// effective search space of each context.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:887-999
/// (`BLAST_GapAlignSetUp`), which calls c++/src/algo/blast/core/blast_setup.c:699-850
/// (`BLAST_CalcEffLengths`).
fn effective_search_spaces<const UNGAPPED: bool>(
    inputs: &ContextInputs<'_>,
    outputs: BatchOutputs<'_, '_>,
    context_karlin: &[Option<ContextKarlin>],
    context_valid: &[bool],
    any_valid_context: bool,
) -> Result<(Vec<ContextKarlin>, Vec<i64>)> {
    let &ContextInputs {
        args,
        query_contexts,
        first_batch,
        stage,
        subjects,
        seq_data,
        megablast,
        shown,
        ..
    } = inputs;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:908-939
    // ```c
    //    if (seq_src) {
    //       total_length = BlastSeqSrcGetTotLenStats(seq_src);
    //       if (total_length <= 0)
    //           total_length = BlastSeqSrcGetTotLen(seq_src);
    //       ...
    //       if (total_length > 0) {
    //           num_seqs = BlastSeqSrcGetNumSeqsStats(seq_src);
    //           if (num_seqs <= 0)
    //               num_seqs = BlastSeqSrcGetNumSeqs(seq_src);
    //       } else {
    //           /* Not a database search; each subject sequence is considered
    //              individually */
    //           Int4 oid = 0;  /* Get length of first sequence. */
    //           if ( (total_length = BlastSeqSrcGetSeqLen(seq_src, (void*) &oid)) < 0) {
    //               total_length = -1;
    //               num_seqs = -1;
    //           }
    //           num_seqs = 1;
    //       }
    //    }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:969-980
    // ```c
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
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_message.c:215-219
    // ```c
    //     case BLASTERR_SUBJECT_LENGTH_INVALID:
    //         new_msg->message = strdup("The average subject length is too short");
    //         new_msg->severity = eBlastSevFatal;
    //         new_msg->context = context;
    //         break;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:310-314
    // ```c
    //             retval = CPrelimSearchRunner(*m_InternalData, opts_memento.get())();
    //             if (retval) {
    //                 NCBI_THROW(CBlastException, eCoreBlastError,
    //                            BlastErrorCode2String(retval));
    //             }
    // ```
    // The subject set has no BLAST database statistics, and BLASTN has no Gumbel block
    // (blast_setup.c:475-480): a batch with a valid context stops NCBI when the subjects
    // have fewer letters than subjects (subjects without letters count with none). A batch
    // without a valid context is not searched (`CheckInternalData`), so it does not.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:448-458
    // ```c++
    // static Int4
    // s_SeqDbGetSeqLen(void* seqdb_handle, void* args)
    // {
    //     Int4* oid = (Int4*) args;
    //
    //     if (!seqdb_handle || !oid)
    //        return BLAST_SEQSRC_ERROR;
    //
    //     CSeqDB & seqdb = **(TSeqDBData *) seqdb_handle;
    //     return seqdb.GetSeqLength(*oid);
    // }
    // ```
    // A database without sequences (makeblastdb builds none; a volume whose index holds
    // none) has a total length of 0, so the setup above asks its reader for the length of
    // ordinal id 0, which it does not have: NCBI's search stage reports the reader's error
    // as a database file that cannot be read (local_blast.cpp:240-251). Shown with the
    // oracle (campaign dbsearch-20261012, `$BUILD_ROOT/dbsearch/db4/smoke/empty`: the
    // nucleotide fixture's index with no sequence): the outfmt 0 prolog, then `BLAST engine
    // error: Database memory map file error`, exit status 3, for outfmt 0, 6 and 7.
    if any_valid_context {
        if let SubjectSource::Database(database) = subjects {
            if database.metadata().db_num_seqs == 0 {
                return Err(database_search_error(
                    database::PreliminaryReadError("Database memory map file error").into(),
                    first_batch,
                    outputs,
                    megablast,
                    shown,
                ));
            }
        }
    }
    if any_valid_context {
        let total_length = seq_data.db_len_total as i64;
        let (total_length, num_seqs) = if total_length > 0 {
            (total_length, seq_data.db_num_seqs as i64)
        } else {
            // A database with a total length of 0 holds sequences without letters only:
            // the length of its first sequence (`BlastSeqSrcGetSeqLen` of ordinal id 0) is 0.
            (
                subjects.prepared().map_or(0, |prepared| {
                    if prepared.records.is_empty() {
                        0
                    } else {
                        prepared.letters().len(0)
                    }
                }) as i64,
                1,
            )
        };
        if total_length / num_seqs <= 0 {
            let error = anyhow::Error::from(crate::cli::NativeError {
                exit: 3,
                message: "BLAST engine error: The average subject length is too short\n"
                    .to_string(),
            });
            // A batch after the first: the batches before had no valid query (one with a
            // valid query would have raised the error), and NCBI has written their reports
            // (`run_in_pool`).
            if !first_batch {
                return Err(ErrorAfterInvalidBatches(error).into());
            }
            let (outputs, output_formats) = outputs.expect(FIRST_BATCH_OUTPUTS);
            write_pairwise_prologs(outputs, output_formats, megablast, shown, false)?;
            return Err(error);
        }
    }
    // Gap costs that the tables accept but LOSAT's greedy extension does not (`scoring.rs`);
    // with no valid context, NCBI extends nothing.
    if any_valid_context && !UNGAPPED {
        check_greedy_gap_costs(args)?;
    }
    let search_karlin: Vec<ContextKarlin> = {
        let fallback = context_karlin
            .iter()
            .flatten()
            .next()
            .copied()
            .unwrap_or(ContextKarlin {
                ungapped: KarlinParams::default(),
                gapped: KarlinParams::default(),
            });
        context_karlin
            .iter()
            .map(|blocks| blocks.unwrap_or(fallback))
            .collect()
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-846
    // ```c
    // BLAST_ComputeLengthAdjustment(..., query_length, db_length, db_num_seqs, &length_adjustment);
    // effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
    // if (effective_db_length <= 0) effective_db_length = 1;
    // effective_search_space = effective_db_length * (query_length - length_adjustment);
    // query_info->contexts[index].eff_searchsp = effective_search_space;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-721
    // ```c
    //    /* use overriding value from effective lengths options or the real value
    //       from effective lengths parameters. */
    //    if (eff_len_options->db_length > 0)
    //       db_length = eff_len_options->db_length;
    //    else
    //       db_length = eff_len_params->real_db_length;
    // ```
    // A `-dbsize` above 0 replaces the total length in the effective lengths (with `-db`
    // and with `-subject`); the number of sequences, the average subject length above, the
    // batch mixer and the reports keep the real numbers.
    let db_len_total_i64 = if args.db_length() > 0 {
        args.db_length()
    } else {
        seq_data.db_len_total as i64
    };
    let db_num_seqs_i64 = seq_data.db_num_seqs as i64;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:778-786
    // ```c
    //       Int8 effective_search_space =
    //           s_GetEffectiveSearchSpaceForContext(eff_len_options, index,
    //                                               blast_message);
    //     ...
    //       if (query_info->contexts[index].is_valid &&
    //           ((query_length = query_info->contexts[index].query_length) > 0) ) {
    // ```
    // An invalid context keeps the search space of the options, which is 0, except in a query
    // chunk: NCBI sets the options to the search spaces of the batch's contexts
    // (`SplitQuery_SetEffectiveSearchSpace`), and a chunk's context takes the one at its own
    // index in the chunk, not its context in the batch (0, for an invalid context of the
    // batch, has the chunk compute its own).
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:177-182
    // ```c
    //     vector<Int8> eff_searchsp;
    //     for (size_t index = 0; index <= (size_t)qinfo->last_context; index++) {
    //         eff_searchsp.push_back(calc.GetEffSearchSpaceForContext(index));
    //     }
    //     options->SetEffectiveSearchSpace(eff_searchsp);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:690-693
    // ```c
    //     } else if (eff_len_options->num_searchspaces > 1) {
    //         ASSERT(context_index < eff_len_options->num_searchspaces);
    //         retval = eff_len_options->searchsp_eff[context_index];
    //     } else {
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:826-846
    // ```c
    //          if (effective_search_space == 0) {
    //          ...
    //              effective_search_space = effective_db_length *
    //                              (query_length - length_adjustment);
    //          }
    //       }
    //       query_info->contexts[index].eff_searchsp = effective_search_space;
    // ```
    let query_eff_searchsp: Vec<i64> = query_contexts
        .iter()
        .zip(&search_karlin)
        .zip(context_valid)
        .enumerate()
        .map(|(index, ((ctx, blocks), &valid))| {
            let options_searchsp = match stage {
                BatchStage::Search { .. } => args.options_searchsp(),
                BatchStage::ChunkPrelim {
                    batch_eff_searchsp, ..
                } => batch_eff_searchsp[index],
            };
            if !valid || options_searchsp != 0 {
                return Ok(options_searchsp);
            }
            context_lengths(
                ctx.len as i64,
                db_len_total_i64,
                db_num_seqs_i64,
                &blocks.gapped,
            )
            .map(|(_, effective_search_space)| effective_search_space)
        })
        .collect::<Result<_>>()?;
    Ok((search_karlin, query_eff_searchsp))
}
