//! Comparison-only diagnostic; no public BLASTX search or report success.
use anyhow::{bail, Result};
use LOSAT::{
    algorithm::blastx::{
        input::read_fasta,
        kappa::KappaState,
        linking::LinkedHspList,
        parameters::subject_parameters,
        search::search_preliminary,
        statistics::{
            evaluate_gapped, identity_and_positive, normalize_scores, reap_final, reap_preliminary,
        },
        traceback::ordinary_traceback,
    },
    cli::{Cli, Commands},
};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1802-1810
// ```c++
//     /* Sort the HSP array by score */
//     Blast_HSPListSortByScore(hsp_list);
//
//     /* Find and fill the best e-value */
//     hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
//     for (index = 1; index < hsp_list->hspcnt; ++index) {
//         if (hsp_list->hsp_array[index]->evalue < hsp_list->best_evalue)
//             hsp_list->best_evalue = hsp_list->hsp_array[index]->evalue;
//     }
// ```
fn dump(stage: &str, call: usize, oid: usize, list: &LinkedHspList) {
    if list.hsps.is_empty() {
        return;
    }
    println!("{stage}_COUNT\t{call}\t{oid}\t{}", list.hsps.len());
    for (i, h) in list.hsps.iter().enumerate() {
        let p = &h.hsp;
        println!(
            "{stage}\t{call}\t{oid}\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:016x}",
            p.context,
            p.score,
            p.frame,
            p.q_start,
            p.q_end,
            p.q_gapped_start,
            p.s_start,
            p.s_end,
            p.s_gapped_start,
            0,
            h.num,
            h.evalue.to_bits()
        );
    }
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
fn main() -> Result<()> {
    let cli: Cli = LOSAT::cli::try_parse_from(std::env::args_os())?;
    let Commands::Blastx(args) = cli.command else {
        bail!("BLASTX diagnostic only")
    };
    let options = args.resolve()?;
    let query = read_fasta(&args.query, false, args.lcase_masking)?;
    let subject = read_fasta(&args.subject, true, false)?;
    let db_length = subject.iter().map(|s| s.sequence.len() as i64).sum();
    for batch in search_preliminary(&query, &subject, &options)? {
        let mut final_call = 0;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:972-989
        // ```c++
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
        //        return status;
        //    }
        //
        //    /* To initialize the gapped alignment structure, we need to know the
        // ```
        let mut full_parameters = batch.full_parameters.clone();
        subject_parameters(
            &batch.restored,
            &mut full_parameters,
            &options,
            db_length,
            subject.len() as i64,
            subject
                .iter()
                .map(|s| s.sequence.len())
                .min()
                .expect("protein subject"),
            None,
        )?;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3360-3367
        // ```c++
        //             matrix_tld[i] = sbp_tld[i]->psi_matrix->pssm->data;
        //         } else {
        //             matrix_tld[i] = sbp_tld[i]->matrix->data;
        //         }
        //         /**** Validate parameters *************/
        //         if (matrix_tld[i] == NULL) {
        //             goto function_cleanup;
        //         }
        // ```
        let mut kappa_state = if options.composition == 2 {
            Some(KappaState::new(
                &batch.restored,
                &full_parameters,
                &options,
            )?)
        } else {
            None
        };
        anyhow::ensure!(
            batch.chunks.len() == 1,
            "BLASTX Kappa checkpoint requires an unsplit query batch"
        );
        for chunk in &batch.chunks {
            for s in &chunk.subjects {
                let mut linked = evaluate_gapped(
                    &s.purged,
                    &chunk.chunk.prepared,
                    &s.parameters,
                    &options,
                    subject[s.oid].sequence.len() as i32,
                    db_length,
                    1.0,
                )?;
                reap_preliminary(&mut linked, &options);
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3541-3543
                // ```c++
                //
                //                 query_index = localMatch->query_index;
                //                 context_index = query_index * numFrames;
                // ```
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:91-96
                // ```c++
                // static int s_SortHSPListByOid(const void *x, const void *y)
                // {
                //         BlastHSPList **xx = (BlastHSPList **)x;
                //             BlastHSPList **yy = (BlastHSPList **)y;
                //                 return (*yy)->oid - (*xx)->oid;
                // }
                // ```
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:173-202
                // ```c++
                //            Int4 alloc = MAX(num_hsplists + hitlist->hsplist_count + 100,
                //                             2 * hsp_stream->num_hsplists_alloc);
                //            hsp_stream->num_hsplists_alloc = alloc;
                //            hsp_stream->sorted_hsplists = (BlastHSPList **)realloc(
                //                                          hsp_stream->sorted_hsplists,
                //                                          alloc * sizeof(BlastHSPList *));
                //        }
                //
                //        for (j = k = 0; j < hitlist->hsplist_count; j++) {
                //
                //            BlastHSPList *hsplist = hitlist->hsplist_array[j];
                //            if (hsplist == NULL)
                //                continue;
                //
                //            hsplist->query_index = i;
                //            hsp_stream->sorted_hsplists[num_hsplists + k] = hsplist;
                //            k++;
                //        }
                //
                //        hitlist->hsplist_count = 0;
                //        num_hsplists += k;
                //    }
                //
                //    /* sort in order of decreasing subject OID. HSPLists will be
                //       read out from the end of hsplist_array later */
                //
                //    hsp_stream->num_hsplists = num_hsplists;
                //    if (num_hsplists > 1) {
                //       qsort(hsp_stream->sorted_hsplists, num_hsplists,
                //                     sizeof(BlastHSPList *), s_SortHSPListByOid);
                // ```
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:320-329
                // ```c++
                //    } else {
                //        /* return the next HSPlist out of the collection stored */
                //
                //        if (!hsp_stream->num_hsplists)
                //           return kBlastHSPStream_Eof;
                //
                //        *hsp_list_out =
                //            hsp_stream->sorted_hsplists[--hsp_stream->num_hsplists];
                //
                //    }
                // ```
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:599-613
                // ```c++
                //    hsplist = hsp_stream->sorted_hsplists[num_hsplists - 1];
                //    target_oid = hsplist->oid;
                //
                //    for (i = 0; i < num_hsplists; i++) {
                //        hsplist = hsp_stream->sorted_hsplists[num_hsplists - 1 - i];
                //        if (hsplist->oid != target_oid)
                //            break;
                //
                //        batch->hsplist_array[i] = hsplist;
                //    }
                //
                //    hsp_stream->num_hsplists = num_hsplists - i;
                //    batch->num_hsplists = i;
                //
                //    return kBlastHSPStream_Success;
                // ```
                let query_indices = (0..batch.restored.original_lengths.len()).rev();
                for query_index in query_indices {
                    let input: Vec<_> = linked
                        .hsps
                        .iter()
                        .filter(|h| h.hsp.context / 6 == query_index)
                        .map(|h| h.hsp.clone())
                        .collect();
                    if input.is_empty() {
                        continue;
                    }
                    let encoded: Vec<_> = subject[s.oid]
                        .sequence
                        .iter()
                        .copied()
                        .map(LOSAT::utils::matrix::aa_char_to_ncbistdaa)
                        .collect();
                    let redone = if options.composition == 2 {
                        kappa_state
                            .as_mut()
                            .expect("Kappa state initialized")
                            .redo_list(&input, &batch.restored, &encoded)?
                    } else {
                        ordinary_traceback(
                            &input,
                            &batch.restored,
                            &s.parameters,
                            &options,
                            &encoded,
                        )?
                    };
                    let hsps: Vec<_> = redone.iter().map(|h| h.hsp.clone()).collect();
                    let mut final_list = evaluate_gapped(
                        &hsps,
                        &batch.restored,
                        &full_parameters,
                        &options,
                        encoded.len() as i32,
                        db_length,
                        if options.composition == 2 { 32.0 } else { 1.0 },
                    )?;
                    let call = final_call;
                    if !final_list.hsps.is_empty() {
                        final_call += 1;
                    }
                    dump("D_NUMERIC_FINAL", call, s.oid, &final_list);
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3675-3691
                    // ```c++
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
                    // ```
                    reap_final(&mut final_list, &options);
                    dump("D_REAP_FINAL", call, s.oid, &final_list);
                    let bits = normalize_scores(&mut final_list, &options);
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3691-3719
                    // ```c++
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
                    //                     }
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
                    // NCBI reference (598d8ae6): c++/include/algo/blast/core/gapinfo.h:44-54
                    // ```c++
                    // typedef enum EGapAlignOpType {
                    //    eGapAlignDel = 0, /**< Deletion: a gap in query */
                    //    eGapAlignDel2 = 1,/**< Frame shift deletion of two nucleotides */
                    //    eGapAlignDel1 = 2,/**< Frame shift deletion of one nucleotide */
                    //    eGapAlignSub = 3, /**< Substitution */
                    //    eGapAlignIns1 = 4,/**< Frame shift insertion of one nucleotide */
                    //    eGapAlignIns2 = 5,/**< Frame shift insertion of two nucleotides */
                    //    eGapAlignIns = 6, /**< Insertion: a gap in subject */
                    //    eGapAlignDecline = 7, /**< Non-aligned region */
                    //    eGapAlignInvalid = 8 /**< Invalid operation */
                    // } EGapAlignOpType;
                    // ```
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:321-342
                    // ```c++
                    //     for (align = *alignments;  NULL != align;  align = align->next) {
                    //         BlastHSP * new_hsp = NULL;
                    //         GapEditScript * editScript = align->context;
                    //         align->context = NULL;
                    //
                    //         status = Blast_HSPInit(align->queryStart, align->queryEnd,
                    //                                align->matchStart, align->matchEnd,
                    //                                unknown_value, unknown_value,
                    //                                align->queryIndex,
                    //                                frame, (Int2) align->frame, align->score,
                    //                                &editScript, &new_hsp);
                    //         switch (align->matrix_adjust_rule) {
                    //         case eDontAdjustMatrix:
                    //             new_hsp->comp_adjustment_method = eNoCompositionBasedStats;
                    //             break;
                    //         case eCompoScaleOldMatrix:
                    //             new_hsp->comp_adjustment_method = eCompositionBasedStats;
                    //             break;
                    //         default:
                    //             new_hsp->comp_adjustment_method = eCompositionMatrixAdjust;
                    //             break;
                    //         }
                    // ```
                    for (i, h) in final_list.hsps.iter().enumerate() {
                        let payload = &redone[h.source_index];
                        print!(
                            "D_EDIT_FINAL\t{}\t{}\t{}\t{}\t{}",
                            call,
                            s.oid,
                            i,
                            payload.composition_method,
                            payload.edit_script.len()
                        );
                        for op in &payload.edit_script {
                            let (kind, length) = match op {
                                LOSAT::common::GapEditOp::Del(n) => (0, *n),
                                LOSAT::common::GapEditOp::Sub(n) => (3, *n),
                                LOSAT::common::GapEditOp::Ins(n) => (6, *n),
                            };
                            print!("\t{kind}:{length}");
                        }
                        println!();
                    }
                    if !final_list.hsps.is_empty() {
                        println!(
                            "D_STANDARD_FINAL_COUNT\t{}\t{}\t{}",
                            call,
                            s.oid,
                            final_list.hsps.len()
                        );
                    }
                    for (i, (h, bit)) in final_list.hsps.iter().zip(bits).enumerate() {
                        let p = &h.hsp;
                        let (ident, positive) = identity_and_positive(
                            &redone[h.source_index],
                            &batch.restored,
                            &encoded,
                        );
                        println!("D_STANDARD_FINAL\t{}\t{}\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:016x}\t{:016x}\t{ident}\t{positive}",call,s.oid,p.context,p.score,p.frame,p.q_start,p.q_end,p.q_gapped_start,p.s_start,p.s_end,p.s_gapped_start,0,h.num,h.evalue.to_bits(),bit.to_bits());
                    }
                }
            }
        }
    }
    Ok(())
}
