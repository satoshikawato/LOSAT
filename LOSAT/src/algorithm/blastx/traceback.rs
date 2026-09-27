//! Ordinary BLASTX composition-0 traceback in context-local amino-acid coordinates.
use super::{
    args::ResolvedOptions,
    kappa::RedoneHsp,
    parameters::ContextParameters,
    preliminary::{purge_endpoints, sort_by_score, PreliminaryHsp},
    query_setup::PreparedQueryBatch,
};
use crate::{
    algorithm::{
        blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp},
        blastp::gapalign::{
            adjust_subject_range, blast_gapped_alignment_with_traceback_with_scratch,
            GapAlignScratch,
        },
    },
    config::{ProteinScoringSpec, ScoringMatrix},
    stats::lookup_protein_params,
};
use anyhow::{ensure, Result};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:212-226
// ```c++
// s_GetQueryStrandOffset(const BlastQueryInfo *query_info,
//                        Int4 context)
// {
//     Int4 c = context;
//
//     while (c) {
//         Int4 frame = query_info->contexts[c].frame;
//         if (frame == 0 || SIGN(frame) !=
//             SIGN(query_info->contexts[c-1].frame)) {
//             break;
//         }
//         c--;
//     }
//
//     return query_info->contexts[c].query_offset;
// ```
fn tree_hsp(h: &PreliminaryHsp, batch: &PreparedQueryBatch) -> TreeHsp {
    let strand_context = h.context / 6 * 6 + if h.context % 6 < 3 { 0 } else { 3 };
    TreeHsp {
        query_offset: h.q_start,
        query_end: h.q_end,
        subject_offset: h.s_start,
        subject_end: h.s_end,
        score: h.score,
        query_frame: h.frame as i32,
        query_length: batch.contexts[h.context].length as i32,
        query_context_offset: batch.contexts[strand_context].offset as i32,
        subject_frame_sign: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:367-403
// ```c++
//    /* set up the tree for HSP containment tests. subject_length
//       is zero only for translated subject sequences, whose maximum
//       length is bounded by the length of the first frame */
//
//    tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
//                                  0, (subject_length > 0 ? subject_length :
//                                  subject_blk->length / CODON_LENGTH) + 1);
//
//    for (index=0; index < num_initial_hsps; index++) {
//       hsp = hsp_array[index];
//       if (program_number == eBlastTypeBlastx && kIsOutOfFrame) {
//           Int4 context = hsp->context - hsp->context % CODON_LENGTH;
//           Int4 context_offset = query_info->contexts[context].query_offset;
//
//           query = query_blk->oof_sequence + CODON_LENGTH + context_offset;
//           query_nomask = query;
//           query_length = query_info->contexts[context+2].query_offset +
//               query_info->contexts[context+2].query_length - context_offset;
//       } else {
//           query = query_blk->sequence +
//               query_info->contexts[hsp->context].query_offset;
//           query_nomask = query_blk->sequence_nomask +
//               query_info->contexts[hsp->context].query_offset;
//           query_length = query_info->contexts[hsp->context].query_length;
//       }
//
//       /* preliminary RPS blast alignments have not had
//          the composition-based correction applied yet, so
//          we cannot reliably check whether an HSP is contained
//          within another */
//
//       /** @todo FIXME Traceback is always performed for rpsblast
//        * because the composition-based correction can change an
//        * HSP. It is optional for RPStblastn since no corrections
//        * are applied there. Such corrections should be added.
//        */
//       if (program_number == eBlastTypeRpsBlast ||
// ```
pub fn ordinary_traceback(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    _parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
) -> Result<Vec<RedoneHsp>> {
    ordinary_traceback_observed(input, batch, _parameters, options, subject, &mut |_, _| {})
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:402-407
// ```c++
//        */
//       if (program_number == eBlastTypeRpsBlast ||
//           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
//                                hit_options->min_diag_separation)) {
//
//          Int4 start_shift = 0;
// ```
pub(crate) fn ordinary_traceback_observed(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    _parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    observe: &mut dyn FnMut(&PreliminaryHsp, bool),
) -> Result<Vec<RedoneHsp>> {
    ensure!(
        options.gapped && options.composition == 0,
        "BLASTX ordinary traceback requires gapped composition 0"
    );
    let ka = lookup_protein_params(&ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    });
    let xdrop = ((options.gap_x_dropoff_final * std::f64::consts::LN_2 / ka.lambda) as i32)
        .max((options.gap_x_dropoff * std::f64::consts::LN_2 / ka.lambda) as i32);
    let mut tree = BlastIntervalTree::new(
        0,
        batch.sequence_start.len() as i32,
        0,
        subject.len() as i32 + 1,
    );
    let mut scratch = GapAlignScratch::new();
    let mut output = Vec::new();
    for original in input {
        let mut h = original.clone();
        let node = tree_hsp(&h, batch);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:402-407
        // ```c++
        //        */
        //       if (program_number == eBlastTypeRpsBlast ||
        //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
        //                                hit_options->min_diag_separation)) {
        //
        //          Int4 start_shift = 0;
        // ```
        let contained = tree.contains_hsp(&node, node.query_context_offset, 0);
        observe(&h, contained);
        if contained {
            continue;
        }
        let c = &batch.contexts[h.context];
        let query = &batch.sequence_start[1 + c.offset..1 + c.offset + c.length];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:437-475
        // ```c++
        //                                hsp->subject.gapped_start == 0) {
        //             Boolean retval =
        //                BlastGetOffsetsForGappedAlignment(query, subject, sbp,
        //                    hsp, &q_start, &s_start);
        //             if (!retval)
        //             {  /* Unable to find start for this HSP */
        //                hsp_array[index] = Blast_HSPFree(hsp);
        //                continue;
        //             }
        //             hsp->query.gapped_start = q_start;
        //             hsp->subject.gapped_start = s_start;
        //          } else {
        //             if(kIsOutOfFrame) {
        //                /* Code below should be investigated for possible
        //                   optimization for OOF */
        //                gap_align->subject_start = 0;
        //                gap_align->query_start = 0;
        //             } else if (program_number == eBlastTypeBlastn ||
        //                        program_number == eBlastTypeMapping) {
        //                /* Find the optimal starting offset */
        //                BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
        //             }
        //             q_start = hsp->query.gapped_start;
        //             s_start = hsp->subject.gapped_start;
        //          }
        //
        //          adjusted_s_length = subject_length;
        //          adjusted_subject = subject;
        //
        //          if (!kTranslateSubject && !kSmithWaterman) {
        //              AdjustSubjectRange(&s_start, &adjusted_s_length, q_start,
        //                                 query_length, &start_shift);
        //              adjusted_subject = subject + start_shift;
        //              /* Shift the gapped start in HSP structure, to compensate for
        //                 a shift in the other direction later. */
        //              hsp->subject.gapped_start = s_start;
        //          }
        //
        //          /* compute the cutoff score to use */
        // ```
        if h.q_gapped_start == 0 && h.s_gapped_start == 0 {
            let Some((q, s)) = get_offsets_for_gapped_alignment(query, subject, &h) else {
                continue;
            };
            h.q_gapped_start = q;
            h.s_gapped_start = s;
        }
        let (s_start, s_length, shift) = adjust_subject_range(
            h.s_gapped_start as usize,
            subject.len(),
            h.q_gapped_start as usize,
            c.length,
        );
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:479-513
        // ```c++
        //                rps_context = hsp_list->query_index * NUM_FRAMES +
        //                         BLAST_FrameToContext(hsp->subject.frame, program_number);
        //             }
        //             cutoff = hit_params->cutoffs[rps_context].cutoff_score;
        //          }
        //          else {
        //             cutoff = hit_params->cutoffs[hsp->context].cutoff_score;
        //          }
        //
        //          /* Perform the gapped extension with traceback */
        //          if (kSmithWaterman) {
        //
        //              /* with Smith-Waterman, 'hsp' is a placeholder that gives
        //                 the context and frames of local alignments. The following
        //                 will compute all of the HSPs that really belong to this
        //                 query-subject pair, then append them to hsp_list. */
        //              SmithWatermanScoreWithTraceback(program_number,
        //                          query, query_length,
        //                          adjusted_subject, adjusted_s_length,
        //                          hsp, hsp_list, score_params,
        //                          hit_params, gap_align, start_shift, cutoff);
        //              /* remove the original HSP unconditionally */
        //              gap_align->score = INT4_MIN;
        //
        //          } else if (kGreedyTraceback) {
        //              BLAST_GreedyGappedAlignment(query, adjusted_subject,
        //                  query_length, adjusted_s_length, gap_align,
        //                  score_params, q_start, s_start, FALSE, TRUE,
        //                  fence_hit);
        //          } else {
        //            BLAST_GappedAlignmentWithTraceback(program_number, query,
        //                  adjusted_subject, gap_align, score_params, q_start, s_start,
        //                  query_length, adjusted_s_length,
        //                  fence_hit);
        //              ASSERT(!(kFullTranslation && *fence_hit));
        // ```
        let Some(aligned) = blast_gapped_alignment_with_traceback_with_scratch(
            query,
            &subject[shift..shift + s_length],
            h.q_gapped_start as usize,
            s_start,
            ScoringMatrix::Blosum62,
            None,
            options.gap_open,
            options.gap_extend,
            xdrop,
            &mut scratch,
            None,
        ) else {
            continue;
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:581-605
        // ```c++
        //          if ( hsp_array[index] )
        //          {
        //            Blast_HSPUpdateWithTraceback(gap_align, hsp);
        //
        //            if (!delete_hsp && !kGreedyTraceback) {
        //                /* Calculate number of identities and check if this HSP meets the
        //                   percent identity and length criteria. */
        //                Int4 align_length = 0;
        //                Blast_HSPGetNumIdentitiesAndPositives(query_nomask,
        //                        							   adjusted_subject,
        //                        							   hsp,
        //                        							   score_options,
        //                        							   &align_length,
        //                        							   sbp);
        //
        //                delete_hsp = Blast_HSPTest(hsp, hit_options, align_length);
        //            }
        //            if (!delete_hsp) {
        //               Blast_HSPAdjustSubjectOffset(hsp, start_shift);
        //               status = BlastIntervalTreeAddHSP(hsp, tree, query_info,
        //                                          eQueryAndSubject);
        //               if (status) return status;
        //            } else {
        //               hsp_array[index] = Blast_HSPFree(hsp);
        //            }
        // ```
        h.score = aligned.score;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583,599
        // ```c
        // Blast_HSPUpdateWithTraceback(gap_align, hsp);
        // Blast_HSPAdjustSubjectOffset(hsp, start_shift);
        // ```
        h.q_start = aligned.query_start;
        h.q_end = aligned.query_stop;
        h.s_start = aligned.subject_start + i32::try_from(shift)?;
        h.s_end = aligned.subject_stop + i32::try_from(shift)?;
        let node = tree_hsp(&h, batch);
        tree.add_hsp(
            node.clone(),
            node.query_context_offset,
            IndexMethod::QueryAndSubject,
        );
        output.push(RedoneHsp {
            hsp: h,
            edit_script: aligned.edit_script,
            composition_method: 0,
            identity_matrix: None,
        });
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2467
    // ```c++
    // Blast_HSPListPurgeHSPsWithCommonEndpoints(EBlastProgramType program,
    //                                           BlastHSPList* hsp_list,
    //                                           Boolean purge)
    //
    // {
    //    BlastHSP** hsp_array;  /* hsp_array to purge. */
    //    BlastHSP* hsp;
    //    Int4 i, j, k;
    //    Int4 hsp_count;
    //    purge |= (program != eBlastTypeBlastn);
    //
    //    /* If HSP list is empty, return immediately. */
    //    if (hsp_list == NULL || hsp_list->hspcnt == 0)
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:634-693
    // ```c++
    //
    //        /* Remove any HSPs that share a starting or ending diagonal
    //           with a higher-scoring HSP. */
    //        Int4 extra_start =
    //            Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
    //
    //        /* Low level greedy algorithm ignores ambiguities, so the score
    //         * needs to be reevaluated. */
    //        if (kGreedyTraceback) {
    //           extra_start = 0;
    //        }
    //        /* Try to make use of the remaining part of the longer hsps that
    //           get purged otherwise */
    //        for (index=extra_start; index < hsp_list->hspcnt; index++) {
    //           Boolean delete_hsp = FALSE;
    //           hsp = hsp_array[index];
    //           if (!hsp) continue;
    //           query = query_blk->sequence +
    //               query_info->contexts[hsp->context].query_offset;
    //           query_nomask = query_blk->sequence_nomask +
    //               query_info->contexts[hsp->context].query_offset;
    //           query_length = query_info->contexts[hsp->context].query_length;
    //           /* the remaining part of the hsp may be extended further */
    //           delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(hsp, query,
    //                        query_length, subject, subject_length, hit_params,
    //                        score_params, sbp);
    //           if (!delete_hsp)
    //               delete_hsp = Blast_HSPTestIdentityAndLength(program_number, hsp, query_nomask,
    //                                                        subject, score_options, hit_options);
    //           if (delete_hsp)
    //               hsp_array[index] = Blast_HSPFree(hsp);
    //        }
    //        Blast_HSPListPurgeNullHSPs(hsp_list);
    //        if(program_number == eBlastTypeBlastn) {
    //     	   Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    //        }
    //
    //        /* Sort HSPs by score again, as the scores might have changed. */
    //        Blast_HSPListSortByScore(hsp_list);
    //
    //        /* Remove any HSPs that are contained within other HSPs.
    //           Since the list is sorted by score already, any HSP
    //           contained by a previous HSP is guaranteed to have a
    //           lower score, and may be purged. */
    //        Blast_IntervalTreeReset(tree);
    //        for (index = 0; index < hsp_list->hspcnt; index++) {
    //            hsp = hsp_array[index];
    //
    //            if (BlastIntervalTreeContainsHSP(tree, hsp, query_info,
    //                                      hit_options->min_diag_separation)) {
    //                hsp_array[index] = Blast_HSPFree(hsp);
    //            }
    //            else {
    //                status = BlastIntervalTreeAddHSP(hsp, tree, query_info,
    //                                        eQueryAndSubject);
    //                if (status)
    //                   return status;
    //            }
    //        }
    //    }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:672-693
    // ```c++
    //        Blast_HSPListSortByScore(hsp_list);
    //
    //        /* Remove any HSPs that are contained within other HSPs.
    //           Since the list is sorted by score already, any HSP
    //           contained by a previous HSP is guaranteed to have a
    //           lower score, and may be purged. */
    //        Blast_IntervalTreeReset(tree);
    //        for (index = 0; index < hsp_list->hspcnt; index++) {
    //            hsp = hsp_array[index];
    //
    //            if (BlastIntervalTreeContainsHSP(tree, hsp, query_info,
    //                                      hit_options->min_diag_separation)) {
    //                hsp_array[index] = Blast_HSPFree(hsp);
    //            }
    //            else {
    //                status = BlastIntervalTreeAddHSP(hsp, tree, query_info,
    //                                        eQueryAndSubject);
    //                if (status)
    //                   return status;
    //            }
    //        }
    //    }
    // ```
    // Preserve each edit script with its original HSP while the source comparator
    // and endpoint purge rearrange the HSP values. Duplicate keys retain input order.
    let key = |h: &PreliminaryHsp| {
        (
            h.context,
            h.frame,
            h.score,
            h.q_start,
            h.q_end,
            h.q_gapped_start,
            h.s_start,
            h.s_end,
            h.s_gapped_start,
        )
    };
    let mut payloads = std::collections::HashMap::<_, std::collections::VecDeque<_>>::new();
    for p in &output {
        payloads
            .entry(key(&p.hsp))
            .or_default()
            .push_back(p.clone());
    }
    let mut hsps: Vec<_> = output.iter().map(|h| h.hsp.clone()).collect();
    purge_endpoints(&mut hsps);
    sort_by_score(&mut hsps);
    tree.reset();
    let mut retained = Vec::new();
    for h in hsps {
        let node = tree_hsp(&h, batch);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:677-690
        // ```c++
        //           lower score, and may be purged. */
        //        Blast_IntervalTreeReset(tree);
        //        for (index = 0; index < hsp_list->hspcnt; index++) {
        //            hsp = hsp_array[index];
        //
        //            if (BlastIntervalTreeContainsHSP(tree, hsp, query_info,
        //                                      hit_options->min_diag_separation)) {
        //                hsp_array[index] = Blast_HSPFree(hsp);
        //            }
        //            else {
        //                status = BlastIntervalTreeAddHSP(hsp, tree, query_info,
        //                                        eQueryAndSubject);
        //                if (status)
        //                   return status;
        // ```
        let contained = tree.contains_hsp(&node, node.query_context_offset, 0);
        observe(&h, contained);
        if contained {
            continue;
        }
        tree.add_hsp(
            node.clone(),
            node.query_context_offset,
            IndexMethod::QueryAndSubject,
        );
        retained.push(
            payloads
                .get_mut(&key(&h))
                .and_then(|queue| queue.pop_front())
                .expect("traceback payload remains attached"),
        );
    }
    Ok(retained)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3248-3330
// ```c++
// BlastGetOffsetsForGappedAlignment (const Uint1* query, const Uint1* subject,
//    const BlastScoreBlk* sbp, BlastHSP* hsp, Int4* q_retval, Int4* s_retval)
// {
//     Int4 index1, max_offset, score, max_score, hsp_end;
//     const Uint1* query_var,* subject_var;
//     Boolean positionBased = (sbp->psi_matrix != NULL);
//     Int4 q_length = hsp->query.end - hsp->query.offset;
//     Int4 s_length = hsp->subject.end - hsp->subject.offset;
//     int q_start = hsp->query.offset;
//     int s_start = hsp->subject.offset;
//
//     if (q_length <= HSP_MAX_WINDOW) {
//         *q_retval = q_start + q_length/2;
//         *s_retval = s_start + q_length/2;
//         return TRUE;
//     }
//
//     hsp_end = q_start + HSP_MAX_WINDOW;
//     query_var = query + q_start;
//     subject_var = subject + s_start;
//     score=0;
//     for (index1=q_start; index1<hsp_end; index1++) {
//         if (!(positionBased))
//             score += sbp->matrix->data[*query_var][*subject_var];
//         else
//             score += sbp->psi_matrix->pssm->data[index1][*subject_var];
//         query_var++; subject_var++;
//     }
//     max_score = score;
//     max_offset = hsp_end - 1;
//     hsp_end = q_start + MIN(q_length, s_length);
//     for (index1=q_start + HSP_MAX_WINDOW; index1<hsp_end; index1++) {
//         if (!(positionBased)) {
//             score -= sbp->matrix->data[*(query_var-HSP_MAX_WINDOW)][*(subject_var-HSP_MAX_WINDOW)];
//             score += sbp->matrix->data[*query_var][*subject_var];
//         } else {
//             score -= sbp->psi_matrix->pssm->data[index1-HSP_MAX_WINDOW][*(subject_var-HSP_MAX_WINDOW)];
//             score += sbp->psi_matrix->pssm->data[index1][*subject_var];
//         }
//         if (score > max_score) {
//             max_score = score;
//             max_offset = index1;
//         }
//         query_var++; subject_var++;
//     }
//
//     if (max_score > 0)
//     {
//         *q_retval = max_offset;
//         *s_retval = (max_offset - q_start) + s_start;
//         return TRUE;
//     }
//     else  /* Test the window around the ends of the HSP. */
//     {
//         score=0;
//         query_var = query + q_start + q_length - HSP_MAX_WINDOW;
//         subject_var = subject + s_start + s_length - HSP_MAX_WINDOW;
//         for (index1=hsp->query.end-HSP_MAX_WINDOW; index1<hsp->query.end; index1++) {
//             if (!(positionBased))
//                 score += sbp->matrix->data[*query_var][*subject_var];
//             else
//                 score += sbp->psi_matrix->pssm->data[index1][*subject_var];
//             query_var++; subject_var++;
//         }
//         if (score > 0)
//         {
//             *q_retval = hsp->query.end - HSP_MAX_WINDOW/2;
//             *s_retval = hsp->subject.end - HSP_MAX_WINDOW/2;
//             return TRUE;
//         }
//     }
//     return FALSE;
// }
//
// void
// BlastGetStartForGappedAlignmentNucl (const Uint1* query, const Uint1* subject,
//    BlastHSP* hsp)
// {
//     /* We will stop when the identity count reaches to this number */
//     int hspMaxIdentRun = 10;
//     const Uint1 *q, *s;
//     Int4 index, max_offset, score, max_score, q_start, s_start, q_len;
//     Boolean match, prev_match;
// ```
fn get_offsets_for_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    h: &PreliminaryHsp,
) -> Option<(i32, i32)> {
    let q_start = h.q_start as usize;
    let s_start = h.s_start as usize;
    let q_length = (h.q_end - h.q_start) as usize;
    let s_length = (h.s_end - h.s_start) as usize;
    if q_length <= 11 {
        return Some((
            h.q_start + (q_length / 2) as i32,
            h.s_start + (q_length / 2) as i32,
        ));
    }
    let score_at = |q: usize, s: usize| crate::utils::matrix::blosum62_score(query[q], subject[s]);
    let mut score: i32 = (0..11).map(|i| score_at(q_start + i, s_start + i)).sum();
    let mut max_score = score;
    let mut max_offset = q_start + 10;
    for i in 11..q_length.min(s_length) {
        score -= score_at(q_start + i - 11, s_start + i - 11);
        score += score_at(q_start + i, s_start + i);
        if score > max_score {
            max_score = score;
            max_offset = q_start + i;
        }
    }
    if max_score > 0 {
        return Some((max_offset as i32, (max_offset - q_start + s_start) as i32));
    }
    let end_score: i32 = (0..11)
        .map(|i| score_at(q_start + q_length - 11 + i, s_start + s_length - 11 + i))
        .sum();
    if end_score > 0 {
        Some((h.q_end - 5, h.s_end - 5))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3248-3330
    // ```c++
    // BlastGetOffsetsForGappedAlignment (const Uint1* query, const Uint1* subject,
    //    const BlastScoreBlk* sbp, BlastHSP* hsp, Int4* q_retval, Int4* s_retval)
    // {
    //     Int4 index1, max_offset, score, max_score, hsp_end;
    //     const Uint1* query_var,* subject_var;
    //     Boolean positionBased = (sbp->psi_matrix != NULL);
    //     Int4 q_length = hsp->query.end - hsp->query.offset;
    //     Int4 s_length = hsp->subject.end - hsp->subject.offset;
    //     int q_start = hsp->query.offset;
    //     int s_start = hsp->subject.offset;
    //
    //     if (q_length <= HSP_MAX_WINDOW) {
    //         *q_retval = q_start + q_length/2;
    //         *s_retval = s_start + q_length/2;
    //         return TRUE;
    //     }
    //
    //     hsp_end = q_start + HSP_MAX_WINDOW;
    //     query_var = query + q_start;
    //     subject_var = subject + s_start;
    //     score=0;
    //     for (index1=q_start; index1<hsp_end; index1++) {
    //         if (!(positionBased))
    //             score += sbp->matrix->data[*query_var][*subject_var];
    //         else
    //             score += sbp->psi_matrix->pssm->data[index1][*subject_var];
    //         query_var++; subject_var++;
    //     }
    //     max_score = score;
    //     max_offset = hsp_end - 1;
    //     hsp_end = q_start + MIN(q_length, s_length);
    //     for (index1=q_start + HSP_MAX_WINDOW; index1<hsp_end; index1++) {
    //         if (!(positionBased)) {
    //             score -= sbp->matrix->data[*(query_var-HSP_MAX_WINDOW)][*(subject_var-HSP_MAX_WINDOW)];
    //             score += sbp->matrix->data[*query_var][*subject_var];
    //         } else {
    //             score -= sbp->psi_matrix->pssm->data[index1-HSP_MAX_WINDOW][*(subject_var-HSP_MAX_WINDOW)];
    //             score += sbp->psi_matrix->pssm->data[index1][*subject_var];
    //         }
    //         if (score > max_score) {
    //             max_score = score;
    //             max_offset = index1;
    //         }
    //         query_var++; subject_var++;
    //     }
    //
    //     if (max_score > 0)
    //     {
    //         *q_retval = max_offset;
    //         *s_retval = (max_offset - q_start) + s_start;
    //         return TRUE;
    //     }
    //     else  /* Test the window around the ends of the HSP. */
    //     {
    //         score=0;
    //         query_var = query + q_start + q_length - HSP_MAX_WINDOW;
    //         subject_var = subject + s_start + s_length - HSP_MAX_WINDOW;
    //         for (index1=hsp->query.end-HSP_MAX_WINDOW; index1<hsp->query.end; index1++) {
    //             if (!(positionBased))
    //                 score += sbp->matrix->data[*query_var][*subject_var];
    //             else
    //                 score += sbp->psi_matrix->pssm->data[index1][*subject_var];
    //             query_var++; subject_var++;
    //         }
    //         if (score > 0)
    //         {
    //             *q_retval = hsp->query.end - HSP_MAX_WINDOW/2;
    //             *s_retval = hsp->subject.end - HSP_MAX_WINDOW/2;
    //             return TRUE;
    //         }
    //     }
    //     return FALSE;
    // }
    //
    // void
    // BlastGetStartForGappedAlignmentNucl (const Uint1* query, const Uint1* subject,
    //    BlastHSP* hsp)
    // {
    //     /* We will stop when the identity count reaches to this number */
    //     int hspMaxIdentRun = 10;
    //     const Uint1 *q, *s;
    //     Int4 index, max_offset, score, max_score, q_start, s_start, q_len;
    //     Boolean match, prev_match;
    // ```
    #[test]
    fn offsets_match_pinned_ncbi_boundary_oracle() {
        let expected =
            include_str!("../../../../docs/evidence/losatx_stage_d/offsets_expected.tsv");
        let mut count = 0;
        for line in expected.lines() {
            let f: Vec<i32> = line.split('\t').map(|v| v.parse().unwrap()).collect();
            let mode = f[0];
            let n = f[1] as usize;
            let q = vec![1u8; 100];
            let mut s = vec![if mode == 0 { 1 } else { 21 }; 100];
            let end = 7 + n + if mode == 3 { 12 } else { 0 };
            if mode == 2 {
                s[7 + n / 2..7 + n].fill(1);
            }
            if mode == 3 && n > 11 {
                s[end - 11..end].fill(1);
            }
            let h = PreliminaryHsp {
                context: 0,
                frame: 1,
                score: 0,
                q_start: 3,
                q_end: 3 + n as i32,
                q_gapped_start: 0,
                s_start: 7,
                s_end: end as i32,
                s_gapped_start: 0,
            };
            let result = get_offsets_for_gapped_alignment(&q, &s, &h);
            assert_eq!(
                result,
                if f[2] != 0 { Some((f[3], f[4])) } else { None },
                "mode {mode} length {n}"
            );
            count += 1;
        }
        assert_eq!(count, 204);
    }
}
