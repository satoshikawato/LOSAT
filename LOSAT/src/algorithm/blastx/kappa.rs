//! BLASTX six-context Kappa redo adapter. No BLASTP searches or external runtime calls.
use super::{
    args::ResolvedOptions, linking::link_parameters, parameters::ContextParameters,
    preliminary::PreliminaryHsp, query_setup::PreparedQueryBatch,
};
use crate::{
    algorithm::blastp::gapalign::{
        blast_gapped_alignment_with_traceback_with_scratch, GapAlignScratch,
    },
    config::{ProteinScoringSpec, ScoringMatrix},
    core::{
        blast_stat::compute_blosum62_ideal_karlin_params,
        composition_adjustment::{
            adjust_scores::{
                build_matrix_info, AdjustedProteinMatrix, BlastAminoAcidComposition,
                BlastCompositionWorkspace,
            },
            redo_alignment::*,
        },
    },
    stats::lookup_protein_params,
    utils::{
        matrix::BLASTAA_SIZE,
        seg::{SegMasker, SegParams},
    },
};
use anyhow::{bail, ensure, Context, Result};
use std::{cell::Cell, ffi::c_void, ptr::NonNull};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2309-2336
// ```c++
// s_GetQueryInfo(Uint1 * query_data, const BlastQueryInfo * blast_query_info, Boolean skip)
// {
//     int i;                   /* loop index */
//     BlastCompo_QueryInfo *
//         compo_query_info;    /* the new array */
//     int num_queries;         /* the number of queries/elements in
//                                 compo_query_info */
//
//     num_queries = blast_query_info->last_context + 1;
//     compo_query_info = calloc(num_queries, sizeof(BlastCompo_QueryInfo));
//     if (compo_query_info != NULL) {
//         for (i = 0;  i < num_queries;  i++) {
//             BlastCompo_QueryInfo * query_info = &compo_query_info[i];
//             const BlastContextInfo * query_context = &blast_query_info->contexts[i];
//
//             query_info->eff_search_space =
//                 (double) query_context->eff_searchsp;
//             query_info->origin = query_context->query_offset;
//             query_info->seq.data = &query_data[query_info->origin];
//             query_info->seq.length = query_context->query_length;
//             query_info->words = NULL;
//
//             s_CreateWordArray(query_info->seq.data, query_info->seq.length,
//                               &query_info->words);
//             if (! skip) {
//                 Blast_ReadAaComposition(&query_info->composition, BLASTAA_SIZE,
//                                         query_info->seq.data,
//                                         query_info->seq.length);
// ```
pub fn query_infos(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
) -> Vec<BlastCompoQueryInfo> {
    batch
        .contexts
        .iter()
        .zip(parameters)
        .map(|(c, p)| {
            let query = &batch.sequence_start[1 + c.offset..1 + c.offset + c.length];
            BlastCompoQueryInfo {
                origin: c.offset as i32,
                seq: BlastCompoSequenceData::from_ncbistdaa(query),
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3244-3248
                // ```c++
                //         query_info_tld[i] = s_GetQueryInfo(
                //                 queryBlk->sequence,
                //                 queryInfo,
                //                 (program_number == eBlastTypeBlastx)
                //         );
                // ```
                composition: BlastAminoAcidComposition::empty(),
                eff_search_space: p.search_space as f64,
                words: Some(build_query_word_hashes(query)),
            }
        })
        .collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3100-3109
// ```c++
//     if (compo_adjust_mode != eNoCompositionBasedStats) {
//         if((0 == strcmp(scoringParams->options->matrix, "BLOSUM62_20"))) {
//             localScalingFactor = SCALING_FACTOR / 10;
//         } else {
//             localScalingFactor = SCALING_FACTOR;
//         }
//     } else {
//         localScalingFactor = 1.0;
//     }
//     s_RescaleSearch(sbp, scoringParams, numContexts, localScalingFactor);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2438-2479
// ```c++
//     Int4 index;
//     for (index = queryInfo->first_context;
//                      index <= queryInfo->last_context; ++index) {
//
//         if ((queryInfo->contexts[index].is_valid)) {
//     		near_identical_cutoff =
//        		 (near_identical_cutoff_bits * NCBIMATH_LN2)
//         		/ context->sbp->kbp_gap[index]->Lambda;
// 		break;
// 	}
//     }
//
//     if (do_link_hsps) {
//         ASSERT(hitParams->link_hsp_params != NULL);
//         cutoff_s =
//             (int) (hitParams->cutoff_score_min * context->localScalingFactor);
//     } else {
//         /* There is no cutoff score; we consider e-values instead */
//         cutoff_s = 1;
//     }
//     cutoff_e = hitParams->options->expect_value;
//     rows = positionBased ? queryInfo->max_length : BLASTAA_SIZE;
//     scaledMatrixInfo = Blast_MatrixInfoNew(rows, BLASTAA_SIZE, positionBased);
//     status = s_MatrixInfoInit(scaledMatrixInfo, queryBlk, context->sbp,
//                               context->localScalingFactor,
//                               context->scoringParams->options->matrix);
//     if (status != 0) {
//         return NULL;
//     }
//     gapping_params = s_GappingParamsNew(context, extendParams,
//                                         queryInfo->last_context + 1);
//     if (gapping_params == NULL) {
//         return NULL;
//     } else {
//         return
//             Blast_RedoAlignParamsNew(&scaledMatrixInfo, &gapping_params,
//                                      compo_adjust_mode, positionBased,
//                                      query_is_translated,
//                                      subject_is_translated,
//                                      queryInfo->max_length, cutoff_s, cutoff_e,
//                                      do_link_hsps, &redo_align_callbacks,
//                                      near_identical_cutoff);
// ```
pub fn redo_params(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
) -> Result<BlastRedoAlignParams> {
    ensure!(
        options.composition == 2 && options.gapped,
        "BLASTX Kappa requires gapped composition 2"
    );
    let scale = 32.0;
    let ka = lookup_protein_params(&ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    });
    let scaled_lambda = ka.lambda / scale;
    let ideal = compute_blosum62_ideal_karlin_params().map_err(anyhow::Error::msg)?;
    let do_link = link_parameters(options).is_some();
    let cutoff = parameters
        .iter()
        .filter(|p| p.valid)
        .map(|p| p.hit_cutoff)
        .min()
        .context("BLASTX Kappa requires a valid context")?;
    let extension_final = (options.gap_x_dropoff_final * std::f64::consts::LN_2 / ka.lambda)
        .max(options.gap_x_dropoff * std::f64::consts::LN_2 / ka.lambda)
        as i32;
    Ok(BlastRedoAlignParams {
        matrix_info: build_matrix_info(ScoringMatrix::Blosum62, ideal.lambda / scale)?,
        gapping_params: BlastCompoGappingParams {
            gap_open: options.gap_open * 32,
            gap_extend: options.gap_extend * 32,
            decline_align: 0,
            x_dropoff: ((options.gap_x_dropoff_final * std::f64::consts::LN_2 / scaled_lambda)
                as i32)
                .max(extension_final),
            context: Cell::new(None),
        },
        compo_adjust_mode: BlastCompoAdjustMode::CompositionMatrixAdjust,
        alphsize: BLASTAA_SIZE as i32,
        composition_test_index: 0,
        unified_p: false,
        log_k: 0.0,
        score_divisor: scale,
        restricted_alignment: false,
        smith_waterman: false,
        is_same_adjustment: false,
        near_identical_cutoff: 1.74 * std::f64::consts::LN_2 / scaled_lambda,
        position_based: false,
        re_matrix_adjustment_pseudocounts: 20,
        ccat_query_length: batch
            .contexts
            .iter()
            .map(|c| c.length as i32)
            .max()
            .unwrap_or(0),
        query_is_translated: true,
        subject_is_translated: false,
        cutoff_score: if do_link {
            (cutoff as f64 * scale) as i32
        } else {
            1
        },
        cutoff_evalue: options.evalue,
        do_link_hsps: do_link,
    })
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1896-1913
// ```c++
//  * @sa redo_one_alignment_type
//  */
// static BlastCompo_Alignment *
// s_RedoOneAlignment(BlastCompo_Alignment * in_align,
//                    EMatrixAdjustRule matrix_adjust_rule,
//                    BlastCompo_SequenceData * query_data,
//                    BlastCompo_SequenceRange * query_range,
//                    int ccat_query_length,
//                    BlastCompo_SequenceData * subject_data,
//                    BlastCompo_SequenceRange * subject_range,
//                    int full_subject_length,
//                    BlastCompo_GappingParams * gapping_params)
// {
//     int status;                /* return code */
//     Int4 q_start, s_start;     /* starting point in query and subject */
//     /* BLAST-specific parameters needed to compute a gapped alignment */
//     BlastKappa_GappingParamsContext * context = gapping_params->context;
//     /* Auxiliary structure for computing gapped alignments */
// ```
struct RedoContext<'a> {
    input: &'a [PreliminaryHsp],
    scratch: NonNull<GapAlignScratch>,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1677-1709
// ```c++
//                    const Uint8* query_words,
//                    const BlastCompo_Alignment *align,
//                    const Boolean shouldTestIdentical,
//                    const ECompoAdjustModes compo_adjust_mode,
//                    const Boolean isSmithWaterman,
//                    Boolean* subject_maybe_biased)
// {
//     Int4 idx;
//     BlastKappa_SequenceInfo * seq_info = self->local_data;
//     Uint1 *origData = query->data + q_range->begin;
//     /* Copy the query sequence (necessary for SEG filtering.) */
//     queryData->length = q_range->end - q_range->begin;
//     queryData->buffer = calloc((queryData->length + 2), sizeof(Uint1));
//     queryData->data   = queryData->buffer + 1;
//
//     for (idx = 0;  idx < queryData->length;  idx++) {
//         /* Copy the sequence data, replacing occurrences of amino acid
//          * number 24 (Selenocysteine) with number 3 (Cysteine). */
//         queryData->data[idx] = (origData[idx] != 24) ? origData[idx] : 3;
//     }
//     if (seq_info && seq_info->prog_number ==  eBlastTypeTblastn) {
//         /* The sequence must be translated. */
//         return s_SequenceGetTranslatedRange(self, s_range, seqData,
//                                             q_range, queryData, query_words,
//                                             align, shouldTestIdentical,
//                                             compo_adjust_mode, isSmithWaterman,
//                                             subject_maybe_biased);
//     } else {
//         return s_SequenceGetProteinRange(self, s_range, seqData,
//                                          q_range, queryData, query_words,
//                                          align, shouldTestIdentical,
//                                          compo_adjust_mode, isSmithWaterman,
//                                          subject_maybe_biased);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1626-1646
// ```c++
//         if (compo_adjust_mode
//             && (!subject_maybe_biased || *subject_maybe_biased)) {
//
//             if ( (!shouldTestIdentical)
//                  || (shouldTestIdentical
//                      && (!s_TestNearIdentical(seqData, 0, queryData,
//                                               q_range->begin, query_words,
//                                               align)))) {
//
//                 status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
//                                              subject_maybe_biased);
//             }
//         }
//     }
//     /* Fit the data to the range. */
//     seqData ->data    = &seqData->data[range->begin - 1];
//     *seqData->data++  = '\0';
//     seqData ->length  = range->end - range->begin;
//
//     if (status != 0) {
//         free(seqData->buffer);
// ```
fn get_range<'a>(
    matching: &BlastCompoMatchingSequence<'a>,
    subject_range: &BlastCompoSequenceRange,
    orig_query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    words: Option<&[u64]>,
    align: &BlastCompoAlignment,
    should_test: bool,
    biased: bool,
    params: &BlastRedoAlignParams,
) -> Result<BlastRedoRangeResult> {
    let query =
        BlastCompoSequenceData::copy_query_range_with_selenocysteine_fix(orig_query, query_range);
    let mut subject = BlastCompoSequenceData::from_ncbistdaa(matching.data);
    let mut biased = biased;
    if params.uses_composition_based_stats()
        && biased
        && (!should_test
            || !test_near_identical(&subject, 0, &query, query_range.begin, words, align))
    {
        let intervals =
            // BLASTX keeps LOSAT's former SEG until SX (plan DW-10).
            SegMasker::with_params(&SegParams::new(10, 1.8, 2.1))
                .keeping_all_left_segments()
                .mask_sequence(subject.data());
        biased = !intervals.is_empty();
        for interval in intervals {
            for residue in &mut subject.buffer[1 + interval.start..1 + interval.end] {
                *residue = 21;
            }
        }
    }
    subject.fit_protein_range_in_place(subject_range);
    Ok(BlastRedoRangeResult {
        query,
        subject,
        subject_maybe_biased: biased,
    })
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1896-1957
// ```c++
//  * @sa redo_one_alignment_type
//  */
// static BlastCompo_Alignment *
// s_RedoOneAlignment(BlastCompo_Alignment * in_align,
//                    EMatrixAdjustRule matrix_adjust_rule,
//                    BlastCompo_SequenceData * query_data,
//                    BlastCompo_SequenceRange * query_range,
//                    int ccat_query_length,
//                    BlastCompo_SequenceData * subject_data,
//                    BlastCompo_SequenceRange * subject_range,
//                    int full_subject_length,
//                    BlastCompo_GappingParams * gapping_params)
// {
//     int status;                /* return code */
//     Int4 q_start, s_start;     /* starting point in query and subject */
//     /* BLAST-specific parameters needed to compute a gapped alignment */
//     BlastKappa_GappingParamsContext * context = gapping_params->context;
//     /* Auxiliary structure for computing gapped alignments */
//     BlastGapAlignStruct* gapAlign = context->gap_align;
//     /* The preliminary gapped HSP that were are recomputing */
//     BlastHSP * hsp = in_align->context;
//     Boolean fence_hit = FALSE;
//
//     /* suppress unused parameter warnings; this is a callback
//        function, so these parameter cannot be deleted */
//     (void) ccat_query_length;
//     (void) full_subject_length;
//
//     /* Use the starting point supplied by the HSP. */
//     q_start = hsp->query.gapped_start - query_range->begin;
//     s_start = hsp->subject.gapped_start - subject_range->begin;
//
//     gapAlign->gap_x_dropoff = gapping_params->x_dropoff;
//
//     /*
//      * Previously, last argument was NULL which could cause problems for
//      * tblastn.
//      */
//     status =
//         BLAST_GappedAlignmentWithTraceback(context->prog_number,
//                                            query_data->data,
//                                            subject_data->data, gapAlign,
//                                            context->scoringParams,
//                                            q_start, s_start,
//                                            query_data->length,
//                                            subject_data->length,
//                                            &fence_hit);
//     if (status == 0) {
//         return s_NewAlignmentFromGapAlign(gapAlign, &gapAlign->edit_script,
//                                           query_range, subject_range,
//                                           matrix_adjust_rule);
//     } else {
//         return NULL;
//     }
// }
//
//
// /**
//  * A BlastKappa_SavedParameters holds the value of certain search
//  * parameters on entry to RedoAlignmentCore.  These values are
//  * restored on exit.
//  */
// ```
fn redo_alignment(
    incoming: &BlastCompoAlignment,
    rule: EMatrixAdjustRule,
    adjusted: Option<&AdjustedProteinMatrix>,
    query: &BlastCompoSequenceData,
    qrange: &BlastCompoSequenceRange,
    _full_query: i32,
    subject: &BlastCompoSequenceData,
    srange: &BlastCompoSequenceRange,
    _full_subject: i32,
    params: &BlastRedoAlignParams,
) -> Result<Option<Box<BlastCompoAlignment>>> {
    let Some(BlastCompoAlignmentContext::PreliminaryHspIndex(index)) = incoming.context.as_ref()
    else {
        bail!("BLASTX Kappa incoming HSP index missing")
    };
    let ptr = params
        .gapping_params
        .context
        .get()
        .context("BLASTX Kappa context missing")?;
    // Safety: redo_context installs this stack-owned state only for the synchronous redo call and restores the pointer afterward.
    let context = unsafe { &*ptr.cast::<RedoContext<'_>>().as_ptr() };
    let h = &context.input[*index];
    // Safety: the calling function retains exclusive ownership of scratch throughout the synchronous callback.
    let scratch = unsafe { &mut *context.scratch.as_ptr() };
    let mut fence = false;
    let Some(mut result) = blast_gapped_alignment_with_traceback_with_scratch(
        query.data(),
        subject.data(),
        usize::try_from(h.q_gapped_start - qrange.begin)?,
        usize::try_from(h.s_gapped_start - srange.begin)?,
        ScoringMatrix::Blosum62,
        adjusted,
        params.gapping_params.gap_open,
        params.gapping_params.gap_extend,
        params.gapping_params.x_dropoff,
        scratch,
        Some(&mut fence),
    ) else {
        return Ok(None);
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1760-1770
    // ```c
    // queryStart = gap_align->query_start + query_range->begin;
    // queryEnd = gap_align->query_stop + query_range->begin;
    // matchStart = gap_align->subject_start + subject_range->begin;
    // matchEnd = gap_align->subject_stop + subject_range->begin;
    // obj = BlastCompo_AlignmentNew(gap_align->score, matrix_adjust_rule,
    //     queryStart, queryEnd, queryIndex, matchStart, matchEnd, frame, *edit_script);
    // ```
    Ok(Some(blast_compo_alignment_new(
        result.score,
        rule,
        result.query_start + qrange.begin,
        result.query_stop + qrange.begin,
        qrange.context,
        result.subject_start + srange.begin,
        result.subject_stop + srange.begin,
        srange.context,
        Some(BlastCompoAlignmentContext::EditScript(std::mem::take(
            &mut result.edit_script,
        ))),
    )))
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2389-2394
// ```c++
// /** Callbacks used by the Blast_RedoOneMatch* routines */
// static const Blast_RedoAlignCallbacks
// redo_align_callbacks = {
//     s_CalcLambda, s_SequenceGetRange, s_RedoOneAlignment,
//     s_NewAlignmentUsingXdrop, s_FreeEditScript
// };
// ```
const CALLBACKS: BlastRedoAlignCallbacks = BlastRedoAlignCallbacks {
    calc_lambda: Some(redo_calc_lambda),
    get_range,
    redo_one_alignment: redo_alignment,
    new_xdrop_align: None,
    free_align_traceback: None,
};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:770-803
// ```c++
// s_ResultHspToDistinctAlign(BlastCompo_Alignment **self,
//                            int *numAligns,
//                            BlastHSP * hsp_array[], Int4 hspcnt,
//                            int init_context,
//                            const BlastQueryInfo* queryInfo,
//                            double localScalingFactor)
// {
//     BlastCompo_Alignment * tail[6];        /* last element in aligns */
//     int hsp_index;                             /* loop index */
//     int frame_index;
//
//     for (frame_index = 0; frame_index < 6; frame_index++) {
//         tail[frame_index] = NULL;
//         numAligns[frame_index] = 0;
//     }
//
//     for (hsp_index = 0;  hsp_index < hspcnt;  hsp_index++) {
//         BlastHSP * hsp = hsp_array[hsp_index]; /* current HSP */
//         BlastCompo_Alignment * new_align;      /* newly-created alignment */
//         frame_index = hsp->context - init_context;
//         ASSERT(frame_index < 6 && frame_index >= 0);
//         /* Incoming alignments will have coordinates of the query
//            portion relative to a particular query context; they must
//            be shifted for used in the composition_adjustment library.
//         */
//         new_align =
//             BlastCompo_AlignmentNew((int) (hsp->score * localScalingFactor),
//                                     eDontAdjustMatrix,
//                                     hsp->query.offset, hsp->query.end, hsp->context,
//                                     hsp->subject.offset, hsp->subject.end,
//                                     hsp->subject.frame, hsp);
//         if (new_align == NULL) /* out of memory */
//             return -1;
//         if (tail[frame_index] == NULL) { /* if the list aligns is empty; */
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3590-3641
// ```c++
//                 hsp_list = Blast_HSPListNew(0);
//                 for (frame_index = 0;
//                         frame_index < numFrames;
//                         frame_index++, context_index++) {
//                     incoming_aligns = incoming_align_set[frame_index];
//                     if (!incoming_aligns) {
//                         continue;
//                     }
//                     /*
//                      * All alignments in thisMatch should be to the same query
//                      */
//                     kbp = sbp->kbp_gap[context_index];
//                     if (smithWaterman) {
//                         *pStatusCode =
//                                 Blast_RedoOneMatchSmithWaterman(
//                                         alignments,
//                                         redo_align_params,
//                                         incoming_aligns,
//                                         numAligns[frame_index],
//                                         kbp->Lambda,
//                                         kbp->logK,
//                                         &matchingSeq,
//                                         query_info,
//                                         numQueries,
//                                         matrix,
//                                         BLASTAA_SIZE,
//                                         NRrecord,
//                                         forbidden,
//                                         redoneMatches,
//                                         &pvalueForThisPair,
//                                         compositionTestIndex,
//                                         &LambdaRatio
//                                 );
//                     } else {
//                         *pStatusCode =
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
// ```
pub(crate) fn redo_context(
    input: &[PreliminaryHsp],
    context_index: usize,
    infos: &[BlastCompoQueryInfo],
    subject: &[u8],
    params: &BlastRedoAlignParams,
    lambda: f64,
    scratch: &mut GapAlignScratch,
    workspace: &mut BlastCompositionWorkspace,
    matrix_state: &mut Option<AdjustedProteinMatrix>,
) -> Result<BlastRedoOneMatchResult> {
    redo_context_observed(
        input,
        context_index,
        infos,
        subject,
        params,
        lambda,
        scratch,
        workspace,
        matrix_state,
        &mut |_| {},
    )
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3641
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
// ```
pub(crate) fn redo_context_observed(
    input: &[PreliminaryHsp],
    context_index: usize,
    infos: &[BlastCompoQueryInfo],
    subject: &[u8],
    params: &BlastRedoAlignParams,
    lambda: f64,
    scratch: &mut GapAlignScratch,
    workspace: &mut BlastCompositionWorkspace,
    matrix_state: &mut Option<AdjustedProteinMatrix>,
    trace: &mut dyn for<'a> FnMut(BlastRedoTraceEvent<'a>),
) -> Result<BlastRedoOneMatchResult> {
    let mut incoming = None;
    let mut tail = &mut incoming;
    for (index, h) in input
        .iter()
        .enumerate()
        .filter(|(_, h)| h.context == context_index)
    {
        let node = blast_compo_alignment_new(
            (h.score as f64 * params.score_divisor) as i32,
            EMatrixAdjustRule::DontAdjustMatrix,
            h.q_start,
            h.q_end,
            h.context as i32,
            h.s_start,
            h.s_end,
            0,
            Some(BlastCompoAlignmentContext::PreliminaryHspIndex(index)),
        );
        *tail = Some(node);
        tail = &mut tail.as_mut().expect("incoming inserted").next;
    }
    let matching = BlastCompoMatchingSequence::new(0, subject);
    let state = RedoContext {
        input,
        scratch: NonNull::from(scratch),
    };
    let old = params
        .gapping_params
        .context
        .replace(Some(NonNull::from(&state).cast::<c_void>()));
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3641
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
    // ```
    trace(BlastRedoTraceEvent::MatchStart {
        context: context_index,
        count: input.iter().filter(|h| h.context == context_index).count(),
        matrix: matrix_state.as_ref(),
    });
    let result = blast_redo_one_match_with_workspace_queries_and_matrix_observed(
        &incoming,
        params,
        &matching,
        infos,
        lambda,
        &CALLBACKS,
        workspace,
        matrix_state,
        trace,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3641-3647
    // ```c++
    //                                 );
    //                     }
    //
    //                     if (*pStatusCode != 0) {
    //                         goto match_loop_cleanup;
    //                     }
    //
    // ```
    trace(BlastRedoTraceEvent::MatchEnd {
        context: context_index,
        matrix: matrix_state.as_ref(),
    });
    params.gapping_params.context.set(old);
    result
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:305-356
// ```c++
// s_HSPListFromDistinctAlignments(BlastHSPList *hsp_list,
//                                 BlastCompo_Alignment ** alignments,
//                                 int oid,
//                                 const BlastQueryInfo* queryInfo,
//                                 int frame)
// {
//     int status = 0;                    /* return code for any routine called */
//     static const int unknown_value = 0;   /* dummy constant to use when a
//                                              parameter value is not known */
//     BlastCompo_Alignment * align;  /* an alignment in the list */
//
//     if (hsp_list == NULL) {
//         return -1;
//     }
//     hsp_list->oid = oid;
//
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
//         if (status != 0)
//             break;
//         /* At this point, the subject and possibly the query sequence have
//          * been filtered; since it is not clear that num_ident of the
//          * filtered sequences, rather than the original, is desired,
//          * explicitly leave num_ident blank. */
//         new_hsp->num_ident = 0;
//
//         status = Blast_HSPListSaveHSP(hsp_list, new_hsp);
//         if (status != 0)
//             break;
//     }
//     if (status == 0) {
//         BlastCompo_AlignmentsFree(alignments, s_FreeEditScript);
// ```
#[derive(Clone, Debug)]
pub struct RedoneHsp {
    pub hsp: PreliminaryHsp,
    pub edit_script: Vec<crate::common::GapEditOp>,
    pub composition_method: i32,
    pub identity_matrix: Option<std::sync::Arc<AdjustedProteinMatrix>>,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:223-278
// ```c++
// static void
// s_HitlistReapContained(BlastHSP * hsp_array[], Int4 * hspcnt)
// {
//     Int4 iread;       /* iteration index used to read the hitlist */
//     Int4 iwrite;      /* iteration index used to write to the hitlist */
//     Int4 old_hspcnt;  /* number of HSPs in the hitlist on entry */
//
//     old_hspcnt = *hspcnt;
//
//     for (iread = 1;  iread < *hspcnt;  iread++) {
//         /* for all HSPs in the hitlist */
//         Int4      ireadBack;  /* iterator over indices less than iread */
//         BlastHSP *hsp1;       /* an HSP that is a candidate for deletion */
//
//         hsp1 = hsp_array[iread];
//         for (ireadBack = 0;  ireadBack < iread && hsp1 != NULL;  ireadBack++) {
//             /* for all HSPs before hsp1 in the hitlist and while hsp1
//              * has not been deleted */
//             BlastHSP *hsp2;    /* an HSP that occurs earlier in hsp_array
//                                 * than hsp1 */
//             hsp2 = hsp_array[ireadBack];
//
//             if( hsp2 == NULL ) {  /* hsp2 was deleted in a prior iteration. */
//                 continue;
//             }
//             if (hsp2->query.frame == hsp1->query.frame &&
//                 hsp2->subject.frame == hsp1->subject.frame) {
//                 /* hsp1 and hsp2 are in the same query/subject frame. */
//                 if (CONTAINED_IN_HSP
//                     (hsp2->query.offset, hsp2->query.end, hsp1->query.offset,
//                      hsp2->subject.offset, hsp2->subject.end,
//                      hsp1->subject.offset) &&
//                     CONTAINED_IN_HSP
//                     (hsp2->query.offset, hsp2->query.end, hsp1->query.end,
//                      hsp2->subject.offset, hsp2->subject.end,
//                      hsp1->subject.end)    &&
//                     hsp1->score <= hsp2->score) {
//                     hsp1 = hsp_array[iread] = Blast_HSPFree(hsp_array[iread]);
//                 }
//             } /* end if hsp1 and hsp2 are in the same query/subject frame */
//         } /* end for all HSPs before hsp1 in the hitlist */
//     } /* end for all HSPs in the hitlist */
//
//     /* Condense the hsp_array, removing any NULL items. */
//     iwrite = 0;
//     for (iread = 0;  iread < *hspcnt;  iread++) {
//         if (hsp_array[iread] != NULL) {
//             hsp_array[iwrite++] = hsp_array[iread];
//         }
//     }
//     *hspcnt = iwrite;
//     /* Fill the remaining memory in hsp_array with NULL pointers. */
//     for ( ;  iwrite < old_hspcnt;  iwrite++) {
//         hsp_array[iwrite] = NULL;
//     }
// }
// ```
// BLASTX's untranslated protein subject has frame 0 throughout this owner;
// the NCBI subject-frame equality predicate is therefore always true here.
fn reap_contained_hsps(
    output: Vec<RedoneHsp>,
    trace: &mut dyn for<'a> FnMut(BlastRedoTraceEvent<'a>),
) -> Vec<RedoneHsp> {
    let mut slots: Vec<Option<RedoneHsp>> = output.into_iter().map(Some).collect();
    for (index, candidate) in slots.iter().enumerate() {
        let h = &candidate.as_ref().expect("live input HSP").hsp;
        trace(BlastRedoTraceEvent::ReapList {
            phase: "IN",
            index,
            context: h.context,
            frame: h.frame,
            score: h.score,
            q_start: h.q_start,
            q_end: h.q_end,
            s_start: h.s_start,
            s_end: h.s_end,
        });
    }
    for index in 1..slots.len() {
        let h = &slots[index].as_ref().expect("live candidate HSP").hsp;
        let mut contained = false;
        for (previous_index, other) in slots[..index].iter().enumerate() {
            let Some(other) = other else { continue };
            let p = &other.hsp;
            let result = p.frame == h.frame
                && p.q_start <= h.q_start
                && h.q_start <= p.q_end
                && p.s_start <= h.s_start
                && h.s_start <= p.s_end
                && p.q_start <= h.q_end
                && h.q_end <= p.q_end
                && p.s_start <= h.s_end
                && h.s_end <= p.s_end
                && h.score <= p.score;
            trace(BlastRedoTraceEvent::ReapCompare {
                index,
                previous_index,
                result,
            });
            if result {
                contained = true;
                break;
            }
        }
        trace(BlastRedoTraceEvent::ReapContained {
            frame: h.frame,
            score: h.score,
            q_start: h.q_start,
            q_end: h.q_end,
            s_start: h.s_start,
            s_end: h.s_end,
            result: contained,
        });
        if contained {
            slots[index] = None;
        }
    }
    let retained: Vec<RedoneHsp> = slots.into_iter().flatten().collect();
    for (index, candidate) in retained.iter().enumerate() {
        let h = &candidate.hsp;
        trace(BlastRedoTraceEvent::ReapList {
            phase: "OUT",
            index,
            context: h.context,
            frame: h.frame,
            score: h.score,
            q_start: h.q_start,
            q_end: h.q_end,
            s_start: h.s_start,
            s_end: h.s_end,
        });
    }
    retained
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3590-3675
// ```c++
//                 hsp_list = Blast_HSPListNew(0);
//                 for (frame_index = 0;
//                         frame_index < numFrames;
//                         frame_index++, context_index++) {
//                     incoming_aligns = incoming_align_set[frame_index];
//                     if (!incoming_aligns) {
//                         continue;
//                     }
//                     /*
//                      * All alignments in thisMatch should be to the same query
//                      */
//                     kbp = sbp->kbp_gap[context_index];
//                     if (smithWaterman) {
//                         *pStatusCode =
//                                 Blast_RedoOneMatchSmithWaterman(
//                                         alignments,
//                                         redo_align_params,
//                                         incoming_aligns,
//                                         numAligns[frame_index],
//                                         kbp->Lambda,
//                                         kbp->logK,
//                                         &matchingSeq,
//                                         query_info,
//                                         numQueries,
//                                         matrix,
//                                         BLASTAA_SIZE,
//                                         NRrecord,
//                                         forbidden,
//                                         redoneMatches,
//                                         &pvalueForThisPair,
//                                         compositionTestIndex,
//                                         &LambdaRatio
//                                 );
//                     } else {
//                         *pStatusCode =
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
//                     if (alignments[context_index] != NULL) {
//                         Int2 qframe = frame_index;
//                         if (program_number == eBlastTypeBlastx) {
//                             if (qframe < 3) {
//                                 qframe++;
//                             } else {
//                                 qframe = 2 - qframe;
//                             }
//                         }
//                         *pStatusCode =
//                                 s_HSPListFromDistinctAlignments(hsp_list,
//                                         &alignments[context_index],
//                                         matchingSeq.index,
//                                         queryInfo, qframe);
//                         if (*pStatusCode) {
//                             goto match_loop_cleanup;
//                         }
//                     }
//                     BlastCompo_AlignmentsFree(&incoming_aligns, NULL);
//                     incoming_align_set[frame_index] = NULL;
//                 }
//
//                 if (hsp_list->hspcnt > 1) {
//                     s_HitlistReapContained(hsp_list->hsp_array,
//                             &hsp_list->hspcnt);
//                 }
//                 *pStatusCode =
//                         s_HitlistEvaluateAndPurge(&best_score, &best_evalue,
// ```
pub fn redo_list(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
) -> Result<Vec<RedoneHsp>> {
    let mut state = KappaState::new(batch, parameters, options)?;
    state.redo_list(input, batch, subject)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3244-3248
// ```c++
//         query_info_tld[i] = s_GetQueryInfo(
//                 queryBlk->sequence,
//                 queryInfo,
//                 (program_number == eBlastTypeBlastx)
//         );
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3340-3367
// ```c++
//         gapping_params_context_tld[i].gap_align = gap_align_tld[i];
//         gapping_params_context_tld[i].scoringParams = score_params_tld[i];
//         gapping_params_context_tld[i].sbp = sbp_tld[i];
//         gapping_params_context_tld[i].localScalingFactor = localScalingFactor;
//         gapping_params_context_tld[i].prog_number = program_number;
//
//         redo_align_params_tld[i] =
//             s_GetAlignParams(
//                     &gapping_params_context_tld[i],
//                     queryBlk,
//                     queryInfo,
//                     hitParams,
//                     extendParams
//             );
//         if (redo_align_params_tld[i] == NULL) {
//             status_code = -1;
//             goto function_cleanup;
//         }
//
//         if (positionBased) {
//             matrix_tld[i] = sbp_tld[i]->psi_matrix->pssm->data;
//         } else {
//             matrix_tld[i] = sbp_tld[i]->matrix->data;
//         }
//         /**** Validate parameters *************/
//         if (matrix_tld[i] == NULL) {
//             goto function_cleanup;
//         }
// ```
pub struct KappaState {
    params: BlastRedoAlignParams,
    // EXPERIMENT (LOSAT_X_BXPAR): read-only, shared with the per-thread copies.
    infos: std::sync::Arc<Vec<BlastCompoQueryInfo>>,
    lambda: f64,
    scratch: GapAlignScratch,
    workspace: BlastCompositionWorkspace,
    matrix: Option<AdjustedProteinMatrix>,
}
impl KappaState {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3244-3248
    // ```c++
    //         query_info_tld[i] = s_GetQueryInfo(
    //                 queryBlk->sequence,
    //                 queryInfo,
    //                 (program_number == eBlastTypeBlastx)
    //         );
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3340-3367
    // ```c++
    //         gapping_params_context_tld[i].gap_align = gap_align_tld[i];
    //         gapping_params_context_tld[i].scoringParams = score_params_tld[i];
    //         gapping_params_context_tld[i].sbp = sbp_tld[i];
    //         gapping_params_context_tld[i].localScalingFactor = localScalingFactor;
    //         gapping_params_context_tld[i].prog_number = program_number;
    //
    //         redo_align_params_tld[i] =
    //             s_GetAlignParams(
    //                     &gapping_params_context_tld[i],
    //                     queryBlk,
    //                     queryInfo,
    //                     hitParams,
    //                     extendParams
    //             );
    //         if (redo_align_params_tld[i] == NULL) {
    //             status_code = -1;
    //             goto function_cleanup;
    //         }
    //
    //         if (positionBased) {
    //             matrix_tld[i] = sbp_tld[i]->psi_matrix->pssm->data;
    //         } else {
    //             matrix_tld[i] = sbp_tld[i]->matrix->data;
    //         }
    //         /**** Validate parameters *************/
    //         if (matrix_tld[i] == NULL) {
    //             goto function_cleanup;
    //         }
    // ```
    pub fn new(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
    ) -> Result<Self> {
        let params = redo_params(batch, parameters, options)?;
        let infos = query_infos(batch, parameters);
        let ka = lookup_protein_params(&ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
        });
        Ok(Self {
            params,
            infos: std::sync::Arc::new(infos),
            lambda: ka.lambda / 32.0,
            scratch: GapAlignScratch::new(),
            workspace: BlastCompositionWorkspace::new_blosum62(),
            matrix: None,
        })
    }
    /// EXPERIMENT (LOSAT_X_BXPAR): a second state for another thread, as NCBI
    /// c++/src/algo/blast/core/blast_kappa.c:3244-3340 allocates per thread
    /// (`redo_align_params_tld`, `gap_align_tld`, `NRrecord_tld`); the query
    /// information is read-only and shared.
    pub(crate) fn x_for_thread(
        &self,
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
    ) -> Result<Self> {
        Ok(Self {
            params: redo_params(batch, parameters, options)?,
            infos: std::sync::Arc::clone(&self.infos),
            lambda: self.lambda,
            scratch: GapAlignScratch::new(),
            workspace: BlastCompositionWorkspace::new_blosum62(),
            matrix: None,
        })
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3590-3641
    // ```c++
    //                 hsp_list = Blast_HSPListNew(0);
    //                 for (frame_index = 0;
    //                         frame_index < numFrames;
    //                         frame_index++, context_index++) {
    //                     incoming_aligns = incoming_align_set[frame_index];
    //                     if (!incoming_aligns) {
    //                         continue;
    //                     }
    //                     /*
    //                      * All alignments in thisMatch should be to the same query
    //                      */
    //                     kbp = sbp->kbp_gap[context_index];
    //                     if (smithWaterman) {
    //                         *pStatusCode =
    //                                 Blast_RedoOneMatchSmithWaterman(
    //                                         alignments,
    //                                         redo_align_params,
    //                                         incoming_aligns,
    //                                         numAligns[frame_index],
    //                                         kbp->Lambda,
    //                                         kbp->logK,
    //                                         &matchingSeq,
    //                                         query_info,
    //                                         numQueries,
    //                                         matrix,
    //                                         BLASTAA_SIZE,
    //                                         NRrecord,
    //                                         forbidden,
    //                                         redoneMatches,
    //                                         &pvalueForThisPair,
    //                                         compositionTestIndex,
    //                                         &LambdaRatio
    //                                 );
    //                     } else {
    //                         *pStatusCode =
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
    // ```
    pub fn redo_list(
        &mut self,
        input: &[PreliminaryHsp],
        batch: &PreparedQueryBatch,
        subject: &[u8],
    ) -> Result<Vec<RedoneHsp>> {
        self.redo_list_observed(input, batch, subject, &mut |_| {})
    }
    /// EXPERIMENT (LOSAT_X_BXPAR): the matrix a match leaves for the next one.
    pub(crate) fn x_take_matrix(&mut self) -> Option<AdjustedProteinMatrix> {
        self.matrix.take()
    }
    /// EXPERIMENT (LOSAT_X_BXPAR): hand a match the matrix of its predecessor.
    pub(crate) fn x_set_matrix(&mut self, matrix: Option<AdjustedProteinMatrix>) {
        self.matrix = matrix;
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
    pub fn redo_list_observed(
        &mut self,
        input: &[PreliminaryHsp],
        batch: &PreparedQueryBatch,
        subject: &[u8],
        trace: &mut dyn for<'a> FnMut(BlastRedoTraceEvent<'a>),
    ) -> Result<Vec<RedoneHsp>> {
        let query_index = input.first().map(|h| h.context / 6);
        ensure!(
            input.iter().all(|h| Some(h.context / 6) == query_index),
            "BLASTX Kappa requires one query/OID HSP list"
        );
        let mut output = Vec::new();

        for context_index in query_index.unwrap_or(0) * 6..query_index.unwrap_or(0) * 6 + 6 {
            if !input.iter().any(|h| h.context == context_index) {
                continue;
            }
            let result = redo_context_observed(
                input,
                context_index,
                &self.infos,
                subject,
                &self.params,
                self.lambda,
                &mut self.scratch,
                &mut self.workspace,
                &mut self.matrix,
                trace,
            )?;
            for (context, alignments) in result.alignments_by_query.into_iter().enumerate() {
                let mut current = alignments;
                while let Some(mut align) = current {
                    current = align.next.take();
                    let Some(BlastCompoAlignmentContext::EditScript(script)) = align.context.take()
                    else {
                        bail!("BLASTX Kappa alignment edit script missing")
                    };
                    let method = match align.matrix_adjust_rule {
                        EMatrixAdjustRule::DontAdjustMatrix => 0,
                        EMatrixAdjustRule::CompoScaleOldMatrix => 1,
                        _ => 2,
                    };
                    output.push(RedoneHsp {
                        hsp: PreliminaryHsp {
                            context,
                            frame: batch.contexts[context].frame,
                            score: align.score,
                            q_start: align.query_start,
                            q_end: align.query_end,
                            q_gapped_start: 0,
                            s_start: align.match_start,
                            s_end: align.match_end,
                            s_gapped_start: 0,
                        },
                        edit_script: script,
                        composition_method: method,
                        identity_matrix: None,
                    });
                }
            }
            output.sort_by(|a, b| super::preliminary::compare_score(&a.hsp, &b.hsp));
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:238-266
        // ```c++
        //         for (ireadBack = 0;  ireadBack < iread && hsp1 != NULL;  ireadBack++) {
        //             /* for all HSPs before hsp1 in the hitlist and while hsp1
        //              * has not been deleted */
        //             BlastHSP *hsp2;    /* an HSP that occurs earlier in hsp_array
        //                                 * than hsp1 */
        //             hsp2 = hsp_array[ireadBack];
        //
        //             if( hsp2 == NULL ) {  /* hsp2 was deleted in a prior iteration. */
        //                 continue;
        //             }
        //             if (hsp2->query.frame == hsp1->query.frame &&
        //                 hsp2->subject.frame == hsp1->subject.frame) {
        //                 /* hsp1 and hsp2 are in the same query/subject frame. */
        //                 if (CONTAINED_IN_HSP
        //                     (hsp2->query.offset, hsp2->query.end, hsp1->query.offset,
        //                      hsp2->subject.offset, hsp2->subject.end,
        //                      hsp1->subject.offset) &&
        //                     CONTAINED_IN_HSP
        //                     (hsp2->query.offset, hsp2->query.end, hsp1->query.end,
        //                      hsp2->subject.offset, hsp2->subject.end,
        //                      hsp1->subject.end)    &&
        //                     hsp1->score <= hsp2->score) {
        //                     hsp1 = hsp_array[iread] = Blast_HSPFree(hsp_array[iread]);
        //                 }
        //             } /* end if hsp1 and hsp2 are in the same query/subject frame */
        //         } /* end for all HSPs before hsp1 in the hitlist */
        //     } /* end for all HSPs in the hitlist */
        //
        //     /* Condense the hsp_array, removing any NULL items. */
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3689-3703
        // ```c++
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
        let final_matrix = self
            .matrix
            .as_ref()
            .map(|matrix| std::sync::Arc::new(matrix.clone()));
        for h in &mut output {
            h.identity_matrix = final_matrix.clone();
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3670-3673
        // ```c++
        // if (hsp_list->hspcnt > 1) {
        //     s_HitlistReapContained(hsp_list->hsp_array, &hsp_list->hspcnt);
        // }
        // ```
        if output.len() > 1 {
            Ok(reap_contained_hsps(output, trace))
        } else {
            Ok(output)
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:248-260,266-276
// ```c++
//             if (hsp2->query.frame == hsp1->query.frame &&
//                 hsp2->subject.frame == hsp1->subject.frame) {
//                 /* hsp1 and hsp2 are in the same query/subject frame. */
//                 if (CONTAINED_IN_HSP
//                     (hsp2->query.offset, hsp2->query.end, hsp1->query.offset,
//                      hsp2->subject.offset, hsp2->subject.end,
//                      hsp1->subject.offset) &&
//                     CONTAINED_IN_HSP
//                     (hsp2->query.offset, hsp2->query.end, hsp1->query.end,
//                      hsp2->subject.offset, hsp2->subject.end,
//                      hsp1->subject.end)    &&
//                     hsp1->score <= hsp2->score) {
//                     hsp1 = hsp_array[iread] = Blast_HSPFree(hsp_array[iread]);
//     /* Condense the hsp_array, removing any NULL items. */
//     iwrite = 0;
//     for (iread = 0;  iread < *hspcnt;  iread++) {
//         if (hsp_array[iread] != NULL) {
//             hsp_array[iwrite++] = hsp_array[iread];
//         }
//     }
//     *hspcnt = iwrite;
//     /* Fill the remaining memory in hsp_array with NULL pointers. */
//     for ( ;  iwrite < old_hspcnt;  iwrite++) {
//         hsp_array[iwrite] = NULL;
// ```
#[cfg(test)]
mod contained_reap_tests {
    use super::*;

    #[test]
    fn final_reap_matches_direct_pinned_c_boundaries() {
        let text = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/unit/blastx_contained_reap_expected_v020.tsv"
        ));
        let mut cases: std::collections::BTreeMap<&str, Vec<(RedoneHsp, i32)>> =
            std::collections::BTreeMap::new();
        for line in text
            .lines()
            .filter(|line| !line.starts_with('#') && !line.is_empty())
        {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 10);
            let n: Vec<i32> = fields[1..]
                .iter()
                .map(|x| x.parse().expect("C integer field"))
                .collect();
            let entries = cases.entry(fields[0]).or_default();
            assert_eq!(n[0] as usize, entries.len(), "C PRE slot order");
            entries.push((
                RedoneHsp {
                    hsp: PreliminaryHsp {
                        context: n[1] as usize,
                        score: n[2],
                        frame: n[3] as i8,
                        q_start: n[4],
                        q_end: n[5],
                        q_gapped_start: 0,
                        s_start: n[6],
                        s_end: n[7],
                        s_gapped_start: 0,
                    },
                    edit_script: Vec::new(),
                    composition_method: 0,
                    identity_matrix: None,
                },
                n[8],
            ));
        }
        assert_eq!(cases.len(), 21);
        for (name, entries) in cases {
            let mut expected: Vec<(i32, PreliminaryHsp)> = entries
                .iter()
                .filter(|(_, slot)| *slot >= 0)
                .map(|(h, slot)| (*slot, h.hsp.clone()))
                .collect();
            expected.sort_by_key(|(slot, _)| *slot);
            let input = entries.into_iter().map(|(h, _)| h).collect();
            let actual = reap_contained_hsps(input, &mut |_| {});
            let actual_hsps: Vec<PreliminaryHsp> = actual.into_iter().map(|h| h.hsp).collect();
            let expected_hsps: Vec<PreliminaryHsp> = expected.into_iter().map(|(_, h)| h).collect();
            assert_eq!(actual_hsps, expected_hsps, "direct C boundary case {name}");
        }
    }
}
