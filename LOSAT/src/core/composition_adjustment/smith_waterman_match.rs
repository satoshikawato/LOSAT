//! The Smith-Waterman redo of one match (`-use_sw_tback`; series RG, items X4 and X4b):
//! NCBI's `Blast_RedoOneMatchSmithWaterman` with a square matrix, one body for BLASTP's protein
//! windows and TBLASTN's translated subject windows, in place of the X-drop redo
//! (`Blast_RedoOneMatch`, `redo_alignment::blast_redo_one_match_with_workspace_queries_and_matrix`).
//! The windows (`s_WindowsFromAligns`), the subject range (`get_range`), the composition adjustment
//! (`Blast_AdjustScores`) and the alignment list it returns are the X-drop redo's; the round loop of a
//! window (the score-only scan, the significance test, the regime of a significant score of 0, the
//! start scan, the X-drop traceback from the Smith-Waterman start and the forbidden ranges;
//! `smith_waterman::smith_waterman_redo_window`) is the one RPS-BLAST's redo runs (line PS, PS11c-3),
//! on the rows that the square matrix gives the query's letters (`matrix[query[queryPos]]`, which
//! `ALIGN_EX` and the kernels read as they read a PSSM's `pssm[queryPos]`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039
//! ```c
//!     Boolean smithWaterman =
//!         (Boolean) (extendParams->options->eTbackExt == eSmithWatermanTbck);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3602-3622
//! ```c
//!                     if (smithWaterman) {
//!                         *pStatusCode =
//!                                 Blast_RedoOneMatchSmithWaterman(
//!                                         alignments,
//!                                         redo_align_params,
//!                                         incoming_aligns,
//!                                         numAligns[frame_index],
//!                                         kbp->Lambda,
//!                                         kbp->logK,
//!                                         &matchingSeq,
//!                                         query_info,
//!                                         numQueries,
//!                                         matrix,
//!                                         BLASTAA_SIZE,
//!                                         NRrecord,
//!                                         forbidden,
//!                                         redoneMatches,
//!                                         &pvalueForThisPair,
//!                                         compositionTestIndex,
//!                                         &LambdaRatio
//!                                 );
//! ```
//! The steps after it (`s_HSPListFromDistinctAlignments`, the reap, `s_HitlistEvaluateAndPurge`,
//! the normalization and the identities) are the X-drop redo's (blast_kappa.c:3654-3707): the
//! callers' (BLASTP's `postprocess_preliminary_hits_smith_waterman`, TBLASTN's stage D).

use anyhow::{bail, Result};

use crate::algorithm::pssm_search::gapalign::{pssm_sw_xdrop_traceback, PssmGapScratch};
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::core::composition_adjustment::adjust_scores::{
    blast_adjust_scores, AdjustedProteinMatrix, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    blast_compo_alignment_new, preliminary_test_near_identical, redo_calc_lambda,
    translated_subject_composition, windows_from_protein_aligns, windows_from_translated_aligns,
    x_whole_subject_composition, BlastCompoAlignment, BlastCompoAlignmentContext,
    BlastCompoMatchingSequence, BlastCompoQueryInfo, BlastRedoAlignCallbacks, BlastRedoAlignParams,
    BlastRedoOneMatchResult, EMatrixAdjustRule, WINDOW_BORDER,
};
use crate::core::composition_adjustment::smith_waterman::{
    smith_waterman_redo_window, SwRedoScratch, SwSignificance, SwWindow, SwWindowRedo,
};

/// The square matrix's memory of its windows (`SwRedoScratch::window`): the score rows of the
/// window's query letters in the window's matrix (`matrix[query[queryPos]]`) and the X-drop
/// traceback's memory (NCBI's `gap_align`).
pub(crate) struct SquareWindow {
    rows: QueryPssm,
    gap: PssmGapScratch,
}

impl Default for SquareWindow {
    fn default() -> Self {
        Self {
            rows: QueryPssm { rows: Vec::new() },
            gap: PssmGapScratch::new(),
        }
    }
}

/// The memory of the square matrix's Smith-Waterman redo from one match to the next: the round
/// loop's (`SwRedoScratch::rounds`, shared with RPS-BLAST's redo) and its windows' (`SquareWindow`).
pub(crate) type SwMatchScratch = SwRedoScratch<SquareWindow>;

/// NCBI's failed Smith-Waterman redo of a match: a significant alignment of score 0 whose X-drop
/// traceback has no operation (`GapEditScriptNew(0)` is NULL, `s_NewAlignmentUsingXdrop` returns
/// -1), in any window of the call. NCBI frees the match's alignments (of every window of the call)
/// and keeps the status -1 (`*pStatusCode`) until the next match that it redoes sets it again; a
/// status left by the last match of a query batch fails the traceback (`Traceback failed`,
/// traceback_stage.cpp:251-252), or crashes NCBI's clean-up when the batch's first query holds a
/// redone match (`sw_batch_failed_error`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1542-1548
/// ```c
/// function_level_cleanup:
///     if (status != 0) {
///         for (query_index = 0;  query_index < numQueries;  query_index++) {
///             BlastCompo_AlignmentsFree(&alignments[query_index],
///                                       callbacks->free_align_traceback);
///         }
///     }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SwRedoFailed;

impl std::fmt::Display for SwRedoFailed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("the Smith-Waterman redo of a match makes an alignment without operations")
    }
}

impl std::error::Error for SwRedoFailed {}

/// NCBI's Smith-Waterman redo of a match that does not end: a score of 0 is significant, the
/// window has two or more HSPs and the first query letter's row scores the window's first subject
/// letter above 0 (`SwWindowRedo::Endless`, the regime of RPS-BLAST's `SwEndless`;
/// PD-LOSAT-NCBI-DEFECTS 1.53, within the class of 1.30). The window that reaches it decides; the
/// windows before it in the call belong to the same match, which NCBI never finishes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SwRedoEndless;

impl std::fmt::Display for SwRedoEndless {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("the Smith-Waterman redo of a match does not end")
    }
}

impl std::error::Error for SwRedoEndless {}

/// `Blast_RedoOneMatchSmithWaterman` with a square matrix, for BLASTP's protein windows
/// (`s_WindowsFromProteinAligns`) and TBLASTN's translated subject windows
/// (`s_WindowsFromTranslatedAligns`, `params.subject_is_translated`): each window's alignments in
/// NCBI's list order (the last one made first) in `alignments_by_query[query_index]`, as the X-drop
/// redo returns them, with the composition p-value and lambda ratio of the last adjustment.
/// `lambda` and `log_k` are the match's scaled block (`kbp_gap[context_index]` of the redo's copy);
/// `would_insert(query_index, evalue, score)` is `BlastCompo_HeapWouldInsert(&significantMatches
/// [query_index], ...)` with the match's subject index, asked before the query's first alignment of
/// the call in e-value mode (an error from it ends the call with that error); `cutoff_e` is
/// `-evalue` (`params->cutoff_e`). Without composition-based statistics no window is adjusted: its
/// rows are those of `unadjusted`, NCBI's `matrix` argument as the redo's set-up leaves it
/// (`sbp->matrix->data`, with a scaling factor of 1), and its rule `eDontAdjustMatrix`; the subject
/// range is not masked by SEG (`get_range`). The errors `SwRedoFailed` and `SwRedoEndless` are
/// NCBI's failed and endless redos (see there).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1345-1357
/// ```c
///     assert((int) compo_adjust_mode < 2 || !positionBased);
///     for (query_index = 0;  query_index < numQueries;  query_index++) {
///         alignments[query_index] = NULL;
///     }
///     /* Find the multiple translation windows used by tblastn queries. */
///     status =
///         s_WindowsFromAligns(incoming_aligns, query_info, hspcnt, numQueries,
///                             kWindowBorder, matchingSeq->length, &windows,
///                             &nWindows, query_is_translated, subject_is_translated, params->positionBased);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1366-1430
/// ```c
///         window = windows[window_index];
///         query_index = window->query_range.context;
///         query_composition = &query_info[query_index].composition;
///
///         nearIdenticalStatus = s_preliminaryTestNearIdentical(query_info,
///                                              window,
///                                              window->align,
///                                              params->near_identical_cutoff);
///
///         status =
///             callbacks->get_range(matchingSeq, &window->subject_range,
///                                  &subject,
///                                  &query_info[query_index].seq,
///                                  &window->query_range,
///                                  &query,
///                                  query_info[query_index].words,
///                                  window->align, nearIdenticalStatus,
///                                  compo_adjust_mode, TRUE,
///                                  NULL);
///         ...
///         searchsp = query_info[query_index].eff_search_space;
///
///         /* For Smith-Waterman alignments, adjust the search using the
///          * composition of the highest scoring alignment in window */
///         if (compo_adjust_mode != eNoCompositionBasedStats) {
///             Blast_AminoAcidComposition subject_composition;
///             s_GetComposition(&subject_composition, alphsize,
///                                     &subject, &window->subject_range,
///                                     window->align, FALSE, subject_is_translated);
///             adjust_search_failed =
///                 Blast_AdjustScores(matrix,
///                                    query_composition, query.length,
///                                    &subject_composition, subject.length,
///                                    scaledMatrixInfo, compo_adjust_mode,
///                                    RE_pseudocounts, NRrecord,
///                                    &matrix_adjust_rule, callbacks->calc_lambda,
///                                    pvalueForThisPair,
///                                    compositionTestIndex,
///                                    LambdaRatio);
///             ...
///         }
///         if ( !adjust_search_failed ) {
/// ```
/// The subject range's last argument is `NULL`: SEG runs (with composition-based statistics) unless
/// the window's first alignment is near identical, as for the X-drop redo's first range
/// (`subject_maybe_biased` starts true there and this redo does not read its answer).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3359-3362 (`matrix`, which
/// `Blast_AdjustScores` overwrites, is the thread's score block's matrix)
/// ```c
///         if (positionBased) {
///             matrix_tld[i] = sbp_tld[i]->psi_matrix->pssm->data;
///         } else {
///             matrix_tld[i] = sbp_tld[i]->matrix->data;
/// ```
#[allow(clippy::too_many_arguments)]
pub(crate) fn redo_one_match_smith_waterman<'seq>(
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    params: &BlastRedoAlignParams,
    lambda: f64,
    log_k: f64,
    matching_seq: &BlastCompoMatchingSequence<'seq>,
    query_infos: &[BlastCompoQueryInfo],
    callbacks: &BlastRedoAlignCallbacks,
    composition_workspace: &mut BlastCompositionWorkspace,
    matrix_state: &mut Option<AdjustedProteinMatrix>,
    unadjusted: Option<&AdjustedProteinMatrix>,
    cutoff_e: f64,
    would_insert: &mut dyn FnMut(usize, f64, i32) -> Result<bool>,
    scratch: &mut SwMatchScratch,
) -> Result<BlastRedoOneMatchResult> {
    // BLASTX's translated query windows (item X4b-2) and a PSSM's rows (psiblast, tblastn
    // `-in_pssm`) do not reach this redo: their options refuse `-use_sw_tback`.
    if params.query_is_translated || params.position_based {
        bail!(
            "the Smith-Waterman redo of translated query or position-based windows is not ported"
        );
    }
    let uses_composition = params.uses_composition_based_stats();
    let unadjusted = match unadjusted {
        Some(matrix) => Some(matrix),
        None if uses_composition => None,
        None => bail!("the Smith-Waterman redo without composition-based statistics needs the search's matrix"),
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:879-888
    // ```c++
    //     if (subject_is_translated || query_is_translated) {
    //         return s_WindowsFromTranslatedAligns(alignments, query_info,
    //                                              hspcnt, border,
    //                                              sequence_length,
    //                                              pwindows, nWindows,
    //                                              subject_is_translated, is_pos_based);
    //     } else {
    //         return s_WindowsFromProteinAligns(alignments, query_info,
    //                                           numQueries, sequence_length,
    //                                           pwindows, nWindows);
    // ```
    let windows = if params.subject_is_translated {
        if matching_seq.genetic_code_id.is_none() {
            bail!("the Smith-Waterman redo of a translated subject needs its genetic code");
        }
        let query_lengths: Vec<_> = query_infos.iter().map(|info| info.seq.length).collect();
        windows_from_translated_aligns(
            incoming_aligns,
            &query_lengths,
            matching_seq.length,
            WINDOW_BORDER,
            params.position_based,
        )?
    } else {
        let query_info = query_infos
            .first()
            .ok_or_else(|| anyhow::anyhow!("BLASTP redo requires query information"))?;
        windows_from_protein_aligns(incoming_aligns, query_info, matching_seq.length)?
    };
    let mut alignments_by_query: Vec<Option<Box<BlastCompoAlignment>>> = Vec::new();
    let mut pvalue_for_this_pair = None;
    let mut lambda_ratio = None;
    let mut adjusted_matrix = matrix_state.take();
    let gap_open = params.gapping_params.gap_open;
    let gap_extend = params.gapping_params.gap_extend;

    for window in &windows {
        let query_index = usize::try_from(window.query_range.context)?;
        // The translated windows name each query of the call (`query_info[query_index]`); BLASTP
        // passes its one query (`windows_from_protein_aligns`).
        let query_info = if params.subject_is_translated {
            query_infos.get(query_index).ok_or_else(|| {
                anyhow::anyhow!("the Smith-Waterman redo's query {query_index} is missing")
            })?
        } else {
            &query_infos[0]
        };
        if alignments_by_query.len() <= query_index {
            alignments_by_query.resize_with(query_index + 1, || None);
        }
        let Some(window_align) = window.align.as_deref() else {
            continue;
        };
        let near_identical = preliminary_test_near_identical(
            query_info,
            window,
            window_align,
            params.near_identical_cutoff,
        );
        let range = (callbacks.get_range)(
            matching_seq,
            &window.subject_range,
            &query_info.seq,
            &window.query_range,
            query_info.words.as_deref(),
            window_align,
            near_identical,
            true,
            params,
        )?;
        let searchsp = query_info.eff_search_space;
        // NCBI's `matrix[query[queryPos]]` for each position of the window's query: the rows that
        // the kernels and `ALIGN_EX` read (`positionBased` false, blast_gapalign.c:543-554).
        let query = range.query.data();
        let SquareWindow { rows, gap } = &mut scratch.window;
        let fill_rows = |rows: &mut QueryPssm, matrix: &AdjustedProteinMatrix| {
            rows.rows.clear();
            rows.rows.extend(
                query
                    .iter()
                    .map(|&letter| matrix.scores[usize::from(letter)]),
            );
        };
        let matrix_adjust_rule = match unadjusted {
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1332
            // ```c
            //     EMatrixAdjustRule matrix_adjust_rule = eDontAdjustMatrix;
            // ```
            Some(matrix) if !uses_composition => {
                fill_rows(rows, matrix);
                EMatrixAdjustRule::DontAdjustMatrix
            }
            _ => {
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:923-936
                // ```c
                //     if (query_is_translated || subject_is_translated) {
                //         int start;
                //         int end;
                //         start = ((query_is_translated) ?
                //                  align->queryStart : align->matchStart) - range->begin;
                //         end   = ((query_is_translated) ?
                //                  align->queryEnd   : align->matchEnd  ) - range->begin;
                //         Blast_GetCompositionRange(&left, &right, data, length, start, end);
                //     } else {
                //         /* Use the whole subject to compute the composition */
                //         left = 0;
                //         right = length;
                //     }
                //     Blast_ReadAaComposition(composition, alphsize, &data[left], right-left);
                // ```
                let subject_composition = if params.subject_is_translated {
                    translated_subject_composition(
                        &range.subject.buffer,
                        &window.subject_range,
                        window_align,
                    )?
                } else {
                    x_whole_subject_composition(&range)
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1439-1444
                // ```c++
                //     if (query_composition->numTrueAminoAcids == 0 ||
                //         subject_composition->numTrueAminoAcids == 0) {
                //         /* Either the query or subject contains only ambiguity
                //            characters, most likely because the entire subject has been
                //            SEGed.  Compositional adjustment is meaningless. */
                //         return 1;
                // ```
                if query_info.composition.num_true_amino_acids == 0
                    || subject_composition.num_true_amino_acids == 0
                {
                    continue;
                }
                let Some(adjusted) = blast_adjust_scores(
                    &params.matrix_info,
                    &query_info.composition,
                    range.query.length,
                    &subject_composition,
                    range.subject.length,
                    params.compo_adjust_mode,
                    params.composition_test_index,
                    composition_workspace,
                    callbacks.calc_lambda.unwrap_or(redo_calc_lambda),
                )?
                else {
                    continue;
                };
                pvalue_for_this_pair = adjusted.pvalue_for_this_pair;
                lambda_ratio = Some(adjusted.lambda_ratio);
                fill_rows(rows, adjusted_matrix.insert(adjusted.adjusted_matrix));
                adjusted.matrix_adjust_rule
            }
        };
        let rows = &*rows;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3261-3264
        // ```c
        //         if(smithWaterman) {
        // 	   forbidden_tld[i] = calloc(1,sizeof(Blast_ForbiddenRanges));
        //            status_code =
        //                Blast_ForbiddenRangesInitialize(forbidden_tld[i], queryInfo->max_length);
        // ```
        // The round loop (`smith_waterman_redo_window`, redo_alignment.c:1438-1529) on the
        // subject's range from NCBI's `subject->data - 1`, which `ALIGN_EX` reads from
        // (`s_SWFindFinalEndsUsingXdrop`): the sentinel before the window's letters
        // (`s_SequenceGetProteinRange`, `s_SequenceGetTranslatedRange`). `Lambda` and `logK` are the
        // redo's scaled block (`kbp_gap` of NCBI's copy), `searchsp` the query's
        // `eff_search_space`, `cutoff_s` the redo's cutoff score.
        let sw_window = SwWindow {
            subject_from_sentinel: &range.subject.buffer[range.subject.data_offset - 1..],
            subject_length: range.subject.data().len(),
            query_length: query.len(),
            rows: |pos| &rows.rows[pos],
            gap_open,
            gap_extend,
            // NCBI's `int` count, never negative (the window functions count HSPs up).
            hspcnt: usize::try_from(window.hspcnt).expect("a window's HSP count"),
        };
        let significance = if params.do_link_hsps {
            SwSignificance::CutoffScore(params.cutoff_score)
        } else {
            SwSignificance::Evalue {
                lambda,
                log_k,
                searchsp,
                cutoff_e,
            }
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1465-1474
        // ```c
        //                     if (alignments[query_index] == NULL) {
        //                         /* this is the most significant alignment; if
        //                          * it will not be accepted, no alignments from
        //                          * this match will */
        //                         alignment_is_significant =
        //                             alignment_is_significant &&
        //                             BlastCompo_HeapWouldInsert(
        //                                 &significantMatches[query_index],
        //                                 newSwEvalue, aSwScore, matchingSeq->index);
        //                     }
        // ```
        // The query's list of the whole call: the windows of one query (TBLASTN's frames and
        // regions) follow each other in it.
        let query_had_none = alignments_by_query[query_index].is_none();
        let mut heap_error = None;
        let x_dropoff = params.gapping_params.x_dropoff;
        let redone = smith_waterman_redo_window(
            &sw_window,
            significance,
            query_had_none,
            |evalue, score| match would_insert(query_index, evalue, score) {
                Ok(inserts) => inserts,
                Err(error) => {
                    heap_error = Some(error);
                    false
                }
            },
            |query_span, subject_span, score| {
                pssm_sw_xdrop_traceback(
                    rows,
                    sw_window.subject_from_sentinel,
                    query_span,
                    subject_span,
                    gap_open,
                    gap_extend,
                    x_dropoff,
                    score,
                    gap,
                )
            },
            &mut scratch.rounds,
        );
        if let Some(error) = heap_error {
            return Err(error);
        }
        let made = match redone {
            SwWindowRedo::Alignments(made) => made,
            SwWindowRedo::Failed => return Err(SwRedoFailed.into()),
            SwWindowRedo::Endless => return Err(SwRedoEndless.into()),
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1856-1873
        // ```c
        //     if (editScript != NULL) {
        //         /* Shifted values of the endpoints */
        //         Int4 aqueryStart =  queryStart + query_range->begin;
        //         Int4 aqueryEnd   = *pqueryEnd  + query_range->begin;
        //         Int4 amatchStart =  matchStart + subject_range->begin;
        //         Int4 amatchEnd   = *pmatchEnd  + subject_range->begin;
        //
        //         obj = BlastCompo_AlignmentNew(newScore, matrix_adjust_rule,
        //                                       aqueryStart, aqueryEnd,
        //                                       query_range->context,
        //                                       amatchStart, amatchEnd,
        //                                       subject_range->context, editScript);
        // ```
        // Each alignment prepended in the order made (`newAlign->next = alignments[query_index]`).
        for align in made {
            let mut new_align = blast_compo_alignment_new(
                align.score,
                matrix_adjust_rule,
                align.query_start + window.query_range.begin,
                align.query_end + window.query_range.begin,
                window.query_range.context,
                align.match_start + window.subject_range.begin,
                align.match_end + window.subject_range.begin,
                window.subject_range.context,
                Some(BlastCompoAlignmentContext::EditScript(align.edit_script)),
            );
            new_align.next = alignments_by_query[query_index].take();
            alignments_by_query[query_index] = Some(new_align);
        }
    }
    *matrix_state = adjusted_matrix;
    Ok(BlastRedoOneMatchResult {
        alignments_by_query,
        pvalue_for_this_pair,
        lambda_ratio,
    })
}

/// The redo's `matrix` without composition-based statistics: the score block's matrix
/// (`sbp->matrix->data`, every pair of NCBIstdaa letters, `utils::matrix::standard_score`), which
/// the redo's set-up does not scale (a scaling factor of 1) and no window adjusts.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3100-3109
/// ```c
///     if (compo_adjust_mode != eNoCompositionBasedStats) {
///         if((0 == strcmp(scoringParams->options->matrix, "BLOSUM62_20"))) {
///             localScalingFactor = SCALING_FACTOR / 10;
///         } else {
///             localScalingFactor = SCALING_FACTOR;
///         }
///     } else {
///         localScalingFactor = 1.0;
///     }
///     s_RescaleSearch(sbp, scoringParams, numContexts, localScalingFactor);
/// ```
pub(crate) fn unadjusted_square_matrix(
    matrix: crate::config::ScoringMatrix,
) -> AdjustedProteinMatrix {
    let mut scores = [[0; crate::utils::matrix::BLASTAA_SIZE]; crate::utils::matrix::BLASTAA_SIZE];
    for (query, row) in (0u8..).zip(scores.iter_mut()) {
        for (subject, score) in (0u8..).zip(row.iter_mut()) {
            *score = crate::utils::matrix::standard_score(matrix, query, subject);
        }
    }
    AdjustedProteinMatrix { scores }
}

/// The error of a query batch whose last redone match failed its Smith-Waterman redo
/// (`SwRedoFailed`), for `program` (NCBI's name, `blastp`): NCBI's `Traceback failed` (exit 3),
/// or, when the batch's first query holds a redone match (`first_query_holds_match`), LOSAT's
/// refusal of the batch where NCBI's clean-up crashes (PD-LOSAT-NCBI-DEFECTS 1.53, rule 3; exit 1).
/// The caller puts it after the reports of the batches before.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3775-3779
/// ```c
///     for (i = 0; i < actual_num_threads; ++i) {
///         if (status_code_tld[i] != 0) {
///             status_code = status_code_tld[i];
///         }
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:251-253
/// ```c++
///     if (status) {
///         NCBI_THROW(CBlastException, eCoreBlastError, "Traceback failed");
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3795-3802
/// ```c
///         } else {
///             if (redoneMatches_tld[i] != NULL) {
///                 int qi;
///                 for (qi = 0; qi < numQueries; ++qi) {
///                     sfree(redoneMatches_tld[i][qi].array);
///                     sfree(redoneMatches_tld[i][qi].heapArray);
///                 }
///                 s_ClearHeap(redoneMatches_tld[i]);
/// ```
#[cold]
#[inline(never)]
pub(crate) fn sw_batch_failed_error(program: &str, first_query_holds_match: bool) -> anyhow::Error {
    if first_query_holds_match {
        anyhow::anyhow!(
            "a query batch whose last match fails its Smith-Waterman redo (-use_sw_tback) after a match of the batch's first query was redone and kept, with which NCBI BLAST+ 2.17.0's {program} crashes, is not supported by LOSAT's {}",
            program.to_ascii_uppercase()
        )
    } else {
        crate::blastinput::app::engine_error("Traceback failed")
    }
}

/// LOSAT's refusal of a match whose Smith-Waterman redo does not end in NCBI (`SwRedoEndless`;
/// PD-LOSAT-NCBI-DEFECTS 1.53, rule 3, within the class of RPS-BLAST's refusal of 1.30), for
/// `program` (NCBI's name, `blastp`) and what makes the score of 0 significant there (`why`). The
/// caller puts it after the reports of the batches before. Exit status 1, its message.
#[cold]
#[inline(never)]
pub(crate) fn sw_endless_error(program: &str, why: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "a match whose Smith-Waterman redo (-use_sw_tback) does not end in NCBI BLAST+ 2.17.0 ({why}) is not supported by LOSAT's {}",
        program.to_ascii_uppercase()
    )
}

#[cfg(test)]
#[path = "smith_waterman_match_tests.rs"]
mod tests;
