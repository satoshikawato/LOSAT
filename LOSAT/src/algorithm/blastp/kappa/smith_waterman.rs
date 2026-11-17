//! The Smith-Waterman redo of one match with a square matrix (`-use_sw_tback`; series RG, item
//! X4): NCBI's `Blast_RedoOneMatchSmithWaterman` for BLASTP's protein window, in place of the
//! X-drop redo (`Blast_RedoOneMatch`, `blast_redo_one_match_with_workspace_queries_and_matrix`).
//! The window, the subject range (`get_range`), the composition adjustment (`Blast_AdjustScores`)
//! and the alignment list it returns are the X-drop redo's; the kernels (the score-only scan, the
//! start and the forbidden ranges) and the X-drop traceback from the Smith-Waterman start
//! (`s_SWFindFinalEndsUsingXdrop`) are RPS-BLAST's (line PS, PS11c-3), on the rows that the square
//! matrix gives the query's letters (`matrix[query[queryPos]]`, which `ALIGN_EX` and the kernels
//! read as they read a PSSM's `pssm[queryPos]`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039
//! ```c
//!     Boolean smithWaterman =
//!         (Boolean) (extendParams->options->eTbackExt == eSmithWatermanTbck);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3602-3621
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
//! caller's (`postprocess_preliminary_hits_smith_waterman`).

use anyhow::{bail, Result};

use crate::algorithm::pssm_search::gapalign::{pssm_sw_xdrop_traceback, PssmGapScratch};
use crate::algorithm::pssm_search::redo1::smith_waterman::evalue_from_score;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::core::composition_adjustment::adjust_scores::{
    blast_adjust_scores, AdjustedProteinMatrix, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    blast_compo_alignment_new, preliminary_test_near_identical, redo_calc_lambda,
    windows_from_protein_aligns, x_whole_subject_composition, BlastCompoAlignment,
    BlastCompoAlignmentContext, BlastCompoMatchingSequence, BlastCompoQueryInfo,
    BlastRedoAlignCallbacks, BlastRedoAlignParams, BlastRedoOneMatchResult,
};
use crate::core::composition_adjustment::smith_waterman::{
    row_column, smith_waterman_find_start, smith_waterman_score_only, ForbiddenRanges, SwScratch,
};

/// What the Smith-Waterman redo keeps from one match to the next: the forbidden ranges (NCBI's
/// thread-local `forbidden_tld`), the scans' memory, the score rows of the window's query letters
/// and the X-drop traceback's memory (NCBI's `gap_align`).
pub(crate) struct SwMatchScratch {
    forbidden: ForbiddenRanges,
    scan: SwScratch,
    rows: QueryPssm,
    gap: PssmGapScratch,
}

impl SwMatchScratch {
    pub(crate) fn new() -> Self {
        Self {
            forbidden: ForbiddenRanges::new(),
            scan: SwScratch::new(),
            rows: QueryPssm { rows: Vec::new() },
            gap: PssmGapScratch::new(),
        }
    }
}

impl Default for SwMatchScratch {
    fn default() -> Self {
        Self::new()
    }
}

/// NCBI's failed Smith-Waterman redo of a match: a significant alignment of score 0 whose X-drop
/// traceback has no operation (`GapEditScriptNew(0)` is NULL, `s_NewAlignmentUsingXdrop` returns
/// -1). NCBI frees the match's alignments and keeps the status -1 (`*pStatusCode`) until the next
/// match that it redoes sets it again; a status left by the last match of a query batch fails the
/// traceback (`Traceback failed`, traceback_stage.cpp:251-252).
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

/// NCBI's Smith-Waterman redo of a match that does not end: the score of 0 is significant, the
/// window has two or more HSPs and the first query letter's row scores the subject's first letter
/// above 0 (the same regime as RPS-BLAST's `SwEndless`, PD-LOSAT-NCBI-DEFECTS 1.29).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SwRedoEndless;

impl std::fmt::Display for SwRedoEndless {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("the Smith-Waterman redo of a match does not end")
    }
}

impl std::error::Error for SwRedoEndless {}

/// `Blast_RedoOneMatchSmithWaterman` with a square matrix for the protein windows of BLASTP
/// (`s_WindowsFromProteinAligns`): each window's alignments in NCBI's list order (the last one
/// made first) in `alignments_by_query[query_index]`, as the X-drop redo returns them, with the
/// composition p-value and lambda ratio of the last adjustment. `would_insert(query_index, evalue,
/// score)` is `BlastCompo_HeapWouldInsert(&significantMatches[query_index], ...)` with the match's
/// subject index; `cutoff_e` is `-evalue` (`params->cutoff_e`). The errors `SwRedoFailed` and
/// `SwRedoEndless` are NCBI's failed and endless redos (see there).
///
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
/// The subject range's last argument is `NULL`: SEG runs unless the window's first alignment is
/// near identical, as for the X-drop redo's first range (`subject_maybe_biased` starts true there
/// and this redo does not read its answer).
#[allow(clippy::too_many_arguments)]
pub(crate) fn redo_one_match_smith_waterman<'seq>(
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    params: &BlastRedoAlignParams,
    matching_seq: &BlastCompoMatchingSequence<'seq>,
    query_infos: &[BlastCompoQueryInfo],
    lambda: f64,
    callbacks: &BlastRedoAlignCallbacks,
    composition_workspace: &mut BlastCompositionWorkspace,
    matrix_state: &mut Option<AdjustedProteinMatrix>,
    cutoff_e: f64,
    would_insert: &mut dyn FnMut(usize, f64, i32) -> bool,
    scratch: &mut SwMatchScratch,
) -> Result<BlastRedoOneMatchResult> {
    // The translated windows (BLASTX, TBLASTN) and a PSSM's rows (psiblast) do not reach this
    // redo: their options refuse `-use_sw_tback` (item X4).
    if params.query_is_translated || params.subject_is_translated || params.position_based {
        bail!("the Smith-Waterman redo of translated or position-based windows is not ported");
    }
    // `-comp_based_stats 0` with `-use_sw_tback` is refused by BLASTP's options (item X4).
    if !params.uses_composition_based_stats() {
        bail!("the Smith-Waterman redo without composition-based statistics is not ported");
    }
    let query_info = query_infos
        .first()
        .ok_or_else(|| anyhow::anyhow!("BLASTP redo requires query information"))?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1352-1357
    // ```c
    //     status =
    //         s_WindowsFromAligns(incoming_aligns, query_info, hspcnt, numQueries,
    //                             kWindowBorder, matchingSeq->length, &windows,
    //                             &nWindows, query_is_translated, subject_is_translated, params->positionBased);
    // ```
    let windows = windows_from_protein_aligns(incoming_aligns, query_info, matching_seq.length)?;
    let mut alignments_by_query: Vec<Option<Box<BlastCompoAlignment>>> = Vec::new();
    let mut pvalue_for_this_pair = None;
    let mut lambda_ratio = None;
    let mut adjusted_matrix = matrix_state.take();
    let gap_open = params.gapping_params.gap_open;
    let gap_extend = params.gapping_params.gap_extend;

    for window in &windows {
        let query_index = usize::try_from(window.query_range.context)?;
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1439-1444
        // ```c++
        //     if (query_composition->numTrueAminoAcids == 0 ||
        //         subject_composition->numTrueAminoAcids == 0) {
        //         /* Either the query or subject contains only ambiguity
        //            characters, most likely because the entire subject has been
        //            SEGed.  Compositional adjustment is meaningless. */
        //         return 1;
        // ```
        let subject_composition = x_whole_subject_composition(&range);
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
        let matrix_adjust_rule = adjusted.matrix_adjust_rule;
        pvalue_for_this_pair = adjusted.pvalue_for_this_pair;
        lambda_ratio = Some(adjusted.lambda_ratio);
        let matrix = &*adjusted_matrix.insert(adjusted.adjusted_matrix);

        // NCBI's `matrix[query[queryPos]]` for each position of the window's query: the rows that
        // the kernels and `ALIGN_EX` read (`positionBased` false, blast_gapalign.c:543-554).
        let query = range.query.data();
        scratch.rows.rows.clear();
        scratch.rows.rows.extend(
            query
                .iter()
                .map(|&letter| matrix.scores[usize::from(letter)]),
        );
        let rows = &scratch.rows.rows;
        let letters = range.subject.data();
        // NCBI's `subject->data - 1`, which `ALIGN_EX` reads from (`s_SWFindFinalEndsUsingXdrop`):
        // the sentinel before the window's letters (`s_SequenceGetProteinRange`).
        let subject_from_sentinel = &range.subject.buffer[range.subject.data_offset - 1..];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3261-3264
        // ```c
        //         if(smithWaterman) {
        // 	   forbidden_tld[i] = calloc(1,sizeof(Blast_ForbiddenRanges));
        //            status_code =
        //                Blast_ForbiddenRangesInitialize(forbidden_tld[i], queryInfo->max_length);
        // ```
        // Here the window's query length, which every range it pushes is inside of.
        scratch.forbidden.initialize(query.len());

        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1431-1525
        // ```c
        //             Blast_ForbiddenRangesClear(forbidden);
        //             do {
        //                 ...
        //                 status =
        //                     Blast_SmithWatermanScoreOnly(&aSwScore, &matchEnd,
        //                                                  &queryEnd,
        //                                                  subject.data,
        //                                                  subject.length,
        //                                                  query.data,
        //                                                  query.length, matrix,
        //                                                  gap_open, gap_extend,
        //                                                  positionBased,
        //                                                  forbidden);
        //                 ...
        //                 if (do_link_hsps) {
        //                     alignment_is_significant = aSwScore >= params->cutoff_s;
        //                 } else {
        //                     double newSwEvalue;     /* evalue as computed by the
        //                                              * Smith-Waterman algorithm */
        //                     newSwEvalue =
        //                         s_EvalueFromScore(aSwScore, Lambda, logK, searchsp);
        //
        //                     alignment_is_significant = newSwEvalue < params->cutoff_e;
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
        //                 }
        //                 if (alignment_is_significant) {
        //                     ...
        //                     status =
        //                         Blast_SmithWatermanFindStart(&updatedScore,
        //                                                      &matchStart,
        //                                                      &queryStart,
        //                                                      subject.data,
        //                                                      subject.length,
        //                                                      query.data,
        //                                                      matrix, gap_open,
        //                                                      gap_extend,
        //                                                      matchEnd,
        //                                                      queryEnd,
        //                                                      aSwScore,
        //                                                      positionBased,
        //                                                      forbidden);
        //                     ...
        //                     status =
        //                         callbacks->
        //                         new_xdrop_align(&newAlign, &queryEnd, &matchEnd,
        //                                         queryStart, matchStart, aSwScore,
        //                                         &query, &window->query_range,
        //                                         ccat_query_length,
        //                                         &subject, &window->subject_range,
        //                                         matchingSeq->length,
        //                                         gapping_params, matrix_adjust_rule);
        //                     if (status != 0) {
        //                         goto window_index_loop_cleanup;
        //                     }
        //                     newAlign->next = alignments[query_index];
        //                     alignments[query_index] = newAlign;
        //
        //                     if (window->hspcnt > 1) {
        //                         /* We may compute more alignments; make the range
        //                            of the current alignment forbidden */
        //                         status =
        //                             Blast_ForbiddenRangesPush(forbidden,
        //                                                       queryStart, queryEnd,
        //                                                       matchStart, matchEnd);
        //                     }
        //                     ...
        //                 }
        //                 /* end if the next local alignment is significant */
        //             } while (alignment_is_significant && window->hspcnt > 1);
        // ```
        // `Lambda` and `logK` are the redo's scaled block (`kbp_gap` of NCBI's copy), `searchsp`
        // the query's `eff_search_space`, `cutoff_s` the redo's cutoff score.
        scratch.forbidden.clear();
        let mut first_round = true;
        loop {
            let end = smith_waterman_score_only(
                letters,
                query.len(),
                |pos| &rows[pos],
                gap_open,
                gap_extend,
                &scratch.forbidden,
                &mut scratch.scan,
            );
            // Whether this round is significant, and whether a round of score 0 would be.
            let (alignment_is_significant, zero_is_significant) = if params.do_link_hsps {
                (end.score >= params.cutoff_score, 0 >= params.cutoff_score)
            } else {
                let evalue = evalue_from_score(end.score, lambda, params.log_k, searchsp);
                let mut significant = evalue < cutoff_e;
                if alignments_by_query[query_index].is_none() {
                    significant = significant && would_insert(query_index, evalue, end.score);
                }
                (
                    significant,
                    evalue_from_score(0, lambda, params.log_k, searchsp) < cutoff_e,
                )
            };
            if !alignment_is_significant {
                break;
            }
            // The regime of a significant score 0, decided at the window's first round as
            // RPS-BLAST's redo decides it (`RedoMatch::redo_one_match_smith_waterman`, whose
            // comment gives the proof): with two or more HSPs every later round is significant (a
            // later score is at least 0, and the heap is asked only before the query's first
            // alignment), every round of a positive score forbids its free start cell, so the loop
            // comes to a round of score 0, whose find-start stops at (0, 0) and whose X-drop
            // traceback aligns the 1 x 1 rectangle there: the first query letter's row against the
            // subject's first letter after SEG. Above 0 that alignment is made again and again
            // without end (`SwRedoEndless`); otherwise it has no operation and the redo fails
            // (`SwRedoFailed`). With one HSP the loop makes one round, whose score 0 (all cells 0
            // or less) also fails, below.
            if first_round && window.hspcnt > 1 && zero_is_significant {
                // `ALIGN_EX`'s cell (1, 1): row 0 against `B[1]`, the subject's first letter (the
                // sentinel when it has none).
                let first_letter = subject_from_sentinel.get(1).copied().unwrap_or(0);
                let first_cell = rows.first().map_or(0, |row| row[row_column(first_letter)]);
                return Err(if first_cell > 0 {
                    SwRedoEndless.into()
                } else {
                    SwRedoFailed.into()
                });
            }
            let start = smith_waterman_find_start(
                letters,
                |pos| &rows[pos],
                gap_open,
                gap_extend,
                end.match_end,
                end.query_end,
                end.score,
                &scratch.forbidden,
                &mut scratch.scan,
            );
            let traced = pssm_sw_xdrop_traceback(
                &scratch.rows,
                subject_from_sentinel,
                (start.query_start, end.query_end),
                (start.match_start, end.match_end),
                gap_open,
                gap_extend,
                params.gapping_params.x_dropoff,
                end.score,
                &mut scratch.gap,
            );
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
            //         ...
            //     }
            //     *pnewAlign = obj;
            //
            //     return obj != NULL ? 0 : -1;
            // ```
            // An X-drop alignment without operations has no edit script (`GapEditScriptNew(0)`
            // is NULL): the redo of the match fails.
            if traced.edit_script.is_empty() {
                return Err(SwRedoFailed.into());
            }
            let query_start = i32::try_from(start.query_start)?;
            let match_start = i32::try_from(start.match_start)?;
            let query_end = query_start + traced.query_extent;
            let match_end = match_start + traced.match_extent;
            let mut new_align = blast_compo_alignment_new(
                traced.score,
                matrix_adjust_rule,
                query_start + window.query_range.begin,
                query_end + window.query_range.begin,
                window.query_range.context,
                match_start + window.subject_range.begin,
                match_end + window.subject_range.begin,
                window.subject_range.context,
                Some(BlastCompoAlignmentContext::EditScript(traced.edit_script)),
            );
            new_align.next = alignments_by_query[query_index].take();
            alignments_by_query[query_index] = Some(new_align);
            if window.hspcnt <= 1 {
                break;
            }
            // A significant score of 0 is decided above, so this round's score is positive and
            // its range forbids its free start cell: the loop ends.
            debug_assert!(end.score > 0);
            scratch
                .forbidden
                .push(query_start, query_end, match_start, match_end);
            first_round = false;
        }
    }
    *matrix_state = adjusted_matrix;
    Ok(BlastRedoOneMatchResult {
        alignments_by_query,
        pvalue_for_this_pair,
        lambda_ratio,
    })
}

#[cfg(test)]
#[path = "smith_waterman_tests.rs"]
mod tests;
