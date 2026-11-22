//! The Smith-Waterman redo of one match with a square matrix (`-use_sw_tback`; series RG, item
//! X4): NCBI's `Blast_RedoOneMatchSmithWaterman` for BLASTP's protein window, in place of the
//! X-drop redo (`Blast_RedoOneMatch`, `blast_redo_one_match_with_workspace_queries_and_matrix`).
//! The window, the subject range (`get_range`), the composition adjustment (`Blast_AdjustScores`)
//! and the alignment list it returns are the X-drop redo's; the round loop of a window (the
//! score-only scan, the significance test, the regime of a significant score of 0, the start scan,
//! the X-drop traceback from the Smith-Waterman start and the forbidden ranges;
//! `core::composition_adjustment::smith_waterman::smith_waterman_redo_window`) is the one RPS-BLAST's
//! redo runs (line PS, PS11c-3), on the rows that the square matrix gives the query's letters
//! (`matrix[query[queryPos]]`, which `ALIGN_EX` and the kernels read as they read a PSSM's
//! `pssm[queryPos]`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3602-3621
//! The steps after it (`s_HSPListFromDistinctAlignments`, the reap, `s_HitlistEvaluateAndPurge`,
//! the normalization and the identities) are the X-drop redo's (blast_kappa.c:3654-3707): the
//! caller's (`postprocess_preliminary_hits_smith_waterman`).

use anyhow::{bail, Result};

use crate::algorithm::pssm_search::gapalign::{pssm_sw_xdrop_traceback, PssmGapScratch};
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
    smith_waterman_redo_window, SwRedoScratch, SwSignificance, SwWindow, SwWindowRedo,
};

/// BLASTP's memory of its windows (`SwRedoScratch::window`): the score rows of the window's query
/// letters in the adjusted matrix (`matrix[query[queryPos]]`) and the X-drop traceback's memory
/// (NCBI's `gap_align`).
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

/// BLASTP's memory of the Smith-Waterman redo from one match to the next: the round loop's
/// (`SwRedoScratch::rounds`, shared with RPS-BLAST's redo) and its windows' (`SquareWindow`).
pub(crate) type SwMatchScratch = SwRedoScratch<SquareWindow>;

/// NCBI's failed Smith-Waterman redo of a match: a significant alignment of score 0 whose X-drop
/// traceback has no operation (`GapEditScriptNew` of 0 is NULL, `s_NewAlignmentUsingXdrop` returns
/// -1). NCBI frees the match's alignments and keeps the status -1 (`*pStatusCode`) until the next
/// match that it redoes sets it again; a status left by the last match of a query batch fails the
/// traceback (`Traceback failed`, traceback_stage.cpp:251-252), or crashes NCBI's clean-up when
/// the batch's first query holds a redone match (`blast_engine::redo_run::sw_batch_failed`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1542-1548
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
/// above 0 (`SwWindowRedo::Endless`, the regime of RPS-BLAST's `SwEndless`; PD-LOSAT-NCBI-DEFECTS
/// 1.53, within the class of 1.30).
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
/// score)` is `BlastCompo_HeapWouldInsert` on the heap `significantMatches[query_index]` with the match's
/// subject index; `cutoff_e` is `-evalue` (`params->cutoff_e`). The errors `SwRedoFailed` and
/// `SwRedoEndless` are NCBI's failed and endless redos (see there).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1366-1430
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
        let SquareWindow { rows, gap } = &mut scratch.window;
        rows.rows.clear();
        rows.rows.extend(
            query
                .iter()
                .map(|&letter| matrix.scores[usize::from(letter)]),
        );
        let rows = &*rows;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3261-3264
        // The round loop (`smith_waterman_redo_window`, redo_alignment.c:1438-1529) on the
        // subject's range from NCBI's position one before `subject->data`, which `ALIGN_EX` reads from
        // (`s_SWFindFinalEndsUsingXdrop`): the sentinel before the window's letters
        // (`s_SequenceGetProteinRange`). `Lambda` and `logK` are the redo's scaled block (`kbp_gap`
        // of NCBI's copy), `searchsp` the query's `eff_search_space`, `cutoff_s` the redo's cutoff
        // score.
        let sw_window = SwWindow {
            subject_from_sentinel: &range.subject.buffer[range.subject.data_offset - 1..],
            subject_length: range.subject.data().len(),
            query_length: query.len(),
            rows: |pos| &rows.rows[pos],
            gap_open,
            gap_extend,
            // NCBI's `int` count, never negative (`s_WindowsFromProteinAligns` counts HSPs up).
            hspcnt: usize::try_from(window.hspcnt).expect("a window's HSP count"),
        };
        // BLASTP's protein windows are one per query (`s_WindowsFromProteinAligns`), so the
        // query's alignments (`alignments[query_index]`, whose emptiness the heap test reads) are
        // this window's.
        let significance = if params.do_link_hsps {
            SwSignificance::CutoffScore(params.cutoff_score)
        } else {
            SwSignificance::Evalue {
                lambda,
                log_k: params.log_k,
                searchsp,
                cutoff_e,
            }
        };
        let x_dropoff = params.gapping_params.x_dropoff;
        let redone = smith_waterman_redo_window(
            &sw_window,
            significance,
            |evalue, score| would_insert(query_index, evalue, score),
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
        let made = match redone {
            SwWindowRedo::Alignments(made) => made,
            SwWindowRedo::Failed => return Err(SwRedoFailed.into()),
            SwWindowRedo::Endless => return Err(SwRedoEndless.into()),
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1856-1873
        // Each alignment prepended in the order made (`newAlign->next` is set to `alignments[query_index]`).
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

#[cfg(test)]
#[path = "smith_waterman_tests.rs"]
mod tests;
