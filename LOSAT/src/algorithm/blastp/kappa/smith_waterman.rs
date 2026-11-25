//! BLASTP's Smith-Waterman redo of one match (`-use_sw_tback`; series RG, item X4): core's match
//! body with a square matrix (`core::composition_adjustment::smith_waterman_match`, NCBI's
//! `Blast_RedoOneMatchSmithWaterman`), called with BLASTP's protein windows, its scaled block and its
//! heaps. The steps after it (`s_HSPListFromDistinctAlignments`, the reap,
//! `s_HitlistEvaluateAndPurge`, the normalization and the identities) are the X-drop redo's
//! (blast_kappa.c:3654-3707): the caller's (`postprocess_preliminary_hits_smith_waterman`).

use anyhow::Result;

use crate::core::composition_adjustment::adjust_scores::{
    AdjustedProteinMatrix, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    BlastCompoAlignment, BlastCompoMatchingSequence, BlastCompoQueryInfo, BlastRedoAlignCallbacks,
    BlastRedoAlignParams, BlastRedoOneMatchResult,
};
// The match body, its memory and its two ends are core's, one owner for BLASTP's and TBLASTN's
// redo (`core::composition_adjustment::smith_waterman_match`), as are the answers of its heap test
// and the state from which a match whose test was left open resumes.
pub(crate) use crate::core::composition_adjustment::smith_waterman_match::{
    SwHeapAnswer, SwMatchResume, SwMatchScratch,
};

/// `Blast_RedoOneMatchSmithWaterman` with a square matrix for the protein windows of BLASTP
/// (`s_WindowsFromProteinAligns`): core's match body
/// (`smith_waterman_match::redo_one_match_smith_waterman_from`) with BLASTP's scaled block
/// (`lambda`, `params.log_k`, `kbp_gap[0]` of the redo's copy) and BLASTP's heaps. Without
/// composition-based statistics (`-comp_based_stats 0`, item X4-c0) `unadjusted` holds the
/// unscaled search matrix and no window is adjusted; with modes 1 to 3 `unadjusted` is `None`.
/// `would_insert` answers `SwHeapAnswer` (the one-thread loop always `Known`; the parallel redo's
/// jobs leave a test open with `Deferred`, which keeps the test's state in `kept`), and `resume`
/// is such a state, from which the call finishes the match.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3602-3613
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
    unadjusted: Option<&AdjustedProteinMatrix>,
    cutoff_e: f64,
    resume: Option<&SwMatchResume>,
    would_insert: &mut dyn FnMut(usize, f64, i32) -> SwHeapAnswer,
    kept: &mut Option<Box<SwMatchResume>>,
    scratch: &mut SwMatchScratch,
) -> Result<BlastRedoOneMatchResult> {
    crate::core::composition_adjustment::smith_waterman_match::redo_one_match_smith_waterman_from(
        incoming_aligns,
        params,
        lambda,
        params.log_k,
        matching_seq,
        query_infos,
        callbacks,
        composition_workspace,
        matrix_state,
        unadjusted,
        cutoff_e,
        resume,
        &mut |query_index, evalue, score| Ok(would_insert(query_index, evalue, score)),
        kept,
        scratch,
    )
}

#[cfg(test)]
#[path = "smith_waterman_tests.rs"]
mod tests;
