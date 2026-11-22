//! The redo of one TBLASTN match in stage D: NCBI's `Blast_RedoOneMatch` (the X-drop redo of
//! composition-based statistics, `kappa::redo_preliminary_match`) or, with `-use_sw_tback`,
//! `Blast_RedoOneMatchSmithWaterman` (core's match body,
//! `core::composition_adjustment::smith_waterman_match`; series RG, item X4b-1). The search
//! chooses one of them once per query batch (`stage_d_pipeline::run_local_search_with_pool`), and
//! the steps around the redo (early termination, the reap, the linking or the e-values, the heaps)
//! are the same for both (`stage_d_pipeline::run_local_search_with_redo`).
//!
//! NCBI's traceback runs the composition-based redo when composition-based statistics are on or
//! the traceback option is Smith-Waterman, and the redo takes the Smith-Waterman form from that
//! option alone:
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1498
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039

use anyhow::Result;

use super::kappa::{incoming_alignments, redo_preliminary_match, TBLASTN_REDO_CALLBACKS};
use super::kappa_heap::{CompoHeap, CompoHeapRecord};
use super::search_gapped::GappedHsp;
use crate::algorithm::blastp::gapalign::GapAlignScratch;
use crate::blastinput::app::ErrorAfterBatchReports;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::{
    AdjustedProteinMatrix, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    BlastCompoAdjustMode, BlastCompoMatchingSequence, BlastCompoQueryInfo, BlastRedoAlignParams,
    BlastRedoOneMatchResult,
};
use crate::core::composition_adjustment::smith_waterman_match::{
    redo_one_match_smith_waterman, sw_batch_failed_error, sw_endless_error,
    unadjusted_square_matrix, SwBatchStatus, SwMatchScratch,
};
use crate::sequence::packed_subject::SubjectNt;

/// The redo of one match (`Blast_RedoOneMatch` or `Blast_RedoOneMatchSmithWaterman`), chosen once
/// per query batch: the batch's loop (`stage_d_pipeline::run_local_search_with_redo`) is generic
/// over it, so that each body is its own function and the X-drop redo keeps its code.
pub(super) trait MatchRedo: Sync {
    /// NCBI's `smithWaterman`: the traceback option is Smith-Waterman (`eSmithWatermanTbck`).
    const SMITH_WATERMAN: bool;

    /// The memory one thread keeps from one match to the next.
    type Scratch: Send;

    /// A thread's memory for the redo.
    fn scratch(&self) -> Self::Scratch;

    /// Whether a match's redo reads the heaps of the redone matches (the Smith-Waterman
    /// significance test of e-value mode): the matches are then redone in stream order on one
    /// thread, as NCBI's `-num_threads 1` redoes them.
    fn reads_heaps(&self, params: &BlastRedoAlignParams) -> bool;

    /// The redo of the match of `preliminary_hits` (query `query_context`, subject `oid`), its
    /// alignments by query. `lambda` and `log_k` are the redo's scaled block of the query
    /// (`kbp_gap[context_index]`); `heaps` are the batch's heaps of redone matches where the redo
    /// reads them (`reads_heaps`).
    #[allow(clippy::too_many_arguments)]
    fn redo(
        &self,
        preliminary_hits: &[GappedHsp],
        query_context: i32,
        query_infos: &[BlastCompoQueryInfo],
        subject: SubjectNt<'_>,
        genetic_code: u8,
        params: &BlastRedoAlignParams,
        lambda: f64,
        log_k: f64,
        matrix: ScoringMatrix,
        oid: usize,
        heaps: Option<&mut [CompoHeap]>,
        scratch: &mut Self::Scratch,
        workspace: &mut BlastCompositionWorkspace,
        matrix_state: &mut Option<AdjustedProteinMatrix>,
    ) -> Result<BlastRedoOneMatchResult>;
}

/// NCBI's X-drop redo of composition-based statistics (`Blast_RedoOneMatch`): the search without
/// `-use_sw_tback`, whose batches with `-comp_based_stats 0` never reach the redo.
pub(super) struct XDropRedo;

impl MatchRedo for XDropRedo {
    const SMITH_WATERMAN: bool = false;
    type Scratch = GapAlignScratch;

    #[inline(always)]
    fn scratch(&self) -> GapAlignScratch {
        GapAlignScratch::new()
    }

    #[inline(always)]
    fn reads_heaps(&self, _params: &BlastRedoAlignParams) -> bool {
        false
    }

    #[inline(always)]
    fn redo(
        &self,
        preliminary_hits: &[GappedHsp],
        query_context: i32,
        query_infos: &[BlastCompoQueryInfo],
        subject: SubjectNt<'_>,
        genetic_code: u8,
        params: &BlastRedoAlignParams,
        lambda: f64,
        _log_k: f64,
        matrix: ScoringMatrix,
        _oid: usize,
        _heaps: Option<&mut [CompoHeap]>,
        scratch: &mut GapAlignScratch,
        workspace: &mut BlastCompositionWorkspace,
        matrix_state: &mut Option<AdjustedProteinMatrix>,
    ) -> Result<BlastRedoOneMatchResult> {
        redo_preliminary_match(
            preliminary_hits,
            query_context,
            query_infos,
            subject,
            genetic_code,
            params,
            lambda,
            matrix,
            scratch,
            workspace,
            matrix_state,
        )
    }
}

/// NCBI's Smith-Waterman redo (`-use_sw_tback`): every window of the match's subject frames
/// (`s_WindowsFromTranslatedAligns`) with the round loop of core's match body, with
/// composition-based statistics 0 to 3 (0: the unadjusted matrix, `unadjusted`), in cutoff-score
/// mode (sum statistics: the redo's `cutoff_s`) or e-value mode (`-evalue`, `cutoff_e`, and the
/// heap of the match's query before its first alignment of the call).
pub(super) struct SmithWatermanRedo {
    /// The score block's matrix, which the redo's windows read without composition-based
    /// statistics (`unadjusted_square_matrix`); `None` with them.
    unadjusted: Option<AdjustedProteinMatrix>,
    /// `-evalue` (`params->cutoff_e`, the hit-saving parameters' `expect_value`).
    cutoff_e: f64,
}

impl SmithWatermanRedo {
    pub(super) fn new(matrix: ScoringMatrix, mode: BlastCompoAdjustMode, evalue: f64) -> Self {
        Self {
            unadjusted: (!mode.is_enabled()).then(|| unadjusted_square_matrix(matrix)),
            cutoff_e: evalue,
        }
    }
}

impl MatchRedo for SmithWatermanRedo {
    const SMITH_WATERMAN: bool = true;
    type Scratch = SwMatchScratch;

    fn scratch(&self) -> SwMatchScratch {
        SwMatchScratch::new()
    }

    // The heap test of e-value mode; a cutoff score (sum statistics) reads no heap.
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1456-1473
    fn reads_heaps(&self, params: &BlastRedoAlignParams) -> bool {
        !params.do_link_hsps
    }

    fn redo(
        &self,
        preliminary_hits: &[GappedHsp],
        query_context: i32,
        query_infos: &[BlastCompoQueryInfo],
        subject: SubjectNt<'_>,
        genetic_code: u8,
        params: &BlastRedoAlignParams,
        lambda: f64,
        log_k: f64,
        _matrix: ScoringMatrix,
        oid: usize,
        mut heaps: Option<&mut [CompoHeap]>,
        scratch: &mut SwMatchScratch,
        workspace: &mut BlastCompositionWorkspace,
        matrix_state: &mut Option<AdjustedProteinMatrix>,
    ) -> Result<BlastRedoOneMatchResult> {
        let incoming = incoming_alignments(preliminary_hits, query_context, params.score_divisor);
        let matching_seq = BlastCompoMatchingSequence::new_translated(0, subject, genetic_code)?;
        // NCBI's `matchingSeq->index`, the subject's OID (`s_MatchingSequenceInitialize`).
        let subject_index = i32::try_from(oid)?;
        redo_one_match_smith_waterman(
            &incoming,
            params,
            lambda,
            log_k,
            &matching_seq,
            query_infos,
            // TBLASTN's callbacks: the subject window's translation (`s_SequenceGetTranslatedRange`,
            // with SEG under composition-based statistics); the body never redoes one alignment.
            &TBLASTN_REDO_CALLBACKS,
            workspace,
            matrix_state,
            self.unadjusted.as_ref(),
            self.cutoff_e,
            &mut |query_index, evalue, score| {
                let heaps = heaps.as_deref_mut().ok_or_else(|| {
                    anyhow::anyhow!("the Smith-Waterman redo of e-value mode needs the heaps")
                })?;
                let heap = heaps.get_mut(query_index).ok_or_else(|| {
                    anyhow::anyhow!("the Smith-Waterman redo's query {query_index} has no heap")
                })?;
                Ok(heap.would_insert(CompoHeapRecord {
                    best_evalue: evalue,
                    best_score: score,
                    subject_index,
                }))
            },
            scratch,
        )
    }
}

/// The result of one redone match in the batch's loop, and NCBI's status of the last match the
/// redo reached (`status`, core's `SwBatchStatus`, which owns the rule): a Smith-Waterman redo
/// that fails (`SwRedoFailed`) leaves the match without alignments and sets the status, the next
/// match that the redo reaches clears it (a match that early termination skips leaves it; a
/// list whose subject the fetch excludes clears it too, `SwBatchStatus::excluded_at_fetch`); a
/// redo that does not end in NCBI (`SwRedoEndless`) is LOSAT's refusal of the batch
/// (`sw_endless`). The X-drop redo's result as it is.
#[inline(always)]
pub(super) fn match_outcome<R: MatchRedo, T>(
    result: Result<Option<T>>,
    status: &mut SwBatchStatus,
    do_link_hsps: bool,
) -> Result<Option<T>> {
    if !R::SMITH_WATERMAN {
        return result;
    }
    Ok(status
        .redone(result, || sw_endless(do_link_hsps))?
        .flatten())
}

/// The error of a query batch whose last redone match failed its Smith-Waterman redo: NCBI's
/// `Traceback failed` (exit 3), or LOSAT's refusal where NCBI's clean-up crashes, when the
/// batch's first query holds a redone match (`first_query_holds_match`, the heap of the batch's
/// first query is not empty; PD-LOSAT-NCBI-DEFECTS 1.72, rule 3), after the reports of the
/// batches before (`ErrorAfterBatchReports`, `tblastn_app::batch_failure`). NCBI's `-num_threads
/// 1` result at every thread count: with `-db` and more threads NCBI interrupts the batch at the
/// first failure.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3775-3779
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3795-3802
/// (cited at `smith_waterman_match::sw_batch_failed_error`: `s_ClearHeap` pops the first query's
/// heap after the arrays are freed)
#[cold]
#[inline(never)]
pub(super) fn sw_batch_failed(first_query_holds_match: bool) -> anyhow::Error {
    ErrorAfterBatchReports(sw_batch_failed_error("tblastn", first_query_holds_match)).into()
}

/// LOSAT's refusal of a PSSM query's search (`-in_pssm`, the one query batch) whose
/// Smith-Waterman redo reaches a match (a preliminary hit list that the fetch keeps), after the
/// prolog (`ErrorAfterBatchReports`). NCBI's redo of a PSSM query (`positionBased`, the PSSM's rows
/// as its matrix) takes the PSSM's rows by the positions of the window, not of the query, and
/// reads memory it never wrote: it makes alignments whose scores pass 2^30 and frame 0 (survey
/// probe `gdb3.tbn_inpssm_s`), then, by the input, fails the traceback (`Traceback failed`, exit
/// 3), stops with a C++ exception of the Seq-align mapping (exit 255), crashes (`free(): invalid
/// pointer` or `invalid size`, exit 134, or SIGSEGV, with `-subject` and with `-db`), does not end
/// (still running at 120 s), or ends with exit 0 and rows that differ from run to run (NCBI
/// BLAST+ 2.17.0, probes `pssm.*` of X4b-1 fix round 1 and probes 1 and 3 of its review 2: 20
/// PSSM and subject pairs, three runs of each pair that ends). None of these gives a result to
/// match, so the search is refused (rule 3). A search without such a match ends as NCBI's (exit
/// 0). Not ported: the Smith-Waterman redo of a PSSM (`smith_waterman_match` takes a square
/// matrix).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3359-3360
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:86-90
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:866-873
/// (`-use_sw_tback` sets the mode after `SetPSITblastnDefaults`, blast_args.cpp:885)
#[cold]
#[inline(never)]
pub(super) fn sw_pssm_refusal() -> anyhow::Error {
    ErrorAfterBatchReports(anyhow::anyhow!(
        "-use_sw_tback with -in_pssm on a search with hits, on which NCBI BLAST+ 2.17.0's tblastn fails, crashes, does not end or writes rows that differ from run to run, is not supported by LOSAT's TBLASTN"
    ))
    .into()
}

/// LOSAT's refusal of a batch with a match whose Smith-Waterman redo does not end in NCBI (a
/// significant score of 0 in a window of two or more preliminary HSPs whose first query letter
/// scores the window's first subject letter above 0; PD-LOSAT-NCBI-DEFECTS 1.72, rule 3), after
/// the reports of the batches before (`ErrorAfterBatchReports`): the Owner's decision of
/// 2026-11-16 on line RG's Q-X2c, "the results of the query batches before are written". What
/// makes the score of 0 significant: a cutoff score of 0 with sum statistics (a tiny `-dbsize`),
/// or its e-value in e-value mode.
#[cold]
#[inline(never)]
fn sw_endless(do_link_hsps: bool) -> anyhow::Error {
    let why = if do_link_hsps {
        "a score of 0 reaches the cutoff score, a window of the match has several preliminary HSPs and the first query letter scores the window's first subject letter above 0"
    } else {
        "the e-value of a score of 0 passes the threshold, a window of the match has several preliminary HSPs and the first query letter scores the window's first subject letter above 0"
    };
    ErrorAfterBatchReports(sw_endless_error("tblastn", why)).into()
}
