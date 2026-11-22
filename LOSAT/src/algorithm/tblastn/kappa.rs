//! TBLASTN translated-subject Kappa traceback boundary.
//! The public TBLASTN CLI remains gated until stages D and E agree.

use anyhow::{bail, Context, Result};
use std::ffi::c_void;
use std::ptr::NonNull;

use super::search_gapped::{GappedHsp, TargetTranslation};
use super::stage_d_linking::score_compare;
use crate::algorithm::blastp::gapalign::{
    blast_gapped_alignment_with_traceback_with_scratch, protein_identities_from_edit_ops,
    stats_from_edit_ops_protein, GapAlignScratch,
};
use crate::common::GapEditOp;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::{
    AdjustedProteinMatrix, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    blast_compo_alignment_new, blast_redo_one_match_with_workspace_queries_and_matrix,
    redo_calc_lambda, translated_subject_get_range, BlastCompoAdjustMode, BlastCompoAlignment,
    BlastCompoAlignmentContext, BlastCompoMatchingSequence, BlastCompoQueryInfo,
    BlastCompoSequenceData, BlastCompoSequenceRange, BlastRedoAlignCallbacks, BlastRedoAlignParams,
    BlastRedoOneMatchResult, EMatrixAdjustRule,
};

mod x_reap_index;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1738-1773,1896-1957
#[allow(dead_code)] // Called by the TBLASTN D pipeline after the preliminary HSP callback is wired.
pub(super) fn redo_one_alignment_from_local_starts(
    q_start: usize,
    s_start: usize,
    query_data: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_data: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    matrix_adjust_rule: EMatrixAdjustRule,
    adjusted_matrix: Option<&AdjustedProteinMatrix>,
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    scratch: &mut GapAlignScratch,
) -> Result<Option<Box<BlastCompoAlignment>>> {
    let mut fence_hit = false;
    let Some(mut aligned) = blast_gapped_alignment_with_traceback_with_scratch(
        query_data.data(),
        subject_data.data(),
        q_start,
        s_start,
        matrix,
        adjusted_matrix,
        gap_open,
        gap_extend,
        x_dropoff,
        scratch,
        Some(&mut fence_hit),
    ) else {
        return Ok(None);
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1760
    let query_start = aligned
        .query_start
        .checked_add(query_range.begin)
        .context("TBLASTN redo query start overflow")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1761
    let query_end = aligned
        .query_stop
        .checked_add(query_range.begin)
        .context("TBLASTN redo query end overflow")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1763
    let match_start = aligned
        .subject_start
        .checked_add(subject_range.begin)
        .context("TBLASTN redo subject start overflow")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1764
    let match_end = aligned
        .subject_stop
        .checked_add(subject_range.begin)
        .context("TBLASTN redo subject end overflow")?;
    let context = Some(BlastCompoAlignmentContext::EditScript(std::mem::take(
        &mut aligned.edit_script,
    )));
    Ok(Some(blast_compo_alignment_new(
        aligned.score,
        matrix_adjust_rule,
        query_start,
        query_end,
        query_range.context,
        match_start,
        match_end,
        subject_range.context,
        context,
    )))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:224-273;
// c++/src/algo/blast/core/blast_hits_priv.h:65-70
#[allow(dead_code)] // Entered after TBLASTN Kappa alignment conversion is in the public D pipeline.
pub(super) fn reap_contained_postredo_hsps(hsps: &mut Vec<(usize, GappedHsp)>) -> Vec<usize> {
    if hsps.len() < 2 {
        return (0..hsps.len()).collect();
    }
    // No NCBI counterpart: LOSAT_X_TBNREAPIDX gives the same answers from an index by subject start
    // for a list of `X_REAP_INDEX_FROM` HSPs or more (`kappa/x_reap_index.rs`).
    let keep = if x_reap_index::x_indexed(hsps.len()) {
        x_reap_index::x_reap_keep(hsps)
    } else {
        reap_keep(hsps)
    };
    let retained_indices = keep
        .iter()
        .enumerate()
        .filter_map(|(index, &selected)| selected.then_some(index))
        .collect();
    let mut index = 0;
    hsps.retain(|_| {
        let selected = keep[index];
        index += 1;
        selected
    });
    retained_indices
}

/// Whether each HSP of `hsps` is kept: NCBI's loop of `s_HitlistReapContained` (above).
fn reap_keep(hsps: &[(usize, GappedHsp)]) -> Vec<bool> {
    let mut keep = vec![true; hsps.len()];
    for read in 1..hsps.len() {
        let (query_context, inner) = hsps[read];
        for previous in 0..read {
            if !keep[previous] {
                continue;
            }
            let (previous_context, outer) = hsps[previous];
            if query_context != previous_context || inner.frame != outer.frame {
                continue;
            }
            let start_contained = outer.q_start <= inner.q_start
                && outer.q_end >= inner.q_start
                && outer.s_start <= inner.s_start
                && outer.s_end >= inner.s_start;
            let end_contained = outer.q_start <= inner.q_end
                && outer.q_end >= inner.q_end
                && outer.s_start <= inner.s_end
                && outer.s_end >= inner.s_end;
            if start_contained && end_contained && inner.score <= outer.score {
                keep[read] = false;
                break;
            }
        }
    }
    keep
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:305-358:
#[derive(Debug)]
pub(super) struct ConvertedKappaHsp {
    pub context: usize,
    pub hsp: GappedHsp,
    pub edit_script: Vec<crate::common::GapEditOp>,
    pub matrix_adjust_rule: EMatrixAdjustRule,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:305-358:
pub(super) fn convert_distinct_alignments(
    alignments: &mut Option<Box<BlastCompoAlignment>>,
) -> Result<Vec<ConvertedKappaHsp>> {
    let mut converted = Vec::new();
    let mut current = alignments.take();
    while let Some(mut alignment) = current {
        current = alignment.next.take();
        let Some(BlastCompoAlignmentContext::EditScript(edit_script)) = alignment.context.take()
        else {
            bail!("TBLASTN Kappa alignment has no edit script");
        };
        converted.push(ConvertedKappaHsp {
            context: usize::try_from(alignment.query_index)?,
            hsp: GappedHsp {
                frame: i8::try_from(alignment.frame)?,
                score: alignment.score,
                q_start: alignment.query_start,
                q_end: alignment.query_end,
                q_gapped_start: 0,
                s_start: alignment.match_start,
                s_end: alignment.match_end,
                s_gapped_start: 0,
            },
            edit_script,
            matrix_adjust_rule: alignment.matrix_adjust_rule,
        });
    }
    converted.sort_by(|a, b| score_compare(&a.hsp, &b.hsp));
    Ok(converted)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:459-532;
// c++/src/algo/blast/core/blast_hits.c:1147-1228,767-811
#[allow(dead_code)] // Used by the TBLASTN D result pipeline after Kappa relink/reap.
pub(super) fn postredo_num_ident(
    alignment: &BlastCompoAlignment,
    query_sequence: &[u8],
    target: &mut TargetTranslation<'_>,
) -> Result<usize> {
    let Some(BlastCompoAlignmentContext::EditScript(edit_script)) = alignment.context.as_ref()
    else {
        bail!("TBLASTN Kappa identity requires an edit script");
    };
    let (subject, _, subject_base) = target.get(
        i8::try_from(alignment.frame)?,
        alignment.match_start,
        alignment.match_end,
    )?;
    let subject_start = super::search_gapped::window_index(alignment.match_start, subject_base)
        .context("TBLASTN identity subject offset before translation range")?;
    Ok(protein_identities_from_edit_ops(
        query_sequence,
        subject,
        usize::try_from(alignment.query_start)?,
        subject_start,
        edit_script,
        ScoringMatrix::Blosum62,
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:515-526:
// c++/src/algo/blast/core/blast_hits.c:767-811:
#[cfg(test)]
pub(super) fn postredo_converted_stats(
    converted: &ConvertedKappaHsp,
    query_sequence: &[u8],
    target: &mut TargetTranslation<'_>,
) -> Result<(usize, usize, usize, usize, usize, usize)> {
    postredo_converted_stats_with_matrix(converted, query_sequence, target, ScoringMatrix::Blosum62)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:767-811;
// blast_traceback.c:650-662; objtools/align_format/tabular.cpp:1003-1021:
// identity uses aligned residues; positives use the selected scoring matrix.
pub(super) fn postredo_converted_stats_with_matrix(
    converted: &ConvertedKappaHsp,
    query_sequence: &[u8],
    target: &mut TargetTranslation<'_>,
    matrix: ScoringMatrix,
) -> Result<(usize, usize, usize, usize, usize, usize)> {
    let hsp = &converted.hsp;
    let (subject, _, subject_base) = target.get(hsp.frame, hsp.s_start, hsp.s_end)?;
    let subject_start = super::search_gapped::window_index(hsp.s_start, subject_base)
        .context("TBLASTN positive subject offset before translation range")?;
    let (ident, positive, mismatch, gap_opens, gap_letters) = stats_from_edit_ops_protein(
        query_sequence,
        subject,
        usize::try_from(hsp.q_start)?,
        subject_start,
        &converted.edit_script,
        matrix,
    );
    let align_length = converted
        .edit_script
        .iter()
        .map(|op| op.num() as usize)
        .sum();
    Ok((
        ident,
        positive,
        align_length,
        mismatch,
        gap_opens,
        gap_letters,
    ))
}

/// NCBI's engine positives of an HSP of a redone match (`BlastHSP::num_positives`, which the
/// Seq-aligns show): `s_ComputeNumIdentities` counts them on the working (masked) query and
/// the target translation with the score block whose matrix the match's composition
/// adjustments have overwritten (`Blast_AdjustScores` writes its `matrix` argument, the
/// thread-local `sbp->matrix->data`), so a pair is positive by `adjusted`, the matrix of the
/// match's last adjustment. `postredo_converted_stats` counts them with BLOSUM62 for the other
/// reports. With the Seq-aligns only (`LocalStageDProfile::engine_positives`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3500-3502,3625-3635
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3360-3362
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1232-1241
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:525-528
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:800-812
pub(super) fn postredo_engine_positives(
    converted: &ConvertedKappaHsp,
    query_sequence: &[u8],
    target: &mut TargetTranslation<'_>,
    adjusted: &AdjustedProteinMatrix,
) -> Result<usize> {
    let hsp = &converted.hsp;
    let (subject, _, subject_base) = target.get(hsp.frame, hsp.s_start, hsp.s_end)?;
    let mut s_index = super::search_gapped::window_index(hsp.s_start, subject_base)
        .context("TBLASTN positive subject offset before translation range")?;
    let mut q_index = usize::try_from(hsp.q_start)?;
    let mut positives = 0usize;
    for op in &converted.edit_script {
        match *op {
            GapEditOp::Sub(n) => {
                for _ in 0..n {
                    let (q, s) = (query_sequence[q_index], subject[s_index]);
                    if q == s || adjusted.is_positive(q, s) {
                        positives += 1;
                    }
                    q_index += 1;
                    s_index += 1;
                }
            }
            GapEditOp::Del(n) => s_index += n as usize,
            GapEditOp::Ins(n) => q_index += n as usize,
        }
    }
    Ok(positives)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:102-114;
// c++/src/algo/blast/core/ncbi_math.c:437-441
#[allow(dead_code)] // Used by the TBLASTN D result pipeline after Kappa relink/reap.
pub(super) fn normalize_postredo_scores(
    list: &mut super::stage_d_linking::LinkedHspList,
    lambda: f64,
    log_k: f64,
    score_divisor: f64,
) -> Vec<f64> {
    list.hsps
        .iter_mut()
        .map(|linked| {
            let scaled = (linked.hsp.score as f64) / score_divisor;
            linked.hsp.score = (scaled + if scaled >= 0.0 { 0.5 } else { -0.5 }) as i32;
            (linked.hsp.score as f64 * lambda * score_divisor - log_k) / std::f64::consts::LN_2
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:770-803,1910-1948
struct TblastnRedoContext<'a> {
    preliminary_hits: &'a [GappedHsp],
    matrix: ScoringMatrix,
    scratch: NonNull<GapAlignScratch>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1896-1957
fn tblastn_redo_one_alignment_callback(
    incoming_align: &BlastCompoAlignment,
    rule: EMatrixAdjustRule,
    adjusted: Option<&AdjustedProteinMatrix>,
    query_data: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    _ccat_query_length: i32,
    subject_data: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    _full_subject_length: i32,
    params: &BlastRedoAlignParams,
) -> Result<Option<Box<BlastCompoAlignment>>> {
    let Some(BlastCompoAlignmentContext::PreliminaryHspIndex(index)) =
        incoming_align.context.as_ref()
    else {
        bail!("TBLASTN Kappa redo requires a preliminary HSP index");
    };
    let pointer = params
        .gapping_params
        .context
        .get()
        .context("TBLASTN Kappa redo context was not installed")?;
    // Safety: redo_preliminary_match installs this context for the synchronous
    // Blast_RedoOneMatch call and restores the prior pointer immediately after.
    let context = unsafe { &*pointer.cast::<TblastnRedoContext<'_>>().as_ptr() };
    let hit = context
        .preliminary_hits
        .get(*index)
        .context("TBLASTN Kappa preliminary HSP index is missing")?;
    let q_start = usize::try_from(hit.q_gapped_start - query_range.begin)?;
    let s_start = usize::try_from(hit.s_gapped_start - subject_range.begin)?;
    // Safety: the caller owns scratch for the full synchronous redo call.
    let scratch = unsafe { context.scratch.as_ptr().as_mut() }
        .context("TBLASTN Kappa traceback scratch pointer is invalid")?;
    redo_one_alignment_from_local_starts(
        q_start,
        s_start,
        query_data,
        query_range,
        subject_data,
        subject_range,
        rule,
        adjusted,
        context.matrix,
        params.gapping_params.gap_open,
        params.gapping_params.gap_extend,
        params.gapping_params.x_dropoff,
        scratch,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2389-2394
pub(super) const TBLASTN_REDO_CALLBACKS: BlastRedoAlignCallbacks = BlastRedoAlignCallbacks {
    calc_lambda: Some(redo_calc_lambda),
    get_range: translated_subject_get_range,
    redo_one_alignment: tblastn_redo_one_alignment_callback,
};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:795-800
/// The incoming alignments of a match (`s_ResultHspToDistinctAlign`): its preliminary HSPs in
/// order, each score times the redo's scale (`localScalingFactor`, `score_divisor`), each naming
/// its preliminary HSP. The X-drop redo (`redo_preliminary_match`) and the Smith-Waterman redo
/// (`match_redo::SmithWatermanRedo`) start from this list.
#[inline(always)]
pub(super) fn incoming_alignments(
    preliminary_hits: &[GappedHsp],
    query_context: i32,
    score_divisor: f64,
) -> Option<Box<BlastCompoAlignment>> {
    let mut incoming = None;
    let mut tail = &mut incoming;
    for (index, hit) in preliminary_hits.iter().enumerate() {
        let node = blast_compo_alignment_new(
            ((hit.score as f64) * score_divisor) as i32,
            EMatrixAdjustRule::DontAdjustMatrix,
            hit.q_start,
            hit.q_end,
            query_context,
            hit.s_start,
            hit.s_end,
            i32::from(hit.frame),
            Some(BlastCompoAlignmentContext::PreliminaryHspIndex(index)),
        );
        *tail = Some(node);
        tail = &mut tail
            .as_mut()
            .expect("TBLASTN incoming alignment inserted")
            .next;
    }
    incoming
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:770-803,3577-3639
#[allow(dead_code)] // The public TBLASTN path stays gated until D and E complete.
pub(super) fn redo_preliminary_match<'s>(
    preliminary_hits: &[GappedHsp],
    query_context: i32,
    query_infos: &[BlastCompoQueryInfo],
    // EXPERIMENT (LOSAT_X_TBNSTORE): the letters, or a packed subject.
    subject_nt: impl Into<crate::sequence::packed_subject::SubjectNt<'s>>,
    genetic_code_id: u8,
    params: &BlastRedoAlignParams,
    lambda: f64,
    matrix: ScoringMatrix,
    scratch: &mut GapAlignScratch,
    composition_workspace: &mut BlastCompositionWorkspace,
    // The matrix that the match's adjustments leave (NCBI's thread-local `sbp->matrix->data`,
    // which `Blast_AdjustScores` overwrites); the caller passes `None` (no adjustment yet) and
    // reads it after the redo (`postredo_engine_positives`).
    matrix_state: &mut Option<AdjustedProteinMatrix>,
) -> Result<BlastRedoOneMatchResult> {
    if params.compo_adjust_mode == BlastCompoAdjustMode::NoCompositionBasedStats {
        bail!("TBLASTN Kappa redo requires composition mode");
    }
    let incoming = incoming_alignments(preliminary_hits, query_context, params.score_divisor);
    let matching_seq =
        BlastCompoMatchingSequence::new_translated(0, subject_nt.into(), genetic_code_id)?;
    let context = TblastnRedoContext {
        preliminary_hits,
        matrix,
        scratch: NonNull::from(scratch),
    };
    let previous = params
        .gapping_params
        .context
        .replace(Some(NonNull::from(&context).cast::<c_void>()));
    let result = blast_redo_one_match_with_workspace_queries_and_matrix(
        &incoming,
        params,
        &matching_seq,
        query_infos,
        lambda,
        &TBLASTN_REDO_CALLBACKS,
        composition_workspace,
        matrix_state,
    );
    params.gapping_params.context.set(previous);
    result
}

#[cfg(test)]
#[path = "kappa_tests.rs"]
mod tests;
