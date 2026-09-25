//! TBLASTN translated-subject Kappa traceback boundary.
//! The public TBLASTN CLI remains gated until stages D and E agree.

use anyhow::{bail, Context, Result};
use std::ffi::c_void;
use std::ptr::NonNull;

use super::search_gapped::{GappedHsp, TargetTranslation};
use crate::algorithm::blastp::gapalign::{
    blast_gapped_alignment_with_traceback_with_scratch, protein_identities_from_edit_ops,
    GapAlignScratch,
};
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::{
    AdjustedProteinMatrix, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    blast_compo_alignment_new, blast_redo_one_match_with_workspace_queries, redo_calc_lambda,
    translated_subject_get_range, BlastCompoAdjustMode, BlastCompoAlignment,
    BlastCompoAlignmentContext, BlastCompoMatchingSequence, BlastCompoQueryInfo,
    BlastCompoSequenceData, BlastCompoSequenceRange, BlastRedoAlignCallbacks, BlastRedoAlignParams,
    BlastRedoOneMatchResult, EMatrixAdjustRule,
};
use crate::utils::genetic_code::GeneticCode;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1738-1773,1896-1957
// ```c
// q_start = hsp->query.gapped_start - query_range->begin;
// s_start = hsp->subject.gapped_start - subject_range->begin;
// gapAlign->gap_x_dropoff = gapping_params->x_dropoff;
// status = BLAST_GappedAlignmentWithTraceback(context->prog_number,
//     query_data->data, subject_data->data, gapAlign,
//     context->scoringParams, q_start, s_start,
//     query_data->length, subject_data->length, &fence_hit);
// if (status == 0) return s_NewAlignmentFromGapAlign(...);
// ```
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
    let query_start = i32::try_from(aligned.query_start)?
        .checked_add(query_range.begin)
        .context("TBLASTN redo query start overflow")?;
    let query_end = i32::try_from(aligned.query_stop)?
        .checked_add(query_range.begin)
        .context("TBLASTN redo query end overflow")?;
    let match_start = i32::try_from(aligned.subject_start)?
        .checked_add(subject_range.begin)
        .context("TBLASTN redo subject start overflow")?;
    let match_end = i32::try_from(aligned.subject_stop)?
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:224-273;
// c++/src/algo/blast/core/blast_hits_priv.h:65-70
// ```c
// if (hsp2->query.frame == hsp1->query.frame &&
//     hsp2->subject.frame == hsp1->subject.frame &&
//     CONTAINED_IN_HSP(q2start, q2end, q1start, s2start, s2end, s1start) &&
//     CONTAINED_IN_HSP(q2start, q2end, q1end, s2start, s2end, s1end) &&
//     hsp1->score <= hsp2->score)
//     hsp1 = hsp_array[iread] = Blast_HSPFree(hsp_array[iread]);
// #define CONTAINED_IN_HSP(a,b,c,d,e,f) \
//     (((a <= c && b >= c) && (d <= f && e >= f)) ? TRUE : FALSE)
// ```
#[allow(dead_code)] // Entered after TBLASTN Kappa alignment conversion is in the public D pipeline.
pub(super) fn reap_contained_postredo_hsps(hsps: &mut Vec<(usize, GappedHsp)>) {
    if hsps.len() < 2 {
        return;
    }
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
    let mut index = 0;
    hsps.retain(|_| {
        let selected = keep[index];
        index += 1;
        selected
    });
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:459-532;
// c++/src/algo/blast/core/blast_hits.c:1147-1228,767-811
// ```c
// query = query_blk->sequence + query_info->contexts[hsp->context].query_offset;
// const Uint1* target_sequence = Blast_HSPGetTargetTranslation(target_t, hsp, NULL);
// status = Blast_HSPGetNumIdentitiesAndPositives(query, target_sequence,
//                                                hsp, scoring_options, 0, sbp);
// if (*q == *s) { num_ident++; }
// ```
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
    let subject_start = usize::try_from(alignment.match_start)?
        .checked_sub(subject_base)
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:102-114;
// c++/src/algo/blast/core/ncbi_math.c:437-441
// ```c
// hsp->score = (Int4)BLAST_Nint(((double) hsp->score) / scoreDivisor);
// hsp->bit_score = (hsp->score*lambda*scoreDivisor - logK)/NCBIMATH_LN2;
// x += (x >= 0. ? 0.5 : -0.5); return (long)x;
// ```
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:770-803,1910-1948
// ```c
// new_align = BlastCompo_AlignmentNew((int)(hsp->score * localScalingFactor),
//     eDontAdjustMatrix, ..., hsp);
// BlastHSP *hsp = in_align->context;
// q_start = hsp->query.gapped_start - query_range->begin;
// s_start = hsp->subject.gapped_start - subject_range->begin;
// ```
struct TblastnRedoContext<'a> {
    preliminary_hits: &'a [GappedHsp],
    matrix: ScoringMatrix,
    scratch: NonNull<GapAlignScratch>,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1896-1957
// ```c
// BlastHSP * hsp = in_align->context;
// q_start = hsp->query.gapped_start - query_range->begin;
// s_start = hsp->subject.gapped_start - subject_range->begin;
// status = BLAST_GappedAlignmentWithTraceback(..., &fence_hit);
// ```
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2389-2394
// ```c
// static const Blast_RedoAlignCallbacks redo_align_callbacks = {
//     s_CalcLambda, s_SequenceGetRange, s_RedoOneAlignment,
//     s_NewAlignmentUsingXdrop, s_FreeEditScript
// };
// ```
const TBLASTN_REDO_CALLBACKS: BlastRedoAlignCallbacks = BlastRedoAlignCallbacks {
    calc_lambda: Some(redo_calc_lambda),
    get_range: translated_subject_get_range,
    redo_one_alignment: tblastn_redo_one_alignment_callback,
    new_xdrop_align: None,
    free_align_traceback: None,
};

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:770-803,3577-3639
// ```c
// s_ResultHspToDistinctAlign(incoming_align_set, numAligns,
//     localMatch->hsp_array, localMatch->hspcnt, context_index,
//     queryInfo, localScalingFactor);
// Blast_RedoOneMatch(alignments, redo_align_params, incoming_aligns,
//     numAligns[frame_index], kbp->Lambda, &matchingSeq, -1,
//     query_info, numContexts, matrix, BLASTAA_SIZE, NRrecord, ...);
// ```
#[allow(dead_code)] // The public TBLASTN path stays gated until D and E complete.
pub(super) fn redo_preliminary_match(
    preliminary_hits: &[GappedHsp],
    query_context: i32,
    query_infos: &[BlastCompoQueryInfo],
    subject_nt: &[u8],
    genetic_code_id: u8,
    params: &BlastRedoAlignParams,
    lambda: f64,
    matrix: ScoringMatrix,
    scratch: &mut GapAlignScratch,
    composition_workspace: &mut BlastCompositionWorkspace,
) -> Result<BlastRedoOneMatchResult> {
    if params.compo_adjust_mode == BlastCompoAdjustMode::NoCompositionBasedStats {
        bail!("TBLASTN Kappa redo requires composition mode");
    }
    let mut incoming = None;
    let mut tail = &mut incoming;
    for (index, hit) in preliminary_hits.iter().enumerate() {
        let node = blast_compo_alignment_new(
            ((hit.score as f64) * params.score_divisor) as i32,
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
    let matching_seq = BlastCompoMatchingSequence::new_translated(0, subject_nt, genetic_code_id)?;
    let context = TblastnRedoContext {
        preliminary_hits,
        matrix,
        scratch: NonNull::from(scratch),
    };
    let previous = params
        .gapping_params
        .context
        .replace(Some(NonNull::from(&context).cast::<c_void>()));
    let result = blast_redo_one_match_with_workspace_queries(
        &incoming,
        params,
        &matching_seq,
        query_infos,
        lambda,
        &TBLASTN_REDO_CALLBACKS,
        composition_workspace,
    );
    params.gapping_params.context.set(previous);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastn::stage_d_kappa_params::{
        local_extension_final_xdrop, local_kappa_redo_params,
    };
    use crate::common::GapEditOp;
    use crate::utils::matrix::BLASTAA_SIZE;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:770-803,1896-1957,3577-3689
    // ```c
    // s_ResultHspToDistinctAlign(..., localScalingFactor);
    // Blast_RedoOneMatch(..., incoming_aligns, ..., kbp->Lambda, ...);
    // if (alignments[context_index] != NULL)
    //     s_HSPListFromDistinctAlignments(...);
    // ```
    #[test]
    fn hard_seg_preliminary_hsp_redoes_to_ncbi_alignment() {
        use crate::core::composition_adjustment::adjust_scores::{
            build_matrix_info, read_aa_composition,
        };
        use crate::core::composition_adjustment::redo_alignment::{
            build_query_word_hashes, BlastCompoGappingParams,
        };
        use crate::utils::matrix::aa_char_to_ncbistdaa;
        use std::cell::Cell;
        let root = format!(
            "{}/../docs/evidence/tlosan_stage_d",
            env!("CARGO_MANIFEST_DIR")
        );
        let call_trace = std::fs::read_to_string(format!(
            "{root}/kappa_traceback_20260924/seg_hard_query_20260924_default.tsv"
        ))
        .unwrap();
        let call: Vec<_> = call_trace
            .lines()
            .find(|line| line.starts_with("K_TRACE_ENTER\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(call[3], "41");
        let mode_trace = std::fs::read_to_string(format!(
            "{root}/kappa_mode2_20260924/seg_hard_query_20260924_default.tsv"
        ))
        .unwrap();
        let incoming: Vec<_> = mode_trace
            .lines()
            .find(|line| line.starts_with("K_ALIGN\t0\tincoming\t"))
            .unwrap()
            .split('\t')
            .collect();
        let expected: Vec<_> = mode_trace
            .lines()
            .find(|line| line.starts_with("K_ALIGN\t0\tredone\t"))
            .unwrap()
            .split('\t')
            .collect();
        let redo: Vec<_> = mode_trace
            .lines()
            .find(|line| line.contains("\tredo_enter\t"))
            .unwrap()
            .split('\t')
            .collect();
        let query_fasta = std::fs::read_to_string(format!(
            "{}/../docs/evidence/tlosan_stage_c/seg_hard_query_20260924/query.faa",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let mut query: Vec<_> = query_fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.bytes().map(aa_char_to_ncbistdaa))
            .collect();
        query[..40].fill(21);
        let query_info = BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData::from_ncbistdaa(&query),
            composition: read_aa_composition(&query),
            eff_search_space: 15225.0,
            words: Some(build_query_word_hashes(&query)),
        };
        let subject_fasta = std::fs::read_to_string(format!(
            "{}/../docs/evidence/tlosan_stage_c/seg_hard_query_20260924/subjects.fna",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let subject_nt: Vec<_> = subject_fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.bytes())
            .collect();
        let hsp = GappedHsp {
            frame: incoming[12].parse().unwrap(),
            score: 656,
            q_start: incoming[8].parse().unwrap(),
            q_end: incoming[9].parse().unwrap(),
            q_gapped_start: call[4].parse().unwrap(),
            s_start: incoming[10].parse().unwrap(),
            s_end: incoming[11].parse().unwrap(),
            s_gapped_start: call[5].parse().unwrap(),
        };
        assert_eq!(
            (hsp.score as f64 * 32.0) as i32,
            incoming[5].parse().unwrap()
        );
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:870-899;
        // ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3577-3595
        // ```c
        // status = BLAST_LinkHsps(program_number, hsp_list_out, ...);
        // status = s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
        // s_ResultHspToDistinctAlign(incoming_align_set, numAligns,
        //     localMatch->hsp_array, localMatch->hspcnt, ...);
        // ```
        use crate::algorithm::tblastn::search_gapped::{
            preliminary_protein_hsps_in_ncbi_order, PreliminaryProfile,
        };
        use crate::algorithm::tblastn::stage_d_linking::{link_preliminary_hsps, reap_by_evalue};
        use crate::algorithm::tblastn::stage_d_stats::{
            local_parameters_for_call, LocalParameterCall, LocalParameterOptions,
        };
        use crate::config::ProteinScoringSpec;
        use crate::stats::spouge::lookup_protein_gumbel_params;
        use crate::stats::tables::{lookup_protein_params_gapped, lookup_protein_params_ungapped};
        use crate::utils::seg::SegParams;
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            (subject_nt.len() / 3) as i64,
        )
        .unwrap();
        let parameters = local_parameters_for_call(
            &[(160, true)],
            subject_nt.len(),
            &[gapped],
            &[ungapped],
            LocalParameterOptions {
                expect_value: 10.0,
                do_sum_stats: true,
                max_intron_length: 0,
                gap_trigger_bits: 22.0,
                word_xdrop_bits: 7.0,
                scale_factor: 1.0,
                gumbel: Some(&gumbel),
            },
            LocalParameterCall::Initial {
                min_subject_length: (subject_nt.len() / 3) as i32,
                composition_based_stats: 2,
            },
        );
        assert_eq!(parameters.lengths[0].eff_searchsp, 15225);
        let seg = SegParams::default();
        let search_query: Vec<_> = query_fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.bytes())
            .collect();
        let word_xdrop = [parameters.cutoffs[0].word_xdrop];
        let word_cutoff = [parameters.cutoffs[0].word_cutoff];
        let hit_cutoff = [parameters.cutoffs[0].hit_cutoff];
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &word_xdrop,
            word_cutoff: &word_cutoff,
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &hit_cutoff,
            hsp_num_max: i32::MAX as usize,
        };
        let (stage_c_hsps, _) =
            preliminary_protein_hsps_in_ncbi_order(&[&search_query], &subject_nt, 1, profile)
                .unwrap();
        assert_eq!(stage_c_hsps, vec![(0, hsp)]);
        let mut preliminary_linked = link_preliminary_hsps(
            &stage_c_hsps,
            &[160],
            &parameters.lengths,
            subject_nt.len() as i32,
            &[gapped],
            &gumbel,
            parameters.link.as_ref().unwrap(),
        )
        .unwrap();
        reap_by_evalue(&mut preliminary_linked, parameters.prelim_evalue);
        assert_eq!(preliminary_linked.hsps.len(), 1);
        let d_trace = std::fs::read_to_string(format!(
            "{root}/run_20260924/seg_hard_query_20260924_default.trace"
        ))
        .unwrap();
        let linked_row: Vec<_> = d_trace
            .lines()
            .find(|line| line.starts_with("D_HSP\t0\tlink_after\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(preliminary_linked.hsps[0].hsp, hsp);
        assert_eq!(
            preliminary_linked.hsps[0].num,
            linked_row[12].parse().unwrap()
        );
        assert_eq!(
            preliminary_linked.hsps[0].evalue.to_bits(),
            linked_row[13].parse::<f64>().unwrap().to_bits()
        );
        let hsp = preliminary_linked.hsps[0].hsp;
        // NCBI reference: core/blast_kappa.c:2352-2390,2418-2479:
        // s_GetAlignParams reads the initial hit cutoffs and scaled score block.
        let params = local_kappa_redo_params(
            ScoringMatrix::Blosum62,
            11,
            1,
            &[gapped],
            &[true],
            &parameters,
            160,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
            false,
            10.0,
            true,
            25.0,
            local_extension_final_xdrop(15.0, 25.0, gapped.lambda).unwrap(),
        )
        .unwrap();
        assert_eq!(params.gapping_params.x_dropoff, call[8].parse().unwrap());
        assert_eq!(params.gapping_params.gap_open, call[9].parse().unwrap());
        assert_eq!(params.gapping_params.gap_extend, call[10].parse().unwrap());
        assert_eq!(params.cutoff_score, redo[11].parse().unwrap());
        assert_eq!(
            (gapped.lambda / 32.0).to_bits(),
            redo[4].parse::<f64>().unwrap().to_bits()
        );
        let mut scratch = GapAlignScratch::new();
        let mut workspace = BlastCompositionWorkspace::new_blosum62();
        let redone = redo_preliminary_match(
            &[hsp],
            0,
            &[query_info],
            &subject_nt,
            1,
            &params,
            gapped.lambda / 32.0,
            ScoringMatrix::Blosum62,
            &mut scratch,
            &mut workspace,
        )
        .unwrap();
        let align = redone.alignments_by_query[0].as_ref().unwrap();
        assert_eq!(align.score, expected[5].parse().unwrap());
        assert_eq!(
            align.matrix_adjust_rule as i32,
            expected[6].parse().unwrap()
        );
        assert_eq!(align.query_start, expected[8].parse().unwrap());
        assert_eq!(align.query_end, expected[9].parse().unwrap());
        assert_eq!(align.match_start, expected[10].parse().unwrap());
        assert_eq!(align.match_end, expected[11].parse().unwrap());
        assert_eq!(align.frame, expected[12].parse().unwrap());
        assert_eq!(redone.lambda_ratio.unwrap().to_bits(), 1.0f64.to_bits());
        assert!(params.gapping_params.context.get().is_none());
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1460-1466;
    // core/blast_kappa.c:1475-1552,3577-3736;
    // composition_adjustment/compo_heap.c:252-275
    // ```c
    // if (seq_arg.seq->gen_code_string == NULL)
    //     seq_arg.seq->gen_code_string = GenCodeSingletonFind(db_options->genetic_code);
    // Blast_RedoOneMatch(..., incoming_aligns, ..., kbp->Lambda, ...);
    // s_HitlistEvaluateAndPurge(...);
    // s_HSPListNormalizeScores(hsp_list, kbp->Lambda, kbp->logK,
    //                          localScalingFactor);
    // BlastCompo_HeapInsert(..., best_evalue, best_score, ...);
    // ```
    #[test]
    fn code32_local_api_subject_translation_redo_and_scores_match_ncbi() {
        use crate::algorithm::tblastn::stage_d_linking::{link_preliminary_hsps, reap_by_evalue};
        use crate::algorithm::tblastn::stage_d_stats::{
            local_parameters_for_call, LocalParameterCall, LocalParameterOptions,
        };
        use crate::config::ProteinScoringSpec;
        use crate::core::composition_adjustment::adjust_scores::{
            blast_adjust_scores, build_matrix_info, read_aa_composition,
        };
        use crate::core::composition_adjustment::redo_alignment::{
            build_query_word_hashes, BlastCompoGappingParams,
        };
        use crate::stats::spouge::lookup_protein_gumbel_params;
        use crate::stats::tables::{
            lookup_protein_params_gapped, lookup_protein_params_ungapped, KarlinParams,
        };
        use crate::utils::matrix::aa_char_to_ncbistdaa;
        use std::cell::Cell;

        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/code32_local_api_20260925_cli_calibrated"
        );
        let stage_a = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_a/fixtures"
        );
        let query_fasta = std::fs::read_to_string(format!("{stage_a}/query.faa")).unwrap();
        let query: Vec<_> = query_fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.bytes().map(aa_char_to_ncbistdaa))
            .collect();
        let subject_fasta =
            std::fs::read_to_string(format!("{stage_a}/subject_code32.fna")).unwrap();
        let subject_nt: Vec<_> = subject_fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.bytes())
            .collect();
        assert_eq!((query.len(), subject_nt.len()), (120, 360));
        let trace = std::fs::read_to_string(format!("{root}/code32.trace")).unwrap();
        let mode = std::fs::read_to_string(format!("{root}/code32_mode2.trace")).unwrap();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3577-3595
        // ```c
        // s_ResultHspToDistinctAlign(incoming_align_set, numAligns,
        //     localMatch->hsp_array, localMatch->hspcnt, ...);
        // incoming_aligns = incoming_align_set[frame_index];
        // ```
        let prelim_rows: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_PRELIM\t"))
            .collect();
        let ncbi_preliminary: Vec<GappedHsp> = prelim_rows
            .iter()
            .enumerate()
            .map(|(index, row)| {
                let f: Vec<_> = row.split('\t').collect();
                assert_eq!(f[2].parse::<usize>().unwrap(), index);
                assert_eq!(f[4], "0");
                GappedHsp {
                    frame: f[5].parse().unwrap(),
                    score: f[3].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[9].parse().unwrap(),
                    s_end: f[10].parse().unwrap(),
                    s_gapped_start: f[11].parse().unwrap(),
                }
            })
            .collect();
        assert_eq!(ncbi_preliminary.len(), 2);
        assert_eq!(
            ncbi_preliminary.iter().map(|h| h.score).collect::<Vec<_>>(),
            [656, 16]
        );
        let call: Vec<_> = trace
            .lines()
            .find(|line| line.starts_with("K_TRACE_ENTER\t"))
            .unwrap()
            .split('\t')
            .collect();
        let redo: Vec<_> = mode
            .lines()
            .find(|line| line.contains("\tredo_enter\t"))
            .unwrap()
            .split('\t')
            .collect();
        let expected: Vec<_> = mode
            .lines()
            .find(|line| line.contains("\tredone\t"))
            .unwrap()
            .split('\t')
            .collect();
        let heap: Vec<_> = trace
            .lines()
            .find(|line| line.starts_with("K_TRACE_HEAP_HSP\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(redo[3], "2");
        assert_eq!(redo[6], "1");
        assert_eq!(redo[7], "360");
        assert_eq!(redo[8], "2");
        assert_eq!(redo[9], "1");
        assert_eq!(redo[11], "288");
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            120,
        )
        .unwrap();
        let parameters = local_parameters_for_call(
            &[(query.len(), true)],
            subject_nt.len(),
            &[gapped],
            &[ungapped],
            LocalParameterOptions {
                expect_value: 10.0,
                do_sum_stats: true,
                max_intron_length: 0,
                gap_trigger_bits: 22.0,
                word_xdrop_bits: 7.0,
                scale_factor: 1.0,
                gumbel: Some(&gumbel),
            },
            LocalParameterCall::Initial {
                min_subject_length: 120,
                composition_based_stats: 2,
            },
        );
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:870-899;
        // ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3577-3595
        // ```c
        // status = BLAST_LinkHsps(program_number, hsp_list_out, ...);
        // status = s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
        // s_ResultHspToDistinctAlign(incoming_align_set, numAligns,
        //     localMatch->hsp_array, localMatch->hspcnt, ...);
        // ```
        use crate::algorithm::tblastn::search_gapped::{
            preliminary_protein_hsps_in_ncbi_order, PreliminaryProfile,
        };
        use crate::utils::seg::SegParams;
        let search_query: Vec<_> = query_fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.bytes())
            .collect();
        let word_xdrop = [parameters.cutoffs[0].word_xdrop];
        let word_cutoff = [parameters.cutoffs[0].word_cutoff];
        let hit_cutoff = [parameters.cutoffs[0].hit_cutoff];
        let seg = SegParams::default();
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &word_xdrop,
            word_cutoff: &word_cutoff,
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &hit_cutoff,
            hsp_num_max: i32::MAX as usize,
        };
        let (stage_c_hsps, _) =
            preliminary_protein_hsps_in_ncbi_order(&[&search_query], &subject_nt, 32, profile)
                .unwrap();
        assert_eq!(
            stage_c_hsps,
            ncbi_preliminary
                .iter()
                .copied()
                .map(|hsp| (0, hsp))
                .collect::<Vec<_>>()
        );
        let mut preliminary_linked = link_preliminary_hsps(
            &stage_c_hsps,
            &[query.len() as i32],
            &parameters.lengths,
            subject_nt.len() as i32,
            &[gapped],
            &gumbel,
            parameters.link.as_ref().unwrap(),
        )
        .unwrap();
        reap_by_evalue(&mut preliminary_linked, parameters.prelim_evalue);
        assert_eq!(preliminary_linked.hsps.len(), 2);
        let preliminary_hsps: Vec<_> = preliminary_linked
            .hsps
            .iter()
            .map(|linked| linked.hsp)
            .collect();
        assert_eq!(preliminary_hsps, ncbi_preliminary);
        let d_trace = std::fs::read_to_string(format!("{root}/code32_d.trace")).unwrap();
        let linked_rows: Vec<_> = d_trace
            .lines()
            .filter(|line| line.starts_with("D_HSP\t0\tlink_after\t"))
            .collect();
        assert_eq!(linked_rows.len(), preliminary_linked.hsps.len());
        for (linked, row) in preliminary_linked.hsps.iter().zip(linked_rows) {
            let f: Vec<_> = row.split('\t').collect();
            assert_eq!(linked.hsp.score, f[10].parse().unwrap());
            assert_eq!(linked.num, f[12].parse().unwrap());
            assert_eq!(
                linked.evalue.to_bits(),
                f[13].parse::<f64>().unwrap().to_bits()
            );
        }
        let query_info = BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData::from_ncbistdaa(&query),
            composition: read_aa_composition(&query),
            eff_search_space: parameters.lengths[0].eff_searchsp as f64,
            words: Some(build_query_word_hashes(&query)),
        };
        // NCBI reference: core/blast_kappa.c:2352-2390,2418-2479:
        // s_GetAlignParams reads the initial hit cutoffs and scaled score block.
        let params = local_kappa_redo_params(
            ScoringMatrix::Blosum62,
            11,
            1,
            &[gapped],
            &[true],
            &parameters,
            query.len() as i32,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
            false,
            10.0,
            true,
            25.0,
            local_extension_final_xdrop(15.0, 25.0, gapped.lambda).unwrap(),
        )
        .unwrap();
        assert_eq!(params.gapping_params.x_dropoff, call[8].parse().unwrap());
        assert_eq!(params.gapping_params.gap_open, call[9].parse().unwrap());
        assert_eq!(params.gapping_params.gap_extend, call[10].parse().unwrap());
        assert_eq!(params.cutoff_score, redo[11].parse().unwrap());
        assert_eq!(
            (gapped.lambda / 32.0).to_bits(),
            redo[4].parse::<f64>().unwrap().to_bits()
        );
        // NCBI reference: composition_adjustment/redo_alignment.c:1219-1254;
        // composition_adjustment/composition_adjustment.c:1414-1530
        // ```c
        // s_GetComposition(&subject_composition, alphsize, &subject,
        //     &window->subject_range, in_align, FALSE, subject_is_translated);
        // Blast_AdjustScores(matrix, query_composition, query.length,
        //     &subject_composition, subject.length, scaledMatrixInfo,
        //     compo_adjust_mode, RE_pseudocounts, NRrecord, ...);
        // ```
        let comp_trace =
            std::fs::read_to_string(format!("{root}/code32_composition.trace")).unwrap();
        let query_composition = read_aa_composition(&query);
        for which in ["query", "subject"] {
            let prefix = format!("K_COMP\t1\t{which}\t");
            let row: Vec<_> = comp_trace
                .lines()
                .find(|line| line.starts_with(&prefix))
                .unwrap()
                .split('\t')
                .collect();
            assert_eq!(
                row[3].parse::<i32>().unwrap(),
                query_composition.num_true_amino_acids
            );
            for (actual, bits) in query_composition.prob.iter().zip(&row[4..]) {
                assert_eq!(actual.to_bits(), u64::from_str_radix(bits, 16).unwrap());
            }
        }
        let mut matrix_workspace = BlastCompositionWorkspace::new_blosum62();
        let adjusted = blast_adjust_scores(
            &params.matrix_info,
            &query_composition,
            query.len() as i32,
            &query_composition,
            query.len() as i32,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
            0,
            &mut matrix_workspace,
            redo_calc_lambda,
        )
        .unwrap()
        .unwrap();
        let result: Vec<_> = comp_trace
            .lines()
            .find(|line| line.starts_with("K_COMP_RESULT\t1\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(
            adjusted.matrix_adjust_rule as i32,
            result[3].parse().unwrap()
        );
        assert_eq!(
            adjusted.lambda_ratio.to_bits(),
            result[4].parse::<f64>().unwrap().to_bits()
        );
        let matrix: Vec<_> = comp_trace
            .lines()
            .find(|line| line.starts_with("K_ADJUSTED\t1\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(matrix.len(), 2 + BLASTAA_SIZE * BLASTAA_SIZE);
        for (row, cells) in matrix[2..].chunks_exact(BLASTAA_SIZE).enumerate() {
            for (col, value) in cells.iter().enumerate() {
                assert_eq!(
                    adjusted.adjusted_matrix.scores[row][col],
                    value.parse::<i32>().unwrap(),
                    "code32 matrix {row},{col}"
                );
            }
        }
        let mut scratch = GapAlignScratch::new();
        let mut workspace = BlastCompositionWorkspace::new_blosum62();
        let redone = redo_preliminary_match(
            &preliminary_hsps,
            0,
            &[query_info],
            &subject_nt,
            32,
            &params,
            gapped.lambda / 32.0,
            ScoringMatrix::Blosum62,
            &mut scratch,
            &mut workspace,
        )
        .unwrap();
        let alignment = redone.alignments_by_query[0].as_deref().unwrap();
        assert_eq!(alignment.score, expected[5].parse().unwrap());
        assert_eq!(
            alignment.matrix_adjust_rule as i32,
            expected[6].parse().unwrap()
        );
        assert_eq!(alignment.query_start, expected[8].parse().unwrap());
        assert_eq!(alignment.query_end, expected[9].parse().unwrap());
        assert_eq!(alignment.match_start, expected[10].parse().unwrap());
        assert_eq!(alignment.match_end, expected[11].parse().unwrap());
        assert_eq!(alignment.frame, expected[12].parse().unwrap());
        assert_eq!(redone.lambda_ratio.unwrap().to_bits(), 1.0f64.to_bits());
        let code = GeneticCode::try_from_id(32).unwrap();
        let mut target = TargetTranslation::new(&subject_nt, &code);
        assert_eq!(
            postredo_num_ident(alignment, &query, &mut target).unwrap(),
            heap[6].parse().unwrap()
        );
        let scaled = KarlinParams {
            lambda: gapped.lambda / 32.0,
            ..gapped
        };
        let mut link_params = parameters.link.unwrap();
        link_params.cutoff_small_gap = 0;
        let mut linked = link_preliminary_hsps(
            &[(
                0,
                GappedHsp {
                    frame: alignment.frame as i8,
                    score: alignment.score,
                    q_start: alignment.query_start,
                    q_end: alignment.query_end,
                    q_gapped_start: 0,
                    s_start: alignment.match_start,
                    s_end: alignment.match_end,
                    s_gapped_start: 0,
                },
            )],
            &[query.len() as i32],
            &parameters.lengths,
            subject_nt.len() as i32,
            &[scaled],
            &gumbel,
            &link_params,
        )
        .unwrap();
        reap_by_evalue(&mut linked, 10.0);
        assert_eq!(linked.hsps.len(), 1);
        let would: Vec<_> = trace
            .lines()
            .find(|line| line.starts_with("K_TRACE_HEAP_WOULD\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(linked.hsps[0].hsp.score, would[3].parse().unwrap());
        assert_eq!(
            linked.best_evalue.to_bits(),
            would[2].parse::<f64>().unwrap().to_bits()
        );
        let bit =
            normalize_postredo_scores(&mut linked, gapped.lambda / 32.0, gapped.k.ln(), 32.0)[0];
        assert_eq!(linked.hsps[0].hsp.score, heap[3].parse().unwrap());
        assert_eq!(bit.to_bits(), heap[4].parse::<f64>().unwrap().to_bits());
        assert_eq!(
            linked.hsps[0].evalue.to_bits(),
            heap[5].parse::<f64>().unwrap().to_bits()
        );
        assert!(params.gapping_params.context.get().is_none());
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:770-803,3577-3689;
    // composition_adjustment/redo_alignment.c:1101-1292
    // ```c
    // s_ResultHspToDistinctAlign(incoming_align_set, numAligns,
    //     localMatch->hsp_array, localMatch->hspcnt, context_index,
    //     queryInfo, localScalingFactor);
    // for (frame_index = 0; frame_index < numFrames; frame_index++)
    //     Blast_RedoOneMatch(alignments, ..., incoming_aligns,
    //         numAligns[frame_index], kbp->Lambda, ..., query_info,
    //         numContexts, ...);
    // ```
    #[test]
    fn multi_query_preliminary_hsps_redo_to_ncbi_alignment_order() {
        use crate::core::composition_adjustment::adjust_scores::{
            build_matrix_info, read_aa_composition,
        };
        use crate::core::composition_adjustment::redo_alignment::{
            build_query_word_hashes, BlastCompoGappingParams,
        };
        use crate::utils::matrix::aa_char_to_ncbistdaa;
        use std::cell::Cell;
        let root = format!(
            "{}/../docs/evidence/tlosan_stage_d",
            env!("CARGO_MANIFEST_DIR")
        );
        let mode = std::fs::read_to_string(format!(
            "{root}/kappa_mode2_20260924/multi_query_20260924_default.tsv"
        ))
        .unwrap();
        let call_trace = std::fs::read_to_string(format!(
            "{root}/kappa_traceback_20260924/multi_query_20260924_default.tsv"
        ))
        .unwrap();
        let query_fasta = std::fs::read_to_string(format!(
            "{}/../docs/evidence/tlosan_stage_c/multi_query_20260924/query.faa",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let mut queries: Vec<Vec<u8>> = Vec::new();
        for line in query_fasta.lines() {
            if line.starts_with('>') {
                queries.push(Vec::new());
            } else {
                queries
                    .last_mut()
                    .unwrap()
                    .extend(line.bytes().map(aa_char_to_ncbistdaa));
            }
        }
        assert_eq!(
            queries.iter().map(Vec::len).collect::<Vec<_>>(),
            vec![120, 70, 120]
        );
        let query_infos: Vec<_> = queries
            .iter()
            .enumerate()
            .map(|(index, query)| BlastCompoQueryInfo {
                origin: 0,
                seq: BlastCompoSequenceData::from_ncbistdaa(query),
                composition: read_aa_composition(query),
                eff_search_space: [182004.0, 88074.0, 182004.0][index],
                words: Some(build_query_word_hashes(query)),
            })
            .collect();
        let subject_fasta = std::fs::read_to_string(format!(
            "{}/../docs/evidence/tlosan_stage_c/multi_query_20260924/subjects.fna",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let subject_nt: Vec<_> = subject_fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.bytes())
            .collect();
        // NCBI core/blast_kappa.c:3658-3685 converts, score-sorts,
        // relinks, and reaps the generated Kappa alignments in this order.
        use crate::algorithm::tblastn::stage_d_linking::{
            link_preliminary_hsps, reap_by_evalue, score_compare,
        };
        use crate::algorithm::tblastn::stage_d_stats::{
            local_parameters_for_call, LocalParameterCall, LocalParameterOptions,
        };
        use crate::config::ProteinScoringSpec;
        use crate::stats::spouge::lookup_protein_gumbel_params;
        use crate::stats::tables::{
            lookup_protein_params_gapped, lookup_protein_params_ungapped, KarlinParams,
        };
        let d_trace = std::fs::read_to_string(format!(
            "{root}/run_20260924/multi_query_20260924_default.trace"
        ))
        .unwrap();
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let scaled = KarlinParams {
            lambda: gapped.lambda / 32.0,
            ..gapped
        };
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            2_125,
        )
        .unwrap();
        let query_lengths = [120, 70, 120];
        let parameters = local_parameters_for_call(
            &[(120, true), (70, true), (120, false)],
            6_377,
            &[gapped; 3],
            &[ungapped; 3],
            LocalParameterOptions {
                expect_value: 10.0,
                do_sum_stats: true,
                max_intron_length: 0,
                gap_trigger_bits: 22.0,
                word_xdrop_bits: 7.0,
                scale_factor: 1.0,
                gumbel: Some(&gumbel),
            },
            LocalParameterCall::Initial {
                min_subject_length: 2_125,
                composition_based_stats: 2,
            },
        );
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:870-899;
        // ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3577-3595
        // ```c
        // status = BLAST_LinkHsps(program_number, hsp_list_out, ...);
        // status = s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
        // s_ResultHspToDistinctAlign(incoming_align_set, numAligns,
        //     localMatch->hsp_array, localMatch->hspcnt, ...);
        // incoming_aligns = incoming_align_set[frame_index];
        // ```
        // Exercise the natural Stage C output through the same preliminary
        // link/reap before passing its retained HSPs to Kappa.
        use crate::algorithm::tblastn::search_gapped::{
            preliminary_protein_hsps_in_ncbi_order, PreliminaryProfile,
        };
        use crate::utils::seg::SegParams;
        let word_xdrop: Vec<_> = parameters.cutoffs.iter().map(|p| p.word_xdrop).collect();
        let word_cutoff: Vec<_> = parameters.cutoffs.iter().map(|p| p.word_cutoff).collect();
        let hit_cutoff: Vec<_> = parameters.cutoffs.iter().map(|p| p.hit_cutoff).collect();
        let seg = SegParams::default();
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &word_xdrop,
            word_cutoff: &word_cutoff,
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &hit_cutoff,
            hsp_num_max: i32::MAX as usize,
        };
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:614-625;
        // ncbi-blast/c++/src/algo/blast/core/blast_engine.c:484-525
        // ```c
        // BlastSetUp_MaskQuery(query_blk, ...);
        // WordFinder(..., query, ...);
        // ```
        // Stage C encodes FASTA protein bytes for lookup; query_infos above
        // hold the NCBISTDAA bytes consumed by Kappa.
        let mut search_queries: Vec<Vec<u8>> = Vec::new();
        for line in query_fasta.lines() {
            if line.starts_with('>') {
                search_queries.push(Vec::new());
            } else {
                search_queries.last_mut().unwrap().extend(line.bytes());
            }
        }
        let query_refs: Vec<_> = search_queries.iter().map(Vec::as_slice).collect();
        let (stage_c_hsps, _) =
            preliminary_protein_hsps_in_ncbi_order(&query_refs, &subject_nt, 1, profile).unwrap();
        assert_eq!(stage_c_hsps.len(), 21);
        let mut preliminary_linked = link_preliminary_hsps(
            &stage_c_hsps,
            &query_lengths,
            &parameters.lengths,
            subject_nt.len() as i32,
            &[gapped; 3],
            &gumbel,
            parameters.link.as_ref().unwrap(),
        )
        .unwrap();
        reap_by_evalue(&mut preliminary_linked, parameters.prelim_evalue);
        assert_eq!(preliminary_linked.hsps.len(), 20);
        let stage_c_by_query: Vec<Vec<GappedHsp>> = (0..queries.len())
            .map(|context| {
                preliminary_linked
                    .hsps
                    .iter()
                    .filter(|linked| linked.context == context)
                    .map(|linked| linked.hsp)
                    .collect()
            })
            .collect();
        let mut post_link = parameters.link.unwrap();
        post_link.cutoff_small_gap = 0;
        let redo_rows: Vec<_> = mode
            .lines()
            .filter(|line| line.contains("\tredo_enter\t"))
            .collect();
        assert_eq!(redo_rows.len(), 2);
        for (redo_index, redo_row) in redo_rows.iter().enumerate() {
            let redo: Vec<_> = redo_row.split('\t').collect();
            let prelim_rows: Vec<_> = call_trace
                .lines()
                .filter(|line| line.starts_with(&format!("K_TRACE_PRELIM\t{redo_index}\t")))
                .collect();
            assert_eq!(prelim_rows.len(), redo[3].parse::<usize>().unwrap());
            let mut preliminary = Vec::new();
            let mut query_context = None;
            for (index, row) in prelim_rows.iter().enumerate() {
                let f: Vec<_> = row.split('\t').collect();
                assert_eq!(f[2].parse::<usize>().unwrap(), index);
                let context = f[4].parse::<i32>().unwrap();
                assert!(query_context.is_none_or(|current| current == context));
                query_context = Some(context);
                preliminary.push(GappedHsp {
                    frame: f[5].parse().unwrap(),
                    score: f[3].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[9].parse().unwrap(),
                    s_end: f[10].parse().unwrap(),
                    s_gapped_start: f[11].parse().unwrap(),
                });
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3577-3595
            // ```c
            // s_ResultHspToDistinctAlign(incoming_align_set, numAligns,
            //     localMatch->hsp_array, localMatch->hspcnt, ...);
            // incoming_aligns = incoming_align_set[frame_index];
            // ```
            // Compare all saved incoming fields, then redo the HSPs computed
            // by Rust Stage C rather than substituting the trace as input.
            assert_eq!(stage_c_by_query[redo_index], preliminary);
            let preliminary = stage_c_by_query[redo_index].clone();
            let context = query_context.unwrap();
            assert_eq!(context, redo_index as i32);
            let first_call = call_trace
                .lines()
                .find(|line| line.starts_with(&format!("K_TRACE_ENTER\t{redo_index}\t")))
                .unwrap();
            let call: Vec<_> = first_call.split('\t').collect();
            // NCBI reference: core/blast_kappa.c:2352-2390,2418-2479:
            // s_GetAlignParams reads the initial hit cutoffs and scaled score block.
            let params = local_kappa_redo_params(
                ScoringMatrix::Blosum62,
                11,
                1,
                &[gapped; 3],
                &[true, true, false],
                &parameters,
                120,
                BlastCompoAdjustMode::CompositionMatrixAdjust,
                false,
                10.0,
                true,
                25.0,
                local_extension_final_xdrop(15.0, 25.0, gapped.lambda).unwrap(),
            )
            .unwrap();
            assert_eq!(params.gapping_params.x_dropoff, call[8].parse().unwrap());
            assert_eq!(params.gapping_params.gap_open, call[9].parse().unwrap());
            assert_eq!(params.gapping_params.gap_extend, call[10].parse().unwrap());
            assert_eq!(params.cutoff_score, redo[11].parse().unwrap());
            assert_eq!(
                (gapped.lambda / 32.0).to_bits(),
                redo[4].parse::<f64>().unwrap().to_bits()
            );
            let mut scratch = GapAlignScratch::new();
            let mut workspace = BlastCompositionWorkspace::new_blosum62();
            let redone = redo_preliminary_match(
                &preliminary,
                context,
                &query_infos,
                &subject_nt,
                1,
                &params,
                gapped.lambda / 32.0,
                ScoringMatrix::Blosum62,
                &mut scratch,
                &mut workspace,
            )
            .unwrap();
            let expected: Vec<_> = mode
                .lines()
                .filter(|line| line.starts_with(&format!("K_ALIGN\t{}\tredone\t", redo[1])))
                .collect();
            let mut observed = redone.alignments_by_query[context as usize].as_deref();
            for (index, row) in expected.iter().enumerate() {
                let align = observed
                    .unwrap_or_else(|| panic!("redo {redo_index} missing alignment {index}"));
                let f: Vec<_> = row.split('\t').collect();
                assert_eq!(
                    align.score,
                    f[5].parse().unwrap(),
                    "redo {redo_index}, HSP {index} score"
                );
                assert_eq!(
                    align.matrix_adjust_rule as i32,
                    f[6].parse().unwrap(),
                    "redo {redo_index}, HSP {index} rule"
                );
                assert_eq!(align.query_index, f[7].parse().unwrap());
                assert_eq!(
                    align.query_start,
                    f[8].parse().unwrap(),
                    "redo {redo_index}, HSP {index} qstart"
                );
                assert_eq!(
                    align.query_end,
                    f[9].parse().unwrap(),
                    "redo {redo_index}, HSP {index} qend"
                );
                assert_eq!(
                    align.match_start,
                    f[10].parse().unwrap(),
                    "redo {redo_index}, HSP {index} sstart"
                );
                assert_eq!(
                    align.match_end,
                    f[11].parse().unwrap(),
                    "redo {redo_index}, HSP {index} send"
                );
                assert_eq!(
                    align.frame,
                    f[12].parse().unwrap(),
                    "redo {redo_index}, HSP {index} frame"
                );
                observed = align.next.as_deref();
            }
            assert!(observed.is_none(), "redo {redo_index} has extra alignments");
            // NCBI reference: core/blast_kappa.c:305-357,3658-3685:
            // ```c
            // s_HSPListFromDistinctAlignments(...);
            // Blast_HSPListSortByScore(hsp_list);
            // if (hsp_list->hspcnt > 1) s_HitlistReapContained(...);
            // s_HitlistEvaluateAndPurge(...);
            // ```
            let mut linked_input = Vec::new();
            let mut align = redone.alignments_by_query[context as usize].as_deref();
            while let Some(value) = align {
                linked_input.push((
                    context as usize,
                    GappedHsp {
                        frame: value.frame as i8,
                        score: value.score,
                        q_start: value.query_start,
                        q_end: value.query_end,
                        q_gapped_start: 0,
                        s_start: value.match_start,
                        s_end: value.match_end,
                        s_gapped_start: 0,
                    },
                ));
                align = value.next.as_deref();
            }
            linked_input.sort_by(|a, b| score_compare(&a.1, &b.1));
            reap_contained_postredo_hsps(&mut linked_input);
            let (link_event, reap_event) = if redo_index == 0 { (3, 5) } else { (6, 8) };
            let before_prefix = format!("D_HSP\t{link_event}\tlink_before\t");
            let before: Vec<_> = d_trace
                .lines()
                .filter(|line| line.starts_with(&before_prefix))
                .collect();
            assert_eq!(linked_input.len(), before.len());
            for (index, ((ctx, hsp), row)) in linked_input.iter().zip(before).enumerate() {
                let f: Vec<_> = row.split('\t').collect();
                assert_eq!(
                    *ctx,
                    f[4].parse().unwrap(),
                    "redo {redo_index} prelink {index} context"
                );
                assert_eq!(
                    hsp.frame,
                    f[5].parse().unwrap(),
                    "redo {redo_index} prelink {index} frame"
                );
                assert_eq!(
                    hsp.q_start,
                    f[6].parse().unwrap(),
                    "redo {redo_index} prelink {index} qstart"
                );
                assert_eq!(
                    hsp.q_end,
                    f[7].parse().unwrap(),
                    "redo {redo_index} prelink {index} qend"
                );
                assert_eq!(
                    hsp.s_start,
                    f[8].parse().unwrap(),
                    "redo {redo_index} prelink {index} sstart"
                );
                assert_eq!(
                    hsp.s_end,
                    f[9].parse().unwrap(),
                    "redo {redo_index} prelink {index} send"
                );
                assert_eq!(
                    hsp.score,
                    f[10].parse().unwrap(),
                    "redo {redo_index} prelink {index} score"
                );
            }
            let mut linked = link_preliminary_hsps(
                &linked_input,
                &query_lengths,
                &parameters.lengths,
                6_377,
                &[scaled; 3],
                &gumbel,
                &post_link,
            )
            .unwrap();
            let after_prefix = format!("D_HSP\t{link_event}\tlink_after\t");
            let after: Vec<_> = d_trace
                .lines()
                .filter(|line| line.starts_with(&after_prefix))
                .collect();
            assert_eq!(linked.hsps.len(), after.len());
            for (index, (hsp, row)) in linked.hsps.iter().zip(after).enumerate() {
                let f: Vec<_> = row.split('\t').collect();
                assert_eq!(
                    hsp.hsp.score,
                    f[10].parse().unwrap(),
                    "redo {redo_index} linked {index} score"
                );
                assert_eq!(
                    hsp.num,
                    f[12].parse().unwrap(),
                    "redo {redo_index} linked {index} num"
                );
                assert_eq!(
                    hsp.evalue.to_bits(),
                    f[13].parse::<f64>().unwrap().to_bits(),
                    "redo {redo_index} linked {index} E-value"
                );
            }
            reap_by_evalue(&mut linked, 10.0);
            let reap_prefix = format!("D_HSP\t{reap_event}\treap_after\t");
            let reaped: Vec<_> = d_trace
                .lines()
                .filter(|line| line.starts_with(&reap_prefix))
                .collect();
            assert_eq!(
                linked.hsps.len(),
                reaped.len(),
                "redo {redo_index} reap count"
            );
            for (index, (hsp, row)) in linked.hsps.iter().zip(reaped.iter()).enumerate() {
                let f: Vec<_> = row.split('\t').collect();
                assert_eq!(
                    hsp.hsp.score,
                    f[10].parse().unwrap(),
                    "redo {redo_index} reaped {index} score"
                );
                assert_eq!(
                    hsp.evalue.to_bits(),
                    f[13].parse::<f64>().unwrap().to_bits(),
                    "redo {redo_index} reaped {index} E-value"
                );
            }
            // NCBI reference: core/blast_kappa.c:3671-3736
            // ```c
            // s_HitlistReapContained(...);
            // s_HitlistEvaluateAndPurge(&best_score, &best_evalue, ...);
            // s_HSPListNormalizeScores(hsp_list, kbp->Lambda, kbp->logK,
            //         localScalingFactor);
            // s_ComputeNumIdentities(...);
            // BlastCompo_HeapWouldInsert(..., best_evalue, best_score, ...);
            // BlastCompo_HeapInsert(..., hsp_list, best_evalue, best_score, ...);
            // ```
            let heap_trace = std::fs::read_to_string(format!(
                "{root}/kappa_traceback_20260924/multi_query_20260924_default.tsv"
            ))
            .unwrap();
            let would_rows: Vec<_> = heap_trace
                .lines()
                .filter(|line| line.starts_with("K_TRACE_HEAP_WOULD\t"))
                .collect();
            let would: Vec<_> = would_rows[redo_index].split('\t').collect();
            assert_eq!(
                linked.best_evalue.to_bits(),
                would[2].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(linked.hsps[0].hsp.score, would[3].parse().unwrap());
            assert_eq!(would[9], "1");
            let bits =
                normalize_postredo_scores(&mut linked, gapped.lambda / 32.0, gapped.k.ln(), 32.0);
            let heap_rows: Vec<_> = heap_trace
                .lines()
                .filter(|line| line.starts_with("K_TRACE_HEAP_HSP\t"))
                .collect();
            let expected_heap: Vec<_> = heap_rows
                .iter()
                .filter(|row| row.split('\t').nth(7).unwrap() == redo_index.to_string())
                .collect();
            assert_eq!(linked.hsps.len(), expected_heap.len());
            let code = GeneticCode::try_from_id(1).unwrap();
            let mut target = TargetTranslation::new(&subject_nt, &code);
            let mut identity_by_hsp = std::collections::HashMap::new();
            let mut alignment = redone.alignments_by_query[context as usize].as_deref();
            while let Some(value) = alignment {
                let identity = postredo_num_ident(
                    value,
                    query_infos[context as usize].seq.data(),
                    &mut target,
                )
                .unwrap();
                identity_by_hsp.insert(
                    (
                        value.frame,
                        value.query_start,
                        value.query_end,
                        value.match_start,
                        value.match_end,
                        value.score,
                    ),
                    identity,
                );
                alignment = value.next.as_deref();
            }
            for (index, ((hsp, bit_score), row)) in linked
                .hsps
                .iter()
                .zip(bits.iter())
                .zip(expected_heap.iter())
                .enumerate()
            {
                let f: Vec<_> = row.split('\t').collect();
                assert_eq!(
                    hsp.hsp.score,
                    f[3].parse().unwrap(),
                    "normalized HSP {index} score"
                );
                assert_eq!(
                    bit_score.to_bits(),
                    f[4].parse::<f64>().unwrap().to_bits(),
                    "normalized HSP {index} bit score"
                );
                assert_eq!(
                    hsp.evalue.to_bits(),
                    f[5].parse::<f64>().unwrap().to_bits(),
                    "normalized HSP {index} E-value"
                );
                let original_score: i32 =
                    reaped[index].split('\t').nth(10).unwrap().parse().unwrap();
                let key = (
                    i32::from(hsp.hsp.frame),
                    hsp.hsp.q_start,
                    hsp.hsp.q_end,
                    hsp.hsp.s_start,
                    hsp.hsp.s_end,
                    original_score,
                );
                assert_eq!(
                    identity_by_hsp[&key],
                    f[6].parse::<usize>().unwrap(),
                    "normalized HSP {index} identities"
                );
            }
            assert!(params.gapping_params.context.get().is_none());
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:224-273;
    // c++/src/algo/blast/core/blast_hits_priv.h:65-70
    // ```c
    // if (CONTAINED_IN_HSP(...) && CONTAINED_IN_HSP(...) &&
    //     hsp1->score <= hsp2->score) hsp1 = Blast_HSPFree(hsp1);
    // ```
    #[test]
    fn postredo_containment_keeps_ncbi_score_frame_and_endpoint_order() {
        let hsp = |frame, score, q_start, q_end, s_start, s_end| GappedHsp {
            frame,
            score,
            q_start,
            q_end,
            s_start,
            s_end,
            q_gapped_start: 0,
            s_gapped_start: 0,
        };
        let mut rows = vec![
            (0, hsp(1, 100, 10, 50, 20, 60)),
            (0, hsp(1, 99, 10, 50, 20, 60)),
            (0, hsp(-1, 99, 15, 45, 25, 55)),
            (1, hsp(1, 99, 15, 45, 25, 55)),
            (0, hsp(1, 101, 15, 45, 25, 55)),
            (0, hsp(1, 98, 10, 50, 20, 61)),
        ];
        reap_contained_postredo_hsps(&mut rows);
        assert_eq!(
            rows.iter().map(|row| row.1.score).collect::<Vec<_>>(),
            vec![100, 99, 99, 101, 98]
        );
        assert_eq!(rows[0].1.q_start, 10);
        assert_eq!(rows[4].1.s_end, 61);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1924-1948;
    // core/gapinfo.h:44-65
    // ```c
    // status = BLAST_GappedAlignmentWithTraceback(..., &fence_hit);
    // if (status == 0) return s_NewAlignmentFromGapAlign(...);
    // eGapAlignDel = 0; eGapAlignSub = 3; eGapAlignIns = 6;
    // ```
    #[test]
    fn kappa_traceback_matches_ncbi_saved_inputs_and_edit_scripts() {
        fn decode(hex: &str) -> Vec<u8> {
            hex.as_bytes()
                .chunks_exact(2)
                .map(|pair| u8::from_str_radix(std::str::from_utf8(pair).unwrap(), 16).unwrap())
                .collect()
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1896-1957;
        // composition_adjustment/redo_alignment.c:1219-1254
        // ```c
        // status = BLAST_GappedAlignmentWithTraceback(..., &fence_hit);
        // Blast_AdjustScores(matrix, query_composition, ..., &subject_composition,
        //     ..., compo_adjust_mode, ...);
        // ```
        // Include both code-32 calls in the CLI-calibrated local-subject
        // oracle, including the weak alignment removed after redo.
        for (case, trace_path, matrix_path) in [
            (
                "seg_hard_query_20260924_default",
                "kappa_traceback_20260924/seg_hard_query_20260924_default.tsv",
                "kappa_composition_matrix_scores_20260924/seg_hard_query_20260924_default.tsv",
            ),
            (
                "multi_query_20260924_default",
                "kappa_traceback_20260924/multi_query_20260924_default.tsv",
                "kappa_composition_matrix_scores_20260924/multi_query_20260924_default.tsv",
            ),
            (
                "code32_local_subject",
                "code32_local_api_20260925_cli_calibrated/code32.trace",
                "code32_local_api_20260925_cli_calibrated/code32_composition.trace",
            ),
            (
                "heap_replacement_112_subjects",
                "kappa_heap_rejection_20260925/result_order_20260925/ncbi.trace",
                "kappa_heap_rejection_20260925/natural_mode2_20260925/composition.tsv",
            ),
        ] {
            let root = format!(
                "{}/../docs/evidence/tlosan_stage_d",
                env!("CARGO_MANIFEST_DIR")
            );
            let trace = std::fs::read_to_string(format!("{root}/{trace_path}")).unwrap();
            let matrix_trace = std::fs::read_to_string(format!("{root}/{matrix_path}")).unwrap();
            let matrices: Vec<_> = matrix_trace
                .lines()
                .filter(|line| line.starts_with("K_ADJUSTED\t"))
                .collect();
            let rows: Vec<_> = trace
                .lines()
                .filter(|line| {
                    [
                        "K_TRACE_ENTER\t",
                        "K_TRACE_QUERY\t",
                        "K_TRACE_SUBJECT\t",
                        "K_TRACE_RETURN\t",
                        "K_TRACE_EDIT\t",
                    ]
                    .iter()
                    .any(|prefix| line.starts_with(prefix))
                })
                .collect();
            assert_eq!(rows.len() % 5, 0);
            assert_eq!(rows.len() / 5, matrices.len());
            for (event, block) in rows.chunks_exact(5).enumerate() {
                let enter: Vec<_> = block[0].split('\t').collect();
                let query: Vec<_> = block[1].split('\t').collect();
                let subject: Vec<_> = block[2].split('\t').collect();
                let outcome: Vec<_> = block[3].split('\t').collect();
                let edit: Vec<_> = block[4].split('\t').collect();
                assert_eq!(
                    (enter[0], query[0], subject[0], outcome[0], edit[0]),
                    (
                        "K_TRACE_ENTER",
                        "K_TRACE_QUERY",
                        "K_TRACE_SUBJECT",
                        "K_TRACE_RETURN",
                        "K_TRACE_EDIT"
                    )
                );
                assert_eq!(enter[2].parse::<usize>().unwrap(), event);
                assert_eq!(query[1].parse::<usize>().unwrap(), event);
                assert_eq!(subject[1].parse::<usize>().unwrap(), event);
                assert_eq!(outcome[2].parse::<usize>().unwrap(), event);
                assert_eq!(edit[1].parse::<usize>().unwrap(), event);
                assert_eq!(enter[3], "41");
                assert_eq!(enter[11], "1");
                assert_eq!(enter[12], "32");
                let matrix_fields: Vec<_> = matrices[event].split('\t').collect();
                assert_eq!(matrix_fields[1].parse::<usize>().unwrap(), event);
                assert_eq!(matrix_fields.len(), 2 + BLASTAA_SIZE * BLASTAA_SIZE);
                let mut scores = [[0i32; BLASTAA_SIZE]; BLASTAA_SIZE];
                for (row, cells) in matrix_fields[2..].chunks_exact(BLASTAA_SIZE).enumerate() {
                    for (col, value) in cells.iter().enumerate() {
                        scores[row][col] = value.parse().unwrap();
                    }
                }
                let adjusted = AdjustedProteinMatrix { scores };
                let query_data = BlastCompoSequenceData::from_ncbistdaa(&decode(query[2]));
                let subject_data = BlastCompoSequenceData::from_ncbistdaa(&decode(subject[2]));
                assert_eq!(query_data.length, enter[6].parse().unwrap());
                assert_eq!(subject_data.length, enter[7].parse().unwrap());
                let query_range = BlastCompoSequenceRange {
                    begin: 0,
                    end: query_data.length,
                    context: 0,
                };
                let subject_range = BlastCompoSequenceRange {
                    begin: 0,
                    end: subject_data.length,
                    context: 1,
                };
                let mut scratch = GapAlignScratch::new();
                let observed = redo_one_alignment_from_local_starts(
                    enter[4].parse().unwrap(),
                    enter[5].parse().unwrap(),
                    &query_data,
                    &query_range,
                    &subject_data,
                    &subject_range,
                    EMatrixAdjustRule::UserSpecifiedRelEntropy,
                    Some(&adjusted),
                    ScoringMatrix::Blosum62,
                    enter[9].parse().unwrap(),
                    enter[10].parse().unwrap(),
                    enter[8].parse().unwrap(),
                    &mut scratch,
                )
                .unwrap()
                .unwrap();
                assert_eq!(outcome[3], "0", "{case}: event {event}");
                assert_eq!(
                    observed.score,
                    outcome[4].parse().unwrap(),
                    "{case}: event {event}"
                );
                assert_eq!(
                    observed.query_start,
                    outcome[5].parse().unwrap(),
                    "{case}: event {event}"
                );
                assert_eq!(
                    observed.query_end,
                    outcome[6].parse().unwrap(),
                    "{case}: event {event}"
                );
                assert_eq!(
                    observed.match_start,
                    outcome[7].parse().unwrap(),
                    "{case}: event {event}"
                );
                assert_eq!(
                    observed.match_end,
                    outcome[8].parse().unwrap(),
                    "{case}: event {event}"
                );
                let Some(BlastCompoAlignmentContext::EditScript(ref script)) = observed.context
                else {
                    panic!("TBLASTN redo returned no edit script");
                };
                assert_eq!(script.len(), edit[2].parse::<usize>().unwrap());
                let observed_ops: Vec<_> = script
                    .iter()
                    .map(|op| match op {
                        GapEditOp::Sub(n) => format!("3:{n}"),
                        GapEditOp::Del(n) => format!("0:{n}"),
                        GapEditOp::Ins(n) => format!("6:{n}"),
                    })
                    .collect();
                assert_eq!(observed_ops, edit[3..], "{case}: event {event}");
            }
        }
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:870-899;
    // c++/src/algo/blast/core/blast_hits.c:3071-3115,3404-3417;
    // c++/src/algo/blast/core/blast_hspstream.c:144-152,289-319;
    // c++/src/algo/blast/core/blast_kappa.c:3391-3425
    // ```c
    // BLAST_LinkHsps(..., hsp_list_out, ...);
    // s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
    // Blast_HSPResultsReverseSort(hsp_stream->results);
    // *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
    // while (BlastHSPStreamRead(hsp_stream, &localMatch) != kBlastHSPStream_Eof)
    //     entry->match = localMatch;
    // ```
    #[test]
    fn multi_subject_stage_c_lists_enter_kappa_in_ncbi_stream_order() {
        use crate::algorithm::tblastn::search_gapped::{
            preliminary_protein_hsps_in_ncbi_order, PreliminaryProfile,
        };
        use crate::algorithm::tblastn::stage_d_linking::{
            compare_preliminary_lists_for_kappa, link_preliminary_hsps, reap_by_evalue,
            LinkedHspList,
        };
        use crate::algorithm::tblastn::stage_d_stats::{
            local_parameters_for_call, LocalParameterCall, LocalParameterOptions,
        };
        use crate::config::ProteinScoringSpec;
        use crate::stats::spouge::lookup_protein_gumbel_params;
        use crate::stats::tables::{lookup_protein_params_gapped, lookup_protein_params_ungapped};
        use crate::utils::seg::SegParams;
        use std::collections::HashMap;

        let root = format!("{}/../docs/evidence", env!("CARGO_MANIFEST_DIR"));
        let read_fasta = |path: &str| -> Vec<(String, Vec<u8>)> {
            let mut records: Vec<(String, Vec<u8>)> = Vec::new();
            for line in std::fs::read_to_string(path).unwrap().lines() {
                if let Some(id) = line.strip_prefix('>') {
                    records.push((id.to_string(), Vec::new()));
                } else {
                    records.last_mut().unwrap().1.extend(line.bytes());
                }
            }
            records
        };
        let query = read_fasta(&format!("{root}/tlosan_stage_c/run_20260923/query.faa"))
            .remove(0)
            .1;
        let subjects = read_fasta(&format!("{root}/tlosan_stage_c/run_20260923/subjects.fna"));
        let total_nt_length: usize = subjects.iter().map(|(_, seq)| seq.len()).sum();
        assert_eq!(
            (query.len(), subjects.len(), total_nt_length),
            (120, 13, 4694)
        );
        let trace = std::fs::read_to_string(format!(
            "{root}/tlosan_stage_d/run_20260924/run_20260923_default.trace"
        ))
        .unwrap();
        let kappa_trace = std::fs::read_to_string(format!(
            "{root}/tlosan_stage_d/kappa_traceback_20260924/run_20260923_default.tsv"
        ))
        .unwrap();
        let mut event_oid = HashMap::new();
        for line in trace.lines().filter(|line| line.starts_with("D_LIST\t")) {
            let f: Vec<_> = line.split('\t').collect();
            let event: usize = f[1].parse().unwrap();
            if event <= 22 && f[2] == "link_before" && f[3] != "NULL" {
                event_oid.insert(event, f[3].parse::<usize>().unwrap());
            }
        }
        let mut before: HashMap<usize, Vec<Vec<&str>>> = HashMap::new();
        let mut after: HashMap<usize, Vec<Vec<&str>>> = HashMap::new();
        for line in trace.lines().filter(|line| line.starts_with("D_HSP\t")) {
            let f: Vec<_> = line.split('\t').collect();
            let event: usize = f[1].parse().unwrap();
            let Some(&oid) = event_oid.get(&event) else {
                continue;
            };
            match f[2] {
                "link_before" => before.entry(oid).or_default().push(f),
                "link_after" => after.entry(oid).or_default().push(f),
                _ => {}
            }
        }
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            (total_nt_length / 3) as i64,
        )
        .unwrap();
        let parameters = local_parameters_for_call(
            &[(query.len(), true)],
            total_nt_length,
            &[gapped],
            &[ungapped],
            LocalParameterOptions {
                expect_value: 10.0,
                do_sum_stats: true,
                max_intron_length: 0,
                gap_trigger_bits: 22.0,
                word_xdrop_bits: 7.0,
                scale_factor: 1.0,
                gumbel: Some(&gumbel),
            },
            LocalParameterCall::Initial {
                min_subject_length: 120,
                composition_based_stats: 2,
            },
        );
        assert_eq!(parameters.cutoffs[0].word_cutoff, 18);
        assert_eq!(parameters.cutoffs[0].hit_cutoff, 18);
        assert_eq!(parameters.cutoffs[0].word_xdrop, 16);
        let seg = SegParams::default();
        let word_xdrop = [parameters.cutoffs[0].word_xdrop];
        let word_cutoff = [parameters.cutoffs[0].word_cutoff];
        let hit_cutoff = [parameters.cutoffs[0].hit_cutoff];
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &word_xdrop,
            word_cutoff: &word_cutoff,
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &hit_cutoff,
            hsp_num_max: i32::MAX as usize,
        };
        let mut retained: Vec<(usize, LinkedHspList)> = Vec::new();
        for (oid, (_, subject)) in subjects.iter().enumerate() {
            let (preliminary, _) =
                preliminary_protein_hsps_in_ncbi_order(&[&query], subject, 1, profile).unwrap();
            let expected = before.get(&oid).map(Vec::as_slice).unwrap_or(&[]);
            assert_eq!(preliminary.len(), expected.len(), "OID {oid} prelink count");
            for ((context, hsp), row) in preliminary.iter().zip(expected) {
                assert_eq!(*context, row[4].parse().unwrap(), "OID {oid} context");
                assert_eq!(hsp.frame, row[5].parse().unwrap(), "OID {oid} frame");
                assert_eq!(hsp.q_start, row[6].parse().unwrap(), "OID {oid} qstart");
                assert_eq!(hsp.q_end, row[7].parse().unwrap(), "OID {oid} qend");
                assert_eq!(hsp.s_start, row[8].parse().unwrap(), "OID {oid} sstart");
                assert_eq!(hsp.s_end, row[9].parse().unwrap(), "OID {oid} send");
                assert_eq!(hsp.score, row[10].parse().unwrap(), "OID {oid} score");
                assert_eq!(
                    hsp.s_gapped_start,
                    row[16].parse().unwrap(),
                    "OID {oid} gapped start"
                );
            }
            if preliminary.is_empty() {
                continue;
            }
            let mut linked = link_preliminary_hsps(
                &preliminary,
                &[query.len() as i32],
                &parameters.lengths,
                subject.len() as i32,
                &[gapped],
                &gumbel,
                parameters.link.as_ref().unwrap(),
            )
            .unwrap();
            let expected_linked = after.get(&oid).map(Vec::as_slice).unwrap_or(&[]);
            assert_eq!(
                linked.hsps.len(),
                expected_linked.len(),
                "OID {oid} linked count"
            );
            for (hsp, row) in linked.hsps.iter().zip(expected_linked) {
                assert_eq!(
                    hsp.hsp.score,
                    row[10].parse().unwrap(),
                    "OID {oid} linked score"
                );
                assert_eq!(hsp.num, row[12].parse().unwrap(), "OID {oid} linked num");
                assert_eq!(
                    hsp.evalue.to_bits(),
                    row[13].parse::<f64>().unwrap().to_bits(),
                    "OID {oid} linked E-value"
                );
            }
            reap_by_evalue(&mut linked, parameters.prelim_evalue);
            if !linked.hsps.is_empty() {
                retained.push((oid, linked));
            }
        }
        assert_eq!(retained.len(), 11);
        retained.sort_by(|(oid_a, a), (oid_b, b)| {
            compare_preliminary_lists_for_kappa(*oid_a as i32, a, *oid_b as i32, b)
        });
        let expected_oids: Vec<_> = kappa_trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_HEAP_WOULD\t"))
            .map(|line| line.split('\t').nth(1).unwrap().parse::<usize>().unwrap())
            .collect();
        assert_eq!(
            retained.iter().map(|(oid, _)| *oid).collect::<Vec<_>>(),
            expected_oids
        );
        for (redo_index, (_, linked)) in retained.iter().enumerate() {
            let expected_prelim: Vec<_> = kappa_trace
                .lines()
                .filter(|line| line.starts_with(&format!("K_TRACE_PRELIM\t{redo_index}\t")))
                .collect();
            assert_eq!(linked.hsps.len(), expected_prelim.len());
            for (linked, row) in linked.hsps.iter().zip(expected_prelim) {
                let f: Vec<_> = row.split('\t').collect();
                assert_eq!(linked.hsp.score, f[3].parse().unwrap());
                assert_eq!(linked.hsp.frame, f[5].parse().unwrap());
                assert_eq!(linked.hsp.q_start, f[6].parse().unwrap());
                assert_eq!(linked.hsp.q_end, f[7].parse().unwrap());
                assert_eq!(linked.hsp.q_gapped_start, f[8].parse().unwrap());
                assert_eq!(linked.hsp.s_start, f[9].parse().unwrap());
                assert_eq!(linked.hsp.s_end, f[10].parse().unwrap());
                assert_eq!(linked.hsp.s_gapped_start, f[11].parse().unwrap());
            }
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3577-3736;
        // c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275,330-391,439-466
        // ```c
        // Blast_RedoOneMatch(alignments, redo_align_params, incoming_aligns, ...);
        // s_HSPListFromDistinctAlignments(...);
        // s_HitlistReapContained(hsp_list);
        // s_HitlistEvaluateAndPurge(&best_score, &best_evalue, ...);
        // s_HSPListNormalizeScores(hsp_list, kbp->Lambda, kbp->logK, localScalingFactor);
        // s_ComputeNumIdentities(...);
        // if (BlastCompo_HeapWouldInsert(...)) BlastCompo_HeapInsert(...);
        // while ((hsp_list = BlastCompo_HeapPop(heap)) != NULL) ...;
        // ```
        use crate::algorithm::tblastn::kappa_heap::{CompoHeap, CompoHeapRecord};
        use crate::algorithm::tblastn::stage_d_linking::score_compare;
        use crate::algorithm::tblastn::stage_d_results::{
            KappaHspPayload, KappaResultHitList, KappaResultList,
        };
        use crate::core::composition_adjustment::adjust_scores::{
            build_matrix_info, read_aa_composition,
        };
        use crate::core::composition_adjustment::redo_alignment::{
            build_query_word_hashes, BlastCompoGappingParams,
        };
        use crate::stats::tables::KarlinParams;
        use crate::utils::matrix::aa_char_to_ncbistdaa;
        use std::cell::Cell;
        let mode_trace = std::fs::read_to_string(format!(
            "{root}/tlosan_stage_d/kappa_mode2_20260924/run_20260923_default.tsv"
        ))
        .unwrap();
        let redo_rows: Vec<_> = mode_trace
            .lines()
            .filter(|line| line.contains("\tredo_enter\t"))
            .collect();
        let would_rows: Vec<_> = kappa_trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_HEAP_WOULD\t"))
            .collect();
        let heap_hsp_rows: Vec<_> = kappa_trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_HEAP_HSP\t"))
            .collect();
        assert_eq!(
            (redo_rows.len(), would_rows.len(), heap_hsp_rows.len()),
            (11, 11, 11)
        );
        let query_aa: Vec<_> = query.iter().copied().map(aa_char_to_ncbistdaa).collect();
        let query_infos = [BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData::from_ncbistdaa(&query_aa),
            composition: read_aa_composition(&query_aa),
            eff_search_space: parameters.lengths[0].eff_searchsp as f64,
            words: Some(build_query_word_hashes(&query_aa)),
        }];
        let scaled = KarlinParams {
            lambda: gapped.lambda / 32.0,
            ..gapped
        };
        let mut post_link = parameters.link.unwrap();
        post_link.cutoff_small_gap = 0;
        let mut scratch = GapAlignScratch::new();
        let mut workspace = BlastCompositionWorkspace::new_blosum62();
        let mut heap = CompoHeap::new(500, 0.002).unwrap();
        let mut postredo_by_oid = HashMap::new();
        let mut result = KappaResultHitList::new(500).unwrap();
        let result_trace = std::fs::read_to_string(format!(
            "{root}/tlosan_stage_d/kappa_result_order_20260925/run_20260923_default.tsv"
        ))
        .unwrap();
        let result_in: Vec<_> = result_trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_RESULT_UPDATE_IN\t"))
            .collect();
        let result_out: Vec<_> = result_trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_RESULT_UPDATE_OUT\t"))
            .collect();
        let code = GeneticCode::try_from_id(1).unwrap();
        for (redo_index, (oid, preliminary)) in retained.iter().enumerate() {
            let redo: Vec<_> = redo_rows[redo_index].split('\t').collect();
            let call: Vec<_> = kappa_trace
                .lines()
                .find(|line| line.starts_with(&format!("K_TRACE_ENTER\t{redo_index}\t")))
                .unwrap()
                .split('\t')
                .collect();
            // NCBI reference: core/blast_kappa.c:2352-2390,2418-2479:
            // s_GetAlignParams reads the initial hit cutoffs and scaled score block.
            let params = local_kappa_redo_params(
                ScoringMatrix::Blosum62,
                11,
                1,
                &[gapped],
                &[true],
                &parameters,
                query.len() as i32,
                BlastCompoAdjustMode::CompositionMatrixAdjust,
                false,
                10.0,
                true,
                25.0,
                local_extension_final_xdrop(15.0, 25.0, gapped.lambda).unwrap(),
            )
            .unwrap();
            assert_eq!(params.gapping_params.x_dropoff, call[8].parse().unwrap());
            assert_eq!(params.gapping_params.gap_open, call[9].parse().unwrap());
            assert_eq!(params.gapping_params.gap_extend, call[10].parse().unwrap());
            assert_eq!(params.cutoff_score, redo[11].parse().unwrap());
            assert_eq!(
                (gapped.lambda / 32.0).to_bits(),
                redo[4].parse::<f64>().unwrap().to_bits()
            );
            let preliminary_hsps: Vec<_> =
                preliminary.hsps.iter().map(|linked| linked.hsp).collect();
            let subject = &subjects[*oid].1;
            let redone = redo_preliminary_match(
                &preliminary_hsps,
                0,
                &query_infos,
                subject,
                1,
                &params,
                gapped.lambda / 32.0,
                ScoringMatrix::Blosum62,
                &mut scratch,
                &mut workspace,
            )
            .unwrap();
            let align = redone.alignments_by_query[0].as_deref().unwrap();
            assert!(align.next.is_none(), "OID {oid} extra alignment");
            let align_prefix = format!("K_ALIGN\t{}\tredone\t", redo[1]);
            let expected_align: Vec<_> = mode_trace
                .lines()
                .find(|line| line.starts_with(&align_prefix))
                .unwrap()
                .split('\t')
                .collect();
            assert_eq!(
                align.score,
                expected_align[5].parse().unwrap(),
                "OID {oid} redo score"
            );
            assert_eq!(
                align.matrix_adjust_rule as i32,
                expected_align[6].parse().unwrap()
            );
            assert_eq!(align.query_start, expected_align[8].parse().unwrap());
            assert_eq!(align.query_end, expected_align[9].parse().unwrap());
            assert_eq!(align.match_start, expected_align[10].parse().unwrap());
            assert_eq!(align.match_end, expected_align[11].parse().unwrap());
            assert_eq!(align.frame, expected_align[12].parse().unwrap());
            let mut input = vec![(
                0,
                GappedHsp {
                    frame: align.frame as i8,
                    score: align.score,
                    q_start: align.query_start,
                    q_end: align.query_end,
                    q_gapped_start: 0,
                    s_start: align.match_start,
                    s_end: align.match_end,
                    s_gapped_start: 0,
                },
            )];
            input.sort_by(|a, b| score_compare(&a.1, &b.1));
            reap_contained_postredo_hsps(&mut input);
            let mut postredo = link_preliminary_hsps(
                &input,
                &[query.len() as i32],
                &parameters.lengths,
                subject.len() as i32,
                &[scaled],
                &gumbel,
                &post_link,
            )
            .unwrap();
            reap_by_evalue(&mut postredo, 10.0);
            assert_eq!(postredo.hsps.len(), 1);
            let would: Vec<_> = would_rows[redo_index].split('\t').collect();
            assert_eq!(*oid, would[1].parse::<usize>().unwrap());
            assert_eq!(
                postredo.best_evalue.to_bits(),
                would[2].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(postredo.hsps[0].hsp.score, would[3].parse().unwrap());
            let candidate = CompoHeapRecord {
                best_evalue: postredo.best_evalue,
                best_score: postredo.hsps[0].hsp.score,
                subject_index: *oid as i32,
            };
            assert_eq!(heap.len(), would[4].parse().unwrap());
            assert_eq!(heap.capacity(), would[6].parse().unwrap());
            assert_eq!(
                heap.worst_evalue().to_bits(),
                would[8].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(heap.would_insert(candidate), would[9] == "1");
            let bits = normalize_postredo_scores(&mut postredo, scaled.lambda, gapped.k.ln(), 32.0);
            let mut target = TargetTranslation::new(subject, &code);
            let identities = postredo_num_ident(align, &query_aa, &mut target).unwrap();
            let heap_hsp: Vec<_> = heap_hsp_rows[redo_index].split('\t').collect();
            assert_eq!(*oid, heap_hsp[1].parse::<usize>().unwrap());
            assert_eq!(postredo.hsps[0].hsp.score, heap_hsp[3].parse().unwrap());
            assert_eq!(
                bits[0].to_bits(),
                heap_hsp[4].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(
                postredo.hsps[0].evalue.to_bits(),
                heap_hsp[5].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(identities, heap_hsp[6].parse().unwrap());
            assert_eq!(postredo.hsps[0].hsp.frame, heap_hsp[8].parse().unwrap());
            assert!(heap.insert(candidate).is_none());
            // NCBI c++/src/algo/blast/core/blast_kappa.c:305-358,3687-3713:
            // Blast_HSPInit(..., &editScript, &new_hsp);
            // s_HSPListNormalizeScores(...); s_ComputeNumIdentities(...);
            let Some(BlastCompoAlignmentContext::EditScript(script)) = align.context.as_ref()
            else {
                panic!("OID {oid} Kappa HSP lacks an edit script");
            };
            let report = KappaHspPayload {
                bit_score: bits[0],
                num_ident: i32::try_from(identities).unwrap(),
                edit_script: script.clone(),
                matrix_adjust_rule: align.matrix_adjust_rule,
            };
            postredo_by_oid.insert(*oid as i32, (postredo, vec![report]));
            assert!(params.gapping_params.context.get().is_none());
        }
        let pop_rows: Vec<_> = kappa_trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_HEAP_POP\t"))
            .collect();
        let mut result_index = 0;
        for row in pop_rows {
            let f: Vec<_> = row.split('\t').collect();
            let popped = heap.pop();
            if f[1] == "-1" {
                assert!(popped.is_none());
            } else {
                let popped = popped.unwrap();
                assert_eq!(popped.subject_index, f[1].parse().unwrap());
                assert_eq!(heap.len(), f[2].parse::<usize>().unwrap());
                assert_eq!(
                    popped.best_evalue.to_bits(),
                    f[3].parse::<f64>().unwrap().to_bits()
                );
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2494-2515;
                // core/blast_hits.c:3243-3297,3420-3437:
                // while ((hsp_list = BlastCompo_HeapPop(heap)) != NULL)
                //     Blast_HitListUpdate(hitlist, hsp_list);
                // Blast_HSPResultsReverseOrder(results);
                let expected_in: Vec<_> = result_in[result_index].split('\t').collect();
                assert_eq!(popped.subject_index, expected_in[1].parse().unwrap());
                assert_eq!(result.lists().len(), expected_in[2].parse().unwrap());
                assert_eq!(
                    result.worst_evalue().to_bits(),
                    expected_in[4].parse::<f64>().unwrap().to_bits()
                );
                assert_eq!(result.low_score(), expected_in[5].parse().unwrap());
                let (hsps, payloads) = postredo_by_oid.remove(&popped.subject_index).unwrap();
                assert_eq!(hsps.hsps[0].hsp.score, expected_in[8].parse().unwrap());
                let saved_hsp: Vec<_> = heap_hsp_rows
                    .iter()
                    .find(|row| row.split('\t').nth(1) == Some(expected_in[1]))
                    .unwrap()
                    .split('\t')
                    .collect();
                assert_eq!(
                    payloads[0].bit_score.to_bits(),
                    saved_hsp[4].parse::<f64>().unwrap().to_bits()
                );
                assert_eq!(payloads[0].num_ident, saved_hsp[6].parse().unwrap());
                result
                    .update(KappaResultList {
                        oid: popped.subject_index,
                        hsps,
                        payloads,
                    })
                    .unwrap();
                let expected_out: Vec<_> = result_out[result_index].split('\t').collect();
                assert_eq!(result.lists().len(), expected_out[3].parse().unwrap());
                assert_eq!(
                    result.worst_evalue().to_bits(),
                    expected_out[4].parse::<f64>().unwrap().to_bits()
                );
                assert_eq!(result.low_score(), expected_out[5].parse().unwrap());
                assert_eq!(
                    result
                        .lists()
                        .iter()
                        .map(|list| list.oid)
                        .collect::<Vec<_>>(),
                    expected_out[7..]
                        .iter()
                        .map(|value| value.parse::<i32>().unwrap())
                        .collect::<Vec<_>>()
                );
                result_index += 1;
            }
        }
        assert_eq!(result_index, result_in.len());
        let reverse_out: Vec<_> = result_trace
            .lines()
            .find(|line| line.starts_with("K_TRACE_RESULT_REVERSE_OUT\t0\t"))
            .unwrap()
            .split('\t')
            .collect();
        result.reverse_order();
        assert_eq!(
            result
                .lists()
                .iter()
                .map(|list| list.oid)
                .collect::<Vec<_>>(),
            reverse_out[3..]
                .iter()
                .map(|value| value.parse::<i32>().unwrap())
                .collect::<Vec<_>>()
        );
    }
}
