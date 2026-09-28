//! Source-derived TBLASTN Kappa redo parameters for local subjects.

use std::cell::Cell;

use anyhow::{ensure, Context, Result};

use super::stage_d_stats::LocalSubjectParameters;
use crate::config::ScoringMatrix;
use crate::core::blast_stat::compute_blosum62_ideal_karlin_params;
use crate::core::composition_adjustment::adjust_scores::build_matrix_info;
use crate::core::composition_adjustment::redo_alignment::{
    BlastCompoAdjustMode, BlastCompoGappingParams, BlastRedoAlignParams,
};
use crate::stats::tables::KarlinParams;
use crate::utils::matrix::BLASTAA_SIZE;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2352-2384,
// 2433-2479,3099-3110; composition_adjustment/redo_alignment.c:1014-1050:
// gapping_params->gap_open = scoring->gap_open;
// gapping_params->gap_extend = scoring->gap_extend;
// gapping_params->x_dropoff = (Int4)
//     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
//         extendParams->gap_x_dropoff_final);
// near_identical_cutoff = (1.74 * NCBIMATH_LN2)
//                         / context->sbp->kbp_gap[index]->Lambda;
// cutoff_s = do_link_hsps ?
//     (int)(hitParams->cutoff_score_min * context->localScalingFactor) : 1;
// localScalingFactor = (compo_adjust_mode != eNoCompositionBasedStats)
//                      ? SCALING_FACTOR : 1.0;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:457-463:
// params->gap_x_dropoff = (Int4)(options->gap_x_dropoff*NCBIMATH_LN2/min_lambda);
// params->gap_x_dropoff_final = (Int4)
//     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2/min_lambda,
//         params->gap_x_dropoff);
#[allow(dead_code)] // Calculated at extension-parameter setup before Kappa.
pub(super) fn local_extension_final_xdrop(
    prelim_bits: f64,
    final_bits: f64,
    min_lambda: f64,
) -> Result<i32> {
    ensure!(min_lambda > 0.0, "TBLASTN gapped lambda must be positive");
    let prelim = ((prelim_bits * std::f64::consts::LN_2) / min_lambda) as i32;
    let final_drop = ((final_bits * std::f64::consts::LN_2) / min_lambda) as i32;
    Ok(final_drop.max(prelim))
}

#[allow(dead_code)] // The public local-subject path stays gated through Stage E.
pub(super) fn local_kappa_redo_params(
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    gapped: &[KarlinParams],
    valid_contexts: &[bool],
    initial: &LocalSubjectParameters,
    max_query_length: i32,
    mode: BlastCompoAdjustMode,
    unified_p: bool,
    expect_value: f64,
    do_sum_stats: bool,
    final_xdrop_bits: f64,
    extension_final_xdrop: i32,
) -> Result<BlastRedoAlignParams> {
    ensure!(
        !gapped.is_empty(),
        "TBLASTN Kappa needs a valid query context"
    );
    ensure!(
        max_query_length > 0,
        "TBLASTN Kappa query length must be positive"
    );
    ensure!(
        gapped.len() == valid_contexts.len(),
        "TBLASTN Kappa context state mismatch"
    );
    ensure!(
        matrix == ScoringMatrix::Blosum62,
        "TBLASTN Kappa ideal score block is only ported for BLOSUM62"
    );
    ensure!(
        valid_contexts.iter().any(|&valid| valid),
        "TBLASTN Kappa needs a valid query context"
    );
    ensure!(
        gapped
            .iter()
            .zip(valid_contexts)
            .all(|(p, &valid)| !valid || p.lambda > 0.0),
        "TBLASTN Kappa lambda must be positive"
    );
    ensure!(
        !do_sum_stats || initial.link.is_some(),
        "TBLASTN Kappa linking parameters are missing"
    );
    let scale = if mode == BlastCompoAdjustMode::NoCompositionBasedStats {
        1.0
    } else {
        32.0
    };
    // NCBI c++/src/algo/blast/core/blast_kappa.c:2369-2378,2439-2448:
    // for (i=0; i<num_queries; i++) if (kbp_gap[i] != NULL &&
    //     kbp_gap[i]->Lambda < min_lambda) min_lambda = kbp_gap[i]->Lambda;
    // for (index=first_context; index<=last_context; ++index)
    //     if (contexts[index].is_valid) { near_identical_cutoff =
    //         (1.74*NCBIMATH_LN2)/kbp_gap[index]->Lambda; break; }
    let active = gapped
        .iter()
        .zip(valid_contexts)
        .filter(|(_, &valid)| valid);
    let min_lambda = active
        .clone()
        .fold(f64::MAX, |min, (p, _)| min.min(p.lambda));
    let first_valid_lambda = active.map(|(p, _)| p.lambda).next().unwrap();
    let scaled_lambda = min_lambda / scale;
    let first_valid_scaled_lambda = first_valid_lambda / scale;
    // NCBI c++/src/algo/blast/core/blast_kappa.c:2372-2384:
    // gapping_params->x_dropoff = (Int4)
    //     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2/min_lambda,
    //         extendParams->gap_x_dropoff_final);
    let redo_xdrop = ((final_xdrop_bits * std::f64::consts::LN_2) / scaled_lambda) as i32;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2833-2852;
    // blast_kappa.c:2216-2231:
    // Blast_ScoreBlkKbpIdealCalc(sbp);
    // self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
    let ideal = compute_blosum62_ideal_karlin_params()
        .map_err(anyhow::Error::msg)
        .context("TBLASTN Kappa ideal score block")?;
    let matrix_info = build_matrix_info(matrix, ideal.lambda / scale)?;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2438-2479;
    // composition_adjustment/redo_alignment.c:1028-1043:
    // near_identical_cutoff = (1.74 * NCBIMATH_LN2) / kbp_gap[index]->Lambda;
    // cutoff_s = do_link_hsps ?
    //     (int)(hitParams->cutoff_score_min * localScalingFactor) : 1;
    // params->ccat_query_length = ccat_query_length;
    Ok(BlastRedoAlignParams {
        matrix_info,
        gapping_params: BlastCompoGappingParams {
            gap_open: ((gap_open as f64) * scale).round() as i32,
            gap_extend: ((gap_extend as f64) * scale).round() as i32,
            decline_align: 0,
            x_dropoff: redo_xdrop.max(extension_final_xdrop),
            context: Cell::new(None),
        },
        compo_adjust_mode: mode,
        alphsize: BLASTAA_SIZE as i32,
        composition_test_index: i32::from(unified_p),
        unified_p,
        log_k: 0.0,
        score_divisor: scale,
        restricted_alignment: false,
        smith_waterman: false,
        is_same_adjustment: false,
        near_identical_cutoff: (1.74 * std::f64::consts::LN_2) / first_valid_scaled_lambda,
        position_based: false,
        re_matrix_adjustment_pseudocounts: 20,
        ccat_query_length: max_query_length,
        query_is_translated: false,
        subject_is_translated: true,
        cutoff_score: if do_sum_stats {
            ((initial.hit_cutoff_min as f64) * scale) as i32
        } else {
            1
        },
        cutoff_evalue: expect_value,
        do_link_hsps: do_sum_stats,
    })
}

#[cfg(test)]
mod tests {
    use super::super::stage_d_stats::LocalLinkParameters;
    use super::*;

    // NCBI reference: c++/src/algo/blast/core/blast_kappa.c:2369-2382,
    // 2439-2448:
    // for (i=0; i<num_queries; i++) if (kbp_gap[i] != NULL &&
    //     kbp_gap[i]->Lambda < min_lambda) min_lambda = kbp_gap[i]->Lambda;
    // for (index=first_context; index<=last_context; ++index)
    //     if (contexts[index].is_valid) { near_identical_cutoff =
    //         (1.74*NCBIMATH_LN2)/kbp_gap[index]->Lambda; break; }
    #[test]
    fn first_valid_lambda_differs_from_minimum_redo_lambda() {
        let gapped = [
            KarlinParams {
                lambda: 0.3,
                ..KarlinParams::default()
            },
            KarlinParams {
                lambda: 0.267,
                ..KarlinParams::default()
            },
            KarlinParams {
                lambda: 0.2,
                ..KarlinParams::default()
            },
        ];
        let initial = LocalSubjectParameters {
            lengths: Vec::new(),
            cutoffs: Vec::new(),
            link: Some(LocalLinkParameters {
                gap_decay_rate: 0.1,
                gap_size: 40,
                overlap_size: 9,
                longest_intron: 40,
                cutoff_small_gap: 0,
            }),
            hit_cutoff_min: 25,
            word_cutoff_min: 25,
            prelim_evalue: 50.0,
        };
        let params = local_kappa_redo_params(
            ScoringMatrix::Blosum62,
            11,
            1,
            &gapped,
            &[false, true, true],
            &initial,
            120,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
            false,
            10.0,
            true,
            30.0,
            70,
        )
        .unwrap();
        assert_eq!(
            params.near_identical_cutoff.to_bits(),
            ((1.74 * std::f64::consts::LN_2) / (0.267 / 32.0)).to_bits()
        );
        assert_eq!(
            params.gapping_params.x_dropoff,
            (((30.0 * std::f64::consts::LN_2) / (0.2 / 32.0)) as i32).max(70)
        );
        assert_eq!(params.cutoff_score, 800);
    }

    // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:457-463:
    // params->gap_x_dropoff_final = (Int4)
    //     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2/min_lambda,
    //         params->gap_x_dropoff);
    #[test]
    fn extension_final_xdrop_uses_option_bits_and_previous_cutoff() {
        let lambda = 0.267;
        let prelim = ((35.0 * std::f64::consts::LN_2) / lambda) as i32;
        assert_eq!(
            local_extension_final_xdrop(35.0, 25.0, lambda).unwrap(),
            prelim
        );
    }
}
