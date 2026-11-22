//! Source-derived TBLASTN Kappa redo parameters for local subjects.

use std::cell::Cell;

use anyhow::{ensure, Context, Result};

use super::stage_d_stats::LocalSubjectParameters;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::{
    build_matrix_info, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    BlastCompoAdjustMode, BlastCompoGappingParams, BlastRedoAlignParams,
};
use crate::stats::compute_blosum62_ideal_karlin_params;
use crate::stats::tables::KarlinParams;
use crate::utils::matrix::BLASTAA_SIZE;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2352-2384,2445,
// 2433-2479,3099-3110; composition_adjustment/redo_alignment.c:1014-1050:
// near_identical_cutoff = (1.74 * NCBIMATH_LN2)
// cutoff_s = do_link_hsps ?
// localScalingFactor = (compo_adjust_mode != eNoCompositionBasedStats)
//                      ? SCALING_FACTOR : 1.0;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:457-463:
pub(super) fn local_extension_final_xdrop(
    prelim_bits: f64,
    final_bits: f64,
    min_lambda: f64,
) -> Result<i32> {
    ensure!(min_lambda > 0.0, "TBLASTN gapped lambda must be positive");
    // The MAX is taken of doubles and then converted as the x86-64 build does (INT_MIN beyond
    // Int4), as every protein program converts it (`core::x_drop::gapped_x_dropoffs`).
    Ok(crate::core::x_drop::gapped_x_dropoffs(prelim_bits, final_bits, min_lambda).1)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1316-1324
/// Makes `workspace`, BLOSUM62's (`new_blosum62`), the composition workspace of `matrix` when the
/// search matrix is another one: from its joint probabilities (`BlastCompositionWorkspace::new`;
/// IDENTITY has none and is searched in mode 0, `args::search_settings`). Not inlined, so that the
/// caller keeps BLOSUM62's workspace in place.
#[cold]
#[inline(never)]
pub(super) fn use_matrix_workspace(
    workspace: &mut BlastCompositionWorkspace,
    matrix: ScoringMatrix,
) -> Result<()> {
    *workspace = BlastCompositionWorkspace::new(matrix)?;
    Ok(())
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2428,2450-2457
    // `hitParams->do_sum_stats` is the flag of the hit-saving parameters, which a
    // `-max_intron_length` of 1 to 4 clears while the option `-sum_stats` stays on
    // (blast_parameters.c:802-806, `core::link_hsps::translated_longest_intron`): the redo
    // links exactly when the search has linking parameters.
    let do_link_hsps = initial.link.is_some();
    let scale = if mode == BlastCompoAdjustMode::NoCompositionBasedStats {
        1.0
    } else {
        32.0
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2369-2378,2439-2448:
    // Scans all valid query contexts for the minimum gapped Lambda, whereas
    // `near_identical_cutoff` takes the first valid context's gapped Lambda
    // (1.74 times `NCBIMATH_LN2` divided by Lambda).
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2372-2384:
    let redo_xdrop =
        crate::core::x_drop::redo_x_dropoff(final_xdrop_bits, scaled_lambda, extension_final_xdrop);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2833-2852;
    // blast_kappa.c:2216-2231:
    // BLOSUM62's ideal block is computed once per process; another matrix's is computed as the
    // lookup table's query preparation computes it (`ideal_karlin_params_of`).
    let ideal_lambda = if matrix == ScoringMatrix::Blosum62 {
        compute_blosum62_ideal_karlin_params()
            .map_err(anyhow::Error::msg)
            .context("TBLASTN Kappa ideal score block")?
            .lambda
    } else {
        crate::algorithm::tblastx::lookup::ideal_karlin_params_of(matrix).lambda
    };
    let matrix_info = build_matrix_info(matrix, ideal_lambda / scale)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2438-2479;
    // composition_adjustment/redo_alignment.c:1028-1043:
    // Computes `near_identical_cutoff` from the first valid context's gapped Lambda
    // (1.74 times `NCBIMATH_LN2` divided by Lambda) and chooses `cutoff_s` according to
    // whether HSP linking is enabled.
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1037-1037
    Ok(BlastRedoAlignParams {
        matrix_info,
        gapping_params: BlastCompoGappingParams {
            gap_open: ((gap_open as f64) * scale).round() as i32,
            gap_extend: ((gap_extend as f64) * scale).round() as i32,
            x_dropoff: redo_xdrop,
            context: Cell::new(None),
        },
        compo_adjust_mode: mode,
        alphsize: BLASTAA_SIZE as i32,
        composition_test_index: i32::from(unified_p),
        log_k: 0.0,
        score_divisor: scale,
        smith_waterman: false,
        near_identical_cutoff: (1.74 * std::f64::consts::LN_2) / first_valid_scaled_lambda,
        position_based: false,
        ccat_query_length: max_query_length,
        query_is_translated: false,
        subject_is_translated: true,
        cutoff_score: if do_link_hsps {
            ((initial.hit_cutoff_min as f64) * scale) as i32
        } else {
            1
        },
        do_link_hsps,
    })
}

#[cfg(test)]
#[path = "stage_d_kappa_params_tests.rs"]
mod tests;
