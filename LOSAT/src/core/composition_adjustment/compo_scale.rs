//! The core of the composition-based statistics of mode 1 ("standard composition-based
//! statistics", `eCompoScaleOldMatrix`) for a position-specific score matrix: NCBI's
//! position-based `Blast_MatrixInfo` (`s_MatrixInfoInit` with `positionBased`) and
//! `Blast_CompositionBasedStats` for a PSSM (series PS, item PS4b; DECISIONS 11 and 13 of line PS).
//!
//! The parts that do not depend on where the rows come from (the score probabilities' range, the
//! lambda ratio, the rounding of a scaled row) are separate functions, so that the square-matrix
//! entry of NCBI (`s_GetMatrixScoreProbs`, `s_ScaleSquareMatrix`) can sit beside the PSSM's
//! (`pssm_score_probs`, `scale_pssm`) when blastp's mode 1 (series SC, item SC5) moves here.
//! BLASTP's square-matrix scaling stays in `adjust_scores.rs` until then.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1104-1158

use crate::core::psi::PsiScoreBlk;
use crate::stats::karlin_calc::calc_lambda_of_score_probs;
use crate::utils::glibc_fma_libm::GlibcLibm;
use crate::utils::matrix::{ncbistdaa, BLASTAA_SIZE};

/// `COMPO_SCORE_MIN` (composition_constants.h:43: `#define COMPO_SCORE_MIN INT2_MIN`).
pub(crate) const COMPO_SCORE_MIN: i32 = i16::MIN as i32;

/// Which redo runs the steps that the X-drop redo and the Smith-Waterman redo share (this
/// module's scaling, `RedoMatch`'s subject range, reap and evaluation): one source, one instance
/// for each, so that each instance has one caller and stays inlined where it is, as before the
/// second redo (line PS `DECISIONS.md` 11 and 63; series PR, item PS11c-3). `XDropRedo` is
/// psiblast's and rpsblast's redo (`Blast_RedoOneMatch`), `SwRedo` rpsblast's `-use_sw_tback`
/// (`Blast_RedoOneMatchSmithWaterman`), `FrameRedo` rpstblastn's (`Blast_RedoOneMatch` with a
/// translated subject: series PR, item PS12b).
pub(crate) trait RedoInstance {
    /// NCBI's `subject_is_translated` of the redo (`s_GetAlignParams`): the subject is the six
    /// frames of a translated sequence, each alignment in its frame.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2422-2423
    const SUBJECT_IS_TRANSLATED: bool = false;
}

/// The X-drop redo (`RedoInstance`).
pub(crate) enum XDropRedo {}

/// The Smith-Waterman redo (`RedoInstance`).
pub(crate) enum SwRedo {}

/// The X-drop redo of a translated subject (`RedoInstance`): rpstblastn's, whose subject is the
/// user query's frames.
pub(crate) enum FrameRedo {}

impl RedoInstance for XDropRedo {}
impl RedoInstance for SwRedo {}
impl RedoInstance for FrameRedo {
    const SUBJECT_IS_TRANSLATED: bool = true;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:43-54
const TRUE_CHAR_POSITIONS: [usize; 20] = [
    1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
];
/// Whether a letter is one of the 20 true amino acids (`alphaConvert[j] >= 0`).
const TRUE_AMINO_ACID: [bool; BLASTAA_SIZE] = {
    let mut table = [false; BLASTAA_SIZE];
    let mut i = 0;
    while i < TRUE_CHAR_POSITIONS.len() {
        table[TRUE_CHAR_POSITIONS[i]] = true;
        i += 1;
    }
    table
};

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:80
const MAXIMUM_X_SCORE: f64 = -1.0;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1099-1101
const LAMBDA_RATIO_LOWER_BOUND: f64 = 0.5;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:603-604,607-608
const POS_EPSILON: f64 = 0.0001;

/// The position-based `Blast_MatrixInfo` of a PSSM query (`positionBased`, `rows` = the query
/// length, 28 columns).
#[derive(Debug, Clone)]
pub(crate) struct PosMatrixInfo {
    /// `startFreqRatios`: the frequency ratios each row of the redo's PSSM is scaled from.
    pub(crate) start_freq_ratios: Vec<[f64; BLASTAA_SIZE]>,
    /// `startMatrix`: the PSSM scaled by the redo's factor; its score range sets the score
    /// probabilities and its `*` column stays in every adjusted PSSM.
    pub(crate) start_matrix: Vec<[i32; BLASTAA_SIZE]>,
    /// `ungappedLambda`: `kbp_psi[0]->Lambda` of the factor-1 PSSM over the scaling factor.
    pub(crate) ungapped_lambda: f64,
    /// `kbp_gap_psi[0]->K` that `_PSIUpdateLambdaK` leaves in the redo's copy of the score
    /// block: the e-values and bit scores of the redo read it (design question Q1, settled by an
    /// LD_PRELOAD probe of NCBI 2.17.0: `BLAST_SpougeStoE` gets this K and the gapped Lambda of
    /// the PSSM file over the scaling factor).
    pub(crate) gapped_k: f64,
}

impl PosMatrixInfo {
    /// `s_MatrixInfoInit` for a position-based search: the query `query` (NCBIstdaa, as the
    /// search holds it: O as X), the PSSM's frequency ratios `freq_ratios` (28 per position, 0.0
    /// where the PSSM has none), the score block `sbp` (the search matrix's), the scaling factor
    /// `scale_factor` (`SCALING_FACTOR`, 32) and the search's `kbp_psi[0]` (`kbp_psi`: the PSSM
    /// file's ungapped values, or `kbp_std`'s), which stays where the factor-1 PSSM has no score
    /// frequencies (`Blast_ScoreFreqNew` NULL, review S1 of PS4b). `None` when the factor-1 PSSM
    /// has score frequencies but no Karlin-Altschul parameters (NCBI goes on with a Lambda of -1).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2215-2227
    pub(crate) fn new(
        query: &[u8],
        freq_ratios: &[[f64; BLASTAA_SIZE]],
        sbp: &PsiScoreBlk,
        scale_factor: f64,
        kbp_psi: crate::core::psi::KarlinBlk,
    ) -> Option<Self> {
        let start_freq_ratios = pos_based_start_freq_ratios(query, freq_ratios, sbp);
        let scaled =
            crate::core::psi::scale_pos_matrix(query, freq_ratios, sbp, scale_factor, kbp_psi);
        if !(scaled.kbp_psi.lambda > 0.0 && scaled.kbp_psi.k > 0.0) {
            return None;
        }
        Some(Self {
            start_freq_ratios,
            start_matrix: scaled.start_matrix,
            ungapped_lambda: scaled.kbp_psi.lambda / scale_factor,
            gapped_k: scaled.gapped_k,
        })
    }
}

/// `s_GetPosBasedStartFreqRatios`: each row starts as the matrix's frequency ratios of the query
/// letter; a cell takes the PSSM's ratio over the letter's standard probability where both
/// letters are standard, the column is neither `*` nor X, and the PSSM's ratio is above 1e-4.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:606-631
fn pos_based_start_freq_ratios(
    query: &[u8],
    freq_ratios: &[[f64; BLASTAA_SIZE]],
    sbp: &PsiScoreBlk,
) -> Vec<[f64; BLASTAA_SIZE]> {
    let standard_prob = crate::stats::karlin_calc::compute_std_aa_composition();
    query
        .iter()
        .zip(freq_ratios)
        .map(|(&residue, numerator)| {
            let mut row = sbp.freq_ratios[residue as usize];
            for (j, ratio) in row.iter_mut().enumerate() {
                if standard_prob[residue as usize] > POS_EPSILON
                    && standard_prob[j] > POS_EPSILON
                    && j != ncbistdaa::STOP as usize
                    && j != ncbistdaa::X as usize
                    && numerator[j] > POS_EPSILON
                {
                    *ratio = numerator[j] / standard_prob[j];
                }
            }
            row
        })
        .collect()
}

/// The score probabilities of a matrix (`scoreArray`): `probs[s - obs_min]` is the probability
/// of score `s`, for `obs_min..=obs_max`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ScoreProbs {
    pub(crate) probs: Vec<f64>,
    pub(crate) obs_min: i32,
    pub(crate) obs_max: i32,
}

/// `s_GetScoreRange`: the smallest score above `COMPO_SCORE_MIN` (at most 0) and the largest
/// score (at least 0) of the true amino-acid columns of `rows`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:806-828
pub(crate) fn score_range<C: RedoInstance>(rows: &[[i32; BLASTAA_SIZE]]) -> (i32, i32) {
    let (mut min_score, mut max_score) = (0, 0);
    for row in rows {
        for &jcol in &TRUE_CHAR_POSITIONS {
            let score = row[jcol];
            if score < min_score && score > COMPO_SCORE_MIN {
                min_score = score;
            }
            if score > max_score {
                max_score = score;
            }
        }
    }
    (min_score, max_score)
}

/// `s_GetPssmScoreProbs`: the score probabilities of the PSSM `matrix` against a subject of the
/// letter probabilities `subject_prob`, each position weighing `1/rows`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:921-940
pub(crate) fn pssm_score_probs<C: RedoInstance>(
    matrix: &[[i32; BLASTAA_SIZE]],
    subject_prob: &[f64; BLASTAA_SIZE],
) -> ScoreProbs {
    let (obs_min, obs_max) = score_range::<C>(matrix);
    let mut probs = vec![0.0; (obs_max - obs_min + 1) as usize];
    let one_pos_frac = 1.0 / (matrix.len() as f64);
    for row in matrix {
        for &jcol in &TRUE_CHAR_POSITIONS {
            if row[jcol] >= obs_min {
                probs[(row[jcol] - obs_min) as usize] += one_pos_frac * subject_prob[jcol];
            }
        }
    }
    ScoreProbs {
        probs,
        obs_min,
        obs_max,
    }
}

/// The lambda ratio of `Blast_CompositionBasedStats`: the Lambda of the score probabilities
/// (`s_CalcLambda`, started from `ungapped_lambda`) over `ungapped_lambda`, at most 1 without the
/// composition p-value (`pValueAdjustment` 0, psiblast's case) and at least 0.5.
pub(crate) fn lambda_ratio<C: RedoInstance>(
    probs: ScoreProbs,
    ungapped_lambda: f64,
    p_value_adjustment: bool,
) -> f64 {
    let correct_ungapped_lambda =
        calc_lambda_of_score_probs(probs.probs, probs.obs_min, probs.obs_max, ungapped_lambda);
    let mut ratio = correct_ungapped_lambda / ungapped_lambda;
    if !p_value_adjustment {
        ratio = ratio.min(1.0);
    }
    ratio.max(LAMBDA_RATIO_LOWER_BOUND)
}

/// `Blast_FreqRatioToScore` of one ratio: `COMPO_SCORE_MIN` for 0, `log(ratio)/Lambda` otherwise.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:403-415
fn freq_ratio_to_score<C: RedoInstance>(ratio: f64, lambda: f64) -> f64 {
    if ratio == 0.0 {
        COMPO_SCORE_MIN as f64
    } else {
        ratio.glibc_ln() / lambda
    }
}

/// `s_CalcXScore` of one row: the row's average score against the true amino acids of the
/// probabilities `probs`, at most -1.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:616-645
fn calc_x_score<C: RedoInstance>(row: &[f64; BLASTAA_SIZE], probs: &[f64; BLASTAA_SIZE]) -> f64 {
    let mut score_ix = 0.0;
    for j in 0..BLASTAA_SIZE {
        if TRUE_AMINO_ACID[j] {
            score_ix += row[j] * probs[j];
        }
    }
    score_ix.min(MAXIMUM_X_SCORE)
}

/// `s_RoundScoreMatrix` of one score: `INT_MIN` below it, the nearest integer (`Nint`) otherwise.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:707-738
pub(crate) fn round_score<C: RedoInstance>(score: f64) -> i32 {
    if score < i32::MIN as f64 {
        i32::MIN
    } else {
        (score + if score >= 0.0 { 0.5 } else { -0.5 }) as i64 as i32
    }
}

/// `s_ScalePSSM`: each row of `matrix` from the row of `start_freq_ratios` at the scale
/// `lambda`, X (and O) the row's average against the subject's letter probabilities `col_prob`,
/// U the row's C, `*` the start matrix's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1025-1050
pub(crate) fn scale_pssm<C: RedoInstance>(
    matrix: &mut [[i32; BLASTAA_SIZE]],
    start_freq_ratios: &[[f64; BLASTAA_SIZE]],
    start_matrix: &[[i32; BLASTAA_SIZE]],
    col_prob: &[f64; BLASTAA_SIZE],
    lambda: f64,
) {
    for ((out, ratios), start) in matrix.iter_mut().zip(start_freq_ratios).zip(start_matrix) {
        let mut row = [0.0; BLASTAA_SIZE];
        for (score, &ratio) in row.iter_mut().zip(ratios) {
            *score = freq_ratio_to_score::<C>(ratio, lambda);
        }
        let x_score = calc_x_score::<C>(&row, col_prob);
        row[ncbistdaa::X as usize] = x_score;
        row[ncbistdaa::U as usize] = row[ncbistdaa::C as usize];
        row[ncbistdaa::O as usize] = x_score;
        for (cell, &score) in out.iter_mut().zip(&row) {
            *cell = round_score::<C>(score);
        }
        out[ncbistdaa::STOP as usize] = start[ncbistdaa::STOP as usize];
    }
}

/// `Blast_CompositionBasedStats` for a PSSM (`positionBased`, no composition p-value): `matrix`
/// becomes the PSSM of `info` scaled for a subject of the letter probabilities `subject_prob`;
/// returns the lambda ratio (`LambdaRatio`).
pub(crate) fn pssm_composition_based_stats<C: RedoInstance>(
    matrix: &mut [[i32; BLASTAA_SIZE]],
    info: &PosMatrixInfo,
    subject_prob: &[f64; BLASTAA_SIZE],
) -> f64 {
    let probs = pssm_score_probs::<C>(&info.start_matrix, subject_prob);
    let ratio = lambda_ratio::<C>(probs, info.ungapped_lambda, false);
    // The ratio is at least 0.5 (`LambdaRatioLowerBound`), so NCBI's `*LambdaRatio > 0` holds.
    let scaled_lambda = info.ungapped_lambda / ratio;
    scale_pssm::<C>(
        matrix,
        &info.start_freq_ratios,
        &info.start_matrix,
        subject_prob,
        scaled_lambda,
    );
    ratio
}

/// Whether every frequency ratio is at most `kEpsilon` (1e-4), which `PsiBlastSetupScoreBlock`
/// warns of when the composition-based statistics are on (a PSSM without ratios counts as all
/// zeros: its flag starts true).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:151-169
pub(crate) fn freq_ratios_all_zeros(freq_ratios: &[[f64; BLASTAA_SIZE]]) -> bool {
    !freq_ratios
        .iter()
        .flatten()
        .any(|&ratio| ratio > POS_EPSILON)
}

#[cfg(test)]
#[path = "compo_scale_tests.rs"]
mod tests;
