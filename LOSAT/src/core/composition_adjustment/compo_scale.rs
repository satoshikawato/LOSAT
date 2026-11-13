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
//! ```c
//! int
//! Blast_CompositionBasedStats(int ** matrix, double * LambdaRatio,
//!                             const Blast_MatrixInfo * ss,
//!                             const double queryProb[], const double resProb[],
//!                             double (*calc_lambda)(double*,int,int,double),
//!                             int pValueAdjustment)
//! {
//!     ...
//!     if (ss->positionBased) {
//!         out_of_memory =
//!             s_GetPssmScoreProbs(&scoreArray, &obs_min, &obs_max,
//!                                 ss->startMatrix, ss->rows, resProb);
//!     } else {
//!         ...
//!     }
//!     ...
//!     correctUngappedLambda =
//!         calc_lambda(scoreArray, obs_min, obs_max, ss->ungappedLambda);
//!     ...
//!     *LambdaRatio = correctUngappedLambda / ss->ungappedLambda;
//!     if (0 == pValueAdjustment)
//!       *LambdaRatio = MIN(1, *LambdaRatio);
//!     *LambdaRatio = MAX(*LambdaRatio, LambdaRatioLowerBound);
//!
//!     if (*LambdaRatio > 0) {
//!         double scaledLambda = ss->ungappedLambda/(*LambdaRatio);
//!         if (ss->positionBased) {
//!             s_ScalePSSM(matrix, ss->rows, ss->cols, ss->startFreqRatios,
//!                         ss->startMatrix, resProb, scaledLambda);
//!         } else {
//!             ...
//!         }
//!     }
//!     free(scoreArray);
//!
//!     return 0;
//! }
//! ```

use crate::core::psi::PsiScoreBlk;
use crate::stats::karlin_calc::calc_lambda_of_score_probs;
use crate::utils::matrix::{ncbistdaa, BLASTAA_SIZE};

/// `COMPO_SCORE_MIN` (composition_constants.h:43: `#define COMPO_SCORE_MIN INT2_MIN`).
pub(crate) const COMPO_SCORE_MIN: i32 = i16::MIN as i32;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:43-54
// ```c
// /**positions of true characters in protein alphabet*/
// static int trueCharPositions[COMPO_NUM_TRUE_AA] =
// {1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22};
//
// static int alphaConvert[COMPO_LARGEST_ALPHABET] =
//   {(-1), 0, (-1),  4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15,
//    16, 19, 17, (-1), 18, (-1), (-1), (-1), (-1), (-1)};
// ```
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
// ```c
// static const double kMaximumXscore = -1.0;
// ```
const MAXIMUM_X_SCORE: f64 = -1.0;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1099-1101
// ```c
// #define LambdaRatioLowerBound 0.5
// ```
const LAMBDA_RATIO_LOWER_BOUND: f64 = 0.5;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:607-608
// ```c
//     const double kPosEpsilon = 0.0001;   /* values below this cutoff
//                                             are treated specially */
// ```
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2216-2227
    /// ```c
    ///     if (self->positionBased) {
    ///         status = s_GetPosBasedStartFreqRatios(self->startFreqRatios,
    ///                                               queryBlk->length,
    ///                                               queryBlk->sequence,
    ///                                               matrixName,
    ///                                               sbp->psi_matrix->freq_ratios);
    ///         if (status == 0) {
    ///             status = s_ScalePosMatrix(self->startMatrix, matrixName,
    ///                                       sbp->psi_matrix->freq_ratios,
    ///                                       queryBlk->sequence,
    ///                                       queryBlk->length, sbp, scale_factor);
    ///             self->ungappedLambda = sbp->kbp_psi[0]->Lambda / scale_factor;
    ///         }
    /// ```
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:610-631
/// ```c
///     stdFreqRatios = _PSIMatrixFrequencyRatiosNew(matrixName);
///     ...
///     for (i = 0;  i < numPositions;  i++) {
///         for (j = 0;  j < BLASTAA_SIZE;  j++) {
///             returnRatios[i][j] = stdFreqRatios->data[query[i]][j];
///         }
///     }
///     stdFreqRatios = _PSIMatrixFrequencyRatiosFree(stdFreqRatios);
///
///     standardProb = BLAST_GetStandardAaProbabilities();
///     ...
///     /*reverse multiplication done in posit.c*/
///     for (i = 0;  i < numPositions;  i++) {
///         for (j = 0;  j < BLASTAA_SIZE;  j++) {
///             if ((standardProb[query[i]] > kPosEpsilon) &&
///                 (standardProb[j] > kPosEpsilon) &&
///                 (j != eStopChar) && (j != eXchar) &&
///                 (startNumerator[i][j] > kPosEpsilon)) {
///                 returnRatios[i][j] = startNumerator[i][j] / standardProb[j];
///             }
///         }
///     }
/// ```
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:806-827
/// ```c
///     minScore = maxScore = 0;
///     for (irow = 0;  irow < rows;  irow++) {
///         for (aa = 0;  aa < COMPO_NUM_TRUE_AA;  aa++) {
///             jcol = trueCharPositions[aa];
///             if (matrix[irow][jcol] < minScore &&
///                 matrix[irow][jcol] > COMPO_SCORE_MIN)
///                 minScore = matrix[irow][jcol];
///             if (matrix[irow][jcol] > maxScore)
///                 maxScore = matrix[irow][jcol];
///         }
///     }
///     *obs_min = minScore;
///     *obs_max = maxScore;
/// ```
pub(crate) fn score_range(rows: &[[i32; BLASTAA_SIZE]]) -> (i32, i32) {
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:922-940
/// ```c
///     s_GetScoreRange(obs_min, obs_max, matrix, rows);
///     minScore = *obs_min;
///     range = *obs_max - *obs_min + 1;
///     *scoreProb = calloc(range, sizeof(double));
///     ...
///     sprob = &((*scoreProb)[-(*obs_min)]); /*center around 0*/
///     onePosFrac = 1.0/ ((double) rows);
///     for (irow = 0;  irow < rows;  irow++) {
///         for (aa = 0;  aa < COMPO_NUM_TRUE_AA;  aa++) {
///             jcol = trueCharPositions[aa];
///             if (matrix[irow][jcol] >= minScore) {
///                 sprob[matrix[irow][jcol]] +=
///                     onePosFrac * subjectProbArray[jcol];
///             }
///         }
///     }
/// ```
pub(crate) fn pssm_score_probs(
    matrix: &[[i32; BLASTAA_SIZE]],
    subject_prob: &[f64; BLASTAA_SIZE],
) -> ScoreProbs {
    let (obs_min, obs_max) = score_range(matrix);
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
pub(crate) fn lambda_ratio(
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
/// ```c
///             if (0.0 == matrix[i][j]) {
///                 matrix[i][j] = COMPO_SCORE_MIN;
///             } else {
///                 matrix[i][j] = log(matrix[i][j])/Lambda;
///             }
/// ```
fn freq_ratio_to_score(ratio: f64, lambda: f64) -> f64 {
    if ratio == 0.0 {
        COMPO_SCORE_MIN as f64
    } else {
        ratio.ln() / lambda
    }
}

/// `s_CalcXScore` of one row: the row's average score against the true amino acids of the
/// probabilities `probs`, at most -1.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:616-645
/// ```c
///     for (j = 0;  j < alphsize;  j++) {
///         if (alphaConvert[j] >= 0) {
///             /* If the column corresponds to a true amino acid */
///             score_iX += M[j * incM] * probs[j];
///         }
///     }
///     return score_iX;
/// ...
///     return MIN( s_CalcAvgScore(M, alphsize, incM, probs), kMaximumXscore);
/// ```
fn calc_x_score(row: &[f64; BLASTAA_SIZE], probs: &[f64; BLASTAA_SIZE]) -> f64 {
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
/// ```c
/// static long Nint(double x)
/// {
///     x += (x >= 0. ? 0.5 : -0.5);
///     return (long)x;
/// }
/// ...
///             if (floatScoreMatrix[p][c] < INT_MIN) {
///                 matrix[p][c] = INT_MIN;
///             } else {
///                 matrix[p][c] = Nint(floatScoreMatrix[p][c]);
///             }
/// ```
pub(crate) fn round_score(score: f64) -> i32 {
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
/// ```c
///     for (p = 0;  p < rows;  p++) {
///         double Xscore;
///         memcpy(row, freq_ratios[p], cols * sizeof(double));
///
///         Blast_FreqRatioToScore(row_matrix, 1, cols, Lambda);
///         row[eXchar] = Xscore = s_CalcXScore(row, cols, 1, col_prob);
///         /* use Cysteine score for Selenocysteine */
///         row[eSelenocysteine] = row[eCchar];
///         if (cols > eOchar) {
///             row[eOchar] = Xscore;
///         }
///         s_RoundScoreMatrix(&matrix[p], 1, cols, row_matrix);
///
///         matrix[p][eStopChar] = start_matrix[p][eStopChar];
///     }
/// ```
pub(crate) fn scale_pssm(
    matrix: &mut [[i32; BLASTAA_SIZE]],
    start_freq_ratios: &[[f64; BLASTAA_SIZE]],
    start_matrix: &[[i32; BLASTAA_SIZE]],
    col_prob: &[f64; BLASTAA_SIZE],
    lambda: f64,
) {
    for ((out, ratios), start) in matrix.iter_mut().zip(start_freq_ratios).zip(start_matrix) {
        let mut row = [0.0; BLASTAA_SIZE];
        for (score, &ratio) in row.iter_mut().zip(ratios) {
            *score = freq_ratio_to_score(ratio, lambda);
        }
        let x_score = calc_x_score(&row, col_prob);
        row[ncbistdaa::X as usize] = x_score;
        row[ncbistdaa::U as usize] = row[ncbistdaa::C as usize];
        row[ncbistdaa::O as usize] = x_score;
        for (cell, &score) in out.iter_mut().zip(&row) {
            *cell = round_score(score);
        }
        out[ncbistdaa::STOP as usize] = start[ncbistdaa::STOP as usize];
    }
}

/// `Blast_CompositionBasedStats` for a PSSM (`positionBased`, no composition p-value): `matrix`
/// becomes the PSSM of `info` scaled for a subject of the letter probabilities `subject_prob`;
/// returns the lambda ratio (`LambdaRatio`).
pub(crate) fn pssm_composition_based_stats(
    matrix: &mut [[i32; BLASTAA_SIZE]],
    info: &PosMatrixInfo,
    subject_prob: &[f64; BLASTAA_SIZE],
) -> f64 {
    let probs = pssm_score_probs(&info.start_matrix, subject_prob);
    let ratio = lambda_ratio(probs, info.ungapped_lambda, false);
    // The ratio is at least 0.5 (`LambdaRatioLowerBound`), so NCBI's `*LambdaRatio > 0` holds.
    let scaled_lambda = info.ungapped_lambda / ratio;
    scale_pssm(
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
/// ```c++
///     bool freq_ratios_all_zeros = true;
///     ...
///                 if ((*freq_ratios)(r,c) > kEpsilon) {
///                     freq_ratios_all_zeros = false;
///                 }
/// ```
pub(crate) fn freq_ratios_all_zeros(freq_ratios: &[[f64; BLASTAA_SIZE]]) -> bool {
    !freq_ratios
        .iter()
        .flatten()
        .any(|&ratio| ratio > POS_EPSILON)
}

#[cfg(test)]
mod tests {
    use super::*;

    // `s_GetScoreRange` keeps 0 in the range and skips `COMPO_SCORE_MIN`.
    #[test]
    fn the_score_range_holds_zero_and_skips_the_minimum() {
        let mut row = [5; BLASTAA_SIZE];
        row[1] = COMPO_SCORE_MIN;
        row[3] = -7;
        // Not a true amino acid: not read.
        row[21] = -100;
        assert_eq!(score_range(&[row]), (-7, 5));
        assert_eq!(score_range(&[[2; BLASTAA_SIZE]]), (0, 2));
    }

    // `s_RoundScoreMatrix`: `Nint` rounds half away from zero; below `INT_MIN` is `INT_MIN`.
    #[test]
    fn scores_round_as_ncbis_nint() {
        assert_eq!(round_score(2.5), 3);
        assert_eq!(round_score(-2.5), -3);
        assert_eq!(round_score(-2.49), -2);
        assert_eq!(round_score(COMPO_SCORE_MIN as f64), COMPO_SCORE_MIN);
        assert_eq!(round_score(-1.0e12), i32::MIN);
    }

    // `s_ScalePSSM` on one hand-checked row: ratio 0 gives `COMPO_SCORE_MIN`, X and O the
    // average (capped at -1), U the C score, `*` the start matrix's.
    #[test]
    fn a_scaled_row_is_ncbis() {
        let mut ratios = [1.0; BLASTAA_SIZE];
        ratios[0] = 0.0;
        ratios[1] = std::f64::consts::E; // A: 1/lambda
        ratios[3] = std::f64::consts::E * std::f64::consts::E; // C: 2/lambda
        let mut start = [0; BLASTAA_SIZE];
        start[ncbistdaa::STOP as usize] = -77;
        let mut probs = [0.0; BLASTAA_SIZE];
        probs[1] = 0.5;
        probs[3] = 0.5;
        let mut matrix = [[0; BLASTAA_SIZE]];
        scale_pssm(&mut matrix, &[ratios], &[start], &probs, 0.5);
        let row = matrix[0];
        assert_eq!(row[0], COMPO_SCORE_MIN);
        assert_eq!(row[1], 2);
        assert_eq!(row[3], 4);
        assert_eq!(row[ncbistdaa::U as usize], 4);
        // The average is 3 > -1: X and O are -1.
        assert_eq!(row[ncbistdaa::X as usize], -1);
        assert_eq!(row[ncbistdaa::O as usize], -1);
        assert_eq!(row[ncbistdaa::STOP as usize], -77);
        assert_eq!(row[5], 0);
    }

    fn gunzip(name: &str) -> Vec<u8> {
        use std::io::Read;
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(name);
        let mut out = Vec::new();
        flate2::read::GzDecoder::new(std::fs::read(path).unwrap().as_slice())
            .read_to_end(&mut out)
            .unwrap();
        out
    }

    // Design question Q1 of PS4b, by an LD_PRELOAD probe of NCBI 2.17.0 (`psiblast -in_pssm
    // tk.asn -db paja -comp_based_stats 1 -num_threads 1`, `PS/tmp/q1/`): in
    // `Blast_RedoAlignmentCore_MT`, `Blast_KarlinBlkUngappedCalc` of the factor-1 PSSM gives
    // Lambda 0.32206865605351226 and K 0.13750086820961227, and every `BLAST_SpougeStoE` of the
    // redo gets Lambda 0.267/32 and K 0.042084934723385127 (= K * 0.041 / 0.13395614448848181).
    #[test]
    fn the_matrix_information_of_the_tk_pssm_is_ncbis() {
        let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(
            "pssm/ps4_tk_paja.asn.gz",
        ))
        .unwrap();
        let search = crate::algorithm::pssm_search::in_pssm::psi_search_of(pssm).unwrap();
        let query: Vec<u8> = search
            .query
            .seq()
            .iter()
            .map(|&c| crate::utils::matrix::aa_char_to_ncbistdaa(c))
            .collect();
        let sbp = crate::algorithm::pssm_search::in_pssm::score_block(
            crate::config::ScoringMatrix::Blosum62,
            crate::core::psi::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                h: 0.14,
            },
        )
        .unwrap();
        assert_eq!(sbp.ideal_k, 0.13395614448848181);
        let info = PosMatrixInfo::new(
            &query,
            &search.freq_ratios,
            &sbp,
            32.0,
            file_ungapped(&search),
        )
        .unwrap();
        assert_eq!(info.ungapped_lambda, 0.32206865605351226 / 32.0);
        assert_eq!(info.gapped_k, 0.042084934723385127);
        assert_eq!(info.start_matrix.len(), 398);
    }

    /// The ungapped block of a PSSM file (all three values present).
    fn file_ungapped(
        search: &crate::algorithm::pssm_search::PsiSearch,
    ) -> crate::core::psi::KarlinBlk {
        crate::core::psi::KarlinBlk {
            lambda: search.stats.lambda_ungapped.unwrap(),
            k: search.stats.kappa_ungapped.unwrap(),
            h: search.stats.h_ungapped.unwrap(),
        }
    }

    // Review S1 of PS4b: the tk PSSM with every frequency ratio 5e-5 (`ps4_tk_eps`) has
    // factor-1 scores of -24 to -17 for the 20 letters, so NCBI's `Blast_ScoreFreqNew` gives
    // NULL, `Blast_KarlinBlkUngappedCalc` returns 1 and the redo keeps the file's ungapped block:
    // an LD_PRELOAD probe of NCBI 2.17.0 (`psiblast -in_pssm tk_eps.asn -db paja`) saw Lambda
    // 0.31673860011845001 and K 0.13910107343338302 in the redo and K 0.042574710048176159 in
    // every redo e-value.
    #[test]
    fn a_factor_1_pssm_without_score_frequencies_keeps_the_files_block() {
        let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(
            "pssm/ps4_tk_eps.asn.gz",
        ))
        .unwrap();
        let search = crate::algorithm::pssm_search::in_pssm::psi_search_of(pssm).unwrap();
        let query: Vec<u8> = search
            .query
            .seq()
            .iter()
            .map(|&c| crate::utils::matrix::aa_char_to_ncbistdaa(c))
            .collect();
        let sbp = crate::algorithm::pssm_search::in_pssm::score_block(
            crate::config::ScoringMatrix::Blosum62,
            crate::core::psi::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                h: 0.14,
            },
        )
        .unwrap();
        let kbp_psi = file_ungapped(&search);
        assert_eq!(kbp_psi.lambda, 0.31673860011845001);
        assert_eq!(kbp_psi.k, 0.13910107343338302);
        let info = PosMatrixInfo::new(&query, &search.freq_ratios, &sbp, 32.0, kbp_psi).unwrap();
        assert_eq!(info.ungapped_lambda, 0.31673860011845001 / 32.0);
        assert_eq!(info.gapped_k, 0.042574710048176159);
    }

    #[test]
    fn all_zero_ratios_are_those_at_most_epsilon() {
        assert!(freq_ratios_all_zeros(&[]));
        assert!(freq_ratios_all_zeros(&[[0.0001; BLASTAA_SIZE]]));
        let mut row = [0.0; BLASTAA_SIZE];
        row[7] = 0.00011;
        assert!(!freq_ratios_all_zeros(&[row]));
    }
}
