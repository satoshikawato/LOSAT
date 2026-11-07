//! `RPSRescalePssm`: the traceback's PSSM of a profile for one user query without
//! composition-based statistics, the profile's scores scaled by the ratio of the ungapped Lambda
//! of the query's composition on the profile to BLOSUM62's (both at the database's scale).

use anyhow::{bail, Result};

use super::blast_nint;
use crate::core::blast_stat::composition::calc_lambda_of_score_probs;
use crate::core::blast_util::BLASTAA_SIZE;

/// `BLAST_SCORE_MIN` (`INT2_MIN`).
const BLAST_SCORE_MIN: i32 = i16::MIN as i32;
/// `BLAST_SCORE_RANGE_MAX`: the size of `RPSRescalePssm`'s array of score probabilities.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:56
/// ```c
/// #define BLAST_SCORE_RANGE_MAX   (BLAST_SCORE_MAX - BLAST_SCORE_MIN) /**< maximum allowed range of BLAST scores. */
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4710
/// ```c
///     scoreArray = (double *)malloc(BLAST_SCORE_RANGE_MAX * sizeof(double));
/// ```
const BLAST_SCORE_RANGE_MAX: i64 = i16::MAX as i64 - i16::MIN as i64;
/// `AMINOACID_TO_NCBISTDAA['X']`.
const X_LETTER: usize = 21;
/// BLOSUM62's ungapped Lambda (`RPSfindUngappedLambda`: the first row of NCBI's table).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:258-259
/// ```c
/// static array_of_8 blosum62_values[BLOSUM62_VALUES_MAX] = {
///     {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2, 0.623757, 4.964660, 4.964660},
/// ```
pub(crate) const BLOSUM62_UNGAPPED_LAMBDA: f64 = 0.3176;

/// `Blast_FillResidueProbability`: each letter's share of the query, X left out of both.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4577-4600
/// ```c
///     denominator = length;
///     for(i = 0; i < BLASTAA_SIZE; i++)
///         frequency[i] = 0;
///
///     for(i = 0; i < length; i++) {
///         if (sequence[i] != AMINOACID_TO_NCBISTDAA['X'])
///             frequency[sequence[i]]++;
///         else
///             denominator--;
///     }
///
///     for(i = 0; i < BLASTAA_SIZE; i++) {
///         if (frequency[i] == 0)
///             resProb[i] = 0.0;
///         else
///             resProb[i] = ((double) frequency[i]) /((double) denominator);
///     }
/// ```
pub(crate) fn residue_probabilities(sequence: &[u8]) -> [f64; BLASTAA_SIZE] {
    let mut frequency = [0i32; BLASTAA_SIZE];
    let mut denominator = sequence.len() as i32;
    for &letter in sequence {
        if usize::from(letter) != X_LETTER {
            frequency[usize::from(letter)] += 1;
        } else {
            denominator -= 1;
        }
    }
    let mut probabilities = [0.0f64; BLASTAA_SIZE];
    for (probability, &count) in probabilities.iter_mut().zip(&frequency) {
        if count != 0 {
            *probability = f64::from(count) / f64::from(denominator);
        }
    }
    probabilities
}

/// `RPSRescalePssm` of the profile `rows` (`db_seq_length` rows of `alphabet` scores, the
/// database's columns) for a query of residue probabilities `query_probabilities`: the rescaled
/// rows into `out`, or `false` where NCBI returns NULL (no ungapped Lambda; the list is skipped).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4649-4691 (`RPSFillScores`)
/// ```c
///     minScore = maxScore = 0;
///     for (i = 0; i < matrixLength; i++) {
///         for (j = 0 ; j < alphabet_size; j++) {
///             if (j == AMINOACID_TO_NCBISTDAA['X'])
///                 continue;
///             if ((matrix[i][j] > BLAST_SCORE_MIN) &&
///                 (matrix[i][j] < minScore))
///                 minScore = matrix[i][j];
///             if (matrix[i][j] > maxScore)
///                 maxScore = matrix[i][j];
///         }
///     }
///
///     return_sfp->obs_min = minScore;
///     return_sfp->obs_max = maxScore;
///     memset(scoreArray, 0, (maxScore - minScore + 1) * sizeof(double));
///
///     return_sfp->sprob = &(scoreArray[-minScore]); /*center around 0*/
///     recipLength = 1.0 / (double) matrixLength;
///     for(i = 0; i < matrixLength; i++) {
///         for (j = 0; j < alphabet_size; j++) {
///             if (j == AMINOACID_TO_NCBISTDAA['X'])
///                 continue;
///             if(matrix[i][j] >= minScore)
///                 return_sfp->sprob[matrix[i][j]] += recipLength *
///                                                 queryProbArray[j];
///         }
///     }
///
///     return_sfp->score_avg = 0;
///     for(i = minScore; i <= maxScore; i++)
///         return_sfp->score_avg += i * return_sfp->sprob[i];
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4719-4757 (`RPSRescalePssm`)
/// ```c
///     initialUngappedLambda = RPSfindUngappedLambda(sbp->name);
///     ASSERT(initialUngappedLambda > 0.0);
///     scaledInitialUngappedLambda = initialUngappedLambda / scalingFactor;
///     correctUngappedLambda = Blast_KarlinLambdaNR(return_sfp,
///                                               scaledInitialUngappedLambda);
///     ...
///     if(correctUngappedLambda == -1.0)
///         return NULL;
///
///     finalLambda = correctUngappedLambda/scaledInitialUngappedLambda;
///     ...
///     for (index = 0; index < db_seq_length; index++) {
///         for (inner_index = 0; inner_index < alphabet_size; inner_index++) {
///             if (posMatrix[index][inner_index] <= BLAST_SCORE_MIN ||
///                 inner_index == AMINOACID_TO_NCBISTDAA['X']) {
///                 returnMatrix[index][inner_index] =
///                     posMatrix[index][inner_index];
///             }
///             else {
///                temp = ((double)(posMatrix[index][inner_index])) * finalLambda;
///                returnMatrix[index][inner_index] = (Int4)BLAST_Nint(temp);
///            }
///         }
///         for (; inner_index < BLASTAA_SIZE; inner_index++) {
///            returnMatrix[index][inner_index] = BLAST_SCORE_MIN;
///         }
///     }
/// ```
/// A profile whose scores span more than NCBI's array (`BLAST_SCORE_RANGE_MAX`, which NCBI's
/// `memset` overruns) is LOSAT's refusal. `calc_lambda_of_score_probs` is NCBI's
/// `Blast_KarlinLambdaNR` on the observed range as given, with this function's average.
pub(crate) fn rps_rescale_pssm(
    scale: f64,
    query_probabilities: &[f64; BLASTAA_SIZE],
    rows: &[[i32; BLASTAA_SIZE]],
    alphabet: usize,
    out: &mut Vec<[i32; BLASTAA_SIZE]>,
) -> Result<bool> {
    let (mut min_score, mut max_score) = (0i32, 0i32);
    for row in rows {
        for (j, &score) in row.iter().enumerate().take(alphabet) {
            if j == X_LETTER {
                continue;
            }
            if score > BLAST_SCORE_MIN && score < min_score {
                min_score = score;
            }
            if score > max_score {
                max_score = score;
            }
        }
    }
    if i64::from(max_score) - i64::from(min_score) + 1 > BLAST_SCORE_RANGE_MAX {
        bail!(
            "a profile whose scores span more than {BLAST_SCORE_RANGE_MAX} values (NCBI BLAST+ writes past its array) is not supported by LOSAT's RPSBLAST"
        );
    }
    let mut probs = vec![0.0f64; (max_score - min_score + 1) as usize];
    let recip_length = 1.0 / rows.len() as f64;
    for row in rows {
        for (j, &score) in row.iter().enumerate().take(alphabet) {
            if j == X_LETTER {
                continue;
            }
            if score >= min_score {
                probs[(score - min_score) as usize] += recip_length * query_probabilities[j];
            }
        }
    }
    let scaled_initial_lambda = BLOSUM62_UNGAPPED_LAMBDA / scale;
    let correct_lambda =
        calc_lambda_of_score_probs(probs, min_score, max_score, scaled_initial_lambda);
    if correct_lambda == -1.0 {
        return Ok(false);
    }
    let final_lambda = correct_lambda / scaled_initial_lambda;
    out.clear();
    out.extend(rows.iter().map(|row| {
        let mut rescaled = [BLAST_SCORE_MIN; BLASTAA_SIZE];
        for (j, &score) in row.iter().enumerate().take(alphabet) {
            rescaled[j] = if score <= BLAST_SCORE_MIN || j == X_LETTER {
                score
            } else {
                blast_nint(f64::from(score) * final_lambda) as i32
            };
        }
        rescaled
    }));
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residue_probabilities_leave_x_out() {
        // A A C X: X is not counted, A is 2/3 and C 1/3.
        let probabilities = residue_probabilities(&[1, 1, 3, 21]);
        assert_eq!(probabilities[1], 2.0 / 3.0);
        assert_eq!(probabilities[3], 1.0 / 3.0);
        assert_eq!(probabilities[21], 0.0);
        assert_eq!(residue_probabilities(&[21, 21]), [0.0; BLASTAA_SIZE]);
    }

    // A profile scaled by 2 whose scores read on the query's composition have BLOSUM62's
    // Lambda at scale 2 rescales to itself; the X column and the columns past the alphabet
    // follow NCBI's rules.
    #[test]
    fn rescaling_keeps_x_and_marks_the_columns_past_the_alphabet() {
        let mut row = [0i32; BLASTAA_SIZE];
        for (j, score) in row.iter_mut().enumerate() {
            *score = if j % 3 == 0 { 8 } else { -6 };
        }
        row[0] = BLAST_SCORE_MIN;
        row[21] = 99;
        let rows = vec![row; 4];
        let mut probabilities = [0.0; BLASTAA_SIZE];
        for (j, p) in probabilities.iter_mut().enumerate().skip(1).take(20) {
            *p = if j % 2 == 0 { 0.06 } else { 0.04 };
        }
        let mut out = Vec::new();
        assert!(rps_rescale_pssm(100.0, &probabilities, &rows, 26, &mut out).unwrap());
        assert_eq!(out.len(), 4);
        assert_eq!(out[0][0], BLAST_SCORE_MIN);
        assert_eq!(out[0][21], 99);
        assert_eq!((out[0][26], out[0][27]), (BLAST_SCORE_MIN, BLAST_SCORE_MIN));
        // Every scaled score has the sign of its source.
        assert!(out[0][3] > 0 && out[0][1] < 0);
        // No negative expected score: no Lambda, NULL.
        let positive = vec![[5i32; BLASTAA_SIZE]; 2];
        assert!(!rps_rescale_pssm(1.0, &probabilities, &positive, 28, &mut out).unwrap());
    }
}
