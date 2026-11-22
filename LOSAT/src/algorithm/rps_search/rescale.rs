//! `RPSRescalePssm`: the traceback's PSSM of a profile for one user query without
//! composition-based statistics, the profile's scores scaled by the ratio of the ungapped Lambda
//! of the query's composition on the profile to BLOSUM62's (both at the database's scale).

use anyhow::{bail, Result};

use super::blast_nint;
use crate::core::blast_util::BLASTAA_SIZE;
use crate::stats::karlin_calc::calc_lambda_of_score_probs;

/// `BLAST_SCORE_MIN` (`INT2_MIN`).
const BLAST_SCORE_MIN: i32 = i16::MIN as i32;
/// `BLAST_SCORE_RANGE_MAX`: the size of `RPSRescalePssm`'s array of score probabilities.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:56
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4710
const BLAST_SCORE_RANGE_MAX: i64 = i16::MAX as i64 - i16::MIN as i64;
/// `AMINOACID_TO_NCBISTDAA['X']`.
const X_LETTER: usize = 21;
/// BLOSUM62's ungapped Lambda (`RPSfindUngappedLambda`: the first row of NCBI's table).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:258-259
pub(crate) const BLOSUM62_UNGAPPED_LAMBDA: f64 = 0.3176;

/// `Blast_FillResidueProbability`: each letter's share of the query, X left out of both.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4577-4600
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4719-4757 (`RPSRescalePssm`)
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
    // An alignment uses each row once: the best scores of the rows sum to the most it can score.
    // A sum above `INT4_MAX` overflows NCBI's `Int4` scores (undefined behavior in C).
    let bound: i64 = out
        .iter()
        .map(|row| i64::from(row.iter().copied().max().unwrap_or(0).max(0)))
        .sum();
    if bound > i64::from(i32::MAX) {
        bail!(
            "a profile whose rescaled scores can add up to {bound} (above 2147483647, where NCBI BLAST+'s 32-bit scores overflow) is not supported by LOSAT's RPSBLAST"
        );
    }
    Ok(true)
}

#[cfg(test)]
#[path = "rescale_tests.rs"]
mod tests;
