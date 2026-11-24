//! The stages of a PSSM made from conserved domains (DELTA-BLAST's CDD input): the port of NCBI's
//! `PSICreatePssmFromCDD` (`algo/blast/core/blast_psi.c`) and of its stages in
//! `blast_psi_priv.c` (`_PSIValidateCdMSA`, `_PSIComputeFrequenciesFromCDs`,
//! `_PSIComputeFreqRatiosFromCDs`, `_PSISaveCDDiagnostics`). The stages after the frequency
//! ratios (`_PSICreateAndScalePssmFromFrequencyRatios`) and the per-column arithmetic of the
//! frequency ratios are the alignment path's (`stages.rs`), in the same order.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:229-321

use super::stages::{
    column_freq_ratios, column_specific_pseudocounts, convert_freq_ratios_to_pssm, gapped_psi_k,
    scale_matrix, ColumnPseudocounts, InternalPssm, GAP, K_EPSILON, PSEUDO_MAX, X,
};
use super::{KarlinBlk, PsiDiagnostics, PsiError, PsiPssm, PsiScoreBlk};
use crate::utils::matrix::BLASTAA_SIZE;

/// The data of one aligned column of a conserved domain (`PSICdMsaCellData`): the domain's
/// weighted residue frequencies at the column (`wfreqs`, 28 letters) and the number of
/// independent observations there (`iobsr`).
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_psi.h:226-231
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CdColumn {
    pub(crate) wfreqs: [f64; BLASTAA_SIZE],
    pub(crate) iobsr: f64,
}

/// A cell of the conserved-domain alignment: not aligned, or aligned with the data of `columns`
/// at the index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CdCell(u32);

impl CdCell {
    /// A cell that no domain column is aligned to (`is_aligned` false).
    pub(crate) const UNALIGNED: Self = Self(u32::MAX);

    /// The cell aligned with `columns[index]`.
    pub(crate) fn aligned(index: usize) -> Self {
        Self(u32::try_from(index).expect("fewer than 2^32 aligned columns"))
    }

    fn index(self) -> Option<usize> {
        (self.0 != u32::MAX).then_some(self.0 as usize)
    }
}

/// The input of the CD engine (`PSICdMsa`): the query, and for each domain a row of cells, one per
/// query position, that name the domain's aligned columns.
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_psi.h:233-251
#[derive(Debug, Clone)]
pub(crate) struct PsiCdMsa {
    /// The query's letters in NCBIstdaa (`query`, `dimensions->query_length` of them).
    pub(crate) query: Vec<u8>,
    /// `msa[domain][position]`.
    pub(crate) rows: Vec<Vec<CdCell>>,
    /// The aligned columns that the cells name.
    pub(crate) columns: Vec<CdColumn>,
}

impl PsiCdMsa {
    fn column(&self, domain: usize, position: usize) -> Option<&CdColumn> {
        self.rows[domain][position]
            .index()
            .map(|index| &self.columns[index])
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:862-907
/// `_PSIValidateCdMSA`: no gap in the query, and every aligned column has data with at least
/// `kEpsilon` independent observations, no negative frequency and frequencies that sum to 1
/// within `kEpsilon`.
pub(crate) fn validate_cd_msa(msa: &PsiCdMsa) -> Result<(), PsiError> {
    if msa.query.contains(&GAP) {
        return Err(PsiError::GapInQuery);
    }
    for domain in 0..msa.rows.len() {
        for position in 0..msa.query.len() {
            let Some(data) = msa.column(domain, position) else {
                continue;
            };
            // `!(iobsr >= kEpsilon)` is not `iobsr < kEpsilon` for a NaN: C's test lets it pass.
            if data.iobsr < K_EPSILON {
                return Err(PsiError::BadProfile);
            }
            let mut sum = 0.0;
            for &frequency in &data.wfreqs {
                if frequency < 0.0 {
                    return Err(PsiError::BadProfile);
                }
                sum += frequency;
            }
            if (sum - 1.0).abs() > K_EPSILON {
                return Err(PsiError::BadProfile);
            }
        }
    }
    Ok(())
}

/// The weights of a CD alignment's columns (`_PSISequenceWeights`' `match_weights` and
/// `independent_observations`) and the standard residue probabilities.
struct CdWeights {
    match_weights: Vec<[f64; BLASTAA_SIZE]>,
    independent_observations: Vec<f64>,
    std_prob: [f64; BLASTAA_SIZE],
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1651-1745
/// `_PSIComputeFrequenciesFromCDs`: the weighted residue frequencies of each query position,
/// the domains' by their independent observations, plus the query's own letter when no domain
/// shows it, and the number of independent observations (at most 400).
fn compute_frequencies_from_cds(msa: &PsiCdMsa) -> CdWeights {
    let query_length = msa.query.len();
    let mut weights = CdWeights {
        match_weights: vec![[0.0; BLASTAA_SIZE]; query_length],
        independent_observations: vec![0.0; query_length],
        std_prob: crate::stats::karlin_calc::compute_std_aa_composition(),
    };
    // No domain aligned to the query: the weights stay zero.
    if msa.rows.is_empty() {
        return weights;
    }
    for position in 0..query_length {
        let mut total_observations = 0.0;
        let query_residue = msa.query[position] as usize;
        let mut sum_weights = [0.0f64; BLASTAA_SIZE];
        for domain in 0..msa.rows.len() {
            let Some(data) = msa.column(domain, position) else {
                continue;
            };
            total_observations += data.iobsr;
            for residue in 0..BLASTAA_SIZE {
                sum_weights[residue] += data.wfreqs[residue] * data.iobsr;
            }
        }
        if total_observations > 0.0
            && query_residue != X as usize
            && sum_weights[query_residue] == 0.0
        {
            sum_weights[query_residue] = 1.0;
            total_observations += 1.0;
        }
        if total_observations > 0.0 {
            for residue in 0..BLASTAA_SIZE {
                weights.match_weights[position][residue] =
                    sum_weights[residue] / total_observations;
            }
        }
        // MIN(400.0, total_observations)
        weights.independent_observations[position] = if 400.0 < total_observations {
            400.0
        } else {
            total_observations
        };
    }
    weights
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2185-2276
/// `_PSIComputeFreqRatiosFromCDs` with `pseudo_count` (the options' 0: the column-specific
/// pseudocounts).
fn compute_freq_ratios_from_cds(
    msa: &PsiCdMsa,
    weights: &CdWeights,
    sbp: &PsiScoreBlk,
    pseudo_count: i32,
    internal: &mut InternalPssm,
) -> Result<(), PsiError> {
    const K_ZERO_OBS_PSEUDO: f64 = 30.0;
    for position in 0..msa.query.len() {
        let query_is_x = msa.query[position] == X;
        let mut column_counts = 0.0;
        let mut observations = 0.0;
        if !query_is_x {
            // MAX(0.0, independent_observations[p] - 1.0)
            let observed = weights.independent_observations[position] - 1.0;
            observations = if 0.0 > observed { 0.0 } else { observed };
            column_counts = if pseudo_count == 0 {
                column_specific_pseudocounts(
                    &weights.match_weights[position],
                    sbp.background,
                    observations,
                )
            } else {
                f64::from(pseudo_count)
            };
        }
        let pseudo_weight = if column_counts >= PSEUDO_MAX {
            observations = 0.0;
            K_ZERO_OBS_PSEUDO
        } else {
            column_counts
        };
        column_freq_ratios(
            ColumnPseudocounts {
                query_is_x,
                pseudo_weight,
                observations,
            },
            &weights.match_weights[position],
            &weights.std_prob,
            sbp,
            false,
            &mut internal.freq_ratios[position],
            None,
        )?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:229-321
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2912-2961
/// `PSICreatePssmFromCDD` with DELTA-BLAST's options (`PSIBlastOptionsNew`: the pseudocount 0, no
/// IMPALA scaling) and its diagnostics request (`PSIDiagnosticsRequestNewEx(false)`: the
/// frequency ratios).
pub(crate) fn create_pssm_from_cdd(
    msa: &PsiCdMsa,
    pseudo_count: i32,
    sbp: &PsiScoreBlk,
) -> Result<PsiPssm, PsiError> {
    validate_cd_msa(msa)?;
    let weights = compute_frequencies_from_cds(msa);
    let mut internal = InternalPssm::new(msa.query.len());
    compute_freq_ratios_from_cds(msa, &weights, sbp, pseudo_count, &mut internal)?;
    convert_freq_ratios_to_pssm(&mut internal, &msa.query, sbp, &weights.std_prob);
    let ungapped = scale_matrix(&msa.query, &weights.std_prob, &mut internal, sbp, None)?;
    let gapped = KarlinBlk {
        k: gapped_psi_k(ungapped, sbp),
        ..sbp.gap_psi
    };
    let frequency_ratios = Some(internal.freq_ratios.clone());
    Ok(PsiPssm {
        scores: internal.into_scores(),
        gapped,
        ungapped,
        diagnostics: PsiDiagnostics {
            frequency_ratios,
            ..PsiDiagnostics::default()
        },
    })
}

#[cfg(test)]
#[path = "cd_stages_tests.rs"]
mod tests;
