//! NCBI BLAST cutoff score calculation for blastn ungapped extensions
//!
//! This module implements NCBI BLAST's cutoff score calculation for blastn ungapped extensions.
//! All functions are direct ports from NCBI C code with original comments preserved.
//!
//! Reference files:
//! - c++/src/algo/blast/core/blast_parameters.c
//! - c++/include/algo/blast/core/blast_options.h

use crate::stats::cutoffs::{cutoff_score, CutoffConversion};
use crate::stats::KarlinParams;

/// Default gap trigger bit score for nucleotide searches
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:140
pub const GAP_TRIGGER_BIT_SCORE_NUCL: f64 = 27.0;

/// Default cutoff E-value for blastn ungapped extension
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:76
pub const CUTOFF_E_BLASTN: f64 = 0.05;

/// The gap trigger's raw score from its bits and the UNGAPPED Karlin params (kbp_std), by every
/// program's one conversion (`core::x_drop::gap_trigger_raw_score`, blast_parameters.c:343-344;
/// NCBI's default for blastn is 27 bits).
pub fn gap_trigger_raw_score(gap_trigger_bits: f64, ungapped_params: &KarlinParams) -> i32 {
    crate::core::x_drop::gap_trigger_raw_score(
        gap_trigger_bits,
        ungapped_params.lambda,
        ungapped_params.k,
    )
}

/// Calculate cutoff_score_max from E-value using GAPPED Karlin params.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:943
///
/// CRITICAL: Uses GAPPED params (kbp_gap) for gapped mode.
///
/// # Arguments
/// * `evalue` - Expected value threshold
/// * `eff_searchsp` - Effective search space
/// * `gapped_params` - GAPPED Karlin-Altschul parameters
///
/// # Returns
/// cutoff_score_max (raw score)
pub fn cutoff_score_max_from_evalue(
    evalue: f64,
    eff_searchsp: i64,
    gapped_params: &KarlinParams,
) -> i32 {
    // `BLAST_Cutoffs` without gap decay from the cutoff of 1 (blast_parameters.c:916 and 943):
    // an invalid block or an e-value not above 0 (`-evalue -nan`) keeps 1.
    cutoff_score(
        evalue,
        gapped_params,
        eff_searchsp,
        None,
        CutoffConversion::BLASTN,
    )
}

/// Calculate cutoff_score for blastn ungapped extension.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:368-374
///
/// For blastn in gapped mode:
/// 1. Start with gap_trigger
/// 2. Apply scale_factor (typically 1.0 for standard scoring)
/// 3. Take MIN with cutoff_score_max
///
/// # Arguments
/// * `gap_trigger` - Gap trigger raw score
/// * `cutoff_score_max` - Maximum cutoff score from E-value
/// * `scale_factor` - Scale factor (typically 1.0 for standard scoring)
///
/// # Returns
/// Final cutoff_score for ungapped extension
pub fn cutoff_score_for_ungapped_extension(
    gap_trigger: i32,
    cutoff_score_max: i32,
    scale_factor: f64,
) -> i32 {
    // Start from `gap_trigger` in gapped mode.
    let mut new_cutoff = gap_trigger;

    // Multiply the cutoff by `sbp->scale_factor` truncated to `Int4` (NCBI truncates the factor first; this code truncates the product; every caller passes 1.0, so the cutoff is NCBI's).
    new_cutoff = (new_cutoff as f64 * scale_factor) as i32;

    // Cap the cutoff at `cutoff_score_max`.
    new_cutoff.min(cutoff_score_max)
}

/// Calculate cutoff_score for ungapped extension when gapped calculation is disabled or matrix-only scoring is active.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:348-367
///
/// This function is used for ungapped blastn or matrix_only_scoring mode.
/// Uses CUTOFF_E_BLASTN = 0.05 (fixed) and searchsp without length adjustment.
///
/// # Arguments
/// * `query_len` - Query length (for both strands if applicable, already * 2)
/// * `subject_len` - Subject length
/// * `_gap_trigger` - Gap trigger raw score (not used for blastn, but kept for compatibility with NCBI code structure)
/// * `ungapped_params` - UNGAPPED Karlin-Altschul parameters (kbp_std)
/// * `scale_factor` - Scale factor (typically 1.0 for standard scoring)
/// * `gap_decay_rate` - Gap decay rate (typically 0.0 for ungapped)
/// * `cutoff_score_max` - Maximum cutoff score from E-value (from `hit_params->cutoffs[context].cutoff_score_max`)
///
/// # Returns
/// Cutoff score for ungapped extension (raw score)
pub fn compute_blastn_cutoff_score_ungapped(
    query_len: i64,
    subject_len: i64,
    _gap_trigger: i32,
    ungapped_params: &KarlinParams,
    scale_factor: f64,
    gap_decay_rate: f64,
    cutoff_score_max: i32,
) -> i32 {
    // Use the fixed blastn ungapped cutoff E-value from `s_GetCutoffEvalue` (`CUTOFF_E_BLASTN`, 0.05).
    let cutoff_e = CUTOFF_E_BLASTN;

    // Compute the unadjusted search space as the product of the shorter sequence length and the subject length.
    let min_len = subject_len.min(query_len);
    let searchsp = min_len * subject_len;

    // `BLAST_Cutoffs` with gap decay from the cutoff of 1 (blast_parameters.c:321 and 360):
    // an invalid block keeps 1.
    let mut new_cutoff = cutoff_score(
        cutoff_e,
        ungapped_params,
        searchsp,
        Some(gap_decay_rate),
        CutoffConversion::BLASTN,
    );

    // Non-blastn programs cap the cutoff at `gap_trigger` here; blastn skips that cap.

    // Multiply the cutoff by `sbp->scale_factor` truncated to `Int4` (NCBI truncates the factor first; this code truncates the product; every caller passes 1.0, so the cutoff is NCBI's).
    new_cutoff = (new_cutoff as f64 * scale_factor) as i32;

    // Cap the cutoff at each context's `cutoff_score_max` from `hit_params->cutoffs`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:372-373
    new_cutoff = new_cutoff.min(cutoff_score_max);

    new_cutoff
}

#[cfg(test)]
#[path = "ncbi_cutoffs_tests.rs"]
mod tests;
