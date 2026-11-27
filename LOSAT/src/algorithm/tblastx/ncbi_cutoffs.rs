//! NCBI BLAST cutoff score calculation - verbatim port
//!
//! This module implements NCBI BLAST's cutoff score calculation for tblastx ungapped extensions:
//! the hit-saving cutoff and the word cutoff of each context. `BLAST_Cutoffs` itself, which they
//! call, is `stats::cutoffs` (every program's one owner).
//!
//! Reference files:
//! - c++/src/algo/blast/core/blast_parameters.c
//! - c++/src/algo/blast/core/blast_setup.c
//! - c++/src/algo/blast/core/blast_stat.c

use crate::stats::cutoffs::{cutoff_score, CutoffConversion};
use crate::stats::KarlinParams;

/// Default E-value cutoff for tblastx in BlastInitialWordParametersUpdate
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:80
pub const CUTOFF_E_TBLASTX: f64 = 1e-300;

/// Default E-value cutoff for ungapped blastp in BlastInitialWordParametersUpdate
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:77
pub const CUTOFF_E_BLASTP: f64 = 1e-300;

/// Default gap decay rate for sum statistics (`core::link_hsps` owns the constant)
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:107
pub(crate) use crate::core::link_hsps::BLAST_GAP_DECAY_RATE;

/// Calculate x_dropoff raw score from bit score using UNGAPPED Karlin params.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:219-221
///
/// CRITICAL: Uses UNGAPPED params (kbp_std), NOT gapped params (kbp_gap).
/// CRITICAL: Uses ceil() rounding, NOT trunc.
///
/// # Arguments
/// * `x_drop_bits` - X-drop in bits (NCBI default: 7.0 for protein)
/// * `ungapped_params` - UNGAPPED Karlin-Altschul parameters (kbp_std)
/// * `scale_factor` - Score scale factor (typically 1.0 for standard BLOSUM62)
///
/// # Returns
/// x_dropoff raw score for ungapped extension (ceiling rounded)
pub fn x_drop_raw_score(
    x_drop_bits: f64,
    ungapped_params: &KarlinParams,
    scale_factor: f64,
) -> i32 {
    // NCBI scales `ceil` of `x_dropoff` in bits times `NCBIMATH_LN2` over `Lambda` by `scale_factor` and truncates to `Int4`,
    // converted as the x86-64 build does past Int4 (`core::x_drop`, shared by the protein programs).
    crate::core::x_drop::ungapped_x_dropoff(x_drop_bits, ungapped_params.lambda, scale_factor)
}

/// Calculate gap_trigger raw score from bit score using UNGAPPED Karlin params.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:340-345
///
/// CRITICAL: Uses UNGAPPED params (kbp_std), NOT gapped params (kbp_gap).
/// CRITICAL: Uses trunc/floor via C's (Int4) cast, NOT ceil.
///
/// # Arguments
/// * `bit_trigger` - Gap trigger in bits (NCBI default: 22.0 for protein)
/// * `ungapped_params` - UNGAPPED Karlin-Altschul parameters (kbp_std)
///
/// # Returns
/// Raw score threshold for gap trigger (truncated to integer; a value past `Int4`, such as one of
/// psiblast's `-gap_trigger 1e300`, `-1e300` or `-inf`, is converted as the x86-64 build of NCBI
/// converts it, to `INT_MIN`: `core::blast_util::ncbi_int4_from_double`, as the X-drops are in
/// `core::x_drop`)
pub fn gap_trigger_raw_score(bit_trigger: f64, ungapped_params: &KarlinParams) -> i32 {
    // NCBI computes `gap_trigger` in bits times `NCBIMATH_LN2`, plus `logK`, all divided by `Lambda`, truncated to `Int4`;
    // every program's one conversion (`core::x_drop`).
    crate::core::x_drop::gap_trigger_raw_score(
        bit_trigger,
        ungapped_params.lambda,
        ungapped_params.k,
    )
}

/// Result of effective length calculation for TBLASTX -subject mode.
///
/// This struct mirrors NCBI's `query_info->contexts[ctx].length_adjustment` and
/// `query_info->contexts[ctx].eff_searchsp` fields which are computed once in
/// `BLAST_CalcEffLengths` and then referenced by cutoff and sum-stats logic.
#[derive(Debug, Clone, Copy)]
pub struct EffLengthsResult {
    /// Length adjustment (Int4 in NCBI: `query_info->contexts[ctx].length_adjustment`)
    pub length_adjustment: i64,
    /// Effective search space (Int8 in NCBI: `query_info->contexts[ctx].eff_searchsp`)
    pub eff_searchsp: i64,
}

/// Calculate effective lengths for the complete local-subject database.
///
/// This is the NCBI-parity function that computes BOTH length_adjustment and eff_searchsp
/// in a single call, matching `BLAST_CalcEffLengths` from blast_setup.c:700-850.
/// These values are then stored and referenced by cutoff and sum-stats calculations.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:734-847
///
/// # Arguments
/// * `query_len_aa` - Query length in amino acids (for this context/frame)
/// * `subject_len_nucl` - Subject length in nucleotides (will be divided by 3)
/// * `karlin_params` - Karlin-Altschul parameters
///
/// # Returns
/// `EffLengthsResult` containing both length_adjustment and eff_searchsp, or `None` when the
/// effective search space overflows NCBI's `Int8` (a `-dbsize` near the `Int8` limit), where
/// NCBI's arithmetic overflows.
#[cfg(test)]
pub fn try_compute_eff_lengths_tblastx(
    query_len_aa: i64,
    subject_len_nucl: i64,
    db_num_seqs: i64,
    karlin_params: &KarlinParams,
) -> Option<EffLengthsResult> {
    if query_len_aa <= 0 {
        return Some(EffLengthsResult {
            length_adjustment: 0,
            eff_searchsp: 0,
        });
    }
    // `stats::setup::calc_eff_lengths` (one owner for the six programs, line RG D1a-3), with the
    // translated subject's length in residues and NCBI's `Int4` count.
    let lengths = crate::stats::setup::calc_eff_lengths(
        query_len_aa,
        subject_len_nucl / 3,
        db_num_seqs as i32,
        karlin_params,
        0,
        "TBLASTX",
    )
    .ok()?;
    Some(EffLengthsResult {
        length_adjustment: i64::from(lengths.length_adjustment),
        eff_searchsp: lengths.eff_searchsp,
    })
}

/// Calculate cutoff_score for BlastInitialWordParametersUpdate (per-subject update).
///
/// This is the NCBI tblastx ungapped path in BlastInitialWordParametersUpdate.
/// It uses a FIXED E-value cutoff (1e-300) and a simple searchsp formula
/// WITHOUT length adjustment.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:348-374
///
/// CRITICAL: subj_length is in NUCLEOTIDES (NOT divided by 3!)
/// CRITICAL: Uses CUTOFF_E_TBLASTX = 1e-300 (NOT user's E-value threshold)
/// CRITICAL: `searchsp` is the minimum of `query_len_aa` and `subj_len_nucl`, times `subj_len_nucl` (mixes AA and nucleotide!)
/// CRITICAL: `dodecay` enabled with `gap_decay_rate`
///
/// # Arguments
/// * `query_len_aa` - Query length in amino acids (translated frame)
/// * `subject_len_nucl` - Subject length in NUCLEOTIDES (NOT divided by 3!)
/// * `gap_trigger` - Gap trigger raw score (from gap_trigger_raw_score)
/// * `cutoff_score_max` - Maximum cutoff from hit_params (from cutoff_score_max_for_tblastx)
/// * `gap_decay_rate` - Gap decay rate (typically BLAST_GAP_DECAY_RATE = 0.5)
/// * `ungapped_params` - UNGAPPED Karlin-Altschul parameters (kbp_std)
/// * `scale_factor` - Score scale factor (typically 1.0 for standard BLOSUM62)
///
/// # Returns
/// Per-subject cutoff_score for ungapped extension (raw score)
pub fn cutoff_score_for_update_tblastx(
    query_len_aa: i64,
    subject_len_nucl: i64,
    gap_trigger: i32,
    cutoff_score_max: i32,
    gap_decay_rate: f64,
    ungapped_params: &KarlinParams,
    scale_factor: f64,
) -> i32 {
    word_cutoff_score_for_update(
        CUTOFF_E_TBLASTX,
        query_len_aa,
        subject_len_nucl,
        gap_trigger,
        cutoff_score_max,
        gap_decay_rate,
        ungapped_params,
        scale_factor,
    )
}

/// The word cutoff of `BlastInitialWordParametersUpdate` for a program without gapped extension of
/// its initial words (TBLASTX; BLASTP with `-ungapped`), for the program's cutoff
/// E-value `cutoff_e` (`s_GetCutoffEvalue`: `CUTOFF_E_TBLASTX` or `CUTOFF_E_BLASTP`, both 1e-300).
/// The body of `cutoff_score_for_update_tblastx`, whose arguments are named for TBLASTX
/// (`subject_len_nucl` is the subject length that NCBI's `subj_length` holds, whatever the
/// program's alphabet); the same NCBI lines (blast_parameters.c:346-374).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:346-374
/// BLASTX's word cutoff (`blastx::parameters::word_cutoff`, the E-value 1.0 of `CUTOFF_E_BLASTX`)
/// and TBLASTN's (`tblastn::stage_d_stats::ungapped_word_cutoff`) are the same lines with their
/// programs' lengths; all three call `stats::cutoffs::cutoff_score` for `BLAST_Cutoffs`.
#[allow(clippy::too_many_arguments)]
pub fn word_cutoff_score_for_update(
    cutoff_e: f64,
    query_len_aa: i64,
    subject_len_nucl: i64,
    gap_trigger: i32,
    cutoff_score_max: i32,
    gap_decay_rate: f64,
    ungapped_params: &KarlinParams,
    scale_factor: f64,
) -> i32 {
    // Search space is the smaller of `subj_length` and `query_length` multiplied by `subj_length`.
    // NOTE: This intentionally mixes AA (query_len_aa) and nucleotide (subject_len_nucl) lengths!
    // This is exactly what NCBI does - see blast_parameters.c:360-362
    let min_len = (query_len_aa as u64).min(subject_len_nucl as u64);
    // NCBI multiplies two `Uint8` values (unsigned 64-bit wrap when `subj_length` is negative).
    let searchsp = min_len.wrapping_mul(subject_len_nucl as u64) as i64;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:348-374 (BlastInitialWordParametersUpdate)
    if crate::utils::xenv::debug_cutoffs_is_some() && subject_len_nucl > 600_000 {
        eprintln!(
            "[DEBUG CUTOFF_UPDATE] query_len_aa={}, subject_len_nucl={}",
            query_len_aa, subject_len_nucl
        );
        eprintln!(
            "[DEBUG CUTOFF_UPDATE] searchsp={} (MIN({}, {}) * {})",
            searchsp, query_len_aa, subject_len_nucl, subject_len_nucl
        );
    }

    // Convert the program cutoff E-value from `s_GetCutoffEvalue` (`CUTOFF_E_TBLASTX`, 1e-300)
    // using `BLAST_Cutoffs` with gap decay enabled (`gap_decay_rate`).
    // Gap decay enabled (`dodecay`) means we apply the gap decay divisor.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:288,310-311
    // Without sum statistics (`-sum_stats false`, no `link_hsp_params`) NCBI passes 0, and the
    // cutoff is the same: 1e-300 and 1e-300 * (1 - 0.5) are both below the 1e-297 that
    // `BlastKarlinEtoS_simple` takes at least (blast_stat.c:4049,4059, `stats::cutoffs::K_SMALL_FLOAT`),
    // so the callers pass `BLAST_GAP_DECAY_RATE` in both modes. The cutoff starts at 1
    // (blast_parameters.c:321): `stats::cutoffs::cutoff_score`, with exception 8's conversion.
    let mut new_cutoff = cutoff_score(
        cutoff_e,
        ungapped_params,
        searchsp,
        Some(gap_decay_rate),
        CutoffConversion::FiniteQuotient,
    );

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:365-373
    if crate::utils::xenv::debug_cutoffs_is_some() && subject_len_nucl > 600_000 {
        eprintln!("[DEBUG CUTOFF_UPDATE] update_cutoff={} (from CUTOFF_E_TBLASTX=1e-300, searchsp={}, gap_decay={})", 
                  new_cutoff, searchsp, gap_decay_rate);
        eprintln!(
            "[DEBUG CUTOFF_UPDATE] gap_trigger={}, cutoff_score_max={}",
            gap_trigger, cutoff_score_max
        );
    }

    // Cap the cutoff score at `gap_trigger`.
    // (Blastn exception at line 366 does not apply to tblastx)
    new_cutoff = new_cutoff.min(gap_trigger);

    // Multiply the cutoff score by `sbp->scale_factor` truncated to `Int4` (NCBI truncates the factor first; this code truncates the product; every caller passes 1.0, so the cutoff is NCBI's).
    new_cutoff = (new_cutoff as f64 * scale_factor) as i32;

    // Cap the cutoff score at the context's `cutoff_score_max` from `hit_params`.
    let final_cutoff = new_cutoff.min(cutoff_score_max);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:371-374
    if crate::utils::xenv::debug_cutoffs_is_some() && subject_len_nucl > 600_000 {
        eprintln!("[DEBUG CUTOFF_UPDATE] final_cutoff={} (MIN(update_cutoff={}, gap_trigger={}, cutoff_score_max={}))", 
                  final_cutoff, new_cutoff, gap_trigger, cutoff_score_max);
    }

    final_cutoff
}

/// Calculate cutoff_score_max for BlastHitSavingParametersNew.
///
/// This is the initial cutoff_score_max set during parameter setup, using
/// the user's E-value threshold and the effective search space (WITH length adjustment).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:942-946
///
/// CRITICAL: Uses user's E-value threshold (NOT CUTOFF_E_TBLASTX)
/// CRITICAL: Uses eff_searchsp WITH length adjustment (from BLAST_CalcEffLengths)
/// CRITICAL: gap decay disabled (`dodecay` false, `gap_decay_rate` 0; no decay)
/// CRITICAL: For tblastx, uses UNGAPPED params (kbp) since kbp_gap is NULL
///
/// # Arguments
/// * `eff_searchsp` - Effective search space (from `try_compute_eff_lengths_tblastx`)
/// * `evalue_threshold` - User's E-value threshold (typically 10.0)
/// * `ungapped_params` - UNGAPPED Karlin-Altschul parameters (kbp for tblastx)
///
/// # Returns
/// cutoff_score_max (raw score)
pub fn cutoff_score_max_for_tblastx(
    eff_searchsp: i64,
    evalue_threshold: f64,
    ungapped_params: &KarlinParams,
) -> i32 {
    // `BLAST_Cutoffs` without gap decay from the cutoff of 1 (blast_parameters.c:916 and 943):
    // a NaN e-value (`-evalue -nan`) is not greater than 0, so the cutoff stays 1.
    cutoff_score(
        evalue_threshold,
        ungapped_params,
        eff_searchsp,
        None,
        CutoffConversion::FiniteQuotient,
    )
}

#[cfg(test)]
#[path = "ncbi_cutoffs_tests.rs"]
mod tests;
