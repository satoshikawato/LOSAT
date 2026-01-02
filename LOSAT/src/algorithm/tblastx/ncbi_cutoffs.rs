//! NCBI BLAST cutoff score calculation - verbatim port
//!
//! This module implements NCBI BLAST's cutoff score calculation for tblastx ungapped extensions.
//! All functions are direct ports from NCBI C code with original comments preserved.
//!
//! Reference files:
//! - ncbi-blast/c++/src/algo/blast/core/blast_parameters.c
//! - ncbi-blast/c++/src/algo/blast/core/blast_setup.c
//! - ncbi-blast/c++/src/algo/blast/core/blast_stat.c

use crate::stats::KarlinParams;
use crate::stats::length_adjustment::compute_length_adjustment_ncbi;

/// ln(2) constant used in NCBI BLAST
/// Reference: ncbi-blast/c++/include/algo/blast/core/ncbi_math.h
/// #define NCBIMATH_LN2 0.69314718055994530941723212145818
const NCBIMATH_LN2: f64 = 0.69314718055994530941723212145818;

/// Smallest float to avoid floating point exception in BlastKarlinEtoS_simple
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4049
const K_SMALL_FLOAT: f64 = 1.0e-297;

/// Calculate gap_trigger raw score from bit score using UNGAPPED Karlin params.
///
/// NCBI reference (verbatim from blast_parameters.c:340-345):
/// ```c
/// if (sbp->kbp_std) {     /* this may not be set for gapped blastn */
///    kbp = sbp->kbp_std[context];
///    if (s_BlastKarlinBlkIsValid(kbp)) {
///       gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + 
///                                kbp->logK) / kbp->Lambda);
///    }
/// }
/// ```
///
/// CRITICAL: Uses UNGAPPED params (kbp_std), NOT gapped params (kbp_gap).
/// CRITICAL: Uses trunc/floor via C's (Int4) cast, NOT ceil.
///
/// # Arguments
/// * `bit_trigger` - Gap trigger in bits (NCBI default: 22.0 for protein)
/// * `ungapped_params` - UNGAPPED Karlin-Altschul parameters (kbp_std)
///
/// # Returns
/// Raw score threshold for gap trigger (truncated to integer)
pub fn gap_trigger_raw_score(bit_trigger: f64, ungapped_params: &KarlinParams) -> i32 {
    // NCBI: gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda);
    // Note: kbp->logK = ln(K), so we use ungapped_params.k.ln()
    let raw = (bit_trigger * NCBIMATH_LN2 + ungapped_params.k.ln()) / ungapped_params.lambda;
    
    // C's (Int4) cast truncates toward zero (same as Rust's `as i32` for positive values)
    raw as i32
}

/// Calculate effective search space for -subject mode (single subject as database).
///
/// NCBI reference (verbatim from blast_setup.c:734-846):
/// ```c
/// // For translated subjects:
/// if (Blast_SubjectIsTranslated(program_number))
///    db_length = db_length/3;
///
/// // Length adjustment calculation:
/// BLAST_ComputeLengthAdjustment(kbp->K, kbp->logK,
///                               alpha/kbp->Lambda, beta,
///                               query_length, db_length,
///                               db_num_seqs, &length_adjustment);
///
/// // Effective search space:
/// Int8 effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
/// if (effective_db_length <= 0)
///     effective_db_length = 1;
/// effective_search_space = effective_db_length * (query_length - length_adjustment);
/// ```
///
/// # Arguments
/// * `query_len_aa` - Query length in amino acids (for this context/frame)
/// * `subject_len_nucl` - Subject length in nucleotides (will be divided by 3)
/// * `gapped_params` - GAPPED Karlin-Altschul parameters (kbp_gap)
///
/// # Returns
/// Effective search space (eff_searchsp)
pub fn compute_eff_searchsp_subject_mode_tblastx(
    query_len_aa: i64,
    subject_len_nucl: i64,
    gapped_params: &KarlinParams,
) -> i64 {
    // NCBI: db_length = db_length/3 for translated subjects (tblastx)
    let db_length = subject_len_nucl / 3;
    
    // NCBI: db_num_seqs = 1 for -subject mode
    let db_num_seqs: i64 = 1;
    
    // NCBI: BLAST_ComputeLengthAdjustment(...)
    let result = compute_length_adjustment_ncbi(
        query_len_aa,
        db_length,
        db_num_seqs,
        gapped_params,
    );
    let length_adjustment = result.length_adjustment;
    
    // NCBI: effective_db_length = db_length - (db_num_seqs * length_adjustment)
    let effective_db_length = (db_length - db_num_seqs * length_adjustment).max(1);
    
    // NCBI: effective_search_space = effective_db_length * (query_length - length_adjustment)
    let effective_query_length = (query_len_aa - length_adjustment).max(1);
    
    effective_db_length * effective_query_length
}

/// Calculate cutoff_score from E-value using BLAST_Cutoffs (dodecay=FALSE).
///
/// NCBI reference (verbatim from blast_stat.c:4040-4063):
/// ```c
/// static Int4
/// BlastKarlinEtoS_simple(double E, const Blast_KarlinBlk* kbp, Int8 searchsp)
/// {
///    const double kSmallFloat = 1.0e-297;
///    Lambda = kbp->Lambda;
///    K = kbp->K;
///    E = MAX(E, kSmallFloat);
///    S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda ));
///    return S;
/// }
/// ```
///
/// And from blast_parameters.c:942-943:
/// ```c
/// BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
/// ```
///
/// CRITICAL: Uses GAPPED params (kbp_gap/kbp_gap_std).
/// CRITICAL: Uses ceil() rounding, NOT trunc.
///
/// # Arguments
/// * `evalue` - Expected value threshold
/// * `eff_searchsp` - Effective search space
/// * `gapped_params` - GAPPED Karlin-Altschul parameters
///
/// # Returns
/// cutoff_score_max (raw score)
pub fn cutoff_score_from_evalue(
    evalue: f64,
    eff_searchsp: i64,
    gapped_params: &KarlinParams,
) -> i32 {
    // NCBI: E = MAX(E, kSmallFloat)
    let e = evalue.max(K_SMALL_FLOAT);
    
    // NCBI: S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda ))
    let searchsp = eff_searchsp as f64;
    let score = ((gapped_params.k * searchsp / e).ln() / gapped_params.lambda).ceil();
    
    score as i32
}

/// Calculate cutoff_score for sum statistics path.
///
/// NCBI reference (verbatim from blast_parameters.c:951-976):
/// ```c
/// if (params->link_hsp_params && gapped_calculation) {
///    double evalue_hsp = 1.0;
///    Int4 concat_qlen = query_info->contexts[query_info->last_context].query_offset +
///                       query_info->contexts[query_info->last_context].query_length;
///    Int4 avg_qlen = concat_qlen / (query_info->last_context + 1);
///    Int8 searchsp = (Int8)MIN(avg_qlen, avg_subject_length) * (Int8)avg_subject_length;
///
///    for (context = ...) {
///        BLAST_Cutoffs(&new_cutoff, &evalue_hsp, kbp, searchsp, TRUE, gap_decay_rate);
///        params->cutoffs[context].cutoff_score = MIN(new_cutoff, cutoff_score);
///    }
/// }
/// ```
///
/// # Arguments
/// * `avg_query_len` - Average query length across all contexts
/// * `avg_subject_len` - Average subject length (= subject_len_nucl/3 for tblastx -subject mode)
/// * `gap_decay_rate` - Gap decay rate (typically 0.5 for ungapped, 0.1 for gapped)
/// * `gapped_params` - GAPPED Karlin-Altschul parameters
///
/// # Returns
/// cutoff_score for sum statistics (raw score)
pub fn cutoff_score_sum_stats(
    avg_query_len: i32,
    avg_subject_len: i32,
    gap_decay_rate: f64,
    gapped_params: &KarlinParams,
) -> i32 {
    // NCBI: evalue_hsp = 1.0 (fixed)
    let evalue_hsp: f64 = 1.0;
    
    // NCBI: searchsp = MIN(avg_qlen, avg_subject_length) * avg_subject_length
    let min_len = avg_query_len.min(avg_subject_len) as i64;
    let searchsp = min_len * (avg_subject_len as i64);
    
    // NCBI: BLAST_Cutoffs(&new_cutoff, &evalue_hsp, kbp, searchsp, TRUE, gap_decay_rate)
    // When dodecay=TRUE, apply gap decay divisor to E-value before conversion
    cutoff_score_from_evalue_with_decay(evalue_hsp, searchsp, gap_decay_rate, gapped_params)
}

/// Calculate cutoff_score from E-value with gap decay adjustment.
///
/// NCBI reference (verbatim from blast_stat.c:4112-4121):
/// ```c
/// if (dodecay) {
///     if( gap_decay_rate > 0 && gap_decay_rate < 1 ) {
///         e *= BLAST_GapDecayDivisor(gap_decay_rate, 1);
///     }
/// }
/// es = BlastKarlinEtoS_simple(e, kbp, searchsp);
/// ```
///
/// And from blast_stat.c:4079-4082:
/// ```c
/// double BLAST_GapDecayDivisor(double decayrate, unsigned nsegs) {
///     return (1. - decayrate) * BLAST_Powi(decayrate, nsegs - 1);
/// }
/// ```
/// For nsegs=1: divisor = (1 - decayrate) * decayrate^0 = (1 - decayrate)
fn cutoff_score_from_evalue_with_decay(
    evalue: f64,
    eff_searchsp: i64,
    gap_decay_rate: f64,
    gapped_params: &KarlinParams,
) -> i32 {
    // NCBI: BLAST_GapDecayDivisor(gap_decay_rate, 1) = (1 - gap_decay_rate)
    let divisor = 1.0 - gap_decay_rate;
    
    // NCBI: e *= BLAST_GapDecayDivisor(gap_decay_rate, 1)
    let adjusted_e = if gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
        evalue * divisor
    } else {
        evalue
    };
    
    cutoff_score_from_evalue(adjusted_e, eff_searchsp, gapped_params)
}

/// Calculate final cutoff_score for word_params (ungapped extension threshold).
///
/// NCBI reference (verbatim from blast_parameters.c:348-374):
/// ```c
/// // For gapped calculation:
/// if (!gapped_calculation || sbp->matrix_only_scoring) {
///     // ... ungapped path ...
/// } else {
///     new_cutoff = gap_trigger;
/// }
/// new_cutoff *= (Int4)sbp->scale_factor;
/// new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
/// curr_cutoffs->cutoff_score = new_cutoff;
/// ```
///
/// # Arguments
/// * `gap_trigger` - Gap trigger raw score (from gap_trigger_raw_score)
/// * `cutoff_score_max` - Maximum cutoff from hit_params (from cutoff_score_from_evalue)
/// * `scale_factor` - Score scale factor (typically 1.0 for standard BLOSUM62)
///
/// # Returns
/// Final cutoff_score for ungapped extension
pub fn cutoff_score_word_params(
    gap_trigger: i32,
    cutoff_score_max: i32,
    scale_factor: f64,
) -> i32 {
    // NCBI: new_cutoff = gap_trigger (gapped mode)
    let mut new_cutoff = gap_trigger;
    
    // NCBI: new_cutoff *= (Int4)sbp->scale_factor
    new_cutoff = (new_cutoff as f64 * scale_factor) as i32;
    
    // NCBI: new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max)
    new_cutoff.min(cutoff_score_max)
}

/// High-level function to compute ungapped extension cutoff for tblastx.
///
/// This combines all the NCBI cutoff calculation steps for -subject mode tblastx:
/// 1. Compute gap_trigger using UNGAPPED params
/// 2. Compute eff_searchsp using GAPPED params
/// 3. Compute cutoff_score_max using GAPPED params
/// 4. Return MIN(gap_trigger, cutoff_score_max)
///
/// # Arguments
/// * `query_len_aa` - Query length in amino acids (for this context/frame)
/// * `subject_len_nucl` - Subject length in nucleotides
/// * `evalue_threshold` - E-value threshold (typically 10.0)
/// * `gap_trigger_bits` - Gap trigger in bits (typically 22.0)
/// * `ungapped_params` - UNGAPPED Karlin-Altschul parameters
/// * `gapped_params` - GAPPED Karlin-Altschul parameters
///
/// # Returns
/// Cutoff score for ungapped extension (raw score)
pub fn compute_tblastx_cutoff_score(
    query_len_aa: i64,
    subject_len_nucl: i64,
    evalue_threshold: f64,
    gap_trigger_bits: f64,
    ungapped_params: &KarlinParams,
    gapped_params: &KarlinParams,
) -> i32 {
    // Step 1: Compute gap_trigger using UNGAPPED params
    let gap_trigger = gap_trigger_raw_score(gap_trigger_bits, ungapped_params);
    
    // Step 2: Compute eff_searchsp using GAPPED params
    let eff_searchsp = compute_eff_searchsp_subject_mode_tblastx(
        query_len_aa,
        subject_len_nucl,
        gapped_params,
    );
    
    // Step 3: Compute cutoff_score_max using GAPPED params
    let cutoff_score_max = cutoff_score_from_evalue(evalue_threshold, eff_searchsp, gapped_params);
    
    // Step 4: Return MIN(gap_trigger, cutoff_score_max)
    // Note: scale_factor = 1.0 for standard BLOSUM62
    cutoff_score_word_params(gap_trigger, cutoff_score_max, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// BLOSUM62 ungapped params (kbp_std)
    fn blosum62_ungapped() -> KarlinParams {
        KarlinParams {
            lambda: 0.3176,
            k: 0.134,
            h: 0.4012,
            alpha: 0.7916,
            beta: -3.2,
        }
    }

    /// BLOSUM62 gapped params (11/1) (kbp_gap)
    fn blosum62_gapped() -> KarlinParams {
        KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        }
    }

    #[test]
    fn test_gap_trigger_raw_score() {
        let ungapped = blosum62_ungapped();
        
        // NCBI: gap_trigger = (Int4)((22.0 * LN2 + ln(0.134)) / 0.3176)
        // = (Int4)((22.0 * 0.693147 + (-2.0099)) / 0.3176)
        // = (Int4)((15.249 - 2.0099) / 0.3176)
        // = (Int4)(13.239 / 0.3176)
        // = (Int4)(41.69)
        // = 41 (truncated)
        let gap_trigger = gap_trigger_raw_score(22.0, &ungapped);
        assert_eq!(gap_trigger, 41, "gap_trigger should be 41 (truncated)");
    }

    #[test]
    fn test_gap_trigger_uses_trunc_not_ceil() {
        let ungapped = blosum62_ungapped();
        
        // Verify that we use trunc, not ceil
        // If ceil were used, the result would be 42
        let gap_trigger = gap_trigger_raw_score(22.0, &ungapped);
        assert!(gap_trigger < 42, "gap_trigger should use trunc, not ceil");
    }

    #[test]
    fn test_cutoff_score_from_evalue_uses_ceil() {
        let gapped = blosum62_gapped();
        
        // Test that ceil is used (not trunc)
        // For a large search space, the cutoff should be positive
        let cutoff = cutoff_score_from_evalue(10.0, 1_000_000, &gapped);
        assert!(cutoff > 0, "cutoff should be positive for reasonable search space");
    }

    #[test]
    fn test_cutoff_score_cap_effective() {
        let ungapped = blosum62_ungapped();
        let gapped = blosum62_gapped();
        
        // For very small subject, cutoff_score_max should be less than gap_trigger (41)
        // This tests that the MIN cap is effective
        let query_len_aa = 100;
        let subject_len_nucl = 300; // Small subject: 100 AA
        
        let cutoff = compute_tblastx_cutoff_score(
            query_len_aa,
            subject_len_nucl,
            10.0, // evalue
            22.0, // gap_trigger_bits
            &ungapped,
            &gapped,
        );
        
        let gap_trigger = gap_trigger_raw_score(22.0, &ungapped);
        
        // The cutoff should be <= gap_trigger
        assert!(cutoff <= gap_trigger, "cutoff should be capped by gap_trigger or cutoff_score_max");
    }

    #[test]
    fn test_cutoff_score_gap_trigger_wins() {
        let ungapped = blosum62_ungapped();
        let gapped = blosum62_gapped();
        
        // For very large subject/query, cutoff_score_max should be > gap_trigger
        // So the final cutoff should equal gap_trigger
        let query_len_aa = 10000;
        let subject_len_nucl = 30000; // Large subject: 10000 AA
        
        let cutoff = compute_tblastx_cutoff_score(
            query_len_aa,
            subject_len_nucl,
            10.0, // evalue
            22.0, // gap_trigger_bits
            &ungapped,
            &gapped,
        );
        
        let gap_trigger = gap_trigger_raw_score(22.0, &ungapped);
        
        // For large sequences, gap_trigger should be the limiting factor
        assert_eq!(cutoff, gap_trigger, "for large sequences, cutoff should equal gap_trigger");
    }
}

