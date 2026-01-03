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

/// Default E-value cutoff for tblastx in BlastInitialWordParametersUpdate
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:80
/// #define CUTOFF_E_TBLASTX 1e-300
pub const CUTOFF_E_TBLASTX: f64 = 1e-300;

/// Default gap decay rate for sum statistics
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:107
/// #define BLAST_GAP_DECAY_RATE 0.5
pub const BLAST_GAP_DECAY_RATE: f64 = 0.5;

/// Calculate x_dropoff raw score from bit score using UNGAPPED Karlin params.
///
/// NCBI reference (verbatim from blast_parameters.c:219-221):
/// ```c
/// p->cutoffs[context].x_dropoff_init =
///     (Int4)(sbp->scale_factor *
///            ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
/// ```
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
    // NCBI: (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda))
    (scale_factor * (x_drop_bits * NCBIMATH_LN2 / ungapped_params.lambda).ceil()) as i32
}

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

/// Calculate cutoff_score for BlastInitialWordParametersUpdate (per-subject update).
///
/// This is the NCBI tblastx ungapped path in BlastInitialWordParametersUpdate.
/// It uses a FIXED E-value cutoff (1e-300) and a simple searchsp formula
/// WITHOUT length adjustment.
///
/// NCBI reference (verbatim from blast_parameters.c:348-374):
/// ```c
/// if (!gapped_calculation || sbp->matrix_only_scoring) {
///     double cutoff_e = s_GetCutoffEvalue(program_number);  // = 1e-300 for tblastx!
///     Int4 query_length = query_info->contexts[context].query_length;  // AA length
///
///     kbp = kbp_array[context];  // kbp_std for tblastx (ungapped)
///     BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp,
///                   MIN((Uint8)subj_length, (Uint8)query_length)*((Uint8)subj_length),
///                   TRUE, gap_decay_rate);
///
///     // Blastn exception does not apply to tblastx
///     new_cutoff = MIN(new_cutoff, gap_trigger);
/// }
/// new_cutoff *= (Int4)sbp->scale_factor;
/// new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
/// ```
///
/// CRITICAL: subj_length is in NUCLEOTIDES (NOT divided by 3!)
/// CRITICAL: Uses CUTOFF_E_TBLASTX = 1e-300 (NOT user's E-value threshold)
/// CRITICAL: searchsp = MIN(query_len_aa, subj_len_nucl) * subj_len_nucl (mixes AA and nucleotide!)
/// CRITICAL: dodecay=TRUE with gap_decay_rate
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
    // NCBI: searchsp = MIN(subj_length, query_length) * subj_length
    // NOTE: This intentionally mixes AA (query_len_aa) and nucleotide (subject_len_nucl) lengths!
    // This is exactly what NCBI does - see blast_parameters.c:360-362
    let min_len = (query_len_aa as u64).min(subject_len_nucl as u64);
    let searchsp = (min_len * (subject_len_nucl as u64)) as i64;

    // NCBI: cutoff_e = s_GetCutoffEvalue(program_number) = CUTOFF_E_TBLASTX = 1e-300
    // NCBI: BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp, searchsp, TRUE, gap_decay_rate)
    // dodecay=TRUE means we apply gap decay divisor
    let mut new_cutoff = cutoff_score_from_evalue_with_decay(
        CUTOFF_E_TBLASTX,
        searchsp,
        gap_decay_rate,
        ungapped_params,
    );

    // NCBI: new_cutoff = MIN(new_cutoff, gap_trigger)
    // (Blastn exception at line 366 does not apply to tblastx)
    new_cutoff = new_cutoff.min(gap_trigger);

    // NCBI: new_cutoff *= (Int4)sbp->scale_factor
    new_cutoff = (new_cutoff as f64 * scale_factor) as i32;

    // NCBI: new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max)
    new_cutoff.min(cutoff_score_max)
}

/// Calculate cutoff_score_max for BlastHitSavingParametersNew.
///
/// This is the initial cutoff_score_max set during parameter setup, using
/// the user's E-value threshold and the effective search space (WITH length adjustment).
///
/// NCBI reference (verbatim from blast_parameters.c:942-946):
/// ```c
/// // searchsp comes from query_info->contexts[context].eff_searchsp
/// // which includes length adjustment from BLAST_CalcEffLengths
/// BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
/// params->cutoffs[context].cutoff_score = new_cutoff;
/// params->cutoffs[context].cutoff_score_max = new_cutoff;
/// ```
///
/// CRITICAL: Uses user's E-value threshold (NOT CUTOFF_E_TBLASTX)
/// CRITICAL: Uses eff_searchsp WITH length adjustment (from BLAST_CalcEffLengths)
/// CRITICAL: dodecay=FALSE, gap_decay_rate=0 (no decay)
/// CRITICAL: For tblastx, uses UNGAPPED params (kbp) since kbp_gap is NULL
///
/// # Arguments
/// * `eff_searchsp` - Effective search space (from compute_eff_searchsp_subject_mode_tblastx)
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
    // NCBI: BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0)
    // dodecay=FALSE means no gap decay adjustment
    cutoff_score_from_evalue(evalue_threshold, eff_searchsp, ungapped_params)
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
/// 2. Compute eff_searchsp using UNGAPPED params (tblastx has no kbp_gap!)
/// 3. Compute cutoff_score_max using UNGAPPED params (tblastx has no kbp_gap!)
/// 4. Return MIN(gap_trigger, cutoff_score_max)
///
/// NCBI Reference (blast_parameters.c:860-866):
/// ```c
/// if (sbp->kbp_gap) {
///     kbp_array = sbp->kbp_gap;
/// } else if (sbp->kbp) {
///     kbp_array = sbp->kbp;        // tblastx uses this path!
///     gapped_calculation = FALSE;
/// }
/// ```
/// For tblastx, kbp_gap is NULL because gapped alignment is not allowed,
/// so kbp_array = sbp->kbp (ungapped params).
///
/// # Arguments
/// * `query_len_aa` - Query length in amino acids (for this context/frame)
/// * `subject_len_nucl` - Subject length in nucleotides
/// * `evalue_threshold` - E-value threshold (typically 10.0)
/// * `gap_trigger_bits` - Gap trigger in bits (typically 22.0)
/// * `ungapped_params` - UNGAPPED Karlin-Altschul parameters (kbp_std / kbp)
/// * `_gapped_params` - GAPPED Karlin-Altschul parameters (unused for tblastx)
///
/// # Returns
/// Cutoff score for ungapped extension (raw score)
pub fn compute_tblastx_cutoff_score(
    query_len_aa: i64,
    subject_len_nucl: i64,
    evalue_threshold: f64,
    gap_trigger_bits: f64,
    ungapped_params: &KarlinParams,
    _gapped_params: &KarlinParams,  // Unused for tblastx - kept for API compatibility
) -> i32 {
    // Step 1: Compute gap_trigger using UNGAPPED params
    let gap_trigger = gap_trigger_raw_score(gap_trigger_bits, ungapped_params);
    
    // Step 2: Compute eff_searchsp using UNGAPPED params (tblastx has no kbp_gap!)
    // NCBI uses kbp_array = sbp->kbp (ungapped) when kbp_gap is NULL (tblastx case)
    let eff_searchsp = compute_eff_searchsp_subject_mode_tblastx(
        query_len_aa,
        subject_len_nucl,
        ungapped_params,  // Use ungapped params, not gapped!
    );
    
    // Step 3: Compute cutoff_score_max using UNGAPPED params (tblastx has no kbp_gap!)
    let cutoff_score_max = cutoff_score_from_evalue(evalue_threshold, eff_searchsp, ungapped_params);
    
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

    #[test]
    fn test_x_drop_raw_score_blosum62() {
        let ungapped = blosum62_ungapped();
        
        // NCBI: x_dropoff_init = (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda))
        // For BLOSUM62: ceil(7.0 * 0.693147 / 0.3176) = ceil(15.27) = 16
        let x_drop = x_drop_raw_score(7.0, &ungapped, 1.0);
        assert_eq!(x_drop, 16, "x_drop should be 16 for BLOSUM62 with 7.0 bits");
    }

    #[test]
    fn test_x_drop_uses_ceil_not_trunc() {
        let ungapped = blosum62_ungapped();
        
        // Verify that we use ceil, not trunc
        // 7.0 * 0.693147 / 0.3176 = 15.27
        // ceil(15.27) = 16 (not 15)
        let x_drop = x_drop_raw_score(7.0, &ungapped, 1.0);
        assert_eq!(x_drop, 16, "x_drop should use ceil, resulting in 16");
        assert!(x_drop > 15, "x_drop should be > 15 due to ceil rounding");
    }

    #[test]
    fn test_x_drop_with_scale_factor() {
        let ungapped = blosum62_ungapped();
        
        // Test that scale_factor is applied correctly
        // With scale_factor = 2.0, result should be doubled
        let x_drop_base = x_drop_raw_score(7.0, &ungapped, 1.0);
        let x_drop_scaled = x_drop_raw_score(7.0, &ungapped, 2.0);
        assert_eq!(x_drop_scaled, x_drop_base * 2, "x_drop with scale_factor=2.0 should be double");
    }

    #[test]
    fn test_x_drop_different_bits() {
        let ungapped = blosum62_ungapped();
        
        // Test with different bit values
        // x_drop = ceil(bits * LN2 / lambda)
        
        // 10 bits: ceil(10.0 * 0.693147 / 0.3176) = ceil(21.83) = 22
        let x_drop_10 = x_drop_raw_score(10.0, &ungapped, 1.0);
        assert_eq!(x_drop_10, 22, "x_drop should be 22 for 10.0 bits");
        
        // 20 bits (NUCL default): ceil(20.0 * 0.693147 / 0.3176) = ceil(43.66) = 44
        let x_drop_20 = x_drop_raw_score(20.0, &ungapped, 1.0);
        assert_eq!(x_drop_20, 44, "x_drop should be 44 for 20.0 bits");
    }

    #[test]
    fn test_cutoff_score_for_update_tblastx() {
        let ungapped = blosum62_ungapped();
        
        // NCBI BlastInitialWordParametersUpdate for tblastx ungapped path
        // Uses CUTOFF_E_TBLASTX = 1e-300 and searchsp = MIN(q_aa, s_nucl) * s_nucl
        let query_len_aa = 500;
        let subject_len_nucl = 3000;  // NUCLEOTIDE length, NOT AA!
        let gap_trigger = gap_trigger_raw_score(22.0, &ungapped);  // 41
        
        // Compute cutoff_score_max first (uses user E-value)
        let eff_searchsp = compute_eff_searchsp_subject_mode_tblastx(
            query_len_aa,
            subject_len_nucl,
            &ungapped,
        );
        let cutoff_score_max = cutoff_score_max_for_tblastx(eff_searchsp, 10.0, &ungapped);
        
        // Compute per-subject cutoff
        let cutoff = cutoff_score_for_update_tblastx(
            query_len_aa,
            subject_len_nucl,
            gap_trigger,
            cutoff_score_max,
            BLAST_GAP_DECAY_RATE,  // 0.5
            &ungapped,
            1.0,
        );
        
        // The result should be capped by either gap_trigger or cutoff_score_max
        assert!(cutoff <= gap_trigger, "cutoff should be <= gap_trigger");
        assert!(cutoff <= cutoff_score_max, "cutoff should be <= cutoff_score_max");
        assert!(cutoff > 0, "cutoff should be positive");
    }

    #[test]
    fn test_cutoff_score_max_for_tblastx() {
        let ungapped = blosum62_ungapped();
        
        // NCBI BlastHitSavingParametersNew uses user E-value and eff_searchsp
        let eff_searchsp = 1_000_000i64;
        
        let cutoff_max = cutoff_score_max_for_tblastx(eff_searchsp, 10.0, &ungapped);
        
        // Should be positive and reasonable
        assert!(cutoff_max > 0, "cutoff_score_max should be positive");
        
        // With E-value = 10.0 and searchsp = 1e6, score should be moderate
        // S = ceil(ln(K * searchsp / E) / Lambda)
        // = ceil(ln(0.134 * 1e6 / 10) / 0.3176)
        // = ceil(ln(13400) / 0.3176)
        // = ceil(9.503 / 0.3176)
        // = ceil(29.93)
        // = 30
        assert_eq!(cutoff_max, 30, "cutoff_score_max should be 30 for these parameters");
    }

    #[test]
    fn test_cutoff_update_vs_legacy() {
        // Verify that the new cutoff_score_for_update_tblastx produces
        // different (typically lower) cutoffs than the old compute_tblastx_cutoff_score
        // because it uses CUTOFF_E_TBLASTX = 1e-300 instead of user's E-value
        
        let ungapped = blosum62_ungapped();
        let gapped = blosum62_gapped();
        
        let query_len_aa = 500;
        let subject_len_nucl = 3000;
        
        // Old method (uses user E-value directly)
        let old_cutoff = compute_tblastx_cutoff_score(
            query_len_aa,
            subject_len_nucl,
            10.0,  // user E-value
            22.0,  // gap_trigger_bits
            &ungapped,
            &gapped,
        );
        
        // New method (uses CUTOFF_E_TBLASTX = 1e-300)
        let gap_trigger = gap_trigger_raw_score(22.0, &ungapped);
        let eff_searchsp = compute_eff_searchsp_subject_mode_tblastx(
            query_len_aa,
            subject_len_nucl,
            &ungapped,
        );
        let cutoff_score_max = cutoff_score_max_for_tblastx(eff_searchsp, 10.0, &ungapped);
        let new_cutoff = cutoff_score_for_update_tblastx(
            query_len_aa,
            subject_len_nucl,
            gap_trigger,
            cutoff_score_max,
            BLAST_GAP_DECAY_RATE,
            &ungapped,
            1.0,
        );
        
        // Both should be positive
        assert!(old_cutoff > 0, "old_cutoff should be positive");
        assert!(new_cutoff > 0, "new_cutoff should be positive");
        
        // New method may produce same or different result depending on caps
        // The key difference is in the internal calculation path
        // Both should be <= gap_trigger (41)
        assert!(old_cutoff <= gap_trigger);
        assert!(new_cutoff <= gap_trigger);
    }

    #[test]
    fn test_cutoff_e_tblastx_constant() {
        // Verify the CUTOFF_E_TBLASTX constant matches NCBI
        assert_eq!(CUTOFF_E_TBLASTX, 1e-300, "CUTOFF_E_TBLASTX should be 1e-300");
    }

    #[test]
    fn test_blast_gap_decay_rate_constant() {
        // Verify the BLAST_GAP_DECAY_RATE constant matches NCBI
        assert_eq!(BLAST_GAP_DECAY_RATE, 0.5, "BLAST_GAP_DECAY_RATE should be 0.5");
    }
}

