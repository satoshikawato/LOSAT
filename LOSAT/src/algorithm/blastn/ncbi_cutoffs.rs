//! NCBI BLAST cutoff score calculation for blastn ungapped extensions
//!
//! This module implements NCBI BLAST's cutoff score calculation for blastn ungapped extensions.
//! All functions are direct ports from NCBI C code with original comments preserved.
//!
//! Reference files:
//! - ncbi-blast/c++/src/algo/blast/core/blast_parameters.c
//! - ncbi-blast/c++/include/algo/blast/core/blast_options.h

use crate::stats::KarlinParams;
use crate::stats::length_adjustment::compute_length_adjustment_ncbi;

/// ln(2) constant used in NCBI BLAST
/// Reference: ncbi-blast/c++/include/algo/blast/core/ncbi_math.h
/// #define NCBIMATH_LN2 0.69314718055994530941723212145818
const NCBIMATH_LN2: f64 = 0.69314718055994530941723212145818;

/// Default gap trigger bit score for nucleotide searches
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:140
/// #define BLAST_GAP_TRIGGER_NUCL 27.0
pub const GAP_TRIGGER_BIT_SCORE_NUCL: f64 = 27.0;

/// Default cutoff E-value for blastn ungapped extension
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:76
/// #define CUTOFF_E_BLASTN 0.05
pub const CUTOFF_E_BLASTN: f64 = 0.05;

/// Smallest float to avoid floating point exception
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4049
const K_SMALL_FLOAT: f64 = 1.0e-297;

/// Calculate gap_trigger raw score from bit score using UNGAPPED Karlin params.
///
/// NCBI reference (verbatim from blast_parameters.c:343-344):
/// ```c
/// gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + 
///                          kbp->logK) / kbp->Lambda);
/// ```
///
/// CRITICAL: Uses UNGAPPED params (kbp_std), NOT gapped params (kbp_gap).
///
/// # Arguments
/// * `gap_trigger_bits` - Gap trigger in bits (NCBI default: 27.0 for blastn)
/// * `ungapped_params` - UNGAPPED Karlin-Altschul parameters (kbp_std)
///
/// # Returns
/// gap_trigger (raw score)
pub fn gap_trigger_raw_score(
    gap_trigger_bits: f64,
    ungapped_params: &KarlinParams,
) -> i32 {
    // NCBI: gap_trigger = (Int4)((gap_trigger_bits * NCBIMATH_LN2 + logK) / Lambda)
    let log_k = ungapped_params.k.ln();
    let numerator = gap_trigger_bits * NCBIMATH_LN2 + log_k;
    let gap_trigger = (numerator / ungapped_params.lambda) as i32;
    
    gap_trigger
}

/// Calculate cutoff_score_max from E-value using GAPPED Karlin params.
///
/// NCBI reference (verbatim from blast_parameters.c:943):
/// ```c
/// BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
/// ```
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
    // NCBI: Check for invalid Karlin parameters (blast_stat.c:4101-4102)
    // if (kbp->Lambda == -1. || kbp->K == -1. || kbp->H == -1.) return 1;
    if gapped_params.lambda < 0.0 || gapped_params.k < 0.0 || gapped_params.h < 0.0 {
        return 1; // Return error value (matches NCBI behavior)
    }
    
    // NCBI: es = 1 (default), then calculate if e > 0.
    // NCBI: if (es > s) *S = es; (pick larger between initial S and calculated es)
    // Reference: blast_stat.c:4108-4129
    let mut s = 1; // Initial value (matches NCBI: es = 1)
    
    // NCBI: if (e > 0.)
    if evalue > 0.0 {
        // NCBI: E = MAX(E, kSmallFloat)
        let e = evalue.max(K_SMALL_FLOAT);
        
        // NCBI: es = BlastKarlinEtoS_simple(e, kbp, searchsp)
        // NCBI: S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda ))
        let searchsp = eff_searchsp as f64;
        let k_times_searchsp = gapped_params.k * searchsp;
        let k_times_searchsp_over_e = k_times_searchsp / e;
        let log_value = k_times_searchsp_over_e.ln();
        let score_before_ceil = log_value / gapped_params.lambda;
        let es = score_before_ceil.ceil() as i32;
        
        // NCBI: if (es > s) *S = es; (pick larger)
        if es > s {
            s = es;
        }
    }
    
    s
}

/// Calculate cutoff_score for blastn ungapped extension.
///
/// NCBI reference (verbatim from blast_parameters.c:368-373):
/// ```c
/// } else {
///     new_cutoff = gap_trigger;
/// }
/// new_cutoff *= (Int4)sbp->scale_factor;
/// new_cutoff = MIN(new_cutoff, 
///                  hit_params->cutoffs[context].cutoff_score_max);
/// curr_cutoffs->cutoff_score = new_cutoff;
/// ```
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
    // NCBI: new_cutoff = gap_trigger (gapped mode)
    let mut new_cutoff = gap_trigger;
    
    // NCBI: new_cutoff *= (Int4)sbp->scale_factor
    new_cutoff = (new_cutoff as f64 * scale_factor) as i32;
    
    // NCBI: new_cutoff = MIN(new_cutoff, cutoff_score_max)
    new_cutoff.min(cutoff_score_max)
}

/// Calculate effective lengths for -subject mode (single subject as database).
///
/// This computes BOTH length_adjustment and eff_searchsp in a single call,
/// matching `BLAST_CalcEffLengths` from blast_setup.c:821-847.
///
/// NCBI reference (verbatim from blast_setup.c:836-843):
/// ```c
/// Int8 effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
/// if (effective_db_length <= 0)
///     effective_db_length = 1;
/// effective_search_space = effective_db_length * (query_length - length_adjustment);
/// ```
///
/// # Arguments
/// * `query_len` - Query length (for both strands if applicable, already * 2)
/// * `subject_len` - Subject length
/// * `gapped_params` - GAPPED Karlin-Altschul parameters (kbp_gap)
///
/// # Returns
/// `(length_adjustment, eff_searchsp)` tuple
pub fn compute_eff_lengths_subject_mode_blastn(
    query_len: i64,
    subject_len: i64,
    gapped_params: &KarlinParams,
) -> (i64, i64) {
    // NCBI: db_num_seqs = 1 for -subject mode
    let db_num_seqs: i64 = 1;
    
    // NCBI blast_setup.c:821-824: BLAST_ComputeLengthAdjustment(...)
    let result = compute_length_adjustment_ncbi(
        query_len,
        subject_len,
        db_num_seqs,
        gapped_params,
    );
    let length_adjustment = result.length_adjustment;
    
    // NCBI blast_setup.c:836: effective_db_length = db_length - (db_num_seqs * length_adjustment)
    let effective_db_length = (subject_len - db_num_seqs * length_adjustment).max(1);
    
    // NCBI blast_setup.c:842-843: effective_search_space = effective_db_length * (query_length - length_adjustment)
    let effective_query_length = (query_len - length_adjustment).max(1);
    let eff_searchsp = effective_db_length * effective_query_length;
    
    (length_adjustment, eff_searchsp)
}

/// Calculate cutoff_score for ungapped extension when !gapped_calculation || matrix_only_scoring.
///
/// NCBI reference (verbatim from blast_parameters.c:348-367):
/// ```c
/// if (!gapped_calculation || sbp->matrix_only_scoring) {
///     double cutoff_e = s_GetCutoffEvalue(program_number);
///     Int4 query_length = query_info->contexts[context].query_length;
///     if (program_number == eBlastTypeBlastn ||
///         program_number == eBlastTypeMapping)
///        query_length *= 2;
///     kbp = kbp_array[context];
///     BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp, 
///                   MIN((Uint8)subj_length, 
///                       (Uint8)query_length)*((Uint8)subj_length), 
///                   TRUE, gap_decay_rate);
///     if (program_number != eBlastTypeBlastn)  
///        new_cutoff = MIN(new_cutoff, gap_trigger);
/// }
/// ```
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
/// * `cutoff_score_max` - Maximum cutoff score from E-value (from hit_params->cutoffs[context].cutoff_score_max)
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
    // NCBI: cutoff_e = s_GetCutoffEvalue(program_number) = CUTOFF_E_BLASTN = 0.05
    let cutoff_e = CUTOFF_E_BLASTN;
    
    // NCBI: searchsp = MIN(subj_length, query_length) * subj_length (length adjustmentなし)
    let min_len = subject_len.min(query_len);
    let searchsp = min_len * subject_len;
    
    // NCBI: Check for invalid Karlin parameters (blast_stat.c:4101-4102)
    // if (kbp->Lambda == -1. || kbp->K == -1. || kbp->H == -1.) return 1;
    if ungapped_params.lambda < 0.0 || ungapped_params.k < 0.0 || ungapped_params.h < 0.0 {
        return 1; // Return error value (matches NCBI behavior)
    }
    
    // NCBI: BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp, searchsp, TRUE, gap_decay_rate)
    // NCBI: es = 1 (default), then calculate if e > 0.
    // NCBI: if (es > s) *S = es; (pick larger between initial S and calculated es)
    // Reference: blast_stat.c:4108-4129
    let mut new_cutoff = 1; // Initial value (matches NCBI: es = 1, s = *S = 1)
    
    // NCBI: if (e > 0.)
    if cutoff_e > 0.0 {
        // NCBI: E = MAX(E, kSmallFloat)
        let e = cutoff_e.max(K_SMALL_FLOAT);
        
        // NCBI: When dodecay=TRUE, apply gap decay divisor to E-value before conversion
        let e_adjusted = if gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
            // NCBI: e *= BLAST_GapDecayDivisor(gap_decay_rate, 1)
            // BLAST_GapDecayDivisor(decayrate, 1) = (1 - decayrate) * decayrate^0 = (1 - decayrate)
            let divisor = 1.0 - gap_decay_rate;
            e * divisor
        } else {
            e
        };
        
        // NCBI: es = BlastKarlinEtoS_simple(e, kbp, searchsp)
        // NCBI: S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda ))
        let searchsp_f64 = searchsp as f64;
        let k_times_searchsp = ungapped_params.k * searchsp_f64;
        let k_times_searchsp_over_e = k_times_searchsp / e_adjusted;
        let log_value = k_times_searchsp_over_e.ln();
        let score_before_ceil = log_value / ungapped_params.lambda;
        let es = score_before_ceil.ceil() as i32;
        
        // NCBI: if (es > s) *S = es; (pick larger)
        if es > new_cutoff {
            new_cutoff = es;
        }
    }
    
    // NCBI: if (program_number != eBlastTypeBlastn) new_cutoff = MIN(new_cutoff, gap_trigger)
    // For blastn, this check is skipped (always true for blastn)
    
    // NCBI: new_cutoff *= (Int4)sbp->scale_factor
    new_cutoff = (new_cutoff as f64 * scale_factor) as i32;
    
    // NCBI: new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max)
    // Reference: blast_parameters.c:372-373
    new_cutoff = new_cutoff.min(cutoff_score_max);
    
    new_cutoff
}

/// High-level function to compute ungapped extension cutoff for blastn (gapped mode).
///
/// This combines all the NCBI cutoff calculation steps for blastn in gapped mode:
/// 1. Compute gap_trigger using UNGAPPED params
/// 2. Compute eff_searchsp with length adjustment using GAPPED params
/// 3. Compute cutoff_score_max using GAPPED params and eff_searchsp
/// 4. Return MIN(gap_trigger * scale_factor, cutoff_score_max)
///
/// NCBI Reference (blast_parameters.c:368-374):
/// For blastn in gapped mode (gapped_calculation = TRUE):
/// ```c
/// } else {
///     new_cutoff = gap_trigger;
/// }
/// new_cutoff *= (Int4)sbp->scale_factor;
/// new_cutoff = MIN(new_cutoff, 
///                  hit_params->cutoffs[context].cutoff_score_max);
/// curr_cutoffs->cutoff_score = new_cutoff;
/// ```
///
/// NCBI Reference (blast_parameters.c:926):
/// searchsp = query_info->contexts[context].eff_searchsp;  // length adjustment included
///
/// # Arguments
/// * `query_len` - Query length (for both strands if applicable, already * 2)
/// * `subject_len` - Subject length
/// * `evalue_threshold` - E-value threshold (typically 10.0, for cutoff_score_max calculation)
/// * `gap_trigger_bits` - Gap trigger in bits (typically 27.0)
/// * `ungapped_params` - UNGAPPED Karlin-Altschul parameters (kbp_std)
/// * `gapped_params` - GAPPED Karlin-Altschul parameters (kbp_gap)
/// * `scale_factor` - Scale factor (typically 1.0 for standard scoring)
///
/// # Returns
/// Cutoff score for ungapped extension (raw score)
pub fn compute_blastn_cutoff_score(
    query_len: i64,
    subject_len: i64,
    evalue_threshold: f64,
    gap_trigger_bits: f64,
    ungapped_params: &KarlinParams,
    gapped_params: &KarlinParams,
    scale_factor: f64,
) -> i32 {
    // Debug mode: LOSAT_DEBUG_CUTOFFS=1 to show detailed calculation
    let debug_cutoffs = std::env::var("LOSAT_DEBUG_CUTOFFS").is_ok();

    // Step 1: Compute gap_trigger using UNGAPPED params
    let gap_trigger = gap_trigger_raw_score(gap_trigger_bits, ungapped_params);

    // Step 2: Compute eff_searchsp with length adjustment using GAPPED params
    // NCBI reference: blast_setup.c:836-843
    // eff_searchsp = (subject_len - length_adj) * (query_len - length_adj)
    let (length_adj, eff_searchsp) = compute_eff_lengths_subject_mode_blastn(
        query_len,
        subject_len,
        gapped_params,
    );

    // Step 3: Compute cutoff_score_max using GAPPED params
    // NCBI reference: blast_parameters.c:943
    // BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
    let cutoff_score_max = cutoff_score_max_from_evalue(
        evalue_threshold,
        eff_searchsp,
        gapped_params,
    );

    // Step 4: Return MIN(gap_trigger * scale_factor, cutoff_score_max)
    // NCBI reference: blast_parameters.c:368-373
    // For blastn in gapped mode: new_cutoff = gap_trigger, then MIN with cutoff_score_max
    let final_cutoff = cutoff_score_for_ungapped_extension(gap_trigger, cutoff_score_max, scale_factor);

    if debug_cutoffs {
        eprintln!(
            "[CUTOFF] q_len={}, s_len={}, len_adj={}, eff_searchsp={}, gap_trigger={}, cutoff_score_max={}, final={}",
            query_len, subject_len, length_adj, eff_searchsp, gap_trigger, cutoff_score_max, final_cutoff
        );
    }

    final_cutoff
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gap_trigger_raw_score() {
        // Test with blastn default params (reward=2, penalty=-3, gap_open=5, gap_extend=2)
        // Ungapped params: lambda=0.625, k=0.41
        let ungapped_params = KarlinParams {
            lambda: 0.625,
            k: 0.41,
            h: 0.78,
            alpha: 0.8,
            beta: -2.0,
        };
        
        let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE_NUCL, &ungapped_params);
        // Expected: (27.0 * 0.6931471805599453 + ln(0.41)) / 0.625
        // = (18.714973875118523 - 0.8915981196417836) / 0.625
        // = 17.82337575547674 / 0.625
        // = 28.517401208762784
        // ≈ 28 or 29
        assert!(gap_trigger >= 28 && gap_trigger <= 29);
    }

    #[test]
    fn test_cutoff_score_max_from_evalue() {
        // Test with blastn gapped params (reward=2, penalty=-3, gap_open=5, gap_extend=2)
        // Gapped params: lambda=0.625, k=0.41 (same as ungapped for blastn)
        let gapped_params = KarlinParams {
            lambda: 0.625,
            k: 0.41,
            h: 0.78,
            alpha: 0.8,
            beta: -2.0,
        };
        
        let query_len = 1000;
        let subject_len = 10000;
        let eff_searchsp = query_len * subject_len;
        let evalue = 10.0;
        
        let cutoff_score_max = cutoff_score_max_from_evalue(evalue, eff_searchsp, &gapped_params);
        // Should be a positive value
        assert!(cutoff_score_max > 0);
    }

    #[test]
    fn test_compute_blastn_cutoff_score_ungapped() {
        // Test with blastn ungapped params (reward=2, penalty=-3)
        // Ungapped params: lambda=0.625, k=0.41
        let ungapped_params = KarlinParams {
            lambda: 0.625,
            k: 0.41,
            h: 0.78,
            alpha: 0.8,
            beta: -2.0,
        };
        
        let query_len = 1000i64 * 2; // Both strands
        let subject_len = 10000i64;
        let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE_NUCL, &ungapped_params);
        let scale_factor = 1.0;
        let gap_decay_rate = 0.0; // No gap decay for ungapped
        let cutoff_score_max = 10000; // Dummy value for testing
        
        let cutoff_score = compute_blastn_cutoff_score_ungapped(
            query_len,
            subject_len,
            gap_trigger,
            &ungapped_params,
            scale_factor,
            gap_decay_rate,
            cutoff_score_max,
        );
        
        // Should be a positive value
        assert!(cutoff_score > 0);
        // Should use CUTOFF_E_BLASTN = 0.05
        // searchsp = MIN(2000, 10000) * 10000 = 2000 * 10000 = 20,000,000
        // Expected score should be calculated from CUTOFF_E_BLASTN = 0.05
        // Should be MIN(new_cutoff, cutoff_score_max)
        assert!(cutoff_score <= cutoff_score_max);
    }

    #[test]
    fn test_compute_blastn_cutoff_score_gapped() {
        // Test with blastn gapped mode (default)
        let ungapped_params = KarlinParams {
            lambda: 0.625,
            k: 0.41,
            h: 0.78,
            alpha: 0.8,
            beta: -2.0,
        };
        let gapped_params = ungapped_params.clone(); // Same for blastn
        
        let query_len = 1000i64 * 2; // Both strands
        let subject_len = 10000i64;
        let evalue_threshold = 10.0;
        let scale_factor = 1.0;
        
        let cutoff_score = compute_blastn_cutoff_score(
            query_len,
            subject_len,
            evalue_threshold,
            GAP_TRIGGER_BIT_SCORE_NUCL,
            &ungapped_params,
            &gapped_params,
            scale_factor,
        );
        
        // Should be a positive value
        assert!(cutoff_score > 0);
        // Should be MIN(gap_trigger, cutoff_score_max)
        let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE_NUCL, &ungapped_params);
        assert!(cutoff_score <= gap_trigger);
    }
}

