//! Cutoff calculations for HSP linking.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c

use crate::stats::KarlinParams;
use crate::stats::sum_statistics::defaults::{GAP_SIZE, OVERLAP_SIZE};

use super::params::{LinkHspCutoffs, BLAST_GAP_DECAY_RATE};

const WINDOW_SIZE: i32 = GAP_SIZE + OVERLAP_SIZE + 1;

/// NCBI BLAST_GAP_PROB (ungapped search default)
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:66
const BLAST_GAP_PROB: f64 = 0.5;

/// NCBI BLAST_Nint - round to nearest integer
/// Reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:437-441
#[inline]
fn blast_nint(x: f64) -> i32 {
    let rounded = if x >= 0.0 { x + 0.5 } else { x - 0.5 };
    rounded as i32
}

/// NCBI CalculateLinkHSPCutoffs - verbatim port
///
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:998-1082
///
/// ```c
/// void
/// CalculateLinkHSPCutoffs(EBlastProgramType program, BlastQueryInfo* query_info,
///    const BlastScoreBlk* sbp, BlastLinkHSPParameters* link_hsp_params,
///    const BlastInitialWordParameters* word_params,
///    Int8 db_length, Int4 subject_length)
/// {
///     Blast_KarlinBlk* kbp = NULL;
///     double gap_prob, gap_decay_rate, x_variable, y_variable;
///     Int4 expected_length, window_size, query_length;
///     Int8 search_sp;
///     const double kEpsilon = 1.0e-9;
///
///     if (!link_hsp_params)
///         return;
///
///     /* Get KarlinBlk for context with smallest lambda (still greater than zero) */
///     s_BlastFindSmallestLambda(sbp->kbp, query_info, &kbp);
///     if (!kbp)
///         return;
///
///     window_size
///         = link_hsp_params->gap_size + link_hsp_params->overlap_size + 1;
///     gap_prob = link_hsp_params->gap_prob = BLAST_GAP_PROB;
///     gap_decay_rate = link_hsp_params->gap_decay_rate;
///     /* Use average query length */
///
///     query_length =
///         (query_info->contexts[query_info->last_context].query_offset +
///         query_info->contexts[query_info->last_context].query_length - 1)
///         / (query_info->last_context + 1);
///
///     if (Blast_SubjectIsTranslated(program) || program == eBlastTypeRpsTblastn) {
///         /* Lengths in subsequent calculations should be on the protein scale */
///         subject_length /= CODON_LENGTH;
///         db_length /= CODON_LENGTH;
///     }
///
///
///     /* Subtract off the expected score. */
///    expected_length = (Int4)BLAST_Nint(log(kbp->K*((double) query_length)*
///                                     ((double) subject_length))/(kbp->H));
///    query_length = query_length - expected_length;
///
///    subject_length = subject_length - expected_length;
///    query_length = MAX(query_length, 1);
///    subject_length = MAX(subject_length, 1);
///
///    /* If this is a database search, use database length, else the single
///       subject sequence length */
///    if (db_length > subject_length) {
///       y_variable = log((double) (db_length)/(double) subject_length)*(kbp->K)/
///          (gap_decay_rate);
///    } else {
///       y_variable = log((double) (subject_length + expected_length)/
///                        (double) subject_length)*(kbp->K)/(gap_decay_rate);
///    }
///
///    search_sp = ((Int8) query_length)* ((Int8) subject_length);
///    x_variable = 0.25*y_variable*((double) search_sp);
///
///    /* To use "small" gaps the query and subject must be "large" compared to
///       the gap size. If small gaps may be used, then the cutoff values must be
///       adjusted for the "bayesian" possibility that both large and small gaps
///       are being checked for. */
///
///    if (search_sp > 8*window_size*window_size) {
///       x_variable /= (1.0 - gap_prob + kEpsilon);
///       link_hsp_params->cutoff_big_gap =
///          (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
///       x_variable = y_variable*(window_size*window_size);
///       x_variable /= (gap_prob + kEpsilon);
///       link_hsp_params->cutoff_small_gap =
///          MAX(word_params->cutoff_score_min,
///              (Int4) floor((log(x_variable)/kbp->Lambda)) + 1);
///    } else {
///       link_hsp_params->cutoff_big_gap =
///          (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
///       /* The following is equivalent to forcing small gap rule to be ignored
///          when linking HSPs. */
///       link_hsp_params->gap_prob = 0;
///       link_hsp_params->cutoff_small_gap = 0;
///    }
///
///    link_hsp_params->cutoff_big_gap *= (Int4)sbp->scale_factor;
///    link_hsp_params->cutoff_small_gap *= (Int4)sbp->scale_factor;
/// }
/// ```
pub fn calculate_link_hsp_cutoffs_ncbi(
    avg_query_length: i32,
    subject_len_nucl: i64,
    db_length: i64, // 0 for -subject mode
    cutoff_score_min: i32,
    scale_factor: f64,
    gap_decay_rate: f64,
    params: &KarlinParams,
) -> LinkHspCutoffs {
    const K_EPSILON: f64 = 1.0e-9;

    // NCBI: gap_prob = link_hsp_params->gap_prob = BLAST_GAP_PROB;
    let mut gap_prob = BLAST_GAP_PROB;

    // NCBI: window_size = link_hsp_params->gap_size + link_hsp_params->overlap_size + 1;
    let window_size = WINDOW_SIZE;

    // NCBI: query_length already provided as avg_query_length
    let mut query_length = avg_query_length;

    // NCBI: if (Blast_SubjectIsTranslated(program)) subject_length /= CODON_LENGTH;
    // For tblastx, subject is translated, so divide by 3
    let mut subject_length = (subject_len_nucl / 3) as i32;
    let db_length = db_length / 3; // Also scale db_length for tblastx

    // NCBI: expected_length = (Int4)BLAST_Nint(log(kbp->K*q*s)/(kbp->H));
    let expected_length = blast_nint(
        (params.k * (query_length as f64) * (subject_length as f64)).ln() / params.h
    );

    // NCBI: query_length = query_length - expected_length;
    query_length -= expected_length;
    // NCBI: subject_length = subject_length - expected_length;
    subject_length -= expected_length;
    // NCBI: query_length = MAX(query_length, 1);
    query_length = query_length.max(1);
    // NCBI: subject_length = MAX(subject_length, 1);
    subject_length = subject_length.max(1);

    // NCBI: y_variable calculation with db_length > subject_length branch
    // For -subject mode, db_length = 0, so we always use the else branch
    let y_variable = if db_length > subject_length as i64 {
        // Database search mode
        ((db_length as f64) / (subject_length as f64)).ln() * params.k / gap_decay_rate
    } else {
        // Subject mode (db_length == 0 or single sequence)
        (((subject_length + expected_length) as f64) / (subject_length as f64)).ln()
            * params.k / gap_decay_rate
    };

    // NCBI: search_sp = ((Int8) query_length)* ((Int8) subject_length);
    let search_sp = (query_length as i64) * (subject_length as i64);

    // NCBI: x_variable = 0.25*y_variable*((double) search_sp);
    let mut x_variable = 0.25 * y_variable * (search_sp as f64);

    let window_sq = (window_size * window_size) as i64;

    let (cutoff_big_gap, cutoff_small_gap, ignore_small_gaps) = if search_sp > 8 * window_sq {
        // NCBI: Large search space - use both small and big gap rules
        // x_variable /= (1.0 - gap_prob + kEpsilon);
        x_variable /= 1.0 - gap_prob + K_EPSILON;

        // cutoff_big_gap = (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
        let cutoff_big = (x_variable.ln() / params.lambda).floor() as i32 + 1;

        // x_variable = y_variable*(window_size*window_size);
        let x_small = y_variable * (window_sq as f64);
        // x_variable /= (gap_prob + kEpsilon);
        let x_small = x_small / (gap_prob + K_EPSILON);

        // cutoff_small_gap = MAX(word_params->cutoff_score_min, floor(log(x)/Lambda)+1);
        let cutoff_small_raw = (x_small.ln() / params.lambda).floor() as i32 + 1;
        let cutoff_small = cutoff_score_min.max(cutoff_small_raw);

        (cutoff_big, cutoff_small, false)
    } else {
        // NCBI: Small search space - disable small gap rule
        // cutoff_big_gap = (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
        let cutoff_big = (x_variable.ln() / params.lambda).floor() as i32 + 1;

        // link_hsp_params->gap_prob = 0;
        gap_prob = 0.0;
        // link_hsp_params->cutoff_small_gap = 0;
        let cutoff_small = 0;

        (cutoff_big, cutoff_small, true)
    };

    // NCBI: cutoff_big_gap *= (Int4)sbp->scale_factor;
    // NCBI: cutoff_small_gap *= (Int4)sbp->scale_factor;
    let scale = scale_factor as i32;

    let final_cutoff_small = cutoff_small_gap * scale;
    let final_cutoff_big = cutoff_big_gap * scale;

    // Debug output for long sequences (600kb+)
    if subject_len_nucl > 600_000 {
        eprintln!("[DEBUG LINKING_CUTOFF] avg_query_length={}, subject_len_nucl={}", avg_query_length, subject_len_nucl);
        eprintln!("[DEBUG LINKING_CUTOFF] query_length_after_adj={}, subject_length_after_adj={}", query_length, subject_length);
        eprintln!("[DEBUG LINKING_CUTOFF] expected_length={}", expected_length);
        eprintln!("[DEBUG LINKING_CUTOFF] search_sp={}, window_sq={}, 8*window_sq={}", search_sp, window_sq, 8 * window_sq);
        eprintln!("[DEBUG LINKING_CUTOFF] search_sp > 8*window_sq: {}", search_sp > 8 * window_sq);
        eprintln!("[DEBUG LINKING_CUTOFF] y_variable={:.6e}", y_variable);
        eprintln!("[DEBUG LINKING_CUTOFF] x_variable_before_div={:.6e}", 0.25 * y_variable * (search_sp as f64));
        if search_sp > 8 * window_sq {
            eprintln!("[DEBUG LINKING_CUTOFF] x_variable_after_div={:.6e}", x_variable);
            eprintln!("[DEBUG LINKING_CUTOFF] x_small={:.6e}", y_variable * (window_sq as f64) / (gap_prob + K_EPSILON));
        }
        eprintln!("[DEBUG LINKING_CUTOFF] cutoff_small_gap_raw={}, cutoff_big_gap_raw={}", cutoff_small_gap, cutoff_big_gap);
        eprintln!("[DEBUG LINKING_CUTOFF] scale_factor={}, final_cutoff_small={}, final_cutoff_big={}", scale_factor, final_cutoff_small, final_cutoff_big);
        eprintln!("[DEBUG LINKING_CUTOFF] gap_prob={}, ignore_small_gaps={}", gap_prob, ignore_small_gaps);
    }

    LinkHspCutoffs {
        cutoff_small_gap: final_cutoff_small,
        cutoff_big_gap: final_cutoff_big,
        gap_prob,
        ignore_small_gaps,
    }
}
