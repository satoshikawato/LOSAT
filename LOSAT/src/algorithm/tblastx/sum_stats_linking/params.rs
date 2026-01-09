//! Linking parameters and helper functions.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c

use crate::stats::KarlinParams;

/// NCBI BLAST_GAP_DECAY_RATE (ungapped search default)
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:68
pub const BLAST_GAP_DECAY_RATE: f64 = 0.5;

/// Parameters for HSP linking cutoffs
/// These are the outputs from calculate_link_hsp_cutoffs_ncbi()
#[derive(Debug, Clone, Copy)]
pub struct LinkHspCutoffs {
    /// Cutoff for small gap linking (index 0)
    pub cutoff_small_gap: i32,
    /// Cutoff for big gap linking (index 1)
    pub cutoff_big_gap: i32,
    /// Gap probability (may be set to 0 for small search spaces)
    pub gap_prob: f64,
    /// Whether to ignore small gaps (when gap_prob=0)
    pub ignore_small_gaps: bool,
}

/// Parameters required for NCBI-style linking
/// These are computed once per subject and passed to linking
#[derive(Debug, Clone)]
pub struct LinkingParams {
    /// Average query length in AA (NCBI CalculateLinkHSPCutoffs formula)
    pub avg_query_length: i32,
    /// Subject length in nucleotides
    pub subject_len_nucl: i64,
    /// Minimum cutoff score across all contexts for this subject
    pub cutoff_score_min: i32,
    /// Scale factor (typically 1.0 for standard BLOSUM62)
    pub scale_factor: f64,
    /// Gap decay rate
    pub gap_decay_rate: f64,
}

impl Default for LinkingParams {
    fn default() -> Self {
        Self {
            avg_query_length: 100,
            subject_len_nucl: 300,
            cutoff_score_min: 0,
            scale_factor: 1.0,
            gap_decay_rate: BLAST_GAP_DECAY_RATE,
        }
    }
}

/// Find Karlin parameters with smallest Lambda from a list of context parameters
///
/// This replicates NCBI's s_BlastFindSmallestLambda (blast_parameters.c:92-112)
/// which finds the Karlin block with the smallest lambda value.
///
/// For tblastx, all contexts typically use kbp_ideal (same lambda), but we
/// maintain this logic for NCBI structural parity.
pub fn find_smallest_lambda_params(params_list: &[KarlinParams]) -> Option<KarlinParams> {
    params_list
        .iter()
        .filter(|p| p.lambda > 0.0)
        .min_by(|a, b| a.lambda.partial_cmp(&b.lambda).unwrap_or(std::cmp::Ordering::Equal))
        .cloned()
}

/// Find smallest Lambda value from a list of context parameters
/// Returns the minimum lambda value (f64::MAX if list is empty or all invalid)
pub fn find_smallest_lambda(params_list: &[KarlinParams]) -> f64 {
    params_list
        .iter()
        .filter(|p| p.lambda > 0.0)
        .map(|p| p.lambda)
        .fold(f64::MAX, f64::min)
}

/// Compute NCBI-style average query length from nucleotide lengths
///
/// This replicates NCBI's SetupQueryInfo_OMF + CalculateLinkHSPCutoffs logic:
/// - For each query, compute 6 context lengths using BLAST_GetTranslatedProteinLength
/// - Build context offsets using s_QueryInfo_SetContext pattern
/// - Apply NCBI's average formula: (last_offset + last_length - 1) / num_contexts
///
/// Reference:
/// - blast_util.c:BLAST_GetTranslatedProteinLength
/// - blast_setup_cxx.cpp:s_QueryInfo_SetContext
/// - blast_parameters.c:1023-1026
pub fn compute_avg_query_length_ncbi(query_nucl_lengths: &[usize]) -> i32 {
    if query_nucl_lengths.is_empty() {
        return 1;
    }

    // Build context offsets and lengths like NCBI
    // Each query has 6 contexts (frames 0-5 in NCBI, which map to frames 1,2,3,-1,-2,-3)
    let mut contexts: Vec<(i32, i32)> = Vec::new(); // (query_offset, query_length)

    for &nucl_len in query_nucl_lengths {
        for context in 0..6u32 {
            // NCBI BLAST_GetTranslatedProteinLength:
            // return (nucleotide_length - context % CODON_LENGTH) / CODON_LENGTH;
            let prot_len = if nucl_len == 0 || nucl_len <= (context % 3) as usize {
                0
            } else {
                ((nucl_len - (context % 3) as usize) / 3) as i32
            };

            // NCBI s_QueryInfo_SetContext offset calculation:
            // Uint4 shift = prev_len ? prev_len + 1 : 0;
            // qinfo->contexts[index].query_offset = prev_loc + shift;
            let (prev_offset, prev_len) = contexts.last().copied().unwrap_or((0, 0));
            let shift = if prev_len > 0 { prev_len + 1 } else { 0 };
            let new_offset = prev_offset + shift;

            contexts.push((new_offset, prot_len));
        }
    }

    if contexts.is_empty() {
        return 1;
    }

    // NCBI average formula: (last_offset + last_length - 1) / num_contexts
    let (last_offset, last_length) = contexts.last().copied().unwrap();
    let num_contexts = contexts.len() as i32;

    let avg = (last_offset + last_length - 1) / num_contexts;
    avg.max(1)
}
