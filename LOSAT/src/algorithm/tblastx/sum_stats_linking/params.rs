//! Linking parameters and helper functions.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c

use crate::stats::KarlinParams;

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
}

impl LinkHspCutoffs {
    /// Whether the linking ignores small gaps: a small-gap cutoff of 0, the rule of
    /// `core::link_hsps`. In TBLASTX that is the search space of at most 8 windows, where
    /// gap_prob is 0 too (`calculate_link_hsp_cutoffs_ncbi`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:495
    /// ```c
    ///    ignore_small_gaps = (cutoff[0] == 0);
    /// ```
    pub fn ignore_small_gaps(&self) -> bool {
        self.cutoff_small_gap == 0
    }
}

/// NCBI reference: c++/src/algo/blast/core/blast_engine.c:1448-1455;
/// c++/src/algo/blast/core/link_hsps.c:468-493
/// CalculateLinkHSPCutoffs(..., db_length, subject->length);
/// `cutoff[0] = link_hsp_params->cutoff_small_gap;`
/// `cutoff[1] = link_hsp_params->cutoff_big_gap;`
/// Per-subject state computed before search and shared by both linking passes.
#[derive(Debug, Clone)]
pub struct LinkingParams {
    pub subject_len_nucl: i64,
    pub gap_decay_rate: f64,
    pub cutoffs: LinkHspCutoffs,
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
        .min_by(|a, b| {
            a.lambda
                .partial_cmp(&b.lambda)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .cloned()
}

/// Compute NCBI-style average query length from nucleotide lengths
///
/// This replicates NCBI's SetupQueryInfo_OMF + CalculateLinkHSPCutoffs logic:
/// - For each query, compute 6 context lengths using BLAST_GetTranslatedProteinLength, 0 for the
///   contexts of the strand that `-strand` does not search (`searched`: plus, minus)
/// - Build context offsets using s_QueryInfo_SetContext pattern
/// - Apply NCBI's average formula: (last_offset + last_length - 1) / num_contexts
///
/// The contexts of the strand not searched stay in the count, so with `-strand plus` or `minus`
/// the average is about half that of both strands.
///
/// Reference:
/// - blast_util.c:BLAST_GetTranslatedProteinLength
/// - blast_setup_cxx.cpp:s_QueryInfo_SetContext
/// - NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:200-215 (the
///   context lengths of `-strand plus` and `minus`)
/// - blast_parameters.c:1023-1026
pub fn compute_avg_query_length_ncbi(query_nucl_lengths: &[usize], searched: (bool, bool)) -> i32 {
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
            let strand_searched = if context < 3 { searched.0 } else { searched.1 };
            let prot_len = if strand_searched { prot_len } else { 0 };

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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:1023-1026
    // No floor: a batch whose contexts hold fewer letters than it has contexts averages 0 (C's
    // division truncates toward zero, as Rust's), and `calculate_link_hsp_cutoffs` takes NCBI's
    // expected length of such a query (0).
    (last_offset + last_length - 1) / num_contexts
}

#[cfg(test)]
mod tests {
    use super::compute_avg_query_length_ncbi;

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:200-215 (the contexts
    // of the strand not searched have length 0) and c++/src/algo/blast/core/blast_parameters.c:
    // 1023-1026 (the average over every context)
    #[test]
    fn average_query_length_counts_the_contexts_of_the_strand_not_searched() {
        // 300 bases: 100, 99 and 99 residues a strand. Both strands: offsets 0, 101, 201, 301, 402
        // and 502, (502 + 99 - 1) / 6 = 100.
        assert_eq!(compute_avg_query_length_ncbi(&[300], (true, true)), 100);
        // Plus: contexts 3 to 5 have no letters, at offset 301: (301 + 0 - 1) / 6 = 50.
        assert_eq!(compute_avg_query_length_ncbi(&[300], (true, false)), 50);
        // Minus: contexts 0 to 2 at offset 0, then 0, 101 and 201: (201 + 99 - 1) / 6 = 49.
        assert_eq!(compute_avg_query_length_ncbi(&[300], (false, true)), 49);
        // A second query of 30 bases: its empty contexts start after the first query's last.
        assert_eq!(compute_avg_query_length_ncbi(&[300, 30], (true, true)), 55);
        assert_eq!(compute_avg_query_length_ncbi(&[300, 30], (true, false)), 27);
        assert_eq!(compute_avg_query_length_ncbi(&[300, 30], (false, true)), 27);
        // No letters in any context: (0 + 0 - 1) / 6, 0 in C and in Rust.
        assert_eq!(compute_avg_query_length_ncbi(&[2], (false, true)), 0);
        // One query of 1,200 bases and 300 of 3 bases, one strand: 1,200 + 300 letters and the
        // separators in 1,806 contexts, an average of 0 (both strands: 1).
        let mut lengths = vec![1200];
        lengths.extend([3; 300]);
        assert_eq!(compute_avg_query_length_ncbi(&lengths, (true, false)), 0);
        assert_eq!(compute_avg_query_length_ncbi(&lengths, (false, true)), 0);
        assert_eq!(compute_avg_query_length_ncbi(&lengths, (true, true)), 1);
    }
}
