//! Cutoff calculations for HSP linking: TBLASTX's view of `CalculateLinkHSPCutoffs`, whose one
//! port is `core::link_hsps::calculate_link_hsp_cutoffs` (BLASTX's ungapped search uses it too).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c

use crate::core::link_hsps::calculate_link_hsp_cutoffs;
use crate::stats::KarlinParams;

use super::params::LinkHspCutoffs;

/// `CalculateLinkHSPCutoffs` of a TBLASTX search with one subject (`calculate_link_hsp_cutoffs`
/// with a translated subject: the subject's and the database's lengths in nucleotides
/// are divided by 3). The small-gap cutoff is 0 (`LinkHspCutoffs::ignore_small_gaps`) only in
/// the branch of a search space of at most 8 windows, the one that sets `gap_prob` to 0: otherwise
/// it is at least `cutoff_score_min`, the smallest context cutoff (the least of the e-value cutoff, `gap_trigger`
/// and `cutoff_score_max`), each at least 1 (BLAST_Cutoffs raises a cutoff from 1 only,
/// blast_stat.c:4090-4130; the 22-bit gap trigger of TBLASTX's matrices is positive).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1451-1455
pub fn calculate_link_hsp_cutoffs_ncbi(
    avg_query_length: i32,
    subject_len_nucl: i64,
    db_length: i64, // Full nucleotide total, including local subjects
    cutoff_score_min: i32,
    params: &KarlinParams,
) -> LinkHspCutoffs {
    let link = calculate_link_hsp_cutoffs(
        params,
        avg_query_length,
        subject_len_nucl,
        db_length,
        true,
        || cutoff_score_min,
    );
    let cutoffs = LinkHspCutoffs {
        cutoff_small_gap: link.cutoff_small_gap,
        cutoff_big_gap: link.cutoff_big_gap,
        gap_prob: link.gap_prob,
    };
    // NCBI computes these values without unconditional stderr diagnostics.
    if crate::utils::xenv::debug_cutoffs_is_some() && subject_len_nucl > 600_000 {
        eprintln!(
            "[DEBUG LINKING_CUTOFF] avg_query_length={}, subject_len_nucl={}, db_length={}, cutoff_score_min={}",
            avg_query_length, subject_len_nucl, db_length, cutoff_score_min
        );
        eprintln!(
            "[DEBUG LINKING_CUTOFF] cutoff_small_gap={}, cutoff_big_gap={}, gap_prob={}, ignore_small_gaps={}",
            cutoffs.cutoff_small_gap, cutoffs.cutoff_big_gap, cutoffs.gap_prob, cutoffs.ignore_small_gaps()
        );
    }
    cutoffs
}

#[cfg(test)]
#[path = "cutoffs_database_cutoff_tests.rs"]
mod database_cutoff_tests;
