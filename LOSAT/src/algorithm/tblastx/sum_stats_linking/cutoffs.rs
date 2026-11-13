//! Cutoff calculations for HSP linking: TBLASTX's view of `CalculateLinkHSPCutoffs`, whose one
//! port is `core::link_hsps::EvenGapLinkParams::ungapped` (BLASTX's ungapped search uses it too).
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c

use crate::core::link_hsps::EvenGapLinkParams;
use crate::stats::KarlinParams;

use super::params::LinkHspCutoffs;

/// `CalculateLinkHSPCutoffs` of a TBLASTX search with one subject (`EvenGapLinkParams::ungapped`
/// with a translated query and subject: the subject's and the database's lengths in nucleotides
/// are divided by 3). `ignore_small_gaps` is the branch of a search space of at most 8 windows,
/// the one that sets `gap_prob` to 0.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1451-1455
/// ```c
///       if (hit_params->link_hsp_params && !kNucleotide &&
///           !gapped_calculation) {
///           CalculateLinkHSPCutoffs(program_number, query_info, sbp,
///             hit_params->link_hsp_params, word_params, db_length,
///             seq_arg.seq->length);
/// ```
pub fn calculate_link_hsp_cutoffs_ncbi(
    avg_query_length: i32,
    subject_len_nucl: i64,
    db_length: i64, // Full nucleotide total, including local subjects
    cutoff_score_min: i32,
    params: &KarlinParams,
) -> LinkHspCutoffs {
    let link = EvenGapLinkParams::ungapped(
        params,
        avg_query_length,
        subject_len_nucl,
        db_length,
        true,
        cutoff_score_min,
        true,
    );
    let cutoffs = LinkHspCutoffs {
        cutoff_small_gap: link.cutoff_small_gap,
        cutoff_big_gap: link.cutoff_big_gap,
        gap_prob: link.gap_prob,
        ignore_small_gaps: link.gap_prob == 0.0,
    };
    // NCBI computes these values without unconditional stderr diagnostics.
    if crate::utils::xenv::debug_cutoffs_is_some() && subject_len_nucl > 600_000 {
        eprintln!(
            "[DEBUG LINKING_CUTOFF] avg_query_length={}, subject_len_nucl={}, db_length={}, cutoff_score_min={}",
            avg_query_length, subject_len_nucl, db_length, cutoff_score_min
        );
        eprintln!(
            "[DEBUG LINKING_CUTOFF] cutoff_small_gap={}, cutoff_big_gap={}, gap_prob={}, ignore_small_gaps={}",
            cutoffs.cutoff_small_gap, cutoffs.cutoff_big_gap, cutoffs.gap_prob, cutoffs.ignore_small_gaps
        );
    }
    cutoffs
}

#[cfg(test)]
mod database_cutoff_tests {
    use super::*;
    // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:998-1082
    // if (db_length > subject_length) { y_variable = log((double)db_length / subject_length)*kbp->K/gap_decay_rate; }
    // Expected values evaluated from the extracted NCBI C function.
    #[test]
    fn complete_database_changes_linking_cutoff() {
        let params = KarlinParams {
            lambda: 0.3176,
            k: 0.134,
            h: 0.4012,
            alpha: 0.7916,
            beta: -3.2,
        };
        for (q, subject, db, small, big, probability) in [
            (300, 900, 0, 30, 22, 0.5),
            (300, 900, 2700, 30, 30, 0.5),
            (30, 90, 90, 0, 8, 0.0),
            (30, 90, 2700, 0, 14, 0.0),
        ] {
            let c = calculate_link_hsp_cutoffs_ncbi(q, subject, db, 30, &params);
            assert_eq!(
                (c.cutoff_small_gap, c.cutoff_big_gap, c.gap_prob),
                (small, big, probability)
            );
        }
    }
}
