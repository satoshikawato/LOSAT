//! TBLASTN local-subject statistical inputs before preliminary search.

use crate::stats::length_adjustment::compute_length_adjustment_ncbi;
use crate::stats::tables::KarlinParams;

// NCBI c++/include/algo/blast/core/blast_query_info.h:60-76:
// Int8 eff_searchsp; Int4 length_adjustment; per query context.
// NCBI c++/src/algo/blast/core/blast_engine.c:1446-1467:
// stat_length = seq_arg.seq->length; translated subjects use stat_length /= 3.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct LocalContextLength {
    pub length_adjustment: i64,
    pub eff_searchsp: i64,
    pub subject_stat_length: i64,
}

// NCBI c++/src/algo/blast/core/blast_setup.c:729-735,770-847:
// if (Blast_SubjectIsTranslated(program_number)) db_length = db_length/3;
// BLAST_ComputeLengthAdjustment(kbp->K, kbp->logK, alpha/kbp->Lambda,
//                               beta, query_length, db_length, db_num_seqs,
//                               &length_adjustment);
// Int8 effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
// if (effective_db_length <= 0) effective_db_length = 1;
// effective_search_space = effective_db_length * (query_length - length_adjustment);
// NCBI c++/src/algo/blast/core/blast_engine.c:1434-1442,1446-1467:
// BLAST_OneSubjectUpdateParameters(..., seq_arg.seq->length, ...);
// stat_length = seq_arg.seq->length;
// if (Blast_SubjectIsTranslated(program_number)) stat_length /= CODON_LENGTH;
#[allow(dead_code)] // Stage D's public path remains gated until C/D/E are complete.
pub(super) fn local_subject_effective_lengths(
    query_contexts: &[(usize, bool)],
    subject_nt_length: usize,
    gapped_params: &[KarlinParams],
) -> Vec<LocalContextLength> {
    assert_eq!(query_contexts.len(), gapped_params.len());
    let db_length = (subject_nt_length / 3) as i64;
    query_contexts
        .iter()
        .zip(gapped_params)
        .map(|(&(query_length, is_valid), params)| {
            // NCBI blast_setup.c:802-847:
            // if (query_info->contexts[index].is_valid &&
            //     ((query_length = query_info->contexts[index].query_length) > 0))
            //     BLAST_ComputeLengthAdjustment(...);
            // query_info->contexts[index].eff_searchsp = effective_search_space;
            if !is_valid || query_length == 0 {
                return LocalContextLength {
                    length_adjustment: 0,
                    eff_searchsp: 0,
                    subject_stat_length: db_length,
                };
            }
            let query_length = query_length as i64;
            let adjustment = compute_length_adjustment_ncbi(query_length, db_length, 1, params)
                .length_adjustment;
            let effective_db_length = (db_length - adjustment).max(1);
            LocalContextLength {
                length_adjustment: adjustment,
                eff_searchsp: effective_db_length * (query_length - adjustment),
                subject_stat_length: db_length,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::ScoringMatrix;
    use crate::stats::tables::lookup_protein_params_gapped;

    // Pinned NCBI blast_setup.c:729-735,770-847 and the actual
    // BLAST_LinkHsps/Blast_HSPListGetEvalues query contexts in
    // docs/evidence/tlosan_stage_d/run_20260924/*.trace.
    #[test]
    fn local_code1_context_lengths_match_ncbi_function_inputs() {
        let params = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        assert_eq!(
            local_subject_effective_lengths(&[(160, true)], 362, &[params]),
            vec![LocalContextLength {
                length_adjustment: 15,
                eff_searchsp: 15225,
                subject_stat_length: 120,
            }]
        );
        assert_eq!(
            local_subject_effective_lengths(
                &[(120, true), (70, true), (120, false)],
                6377,
                &[params; 3],
            ),
            vec![
                LocalContextLength {
                    length_adjustment: 33,
                    eff_searchsp: 182004,
                    subject_stat_length: 2125,
                },
                LocalContextLength {
                    length_adjustment: 28,
                    eff_searchsp: 88074,
                    subject_stat_length: 2125,
                },
                LocalContextLength {
                    length_adjustment: 0,
                    eff_searchsp: 0,
                    subject_stat_length: 2125,
                },
            ]
        );
    }
}
