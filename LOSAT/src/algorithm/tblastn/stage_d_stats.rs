//! TBLASTN local-subject statistical inputs before preliminary search.

use crate::algorithm::tblastx::ncbi_cutoffs::{
    cutoff_score_from_evalue, cutoff_score_sum_stats, gap_trigger_raw_score, x_drop_raw_score,
};
use crate::stats::length_adjustment::compute_length_adjustment_ncbi;
use crate::stats::spouge::{blast_spouge_etos, BlastGumbelBlk};
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

// NCBI c++/include/algo/blast/core/blast_parameters.h:69-71;
// c++/src/algo/blast/core/blast_parameters.c:604-625,793-805:
// #define BLAST_GAP_DECAY_RATE_GAPPED 0.1
// #define BLAST_GAP_SIZE 40
// #define BLAST_OVERLAP_SIZE 9
// params->longest_intron = (DEFAULT_LONGEST_INTRON - 2) / 3;
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct LocalLinkParameters {
    pub gap_decay_rate: f64,
    pub gap_size: i32,
    pub overlap_size: i32,
    pub longest_intron: i32,
    pub cutoff_small_gap: i32,
}

// NCBI c++/src/algo/blast/core/blast_parameters.c:902-999,302-383:
// params->cutoffs[context].cutoff_score = new_cutoff;
// params->cutoffs[context].cutoff_score_max = new_cutoff;
// curr_cutoffs->cutoff_score = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct LocalContextCutoffs {
    pub hit_cutoff: i32,
    pub hit_cutoff_max: i32,
    pub word_cutoff: i32,
    pub word_xdrop: i32,
}

// NCBI c++/src/algo/blast/core/blast_setup.c:1011-1024:
// BLAST_CalcEffLengths(...); BlastHitSavingParametersUpdate(...);
// BlastInitialWordParametersUpdate(...); BlastLinkHSPParametersUpdate(...);
#[derive(Clone, Debug, PartialEq)]
pub(super) struct LocalSubjectParameters {
    pub lengths: Vec<LocalContextLength>,
    pub cutoffs: Vec<LocalContextCutoffs>,
    pub link: Option<LocalLinkParameters>,
    pub hit_cutoff_min: i32,
    pub word_cutoff_min: i32,
    pub prelim_evalue: f64,
}

// NCBI c++/src/algo/blast/core/blast_parameters.c:774-815,831-856;
// c++/src/algo/blast/core/blast_setup.c:1016-1024:
// if (params->do_sum_stats) BlastLinkHSPParametersNew(...);
// BlastHitSavingParametersUpdate(..., subject_length, 0, hit_params);
// The latter call passes zero compositionBasedStats even when the search uses
// composition mode 2; the score block's Gumbel state is still consulted.
#[derive(Clone, Copy)]
pub(super) struct LocalParameterOptions<'a> {
    pub expect_value: f64,
    pub do_sum_stats: bool,
    pub max_intron_length: i32,
    pub gap_trigger_bits: f64,
    pub word_xdrop_bits: f64,
    pub scale_factor: f64,
    pub gumbel: Option<&'a BlastGumbelBlk>,
}

// NCBI c++/src/algo/blast/core/blast_setup.c:964-985,1001-1024:
// min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
// if (Blast_SubjectIsTranslated(program_number)) min_subject_length/=3;
// BlastHitSavingParametersNew(..., min_subject_length,
//                              (*ext_params)->options->compositionBasedStats, ...);
// BlastHitSavingParametersUpdate(..., subject_length, 0, hit_params);
#[derive(Clone, Copy)]
pub(super) enum LocalParameterCall {
    Initial {
        min_subject_length: i32,
        composition_based_stats: i32,
    },
    OneSubjectUpdate,
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

// NCBI c++/src/algo/blast/core/blast_setup.c:1011-1024:
// BlastHitSavingParametersNew(..., min_subject_length,
//                              (*ext_params)->options->compositionBasedStats, ...);
// eff_len_params->real_db_length = subject_length;
// BLAST_CalcEffLengths(...);
// BlastHitSavingParametersUpdate(..., subject_length, 0, hit_params);
// BlastInitialWordParametersUpdate(..., subject_length, word_params);
// BlastLinkHSPParametersUpdate(word_params, hit_params, TRUE);
#[allow(dead_code)] // The public TBLASTN path remains gated until C/D/E pass.
pub(super) fn local_parameters_for_call(
    query_contexts: &[(usize, bool)],
    subject_nt_length: usize,
    gapped_params: &[KarlinParams],
    ungapped_params: &[KarlinParams],
    options: LocalParameterOptions<'_>,
    call: LocalParameterCall,
) -> LocalSubjectParameters {
    assert!(!query_contexts.is_empty());
    assert_eq!(query_contexts.len(), gapped_params.len());
    assert_eq!(query_contexts.len(), ungapped_params.len());
    let lengths = local_subject_effective_lengths(query_contexts, subject_nt_length, gapped_params);
    // NCBI c++/src/algo/blast/core/blast_setup.c:964-985,1011-1024:
    // initial creation passes min_subject_length and compositionBasedStats;
    // the conditional update passes subject_length and literal zero.
    let (hit_subject_length, composition_based_stats) = match call {
        LocalParameterCall::Initial {
            min_subject_length,
            composition_based_stats,
        } => (min_subject_length, composition_based_stats),
        LocalParameterCall::OneSubjectUpdate => (subject_nt_length as i32, 0),
    };
    let mut link = if options.do_sum_stats {
        // NCBI c++/src/algo/blast/core/blast_parameters.c:774-815;
        // c++/include/algo/blast/core/blast_def.h:77-78:
        // if (options->longest_intron == 0)
        //     longest_intron = (DEFAULT_LONGEST_INTRON - 2) / 3;
        // else if ((options->longest_intron - 2)/3 <= 0) linking is disabled.
        let intron = if options.max_intron_length == 0 {
            (122 - 2) / 3
        } else {
            (options.max_intron_length - 2) / 3
        };
        (intron > 0).then_some(LocalLinkParameters {
            gap_decay_rate: 0.1,
            gap_size: 40,
            overlap_size: 9,
            longest_intron: intron,
            cutoff_small_gap: 0,
        })
    } else {
        None
    };
    let mut cutoffs = Vec::with_capacity(query_contexts.len());
    let mut hit_cutoff_min = i32::MAX;
    let mut word_cutoff_min = i32::MAX;
    let mut prelim_evalue = options.expect_value;

    // NCBI c++/src/algo/blast/core/blast_parameters.c:902-999:
    // if (sbp->gbp && sbp->gbp->filled)
    //     new_cutoff = BLAST_SpougeEtoS(evalue, kbp, sbp->gbp,
    //                                   query_length, avg_subject_length);
    // else BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
    // NCBI blast_parameters.c:931-946:
    // cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
    // params->prelim_evalue = cbs_stretch * evalue;
    let mut hit = Vec::with_capacity(query_contexts.len());
    for (context, (&(query_length, valid), params)) in
        query_contexts.iter().zip(gapped_params).enumerate()
    {
        if !valid || query_length == 0 {
            hit.push((i32::MAX, 0));
            continue;
        }
        let max_cutoff = if let Some(gumbel) = options.gumbel {
            // NCBI c++/src/algo/blast/core/blast_parameters.c:931-945:
            // cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
            // params->prelim_evalue = cbs_stretch * evalue;
            let cbs_stretch = if composition_based_stats > 1 {
                5.0
            } else {
                1.0
            };
            prelim_evalue = cbs_stretch * options.expect_value;
            blast_spouge_etos(
                prelim_evalue,
                params,
                gumbel,
                query_length as i32,
                hit_subject_length,
            )
        } else {
            cutoff_score_from_evalue(options.expect_value, lengths[context].eff_searchsp, params)
                .max(1)
        };
        hit.push((max_cutoff, max_cutoff));
    }

    if let Some(link_params) = link {
        // NCBI c++/src/algo/blast/core/blast_parameters.c:950-976:
        // concat_qlen = last.query_offset + last.query_length;
        // avg_qlen = concat_qlen / (last_context + 1);
        // searchsp = MIN(avg_qlen, avg_subject_length) * avg_subject_length;
        // BLAST_Cutoffs(&new_cutoff, &evalue_hsp, kbp, searchsp, TRUE, gap_decay_rate);
        let concat_qlen: usize = query_contexts
            .iter()
            .map(|(length, _)| length + 1)
            .sum::<usize>()
            - 1;
        let avg_qlen = (concat_qlen / query_contexts.len()) as i32;
        for (context, (&(_, valid), params)) in query_contexts.iter().zip(gapped_params).enumerate()
        {
            if valid {
                let cutoff = cutoff_score_sum_stats(
                    avg_qlen,
                    hit_subject_length,
                    link_params.gap_decay_rate,
                    params,
                );
                hit[context].0 = hit[context].0.min(cutoff);
            }
        }
    }

    // NCBI c++/src/algo/blast/core/blast_parameters.c:983-999,302-383:
    // hit cutoff and maximum are scaled before the word cutoff is computed;
    // gap_trigger uses kbp_std and the word cutoff is capped by hit cutoff_max.
    for (context, (&(query_length, valid), ungapped)) in
        query_contexts.iter().zip(ungapped_params).enumerate()
    {
        if !valid || query_length == 0 {
            cutoffs.push(LocalContextCutoffs {
                hit_cutoff: i32::MAX,
                hit_cutoff_max: 0,
                word_cutoff: i32::MAX,
                word_xdrop: 0,
            });
            continue;
        }
        let hit_cutoff = hit[context].0 * options.scale_factor as i32;
        let hit_cutoff_max = hit[context].1 * options.scale_factor as i32;
        hit_cutoff_min = hit_cutoff_min.min(hit_cutoff);
        let gap_trigger = gap_trigger_raw_score(options.gap_trigger_bits, ungapped);
        let word_cutoff = ((gap_trigger as f64 * options.scale_factor) as i32).min(hit_cutoff_max);
        let xdrop_init = x_drop_raw_score(options.word_xdrop_bits, ungapped, options.scale_factor);
        let word_xdrop = if xdrop_init == 0 {
            word_cutoff
        } else {
            xdrop_init
        };
        word_cutoff_min = word_cutoff_min.min(word_cutoff);
        cutoffs.push(LocalContextCutoffs {
            hit_cutoff,
            hit_cutoff_max,
            word_cutoff,
            word_xdrop,
        });
    }
    if let Some(link_params) = &mut link {
        // NCBI c++/src/algo/blast/core/blast_parameters.c:625-649:
        // hit_params->link_hsp_params->cutoff_small_gap =
        //     word_params->cutoff_score_min;
        link_params.cutoff_small_gap = word_cutoff_min;
    }
    LocalSubjectParameters {
        lengths,
        cutoffs,
        link,
        hit_cutoff_min,
        word_cutoff_min,
        prelim_evalue,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{ProteinScoringSpec, ScoringMatrix};
    use crate::stats::spouge::lookup_protein_gumbel_params;
    use crate::stats::tables::{lookup_protein_params_gapped, lookup_protein_params_ungapped};

    // Pinned NCBI blast_setup.c:729-735,770-847 and the actual
    // BLAST_LinkHsps/Blast_HSPListGetEvalues query contexts in
    // docs/evidence/tlosan_stage_d/run_20260924/*.trace.
    // NCBI c++/src/algo/blast/core/blast_setup.c:1011-1024;
    // c++/src/algo/blast/core/blast_parameters.c:902-999,302-383,625-649:
    // update effective lengths, hit cutoffs, word cutoffs, then link cutoff.
    // The saved C fixture records both context cutoffs and word x-drop on every call.
    #[test]
    fn computed_code1_control_parameters_match_stage_c_inputs() {
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let subject_length = 15_000_962usize;
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            (subject_length / 3) as i64,
        )
        .unwrap();
        let result = local_parameters_for_call(
            &[(120, true), (70, true), (120, false)],
            subject_length,
            &[gapped; 3],
            &[ungapped; 3],
            LocalParameterOptions {
                expect_value: 10_000.0,
                do_sum_stats: false,
                max_intron_length: 0,
                gap_trigger_bits: 22.0,
                word_xdrop_bits: 7.0,
                scale_factor: 1.0,
                gumbel: Some(&gumbel),
            },
            LocalParameterCall::Initial {
                min_subject_length: (subject_length / 3) as i32,
                composition_based_stats: 0,
            },
        );
        assert_eq!(result.link, None);
        assert_eq!(
            result
                .cutoffs
                .iter()
                .map(|c| c.hit_cutoff)
                .collect::<Vec<_>>(),
            [28, 25, i32::MAX]
        );
        assert_eq!(
            result
                .cutoffs
                .iter()
                .map(|c| c.hit_cutoff_max)
                .collect::<Vec<_>>(),
            [28, 25, 0]
        );
        assert_eq!(
            result
                .cutoffs
                .iter()
                .map(|c| c.word_cutoff)
                .collect::<Vec<_>>(),
            [28, 25, i32::MAX]
        );
        assert_eq!(
            result
                .cutoffs
                .iter()
                .map(|c| c.word_xdrop)
                .collect::<Vec<_>>(),
            [16, 16, 0]
        );
    }

    // NCBI c++/src/algo/blast/core/link_hsps.c:1788-1802;
    // c++/src/algo/blast/core/blast_hits.c:1811-1926:
    // Blast_HSPListGetEvalues receives translated subject_length/3 and
    // divides each Spouge E-value by BLAST_GapDecayDivisor(rate, 1).
    #[test]
    fn preliminary_individual_evalues_match_ncbi_double_bits() {
        let trace = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/run_20260924/multi_query_20260924_default.trace"
        ));
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            2_125,
        )
        .unwrap();
        let mut count = 0;
        for line in trace
            .lines()
            .filter(|line| line.starts_with("D_HSP\t1\tevalue_after\t"))
        {
            let f: Vec<_> = line.split('\t').collect();
            let context: usize = f[4].parse().unwrap();
            let score: i32 = f[10].parse().unwrap();
            let expected: f64 = f[13].parse().unwrap();
            let actual = crate::stats::spouge::blast_spouge_stoe(
                score,
                &gapped,
                &gumbel,
                [120, 70][context],
                2_125,
            ) / 0.9;
            assert_eq!(actual.to_bits(), expected.to_bits(), "HSP {count}: {line}");
            count += 1;
        }
        assert_eq!(count, 21);
    }

    // NCBI c++/src/algo/blast/core/link_hsps.c:1117-1155;
    // c++/src/algo/blast/core/blast_stat.c:4491-4532:
    // sum_score = new_hsp->sum_score + head_hsp->sum_score;
    // BLAST_UnevenGapSumE(50, 50, 2, sum_score, 105, 135,
    //                     14175, BLAST_GapDecayDivisor(0.1, 2));
    #[test]
    fn positive_uneven_gap_pair_matches_ncbi_double_bits() {
        let params = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let first = params.lambda * 331.0 - params.k.ln();
        let second = params.lambda * 329.0 - params.k.ln();
        let divisor = crate::stats::sum_statistics::gap_decay_divisor(0.1, 2);
        let actual = crate::stats::sum_statistics::uneven_gap_sum_e(
            50,
            50,
            2,
            second + first,
            105,
            135,
            14_175,
            divisor,
        );
        let expected = 3.2037962798257965e-69_f64;
        assert_eq!(actual.to_bits(), expected.to_bits(), "actual={actual:.17e}");
    }

    // NCBI c++/src/algo/blast/core/blast_parameters.c:774-815,950-976,625-649:
    // default translated/gapped searches keep uneven-gap linking and set its
    // small-gap cutoff from word_params->cutoff_score_min after the update.
    #[test]
    fn computed_code1_default_link_parameters_match_stage_d_call() {
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let subject_length = 6_377usize;
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            (subject_length / 3) as i64,
        )
        .unwrap();
        let result = local_parameters_for_call(
            &[(120, true), (70, true), (120, false)],
            subject_length,
            &[gapped; 3],
            &[ungapped; 3],
            LocalParameterOptions {
                expect_value: 10.0,
                do_sum_stats: true,
                max_intron_length: 0,
                gap_trigger_bits: 22.0,
                word_xdrop_bits: 7.0,
                scale_factor: 1.0,
                gumbel: Some(&gumbel),
            },
            LocalParameterCall::Initial {
                min_subject_length: (subject_length / 3) as i32,
                composition_based_stats: 2,
            },
        );
        assert_eq!(result.link.unwrap().longest_intron, 40);
        assert_eq!(result.link.unwrap().cutoff_small_gap, 17, "{result:?}");
    }

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
