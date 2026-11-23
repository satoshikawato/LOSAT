use super::*;
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::stats::spouge::lookup_protein_gumbel_params;
use crate::stats::tables::{lookup_protein_params_gapped, lookup_protein_params_ungapped};
use crate::utils::glibc_fma_libm::GlibcLibm;

// Pinned NCBI blast_setup.c:729-735,770-847 and the actual
// BLAST_LinkHsps/Blast_HSPListGetEvalues query contexts in
// docs/evidence/tlosan_stage_d/run_20260924/*.trace.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1011-1024;
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:729-847:
    let result = local_parameters_for_call(
        &[(120, true), (70, true), (120, false)],
        subject_length as i64,
        1,
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
            eff_searchsp_override: None,
            ungapped: None,
        },
        LocalParameterCall::Initial {
            min_subject_length: (subject_length / 3) as i32,
            composition_based_stats: 0,
        },
    )
    .unwrap();
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:766-768 (the effective lengths from `kbp_gap_std`), c++/src/algo/blast/core/blast_parameters.c:
// 204-221 (the word X-drop from `kbp`), 338-345 (the gap trigger from `kbp_std`) and 902-999
// (the cutoffs from `kbp_gap`): each block of a PSSM query is read where NCBI reads it (tblastn
// `-in_pssm`, line PS item PS7).
#[test]
fn each_karlin_block_is_read_where_ncbi_reads_it() {
    let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let options = LocalParameterOptions {
        expect_value: 10.0,
        do_sum_stats: false,
        max_intron_length: 0,
        gap_trigger_bits: 22.0,
        word_xdrop_bits: 7.0,
        scale_factor: 1.0,
        gumbel: None,
        eff_searchsp_override: None,
        ungapped: None,
    };
    let call = LocalParameterCall::Initial {
        min_subject_length: 1000,
        composition_based_stats: 0,
    };
    let run = |blocks: LocalKarlinBlocks<'_>| {
        local_parameters_for_blocks(&[(300, true)], 3_000_000, 10, blocks, options, call).unwrap()
    };
    let base = local_parameters_for_call(
        &[(300, true)],
        3_000_000,
        10,
        &[gapped],
        &[ungapped],
        options,
        call,
    )
    .unwrap();
    let same = run(LocalKarlinBlocks {
        gapped: &[gapped],
        lengths: &[gapped],
        ungapped: &[ungapped],
        x_dropoff: &[ungapped],
    });
    assert_eq!(base, same);
    // Another gapped block for the cutoffs only: the lengths stay the matrix's.
    let psi_gapped = KarlinParams {
        lambda: gapped.lambda * 0.8,
        k: gapped.k * 2.0,
        ..gapped
    };
    let cutoffs = run(LocalKarlinBlocks {
        gapped: &[psi_gapped],
        lengths: &[gapped],
        ungapped: &[ungapped],
        x_dropoff: &[ungapped],
    });
    assert_eq!(cutoffs.lengths, base.lengths);
    assert_ne!(cutoffs.cutoffs[0].hit_cutoff, base.cutoffs[0].hit_cutoff);
    let lengths = run(LocalKarlinBlocks {
        gapped: &[gapped],
        lengths: &[psi_gapped],
        ungapped: &[ungapped],
        x_dropoff: &[ungapped],
    });
    assert_ne!(lengths.lengths, base.lengths);
    // Another ungapped block for the X-drop only: the gap trigger stays `kbp_std`'s.
    let psi_ungapped = KarlinParams {
        lambda: ungapped.lambda * 0.5,
        ..ungapped
    };
    let xdrop = run(LocalKarlinBlocks {
        gapped: &[gapped],
        lengths: &[gapped],
        ungapped: &[ungapped],
        x_dropoff: &[psi_ungapped],
    });
    assert_eq!(xdrop.cutoffs[0].word_cutoff, base.cutoffs[0].word_cutoff);
    assert_ne!(xdrop.cutoffs[0].word_xdrop, base.cutoffs[0].word_xdrop);
    let trigger = run(LocalKarlinBlocks {
        gapped: &[gapped],
        lengths: &[gapped],
        ungapped: &[psi_ungapped],
        x_dropoff: &[ungapped],
    });
    assert_eq!(trigger.cutoffs[0].word_xdrop, base.cutoffs[0].word_xdrop);
    assert_ne!(trigger.cutoffs[0].word_cutoff, base.cutoffs[0].word_cutoff);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1788-1802;
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1117-1155;
// c++/src/algo/blast/core/blast_stat.c:4491-4532:
#[test]
fn positive_uneven_gap_pair_matches_ncbi_double_bits() {
    let params = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
    let first = params.lambda * 331.0 - params.k.glibc_ln();
    let second = params.lambda * 329.0 - params.k.glibc_ln();
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:774-815,950-976,625-649:
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:729-847:
    let result = local_parameters_for_call(
        &[(120, true), (70, true), (120, false)],
        subject_length as i64,
        1,
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
            eff_searchsp_override: None,
            ungapped: None,
        },
        LocalParameterCall::Initial {
            min_subject_length: (subject_length / 3) as i32,
            composition_based_stats: 2,
        },
    )
    .unwrap();
    assert_eq!(result.link.unwrap().longest_intron, 40);
    assert_eq!(result.link.unwrap().cutoff_small_gap, 17, "{result:?}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:729-847:
// The pinned natural two-subject trace records D_CONTEXT 0, 18, 22848.
#[test]
fn two_local_subjects_use_both_sequences_in_effective_search_space() {
    let params = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
    assert_eq!(
        local_subject_effective_lengths(&[(120, true)], 781, 2, &[params]).unwrap(),
        vec![LocalContextLength {
            length_adjustment: 18,
            eff_searchsp: 22_848,
            subject_stat_length: 260,
        }]
    );
}

#[test]
fn local_code1_context_lengths_match_ncbi_function_inputs() {
    let params = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
    assert_eq!(
        local_subject_effective_lengths(&[(160, true)], 362, 1, &[params]).unwrap(),
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
            1,
            &[params; 3],
        )
        .unwrap(),
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

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:217
// A database of 2^32 letters or more (an alias file's `LENGTH`: 68,719,476,736 letters in 16
// sequences, row db.cnt.avg2p32; 4,294,967,296,000,000,000 with `-searchsp 1000000`) reaches
// the effective lengths whole, as NCBI's `Int8`, on every target: its third is the statistics'
// length (`stats::setup::calc_eff_lengths`, which these lengths are).
#[test]
fn effective_lengths_read_a_database_of_2_to_the_32_letters_or_more() {
    let params = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
    for (total, search_space) in [
        (68_719_476_736i64, 0),
        (4_294_967_296_000_000_000, 1_000_000),
        (4_294_967_295_999_999_999, 1_000_000),
    ] {
        assert!(total > i64::from(u32::MAX));
        let expected = crate::stats::setup::calc_eff_lengths(
            120,
            total / 3,
            16,
            &params,
            search_space,
            "TBLASTN",
        )
        .unwrap();
        let spaces = [search_space];
        assert_eq!(
            local_subject_effective_lengths_with_search_spaces(
                &[(120, true)],
                total,
                16,
                &[params],
                (search_space != 0).then_some(&spaces[..]),
            )
            .unwrap(),
            vec![LocalContextLength {
                length_adjustment: i64::from(expected.length_adjustment),
                eff_searchsp: expected.eff_searchsp,
                subject_stat_length: total / 3,
            }],
            "{total}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4059-4061
// The score of an e-value E in a search space N is ceil(ln(K * N / E) / Lambda).
// The ungapped search's word cutoff: e-value 1 divided by the gap decay divisor (1 - 0.5) in
// the space of the shorter of subject and query times the subject (BLOSUM62's ungapped block,
// Lambda 0.3176, K 0.134): ceil(ln(0.134 * 30000 / 0.5) / 0.3176) = 29, with a query longer
// than the subject
// ceil(ln(0.134 * 90000 / 0.5) / 0.3176) = 32, and without the divisor (a decay outside
// (0, 1)) ceil(ln(0.134 * 30000) / 0.3176) = 27.
#[test]
fn ungapped_word_cutoff_takes_the_shorter_length_and_the_decay() {
    let kbp = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    assert_eq!(ungapped_word_cutoff(100, 300, 0.5, &kbp), 29);
    assert_eq!(ungapped_word_cutoff(1000, 300, 0.5, &kbp), 32);
    assert_eq!(ungapped_word_cutoff(100, 300, 0.0, &kbp), 27);
    // A space below one letter's score keeps NCBI's cutoff of 1.
    assert_eq!(ungapped_word_cutoff(1, 1, 0.5, &kbp), 1);
}
