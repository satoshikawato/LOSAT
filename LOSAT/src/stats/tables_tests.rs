use super::*;

fn spec(reward: i32, penalty: i32, gap_open: i32, gap_extend: i32) -> NuclScoringSpec {
    NuclScoringSpec {
        reward,
        penalty,
        gap_open,
        gap_extend,
    }
}

const UNGAPPED: KarlinParams = KarlinParams {
    lambda: 1.1,
    k: 0.3,
    h: 0.5,
    alpha: 0.0,
    beta: 0.0,
};

fn values(spec: &NuclScoringSpec) -> Result<(KarlinParams, bool), String> {
    let values = NuclValues::new(spec.reward, spec.penalty)?;
    Ok((values.gapped(spec, &UNGAPPED)?, values.round_down))
}

#[test]
fn task_defaults_use_the_table_rows() {
    // blastn_values_1_2[0] (linear) and blastn_values_2_3 row { 5, 2, ... }.
    let (megablast, round_down) = values(&spec(1, -2, 0, 0)).unwrap();
    assert_eq!(
        (
            megablast.lambda,
            megablast.k,
            megablast.h,
            megablast.alpha,
            megablast.beta
        ),
        (1.28, 0.46, 0.85, 1.5, -2.0)
    );
    assert!(!round_down);
    let (blastn, round_down) = values(&spec(2, -3, 5, 2)).unwrap();
    assert_eq!(
        (blastn.lambda, blastn.k, blastn.h, blastn.alpha, blastn.beta),
        (0.625, 0.41, 0.78, 0.8, -2.0)
    );
    assert!(round_down);
}

#[test]
fn scores_with_a_common_divisor_use_the_divided_table() {
    // 2/-4 is 1/-2: the gap costs double, Lambda and alpha halve.
    let (params, round_down) = values(&spec(2, -4, 4, 4)).unwrap();
    assert_eq!(
        (params.lambda, params.k, params.h, params.alpha, params.beta),
        (1.33 / 2.0, 0.62, 1.1, 1.2 / 2.0, 0.0)
    );
    assert!(!round_down);
    // 4/-6 is 2/-3, whose scores are rounded down.
    assert!(values(&spec(4, -6, 10, 4)).unwrap().1);
    assert!(values(&spec(2, -4, 3, 3))
        .unwrap_err()
        .contains("\n4 and 4 are supported"));
}

#[test]
fn gap_costs_beyond_the_table_use_the_ungapped_block() {
    let (params, _) = values(&spec(1, -2, 5, 2)).unwrap();
    assert_eq!(
        (params.lambda, params.k, params.h, params.alpha, params.beta),
        (1.1, 0.3, 0.5, 1.1 / 0.5, 0.0)
    );
    assert_eq!(values(&spec(2, -3, 7, 4)).unwrap().0.beta, -2.0);
}

#[test]
fn unsupported_scores_have_ncbi_messages() {
    assert_eq!(
        values(&spec(2, -5, 5, 2)).unwrap_err(),
        "Gap existence and extension values 5 and 2 are not supported for substitution scores 2 and -5\n\
             2 and 4 are supported existence and extension values\n\
             0 and 4 are supported existence and extension values\n\
             4 and 2 are supported existence and extension values\n\
             2 and 2 are supported existence and extension values\n\
             4 and 4 are supported existence and extension values\n\
             Any values more stringent than 4 and 4 are supported\n"
    );
    assert_eq!(
        values(&spec(1, -6, 0, 0)).unwrap_err(),
        "Substitution scores 1 and -6 are not supported"
    );
    // The message names the divided scores; 5/-4 has no linear row.
    assert_eq!(
        values(&spec(2, -12, 0, 0)).unwrap_err(),
        "Substitution scores 1 and -6 are not supported"
    );
    assert!(values(&spec(5, -4, 0, 0))
        .unwrap_err()
        .starts_with("Gap existence and extension values 0 and 0 are not supported"));
}

// The protein rows that this file carried as literals before the tables were read from
// protein_tables.rs (dev fb8b10ee6): (gap open, gap extend, Lambda, K, H, alpha, beta),
// i32::MAX marking the ungapped row.
const LITERAL_BLOSUM45: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.2291, 0.0924, 0.2514, 0.9113, -5.7),
    ParamEntry::new(13, 3, 0.207, 0.049, 0.14, 1.5, -22.0),
    ParamEntry::new(12, 3, 0.199, 0.039, 0.11, 1.8, -34.0),
    ParamEntry::new(11, 3, 0.190, 0.031, 0.095, 2.0, -38.0),
    ParamEntry::new(10, 3, 0.179, 0.023, 0.075, 2.4, -51.0),
    ParamEntry::new(16, 2, 0.210, 0.051, 0.14, 1.5, -24.0),
    ParamEntry::new(15, 2, 0.203, 0.041, 0.12, 1.7, -31.0),
    ParamEntry::new(14, 2, 0.195, 0.032, 0.10, 1.9, -36.0),
    ParamEntry::new(13, 2, 0.185, 0.024, 0.084, 2.2, -45.0),
    ParamEntry::new(12, 2, 0.171, 0.016, 0.061, 2.8, -65.0),
    ParamEntry::new(19, 1, 0.205, 0.040, 0.11, 1.9, -43.0),
    ParamEntry::new(18, 1, 0.198, 0.032, 0.10, 2.0, -43.0),
    ParamEntry::new(17, 1, 0.189, 0.024, 0.079, 2.4, -57.0),
    ParamEntry::new(16, 1, 0.176, 0.016, 0.063, 2.8, -67.0),
];

const LITERAL_BLOSUM50: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.2318, 0.112, 0.3362, 0.6895, -4.0),
    ParamEntry::new(13, 3, 0.212, 0.063, 0.19, 1.1, -16.0),
    ParamEntry::new(12, 3, 0.206, 0.055, 0.17, 1.2, -18.0),
    ParamEntry::new(11, 3, 0.197, 0.042, 0.14, 1.4, -25.0),
    ParamEntry::new(10, 3, 0.186, 0.031, 0.11, 1.7, -34.0),
    ParamEntry::new(9, 3, 0.172, 0.022, 0.082, 2.1, -48.0),
    ParamEntry::new(16, 2, 0.215, 0.066, 0.20, 1.05, -15.0),
    ParamEntry::new(15, 2, 0.210, 0.058, 0.17, 1.2, -20.0),
    ParamEntry::new(14, 2, 0.202, 0.045, 0.14, 1.4, -27.0),
    ParamEntry::new(13, 2, 0.193, 0.035, 0.12, 1.6, -32.0),
    ParamEntry::new(12, 2, 0.181, 0.025, 0.095, 1.9, -41.0),
    ParamEntry::new(19, 1, 0.212, 0.057, 0.18, 1.2, -21.0),
    ParamEntry::new(18, 1, 0.207, 0.050, 0.15, 1.4, -28.0),
    ParamEntry::new(17, 1, 0.198, 0.037, 0.12, 1.6, -33.0),
    ParamEntry::new(16, 1, 0.186, 0.025, 0.10, 1.9, -42.0),
    ParamEntry::new(15, 1, 0.171, 0.015, 0.063, 2.7, -76.0),
];

const LITERAL_BLOSUM62: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2),
    ParamEntry::new(11, 2, 0.297, 0.082, 0.27, 1.1, -10.0),
    ParamEntry::new(10, 2, 0.291, 0.075, 0.23, 1.3, -15.0),
    ParamEntry::new(9, 2, 0.279, 0.058, 0.19, 1.5, -19.0),
    ParamEntry::new(8, 2, 0.264, 0.045, 0.15, 1.8, -26.0),
    ParamEntry::new(7, 2, 0.239, 0.027, 0.10, 2.5, -46.0),
    ParamEntry::new(6, 2, 0.201, 0.012, 0.061, 3.3, -58.0),
    ParamEntry::new(13, 1, 0.292, 0.071, 0.23, 1.2, -11.0),
    ParamEntry::new(12, 1, 0.283, 0.059, 0.19, 1.5, -19.0),
    ParamEntry::new(11, 1, 0.267, 0.041, 0.14, 1.9, -30.0),
    ParamEntry::new(10, 1, 0.243, 0.024, 0.10, 2.5, -44.0),
    ParamEntry::new(9, 1, 0.206, 0.010, 0.052, 4.0, -87.0),
];

const LITERAL_BLOSUM80: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3430, 0.177, 0.6568, 0.5222, -1.6),
    ParamEntry::new(25, 2, 0.342, 0.17, 0.66, 0.52, -1.6),
    ParamEntry::new(13, 2, 0.336, 0.15, 0.57, 0.59, -3.0),
    ParamEntry::new(9, 2, 0.319, 0.11, 0.42, 0.76, -6.0),
    ParamEntry::new(8, 2, 0.308, 0.090, 0.35, 0.89, -9.0),
    ParamEntry::new(7, 2, 0.293, 0.070, 0.27, 1.1, -14.0),
    ParamEntry::new(6, 2, 0.268, 0.045, 0.19, 1.4, -19.0),
    ParamEntry::new(11, 1, 0.314, 0.095, 0.35, 0.90, -9.0),
    ParamEntry::new(10, 1, 0.299, 0.071, 0.27, 1.1, -14.0),
    ParamEntry::new(9, 1, 0.279, 0.048, 0.20, 1.4, -19.0),
];

const LITERAL_BLOSUM90: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3346, 0.190, 0.7547, 0.4434, -1.4),
    ParamEntry::new(9, 2, 0.310, 0.12, 0.46, 0.67, -6.0),
    ParamEntry::new(8, 2, 0.300, 0.099, 0.39, 0.76, -7.0),
    ParamEntry::new(7, 2, 0.283, 0.072, 0.30, 0.93, -11.0),
    ParamEntry::new(6, 2, 0.259, 0.048, 0.22, 1.2, -16.0),
    ParamEntry::new(11, 1, 0.302, 0.093, 0.39, 0.78, -8.0),
    ParamEntry::new(10, 1, 0.290, 0.075, 0.28, 1.04, -15.0),
    ParamEntry::new(9, 1, 0.265, 0.044, 0.20, 1.3, -19.0),
];

const LITERAL_PAM30: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3400, 0.283, 1.754, 0.1938, -0.3),
    ParamEntry::new(7, 2, 0.305, 0.15, 0.87, 0.35, -3.0),
    ParamEntry::new(6, 2, 0.287, 0.11, 0.68, 0.42, -4.0),
    ParamEntry::new(5, 2, 0.264, 0.079, 0.45, 0.59, -7.0),
    ParamEntry::new(10, 1, 0.309, 0.15, 0.88, 0.34, -3.0),
    ParamEntry::new(9, 1, 0.294, 0.11, 0.61, 0.48, -6.0),
    ParamEntry::new(8, 1, 0.270, 0.072, 0.40, 0.68, -10.0),
];

const LITERAL_PAM70: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.3345, 0.229, 1.029, 0.3250, -0.9),
    ParamEntry::new(8, 2, 0.301, 0.12, 0.54, 0.56, -5.0),
    ParamEntry::new(7, 2, 0.286, 0.093, 0.43, 0.67, -7.0),
    ParamEntry::new(6, 2, 0.264, 0.064, 0.29, 0.90, -12.0),
    ParamEntry::new(11, 1, 0.305, 0.12, 0.52, 0.59, -6.0),
    ParamEntry::new(10, 1, 0.291, 0.091, 0.41, 0.71, -9.0),
    ParamEntry::new(9, 1, 0.270, 0.060, 0.28, 0.97, -14.0),
];

const LITERAL_PAM250: &[ParamEntry] = &[
    ParamEntry::new(i32::MAX, i32::MAX, 0.2252, 0.0868, 0.2223, 0.98, -5.0),
    ParamEntry::new(15, 3, 0.205, 0.049, 0.13, 1.6, -23.0),
    ParamEntry::new(14, 3, 0.200, 0.043, 0.12, 1.7, -26.0),
    ParamEntry::new(13, 3, 0.194, 0.036, 0.10, 1.9, -31.0),
    ParamEntry::new(12, 3, 0.186, 0.029, 0.085, 2.2, -41.0),
    ParamEntry::new(11, 3, 0.174, 0.020, 0.070, 2.5, -48.0),
    ParamEntry::new(17, 2, 0.204, 0.047, 0.12, 1.7, -28.0),
    ParamEntry::new(16, 2, 0.198, 0.038, 0.11, 1.8, -29.0),
    ParamEntry::new(15, 2, 0.191, 0.031, 0.087, 2.2, -44.0),
    ParamEntry::new(14, 2, 0.182, 0.024, 0.073, 2.5, -53.0),
    ParamEntry::new(13, 2, 0.171, 0.017, 0.059, 2.9, -64.0),
    ParamEntry::new(21, 1, 0.205, 0.045, 0.11, 1.8, -34.0),
    ParamEntry::new(20, 1, 0.199, 0.037, 0.10, 1.9, -35.0),
    ParamEntry::new(19, 1, 0.192, 0.029, 0.083, 2.3, -52.0),
    ParamEntry::new(18, 1, 0.183, 0.021, 0.070, 2.6, -60.0),
    ParamEntry::new(17, 1, 0.171, 0.014, 0.052, 3.3, -86.0),
];

#[test]
fn protein_table_rows_equal_the_former_literals_bit_for_bit() {
    let bits = |p: KarlinParams| [p.lambda, p.k, p.h, p.alpha, p.beta].map(f64::to_bits);
    let literal_tables: [(ScoringMatrix, &[ParamEntry]); 8] = [
        (ScoringMatrix::Blosum45, LITERAL_BLOSUM45),
        (ScoringMatrix::Blosum50, LITERAL_BLOSUM50),
        (ScoringMatrix::Blosum62, LITERAL_BLOSUM62),
        (ScoringMatrix::Blosum80, LITERAL_BLOSUM80),
        (ScoringMatrix::Blosum90, LITERAL_BLOSUM90),
        (ScoringMatrix::Pam30, LITERAL_PAM30),
        (ScoringMatrix::Pam70, LITERAL_PAM70),
        (ScoringMatrix::Pam250, LITERAL_PAM250),
    ];
    let mut rows = 0;
    for (matrix, entries) in literal_tables {
        for entry in entries {
            let mut expected = entry.to_karlin_params();
            // Two literals were not NCBI's: blast_stat.c:410, the PAM70 ungapped beta is
            // -0.7 (the literal said -0.9); :384, the PAM30 10/1 alpha is 0.35 (0.34).
            if matrix == ScoringMatrix::Pam70 && entry.gap_open == i32::MAX {
                expected.beta = -0.7;
            }
            if matrix == ScoringMatrix::Pam30 && (entry.gap_open, entry.gap_extend) == (10, 1) {
                expected.alpha = 0.35;
            }
            let got = if entry.gap_open == i32::MAX {
                lookup_protein_params_ungapped(matrix)
            } else {
                let spec = ProteinScoringSpec {
                    matrix,
                    gap_open: entry.gap_open,
                    gap_extend: entry.gap_extend,
                };
                assert!(protein_scoring_supported(&spec), "{matrix:?} {entry:?}");
                lookup_protein_params(&spec)
            };
            assert_eq!(bits(got), bits(expected), "{matrix:?} {entry:?}");
            rows += 1;
        }
    }
    assert_eq!(rows, 90);
}

#[test]
fn every_gap_cost_row_of_ncbis_tables_is_found() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:183-587: 98 rows, 9 of them ungapped. The literal tables had 82
    // gapped rows; PAM30 15/3, 14/2, 14/1, 13/3, PAM70 11/2, 12/3 and IDENTITY 15/2 are
    // the others.
    let mut gapped = 0;
    for matrix in ScoringMatrix::ALL {
        let table = scoring_matrix_table(matrix);
        assert_eq!(table.name, matrix.ncbi_name());
        for row in &table.values[1..] {
            let spec = ProteinScoringSpec {
                matrix,
                gap_open: row.0[0] as i32,
                gap_extend: row.0[1] as i32,
            };
            let params = lookup_protein_params(&spec);
            assert_eq!(
                [params.lambda, params.k, params.h, params.alpha, params.beta],
                [row.0[3], row.0[4], row.0[5], row.0[6], row.0[7]]
            );
            gapped += 1;
        }
    }
    assert_eq!(gapped, 89);
    let pam30 = lookup_protein_params(&ProteinScoringSpec {
        matrix: ScoringMatrix::Pam30,
        gap_open: 15,
        gap_extend: 3,
    });
    assert_eq!((pam30.lambda, pam30.k, pam30.h), (0.339, 0.28, 1.70));
    let identity = lookup_protein_params(&ProteinScoringSpec {
        matrix: ScoringMatrix::Identity,
        gap_open: 15,
        gap_extend: 2,
    });
    assert_eq!(
        (identity.lambda, identity.k, identity.h),
        (0.2835, 0.255, 1.49)
    );
    assert!(!protein_scoring_supported(&ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: 5,
        gap_extend: 5,
    }));
}

#[test]
fn the_ungapped_rows_costs_get_no_alpha_and_beta() {
    // BLAST_GetAlphaBeta (blast_stat.c:3117-3126) searches from row 1; BLAST_CalcEffLengths
    // starts alpha and beta at 0 (blast_setup.c:705). Lambda, K and H are row 0's
    // (Blast_KarlinBlkGappedLoadFromTables searches from row 0).
    for matrix in ScoringMatrix::ALL {
        let spec = ProteinScoringSpec {
            matrix,
            gap_open: 32767,
            gap_extend: 32767,
        };
        assert!(protein_scoring_supported(&spec));
        let params = lookup_protein_params(&spec);
        let ungapped = lookup_protein_params_ungapped(matrix);
        assert_eq!(
            [params.lambda, params.k, params.h, params.alpha, params.beta],
            [ungapped.lambda, ungapped.k, ungapped.h, 0.0, 0.0],
            "{matrix:?}"
        );
        for row in &scoring_matrix_table(matrix).values[1..] {
            assert!(row.0[0] != 32767.0 && row.0[1] != 32767.0);
        }
    }
}

#[test]
fn gapped_lookup_reads_the_best_row() {
    let blosum62 = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
    assert_eq!(
        (blosum62.lambda, blosum62.k, blosum62.alpha, blosum62.beta),
        (0.267, 0.041, 1.9, -30.0)
    );
    let blosum45 = lookup_protein_params_gapped(ScoringMatrix::Blosum45);
    assert_eq!((blosum45.lambda, blosum45.k), (0.195, 0.032));
}

#[test]
#[should_panic(expected = "gap costs 5/5 are not in NCBI's BLOSUM62 table")]
fn gap_costs_outside_the_table_are_a_defect() {
    lookup_protein_params(&ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: 5,
        gap_extend: 5,
    });
}

#[test]
fn test_lookup_protein_params_blosum62() {
    let spec = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: 11,
        gap_extend: 1,
    };
    let params = lookup_protein_params(&spec);
    assert!((params.lambda - 0.267).abs() < 0.01);
}
