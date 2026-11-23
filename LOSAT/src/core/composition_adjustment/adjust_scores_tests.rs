use super::*;
use crate::core::composition_adjustment::redo_alignment::redo_calc_lambda;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1204-1248:
#[test]
fn translated_composition_range_honors_stop_margin_strict_boundaries() {
    let mut data = vec![1u8; 100];
    assert_eq!(
        blast_get_composition_range(&data, 50, 60).unwrap(),
        (0, 100)
    );
    data[20] = 25;
    data[80] = 25;
    assert_eq!(
        blast_get_composition_range(&data, 50, 60).unwrap(),
        (41, 60)
    );
    assert_eq!(
        blast_get_composition_range(&data, 50, 50).unwrap(),
        (41, 60)
    );
    data[29] = 25;
    assert_eq!(
        blast_get_composition_range(&data, 50, 50).unwrap(),
        (50, 60)
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1234-1254;
// composition_adjustment/composition_adjustment.c:1414-1530:
#[test]
fn pinned_tblastn_mode2_adjustment_rules_match_ncbi_call_inputs() {
    for case in [
        "seg_hard_query_20260924_default",
        "multi_query_20260924_default",
    ] {
        let trace = std::fs::read_to_string(format!(
            "{}/../docs/evidence/tlosan_stage_d/kappa_composition_matrix_scores_20260924/{case}.tsv",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let rows: Vec<_> = trace.lines().collect();
        assert_eq!(rows.len() % 6, 0);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2229-2236:
        // Initializes `Blast_MatrixInfo` from `startFreqRatios` and `ungappedLambda`.
        let matrix_row: Vec<_> = rows[1].split('\t').collect();
        assert_eq!(matrix_row[0], "K_MATRIX");
        let matrix_lambda: f64 = matrix_row[5].parse().unwrap();
        let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, matrix_lambda).unwrap();
        let mut workspace = BlastCompositionWorkspace::new_blosum62();
        for block in rows.chunks_exact(6) {
            let call: Vec<_> = block[0].split('\t').collect();
            let matrix: Vec<_> = block[1].split('\t').collect();
            let result: Vec<_> = block[4].split('\t').collect();
            assert_eq!(matrix[0], "K_MATRIX");
            assert_eq!(matrix[1], call[1]);
            assert_eq!(matrix[2], "28");
            assert_eq!(matrix[3], "28");
            assert_eq!(matrix[4], "0");
            assert_eq!(
                matrix[5].parse::<f64>().unwrap().to_bits(),
                matrix_lambda.to_bits()
            );
            assert_eq!(matrix_info.start_matrix[1][1], matrix[6].parse().unwrap());
            assert_eq!(matrix_info.start_matrix[1][16], matrix[7].parse().unwrap());
            assert_eq!(matrix_info.start_matrix[22][22], matrix[8].parse().unwrap());
            assert_eq!(call[0], "K_COMP_CALL");
            assert_eq!(call[4], "2");
            assert_eq!(call[5], "20");
            assert_eq!(result[0], "K_COMP_RESULT");
            assert_eq!(call[1], result[1]);
            assert_eq!(result[2], "0");
            let compositions: Vec<_> = block[2..4]
                .iter()
                .enumerate()
                .map(|(idx, row)| {
                    let fields: Vec<_> = row.split('\t').collect();
                    assert_eq!(fields[0], "K_COMP");
                    assert_eq!(fields[1], call[1]);
                    assert_eq!(fields[2], if idx == 0 { "query" } else { "subject" });
                    let mut comp = BlastAminoAcidComposition::empty();
                    comp.num_true_amino_acids = fields[3].parse().unwrap();
                    assert_eq!(fields.len(), 32);
                    for (value, hex) in comp.prob.iter_mut().zip(&fields[4..]) {
                        *value = f64::from_bits(u64::from_str_radix(hex, 16).unwrap());
                    }
                    comp
                })
                .collect();
            let adjusted = blast_adjust_scores(
                &matrix_info,
                &compositions[0],
                call[2].parse().unwrap(),
                &compositions[1],
                call[3].parse().unwrap(),
                BlastCompoAdjustMode::CompositionMatrixAdjust,
                call[6].parse().unwrap(),
                &mut workspace,
                redo_calc_lambda,
            )
            .unwrap()
            .unwrap();
            assert_eq!(
                adjusted.matrix_adjust_rule as i32,
                result[3].parse().unwrap(),
                "{case}: call {}",
                call[1]
            );
            assert_eq!(
                adjusted.lambda_ratio.to_bits(),
                result[4].parse::<f64>().unwrap().to_bits(),
                "{case}: call {}",
                call[1]
            );
            let adjusted_row: Vec<_> = block[5].split('\t').collect();
            assert_eq!(adjusted_row[0], "K_ADJUSTED");
            assert_eq!(adjusted_row[1], call[1]);
            assert_eq!(adjusted_row.len(), 2 + 28 * 28);
            for (index, (observed, expected)) in adjusted
                .adjusted_matrix
                .scores
                .iter()
                .flatten()
                .zip(&adjusted_row[2..])
                .enumerate()
            {
                assert_eq!(
                    *observed,
                    expected.parse::<i32>().unwrap(),
                    "{case}: call {}, matrix index {index}",
                    call[1]
                );
            }
        }
    }
}

#[test]
fn test_read_aa_composition_ignores_x_and_merges_u_into_c() {
    let seq = [1u8, 3, 21, 24];
    let composition = read_aa_composition(&seq);
    assert_eq!(composition.num_true_amino_acids, 3);
    assert!((composition.prob[1] - (1.0 / 3.0)).abs() < 1.0e-12);
    assert!((composition.prob[3] - (2.0 / 3.0)).abs() < 1.0e-12);
    assert_eq!(composition.prob[21], 0.0);
    assert_eq!(composition.prob[24], 0.0);
}

#[test]
fn test_build_matrix_info_blosum62_rebuilds_start_matrix_from_freq_ratios() {
    let info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176).unwrap();
    assert!((info.start_freq_ratios[1][1] - 3.90294070).abs() < 1.0e-8);
    assert!((info.start_freq_ratios[22][22] - 9.83220341).abs() < 1.0e-8);
    assert_eq!(info.start_matrix[1][1], 4);
    assert_eq!(info.start_matrix[1][16], -2);
    assert!(info.start_freq_ratios[1][1] > 0.0);
}

// The literals this file held before the generated tables (SC0b): BLOSUM62's joint
// probabilities and background frequencies (matrix_frequency_data.c:383-537) and its
// frequency ratios (matrix_freq_ratios.c:438-634).
const LITERAL_BLOSUM62_JOINT_PROBS: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA] =
    include!("blosum62_joint_probs.inc.rs");
const LITERAL_BLOSUM62_BG: [f64; COMPO_NUM_TRUE_AA] = [
    7.4216205067993410e-02,
    5.1614486141284638e-02,
    4.4645808512757915e-02,
    5.3626000838554413e-02,
    2.4687457167944848e-02,
    3.4259650591416023e-02,
    5.4311925684587502e-02,
    7.4146941452644999e-02,
    2.6212984805266227e-02,
    6.7917367618953756e-02,
    9.8907868497150955e-02,
    5.8155682303079680e-02,
    2.4990197579643110e-02,
    4.7418459742284751e-02,
    3.8538003320306206e-02,
    5.7229029476494421e-02,
    5.0891364550287033e-02,
    1.3029956129972148e-02,
    3.2281512313758580e-02,
    7.2919098205619245e-02,
];
const LITERAL_BLOSUM62_START_FREQ_RATIOS: [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET] =
    include!("blosum62_start_freq_ratios.inc.rs");

fn bits<const N: usize>(values: &[f64; N]) -> [u64; N] {
    values.map(f64::to_bits)
}

#[test]
fn generated_blosum62_data_equal_the_former_literals_bit_for_bit() {
    let data = compo_matrix_data(ScoringMatrix::Blosum62).unwrap();
    for i in 0..COMPO_NUM_TRUE_AA {
        let literal = bits(&LITERAL_BLOSUM62_JOINT_PROBS[i]);
        assert_eq!(bits(&data.joint_probs[i]), literal, "row {i}");
        assert_eq!(bits(&BLOSUM62_JOINT_PROBS[i]), literal, "row {i}");
    }
    assert_eq!(bits(data.background), bits(&LITERAL_BLOSUM62_BG));
    assert_eq!(bits(&BLOSUM62_BG), bits(&LITERAL_BLOSUM62_BG));
    let info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176).unwrap();
    for i in 0..COMPO_LARGEST_ALPHABET {
        let literal = bits(&LITERAL_BLOSUM62_START_FREQ_RATIOS[i]);
        assert_eq!(bits(&data.freq_ratios[i]), literal, "row {i}");
        assert_eq!(bits(&info.start_freq_ratios[i]), literal, "row {i}");
    }
}

#[test]
fn eight_matrices_have_composition_data_and_identity_has_none() {
    assert_eq!(COMPO_MATRIX_DATA.len(), 8);
    for matrix in ScoringMatrix::ALL {
        match compo_matrix_data(matrix) {
            Some(data) => assert_eq!(data.name, matrix.ncbi_name()),
            None => assert_eq!(matrix, ScoringMatrix::Identity),
        }
    }
    for data in &COMPO_MATRIX_DATA {
        let matrix = data.name.parse().unwrap();
        assert!(build_matrix_info(matrix, 0.3).is_ok(), "{}", data.name);
        assert!(
            BlastCompositionWorkspace::new(matrix).is_ok(),
            "{}",
            data.name
        );
    }
    // blast_kappa.c:657-660 and composition_adjustment.c:1316-1331 return -1.
    assert_eq!(
        build_matrix_info(ScoringMatrix::Identity, 0.3)
            .unwrap_err()
            .to_string(),
        "composition adjustment has no frequency ratios for the IDENTITY matrix"
    );
    assert_eq!(
        BlastCompositionWorkspace::new(ScoringMatrix::Identity)
            .unwrap_err()
            .to_string(),
        "Matrix IDENTITY not currently supported for RE based adjustment"
    );
}

#[test]
fn workspace_standard_frequencies_are_the_sums_of_the_joint_probabilities() {
    let generic = BlastCompositionWorkspace::new(ScoringMatrix::Blosum62).unwrap();
    let blosum62 = BlastCompositionWorkspace::new_blosum62();
    for i in 0..COMPO_NUM_TRUE_AA {
        assert_eq!(bits(&generic.mat_b[i]), bits(&blosum62.mat_b[i]));
        assert_eq!(bits(&generic.mat_final[i]), bits(&blosum62.mat_final[i]));
    }
    assert_eq!(
        bits(&generic.first_standard_freq),
        bits(&blosum62.first_standard_freq)
    );
    assert_eq!(
        bits(&generic.second_standard_freq),
        bits(&blosum62.second_standard_freq)
    );
    // matrix_frequency_data.c:1406-1416 sums the joint probabilities. The matrices are
    // symmetric, so the row and column sums are the same; they equal the background
    // frequencies of the BLOSUM matrices and differ in the last bits for the PAM ones.
    for data in &COMPO_MATRIX_DATA {
        let matrix = data.name.parse().unwrap();
        let workspace = BlastCompositionWorkspace::new(matrix).unwrap();
        assert_eq!(
            bits(&workspace.first_standard_freq),
            bits(&workspace.second_standard_freq),
            "{}",
            data.name
        );
        let equal = (0..COMPO_NUM_TRUE_AA)
            .filter(|&i| workspace.first_standard_freq[i].to_bits() == data.background[i].to_bits())
            .count();
        let expected = match matrix {
            ScoringMatrix::Pam30 | ScoringMatrix::Pam70 => 1,
            ScoringMatrix::Pam250 => 0,
            _ => COMPO_NUM_TRUE_AA,
        };
        assert_eq!(equal, expected, "{}", data.name);
        for i in 0..COMPO_NUM_TRUE_AA {
            assert!((workspace.first_standard_freq[i] - data.background[i]).abs() < 1.0e-15);
        }
    }
}

#[test]
fn matrix_adjust_rule_measures_against_the_matrix_background() {
    // compo_mode_condition.c:188: the distances are to the matrix's background. A query
    // with BLOSUM62's background is at distance 0 from it (no angle: the relative-entropy
    // rule), but not from PAM30's, where the distant, four times longer match is at 87
    // degrees (the scaling rule).
    let mut subject = [0.8 / 19.0; COMPO_NUM_TRUE_AA];
    subject[11] = 0.2;
    let rule = |matrix| {
        blast_choose_matrix_adjust_rule(
            100,
            400,
            &BLOSUM62_BG,
            &subject,
            matrix,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
        )
    };
    assert_eq!(
        rule(ScoringMatrix::Blosum62).unwrap(),
        EMatrixAdjustRule::UserSpecifiedRelEntropy
    );
    assert_eq!(
        rule(ScoringMatrix::Pam30).unwrap(),
        EMatrixAdjustRule::CompoScaleOldMatrix
    );
    assert_eq!(
        rule(ScoringMatrix::Identity).unwrap_err().to_string(),
        "matrix IDENTITY is not supported for RE based adjustment"
    );
}

#[test]
fn matrix_adjust_rule_angle_follows_the_ncbi_expression() {
    // compo_mode_condition.c:200-204. A right angle (cosine 0) and a straight one (cosine -1)
    // are above 70 degrees.
    assert!(matrix_adjust_rule_angle_above_threshold(0.3, 0.4, 0.5));
    assert!(matrix_adjust_rule_angle_above_threshold(0.3, 0.2, 0.5));
    // A cosine below -1 (a collapsed triangle) and a zero distance give NaN, not 180 or 0.
    assert!(!matrix_adjust_rule_angle_above_threshold(
        0.3,
        0.2,
        0.5000000001
    ));
    assert!(!matrix_adjust_rule_angle_above_threshold(0.3, 0.0, 0.3));
    assert!(!matrix_adjust_rule_angle_above_threshold(0.0, 0.0, 0.0));
    // Around the limit: 70.0 degrees at it, 70.00000000000001 one double below.
    let below = f64::from_bits(ANGLE_COSINE_LIMIT.to_bits() - 1);
    assert!((-1.0..ANGLE_COSINE_LIMIT).contains(&below));
    assert!(!(-1.0..ANGLE_COSINE_LIMIT).contains(&ANGLE_COSINE_LIMIT));
}

#[cfg(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64"))]
#[test]
#[expect(
    clippy::disallowed_methods,
    reason = "checks the limit against glibc's own acos, NCBI's"
)]
fn matrix_adjust_rule_limit_is_glibcs() {
    // NCBI's expression with glibc's acos against the limit, for every double within 2^16
    // steps of it and for cosines across [-1.25, 1.25].
    let ncbi = |cosine: f64| cosine.acos() * HALF_CIRCLE_DEGREES / PI > ANGLE_DEGREE_THRESHOLD;
    let limit = ANGLE_COSINE_LIMIT.to_bits();
    for bits in limit - (1 << 16)..limit + (1 << 16) {
        let cosine = f64::from_bits(bits);
        assert_eq!(
            ncbi(cosine),
            (-1.0..ANGLE_COSINE_LIMIT).contains(&cosine),
            "{cosine}"
        );
    }
    for i in -125_000..=125_000 {
        let cosine = f64::from(i) * 1e-5;
        assert_eq!(
            ncbi(cosine),
            (-1.0..ANGLE_COSINE_LIMIT).contains(&cosine),
            "{cosine}"
        );
    }
}

// The query and match distances of the rule, and the cosine of its angle.
fn matrix_adjust_rule_cosine(
    query: &[f64; COMPO_NUM_TRUE_AA],
    subject: &[f64; COMPO_NUM_TRUE_AA],
) -> f64 {
    let d_m_mat = blast_get_relative_entropy(subject, &BLOSUM62_BG);
    let d_q_mat = blast_get_relative_entropy(query, &BLOSUM62_BG);
    let d_m_q = blast_get_relative_entropy(subject, query);
    (d_m_mat * d_m_mat + d_q_mat * d_q_mat - d_m_q * d_m_q) / 2.0 / d_m_mat / d_q_mat
}

#[test]
fn matrix_adjust_rule_takes_ncbi_degrees_at_the_threshold() {
    // A query shifted from BLOSUM62's background (letter 11 up, letter 2 down) against the
    // distant, four times longer match of the test above. At the first shift the angle is
    // 70.0000000079 degrees with the true pi (the cosine is below cos 70 degrees,
    // 0.3420201433256688) but 69.9999999921 with NCBI's PI (not below the limit), so NCBI
    // takes the relative-entropy rule; past NCBI's crossing it takes the scaling rule.
    let mut subject = [0.8 / 19.0; COMPO_NUM_TRUE_AA];
    subject[11] = 0.2;
    let at = |shift: f64| {
        let mut query = BLOSUM62_BG;
        query[11] += shift;
        query[2] -= shift;
        let cosine = matrix_adjust_rule_cosine(&query, &subject);
        let rule = blast_choose_matrix_adjust_rule(
            100,
            400,
            &query,
            &subject,
            ScoringMatrix::Blosum62,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
        )
        .unwrap();
        (cosine, rule)
    };
    let (cosine, rule) = at(0.03708335181724542);
    assert!(
        (ANGLE_COSINE_LIMIT..0.3420201433256688).contains(&cosine),
        "{cosine}"
    );
    assert_eq!(rule, EMatrixAdjustRule::UserSpecifiedRelEntropy);
    let (_, rule) = at(0.03708335187554269);
    assert_eq!(rule, EMatrixAdjustRule::CompoScaleOldMatrix);
}

#[test]
fn matrix_adjust_rule_collapsed_triangle_takes_relative_entropy() {
    // The query lies just past BLOSUM62's background on the line from the match, so the
    // background is between the two and the triangle of distances collapses: the cosine
    // rounds to -1.6. NCBI's acos is then NaN (the relative-entropy rule); a clamp to -1
    // (180 degrees) would take the scaling rule. The -1.6 is the rounding of `a / temp` in
    // the 2e-9 query distance (IEEE division, the same on every platform); a libm's log
    // error there is about 1e-24, six orders below the squared distance.
    let mut subject = [0.8 / 19.0; COMPO_NUM_TRUE_AA];
    subject[11] = 0.2;
    let query: [f64; COMPO_NUM_TRUE_AA] =
        std::array::from_fn(|i| BLOSUM62_BG[i] - 1.4e-8 * (subject[i] - BLOSUM62_BG[i]));
    let cosine = matrix_adjust_rule_cosine(&query, &subject);
    assert!(cosine < -1.0, "{cosine}");
    assert_eq!(
        blast_choose_matrix_adjust_rule(
            100,
            400,
            &query,
            &subject,
            ScoringMatrix::Blosum62,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
        )
        .unwrap(),
        EMatrixAdjustRule::UserSpecifiedRelEntropy
    );
}

#[test]
fn test_blast_composition_pvalue_clamps_left_tail() {
    assert_eq!(
        blast_composition_pvalue(COMPO_MIN_LAMBDA),
        P_LAMBDA_TABLE[LOW_LAMBDA_BIN_CUT]
    );
}

#[test]
fn test_blast_overall_p_value_matches_fisher_style_formula() {
    let combined = blast_overall_p_value(0.2, 0.3);
    let product = 0.06f64;
    assert!((combined - (product * (1.0 - product.glibc_ln()))).abs() < 1.0e-12);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
// The ignored test times the port of this function against `optimize_target_frequencies_exact`
// and checks the last result bit for bit.
/// Timing of the Newton kernels on realistic inputs (run with --ignored --nocapture;
/// LOSAT_X_NEWTONEXACT_LEVEL selects the vector width of the exact one).
#[test]
#[ignore]
fn x_bench_newton_kernels() {
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let workspace = BlastCompositionWorkspace::new_blosum62();
    let mut q = vec![0.0f64; COMPO_NUM_TRUE_AA * COMPO_NUM_TRUE_AA];
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            q[i * COMPO_NUM_TRUE_AA + j] = workspace.mat_b[i][j];
        }
    }
    let cases: Vec<([f64; COMPO_NUM_TRUE_AA], [f64; COMPO_NUM_TRUE_AA])> = (0..256)
        .map(|_| {
            let mut row = [0.0f64; COMPO_NUM_TRUE_AA];
            let mut col = [0.0f64; COMPO_NUM_TRUE_AA];
            let (mut sr, mut sc) = (0.0, 0.0);
            for k in 0..COMPO_NUM_TRUE_AA {
                row[k] = BLOSUM62_BG[k] * (0.3 + 1.4 * next());
                col[k] = BLOSUM62_BG[k] * (0.3 + 1.4 * next());
                sr += row[k];
                sc += col[k];
            }
            for k in 0..COMPO_NUM_TRUE_AA {
                row[k] /= sr;
                col[k] /= sc;
            }
            apply_pseudocounts(&mut row, 300, &BLOSUM62_BG, 20);
            apply_pseudocounts(&mut col, 250, &BLOSUM62_BG, 20);
            (row, col)
        })
        .collect();
    let mut x_ref = vec![0.0f64; q.len()];
    let mut x_new = vec![0.0f64; q.len()];
    let reps = 8;
    let t0 = std::time::Instant::now();
    let mut status_sum = 0;
    for _ in 0..reps {
        for (row, col) in &cases {
            status_sum += optimize_target_frequencies_reference(
                &mut x_ref,
                COMPO_NUM_TRUE_AA,
                &q,
                row,
                col,
                true,
                K_FIXED_RE_BLOSUM62,
            );
        }
    }
    let t_ref = t0.elapsed().as_secs_f64() / (reps * cases.len()) as f64 * 1e6;
    let t0 = std::time::Instant::now();
    let mut status_sum2 = 0;
    for _ in 0..reps {
        for (row, col) in &cases {
            status_sum2 += super::super::x_newton_exact::optimize_target_frequencies_exact(
                &mut x_new,
                &q,
                row,
                col,
                K_FIXED_RE_BLOSUM62,
            );
        }
    }
    let t_new = t0.elapsed().as_secs_f64() / (reps * cases.len()) as f64 * 1e6;
    // bitwise check on the last case
    for (a, b) in x_ref.iter().zip(x_new.iter()) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    eprintln!(
        "[BENCH] reference {:.1} us/call (status sum {}), exact {:.1} us/call (status sum {})",
        t_ref, status_sum, t_new, status_sum2
    );
    super::super::x_newton_exact::prof_print();
}

#[test]
fn test_blast_adjust_scores_returns_adjusted_matrix_for_blosum62() {
    let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176).unwrap();
    let query = read_aa_composition(&[1u8, 1, 1, 4, 4, 4, 17, 17]);
    let subject = read_aa_composition(&[1u8, 4, 17, 1, 4, 17, 1, 4]);
    let mut workspace = BlastCompositionWorkspace::new_blosum62();
    let adjusted = blast_adjust_scores(
        &matrix_info,
        &query,
        8,
        &subject,
        8,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        0,
        &mut workspace,
        redo_calc_lambda,
    )
    .unwrap()
    .expect("adjusted matrix result");
    assert_eq!(
        adjusted.matrix_adjust_rule,
        EMatrixAdjustRule::UserSpecifiedRelEntropy
    );
    assert_eq!(adjusted.lambda_ratio, 1.0);
    assert_ne!(
        adjusted.adjusted_matrix.scores[1][1],
        matrix_info.start_matrix[1][1]
    );
}

#[test]
fn test_blast_adjust_scores_returns_none_for_zero_composition_ncbi_style() {
    let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176).unwrap();
    let zero = read_aa_composition(&[21u8, 21, 21]);
    let subject = read_aa_composition(&[1u8, 4, 17]);
    let mut workspace = BlastCompositionWorkspace::new_blosum62();
    let adjusted = blast_adjust_scores(
        &matrix_info,
        &zero,
        3,
        &subject,
        3,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        0,
        &mut workspace,
        redo_calc_lambda,
    )
    .unwrap();
    assert!(adjusted.is_none());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:682-690
#[test]
fn test_set_xuo_scores_uses_row_average_for_pairwise_ambiguity_to_x() {
    let mut scores = [[0.0; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET];
    let mut row_probs = [0.0; COMPO_LARGEST_ALPHABET];
    let mut col_probs = [0.0; COMPO_LARGEST_ALPHABET];

    scores[E_BCHAR][E_DCHAR] = -10.0;
    scores[E_BCHAR][E_NCHAR] = -30.0;
    scores[E_DCHAR][E_BCHAR] = -7.0;
    scores[E_NCHAR][E_BCHAR] = -11.0;

    col_probs[E_DCHAR] = 0.75;
    col_probs[E_NCHAR] = 0.25;
    row_probs[E_DCHAR] = 0.25;
    row_probs[E_NCHAR] = 0.75;

    set_xuo_scores(&mut scores, &row_probs, &col_probs);

    assert!((scores[E_BCHAR][E_XCHAR] - (-15.0)).abs() < 1.0e-12);
    assert!((scores[E_XCHAR][E_BCHAR] - (-10.0)).abs() < 1.0e-12);
}

/// NCBI's `Blast_CalcLambdaFullPrecision`, statement by statement. The
/// third value counts the Newton steps it refused because f did not fall
/// below 0.9 of the previous iteration's f (`fold`).
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:220-334
fn ncbi_calc_lambda_full_precision(
    score: &[Vec<f64>],
    alphsize: usize,
    row_prob: &[f64],
    col_prob: &[f64],
) -> (f64, usize, usize) {
    let mut f = 4.0;
    let mut left = 0.0;
    let mut right = 1.0;
    let mut x = 0.367879441171f64;
    let mut is_newton = false;
    let mut max_score = COMPO_SCORE_MIN as f64;
    let mut avg_score = 0.0;
    for i in 0..alphsize {
        for j in 0..alphsize {
            if row_prob[i] == 0.0 || col_prob[j] == 0.0 {
                continue;
            }
            if max_score < score[i][j] {
                max_score = score[i][j];
            }
            avg_score += row_prob[i] * col_prob[j] * score[i][j];
        }
    }
    assert!(max_score > 0.0 && avg_score < 0.0, "lambda must exist");
    let mut fold_refusals = 0;
    let mut k = 0;
    while k < K_LAMBDA_ITERATION_LIMIT {
        let fold = f;
        let lambda = -x.glibc_ln();
        let was_newton = is_newton;
        let x_pow_max_score = (-max_score * lambda).glibc_exp();
        f = -x_pow_max_score;
        let mut slope = max_score * f / x;
        for i in 0..alphsize {
            for j in 0..alphsize {
                if row_prob[i] == 0.0 || col_prob[j] == 0.0 {
                    continue;
                }
                let ff = if max_score != score[i][j] {
                    let diff_score = max_score - score[i][j];
                    let ff = row_prob[i] * col_prob[j] * (-lambda * diff_score).glibc_exp();
                    slope += diff_score * ff / x;
                    ff
                } else {
                    row_prob[i] * col_prob[j]
                };
                f += ff;
            }
        }
        if f > 0.0 {
            left = x;
        } else if f < 0.0 {
            right = x;
        } else {
            break;
        }
        if right - left <= 2.0 * left * (1.0 - right) * K_LAMBDA_ERROR_TOLERANCE
            && (f / x_pow_max_score).abs() <= K_LAMBDA_FUNCTION_TOLERANCE
        {
            x = (left + right) / 2.0;
            break;
        }
        let f_fell_too_little = was_newton && f.abs() > 0.9 * fold.abs();
        if f_fell_too_little || slope >= 0.0 {
            if f_fell_too_little && slope < 0.0 {
                fold_refusals += 1;
            }
            x = (left + right) / 2.0;
        } else {
            let p = -f / slope;
            let y = x + p;
            if y <= left || y >= right {
                x = (left + right) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() <= K_LAMBDA_ERROR_TOLERANCE * x * (1.0 - x)
                    && (f / x_pow_max_score).abs() <= K_LAMBDA_FUNCTION_TOLERANCE
                {
                    break;
                }
            }
        }
        k += 1;
    }
    let lambda = if k < K_LAMBDA_ITERATION_LIMIT {
        -x.glibc_ln()
    } else {
        -1.0
    };
    (lambda, k, fold_refusals)
}

// The two matrices make NCBI refuse one Newton step because f did not
// fall below 0.9 of the previous f. Comparing with -x^max_score instead
// (the order LOSAT had before) takes 6 and 5 iterations, and the first
// lambda differs in its last bits.
#[test]
fn calc_lambda_full_precision_compares_f_with_the_previous_iteration_like_ncbi() {
    let cases = [
        (
            "a",
            [[-2.0, -6.0], [-2.0, 4.0]],
            [5.0 / 7.0, 2.0 / 7.0],
            [2.0 / 3.0, 1.0 / 3.0],
            3,
        ),
        ("b", [[-4.0, -5.0], [-4.0, 4.0]], [0.8, 0.2], [0.5, 0.5], 4),
    ];
    for (name, matrix, rows, cols, iterations) in cases {
        let score: Vec<Vec<f64>> = matrix.iter().map(|row| row.to_vec()).collect();
        let mut row_prob = [0.0; COMPO_NUM_TRUE_AA];
        let mut col_prob = [0.0; COMPO_NUM_TRUE_AA];
        row_prob[..2].copy_from_slice(&rows);
        col_prob[..2].copy_from_slice(&cols);

        let (ncbi_lambda, ncbi_iterations, fold_refusals) =
            ncbi_calc_lambda_full_precision(&score, 2, &row_prob, &col_prob);
        assert_eq!(fold_refusals, 1, "case {name} reaches the refusal");
        assert_eq!(ncbi_iterations, iterations, "case {name}");

        let (lambda, iter_count) =
            calc_lambda_full_precision_from_scores(&score, 2, &row_prob, &col_prob);
        assert_eq!(iter_count, ncbi_iterations, "case {name} iterations");
        assert_eq!(
            lambda.to_bits(),
            ncbi_lambda.to_bits(),
            "case {name}: lambda {lambda} vs NCBI {ncbi_lambda}"
        );

        // lambda solves sum_ij p_i q_j exp(lambda s_ij) = 1.
        let mut phi = -1.0;
        for i in 0..2 {
            for j in 0..2 {
                phi += rows[i] * cols[j] * (lambda * matrix[i][j]).glibc_exp();
            }
        }
        assert!(phi.abs() <= 1.0e-12, "case {name}: phi {phi}");
    }
}
