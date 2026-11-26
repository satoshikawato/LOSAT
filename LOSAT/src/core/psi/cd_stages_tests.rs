//! Tests of the conserved-domain stages (`cd_stages.rs`): the validation of a profile, the
//! weighting of its columns and the PSSM that `create_pssm_from_cdd` makes of them.

use super::*;
use crate::config::ScoringMatrix;
use crate::core::psi::KarlinBlk;

/// BLOSUM62 11/1's gapped block (`Blast_KarlinBlkGappedCalc`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:268
const BLOSUM62_11_1: KarlinBlk = KarlinBlk {
    lambda: 0.267,
    k: 0.041,
    h: 0.14,
};

fn blosum62() -> PsiScoreBlk {
    crate::algorithm::pssm_search::in_pssm::score_block(ScoringMatrix::Blosum62, BLOSUM62_11_1)
        .unwrap()
}

/// A column whose frequency of `letter` is 1.
fn one_letter_column(letter: u8, iobsr: f64) -> CdColumn {
    let mut wfreqs = [0.0; BLASTAA_SIZE];
    wfreqs[letter as usize] = 1.0;
    CdColumn { wfreqs, iobsr }
}

/// A query of NCBIstdaa letters A, C, D, E, F.
const QUERY: [u8; 5] = [1, 3, 4, 5, 6];

fn one_domain(columns: Vec<CdColumn>) -> PsiCdMsa {
    PsiCdMsa {
        query: QUERY.to_vec(),
        rows: vec![(0..QUERY.len()).map(CdCell::aligned).collect()],
        columns,
    }
}

/// The scores of NCBI BLAST+ 2.17.0's PSSM of the query ACDEF without domain hits
/// (`deltablast -query q_acdef.faa -db dbs -rpsdb deltatest -out_pssm p.asn`, the file's
/// `finalData scores`, probe `acdef` of PS13a's fix round 2a, `PR/ps13/fix2_probe/`): one row of
/// the 28 NCBIstdaa letters for each query position.
const ACDEF_SCORES: [[i32; BLASTAA_SIZE]; 5] = [
    [
        -32768, 4, -2, 0, -2, -1, -2, 0, -2, -1, -1, -2, -1, -2, -1, -1, -2, 1, 0, 0, -3, -1, -2,
        -1, -1, -4, -1, -2,
    ],
    [
        -32768, 0, -3, 10, -4, -4, -3, -3, -3, -1, -3, -1, -2, -3, -3, -3, -4, -1, -1, -1, -3, -1,
        -3, -4, -1, -4, -1, -1,
    ],
    [
        -32768, -2, 5, -4, 6, 2, -4, -1, -1, -3, -1, -4, -3, 1, -2, 0, -2, 0, -1, -4, -5, -1, -3,
        1, -1, -4, -1, -4,
    ],
    [
        -32768, -1, 1, -4, 2, 5, -4, -2, 0, -4, 1, -3, -2, 0, -1, 2, 0, 0, -1, -3, -3, -1, -2, 5,
        -1, -4, -1, -3,
    ],
    [
        -32768, -2, -4, -3, -4, -4, 7, -3, -1, 0, -3, 0, 0, -3, -4, -4, -3, -3, -2, -1, 1, -1, 3,
        -4, -1, -4, -1, 0,
    ],
];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:385-414
// A query without domain hits has frequency ratios of zero, so its PSSM takes the matrix's scores
// of its letters, which the scaling to the query's ungapped Lambda then changes (A against L is
// -1 in BLOSUM62, -2 in the PSSM): NCBI's scores, not the matrix's.
#[test]
fn a_query_without_domains_gets_the_matrix_scores_scaled_as_ncbi() {
    let sbp = blosum62();
    let msa = PsiCdMsa {
        query: QUERY.to_vec(),
        rows: Vec::new(),
        columns: Vec::new(),
    };
    let pssm = create_pssm_from_cdd(&msa, 0, &sbp).unwrap();
    assert_eq!(pssm.scores.len(), ACDEF_SCORES.len());
    for (position, (row, expected)) in pssm.scores.iter().zip(&ACDEF_SCORES).enumerate() {
        assert_eq!(&row[..], &expected[..], "position {position}");
    }
    assert_ne!(pssm.scores[0], sbp.matrix[QUERY[0] as usize]);
    // The file's lambdaUngapped and kappaUngapped (15 digits).
    assert!((pssm.ungapped.lambda - 0.320004687599798).abs() < 1e-14);
    assert!((pssm.ungapped.k - 0.146268722049472).abs() < 1e-14);
    let ratios = pssm.diagnostics.frequency_ratios.unwrap();
    assert!(ratios.iter().flatten().all(|&ratio| ratio == 0.0));
}

#[test]
fn a_gap_in_the_query_is_rejected() {
    let mut msa = one_domain(QUERY.iter().map(|&l| one_letter_column(l, 5.0)).collect());
    msa.query[2] = GAP;
    assert_eq!(validate_cd_msa(&msa), Err(PsiError::GapInQuery));
}

#[test]
fn a_column_with_too_few_observations_is_rejected() {
    let mut columns: Vec<CdColumn> = QUERY.iter().map(|&l| one_letter_column(l, 5.0)).collect();
    columns[3].iobsr = K_EPSILON / 2.0;
    assert_eq!(
        validate_cd_msa(&one_domain(columns)),
        Err(PsiError::BadProfile)
    );
}

#[test]
fn a_column_whose_frequencies_do_not_sum_to_one_is_rejected() {
    let mut columns: Vec<CdColumn> = QUERY.iter().map(|&l| one_letter_column(l, 5.0)).collect();
    columns[1].wfreqs[3] = 0.5;
    assert_eq!(
        validate_cd_msa(&one_domain(columns)),
        Err(PsiError::BadProfile)
    );
}

#[test]
fn a_negative_frequency_is_rejected() {
    let mut columns: Vec<CdColumn> = QUERY.iter().map(|&l| one_letter_column(l, 5.0)).collect();
    columns[0].wfreqs[1] = 1.5;
    columns[0].wfreqs[2] = -0.5;
    assert_eq!(
        validate_cd_msa(&one_domain(columns)),
        Err(PsiError::BadProfile)
    );
}

#[test]
fn an_unaligned_cell_is_not_validated() {
    let mut msa = one_domain(vec![one_letter_column(1, 5.0)]);
    msa.rows[0] = vec![CdCell::aligned(0)]
        .into_iter()
        .chain(std::iter::repeat_n(CdCell::UNALIGNED, QUERY.len() - 1))
        .collect();
    assert_eq!(validate_cd_msa(&msa), Ok(()));
}

#[test]
fn the_query_letter_joins_a_column_that_lacks_it() {
    // One domain that shows only G at the first position, where the query has A: the weighted
    // frequencies are (G: 5, A: 1) / 6 and the observations 6.
    let mut msa = one_domain(vec![one_letter_column(7, 5.0)]);
    msa.rows[0] = vec![CdCell::aligned(0)]
        .into_iter()
        .chain(std::iter::repeat_n(CdCell::UNALIGNED, QUERY.len() - 1))
        .collect();
    let weights = compute_frequencies_from_cds(&msa);
    assert_eq!(weights.independent_observations[0], 6.0);
    assert_eq!(weights.match_weights[0][7], 5.0 / 6.0);
    assert_eq!(weights.match_weights[0][1], 1.0 / 6.0);
    // The columns without a domain have no weight.
    assert_eq!(weights.independent_observations[1], 0.0);
    assert!(weights.match_weights[1].iter().all(|&weight| weight == 0.0));
}

#[test]
fn the_observations_stop_at_400() {
    let msa = one_domain(
        QUERY
            .iter()
            .map(|&l| one_letter_column(l, 1000.0))
            .collect(),
    );
    let weights = compute_frequencies_from_cds(&msa);
    assert!(weights.independent_observations.iter().all(|&n| n == 400.0));
}

#[test]
fn a_conserved_column_scores_its_letter_above_the_matrix() {
    let sbp = blosum62();
    let msa = one_domain(QUERY.iter().map(|&l| one_letter_column(l, 50.0)).collect());
    let pssm = create_pssm_from_cdd(&msa, 0, &sbp).unwrap();
    for (row, &letter) in pssm.scores.iter().zip(&QUERY) {
        let letter = letter as usize;
        assert!(row[letter] > sbp.matrix[letter][letter], "{row:?}");
    }
    assert!(pssm.ungapped.lambda > 0.0 && pssm.gapped.k > 0.0);
}

#[test]
fn an_x_in_the_query_has_no_ratios() {
    let sbp = blosum62();
    let mut msa = one_domain(QUERY.iter().map(|&l| one_letter_column(l, 50.0)).collect());
    msa.query[2] = X;
    let pssm = create_pssm_from_cdd(&msa, 0, &sbp).unwrap();
    let ratios = pssm.diagnostics.frequency_ratios.unwrap();
    assert!(ratios[2].iter().all(|&ratio| ratio == 0.0));
    assert!(ratios[1].iter().any(|&ratio| ratio != 0.0));
}

#[test]
fn a_zero_count_nan_column_on_x_matches_ordinary_x_column() {
    let sbp = blosum62();
    let mut msa_ok = one_domain(QUERY.iter().map(|&l| one_letter_column(l, 50.0)).collect());
    msa_ok.query[2] = X;
    let mut msa_nan = msa_ok.clone();
    msa_nan.columns[2].wfreqs = [f64::NAN; BLASTAA_SIZE];
    assert!(validate_cd_msa(&msa_nan).is_ok());
    let pssm_ok = create_pssm_from_cdd(&msa_ok, 0, &sbp).unwrap();
    let pssm_nan = create_pssm_from_cdd(&msa_nan, 0, &sbp).unwrap();
    assert_eq!(pssm_ok.scores, pssm_nan.scores);
    assert_eq!(
        pssm_ok.diagnostics.frequency_ratios,
        pssm_nan.diagnostics.frequency_ratios
    );
}
