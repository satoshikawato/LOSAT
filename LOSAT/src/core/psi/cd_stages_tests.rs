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

#[test]
fn a_query_without_domains_keeps_the_matrix_scores() {
    let sbp = blosum62();
    let msa = PsiCdMsa {
        query: QUERY.to_vec(),
        rows: Vec::new(),
        columns: Vec::new(),
    };
    let pssm = create_pssm_from_cdd(&msa, 0, &sbp).unwrap();
    for (row, &letter) in pssm.scores.iter().zip(&QUERY) {
        assert_eq!(*row, sbp.matrix[letter as usize]);
    }
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
        .chain(std::iter::repeat(CdCell::UNALIGNED).take(QUERY.len() - 1))
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
        .chain(std::iter::repeat(CdCell::UNALIGNED).take(QUERY.len() - 1))
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
