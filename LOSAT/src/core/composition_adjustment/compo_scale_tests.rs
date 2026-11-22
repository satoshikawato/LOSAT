use super::*;

// `s_GetScoreRange` keeps 0 in the range and skips `COMPO_SCORE_MIN`.
#[test]
fn the_score_range_holds_zero_and_skips_the_minimum() {
    let mut row = [5; BLASTAA_SIZE];
    row[1] = COMPO_SCORE_MIN;
    row[3] = -7;
    // Not a true amino acid: not read.
    row[21] = -100;
    assert_eq!(score_range::<XDropRedo>(&[row]), (-7, 5));
    assert_eq!(score_range::<XDropRedo>(&[[2; BLASTAA_SIZE]]), (0, 2));
}

// `s_RoundScoreMatrix`: `Nint` rounds half away from zero; below `INT_MIN` is `INT_MIN`.
#[test]
fn scores_round_as_ncbis_nint() {
    assert_eq!(round_score::<XDropRedo>(2.5), 3);
    assert_eq!(round_score::<XDropRedo>(-2.5), -3);
    assert_eq!(round_score::<XDropRedo>(-2.49), -2);
    assert_eq!(
        round_score::<XDropRedo>(COMPO_SCORE_MIN as f64),
        COMPO_SCORE_MIN
    );
    assert_eq!(round_score::<XDropRedo>(-1.0e12), i32::MIN);
}

// `s_ScalePSSM` on one hand-checked row: ratio 0 gives `COMPO_SCORE_MIN`, X and O the
// average (capped at -1), U the C score, `*` the start matrix's.
#[test]
fn a_scaled_row_is_ncbis() {
    let mut ratios = [1.0; BLASTAA_SIZE];
    ratios[0] = 0.0;
    ratios[1] = std::f64::consts::E; // A: 1/lambda
    ratios[3] = std::f64::consts::E * std::f64::consts::E; // C: 2/lambda
    let mut start = [0; BLASTAA_SIZE];
    start[ncbistdaa::STOP as usize] = -77;
    let mut probs = [0.0; BLASTAA_SIZE];
    probs[1] = 0.5;
    probs[3] = 0.5;
    let mut matrix = [[0; BLASTAA_SIZE]];
    scale_pssm::<XDropRedo>(&mut matrix, &[ratios], &[start], &probs, 0.5);
    let row = matrix[0];
    assert_eq!(row[0], COMPO_SCORE_MIN);
    assert_eq!(row[1], 2);
    assert_eq!(row[3], 4);
    assert_eq!(row[ncbistdaa::U as usize], 4);
    // The average is 3 > -1: X and O are -1.
    assert_eq!(row[ncbistdaa::X as usize], -1);
    assert_eq!(row[ncbistdaa::O as usize], -1);
    assert_eq!(row[ncbistdaa::STOP as usize], -77);
    assert_eq!(row[5], 0);
}

fn gunzip(name: &str) -> Vec<u8> {
    use std::io::Read;
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name);
    let mut out = Vec::new();
    flate2::read::GzDecoder::new(std::fs::read(path).unwrap().as_slice())
        .read_to_end(&mut out)
        .unwrap();
    out
}

// Design question Q1 of PS4b, by an LD_PRELOAD probe of NCBI 2.17.0 (`psiblast -in_pssm
// tk.asn -db paja -comp_based_stats 1 -num_threads 1`, `PS/tmp/q1/`): in
// `Blast_RedoAlignmentCore_MT`, `Blast_KarlinBlkUngappedCalc` of the factor-1 PSSM gives
// Lambda 0.32206865605351226 and K 0.13750086820961227, and every `BLAST_SpougeStoE` of the
// redo gets Lambda 0.267/32 and K 0.042084934723385127 (= K * 0.041 / 0.13395614448848181).
#[test]
fn the_matrix_information_of_the_tk_pssm_is_ncbis() {
    let pssm =
        crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip("pssm/ps4_tk_paja.asn.gz"))
            .unwrap();
    let search = crate::algorithm::pssm_search::in_pssm::psi_search_of(pssm).unwrap();
    let query: Vec<u8> = search
        .query
        .seq()
        .iter()
        .map(|&c| crate::utils::matrix::aa_char_to_ncbistdaa(c))
        .collect();
    let sbp = crate::algorithm::pssm_search::in_pssm::score_block(
        crate::config::ScoringMatrix::Blosum62,
        crate::core::psi::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
        },
    )
    .unwrap();
    assert_eq!(sbp.ideal_k, 0.13395614448848181);
    let info = PosMatrixInfo::new(
        &query,
        &search.freq_ratios,
        &sbp,
        32.0,
        file_ungapped(&search),
    )
    .unwrap();
    assert_eq!(info.ungapped_lambda, 0.32206865605351226 / 32.0);
    assert_eq!(info.gapped_k, 0.042084934723385127);
    assert_eq!(info.start_matrix.len(), 398);
}

/// The ungapped block of a PSSM file (all three values present).
fn file_ungapped(search: &crate::algorithm::pssm_search::PsiSearch) -> crate::core::psi::KarlinBlk {
    crate::core::psi::KarlinBlk {
        lambda: search.stats.lambda_ungapped.unwrap(),
        k: search.stats.kappa_ungapped.unwrap(),
        h: search.stats.h_ungapped.unwrap(),
    }
}

// Review S1 of PS4b: the tk PSSM with every frequency ratio 5e-5 (`ps4_tk_eps`) has
// factor-1 scores of -24 to -17 for the 20 letters, so NCBI's `Blast_ScoreFreqNew` gives
// NULL, `Blast_KarlinBlkUngappedCalc` returns 1 and the redo keeps the file's ungapped block:
// an LD_PRELOAD probe of NCBI 2.17.0 (`psiblast -in_pssm tk_eps.asn -db paja`) saw Lambda
// 0.31673860011845001 and K 0.13910107343338302 in the redo and K 0.042574710048176159 in
// every redo e-value.
#[test]
fn a_factor_1_pssm_without_score_frequencies_keeps_the_files_block() {
    let pssm =
        crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip("pssm/ps4_tk_eps.asn.gz"))
            .unwrap();
    let search = crate::algorithm::pssm_search::in_pssm::psi_search_of(pssm).unwrap();
    let query: Vec<u8> = search
        .query
        .seq()
        .iter()
        .map(|&c| crate::utils::matrix::aa_char_to_ncbistdaa(c))
        .collect();
    let sbp = crate::algorithm::pssm_search::in_pssm::score_block(
        crate::config::ScoringMatrix::Blosum62,
        crate::core::psi::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
        },
    )
    .unwrap();
    let kbp_psi = file_ungapped(&search);
    assert_eq!(kbp_psi.lambda, 0.31673860011845001);
    assert_eq!(kbp_psi.k, 0.13910107343338302);
    let info = PosMatrixInfo::new(&query, &search.freq_ratios, &sbp, 32.0, kbp_psi).unwrap();
    assert_eq!(info.ungapped_lambda, 0.31673860011845001 / 32.0);
    assert_eq!(info.gapped_k, 0.042574710048176159);
}

#[test]
fn all_zero_ratios_are_those_at_most_epsilon() {
    assert!(freq_ratios_all_zeros(&[]));
    assert!(freq_ratios_all_zeros(&[[0.0001; BLASTAA_SIZE]]));
    let mut row = [0.0; BLASTAA_SIZE];
    row[7] = 0.00011;
    assert!(!freq_ratios_all_zeros(&[row]));
}
