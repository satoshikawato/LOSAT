use super::super::stage_d_stats::LocalLinkParameters;
use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2369-2382,
// 2439-2448:
// for (i=0; i<num_queries; i++) if (kbp_gap[i] != NULL &&
//     if (contexts[index].is_valid) { near_identical_cutoff =
//         (1.74*NCBIMATH_LN2)/kbp_gap[index]->Lambda; break; }
#[test]
fn first_valid_lambda_differs_from_minimum_redo_lambda() {
    let gapped = [
        KarlinParams {
            lambda: 0.3,
            ..KarlinParams::default()
        },
        KarlinParams {
            lambda: 0.267,
            ..KarlinParams::default()
        },
        KarlinParams {
            lambda: 0.2,
            ..KarlinParams::default()
        },
    ];
    let initial = LocalSubjectParameters {
        lengths: Vec::new(),
        cutoffs: Vec::new(),
        link: Some(LocalLinkParameters {
            gap_decay_rate: 0.1,
            gap_size: 40,
            overlap_size: 9,
            longest_intron: 40,
            cutoff_small_gap: 0,
        }),
        hit_cutoff_min: 25,
        word_cutoff_min: 25,
        prelim_evalue: 50.0,
    };
    let params = local_kappa_redo_params(
        ScoringMatrix::Blosum62,
        11,
        1,
        &gapped,
        &[false, true, true],
        &initial,
        120,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        false,
        30.0,
        70,
    )
    .unwrap();
    assert_eq!(
        params.near_identical_cutoff.to_bits(),
        ((1.74 * std::f64::consts::LN_2) / (0.267 / 32.0)).to_bits()
    );
    assert_eq!(
        params.gapping_params.x_dropoff,
        (((30.0 * std::f64::consts::LN_2) / (0.2 / 32.0)) as i32).max(70)
    );
    assert_eq!(params.cutoff_score, 800);
    assert!(params.do_link_hsps);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2428,2450-2457: without linking parameters (`-sum_stats false`,
    // or a `-max_intron_length` of 1 to 4) the redo does not link and its cutoff is 1.
    let params = local_kappa_redo_params(
        ScoringMatrix::Blosum62,
        11,
        1,
        &gapped,
        &[false, true, true],
        &LocalSubjectParameters {
            link: None,
            ..initial
        },
        120,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        false,
        30.0,
        70,
    )
    .unwrap();
    assert_eq!(params.cutoff_score, 1);
    assert!(!params.do_link_hsps);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:457-463:
#[test]
fn extension_final_xdrop_uses_option_bits_and_previous_cutoff() {
    let lambda = 0.267;
    let prelim = ((35.0 * std::f64::consts::LN_2) / lambda) as i32;
    assert_eq!(
        local_extension_final_xdrop(35.0, 25.0, lambda).unwrap(),
        prelim
    );
}
