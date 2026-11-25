use super::*;
use crate::report::serial::pssm::{Pssm, PssmIntermediateData};

fn sample_pssm(freqs: Option<Vec<f64>>, obsr: Option<Vec<f64>>) -> Pssm {
    Pssm {
        is_protein: Some(true),
        identifier: None,
        num_rows: BLASTAA_SIZE as i32,
        num_columns: obsr.as_ref().map_or(1, |v| v.len().max(1) as i32),
        row_labels: None,
        by_row: Some(false),
        query: None,
        intermediate_data: Some(PssmIntermediateData {
            weighted_res_freqs_per_pos: freqs,
            num_indept_obsr: obsr,
            ..PssmIntermediateData::default()
        }),
        final_data: None,
    }
}

fn unit_column(first: f64) -> Vec<f64> {
    let mut col = vec![0.0; BLASTAA_SIZE];
    col[0] = first;
    col
}

#[test]
fn uint4_of_accepts_defined_range_and_rejects_ub() {
    assert_eq!(uint4_of(0.0), Some(0));
    assert_eq!(uint4_of(-0.5), Some(0));
    assert_eq!(uint4_of(5999.999), Some(5999));
    assert_eq!(uint4_of(4_294_967_295.9), Some(u32::MAX));
    assert_eq!(uint4_of(-1.0), None);
    assert_eq!(uint4_of(-1000.0), None);
    assert_eq!(uint4_of(4_294_967_296.0), None);
    assert_eq!(uint4_of(f64::NAN), None);
    assert_eq!(uint4_of(f64::INFINITY), None);
    assert_eq!(uint4_of(f64::NEG_INFINITY), None);
}

#[test]
fn encode_obsr_merges_relative_to_run_head_and_truncates() {
    // 5.0 and 5.00006 merge (diff 6e-5 < 1e-4); 5.00012 starts a new run (diff from 5.0 = 1.2e-4);
    // 5.00022 merges with 5.00012 because IEEE-754 5.00022 - 5.00012 < 1e-4; 5.0003199 starts new.
    let vals = [5.0, 5.00006, 5.00012, 5.00012, 5.00022, 5.0003199];
    assert_eq!(encode_obsr(&vals).unwrap(), [6000, 2, 6000, 3, 6000, 1]);
    assert_eq!(encode_obsr(&[-0.5, -0.5]).unwrap(), [500, 2]);
    assert!(matches!(
        encode_obsr(&[-2.0]),
        Err(DeltaError::Rps(RpsError::Refused(_)))
    ));
    assert!(matches!(
        encode_obsr(&[f64::NAN]),
        Err(DeltaError::Rps(RpsError::Refused(_)))
    ));
}

#[test]
fn check_delta_thresholds_validation_and_error_order() {
    // Exact boundary: max(obsr) + 1.0 == 6.0 -> included; 5.9999 < 6.0 -> excluded.
    let p_eq = sample_pssm(Some(unit_column(1.0)), Some(vec![5.0]));
    assert_eq!(
        check_delta(&p_eq, 1, b"a.smp", 6.0, true).unwrap(),
        DeltaCheck::Included
    );
    let p_low = sample_pssm(Some(unit_column(1.0)), Some(vec![4.9999]));
    assert_eq!(
        check_delta(&p_low, 1, b"a.smp", 6.0, false).unwrap(),
        DeltaCheck::ExcludedObsr
    );

    // Sum epsilon boundaries: 1.0001 and 0.9999 are within 0.0001 in f64; 1.00011 and 0.99989 fail.
    for ok_sum in [1.0001, 0.9999] {
        let p = sample_pssm(Some(unit_column(ok_sum)), Some(vec![10.0]));
        assert_eq!(
            check_delta(&p, 1, b"a.smp", 6.0, true).unwrap(),
            DeltaCheck::Included,
            "{ok_sum}"
        );
    }
    for bad_sum in [1.00011, 0.99989] {
        let p = sample_pssm(Some(unit_column(bad_sum)), Some(vec![10.0]));
        assert_eq!(
            check_delta(&p, 1, b"a.smp", 6.0, true).unwrap(),
            DeltaCheck::ExcludedInvalidCd,
            "{bad_sum}"
        );
        assert_eq!(
            check_delta(&p, 1, b"a.smp", 6.0, false).unwrap(),
            DeltaCheck::Included,
            "{bad_sum}"
        );
    }

    // Observation < 1.0 is invalid even when obsr_threshold is 0.0.
    let p_obsr_lt1 = sample_pssm(Some(unit_column(1.0)), Some(vec![0.9999]));
    assert_eq!(
        check_delta(&p_obsr_lt1, 1, b"a.smp", 0.0, true).unwrap(),
        DeltaCheck::ExcludedInvalidCd
    );
    assert_eq!(
        check_delta(&p_obsr_lt1, 1, b"a.smp", 0.0, false).unwrap(),
        DeltaCheck::Included
    );

    // Column count mismatch throws even when exclude_invalid is false, unless excluded by threshold first.
    let p_mis = sample_pssm(Some(vec![0.0; 56]), Some(vec![10.0]));
    assert!(matches!(
        check_delta(&p_mis, 1, b"a.smp", 6.0, false),
        Err(DeltaError::Exception(msg)) if msg == b"Number of frequency and observations columns do not match"
    ));
    assert_eq!(
        check_delta(&p_mis, 1, b"a.smp", 100.0, false).unwrap(),
        DeltaCheck::ExcludedObsr
    );
}

#[test]
fn check_delta_refuses_ncbi_ub_shapes() {
    // byRow TRUE
    let mut p_byrow = sample_pssm(Some(unit_column(1.0)), Some(vec![10.0]));
    p_byrow.by_row = Some(true);
    assert!(matches!(
        check_delta(&p_byrow, 1, b"a.smp", 6.0, true),
        Err(DeltaError::Rps(RpsError::Refused(_)))
    ));

    // numRows < 28 with non-empty weightedResFreqsPerPos
    let mut p_rows26 = sample_pssm(Some(vec![0.0; 26]), Some(vec![10.0]));
    p_rows26.num_rows = 26;
    assert!(matches!(
        check_delta(&p_rows26, 1, b"a.smp", 6.0, true),
        Err(DeltaError::Rps(RpsError::Refused(_)))
    ));

    // Empty numIndeptObsr
    let p_empty_obsr = sample_pssm(Some(unit_column(1.0)), Some(Vec::new()));
    assert!(matches!(
        check_delta(&p_empty_obsr, 1, b"a.smp", 6.0, true),
        Err(DeltaError::Rps(RpsError::Refused(_)))
    ));

    // Partial trailing column: refused when every full column passes (even if obsr < 1.0, because
    // NCBI's column loop runs before the obsr < 1.0 check).
    let mut partial = unit_column(1.0);
    partial.push(0.0);
    let p_partial = sample_pssm(Some(partial.clone()), Some(vec![10.0]));
    assert!(matches!(
        check_delta(&p_partial, 1, b"a.smp", 6.0, true),
        Err(DeltaError::Rps(RpsError::Refused(_)))
    ));
    let p_partial_low_obsr = sample_pssm(Some(partial), Some(vec![0.5]));
    assert!(matches!(
        check_delta(&p_partial_low_obsr, 1, b"a.smp", 0.0, true),
        Err(DeltaError::Rps(RpsError::Refused(_)))
    ));
}

#[test]
fn validate_cd_partial_column_with_earlier_bad_column_is_not_refused() {
    // 29 zeros with 1 obsr = 10.0: column 0 sums to 0.0, so NCBI returns false before reaching the
    // trailing 1-element partial column (no UB).
    let p_bad0 = sample_pssm(Some(vec![0.0; BLASTAA_SIZE + 1]), Some(vec![10.0]));
    assert_eq!(
        check_delta(&p_bad0, 1, b"a.smp", 6.0, true).unwrap(),
        DeltaCheck::ExcludedInvalidCd
    );
    assert_eq!(
        check_delta(&p_bad0, 1, b"a.smp", 6.0, false).unwrap(),
        DeltaCheck::Included
    );

    // 2 full columns + 1 trailing entry: column 0 passes, column 1 fails -> still returns false
    // without reaching the trailing partial column.
    let mut two_cols_bad1 = unit_column(1.0);
    two_cols_bad1.extend_from_slice(&unit_column(0.5));
    two_cols_bad1.push(0.0);
    let p_bad1 = sample_pssm(Some(two_cols_bad1), Some(vec![10.0, 10.0]));
    assert_eq!(
        check_delta(&p_bad1, 2, b"a.smp", 6.0, true).unwrap(),
        DeltaCheck::ExcludedInvalidCd
    );
    assert_eq!(
        check_delta(&p_bad1, 2, b"a.smp", 6.0, false).unwrap(),
        DeltaCheck::Included
    );
}

#[test]
fn delta_offset_arithmetic_wraps_modulo_2_pow_32_like_ncbi() {
    // Crossing 2^31 wraps Int4 to negative, and freq_record_offset divides by 28 truncating
    // toward zero, matching NCBI's 32-bit offset counter and column division.
    assert_eq!(advance_offset(0, 56), 56);
    assert_eq!(freq_record_offset(56), 2);
    assert_eq!(freq_record_offset(57), 2);

    let wrapped = advance_offset(i32::MAX - 27, 28);
    assert_eq!(wrapped, i32::MIN);
    assert_eq!(freq_record_offset(wrapped), i32::MIN / 28);
    assert_eq!(freq_record_offset(-27), 0);
    assert_eq!(freq_record_offset(-28), -1);

    // Crossing 2^32 wraps back to non-negative.
    assert_eq!(advance_offset(-1, 1), 0);
    assert_eq!(advance_offset(-10, 38), 28);
    assert_eq!(freq_record_offset(advance_offset(-10, 38)), 1);

    // Adding more than 2^32 on 64-bit targets takes the low 32 bits.
    #[cfg(target_pointer_width = "64")]
    assert_eq!(advance_offset(0, (1usize << 32) + 56), 56);
}
