use super::*;

// C's `%e`: six digits after the point, a two-digit exponent at least, glibc's words for the
// infinities and NaN (survey case S_kappa).
#[test]
fn kappa_is_written_as_c_scientific() {
    assert_eq!(c_scientific(0.0), "0.000000e+00");
    assert_eq!(c_scientific(1.2345675e200), "1.234568e+200");
    assert_eq!(c_scientific(-2.5), "-2.500000e+00");
    assert_eq!(c_scientific(0.12345675), "1.234568e-01");
    assert_eq!(c_scientific(f64::MIN_POSITIVE), "2.225074e-308");
    assert_eq!(c_scientific(f64::INFINITY), "inf");
    assert_eq!(c_scientific(f64::NEG_INFINITY), "-inf");
}

// `(Int4) BLAST_Nint(x)`: halves away from zero, then the low 32 bits of the `long`; no value
// where the conversion to `long` is undefined (survey cases S_fredge, S_frhuge).
#[test]
fn frequency_ratios_round_half_away_and_wrap_to_32_bits() {
    assert_eq!(nint_int4(0.5 * FREQ_RATIO_SCALE), Some(500_000_000));
    assert_eq!(nint_int4(-1.5 * FREQ_RATIO_SCALE), Some(-1_500_000_000));
    assert_eq!(nint_int4(2.5), Some(3));
    assert_eq!(nint_int4(-2.5), Some(-3));
    assert_eq!(nint_int4(7e5 * FREQ_RATIO_SCALE), Some(1_935_130_624));
    assert_eq!(nint_int4(1e300 * FREQ_RATIO_SCALE), None);
    assert_eq!(nint_int4(f64::NAN), None);
}

// `(Int4)` of the scaled word threshold, kept in 64 bits: truncation toward zero, the mathematical
// value past `Int4` (survey cases V25, S_cwstbig), a NaN refused.
#[test]
fn word_thresholds_truncate_and_keep_their_value_past_int4() {
    assert_eq!(word_threshold(982.5).unwrap(), 982);
    assert_eq!(word_threshold(9.82).unwrap(), 9);
    assert_eq!(word_threshold(-5.5).unwrap(), -5);
    assert_eq!(word_threshold(3e9).unwrap(), 3_000_000_000);
    assert_eq!(word_threshold(1e14).unwrap(), 100_000_000_000_000);
    assert!(matches!(
        word_threshold(f64::NAN),
        Err(RpsError::Refused(_))
    ));
}

// The scores by position (`numRows` per position, the rest `INT2_MIN`) or by row, the extra row
// of `-BLAST_SCORE_MAX`, and NCBI's two errors of a list of the wrong length.
#[test]
fn pssm_rows_follow_the_profile_layout() {
    let by_position: Vec<i32> = (0..2 * 26).collect();
    let rows = pssm_rows(&by_position, 26, false, 2).unwrap();
    assert_eq!(rows.rows.len(), 3);
    assert_eq!(rows.rows[1][0], 26);
    assert_eq!(rows.rows[1][25], 51);
    assert_eq!(rows.rows[1][26], INT2_MIN);
    assert_eq!(rows.rows[2], [-BLAST_SCORE_MAX; BLASTAA_SIZE]);

    let by_row: Vec<i32> = (0..2 * 28).collect();
    let rows = pssm_rows(&by_row, 28, true, 2).unwrap();
    assert_eq!(rows.rows[0][1], 2);
    assert_eq!(rows.rows[1][0], 1);

    assert!(matches!(
        pssm_rows(&by_position[..51], 26, false, 2),
        Err(RpsError::Blast("PSSM was truncated early"))
    ));
    assert!(matches!(
        pssm_rows(&[0; 53], 26, false, 2),
        Err(RpsError::Blast("PSSM too large for this sequence"))
    ));
}

// `GetFreqRatios`: 28 rows of `numColumns` ratios stored row by row, filled by row or by column;
// the `runtime_error` of an empty list; a short list refused.
#[test]
fn freq_ratio_matrix_is_stored_row_by_row() {
    let ratios: Vec<f64> = (0..2 * 28).map(f64::from).collect();
    let by_column = freq_ratio_matrix(&ratios, 28, 2, false).unwrap();
    assert_eq!(by_column.len(), 56);
    assert_eq!(by_column[1], 28.0);
    assert_eq!(by_column[2], 1.0);
    let by_row = freq_ratio_matrix(&ratios, 28, 2, true).unwrap();
    assert_eq!(by_row[1], 1.0);
    assert!(matches!(
        freq_ratio_matrix(&[], 28, 2, false),
        Err(RpsError::Unknown)
    ));
    assert!(matches!(
        freq_ratio_matrix(&ratios[..55], 28, 2, false),
        Err(RpsError::Refused(_))
    ));
}
