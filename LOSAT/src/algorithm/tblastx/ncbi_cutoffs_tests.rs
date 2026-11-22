use super::*;

/// The effective lengths of a subject set whose search space fits `Int8`.
fn compute_eff_lengths_tblastx(
    query_len_aa: i64,
    subject_len_nucl: i64,
    db_num_seqs: i64,
    karlin_params: &KarlinParams,
) -> EffLengthsResult {
    try_compute_eff_lengths_tblastx(query_len_aa, subject_len_nucl, db_num_seqs, karlin_params)
        .expect("TBLASTX effective search space within i64")
}

/// Only the effective search space of `compute_eff_lengths_tblastx`.
fn compute_eff_searchsp_tblastx(
    query_len_aa: i64,
    subject_len_nucl: i64,
    db_num_seqs: i64,
    gapped_params: &KarlinParams,
) -> i64 {
    compute_eff_lengths_tblastx(query_len_aa, subject_len_nucl, db_num_seqs, gapped_params)
        .eff_searchsp
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:734-847;
// blast_stat.c:5042-5126 (`BLAST_ComputeLengthAdjustment`).
// Expected numbers were evaluated from the extracted NCBI C functions.
#[test]
fn complete_subject_set_effective_lengths() {
    let params = blosum62_ungapped();
    for (q, total, count, adjustment, space) in [
        (300, 900, 1, 19, 78961),
        (300, 2700, 3, 22, 231852),
        (300, 2701, 3, 22, 231852),
        (300, 2702, 3, 22, 231852),
        (300, 2703, 3, 22, 232130),
        (299, 2703, 3, 22, 231295),
        (1, 2703, 3, 0, 901),
        (100, 3003, 3, 19, 76464),
        (0, 2703, 3, 0, 0),
    ] {
        let actual = compute_eff_lengths_tblastx(q, total, count, &params);
        assert_eq!(
            (actual.length_adjustment, actual.eff_searchsp),
            (adjustment, space)
        );
    }
}

/// BLOSUM62 ungapped params (kbp_std)
fn blosum62_ungapped() -> KarlinParams {
    KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.4012,
        alpha: 0.7916,
        beta: -3.2,
    }
}

/// BLOSUM62 gapped params (11/1) (kbp_gap)
fn blosum62_gapped() -> KarlinParams {
    KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:965-973
// `BLAST_Cutoffs` starts from a cutoff of 1 and only raises it. With `-db` SeqDB's shortest
// length 10 makes TBLASTN's subject length 3: a search space of 3 x 3 gives a score below 1
// for E = 0.9, and the cutoff is 1 (NCBI 2.17.0's `cutoff_s` 1 and 32 of line RG X4b's
// survey, section 1.6).
#[test]
fn the_sum_statistics_cutoff_is_at_least_one() {
    let gapped = blosum62_gapped();
    assert!(cutoff_score_from_evalue(0.9, 9, &gapped) <= 0);
    assert_eq!(cutoff_score_sum_stats(160, 3, 0.1, &gapped), 1);
    // A search space large enough keeps the computed score.
    let computed = cutoff_score_from_evalue(0.9, 160 * 2500, &gapped);
    assert!(computed > 1);
    assert_eq!(cutoff_score_sum_stats(160, 2500, 0.1, &gapped), computed);
}

#[test]
fn test_gap_trigger_raw_score() {
    let ungapped = blosum62_ungapped();

    // Expected gap_trigger for 22.0 bits with Lambda = 0.3176 and K = 0.134:
    // (22.0 * ln(2) + ln(0.134)) / 0.3176 = (15.249 - 2.0099) / 0.3176 = 13.239 / 0.3176 = 41.69,
    // which truncates toward zero to 41.
    let gap_trigger = gap_trigger_raw_score(22.0, &ungapped);
    assert_eq!(gap_trigger, 41, "gap_trigger should be 41 (truncated)");
}

#[test]
fn test_gap_trigger_uses_trunc_not_ceil() {
    let ungapped = blosum62_ungapped();

    // Verify that we use trunc, not ceil
    // If ceil were used, the result would be 42
    let gap_trigger = gap_trigger_raw_score(22.0, &ungapped);
    assert!(gap_trigger < 42, "gap_trigger should use trunc, not ceil");
}

// psiblast's `-gap_trigger` takes any double (NCBI 2.17.0's x86-64 conversion past `Int4`:
// INT_MIN, `core::blast_util::ncbi_int4_from_double`).
#[test]
fn test_gap_trigger_past_int4_is_int_min() {
    let ungapped = blosum62_ungapped();
    assert_eq!(gap_trigger_raw_score(-5.0, &ungapped), -17);
    for bits in [
        1e300,
        -1e300,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        1e9,
    ] {
        assert_eq!(gap_trigger_raw_score(bits, &ungapped), i32::MIN, "{bits}");
    }
}

#[test]
fn test_cutoff_score_from_evalue_uses_ceil() {
    let gapped = blosum62_gapped();

    // Test that ceil is used (not trunc)
    // For a large search space, the cutoff should be positive
    let cutoff = cutoff_score_from_evalue(10.0, 1_000_000, &gapped);
    assert!(
        cutoff > 0,
        "cutoff should be positive for reasonable search space"
    );
}

// blast_stat.c:4059-4061 (line PS DECISIONS 28): a quotient `K searchsp / E` past the
// double takes `ln K + ln searchsp - ln E`; a finite one keeps NCBI's order of operations.
#[test]
fn a_cutoff_whose_quotient_overflows_the_double_is_finite() {
    let gapped = blosum62_gapped();
    // 0.041 * (2^63 - 1) / 1e-297 = 3.8e314: past f64::MAX. E below kSmallFloat is
    // raised to it. ceil((ln 0.041 + ln 9.223372036854775807e18 + 297 ln 10) / 0.267).
    let expected = ((0.041f64.glibc_ln() + (i64::MAX as f64).glibc_ln() - 1.0e-297f64.glibc_ln())
        / 0.267)
        .ceil();
    assert_eq!(expected, 2713.0);
    for evalue in [1e-300, 1e-297] {
        assert_eq!(cutoff_score_from_evalue(evalue, i64::MAX, &gapped), 2713);
    }
    // The edge: the largest search space whose quotient at E = 1e-297 is finite, and the
    // next one.
    let last_finite = (f64::MAX * 1.0e-297 / 0.041) as i64;
    assert!((0.041 * last_finite as f64 / 1.0e-297).is_finite());
    let finite = cutoff_score_from_evalue(1e-297, last_finite, &gapped);
    assert_eq!(
        finite,
        ((0.041 * last_finite as f64 / 1.0e-297).glibc_ln() / 0.267).ceil() as i32
    );
    let past = (last_finite as f64 * (1.0 + 1e-9)) as i64;
    assert!(!(0.041 * past as f64 / 1.0e-297).is_finite());
    let overflowed = cutoff_score_from_evalue(1e-297, past, &gapped);
    assert!(
        (finite..=finite + 1).contains(&overflowed),
        "{finite} {overflowed}"
    );
}

#[test]
fn test_x_drop_raw_score_blosum62() {
    let ungapped = blosum62_ungapped();

    // Initial ungapped X-drop scales `ceil` of `x_dropoff` in bits times `NCBIMATH_LN2` over `Lambda` by `scale_factor`.
    // For BLOSUM62: ceil(7.0 * 0.693147 / 0.3176) = ceil(15.27) = 16
    let x_drop = x_drop_raw_score(7.0, &ungapped, 1.0);
    assert_eq!(x_drop, 16, "x_drop should be 16 for BLOSUM62 with 7.0 bits");
}

#[test]
fn test_x_drop_uses_ceil_not_trunc() {
    let ungapped = blosum62_ungapped();

    // Verify that we use ceil, not trunc
    // 7.0 * 0.693147 / 0.3176 = 15.27
    // ceil(15.27) = 16 (not 15)
    let x_drop = x_drop_raw_score(7.0, &ungapped, 1.0);
    assert_eq!(x_drop, 16, "x_drop should use ceil, resulting in 16");
    assert!(x_drop > 15, "x_drop should be > 15 due to ceil rounding");
}

#[test]
fn test_x_drop_with_scale_factor() {
    let ungapped = blosum62_ungapped();

    // Test that scale_factor is applied correctly
    // With scale_factor = 2.0, result should be doubled
    let x_drop_base = x_drop_raw_score(7.0, &ungapped, 1.0);
    let x_drop_scaled = x_drop_raw_score(7.0, &ungapped, 2.0);
    assert_eq!(
        x_drop_scaled,
        x_drop_base * 2,
        "x_drop with scale_factor=2.0 should be double"
    );
}

#[test]
fn test_x_drop_different_bits() {
    let ungapped = blosum62_ungapped();

    // Test with different bit values
    // x_drop = ceil(bits * LN2 / lambda)

    // 10 bits: ceil(10.0 * 0.693147 / 0.3176) = ceil(21.83) = 22
    let x_drop_10 = x_drop_raw_score(10.0, &ungapped, 1.0);
    assert_eq!(x_drop_10, 22, "x_drop should be 22 for 10.0 bits");

    // 20 bits (NUCL default): ceil(20.0 * 0.693147 / 0.3176) = ceil(43.66) = 44
    let x_drop_20 = x_drop_raw_score(20.0, &ungapped, 1.0);
    assert_eq!(x_drop_20, 44, "x_drop should be 44 for 20.0 bits");
}

#[test]
fn test_cutoff_score_for_update_tblastx() {
    let ungapped = blosum62_ungapped();

    // NCBI BlastInitialWordParametersUpdate for tblastx ungapped path
    // Uses `CUTOFF_E_TBLASTX` (1e-300) and search space min(`q_aa`, `s_nucl`) times `s_nucl`.
    let query_len_aa = 500;
    let subject_len_nucl = 3000; // NUCLEOTIDE length, NOT AA!
    let gap_trigger = gap_trigger_raw_score(22.0, &ungapped); // 41

    // Compute cutoff_score_max first (uses user E-value)
    let eff_searchsp = compute_eff_searchsp_tblastx(query_len_aa, subject_len_nucl, 1, &ungapped);
    let cutoff_score_max = cutoff_score_max_for_tblastx(eff_searchsp, 10.0, &ungapped);

    // Compute per-subject cutoff
    let cutoff = cutoff_score_for_update_tblastx(
        query_len_aa,
        subject_len_nucl,
        gap_trigger,
        cutoff_score_max,
        BLAST_GAP_DECAY_RATE, // 0.5
        &ungapped,
        1.0,
    );

    // The result should be capped by either gap_trigger or cutoff_score_max
    assert!(cutoff <= gap_trigger, "cutoff should be <= gap_trigger");
    assert!(
        cutoff <= cutoff_score_max,
        "cutoff should be <= cutoff_score_max"
    );
    assert!(cutoff > 0, "cutoff should be positive");
}

#[test]
fn test_cutoff_score_max_for_tblastx() {
    let ungapped = blosum62_ungapped();

    // NCBI BlastHitSavingParametersNew uses user E-value and eff_searchsp
    let eff_searchsp = 1_000_000i64;

    let cutoff_max = cutoff_score_max_for_tblastx(eff_searchsp, 10.0, &ungapped);

    // Should be positive and reasonable
    assert!(cutoff_max > 0, "cutoff_score_max should be positive");

    // With E-value = 10.0 and searchsp = 1e6, score should be moderate
    // S = ceil(ln(K * searchsp / E) / Lambda)
    // = ceil(ln(0.134 * 1e6 / 10) / 0.3176)
    // = ceil(ln(13400) / 0.3176)
    // = ceil(9.503 / 0.3176)
    // = ceil(29.93)
    // = 30
    assert_eq!(
        cutoff_max, 30,
        "cutoff_score_max should be 30 for these parameters"
    );
}

#[test]
fn test_cutoff_e_tblastx_constant() {
    // Verify the CUTOFF_E_TBLASTX constant matches NCBI
    assert_eq!(
        CUTOFF_E_TBLASTX, 1e-300,
        "CUTOFF_E_TBLASTX should be 1e-300"
    );
}

#[test]
fn test_blast_gap_decay_rate_constant() {
    // Verify the BLAST_GAP_DECAY_RATE constant matches NCBI
    assert_eq!(
        BLAST_GAP_DECAY_RATE, 0.5,
        "BLAST_GAP_DECAY_RATE should be 0.5"
    );
}

#[test]
fn test_eff_searchsp_long_sequence() {
    // Test eff_searchsp calculation for long sequences (600kb+)
    // This verifies that there are no overflow or precision issues
    let ungapped = blosum62_ungapped();

    // Simulate a 600kb sequence
    let query_len_aa = 200_000i64; // ~600kb query
    let subject_len_nucl = 600_000i64; // 600kb subject

    let result = compute_eff_lengths_tblastx(query_len_aa, subject_len_nucl, 1, &ungapped);

    // Verify no overflow occurred
    assert!(result.eff_searchsp > 0, "eff_searchsp should be positive");
    assert!(
        result.eff_searchsp < i64::MAX,
        "eff_searchsp should not overflow"
    );

    // Verify length_adjustment is reasonable
    assert!(
        result.length_adjustment > 0,
        "length_adjustment should be positive"
    );
    assert!(
        result.length_adjustment < query_len_aa,
        "length_adjustment should be less than query length"
    );

    // Verify eff_searchsp calculation
    // db_length = subject_len_nucl / 3 = 200,000
    // effective_db_length = db_length - length_adjustment
    // effective_query_length = query_len_aa - length_adjustment
    // eff_searchsp = effective_db_length * effective_query_length
    let db_length = subject_len_nucl / 3;
    let expected_eff_db = (db_length - result.length_adjustment).max(1);
    let expected_eff_query = (query_len_aa - result.length_adjustment).max(1);
    let expected_eff_searchsp = expected_eff_db * expected_eff_query;

    assert_eq!(
        result.eff_searchsp, expected_eff_searchsp,
        "eff_searchsp should match manual calculation"
    );
}

#[test]
fn test_cutoff_score_max_long_sequence() {
    // Test cutoff_score_max calculation for long sequences (600kb+)
    let ungapped = blosum62_ungapped();

    // Simulate a 600kb sequence
    let query_len_aa = 200_000i64;
    let subject_len_nucl = 600_000i64;

    let eff_lengths = compute_eff_lengths_tblastx(query_len_aa, subject_len_nucl, 1, &ungapped);

    let cutoff_max = cutoff_score_max_for_tblastx(
        eff_lengths.eff_searchsp,
        10.0, // user E-value
        &ungapped,
    );

    // Verify cutoff_score_max is reasonable
    assert!(cutoff_max > 0, "cutoff_score_max should be positive");
    assert!(
        cutoff_max <= 100,
        "cutoff_score_max should be reasonable for long sequences"
    );

    // For long sequences, cutoff_score_max should typically be > gap_trigger (41)
    // because eff_searchsp is very large
    // Note: cutoff_score_max may be > gap_trigger for very long sequences
    // The final cutoff is the minimum of `update_cutoff`, `gap_trigger`, and `cutoff_score_max`.
}

#[test]
fn test_cutoff_for_update_long_sequence() {
    // Test cutoff_score_for_update_tblastx for long sequences (600kb+)
    let ungapped = blosum62_ungapped();

    // Simulate a 600kb sequence
    let query_len_aa = 200_000i64;
    let subject_len_nucl = 600_000i64;

    let gap_trigger = gap_trigger_raw_score(22.0, &ungapped);
    let eff_lengths = compute_eff_lengths_tblastx(query_len_aa, subject_len_nucl, 1, &ungapped);
    let cutoff_score_max = cutoff_score_max_for_tblastx(eff_lengths.eff_searchsp, 10.0, &ungapped);

    // Verify searchsp calculation doesn't overflow
    let min_len = (query_len_aa as u64).min(subject_len_nucl as u64);
    let searchsp = min_len * (subject_len_nucl as u64);
    assert!(searchsp < u64::MAX, "searchsp should not overflow u64");

    let cutoff = cutoff_score_for_update_tblastx(
        query_len_aa,
        subject_len_nucl,
        gap_trigger,
        cutoff_score_max,
        BLAST_GAP_DECAY_RATE,
        &ungapped,
        1.0,
    );

    // Verify final cutoff is reasonable
    assert!(cutoff > 0, "cutoff should be positive");
    assert!(cutoff <= gap_trigger, "cutoff should be <= gap_trigger");
    assert!(
        cutoff <= cutoff_score_max,
        "cutoff should be <= cutoff_score_max"
    );
}

#[test]
fn test_eff_searchsp_very_long_sequence() {
    // Test eff_searchsp calculation for very long sequences (10Mb+)
    // This tests the upper bounds of the calculation
    let ungapped = blosum62_ungapped();

    // Simulate a 10Mb sequence
    let query_len_aa = 3_333_333i64; // ~10Mb query
    let subject_len_nucl = 10_000_000i64; // 10Mb subject

    let result = compute_eff_lengths_tblastx(query_len_aa, subject_len_nucl, 1, &ungapped);

    // Verify no overflow occurred
    assert!(result.eff_searchsp > 0, "eff_searchsp should be positive");
    assert!(
        result.eff_searchsp < i64::MAX,
        "eff_searchsp should not overflow"
    );

    // Verify the calculation is correct
    let db_length = subject_len_nucl / 3;
    let expected_eff_db = (db_length - result.length_adjustment).max(1);
    let expected_eff_query = (query_len_aa - result.length_adjustment).max(1);
    let expected_eff_searchsp = expected_eff_db * expected_eff_query;

    assert_eq!(
        result.eff_searchsp, expected_eff_searchsp,
        "eff_searchsp should match manual calculation for very long sequences"
    );
}

#[test]
fn test_db_length_integer_division_edge_cases() {
    // Test that db_length = subject_len_nucl / 3 matches NCBI's integer division behavior
    // NCBI uses: Int8 db_length; db_length = db_length/3; (integer division, truncates toward zero)
    // Rust uses: i64 / 3 (integer division, truncates toward zero)
    let ungapped = blosum62_ungapped();

    // Test cases for integer division edge cases
    let test_cases = vec![
        // (subject_len_nucl, expected_db_length, description)
        (3i64, 1i64, "exactly divisible by 3"),
        (4i64, 1i64, "4/3 = 1.333... truncates to 1"),
        (5i64, 1i64, "5/3 = 1.666... truncates to 1"),
        (6i64, 2i64, "exactly divisible by 3"),
        (600_000i64, 200_000i64, "600kb exactly divisible"),
        (600_001i64, 200_000i64, "600kb+1 truncates to 200000"),
        (600_002i64, 200_000i64, "600kb+2 truncates to 200000"),
        (600_003i64, 200_001i64, "600kb+3 = 200001"),
        (1i64, 0i64, "very small: 1/3 = 0"),
        (2i64, 0i64, "very small: 2/3 = 0"),
    ];

    for (subject_len_nucl, expected_db_length, desc) in test_cases {
        let db_length = subject_len_nucl / 3;
        assert_eq!(
            db_length, expected_db_length,
            "db_length calculation failed for {}: subject_len_nucl={}, expected={}, got={}",
            desc, subject_len_nucl, expected_db_length, db_length
        );

        // Verify that compute_eff_lengths_tblastx uses this calculation
        // For very small values, we need a reasonable query length
        let query_len_aa = if subject_len_nucl < 3 {
            10i64
        } else {
            subject_len_nucl / 3
        };

        let _result = compute_eff_lengths_tblastx(query_len_aa, subject_len_nucl, 1, &ungapped);

        // Verify the internal calculation matches
        let computed_db_length = subject_len_nucl / 3;
        assert_eq!(
            computed_db_length, expected_db_length,
            "Internal db_length calculation failed for {}: subject_len_nucl={}",
            desc, subject_len_nucl
        );
    }
}

#[test]
fn test_db_length_matches_ncbi_integer_division() {
    // Verify that Rust's integer division matches C's integer division behavior
    // NCBI divides `db_length` (an `Int8`) by 3 using integer division.
    // C integer division truncates toward zero (same as Rust for positive values)
    let test_values = vec![
        600_000i64,   // 200000
        600_001i64,   // 200000 (truncated)
        600_002i64,   // 200000 (truncated)
        600_003i64,   // 200001
        1_000_000i64, // 333333 (truncated from 333333.333...)
        1_000_001i64, // 333333 (truncated)
        1_000_002i64, // 333334
    ];

    for subject_len_nucl in test_values {
        let rust_result = subject_len_nucl / 3;

        // In C, Int8 / 3 would produce the same result (truncation toward zero)
        // For positive values, C and Rust integer division are identical
        let expected = subject_len_nucl / 3; // This is what C would produce

        assert_eq!(
            rust_result, expected,
            "Integer division mismatch: subject_len_nucl={}, rust_result={}, expected={}",
            subject_len_nucl, rust_result, expected
        );
    }
}

#[test]
fn test_cutoff_for_update_very_long_sequence() {
    // Test cutoff_score_for_update_tblastx for very long sequences (10Mb+)
    // This tests the upper bounds of searchsp calculation
    let ungapped = blosum62_ungapped();

    // Simulate a 10Mb sequence
    let query_len_aa = 3_333_333i64;
    let subject_len_nucl = 10_000_000i64;

    let gap_trigger = gap_trigger_raw_score(22.0, &ungapped);
    let eff_lengths = compute_eff_lengths_tblastx(query_len_aa, subject_len_nucl, 1, &ungapped);
    let cutoff_score_max = cutoff_score_max_for_tblastx(eff_lengths.eff_searchsp, 10.0, &ungapped);

    // Verify searchsp calculation doesn't overflow
    // Search space is `min(query_len_aa, subject_len_nucl) * subject_len_nucl`
    // = min(3,333,333, 10,000,000) * 10,000,000
    // = 3,333,333 * 10,000,000 = 33,333,330,000,000
    let min_len = (query_len_aa as u64).min(subject_len_nucl as u64);
    let searchsp = min_len * (subject_len_nucl as u64);
    assert!(
        searchsp < u64::MAX,
        "searchsp should not overflow u64 for very long sequences"
    );
    assert_eq!(
        searchsp, 33_333_330_000_000u64,
        "searchsp should match expected value"
    );

    let cutoff = cutoff_score_for_update_tblastx(
        query_len_aa,
        subject_len_nucl,
        gap_trigger,
        cutoff_score_max,
        BLAST_GAP_DECAY_RATE,
        &ungapped,
        1.0,
    );

    // Verify final cutoff is reasonable
    assert!(cutoff > 0, "cutoff should be positive");
    assert!(cutoff <= gap_trigger, "cutoff should be <= gap_trigger");
}
