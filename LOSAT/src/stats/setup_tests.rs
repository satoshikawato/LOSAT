use super::*;

// NCBI BLAST+ 2.17.0 probes of line RG D1a-3 (p5_ns_*, p5_avg_*, p5_min_*): which setups
// fail with `The average subject length is too short`.
#[test]
fn too_short_is_the_uint4_of_ncbis_quotient() {
    let tfn = 10_937;
    assert_eq!(setup_subject_length(tfn, 16, None), Some(683));
    // NSEQ -1, -5, 2^32 - 1: negative quotients pass.
    assert_eq!(setup_subject_length(tfn, -1, None), Some(-10_937));
    assert_eq!(setup_subject_length(tfn, -5, None), Some(-2_187));
    // NSEQ 2^31 - 1, 2^31, 2^31 + 1, 3e9: a quotient of 0.
    for n in [
        i32::MAX,
        i32::MIN,
        (2_147_483_649u32) as i32,
        (3_000_000_000u32) as i32,
    ] {
        assert_eq!(setup_subject_length(tfn, n, None), None, "{n}");
    }
    // LENGTH 16 * 2^32 over 16 sequences: the `Int4` of 2^32 is 0; 2^31 and 2^32 + 5 run.
    assert_eq!(setup_subject_length(16 << 32, 16, None), None);
    assert_eq!(setup_subject_length(16 << 31, 16, None), Some(i32::MIN));
    assert_eq!(
        setup_subject_length(16 * ((1 << 32) + 5), 16, None),
        Some(5)
    );
    // blastn LENGTH 4294967296000000000 over 16: 268435456000000000 = 2^32 * 62500000.
    assert_eq!(
        setup_subject_length(4_294_967_296_000_000_000, 16, None),
        None
    );
    // The Gumbel block's MIN_SEQ_LENGTH: 1 and 2 fail for a translated subject, 3 runs.
    assert_eq!(setup_subject_length(tfn, 16, Some((1, true))), None);
    assert_eq!(setup_subject_length(tfn, 16, Some((2, true))), None);
    assert_eq!(setup_subject_length(tfn, 16, Some((3, true))), Some(1));
    assert_eq!(setup_subject_length(tfn, -1, Some((10, false))), Some(10));
    // A negative `int` is a large `Uint4`: 2^32 - 1 over 3.
    assert_eq!(
        setup_subject_length(tfn, 16, Some((-1, true))),
        Some(1_431_655_765)
    );
}

// probe p5_ns_2p32m1#10 of line RG D1a-3 (NSEQ -1 over tfn): the adjustment grows with the
// length the negative count adds; `-searchsp` keeps the option's space.
#[test]
fn effective_lengths_take_ncbis_int4_count() {
    let ka = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    let positive = calc_eff_lengths(142, 3263, 23, &ka, 0, "BLASTP").unwrap();
    let negative = calc_eff_lengths(142, 3263, -1, &ka, 0, "BLASTP").unwrap();
    let adj = i64::from(negative.length_adjustment);
    assert_eq!(negative.eff_searchsp, (3263 + adj) * (142 - adj));
    assert_eq!(
        positive.eff_searchsp,
        (3263 - 23 * i64::from(positive.length_adjustment))
            * (142 - i64::from(positive.length_adjustment))
    );
    let fixed = calc_eff_lengths(142, 3263, -1, &ka, 1_000_000, "BLASTP").unwrap();
    assert_eq!(
        (fixed.length_adjustment, fixed.eff_searchsp),
        (negative.length_adjustment, 1_000_000)
    );
    let overflow = calc_eff_lengths(142, i64::MAX, -1, &ka, 0, "BLASTP").unwrap_err();
    assert!(overflow
        .to_string()
        .contains("64-bit effective search space"));
    assert_eq!(eff_lengths_database(0, 1, 0, false, false), None);
    assert_eq!(eff_lengths_database(0, 1, 0, true, true), Some((0, 1)));
    assert_eq!(eff_lengths_database(30, 2, 0, false, true), Some((10, 2)));
    assert_eq!(eff_lengths_database(30, 2, 99, false, false), Some((99, 2)));
}

// SURVEY_D1a3 §0.3 / impl2 of line RG D1a-3: the preliminary stage's option spaces
// (`SplitQuery_SetEffectiveSearchSpace`) and the statistics of one subject.
#[test]
fn preliminary_option_spaces_follow_ncbis_split_query_setup() {
    let ka = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    // -searchsp is the options' own.
    assert_eq!(
        preliminary_option_space(142, 0, 23, 0, 777, false, &ka, "BLASTP").unwrap(),
        777
    );
    // A total of 0 without -dbsize: CalcEffLengths returns before the contexts.
    assert_eq!(
        preliminary_option_space(142, 0, 23, 0, 0, false, &ka, "BLASTP").unwrap(),
        0
    );
    // With -dbsize the walk's count, not the setup's fallback of 1.
    let dbsize = preliminary_option_space(142, 0, 23, 100_000, 0, false, &ka, "BLASTP").unwrap();
    assert_eq!(
        dbsize,
        calc_eff_lengths(142, 100_000, 23, &ka, 0, "BLASTP")
            .unwrap()
            .eff_searchsp
    );
    // A negative total: an effective length of 1, the context's query length.
    assert_eq!(
        preliminary_option_space(
            142,
            -9_223_372_036_854_775_807,
            23,
            0,
            0,
            false,
            &ka,
            "BLASTP"
        )
        .unwrap(),
        142
    );
    // A positive total: the setup's own space; an empty context has none.
    assert_eq!(
        preliminary_option_space(142, 3263, 23, 0, 0, false, &ka, "BLASTP").unwrap(),
        calc_eff_lengths(142, 3263, 23, &ka, 0, "BLASTP")
            .unwrap()
            .eff_searchsp
    );
    assert_eq!(
        preliminary_option_space(0, 3263, 23, 0, 0, false, &ka, "BLASTP").unwrap(),
        0
    );
    // One subject: its length (/3 translated) with the setup's count, an option space kept.
    let subject = subject_eff_lengths(142, 300, 1, 0, 0, true, &ka, "TBLASTN").unwrap();
    assert_eq!(
        subject,
        calc_eff_lengths(142, 100, 1, &ka, 0, "TBLASTN").unwrap()
    );
    let kept = subject_eff_lengths(142, 300, 7, 0, 481_440, false, &ka, "BLASTP").unwrap();
    assert_eq!(kept.eff_searchsp, 481_440);
    let fixed = subject_eff_lengths(142, 300, 1, 100_000, 0, false, &ka, "BLASTP").unwrap();
    assert_eq!(
        fixed,
        calc_eff_lengths(142, 100_000, 1, &ka, 0, "BLASTP").unwrap()
    );
}
