use super::{RpsUnit, RpsUnits};

/// A unit of `total` letters (before the fallback) whose setup reads `setup` (`T`, `n`).
fn unit(total: i64, setup: (i64, i32), request_count: i32, engine_size: i128) -> RpsUnit {
    RpsUnit {
        total,
        setup_total: setup.0,
        setup_count: setup.1,
        request_count,
        min_length: 10,
        engine_size,
        sorts_results: total > 0,
    }
}

fn units(units: Vec<RpsUnit>, adjusted: Option<(i64, i32)>, dbsize: i64) -> RpsUnits {
    let (db_length, db_num_seqs) = adjusted.unwrap_or((dbsize, 0));
    let order = (0..units.len()).collect();
    RpsUnits {
        units,
        adjusted,
        db_length,
        db_num_seqs,
        order,
    }
}

fn blosum62_gapped() -> crate::stats::KarlinParams {
    crate::stats::KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    }
}

/// Model M2 of line RG D1a-4: `D` and `M` above 0 replace each unit's length and count, at
/// 0 or below each unit reads its own (and, at 0, its own Gumbel length); `-dbsize` replaces
/// the lengths, never the counts.
#[test]
fn each_unit_reads_its_own_numbers_under_the_options() {
    let volume = unit(190, (190, 1), 1, 191);
    for (adjusted, eff, gumbel) in [
        (Some((560, 3)), (560, 3), 560),
        (Some((5000, 7)), (5000, 7), 5000),
        (Some((0, 0)), (190, 1), 190),
        (Some((0, -5)), (190, 1), 190),
        (
            Some((-9_223_372_036_854_775_807, 3)),
            (190, 3),
            -9_223_372_036_854_775_807,
        ),
    ] {
        let all = units(vec![volume; 3], adjusted, 0);
        assert_eq!(all.eff_lengths(&volume), eff, "{adjusted:?}");
        assert_eq!(all.gumbel_length(&volume), gumbel, "{adjusted:?}");
    }
    let dbsize = units(vec![volume; 3], None, 100_000);
    assert_eq!(dbsize.eff_lengths(&volume), (100_000, 1));
    assert_eq!(dbsize.gumbel_length(&volume), 100_000);
    let negative = units(vec![volume; 3], None, -5);
    assert_eq!(negative.eff_lengths(&volume), (190, 1));
    assert_eq!(negative.gumbel_length(&volume), -5);
    // One volume, LENGTH 0 and NSEQ 600: the setup's fallback to OID 0 and one sequence.
    let fallback = unit(0, (190, 1), 600, 600);
    let one = units(vec![fallback], None, 0);
    assert_eq!(one.eff_lengths(&fallback), (190, 1));
    assert_eq!(one.gumbel_length(&fallback), 0);
}

/// Model M6: the request's search space is the last unit's, with its total and count before
/// the fallback: 0 for a total of 0, the query length for a negative one.
#[test]
fn the_request_space_is_the_last_units_before_the_fallback() {
    let gapped = blosum62_gapped();
    let first = unit(190, (190, 1), 1, 191);
    let last = unit(180, (180, 1), 1, 181);
    let two = units(vec![first, last], Some((0, 0)), 0);
    assert_eq!(
        two.request_space(208, &gapped, "RPSBLAST").unwrap(),
        two.search_space(&last, 208, &gapped, 0, "RPSBLAST")
            .unwrap()
    );
    assert_ne!(
        two.request_space(208, &gapped, "RPSBLAST").unwrap(),
        two.search_space(&first, 208, &gapped, 0, "RPSBLAST")
            .unwrap()
    );
    let zero = units(vec![unit(0, (190, 1), 600, 600)], None, 0);
    assert_eq!(zero.request_space(208, &gapped, "RPSBLAST").unwrap(), 0);
    let negative = units(
        vec![unit(-9_223_372_036_854_775_807, (190, 1), 600, 601)],
        None,
        0,
    );
    assert_eq!(
        negative.request_space(208, &gapped, "RPSBLAST").unwrap(),
        208
    );
    assert_eq!(negative.request_space(0, &gapped, "RPSBLAST").unwrap(), 0);
}

// prelim_stage.cpp:219-221 (cited at `RpsUnits::prelim_lengths`): the preliminary stage reads
// SplitQuery's space, whose lengths are the setup's except for a negative total (review A of
// D1a-4, MF1): `LENGTH` 2^63 + 1 and `NSEQ` 600 over one volume, SplitQuery's 208 (1 times the
// query's 208 letters less no adjustment) where the setup reads OID 0's 190 letters and 1.
#[test]
fn the_preliminary_lengths_differ_from_the_setups_for_a_negative_total_only() {
    let negative = unit(-9_223_372_036_854_775_807, (190, 1), 600, 601);
    let one = |unit: RpsUnit| units(vec![unit], None, 0);
    assert_eq!(
        one(negative).prelim_lengths(&negative, 0),
        Some((-9_223_372_036_854_775_807, 600))
    );
    // `-searchsp`: both stages read it.
    assert_eq!(one(negative).prelim_lengths(&negative, 1_000_000), None);
    // `LENGTH 0`: SplitQuery sets no space, the setup computes its own.
    let zero = unit(0, (190, 1), 600, 600);
    assert_eq!(one(zero).prelim_lengths(&zero, 0), None);
    // A positive total: the same lengths (also `STATS_NSEQ`'s count, and under D and M).
    let plain = unit(560, (560, 7), 7, 563);
    assert_eq!(one(plain).prelim_lengths(&plain, 0), None);
    let volume = unit(190, (190, 1), 1, 191);
    let several = units(vec![volume, volume], Some((-5, -3)), 0);
    assert_eq!(several.prelim_lengths(&volume, 0), None);
    let several = units(vec![volume, volume], Some((5000, 7)), 0);
    assert_eq!(several.prelim_lengths(&volume, 0), None);
}

// blast_setup.c:826-846 (cited at `stats::setup::calc_eff_lengths`): with `-dbsize` each
// volume's space reads the given length and its own sequence count; an effective database
// length of 0 or below is 1; a product past `Int8` is refused.
#[test]
fn a_units_search_space_reads_its_own_sequence_count() {
    use crate::stats::length_adjustment::compute_length_adjustment_ncbi;
    let gapped = blosum62_gapped();
    for (query, db_length, seqs) in [(302i64, 5_000i64, 1i32), (302, 5_000, 3), (40, 100_000, 7)] {
        let volume = unit(190, (190, seqs), seqs, 191);
        let adjustment = compute_length_adjustment_ncbi(query, db_length, i64::from(seqs), &gapped)
            .length_adjustment;
        assert_eq!(
            units(vec![volume; 3], None, db_length)
                .search_space(&volume, query, &gapped, 0, "RPSBLAST")
                .unwrap(),
            (db_length - i64::from(seqs) * adjustment) * (query - adjustment),
            "{query} {db_length} {seqs}"
        );
    }
    let volume = unit(190, (190, 1), 1, 191);
    assert_eq!(
        units(vec![volume], None, 0)
            .search_space(&volume, 302, &gapped, 1_000_000, "RPSBLAST")
            .unwrap(),
        1_000_000
    );
    assert!(units(vec![volume; 3], None, i64::MAX)
        .search_space(&volume, 302, &gapped, 0, "RPSBLAST")
        .unwrap_err()
        .to_string()
        .contains("overflows NCBI BLAST+'s 64-bit effective search space"));
}

/// The strategy's `DbLength`: `D` of any sign where `x_AdjustDbSize` runs, -2^63 too, which the
/// requests write as its digits where NCBI writes `-0` (PD-LOSAT-NCBI-DEFECTS approved
/// exception 17, Int8 -2^63 written as digits; NCBI's `big-integer -0` of probe r2_3v_l2p63#2 of line RG D1a-4's review A;
/// the request's bytes: `report::archive::strategy`'s test
/// `a_big_integer_of_minus_two_to_the_63_is_written_as_its_digits`).
#[test]
fn the_strategy_db_length_is_d_of_any_sign() {
    let volume = unit(190, (190, 1), 1, 191);
    for length in [0, -1, 5000, -9_223_372_036_854_775_807, i64::MIN] {
        assert_eq!(
            units(vec![volume; 3], Some((length, 3)), 0).strategy_db_length(),
            Some(length)
        );
    }
    assert_eq!(units(vec![volume], None, 0).strategy_db_length(), None);
}

/// The setup's division by 0 without the Gumbel block, and the concatenated lengths.
#[test]
fn counts_of_zero_and_concatenated_lengths_are_refused() {
    let zero = units(vec![unit(560, (560, 0), 0, 560)], None, 0);
    assert!(zero.first().check_setup_count(true, "RPSBLAST").is_ok());
    assert!(zero.first().check_setup_count(false, "RPSBLAST").is_err());
    // The average of `NSEQ 1000` over 560 letters is 0 without the Gumbel block: NCBI's "too
    // short" (review A of D1a-4, MF5); the shortest length 10 with it.
    let short = units(vec![unit(560, (560, 1000), 1000, 1560)], None, 0);
    assert!(short.subject_length_too_short(false));
    assert!(!short.subject_length_too_short(true));
    assert_eq!(short.first().setup_subject_length(true), Some(10));
    assert!(!units(vec![unit(560, (560, 3), 3, 563)], None, 0).subject_length_too_short(false));
    let check = |engine_size: i128| {
        unit(560, (560, 3), 3, engine_size)
            .check_engine_length(40, "RPSBLAST")
            .map_err(|error| error.to_string())
    };
    assert!(check(563).is_ok());
    assert!(check(1003).is_ok());
    // LENGTH 0 and NSEQ -300 over the 563 rows of three profiles: searched (an `Lc` below the
    // rows, `RpsRun::concat_length`), also a negative one (LENGTH 0, NSEQ -5).
    assert!(check(3).is_ok());
    assert!(check(260).is_ok());
    assert!(check(-2).is_ok());
    // LENGTH 2^63 + 1 and NSEQ 600: `(Int4)` keeps 601.
    assert!(check(-9_223_372_036_854_775_807 + 600).is_ok());
    // The diagonal table: NCBI's loop ends at 2^30 (`Lc` plus the window 40 of exactly 2^30,
    // also with `(Int4)` of a negative sum) and never past it (LENGTH 2^30).
    assert!(check((1 << 30) - 40).is_ok());
    assert!(check((1 << 30) - 40 - (1i128 << 32)).is_ok());
    assert!(check((1 << 30) - 39)
        .unwrap_err()
        .contains("never stops growing"));
    assert!(check(1_073_741_827)
        .unwrap_err()
        .contains("of 1073741827 concatenated positions"));
    // LENGTH 2^31 - 4 and NSEQ 3: `Lc` is INT4_MAX, whose sum with the window overflows NCBI's
    // `Int4`; the widened sum is above 2^30.
    assert!(check(2_147_483_647)
        .unwrap_err()
        .contains("never stops growing"));
    // LENGTH 2^31: above INT4_MAX (the double free comes first in NCBI's engine).
    assert!(check(2_147_483_651)
        .unwrap_err()
        .contains("above 2147483647"));
    // The window of the one-hit finder (`-window_size 0`): the sum is `Lc` alone.
    let one_hit = |engine_size: i128| {
        unit(560, (560, 3), 3, engine_size)
            .check_engine_length(0, "RPSBLAST")
            .is_ok()
    };
    assert!(one_hit(1 << 30));
    assert!(!one_hit((1 << 30) + 1));
}
