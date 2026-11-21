use super::*;

fn extent_of(&(context, frame, offset, end): &(i64, i32, i64, i64)) -> BestHitExtent {
    BestHitExtent {
        context,
        query_frame: frame,
        offset,
        end,
        query_length: 100,
    }
}

#[test]
fn subject_besthit_removes_later_hsps_within_three_letters() {
    // The second lies within 10-3..50+3 of the first; the third starts 4 before it; the
    // fourth is of another context; the fifth lies within the first's range too.
    let mut hsps = vec![
        (0, 1, 10, 50),
        (0, 1, 7, 53),
        (0, 1, 6, 50),
        (1, -1, 10, 50),
        (0, 1, 8, 49),
    ];
    subject_besthit(&mut hsps, false, extent_of);
    assert_eq!(hsps, vec![(0, 1, 10, 50), (0, 1, 6, 50), (1, -1, 10, 50)]);
}

#[test]
fn subject_besthit_bounds_the_start_at_zero_and_a_negative_end() {
    // Start 1 - 3 is kept at 0: an HSP at -1 (outside its context) stays.
    let mut hsps = vec![(0, 1, 1, 20), (0, 1, -1, 20), (0, 1, 0, 22)];
    subject_besthit(&mut hsps, false, extent_of);
    assert_eq!(hsps, vec![(0, 1, 1, 20), (0, 1, -1, 20)]);
    // An end of -4 (below -3) keeps -4 as the bound: the later HSP ending at -3 stays.
    let mut hsps = vec![(0, 1, 0, -4), (0, 1, 0, -3), (0, 1, 0, -5)];
    subject_besthit(&mut hsps, false, extent_of);
    assert_eq!(hsps, vec![(0, 1, 0, -4), (0, 1, 0, -3)]);
}

#[test]
fn subject_besthit_reflects_onto_the_other_strand_for_blastn_only() {
    // On a context of 100 letters, 10..50 of the plus strand is 50..90 of the minus strand,
    // widened to 47..93: the minus strand's 48..92 goes, its 30..92 stays.
    let hsps = vec![(0, 1, 10, 50), (1, -1, 48, 92), (1, -1, 30, 92)];
    let mut blastn = hsps.clone();
    subject_besthit(&mut blastn, true, extent_of);
    assert_eq!(blastn, vec![(0, 1, 10, 50), (1, -1, 30, 92)]);
    let mut protein = hsps.clone();
    subject_besthit(&mut protein, false, extent_of);
    assert_eq!(protein, hsps);
}

/// NCBI's `testSubjectBestHit`: four HSPs of context 0 at query offset 0, ending at 12, 6, 4
/// and 3 (12 / (i + 1)) in score order, with the query frame and length NCBI's test gives
/// them (frame 0 of a new HSP, the 12-letter query): one HSP is left for tblastn and for
/// blastn. The HSP lists are NCBI's two fixtures:
/// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:1386-1402
/// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:1495-1512
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:1514-1542 testSubjectBestHit
#[test]
fn subject_besthit_keeps_one_of_ncbi_four_nested_hsps() {
    for blastn in [false, true] {
        let mut ends: Vec<i64> = (1..=4).map(|factor| 12 / factor).collect();
        subject_besthit(&mut ends, blastn, |&end| BestHitExtent {
            context: 0,
            query_frame: 0,
            offset: 0,
            end,
            query_length: 12,
        });
        assert_eq!(ends, vec![12], "blastn {blastn}");
    }
}

#[test]
fn query_coverage_adds_a_half_below_99() {
    // 95 of 1000 letters: 9.5 % + 0.5.
    assert_eq!(hsp_query_coverage(5, 100, 1000), 10.0);
    // 988 of 1000: 98.8 % + 0.5 = 99.3; 992 of 1000 stays 99.2.
    assert_eq!(hsp_query_coverage(0, 988, 1000), 98.8 + 0.5);
    assert_eq!(hsp_query_coverage(0, 992, 1000), 99.2);
    assert_eq!(hsp_query_coverage(0, 10, 0), 0.0);
}

#[test]
fn reap_keeps_the_order_and_reports_removals() {
    let mut hsps = vec![(0, 50), (0, 95), (10, 100), (0, 9)];
    let extent = |&(offset, end): &(i64, i64)| (offset, end, 1000);
    assert!(!reap_by_query_coverage(&mut hsps, 0.0, extent));
    assert_eq!(hsps.len(), 4);
    // Coverage 5.5, 10.0, 9.5, 1.4: a limit of 9.5 keeps the second and the third.
    assert!(reap_by_query_coverage(&mut hsps, 9.5, extent));
    assert_eq!(hsps, vec![(0, 95), (10, 100)]);
    assert!(!reap_by_query_coverage(&mut hsps, 9.5, extent));
}

#[test]
fn hsp_test_deletes_a_negative_length_at_the_default_options() {
    // BLASTX's options: percent_identity 0 and min_hit_length 0. A negative length (an
    // ungapped HSP whose end lies before its start) is deleted by the length test (its
    // identity test is `0.0 < -0.0`, false); a length of 0 is kept: both tests are false
    // (`0.0 < 0.0`, `0 < 0`).
    for length in [-4, -1] {
        assert!(hsp_test(0, length, 0.0, 0));
    }
    assert!(!hsp_test(0, 0, 0.0, 0));
    assert!(!hsp_test(3, 7, 0.0, 0));
    // `-perc_identity -0` (in NCBI's range): `0.0 < -0.0` and `500.0 < -0.0` are false.
    assert!(!hsp_test(0, 5, -0.0, 0));
    assert!(!hsp_test(5, 5, -0.0, 0));
    // The identity test: 3 of 7 is 42.9 %, below 50 and not below 40.
    assert!(hsp_test(3, 7, 50.0, 0));
    assert!(!hsp_test(3, 7, 40.0, 0));
    // The minimal hit length.
    assert!(hsp_test(7, 7, 0.0, 8));
    assert!(!hsp_test(7, 7, 0.0, 7));
}

/// NCBI's C API (`Blast_HSPTest` called directly on the pinned NCBI build) on the
/// boundaries: identity equal to the limit kept and just below deleted, the length test,
/// and the empty HSP at the default options.
#[test]
fn hsp_test_matches_the_ncbi_api_oracle() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/hsp_test_20260924/oracle.tsv"
    );
    let oracle = std::fs::read_to_string(path).unwrap();
    for line in oracle.lines() {
        let f: Vec<_> = line.split('\t').collect();
        let observed = hsp_test(
            f[1].parse().unwrap(),
            f[2].parse().unwrap(),
            f[3].parse().unwrap(),
            f[4].parse().unwrap(),
        );
        assert_eq!(observed, f[5] == "1", "{line}");
    }
    assert_eq!(oracle.lines().count(), 6);
}

/// No NCBI counterpart: NCBI's `Int4` cannot hold these counts. The `usize` counts of
/// BLASTN, TBLASTN and TBLASTX are NCBI's `Int4` widened to `i64`: a length past `Int4` is
/// tested as itself, not capped.
#[test]
fn hsp_test_takes_lengths_past_int4() {
    let long = i64::from(i32::MAX) + 10;
    assert!(!hsp_test(long, long, 100.0, 0));
    assert!(hsp_test(long - 1, long, 100.0, 0));
    assert!(!hsp_test(0, long, 0.0, long));
    assert!(hsp_test(0, long - 1, 0.0, long));
}

/// NCBI's test `stop < length/3 - 3` with `stop` 0 (an `Int4` division, negative below 9
/// bases) holds from 12 bases on.
#[test]
fn a_target_translation_is_never_made_below_12_bases() {
    for length in 0..60usize {
        let ncbi_translates = 0 < length as i32 / 3 - 3;
        assert_eq!(
            target_translation_never_made(length),
            !ncbi_translates,
            "{length}"
        );
    }
    assert!(target_translation_never_made(11));
    assert!(!target_translation_never_made(12));
    assert_eq!(
        x_untranslated_frame_refusal("TBLASTX", "the ungapped re-evaluation").to_string(),
        "TBLASTX: the ungapped re-evaluation of an alignment on a subject (or its -subject_loc range) of fewer than 12 bases, where NCBI BLAST+ reads the translation of the subject's frame, which it never makes for a subject this short, and crashes; LOSAT does not support it"
    );
}

#[test]
fn trim_keeps_the_first_hsps() {
    let mut hsps = vec![3, 1, 2];
    trim_by_max_hsps(&mut hsps, 0);
    assert_eq!(hsps, vec![3, 1, 2]);
    trim_by_max_hsps(&mut hsps, 3);
    assert_eq!(hsps, vec![3, 1, 2]);
    trim_by_max_hsps(&mut hsps, 2);
    assert_eq!(hsps, vec![3, 1]);
}
