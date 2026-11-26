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

#[derive(Clone, Debug, PartialEq, Eq)]
struct TestChunkHsp {
    context: i32,
    subject_frame: i32,
    q_start: i32,
    q_end: i32,
    s_start: i32,
    s_end: i32,
    q_gapped_start: i32,
    s_gapped_start: i32,
    score: i32,
    tag: i32,
}

impl ChunkMergeHsp for TestChunkHsp {
    fn merge_context(&self) -> i32 {
        self.context
    }
    fn merge_subject_frame(&self) -> i32 {
        self.subject_frame
    }
    fn merge_q_start(&self) -> i32 {
        self.q_start
    }
    fn merge_q_end(&self) -> i32 {
        self.q_end
    }
    fn merge_s_start(&self) -> i32 {
        self.s_start
    }
    fn merge_s_end(&self) -> i32 {
        self.s_end
    }
    fn merge_score(&self) -> i32 {
        self.score
    }
    fn merge_q_gapped_start(&self) -> i32 {
        self.q_gapped_start
    }
    fn merge_s_gapped_start(&self) -> i32 {
        self.s_gapped_start
    }
    fn set_merge_bounds(&mut self, q_start: i32, q_end: i32, s_start: i32, s_end: i32) {
        self.q_start = q_start;
        self.q_end = q_end;
        self.s_start = s_start;
        self.s_end = s_end;
    }
    fn set_merge_gapped_start(&mut self, q_gapped_start: i32, s_gapped_start: i32) {
        self.q_gapped_start = q_gapped_start;
        self.s_gapped_start = s_gapped_start;
    }
    fn set_merge_score(&mut self, score: i32) {
        self.score = score;
    }
}

fn test_score_cmp(a: &TestChunkHsp, b: &TestChunkHsp) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| a.s_start.cmp(&b.s_start))
        .then_with(|| b.s_end.cmp(&a.s_end))
        .then_with(|| a.q_start.cmp(&b.q_start))
        .then_with(|| b.q_end.cmp(&a.q_end))
}

fn hsp(context: i32, q: (i32, i32), s: (i32, i32), gapped: (i32, i32), score: i32) -> TestChunkHsp {
    TestChunkHsp {
        context,
        subject_frame: 0,
        q_start: q.0,
        q_end: q.1,
        s_start: s.0,
        s_end: s.1,
        q_gapped_start: gapped.0,
        s_gapped_start: gapped.1,
        score,
        tag: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1464-1537,2857-3035
#[test]
fn merge_subject_chunk_hsps_joins_gapped_and_ungapped_across_overlap() {
    // Gapped merge with diagonal shift of 2 (< 10) across a 5,000,000-residue boundary:
    // higher-scoring incoming HSP transfers its gapped_start, and score uses density scaling.
    let mut combined = vec![hsp(
        0,
        (0, 180),
        (4_999_840, 5_000_020),
        (40, 4_999_880),
        400,
    )];
    let mut incoming = vec![hsp(
        0,
        (78, 218),
        (4_999_920, 5_000_060),
        (120, 4_999_962),
        500,
    )];
    merge_subject_chunk_hsps(
        &mut combined,
        &mut incoming,
        4_999_920,
        100,
        i32::MAX as usize,
        true,
        test_score_cmp,
    )
    .unwrap();
    assert!(incoming.is_empty());
    assert_eq!(combined.len(), 1);
    assert_eq!(
        (
            combined[0].q_start,
            combined[0].q_end,
            combined[0].s_start,
            combined[0].s_end
        ),
        (0, 218, 4_999_840, 5_000_060)
    );
    assert_eq!(
        (combined[0].q_gapped_start, combined[0].s_gapped_start),
        (120, 4_999_962)
    );
    // density = (400 + 500) / (180 + 140) = 900 / 320 = 2.8125; 2.8125 * 218 = 613.125 -> 613
    assert_eq!(combined[0].score, 613);

    // Same off-diagonal pair in ungapped mode (allow_gap = false): does NOT merge.
    let mut ungapped = vec![hsp(
        0,
        (0, 180),
        (4_999_840, 5_000_020),
        (40, 4_999_880),
        400,
    )];
    let mut incoming = vec![hsp(
        0,
        (78, 218),
        (4_999_920, 5_000_060),
        (120, 4_999_962),
        500,
    )];
    merge_subject_chunk_hsps(
        &mut ungapped,
        &mut incoming,
        4_999_920,
        100,
        i32::MAX as usize,
        false,
        test_score_cmp,
    )
    .unwrap();
    assert!(incoming.is_empty());
    assert_eq!(ungapped.len(), 2);
    assert_eq!(ungapped[0].score, 500);
    assert_eq!(ungapped[1].score, 400);

    // Diagonal difference of 10 (not < OVERLAP_DIAG_CLOSE): does NOT merge even in gapped mode.
    let mut far_diag = vec![hsp(0, (0, 180), (100, 280), (10, 110), 400)];
    let mut incoming = vec![hsp(0, (40, 190), (130, 280), (50, 140), 450)];
    merge_subject_chunk_hsps(
        &mut far_diag,
        &mut incoming,
        130,
        100,
        i32::MAX as usize,
        true,
        test_score_cmp,
    )
    .unwrap();
    assert_eq!(far_diag.len(), 2);
}

#[test]
fn merge_subject_chunk_hsps_respects_zero_overlap_context_frame_and_undefined_zero_sum() {
    // Zero overlap: incoming starting at split_offset (130) is NOT < 130 + 0, so no merge.
    let mut combined = vec![hsp(0, (0, 140), (110, 250), (10, 120), 600)];
    let mut incoming = vec![hsp(0, (20, 150), (130, 260), (30, 140), 550)];
    merge_subject_chunk_hsps(
        &mut combined,
        &mut incoming,
        130,
        0,
        i32::MAX as usize,
        false,
        test_score_cmp,
    )
    .unwrap();
    assert_eq!(combined.len(), 2);

    // Different context or different subject frame: no merge.
    let mut diff_ctx = vec![hsp(0, (0, 140), (110, 250), (10, 120), 600)];
    let mut incoming = vec![hsp(1, (20, 150), (130, 260), (30, 140), 550)];
    merge_subject_chunk_hsps(
        &mut diff_ctx,
        &mut incoming,
        130,
        100,
        i32::MAX as usize,
        true,
        test_score_cmp,
    )
    .unwrap();
    assert_eq!(diff_ctx.len(), 2);

    let mut diff_frame = vec![hsp(0, (0, 140), (110, 250), (10, 120), 600)];
    let mut neg_frame = hsp(0, (20, 150), (130, 260), (30, 140), 550);
    neg_frame.subject_frame = -1;
    let mut incoming = vec![neg_frame];
    merge_subject_chunk_hsps(
        &mut diff_frame,
        &mut incoming,
        130,
        100,
        i32::MAX as usize,
        true,
        test_score_cmp,
    )
    .unwrap();
    assert_eq!(diff_frame.len(), 2);

    // Zero query length sum in ungapped one-hit merge: returns ChunkMergeUndefined.
    let mut zero_len = vec![hsp(0, (10, 10), (150, 150), (10, 150), 15)];
    let mut incoming = vec![hsp(0, (10, 10), (150, 150), (10, 150), 15)];
    assert_eq!(
        merge_subject_chunk_hsps(
            &mut zero_len,
            &mut incoming,
            130,
            100,
            i32::MAX as usize,
            false,
            test_score_cmp,
        ),
        Err(ChunkMergeUndefined)
    );
}

#[test]
fn merge_subject_chunk_hsps_preserves_swap_displacement_and_hsp_num_max_tie_order() {
    // Three equal-score HSPs in `combined`: non-strip N1 (tag 1), N2 (tag 2), strip S (tag 3).
    // Swap moves S to index 0, displacing N1 to index 2 -> [S(3), N2(2), N1(1)].
    // Adding an equal-score non-merging incoming I (tag 4) with identical coordinates keeps
    // NCBI's swap order under stable sort.
    let make = |s_end, tag| TestChunkHsp {
        context: 0,
        subject_frame: 0,
        q_start: 0,
        q_end: 50,
        s_start: 10,
        s_end,
        q_gapped_start: 10,
        s_gapped_start: 20,
        score: 100,
        tag,
    };
    let mut combined = vec![make(80, 1), make(80, 2), make(150, 3)];
    let mut incoming = vec![make(80, 4)];
    let incoming_cap = incoming.capacity();
    merge_subject_chunk_hsps(
        &mut combined,
        &mut incoming,
        100,
        0,
        i32::MAX as usize,
        true,
        test_score_cmp,
    )
    .unwrap();
    assert!(incoming.is_empty());
    assert!(incoming.capacity() >= incoming_cap);
    // S (s_end=150) sorts before s_end=80; among the three s_end=80 ties, N2(2) precedes N1(1)
    // because S swapped with N1 at index 0, and incoming I(4) follows both.
    let tags: Vec<i32> = combined.iter().map(|h| h.tag).collect();
    assert_eq!(tags, vec![3, 2, 1, 4]);

    // With hsp_num_max = 2, ties prefer `combined` over `incoming`.
    let mut capped = vec![make(80, 10), make(80, 11)];
    let mut incoming = vec![make(80, 20), make(80, 21)];
    merge_subject_chunk_hsps(&mut capped, &mut incoming, 100, 0, 2, true, test_score_cmp).unwrap();
    assert!(incoming.is_empty());
    let tags: Vec<i32> = capped.iter().map(|h| h.tag).collect();
    assert_eq!(tags, vec![10, 11]);
}
