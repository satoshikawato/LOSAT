use super::*;
use crate::algorithm::blastn::alignment::gapped::build_blastna_matrix;

// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:471-478
const NCBI_QUERY: [u8; 22] = [
    0, 1, 2, 3, 0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 0, 1, 2, 2, 1, 1, 0, 0,
];
const NCBI_SUBJECT: [u8; 22] = [
    0, 1, 2, 2, 14, 14, 1, 1, 2, 2, 3, 3, 3, 3, 0, 1, 8, 2, 1, 1, 14, 0,
];

/// `testReevaluateWithAmbiguitiesUngappedNucl` (blasthits_unit_test.cpp:500-600): the HSP
/// 0..20 on 0..20 (score 20) of a 1/-2 search is cut to 6..20 on 6..20 with score 11, and not
/// deleted. NCBI's test sets `word_params->cutoff_score_min` to 3; the re-evaluation reads
/// `word_params->cutoffs[context].cutoff_score` (blast_hits.c:687), not that one, so the 3
/// below is only the test's argument: the expected values are the same for a cutoff of 4 to 11.
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:500-585 testReevaluateWithAmbiguitiesUngappedNucl
#[test]
fn ncbi_test_reevaluate_with_ambiguities_ungapped_nucl() {
    let matrix = build_blastna_matrix(1, -2);
    let got =
        reevaluate_with_ambiguities_ungapped(&NCBI_QUERY, &NCBI_SUBJECT, 0, 0, 20, 3, &matrix);
    assert_eq!(
        got,
        Some(ReevaluatedUngapped {
            score: 11,
            q_offset: 6,
            q_end: 20,
            s_offset: 6,
            s_end: 20
        })
    );
}

/// The same HSP with the cutoff of the other NCBI constant (`kHitCutoff` 11) is kept, with 12
/// it is deleted (the best stretch scores 11).
#[test]
fn reevaluate_ungapped_cutoff_boundary() {
    let matrix = build_blastna_matrix(1, -2);
    let kept =
        reevaluate_with_ambiguities_ungapped(&NCBI_QUERY, &NCBI_SUBJECT, 0, 0, 20, 11, &matrix);
    assert_eq!(
        kept.map(|r| (r.score, r.q_offset, r.q_end)),
        Some((11, 6, 20))
    );
    assert_eq!(
        reevaluate_with_ambiguities_ungapped(&NCBI_QUERY, &NCBI_SUBJECT, 0, 0, 20, 12, &matrix),
        None
    );
}

/// A front part that reached the cutoff before the sum went negative stays the best part
/// unless the later stretch scores more; below the cutoff it is discarded at the break.
/// (Expected values worked out outside with the same recurrence; an `N` scores -1 at 1/-2.)
#[test]
fn reevaluate_ungapped_front_part_kept_only_above_cutoff() {
    let matrix = build_blastna_matrix(1, -2);
    // Two Ns inside: the sum stays positive, the whole stretch is one HSP of score 6.
    let q = [0u8, 1, 2, 3, 0, 0, 0, 1, 2, 3];
    let s = [0u8, 1, 2, 3, 0, 14, 14, 1, 2, 3];
    let got = reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 10, 4, &matrix).unwrap();
    assert_eq!(
        (got.score, got.q_offset, got.q_end, got.s_offset, got.s_end),
        (6, 0, 10, 0, 10)
    );
    assert_eq!(
        reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 10, 7, &matrix),
        None
    );
    // Three mismatches take the sum below zero (5 - 6): front 5, back 7.
    let q = [0u8, 1, 2, 3, 0, 0, 0, 0, 1, 2, 3, 0, 1, 2, 3];
    let s = [0u8, 1, 2, 3, 0, 2, 2, 2, 1, 2, 3, 0, 1, 2, 3];
    let got = reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 15, 4, &matrix).unwrap();
    assert_eq!(
        (got.score, got.q_offset, got.q_end, got.s_offset, got.s_end),
        (7, 8, 15, 8, 15)
    );
    // Cutoff 6: the front (5) is discarded at the break, the back (7) is kept; cutoff 8: both fail.
    let got = reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 15, 6, &matrix).unwrap();
    assert_eq!((got.score, got.q_offset, got.q_end), (7, 8, 15));
    assert_eq!(
        reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 15, 8, &matrix),
        None
    );
    // Front 5 kept at cutoff 4 when the back stretch (3) scores less; at cutoff 6 nothing is left.
    let q = [0u8, 1, 2, 3, 0, 0, 0, 0, 1, 2, 3];
    let s = [0u8, 1, 2, 3, 0, 2, 2, 2, 1, 2, 3];
    let got = reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 11, 4, &matrix).unwrap();
    assert_eq!(
        (got.score, got.q_offset, got.q_end, got.s_offset, got.s_end),
        (5, 0, 5, 0, 5)
    );
    assert_eq!(
        reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 11, 6, &matrix),
        None
    );
}

/// Offsets are those of the passed sequences (the HSP need not start at 0 on either).
#[test]
fn reevaluate_ungapped_with_offsets_and_sentinel_mask() {
    let matrix = build_blastna_matrix(2, -3);
    let q = [3u8, 3, 0, 1, 2, 3, 0, 1];
    let s = [9u8, 0, 1, 2, 3, 0, 1, 1, 1];
    // HSP: q 2..8 against s 1..7 (6 matches except the last: q1 vs s1 equals, 6 letters).
    let got = reevaluate_with_ambiguities_ungapped(&q, &s, 2, 1, 6, 5, &matrix).unwrap();
    assert_eq!(
        (got.score, got.q_offset, got.q_end, got.s_offset, got.s_end),
        (12, 2, 8, 1, 7)
    );
    // A query letter of 0x10 or more is masked with 0x0f (the translated mask is 0xff).
    let qm = [3u8, 3, 0x10, 1, 2, 3, 0, 1];
    let got2 = reevaluate_with_ambiguities_ungapped(&qm, &s, 2, 1, 6, 5, &matrix).unwrap();
    assert_eq!(got, got2);
}

#[test]
fn identities_count_blastna_codes_equal() {
    // N (14) against N (14) is an identity; N against A is not.
    let q = [0u8, 14, 14, 3, 2];
    let s = [0u8, 14, 0, 3, 1];
    assert_eq!(count_identities_ungapped(&q, &s, 0, 0, 5), 3);
    assert_eq!(count_identities_ungapped(&q, &s, 3, 3, 2), 1);
    assert_eq!(count_identities_ungapped(&q, &s, 2, 2, 0), 0);
}

/// The word cutoff of NCBI's `BLAST_Cutoffs` for a search space of `avg` times the smaller of
/// `avg` and twice `q`, with `TRUE` and 0 as its last two arguments: worked by hand
/// as `ceil(ln(K * searchsp / 0.05) / Lambda)` (values from the formula, computed outside).
#[test]
fn word_cutoff_is_ceil_of_etos_with_min_of_avg_and_two_q() {
    let kbp = KarlinParams {
        lambda: 1.33,
        k: 0.621,
        h: 1.12,
        alpha: 0.0,
        beta: 0.0,
    };
    // q 377 (2q 754), avg 373: searchsp 373 * 373 -> ceil(ln(0.621*139129/0.05)/1.33) = 11.
    assert_eq!(word_cutoff_score(377, 373, &kbp, i32::MAX), Ok(11));
    // A short query: 2q 600 < avg 20000 -> searchsp 600 * 20000 (not 300 * 20000).
    assert_eq!(word_cutoff_score(300, 20000, &kbp, i32::MAX), Ok(15));
    // The cutoff is capped by cutoff_score_max.
    assert_eq!(word_cutoff_score(377, 373, &kbp, 7), Ok(7));
    // A very large subject average: 100000 * 100000 -> 20.
    assert_eq!(word_cutoff_score(100_000, 100_000, &kbp, i32::MAX), Ok(20));
    let kbp2 = KarlinParams {
        lambda: 0.634,
        k: 0.408,
        h: 0.912,
        alpha: 0.0,
        beta: -2.0,
    };
    assert_eq!(word_cutoff_score(377, 373, &kbp2, i32::MAX), Ok(22));
}

#[test]
fn word_cutoff_refuses_overflowing_arithmetic() {
    let kbp = KarlinParams {
        lambda: 1.33,
        k: 0.621,
        h: 1.12,
        alpha: 0.0,
        beta: 0.0,
    };
    assert_eq!(
        word_cutoff_score(1 << 30, 1000, &kbp, i32::MAX),
        Err(CutoffLimit::QueryLength)
    );
    assert!(word_cutoff_score((1 << 30) - 1, 1000, &kbp, i32::MAX).is_ok());
    // The check the streaming of reports asks of every query before the first output.
    assert!(cutoff_context_length_supported((1 << 30) - 1));
    assert!(!cutoff_context_length_supported(1 << 30));
    assert!(!cutoff_context_length_supported(usize::MAX));
    // min(avg, 2q) * avg > i64::MAX (not reachable with NCBI's 32-bit averages; the function
    // takes the average as a u64).
    assert_eq!(
        word_cutoff_score((1 << 30) - 1, 5_000_000_000, &kbp, i32::MAX),
        Err(CutoffLimit::SearchSpace)
    );
}

#[test]
fn average_subject_lengths_follow_each_source() {
    // A database's average is `SeqDb::average_length` (seqsrc_seqdb.cpp:203-209).
    assert_eq!(average_subject_length_subject_set(1000, 3), 333);
    // The Uint4 cast of a huge average.
    assert_eq!(
        average_subject_length_subject_set(5_000_000_000, 1),
        705_032_704
    );
}
