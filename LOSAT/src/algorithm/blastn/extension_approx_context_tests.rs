use super::{build_nucl_score_table, extend_hit_ungapped_approx_ncbi};
use crate::algorithm::blastn::alignment::build_blastna_matrix;

fn blastna(text: &str) -> Vec<u8> {
    text.bytes()
        .map(|letter| match letter {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => unreachable!(),
        })
        .collect()
}

fn packed(text: &str) -> Vec<u8> {
    blastna(text)
        .chunks(4)
        .map(|chunk| {
            (0..4).fold(0u8, |byte, index| {
                byte << 2 | chunk.get(index).copied().unwrap_or(0)
            })
        })
        .collect()
}

// Two contexts and the sentinel between them (blast_encoding.c:121, kNuclSentinel = 0xF):
// context 0 at 0..16, the sentinel at 16, context 1 at 17..35. The seed is the first letter
// of context 1 and the subject offset is 2 mod 4, so the first left block holds the last
// letter of context 0, the sentinel and the first two letters of context 1.
#[test]
fn the_approximate_extension_crosses_the_sentinel_and_its_score_decides_the_rerun() {
    let mut query = blastna("ACGTACGTACGTACGT");
    query.push(15);
    query.extend(blastna("GACCGGTTAATTTTGGGG"));
    let context_start = 17;
    let subject = packed("GTACTTGACCGGTTAATTAA");
    let subject_len = 20;
    let score_table = build_nucl_score_table(2, -3);
    let matrix = build_blastna_matrix(2, -3);
    let (q_off, s_off, s_match_end, x_dropoff) = (context_start, 6, 14, 16);
    let extend = |reduced_cutoff: i32| {
        extend_hit_ungapped_approx_ncbi(
            &query,
            &subject,
            q_off,
            s_off,
            s_match_end,
            subject_len,
            x_dropoff,
            &score_table,
            reduced_cutoff,
            &matrix,
        )
    };

    // Approximate score 24 below the reduced cutoff 25 (cutoff 32): the approximate result
    // is kept and starts two letters left of context 1. Such a result is rejected, as its
    // score is below the cutoff.
    let kept = extend(25);
    assert!(kept.q_start < context_start);
    assert_eq!(
        (kept.q_start, kept.s_start, kept.length, kept.score),
        (15, 4, 12, 24)
    );

    // Reduced cutoff 19 (cutoff 24): the crossing block lifts the approximate score to 24,
    // so the exact extension reruns; it stops at the sentinel and scores 24, the cutoff.
    let rerun = extend(19);
    assert_eq!(
        (rerun.q_start, rerun.s_start, rerun.length, rerun.score),
        (17, 6, 12, 24)
    );

    // The same seed on context 1 alone (not NCBI's extension): the first two letters are
    // never scored, the approximate score is 16 < 19, no rerun, and 16 < 24 would reject
    // the hit that NCBI saves.
    let clamped = extend_hit_ungapped_approx_ncbi(
        &query[context_start..],
        &subject,
        0,
        s_off,
        s_match_end,
        subject_len,
        x_dropoff,
        &score_table,
        19,
        &matrix,
    );
    assert_eq!(
        (
            clamped.q_start,
            clamped.s_start,
            clamped.length,
            clamped.score
        ),
        (2, 8, 8, 16)
    );
}
