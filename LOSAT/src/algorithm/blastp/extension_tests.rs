use super::*;
use crate::algorithm::common::hard_chunk_scan::OneWord;
use crate::utils::matrix::{ncbistdaa, protein_score_rows};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1078-1127
#[test]
fn test_extend_one_hit_returns_ncbi_blast_ungapped_data_payload() {
    let query = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
    ];
    let subject = query.clone();

    let result = extend_one_hit(ScoringMatrix::Blosum62, &query, &subject, 1, 1, 7, 3)
        .expect("one-hit ungapped extension");

    assert_eq!(result.ungapped_data.q_start, 0);
    assert_eq!(result.ungapped_data.s_start, 0);
    assert_eq!(result.ungapped_data.length, 5);
    assert_eq!(
        result.ungapped_data.score,
        5 * protein_score(ScoringMatrix::Blosum62, ncbistdaa::A, ncbistdaa::A)
    );
    assert_eq!(result.s_last_off, 5);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1124-1158
#[test]
fn test_extend_two_hit_returns_ncbi_payload_when_left_reaches_first_hit() {
    let query = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
    ];
    let subject = query.clone();

    let result = extend_two_hit(ScoringMatrix::Blosum62, &query, &subject, 0, 3, 3, 7, 3)
        .expect("two-hit ungapped extension");

    assert!(result.right_extend);
    assert_eq!(result.ungapped_data.q_start, 0);
    assert_eq!(result.ungapped_data.s_start, 0);
    assert_eq!(result.ungapped_data.length, 6);
    assert_eq!(
        result.ungapped_data.score,
        6 * protein_score(ScoringMatrix::Blosum62, ncbistdaa::A, ncbistdaa::A)
    );
    assert_eq!(result.s_last_off, 6);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1129-1158
#[test]
fn test_extend_two_hit_returns_ncbi_payload_when_left_stops_before_first_hit() {
    let query = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::D,
        ncbistdaa::D,
        ncbistdaa::D,
        ncbistdaa::C,
        ncbistdaa::C,
        ncbistdaa::C,
    ];
    let subject = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::W,
        ncbistdaa::W,
        ncbistdaa::W,
        ncbistdaa::C,
        ncbistdaa::C,
        ncbistdaa::C,
    ];

    let result = extend_two_hit(ScoringMatrix::Blosum62, &query, &subject, 0, 6, 6, 2, 3)
        .expect("two-hit ungapped extension");

    assert!(!result.right_extend);
    assert_eq!(result.ungapped_data.q_start, 6);
    assert_eq!(result.ungapped_data.s_start, 6);
    assert_eq!(result.ungapped_data.length, 3);
    assert_eq!(
        result.ungapped_data.score,
        3 * protein_score(ScoringMatrix::Blosum62, ncbistdaa::C, ncbistdaa::C)
    );
    assert_eq!(result.s_last_off, 9);
}

/// The one word of a hard chunk shorter than the word (`common::hard_chunk_scan::OneWord`): where
/// the chunk holds the word, its extension gives the score, the direction and the last
/// offset of LOSAT's reference two-hit extension of the chunk's letters
/// (`extend_two_hit`) at subject offset 0, for any first hit's end, negative ones
/// included.
#[test]
fn one_word_extension_is_the_reference_two_hit_extension() {
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = |below: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % below
    };
    let mut extended = 0;
    for _ in 0..20_000 {
        let word = 2 + next(4) as i32;
        let length = word as usize + next(8) as usize;
        let letters: Vec<u8> = (0..length).map(|_| 1 + next(20) as u8).collect();
        let query: Vec<u8> = (0..word as usize + next(12) as usize)
            .map(|_| 1 + next(20) as u8)
            .collect();
        let q = next((query.len() - word as usize + 1) as u64) as i32;
        let s_left_off = next(20) as i32 - 16;
        let x_drop = 1 + next(30) as i32;
        let test = OneWord {
            query: &query,
            letters: &letters,
            length: length as i32,
            word,
            window: 40,
            scores: protein_score_rows(ScoringMatrix::Blosum62),
        };
        let got = test.extend(s_left_off, q, x_drop);
        let want = extend_two_hit(
            ScoringMatrix::Blosum62,
            &query,
            &letters,
            s_left_off,
            0,
            q as usize,
            x_drop,
            word as usize,
        )
        .expect("the word lies within both sequences");
        let data = want.ungapped_data;
        assert_eq!(
            (got.score, got.right_extend, got.s_last_off),
            (data.score, want.right_extend, want.s_last_off),
            "{query:?} {letters:?} q={q} s_left_off={s_left_off} x_drop={x_drop}"
        );
        assert_eq!(
            (got.q_start, got.s_start, got.length),
            (data.q_start, data.s_start, data.length),
            "{query:?} {letters:?} q={q} s_left_off={s_left_off} x_drop={x_drop}"
        );
        extended += usize::from(got.right_extend);
    }
    assert!(extended > 1_000, "only {extended} right walks");
}
