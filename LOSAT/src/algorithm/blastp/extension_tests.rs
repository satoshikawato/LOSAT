use super::*;
use crate::utils::matrix::ncbistdaa;

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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 and :1020-1083
// (cited at `extend_two_hit_outside`)
// On a hit inside its chunk the extensions on the letters around the chunk are the chunk's own:
// every two-hit and one-hit extension of a chunk of 40 letters between 7 before it and 5 after it
// gives what `extend_two_hit` and `extend_one_hit` give on the chunk's letters.
#[test]
fn extensions_around_a_chunk_match_the_chunk_s_inside_it() {
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut letters = || -> Vec<u8> {
        (0..52)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                1 + (state % 22) as u8
            })
            .collect()
    };
    let query = letters();
    let around = letters();
    let chunk = ChunkLetters {
        letters: &around,
        lead: 7,
        length: 40,
    };
    let inside = &around[7..47];
    for matrix in [ScoringMatrix::Blosum62, ScoringMatrix::Blosum45] {
        for q in 0..query.len() - 3 {
            for s in 0..inside.len() - 3 {
                for gap in 3..8 {
                    let Some(s_left) = s.checked_sub(gap) else {
                        continue;
                    };
                    let want = extend_two_hit(matrix, &query, inside, s_left + 3, s, q, 12, 3);
                    let got = extend_two_hit_outside(
                        matrix,
                        &query,
                        chunk,
                        (s_left + 3) as i32,
                        s as i32,
                        q as i32,
                        12,
                        3,
                    );
                    assert_eq!(want, Some(got), "{matrix:?} two-hit q {q} s {s} gap {gap}");
                }
                let want = extend_one_hit(matrix, &query, inside, q, s, 12, 3);
                let got = extend_one_hit_outside(matrix, &query, chunk, q as i32, s as i32, 12, 3);
                assert_eq!(want, Some(got), "{matrix:?} one-hit q {q} s {s}");
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1104-1157, :995 and :948
// A two-hit extension whose first hit lies before the chunk (masks stored out of order): the left
// walk stops at the chunk's first letter, and the right walk runs only where the left one reaches
// the end of the first hit, 2 letters before the chunk, which it does not.
#[test]
fn a_two_hit_extension_from_before_the_chunk_stops_its_left_walk_at_the_chunk() {
    let w = ncbistdaa::W;
    let a = ncbistdaa::A;
    let query = vec![w; 6];
    let around = vec![w, w, w, w, w, w, a, a, a];
    let chunk = ChunkLetters {
        letters: &around,
        lead: 3,
        length: 6,
    };
    let result = extend_two_hit_outside(ScoringMatrix::Blosum62, &query, chunk, -2, 0, 3, 7, 3);
    let ww = protein_score(ScoringMatrix::Blosum62, w, w);
    assert_eq!(
        result,
        BlastpTwoHitUngappedResult {
            ungapped_data: BlastpUngappedData {
                q_start: 3,
                s_start: 0,
                length: 3,
                score: 3 * ww,
            },
            right_extend: false,
            s_last_off: 3,
        }
    );
    // The first hit ending at the chunk's first letter: the left walk reaches it, and the right
    // walk runs (the query ends there).
    let result = extend_two_hit_outside(ScoringMatrix::Blosum62, &query, chunk, 0, 0, 3, 7, 3);
    assert!(result.right_extend);
    assert_eq!(result.ungapped_data.s_start, 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1036-1085, :995 and :948
// A one-hit extension of a word 2 letters before the chunk: no left walk (it would start before
// the chunk), and the HSP starts before the chunk, at the word.
#[test]
fn a_one_hit_extension_of_a_word_before_the_chunk_starts_there() {
    let w = ncbistdaa::W;
    let a = ncbistdaa::A;
    let query = vec![w; 3];
    let around = vec![a, w, w, w, a, a];
    let chunk = ChunkLetters {
        letters: &around,
        lead: 3,
        length: 3,
    };
    let result = extend_one_hit_outside(ScoringMatrix::Blosum62, &query, chunk, 0, -2, 7, 3);
    let ww = protein_score(ScoringMatrix::Blosum62, w, w);
    assert_eq!(
        result,
        BlastpOneHitUngappedResult {
            ungapped_data: BlastpUngappedData {
                q_start: 0,
                s_start: -2,
                length: 3,
                score: 3 * ww,
            },
            s_last_off: 1,
        }
    );
}
