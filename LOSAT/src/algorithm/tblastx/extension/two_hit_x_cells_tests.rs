use super::*;

// No NCBI counterpart (a statistics helper). The counting copy returns what the extension
// returns, on random NCBIstdaa sequences, offsets and X-drops, with the scores and the word
// length of the search: BLOSUM62's instance (`extend_hit_two_hit_fitting`), and every matrix at
// words of 2, 3 and 4 letters (`extend_hit_two_hit_reference`).
#[test]
fn the_counting_copy_returns_what_the_extension_returns() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move |n: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % n
    };
    let tables: Vec<crate::algorithm::tblastx::scoring::MatrixScore> = [2, 3, 4]
        .into_iter()
        .flat_map(|word| {
            crate::config::ScoringMatrix::ALL
                .into_iter()
                .map(move |matrix| {
                    crate::algorithm::tblastx::scoring::MatrixScore::with_word(matrix, word)
                })
        })
        .collect();
    for round in 0..40_000 {
        let q_len = 4 + next(120) as usize;
        let s_len = 4 + next(120) as usize;
        // Mostly a few residues, so that runs of positive scores are common.
        let alphabet = 1 + next(24);
        let q: Vec<u8> = (0..q_len).map(|_| 1 + next(alphabet) as u8).collect();
        let s: Vec<u8> = (0..s_len).map(|_| 1 + next(alphabet) as u8).collect();
        let x_drop = 1 + next(40) as i32;
        // Even rounds: BLOSUM62's instance; odd rounds: a matrix and a word length.
        let sc = (round % 2 == 1).then(|| tables[next(tables.len() as u64) as usize]);
        let word = sc.map_or(3, |sc| sc.word());
        let q_right = next((q_len - word) as u64 + 1) as usize;
        let s_right = next((s_len - word) as u64 + 1) as usize;
        let s_left = s_right.saturating_sub(next(50) as usize);
        let (cells, result) = match sc {
            None => (
                x_two_hit_cells(Blosum62Score, &q, &s, s_left, s_right, q_right, x_drop),
                extend_hit_two_hit_fitting(
                    &q,
                    InBatch::ALONE,
                    &s,
                    s_left,
                    s_right,
                    q_right,
                    x_drop,
                    false,
                ),
            ),
            Some(sc) => (
                x_two_hit_cells(sc, &q, &s, s_left, s_right, q_right, x_drop),
                extend_hit_two_hit_reference(
                    sc,
                    &q,
                    InBatch::ALONE,
                    &s,
                    s_left,
                    s_right,
                    q_right,
                    x_drop,
                    false,
                ),
            ),
        };
        assert_eq!(cells.result, result, "round {round} word {word}");
        if !result.5 {
            assert_eq!(cells.right, 0);
        }
        assert!(cells.s_lo <= s_right && cells.s_hi >= s_right + word && cells.s_hi <= s_len);
    }
}
