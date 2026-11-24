use super::super::extension::{
    build_nucl_score_table, extend_hit_ungapped_approx_ncbi, extend_hit_ungapped_exact_ported,
    query_four_base_byte,
};
use super::*;

fn next(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn key(u: UngappedData) -> (usize, usize, usize, i32) {
    (u.q_start, u.s_start, u.length, u.score)
}

/// Every word of four letters below 16 gives the byte of `query_four_base_byte`.
#[test]
fn the_word_gives_the_overlapping_byte_for_every_letter_quadruple() {
    for w in 0u32..1 << 16 {
        let q = [
            (w & 15) as u8,
            ((w >> 4) & 15) as u8,
            ((w >> 8) & 15) as u8,
            (w >> 12) as u8,
        ];
        assert_eq!(
            four_base_byte_of_word(u32::from_le_bytes(q)),
            query_four_base_byte(&q, 0),
            "{q:?}"
        );
    }
}

/// Random query blocks (contexts of 1 to 40 letters between sentinels, a few ambiguity codes)
/// against random and periodic subjects of 1 to 140 bases packed into exactly `len.div_ceil(4)`
/// bytes (the windows' short paths at both ends and their refills), every hit offset.
fn cases(rounds: usize, mut f: impl FnMut(&[u8], &[u8], usize, usize, usize, usize)) {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    for round in 0..rounds {
        let mut query = Vec::new();
        for _ in 0..1 + next(&mut state) % 4 {
            for _ in 0..1 + next(&mut state) % 40 {
                let r = next(&mut state);
                query.push(if r.is_multiple_of(50) {
                    (4 + r % 11) as u8
                } else {
                    (r % 4) as u8
                });
            }
            query.push(15);
        }
        query.pop();
        let subject_len = 1 + (next(&mut state) % 140) as usize;
        let period = 1 + (next(&mut state) % 6) as usize;
        let motif: Vec<u8> = (0..period).map(|_| (next(&mut state) % 4) as u8).collect();
        // Odd rounds copy the query's letters into the subject, so long walks happen.
        let bases: Vec<u8> = (0..subject_len)
            .map(|i| match round % 3 {
                0 => motif[i % period],
                1 => (next(&mut state) % 4) as u8,
                _ => query.get(i).map_or(0, |&l| l & 3),
            })
            .collect();
        let mut packed = vec![0u8; subject_len.div_ceil(4)];
        for (i, b) in bases.iter().enumerate() {
            packed[i / 4] |= b << (2 * (3 - i % 4));
        }
        for q_off in 0..query.len() {
            for s_off in 0..subject_len {
                f(&query, &packed, q_off, s_off, subject_len, round);
            }
        }
    }
}

#[test]
fn the_window_loops_give_the_ported_loops_results() {
    let mut calls = 0u64;
    for (reward, penalty) in [(1, -3), (2, -3), (1, -1)] {
        let matrix = crate::algorithm::blastn::alignment::build_blastna_matrix(reward, penalty);
        cases(120, |query, packed, q_off, s_off, subject_len, round| {
            for x_dropoff in [5, 16, 40, i32::MIN] {
                assert_eq!(
                    key(exact_window(
                        query,
                        packed,
                        q_off,
                        s_off,
                        subject_len,
                        x_dropoff,
                        &matrix
                    )),
                    key(extend_hit_ungapped_exact_ported(
                        query,
                        packed,
                        q_off,
                        s_off,
                        subject_len,
                        x_dropoff,
                        &matrix
                    )),
                    "round {round} {:?}",
                    (q_off, s_off, subject_len, x_dropoff)
                );
                calls += 1;
            }
        });
    }
    assert!(calls > 1_000_000, "{calls}");
}

/// The block loops give today's approximate loops' results (with no probe bit set, the tests'
/// `extend_hit_ungapped_approx_ncbi` runs today's loops), at reduced cutoffs that keep the
/// approximate result and that rerun the exact extension. The query offsets go up to 3 letters
/// past a hit's usual place, as the extension's own blocks do.
#[test]
fn the_block_loops_give_todays_approximate_results() {
    let mut calls = 0u64;
    for (reward, penalty) in [(1, -3), (2, -3), (1, -2)] {
        let matrix = crate::algorithm::blastn::alignment::build_blastna_matrix(reward, penalty);
        let table = build_nucl_score_table(reward, penalty);
        cases(100, |query, packed, q_off, s_off, subject_len, round| {
            // Today's loops read the query to `q_off + 3` before their first block: skip the
            // offsets where that is past the query (both panic there).
            let len = (4 - s_off % 4) % 4;
            if q_off + len > query.len() {
                return;
            }
            for x_dropoff in [5, 16, 40, i32::MIN] {
                for reduced_cutoff in [i32::MAX, 12, 0] {
                    let s_match_end = (s_off + 11).min(subject_len);
                    let run = |f: fn(
                        &[u8],
                        &[u8],
                        usize,
                        usize,
                        usize,
                        usize,
                        i32,
                        &[i32; 256],
                        i32,
                        &[i32; 256],
                    ) -> UngappedData| {
                        key(f(
                            query,
                            packed,
                            q_off,
                            s_off,
                            s_match_end,
                            subject_len,
                            x_dropoff,
                            &table,
                            reduced_cutoff,
                            &matrix,
                        ))
                    };
                    assert_eq!(
                        run(approx_blocks),
                        run(extend_hit_ungapped_approx_ncbi),
                        "round {round} {:?}",
                        (q_off, s_off, subject_len, x_dropoff, reduced_cutoff)
                    );
                    calls += 1;
                }
            }
        });
    }
    assert!(calls > 1_000_000, "{calls}");
}
