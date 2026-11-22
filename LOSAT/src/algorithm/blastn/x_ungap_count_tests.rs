use super::super::extension::extend_hit_ungapped_exact_ported;
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:197 and
// blast_stat.c:1124-1127
/// `-xdrop_ungap 1e10` (X-drop `INT_MIN`, which negates to itself): the left walk from the hit
/// at 3 passes the two sentinels (`INT4_MIN / 2` each: `INT_MIN`) and the mismatch A/C (-3),
/// where the `Int4` sum wraps to 2^31 - 3, above 0, as in NCBI's x86-64 build; the right walk
/// adds the two matches: `i32::MAX`. A build with overflow checks gives the same.
#[test]
fn a_walk_past_two_sentinels_wraps_as_ncbis_int4() {
    let matrix = crate::algorithm::blastn::alignment::build_blastna_matrix(1, -3);
    let query = [0u8, 15, 15, 0, 0];
    // C A A A A, packed (ncbi2na: A 0, C 1).
    let packed = [0b0100_0000u8, 0];
    let expected = (0, 0, 5, i32::MAX);
    let counted = extend_hit_ungapped_exact_counted(&query, &packed, 3, 3, 5, i32::MIN, &matrix);
    assert_eq!(key(counted), expected);
    let ported = extend_hit_ungapped_exact_ported(&query, &packed, 3, 3, 5, i32::MIN, &matrix);
    assert_eq!(key(ported), expected);
}

/// Every hit offset of random query blocks (contexts of 1 to 40 letters between sentinels,
/// a few ambiguity codes) against random and periodic subjects of 1 to 140 bases packed into
/// exactly `len.div_ceil(4)` bytes, at three X-drops: the counted loops give the ported
/// loops' four results.
#[test]
fn every_offset_matches_the_ported_loops() {
    let matrix = crate::algorithm::blastn::alignment::build_blastna_matrix(1, -3);
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut calls = 0u64;
    for round in 0..400 {
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
        let bases: Vec<u8> = (0..subject_len)
            .map(|i| {
                if round % 2 == 0 {
                    motif[i % period]
                } else {
                    (next(&mut state) % 4) as u8
                }
            })
            .collect();
        let mut packed = vec![0u8; subject_len.div_ceil(4)];
        for (i, b) in bases.iter().enumerate() {
            packed[i / 4] |= b << (2 * (3 - i % 4));
        }
        for q_off in 0..query.len() {
            for s_off in 0..subject_len {
                for x_dropoff in [5, 16, 40, i32::MIN] {
                    let args = (q_off, s_off, subject_len, x_dropoff);
                    assert_eq!(
                        key(extend_hit_ungapped_exact_counted(
                            &query,
                            &packed,
                            q_off,
                            s_off,
                            subject_len,
                            x_dropoff,
                            &matrix
                        )),
                        key(extend_hit_ungapped_exact_ported(
                            &query,
                            &packed,
                            q_off,
                            s_off,
                            subject_len,
                            x_dropoff,
                            &matrix
                        )),
                        "round {round} {args:?}"
                    );
                    calls += 1;
                }
            }
        }
    }
    assert!(calls > 100_000, "{calls}");
}

/// A seed of a MegaBLAST database index can start in its sequence's last 3 letters (the seed runs
/// on into the next sequence of the index's local id): for a sequence of 22 letters in 6 bytes
/// and a seed at 21, the approximate extension's blocks start at 24, past the last letter, and
/// NCBI's `Int4` length gives no block to the right (a `usize` one wrapped and read past the 6
/// bytes); its exact pass gives the exact extension's result.
#[test]
fn a_seed_in_the_last_three_letters_takes_no_block_to_the_right() {
    use super::super::extension::{build_nucl_score_table, extend_hit_ungapped_approx_ncbi};
    let matrix = crate::algorithm::blastn::alignment::build_blastna_matrix(1, -2);
    let score_table = build_nucl_score_table(1, -2);
    let mut state = 0x5eed_d5b2_u64;
    let letters: Vec<u8> = (0..22).map(|_| (next(&mut state) & 3) as u8).collect();
    let mut packed = vec![0u8; 6];
    for (p, &letter) in letters.iter().enumerate() {
        packed[p / 4] |= letter << (2 * (3 - p % 4));
    }
    // The query: the sequence's letters, then other letters (the next sequence's).
    let mut query = letters.clone();
    query.extend([1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]);
    for s_off in 17..22 {
        let exact =
            extend_hit_ungapped_exact_counted(&query, &packed, s_off, s_off, 22, 20, &matrix);
        let approx = |reduced_cutoff| {
            extend_hit_ungapped_approx_ncbi(
                &query,
                &packed,
                s_off,
                s_off,
                s_off + 12,
                22,
                20,
                &score_table,
                reduced_cutoff,
                &matrix,
            )
        };
        assert_eq!(key(approx(0)), key(exact));
        let blocks = approx(i32::MAX);
        assert!(blocks.score > 0 && blocks.q_start + blocks.length <= s_off + 12);
    }
}
