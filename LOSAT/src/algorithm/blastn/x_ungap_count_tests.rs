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
/// at 3 passes the two sentinels (half of `INT4_MIN` each: `INT_MIN`) and the mismatch A/C (-3),
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
