use super::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 17
    }
    fn below(&mut self, bound: u64) -> u64 {
        self.next() % bound
    }
}

/// The reference's statements for one hit (module doc): `None` for a flag-1 cell or a
/// candidate, else the cell's next `last_hit`. Written with the reference's subtraction order.
fn reference(cell: (i32, bool), s: i32, offset: i32, window: i32, word_size: i32) -> Option<i32> {
    let (last_hit, flag) = cell;
    if flag {
        return None;
    }
    let diff = s.wrapping_sub(last_hit.wrapping_sub(offset));
    if diff >= window {
        Some(s.wrapping_add(offset))
    } else if diff < word_size {
        Some(last_hit)
    } else {
        None
    }
}

// LOSAT_X_TBNDIAGPACK: the loop on NCBI's four-byte cells takes the same hits and leaves the same
// values as on `(last_hit, flag)` cells, on random tables of the values a packed table holds
// (`search_init::x_packs`: last_hit from 0 to 2^31 - 1, windows from wordsize + 1 to 2^28,
// offsets from the window to INT4_MAX / 4 + 2 x 5,000,000 + 2 x the window), the cells set
// next to each test's boundary as in the test below. LOSAT_X_TBNHITSTEP: `x_hits_steps` on
// either layout takes the same hits and leaves the same values as well.
#[test]
fn x_hits_bf_on_packed_cells_equals_pair_cells() {
    let mut rng = Rng(0x7b1f_5eed_a6a6_0002);
    let mut taken = [0u64; 3];
    for round in 0..4000 {
        let word_size = [0, 1, 2, 3, 3, 5][round % 6];
        let window = match round % 5 {
            0 => word_size + 1,
            1 | 3 => 40,
            2 => 1 << 28,
            _ => word_size + 1 + rng.below(1000) as i32,
        };
        assert!(fits(window, word_size as usize));
        let offset = match round % 4 {
            0 => window,
            1 => i32::MAX / 4 + 10_000_000 + 2 * window - rng.below(1 << 20) as i32,
            _ => window + rng.below((i32::MAX / 4) as u64) as i32,
        };
        let length = 1usize << rng.below(7);
        let pair: Vec<(i32, bool)> = (0..length)
            .map(|_| ((rng.next() as u32 & 0x7fff_ffff) as i32, rng.below(40) == 0))
            .collect();
        let mut packed: Vec<u32> = pair
            .iter()
            .map(|&(last_hit, flag)| {
                let mut cell = 0u32;
                cell.set(last_hit, flag);
                cell
            })
            .collect();
        let mut pair = pair;
        let n = rng.below(300) as usize;
        let mut pairs = Vec::with_capacity(n);
        for _ in 0..n {
            let q = rng.below(1 << 20) as u32;
            let s = rng.below(1 << 20) as u32;
            let cell = (q.wrapping_sub(s) as usize) & (length - 1);
            let s_diag = (s as i32) + offset;
            let diff = match rng.below(8) {
                0 => window,
                1 => window - 1,
                2 => word_size,
                3 => word_size - 1,
                _ => rng.next() as i32,
            };
            let last_hit = s_diag.wrapping_sub(diff);
            if rng.below(3) == 0 && last_hit >= 0 {
                pair[cell].0 = last_hit;
                packed[cell].set(last_hit, pair[cell].1);
            }
            pairs.push(BlastOffsetPair { q_off: q, s_off: s });
        }
        let mut pair_steps = pair.clone();
        let mut packed_steps = packed.clone();
        let mut done = 0;
        while done < pairs.len() {
            let a = x_hits_bf(&mut pair, offset, window, word_size, &pairs[done..]);
            let b = x_hits_bf(&mut packed, offset, window, word_size, &pairs[done..]);
            let c = x_hits_steps(&mut pair_steps, offset, window, word_size, &pairs[done..]);
            let d = x_hits_steps(&mut packed_steps, offset, window, word_size, &pairs[done..]);
            assert_eq!([a, a, a], [b, c, d], "round {round}: hits taken");
            let cells: Vec<(i32, bool)> = packed.iter().map(|cell| cell.get()).collect();
            assert_eq!(cells, pair, "round {round}: cells");
            assert_eq!(pair_steps, pair, "round {round}: cells of the steps");
            assert_eq!(
                packed_steps, packed,
                "round {round}: packed cells of the steps"
            );
            taken[0] += a as u64;
            done += a;
            if let Some(pair_at) = pairs.get(done) {
                let cell = (pair_at.q_off.wrapping_sub(pair_at.s_off) as usize) & (length - 1);
                taken[if pair[cell].1 { 1 } else { 2 }] += 1;
                // The per-hit body's outcome on the stop: a new hit on the cell.
                let next = (pair_at.s_off as i32) + offset;
                pair[cell] = (next, false);
                packed[cell].set(next, false);
                pair_steps[cell] = (next, false);
                packed_steps[cell].set(next, false);
                done += 1;
            }
        }
    }
    assert!(taken.iter().all(|&count| count > 0), "taken {taken:?}");
}

// The loop against the reference's statements on random tables: windows from wordsize + 1 to
// i32::MAX, offsets from 0 to i32::MAX (with wrapping sums), last_hit values at and next to each
// test's boundary (diff = window, window - 1, wordsize, wordsize - 1, and the wrap of
// diff - wordsize), flag-1 cells; every outcome reached, and the count it returns is the first
// hit the reference leaves to the per-hit body. `x_hits_steps` (LOSAT_X_TBNHITSTEP) on the same
// tables gives the same.
#[test]
fn x_hits_bf_takes_the_reference_outcomes() {
    let mut rng = Rng(0x7b1f_5eed_a6a6_0001);
    // too far, overlap, flag 1, candidate; wrapped diff - wordsize.
    let mut taken = [0u64; 5];
    // diff = window, window - 1, wordsize, wordsize - 1 on a flag-0 cell.
    let mut edges = [0u64; 4];
    for round in 0..4000 {
        let word_size = [0, 1, 2, 3, 3, 5][round % 6];
        let window = match round % 7 {
            0 => word_size + 1,
            1 | 4 => 40,
            2 => 64,
            3 => i32::MAX,
            5 => 1 << 20,
            _ => word_size + 1 + rng.below(1000) as i32,
        };
        assert!(fits(window, word_size as usize));
        let offset = match round % 5 {
            0 => window,
            1 => i32::MAX / 4 + rng.below(1 << 20) as i32,
            2 => i32::MAX - rng.below(1000) as i32,
            _ => rng.next() as i32 & i32::MAX,
        };
        let length = 1usize << rng.below(7);
        let mut table: Vec<(i32, bool)> = (0..length)
            .map(|_| (rng.next() as i32, rng.below(40) == 0))
            .collect();
        let n = rng.below(300) as usize;
        let mut pairs = Vec::with_capacity(n);
        for _ in 0..n {
            let q = rng.below(1 << 20) as u32;
            let s = rng.below(1 << 20) as u32;
            let cell = (q.wrapping_sub(s) as usize) & (length - 1);
            let s_diag = (s as i32).wrapping_add(offset);
            // Set the cell so that diff lands on or next to a boundary.
            let diff = match rng.below(10) {
                0 => window,
                1 => window.wrapping_sub(1),
                2 => word_size,
                3 => word_size - 1,
                4 => i32::MIN
                    .wrapping_add(word_size)
                    .wrapping_sub(1 + rng.below(3) as i32),
                5 => i32::MIN,
                _ => rng.next() as i32,
            };
            if rng.below(3) == 0 {
                table[cell].0 = s_diag.wrapping_sub(diff);
            }
            pairs.push(BlastOffsetPair { q_off: q, s_off: s });
        }
        let mut expected = table.clone();
        let mut stop = n;
        for (i, pair) in pairs.iter().enumerate() {
            let cell = (pair.q_off.wrapping_sub(pair.s_off) as usize) & (length - 1);
            let s = pair.s_off as i32;
            let entry = expected[cell];
            let diff = s.wrapping_sub(entry.0.wrapping_sub(offset));
            if !entry.1 {
                for (edge, value) in [window, window - 1, word_size, word_size - 1]
                    .into_iter()
                    .enumerate()
                {
                    edges[edge] += u64::from(diff == value);
                }
            }
            match reference(entry, s, offset, window, word_size) {
                None => {
                    taken[if entry.1 { 2 } else { 3 }] += 1;
                    stop = i;
                    break;
                }
                Some(next) => {
                    taken[if diff >= window { 0 } else { 1 }] += 1;
                    if diff < word_size && diff.wrapping_sub(word_size) >= 0 {
                        taken[4] += 1;
                    }
                    expected[cell].0 = next;
                }
            }
        }
        let mut steps = table.clone();
        let got = x_hits_bf(&mut table, offset, window, word_size, &pairs);
        assert_eq!(got, stop, "round {round}: stop");
        assert_eq!(table, expected, "round {round}: table");
        let got = x_hits_steps(&mut steps, offset, window, word_size, &pairs);
        assert_eq!(got, stop, "round {round}: stop of the steps");
        assert_eq!(steps, expected, "round {round}: table of the steps");
    }
    assert!(taken.iter().all(|&count| count > 0), "outcomes {taken:?}");
    assert!(edges.iter().all(|&count| count > 0), "boundaries {edges:?}");
}

#[test]
fn fits_needs_a_window_above_the_word_size() {
    assert!(fits(40, 3));
    assert!(fits(1, 0));
    assert!(fits(i32::MAX, 3));
    assert!(!fits(3, 3));
    assert!(!fits(2, 3));
    assert!(!fits(40, i32::MAX as usize + 1));
}

#[test]
fn an_empty_table_takes_no_hit() {
    let pairs = [BlastOffsetPair { q_off: 1, s_off: 0 }];
    assert_eq!(x_hits_bf::<(i32, bool)>(&mut [], 40, 40, 3, &pairs), 0);
    assert_eq!(x_hits_bf::<u32>(&mut [], 40, 40, 3, &pairs), 0);
    assert_eq!(x_hits_steps::<(i32, bool)>(&mut [], 40, 40, 3, &pairs), 0);
    assert_eq!(x_hits_steps::<u32>(&mut [], 40, 40, 3, &pairs), 0);
}
