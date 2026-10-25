//! The two-hit loop of tblastn's unit path (`XTwoHit::pairs`)
//! with its two common outcomes, too far and overlap, taken without a branch.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-551
//! ```c
//!             diag_coord = (query_offset - subject_offset) & diag_mask;
//!
//!             /* If the reset bit is set, an extension just happened. */
//!             if (diag_array[diag_coord].flag) {
//!                 /* If we've already extended past this hit, skip it. */
//!                 if ((Int4) (subject_offset + diag_offset) <
//!                     diag_array[diag_coord].last_hit) {
//!                     continue;
//!                 }
//!                 /* Otherwise, start a new hit. */
//!                 else {
//!                     diag_array[diag_coord].last_hit =
//!                         subject_offset + diag_offset;
//!                     diag_array[diag_coord].flag = 0;
//!                 }
//!             }
//!             /* If the reset bit is cleared, try to start an extension. */
//!             else {
//!                 /* find the distance to the last hit on this diagonal */
//!                 last_hit = diag_array[diag_coord].last_hit - diag_offset;
//!                 diff = subject_offset - last_hit;
//!
//!                 if (diff >= window) {
//!                     /* We are beyond the window for this diagonal; start a
//!                        new hit */
//!                     diag_array[diag_coord].last_hit =
//!                         subject_offset + diag_offset;
//!                     continue;
//!                 }
//!
//!                 /* If the difference is less than the wordsize (i.e. last
//!                    hit and this hit overlap), give up */
//!
//!                 if (diff < wordsize) {
//!                     continue;
//!                 }
//! ```
//! On a genome searched with proteins most hits end in one of the two `continue`s of a flag-0 cell:
//! too far (83-87 % of the hits of the four ledger cases, v4-A10 probe 1) or overlap (12-15 %).
//! Which one is hard to predict: taking both without a branch removed 0.12-0.17 mispredicted
//! branches per hit from the whole run (0.25-0.35 before). A flag-1 cell (0.5-0.7 %) or a hit past
//! both tests, a candidate for the extension (1.2-1.6 %), is rare.
//!
//! `x_hits_bf` reads the cell of each hit and computes `d = (s + diag_offset) - last_hit`. That is
//! the reference's `diff = s - (last_hit - diag_offset)` in wrapping `i32` operations, so in a release
//! build (no overflow checks) the same value whatever the operands; and the operands stay far from
//! overflow (`Diagonals::finish_frame` starts again from `window` once the offset reaches
//! INT4_MAX / 4). A candidate (`wordsize <= d < window`) is
//! `(d - wordsize) as u32 < window - wordsize`, one unsigned compare: for `d >= window` the
//! difference is at least `window - wordsize`; for `d < wordsize` it is negative, or, when it wraps
//! (`d < i32::MIN + wordsize`), at least `2^31 - wordsize`, which is not below
//! `window - wordsize` because `window <= i32::MAX`. So the test is exact whenever
//! `0 <= wordsize < window` (`fits`). A flag-1 cell or a candidate stops the loop at that hit; the
//! caller (`XTwoHit::x_pairs_bf`) runs the body of `pairs` on that one hit and calls the loop again
//! after it. Every other hit is a flag-0 cell whose next `last_hit` is selected without a branch:
//! `s + diag_offset` when too far (`d >= window`), the value it read when it overlaps. The cell's
//! flag stays 0. So every cell holds after each hit what the reference leaves there, and the
//! reference's body sees the same table for every hit it runs and saves the same init HSPs in the
//! same order.
//!
//! The cell index is the hit's diagonal masked by the table's length minus one. The table's length
//! is a power of two and its `mask` is that length minus one (`Diagonals::new`), so this is the
//! reference's index, and it is never past the end of the table.
//!
//! The loop runs for every offset array of the unit path (the former switch `LOSAT_X_TBNHITBF`
//! graduated on 2026-10-26). TBLASTX's unit scan has a loop of this form under its own switch,
//! LOSAT_X_TBXHITBF.

use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair;

/// Whether `x_hits_bf` gives the reference's outcomes for this
/// window and word size (module doc): `0 <= wordsize < window`.
pub(super) fn fits(window: i32, word_size: usize) -> bool {
    i32::try_from(word_size).is_ok_and(|word_size| word_size < window)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-552 (module doc).
/// The hits of `pairs` from the first, each on its cell of
/// `entries` (`(last_hit, flag)`), up to the first that is a flag-1 cell or a candidate for the
/// extension; returns the number of hits taken (`pairs.len()` when there is no such hit). The
/// hits taken are too far (`last_hit = s + offset`) or overlap (no change). The outcomes are
/// exact when `fits(window, word_size)` (module doc); `entries.len()` is a power of two.
#[inline(always)]
pub(super) fn x_hits_bf(
    entries: &mut [(i32, bool)],
    offset: i32,
    window: i32,
    word_size: i32,
    pairs: &[BlastOffsetPair],
) -> usize {
    debug_assert!(entries.is_empty() || entries.len().is_power_of_two());
    debug_assert!(0 <= word_size && word_size < window);
    let Some(last) = entries.len().checked_sub(1) else {
        return 0;
    };
    // [C] diag_coord = (query_offset - subject_offset) & diag_mask; (diag_mask = length - 1)
    let mask = last as u32;
    // window - wordsize: the length of the candidate range of diff - wordsize.
    let k = window.wrapping_sub(word_size) as u32;
    let base = entries.as_mut_ptr();
    let start = pairs.as_ptr();
    let mut p = start;
    // SAFETY: one past the last element of `pairs`.
    let end = unsafe { start.add(pairs.len()) };
    while p < end {
        // SAFETY: start <= p < end.
        let pair = unsafe { p.read() };
        let index = (pair.q_off.wrapping_sub(pair.s_off) & mask) as usize;
        // SAFETY: index <= mask = entries.len() - 1.
        let entry = unsafe { base.add(index) };
        // SAFETY: as above.
        let (last_hit, flag) = unsafe { entry.read() };
        // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
        // [C] diff = subject_offset - last_hit;
        let s_diag = (pair.s_off as i32).wrapping_add(offset);
        let diff = s_diag.wrapping_sub(last_hit);
        let candidate = (diff.wrapping_sub(word_size) as u32) < k;
        // [C] if (diag_array[diag_coord].flag), or a hit past both tests.
        if flag | candidate {
            break;
        }
        // [C] if (diff >= window) last_hit = subject_offset + diag_offset;
        // [C] if (diff < wordsize) continue;
        let next = if diff >= window { s_diag } else { last_hit };
        // SAFETY: as above.
        unsafe { (*entry).0 = next };
        // SAFETY: p < end, so p + 1 <= end.
        p = unsafe { p.add(1) };
    }
    // SAFETY: start <= p <= end, both in `pairs`.
    unsafe { p.offset_from(start) as usize }
}

#[cfg(test)]
mod tests {
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
    fn reference(
        cell: (i32, bool),
        s: i32,
        offset: i32,
        window: i32,
        word_size: i32,
    ) -> Option<i32> {
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

    // The loop against the reference's statements on random tables: windows from wordsize + 1 to
    // i32::MAX, offsets from 0 to i32::MAX (with wrapping sums), last_hit values at and next to each
    // test's boundary (diff = window, window - 1, wordsize, wordsize - 1, and the wrap of
    // diff - wordsize), flag-1 cells; every outcome reached, and the count it returns is the first
    // hit the reference leaves to the per-hit body.
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
            let got = x_hits_bf(&mut table, offset, window, word_size, &pairs);
            assert_eq!(got, stop, "round {round}: stop");
            assert_eq!(table, expected, "round {round}: table");
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
        assert_eq!(x_hits_bf(&mut [], 40, 40, 3, &pairs), 0);
    }
}
