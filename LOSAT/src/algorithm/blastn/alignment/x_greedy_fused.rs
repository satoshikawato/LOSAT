//! No NCBI counterpart: the loop over distances of
//! the score-only `BLAST_GreedyAlign` on an unpacked subject in one function, each distance in one
//! pass over its diagonals.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:509-537,680-686
//!
//! Where it runs. `blast_greedy_align` (the port of `BLAST_GreedyAlign`) calls
//! `x_greedy_fused_dispatch` after its set-up: `max_score[0..xdrop_offset]` cleared, the first
//! `s_FindFirstMismatch`, the early return, `last_seq2_off[0][diag_origin]`, `max_score[0]` and the
//! first diagonal bounds (greedy_align.c:463-505). The call takes over the rest of the function
//! when all of these hold: there is no edit block (score only, the preliminary stage; the subject
//! is always unpacked, `x_greedy_prelim.rs`), the fence was not hit by the first scan,
//! `1 <= max_dist`, and `len1` / `len2` lie within the slices. Otherwise `x_greedy_row` runs the
//! distance (a traceback call, or a call this module declines), and the ported loop only for
//! lengths that do not fit the slices. The affine alignment keeps its own path.
//!
//! Why the result is the same. The function is the C loop statement by statement, on the same
//! memory:
//!
//! - The same two rows (`last_seq2_off[0]` / `[1]` of `GreedyNonAffineMem`, origin 0, cell `k` =
//!   diagonal `k`; distance `d` writes row `d & 1` and reads row `(d - 1) & 1`, as
//!   `last_seq2_off[d + 1] = last_seq2_off[d - 1]` makes the C rows roll) and the same
//!   `max_score` array. Every cell the C loop writes is written here with the same value, and no
//!   other cell: the four sentinels of distance `d - 1`, `kInvalidOffset` for a dropped diagonal
//!   that is not the lower bound (a dropped lower bound is not written), the slide end for a kept
//!   diagonal, and `max_score[d]`. Stale cells of earlier calls therefore stay as they are.
//! - The same order. The diagonals of a distance run in increasing `k` with the C tests in the C
//!   order; the state the loop carries from diagonal to diagonal (`diag_lower`, `diag_upper`,
//!   `curr_extent` and `curr_diag`, the seed and `longest_match_run`, the end flags) is kept in
//!   locals and changed where the C code changes it. `curr_seq2_index` is not carried: it is the
//!   value stored in the cell of `curr_diag` when that diagonal raised `curr_extent`, and a
//!   distance writes each of its cells once, so it is read from that cell after the diagonals
//!   (no diagonal raised the extent while `curr_diag` is 0, as no diagonal is 0). The three cells
//!   of distance `d - 1` a diagonal reads are kept in a window that moves by one cell per
//!   diagonal; distance `d` never writes the row of `d - 1`, which is the same as reading the
//!   cells again.
//! - The same slide. `s_FindFirstMismatch` with a `rem` of 4 continues while
//!   `seq1[i] < 4 && seq1[i] == seq2[i]`; per byte that is `(b1 ^ b2) | (b1 & 0xFC) == 0`, so the
//!   first non-zero byte of `(w1 ^ w2) | (w1 & 0xFC..FC)` over eight positions (lowest byte first
//!   to the right, highest byte first to the left) is the first mismatch. The eight-byte
//!   comparison runs when both sequences have at least eight positions left; a slide that stops
//!   within them stops inside both sequences, so neither end test can fire and the stop byte of
//!   `seq2` is the byte the C fence test reads; when that byte is not the sentry and the slide is
//!   not longer than `longest_match_run`, the C code only sets `diag_upper`, the cell and the
//!   extent for the diagonal. Every other kept diagonal (fewer than eight positions left, eight
//!   or more matches, the sentry at the stop, a new seed) takes `x_cell_cold`: the byte loop
//!   (`x_first_mismatch`, the same predicate eight bytes at a time and then byte by byte), the C
//!   fence test with its bound and the C updates in order. The exception is a diagonal past the
//!   end of `seq1` (`x_past_seq1_end`), which is handled in place.
//! - The same integers. `seq2_index`, `seq1_index`, the extent and the scores are the C
//!   expressions on `i32`. `xdrop_score` is the reference's
//!   `((raw as f64) / (match_cost as f64 / 2.0)).ceil() as i32`, except for `match_cost == 2`
//!   (megablast's reward 1, doubled by `BLAST_AffineGreedyAlign`): the divisor is then 1.0, and
//!   an `i32` divided by 1.0 and rounded up is itself (exact in `f64`), so `raw` is used as it is.
//! - The same exits, with the same values in the outputs: converged -> `best_dist`; past
//!   `max_dist` -> -1; a fence at the stop byte of a slide -> `fence_hit` and 0, at once (the cell
//!   of that diagonal and the seed are not written, as in C). The outputs `seq1_align_len`,
//!   `seq2_align_len` and the seed keep the values of their last update in every case.

#![allow(clippy::too_many_arguments)]

use super::{x_first_mismatch, x_past_seq1_end, GreedyNonAffineMem, GreedySeed, FENCE_SENTRY};
use crate::algorithm::blastn::constants::INVALID_OFFSET;
#[cfg(test)]
use std::cell::Cell;
#[cfg(target_arch = "x86_64")]
use std::sync::OnceLock;

#[cfg(test)]
thread_local! {
    /// Lets a test run the plain build on a CPU with AVX2.
    static X_GREEDY_FUSE_TEST_PLAIN: Cell<bool> = const { Cell::new(false) };
    /// Diagonals that took `x_cell_cold` with fewer than eight positions left, and with a slide
    /// of eight or more.
    static X_GREEDY_FUSE_TEST_REST: Cell<(u64, u64)> = const { Cell::new((0, 0)) };
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:380-387
/// The parameters the distance loop reads, and the values of the set-up it starts from
/// (`diag_lower`, `diag_upper`, `xdrop_offset`).
#[derive(Clone, Copy)]
struct XFuseInput {
    seq1: *const u8,
    len1: i32,
    seq2: *const u8,
    len2: i32,
    xdrop_threshold: i32,
    match_cost: i32,
    mismatch_cost: i32,
    max_dist: i32,
    diag_lower: i32,
    diag_upper: i32,
    max_score_offset: usize,
}

/// NCBI reference (598d8ae6): c++/include/algo/blast/core/greedy_align.h:88-99
/// The two rolling rows of `last_seq2_off` (origin 0) and `aux_data->max_score`.
#[derive(Clone, Copy)]
struct XFuseMem {
    rows: [*mut i32; 2],
    max_score: *mut i32,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:463-470
/// The outputs of `BLAST_GreedyAlign` besides its return value: what the call had set when the
/// loop starts, and what the loop leaves. `seed_match_length` is `longest_match_run`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct XFuseOut {
    seq1_align_len: i32,
    seq2_align_len: i32,
    seed_start_q: i32,
    seed_start_s: i32,
    seed_match_length: i32,
    fence_hit: bool,
}

impl XFuseOut {
    fn read(seq1_align_len: i32, seq2_align_len: i32, seed: &GreedySeed, fence_hit: bool) -> Self {
        Self {
            seq1_align_len,
            seq2_align_len,
            seed_start_q: seed.start_q,
            seed_start_s: seed.start_s,
            seed_match_length: seed.match_length,
            fence_hit,
        }
    }

    fn write(
        &self,
        seq1_align_len: &mut i32,
        seq2_align_len: &mut i32,
        seed: &mut GreedySeed,
        fence_hit: &mut bool,
    ) {
        *seq1_align_len = self.seq1_align_len;
        *seq2_align_len = self.seq2_align_len;
        seed.start_q = self.seed_start_q;
        seed.start_s = self.seed_start_s;
        seed.match_length = self.seed_match_length;
        *fence_hit = self.fence_hit;
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:394-404,537-541
/// The variables the C loop reads or changes on most diagonals. `x_fused_loop` keeps them in
/// locals; `x_cell_cold` takes and returns them. `curr_seq2_index` is not among them: it is
/// the cell of `curr_diag` (module doc).
#[derive(Clone, Copy)]
struct XHot {
    diag_lower: i32,
    diag_upper: i32,
    curr_extent: i32,
    curr_diag: i32,
    longest_match_run: i32,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:578-582,607-614
/// The variables only rare diagonals change (`x_cell_cold`): the end flags, the seed start and
/// the fence.
struct XColdState {
    end1_reached: bool,
    end2_reached: bool,
    seed_start_q: i32,
    seed_start_s: i32,
    fence_hit: bool,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:565-614
/// A kept diagonal, as the C loop handles it after its X-drop test: the slide by
/// `s_FindFirstMismatch` for a `rem` of 4 (`x_first_mismatch`: the same predicate, eight bytes at a
/// time, then byte by byte), its fence test with its bound, and every update in the C order. The
/// diagonals whose slide one eight-byte comparison does not settle (fewer than eight positions
/// left, eight or more matches, the fence at the stop, a longer run than the seed's), except a
/// diagonal past the end of `seq1` (`x_past_seq1_end`), which is handled in place. With a
/// fence, `cold.fence_hit` is set and nothing else changes (the C code returns at once).
///
/// SAFETY: `seq1` / `seq2` readable for `len1` / `len2` bytes, `current` writable at `k`.
#[cold]
#[inline(never)]
unsafe fn x_cell_cold<const REVERSE: bool>(
    seq1: *const u8,
    seq2: *const u8,
    len1: i32,
    len2: i32,
    current: *mut i32,
    k: i32,
    seq1_index: i32,
    seq2_index: i32,
    hot: XHot,
    cold: &mut XColdState,
) -> XHot {
    // The reference indexes the slices with these offsets.
    assert!(seq1_index >= 0 && seq2_index >= 0);
    let mut hot = hot;
    hot.diag_upper = k;
    let index = x_first_mismatch::<REVERSE>(seq1, seq2, len1, len2, seq1_index, seq2_index);
    #[cfg(test)]
    X_GREEDY_FUSE_TEST_REST.with(|n| {
        let (short, long) = n.get();
        if (len1 - seq1_index).min(len2 - seq2_index) < 8 {
            n.set((short + 1, long));
        } else if index >= 8 {
            n.set((short, long + 1));
        }
    });
    let stop1 = seq1_index + index;
    let stop2 = seq2_index + index;
    if stop2 < len2 {
        let at = if REVERSE { len2 - 1 - stop2 } else { stop2 };
        if *seq2.add(at as usize) == FENCE_SENTRY {
            cold.fence_hit = true;
            return hot;
        }
    }
    if index > hot.longest_match_run {
        cold.seed_start_q = seq1_index;
        cold.seed_start_s = seq2_index;
        hot.longest_match_run = index;
    }
    *current.add(k as usize) = stop2;
    if stop1 + stop2 > hot.curr_extent {
        hot.curr_extent = stop1 + stop2;
        hot.curr_diag = k;
    }
    if stop2 == len2 {
        hot.diag_lower = k + 1;
        cold.end2_reached = true;
    }
    if stop1 == len1 {
        hot.diag_upper = k - 1;
        cold.end1_reached = true;
    }
    hot
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:509-686
/// The distance loop of `BLAST_GreedyAlign` without an edit block, for a `rem` of 4, from the
/// state the set-up left (module doc). The parts below say which C lines they port. `WIDE`:
/// both sequences have at least eight bytes, so the eight-byte comparison may run; without it
/// every kept diagonal takes `x_cell_cold`.
///
/// SAFETY: `seq1` / `seq2` readable for `len1` / `len2` bytes (at least eight each with `WIDE`);
/// both rows writable for `2 * max_dist + 6` cells and `max_score` for
/// `max_dist + max_score_offset + 1` cells; the rows and `max_score` do not overlap;
/// `diag_lower` / `diag_upper` are `max_dist + 1` / `max_dist + 3`.
#[inline(always)]
unsafe fn x_fused_loop<const REVERSE: bool, const WIDE: bool>(
    input: &XFuseInput,
    mem: XFuseMem,
    out: &mut XFuseOut,
) -> i32 {
    const NOT_A_BASE: u64 = 0xFCFC_FCFC_FCFC_FCFC;
    let XFuseInput {
        seq1,
        len1,
        seq2,
        len2,
        xdrop_threshold,
        match_cost,
        mismatch_cost,
        max_dist,
        max_score_offset,
        ..
    } = *input;
    let diag_origin = max_dist + 2;
    let op_cost = match_cost + mismatch_cost;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:496
    // `max_score_raw[i]` is `aux_data->max_score[i]`, `max_score[d]` the C `max_score[d]`.
    let max_score_raw = mem.max_score;
    let max_score = mem.max_score.add(max_score_offset);
    // `xdrop_score` as an integer when the divisor is 1.0 (module doc).
    let unit_divisor = match_cost == 2;
    let divisor = match_cost as f64 / 2.0;
    // A diagonal has eight positions left in both sequences when these are >= its offsets.
    let last1 = len1 - 8;
    let last2 = len2 - 8;
    // Row of distance d - 1 and row of distance d.
    let (mut previous, mut current) = (mem.rows[0], mem.rows[1]);

    let mut hot = XHot {
        diag_lower: input.diag_lower,
        diag_upper: input.diag_upper,
        curr_extent: 0,
        curr_diag: 0,
        longest_match_run: out.seed_match_length,
    };
    let mut cold = XColdState {
        end1_reached: false,
        end2_reached: false,
        seed_start_q: out.seed_start_q,
        seed_start_s: out.seed_start_s,
        fence_hit: false,
    };
    let mut seq1_align_len = out.seq1_align_len;
    let mut seq2_align_len = out.seq2_align_len;
    let mut best_dist = 0;
    let mut converged = false;

    let mut d = 1;
    'distances: while d <= max_dist {
        debug_assert!(hot.diag_lower >= diag_origin - d && hot.diag_upper <= diag_origin + d);
        debug_assert!(hot.diag_lower <= hot.diag_upper);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:523-526
        *previous.add((hot.diag_lower - 1) as usize) = INVALID_OFFSET;
        *previous.add(hot.diag_lower as usize) = INVALID_OFFSET;
        *previous.add(hot.diag_upper as usize) = INVALID_OFFSET;
        *previous.add((hot.diag_upper + 1) as usize) = INVALID_OFFSET;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:531-533
        let raw = *max_score_raw.add(d as usize) + op_cost * d - xdrop_threshold;
        let xdrop_score = if unit_divisor {
            raw
        } else {
            ((raw as f64) / divisor).ceil() as i32
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:512-518
        let tmp_diag_lower = hot.diag_lower;
        let tmp_diag_upper = hot.diag_upper;
        let mut diag_lower = hot.diag_lower;
        let mut diag_upper = hot.diag_upper;
        let mut curr_extent: i32 = 0;
        let mut curr_diag: i32 = 0;
        let mut longest_match_run = hot.longest_match_run;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-553
        // `left`, `middle` = cells k - 1, k of distance d - 1; cell k + 1 is read per diagonal.
        let mut left = *previous.add((tmp_diag_lower - 1) as usize);
        let mut middle = *previous.add(tmp_diag_lower as usize);
        let mut k = tmp_diag_lower;
        while k <= tmp_diag_upper {
            let right = *previous.add((k + 1) as usize);
            // The second `max` as `a + max(b - a, 0)` with a sign mask, without a branch (which
            // cell wins varies per diagonal); the offsets are far from the `i32` limits.
            let from_upper = right.max(middle) + 1;
            let gain = left - from_upper;
            let seq2_index = from_upper + (gain & !(gain >> 31));
            left = middle;
            middle = right;
            let seq1_index = seq2_index + k - diag_origin;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:553-565
            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == diag_lower {
                    diag_lower += 1;
                } else {
                    *current.add(k as usize) = INVALID_OFFSET;
                }
                k += 1;
                continue;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:565-614
            // One eight-byte comparison (module doc) settles the slide when both sequences have
            // eight positions left and it stops within them. The C code then runs, for this
            // diagonal: `diag_upper = k`; no fence (the stop byte of `seq2` is not the sentry);
            // no seed (`index <= longest_match_run`); the cell; the extent; no end test fires
            // (the stop lies inside both sequences). Every other kept diagonal is `x_cell_cold`,
            // except a diagonal past the end of `seq1` (`x_past_seq1_end`), which is handled in
            // place.
            if WIDE && ((last1 - seq1_index) | (last2 - seq2_index) | seq1_index) >= 0 {
                let (w1, w2) = if REVERSE {
                    (
                        u64::from_le(
                            (seq1.add((last1 - seq1_index) as usize) as *const u64)
                                .read_unaligned(),
                        ),
                        u64::from_le(
                            (seq2.add((last2 - seq2_index) as usize) as *const u64)
                                .read_unaligned(),
                        ),
                    )
                } else {
                    (
                        u64::from_le(
                            (seq1.add(seq1_index as usize) as *const u64).read_unaligned(),
                        ),
                        u64::from_le(
                            (seq2.add(seq2_index as usize) as *const u64).read_unaligned(),
                        ),
                    )
                };
                let diff = (w1 ^ w2) | (w1 & NOT_A_BASE);
                let bits = if REVERSE {
                    diff.leading_zeros()
                } else {
                    diff.trailing_zeros()
                };
                let index = (bits >> 3) as i32;
                let stop_byte = if REVERSE {
                    (w2 >> (56 - (bits & 0x38))) as u8
                } else {
                    (w2 >> (bits & 0x38)) as u8
                };
                if diff != 0 && stop_byte != FENCE_SENTRY && index <= longest_match_run {
                    diag_upper = k;
                    let stop2 = seq2_index + index;
                    *current.add(k as usize) = stop2;
                    let extent = seq1_index + index + stop2;
                    if extent > curr_extent {
                        curr_extent = extent;
                        curr_diag = k;
                    }
                    k += 1;
                    continue;
                }
            }
            // A diagonal past the end of `seq1` (`x_past_seq1_end`).
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:565-600
            if x_past_seq1_end::<REVERSE>(seq2, len1, len2, seq1_index, seq2_index) {
                diag_upper = k;
                *current.add(k as usize) = seq2_index;
                let extent = seq1_index + seq2_index;
                if extent > curr_extent {
                    curr_extent = extent;
                    curr_diag = k;
                }
                k += 1;
                continue;
            }
            let next = x_cell_cold::<REVERSE>(
                seq1,
                seq2,
                len1,
                len2,
                current,
                k,
                seq1_index,
                seq2_index,
                XHot {
                    diag_lower,
                    diag_upper,
                    curr_extent,
                    curr_diag,
                    longest_match_run,
                },
                &mut cold,
            );
            if cold.fence_hit {
                // The C code returns here; the seed it leaves includes this distance's updates.
                hot = next;
                break 'distances;
            }
            diag_lower = next.diag_lower;
            diag_upper = next.diag_upper;
            curr_extent = next.curr_extent;
            curr_diag = next.curr_diag;
            longest_match_run = next.longest_match_run;
            k += 1;
        }
        hot = XHot {
            diag_lower,
            diag_upper,
            curr_extent,
            curr_diag,
            longest_match_run,
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:596-600,619-634
        // `curr_seq2_index` is the value stored in the cell of `curr_diag` when it raised
        // `curr_extent`; this distance writes each cell once, so the cell still holds it. No
        // diagonal raised it when `curr_diag` is still 0 (every diagonal is at least 1).
        let curr_seq2_index = if hot.curr_diag != 0 {
            *current.add(hot.curr_diag as usize)
        } else {
            0
        };
        let curr_score = hot.curr_extent * (match_cost / 2) - d * op_cost;
        let previous_max = *max_score.add((d - 1) as usize);
        if curr_score >= previous_max {
            *max_score.add(d as usize) = curr_score;
            best_dist = d;
            seq2_align_len = curr_seq2_index;
            seq1_align_len = curr_seq2_index + hot.curr_diag - diag_origin;
        } else {
            *max_score.add(d as usize) = previous_max;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:639-662
        if hot.diag_lower > hot.diag_upper {
            converged = true;
            break;
        }
        if !cold.end2_reached {
            hot.diag_lower -= 1;
        }
        if !cold.end1_reached {
            hot.diag_upper += 1;
        }
        std::mem::swap(&mut previous, &mut current);
        d += 1;
    }

    out.seq1_align_len = seq1_align_len;
    out.seq2_align_len = seq2_align_len;
    out.seed_start_q = cold.seed_start_q;
    out.seed_start_s = cold.seed_start_s;
    out.seed_match_length = hot.longest_match_run;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:574-576,682-686
    if cold.fence_hit {
        out.fence_hit = true;
        return 0;
    }
    if !converged {
        return -1;
    }
    best_dist
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:509
/// The same integer code as `x_fused_plain`, so the same result as the C loop; the compiler may
/// only use other instructions (`tzcnt` / `lzcnt`).
///
/// SAFETY: as `x_fused_loop`, and the CPU must support the features.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,bmi1,bmi2,lzcnt")]
unsafe fn x_fused_avx2(
    reverse: bool,
    input: &XFuseInput,
    mem: XFuseMem,
    out: &mut XFuseOut,
) -> i32 {
    match (reverse, input.len1 >= 8 && input.len2 >= 8) {
        (true, true) => x_fused_loop::<true, true>(input, mem, out),
        (true, false) => x_fused_loop::<true, false>(input, mem, out),
        (false, true) => x_fused_loop::<false, true>(input, mem, out),
        (false, false) => x_fused_loop::<false, false>(input, mem, out),
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:509
/// The build for any CPU; `reverse` chooses the instance (the C `reverse` argument).
///
/// SAFETY: as `x_fused_loop`.
#[inline(never)]
unsafe fn x_fused_plain(
    reverse: bool,
    input: &XFuseInput,
    mem: XFuseMem,
    out: &mut XFuseOut,
) -> i32 {
    match (reverse, input.len1 >= 8 && input.len2 >= 8) {
        (true, true) => x_fused_loop::<true, true>(input, mem, out),
        (true, false) => x_fused_loop::<true, false>(input, mem, out),
        (false, true) => x_fused_loop::<false, true>(input, mem, out),
        (false, false) => x_fused_loop::<false, false>(input, mem, out),
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:509
/// Picks the AVX2 or the plain build of the same loop (the CPU is checked once, the choice is
/// made once per call). Both give the result of the C loop; this is a choice of machine code.
///
/// SAFETY: as `x_fused_loop`.
unsafe fn x_fused_run(reverse: bool, input: &XFuseInput, mem: XFuseMem, out: &mut XFuseOut) -> i32 {
    #[cfg(target_arch = "x86_64")]
    {
        static AVX2: OnceLock<bool> = OnceLock::new();
        #[cfg(test)]
        let plain = X_GREEDY_FUSE_TEST_PLAIN.with(Cell::get);
        #[cfg(not(test))]
        let plain = false;
        if !plain
            && *AVX2.get_or_init(|| {
                std::is_x86_feature_detected!("avx2")
                    && std::is_x86_feature_detected!("bmi1")
                    && std::is_x86_feature_detected!("bmi2")
                    && std::is_x86_feature_detected!("lzcnt")
            })
        {
            return x_fused_avx2(reverse, input, mem, out);
        }
    }
    x_fused_plain(reverse, input, mem, out)
}

/// Runs the distance loop of this module on `non_affine_mem` (checked sizes, see `x_fused_loop`).
fn x_fused_on(
    reverse: bool,
    input: &XFuseInput,
    last_seq2_off: &mut [Vec<i32>; 2],
    max_score: &mut [i32],
    out: &mut XFuseOut,
) -> i32 {
    let base_row_len = (2 * input.max_dist + 6) as usize;
    let max_score_len = input.max_dist as usize + input.max_score_offset + 1;
    let diag_origin = input.max_dist + 2;
    assert!(
        last_seq2_off[0].len() >= base_row_len
            && last_seq2_off[1].len() >= base_row_len
            && max_score.len() >= max_score_len
            && input.diag_lower == diag_origin - 1
            && input.diag_upper == diag_origin + 1,
        "greedy scratch smaller than the set-up allocated"
    );
    let (row0, row1) = last_seq2_off.split_at_mut(1);
    let mem = XFuseMem {
        rows: [row0[0].as_mut_ptr(), row1[0].as_mut_ptr()],
        max_score: max_score.as_mut_ptr(),
    };
    // SAFETY: the dispatch checked `len1` / `len2` against the slices; the sizes and the first
    // bounds are checked above; the band then stays within `diag_origin - d ..= diag_origin + d`
    // (each bound moves outwards by at most one per distance), whose cells and neighbours
    // `1 ..= 2 * max_dist + 3` lie in the rows.
    unsafe { x_fused_run(reverse, input, mem, out) }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:498-509
/// The fused loop in `blast_greedy_align`, after the lines above. `None`: the loop of
/// `blast_greedy_align` runs (a case this module does not take, module doc). `Some(value)`: the call's return value; the
/// outputs, rows and `max_score` hold what the ported loop leaves.
#[inline]
pub(super) fn x_greedy_fused_dispatch(
    seq1: &[u8],
    len1: i32,
    seq2: &[u8],
    len2: i32,
    reverse: bool,
    xdrop_threshold: i32,
    match_cost: i32,
    mismatch_cost: i32,
    seq1_align_len: &mut i32,
    seq2_align_len: &mut i32,
    max_dist: i32,
    non_affine_mem: &mut GreedyNonAffineMem,
    score_only: bool,
    fence_hit: &mut bool,
    seed: &mut GreedySeed,
    diag_lower: i32,
    diag_upper: i32,
    max_score_offset: usize,
) -> Option<i32> {
    if !score_only
        || *fence_hit
        || max_dist < 1
        || len1 < 0
        || len2 < 0
        || len1 as usize > seq1.len()
        || len2 as usize > seq2.len()
    {
        return None;
    }
    let input = XFuseInput {
        seq1: seq1.as_ptr(),
        len1,
        seq2: seq2.as_ptr(),
        len2,
        xdrop_threshold,
        match_cost,
        mismatch_cost,
        max_dist,
        diag_lower,
        diag_upper,
        max_score_offset,
    };
    let mut out = XFuseOut::read(*seq1_align_len, *seq2_align_len, seed, *fence_hit);
    let result = x_fused_on(
        reverse,
        &input,
        &mut non_affine_mem.last_seq2_off,
        &mut non_affine_mem.max_score,
        &mut out,
    );
    out.write(seq1_align_len, seq2_align_len, seed, fence_hit);
    Some(result)
}

#[cfg(test)]
#[path = "x_greedy_fused_tests.rs"]
mod tests;
