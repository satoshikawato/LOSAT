//! No NCBI counterpart: one distance of the megablast greedy alignment (`x_greedy_row`) in two
//! passes per block of the band instead of three, for x86-64 CPUs with AVX2/BMI/LZCNT
//! (`LOSAT_X_GREEDYROWFUSE`; `=0` runs `x_greedy_row_avx2`, the three passes).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-553,571-614
//! ```c
//!         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
//!             seq2_index = MAX(last_seq2_off[d - 1][k + 1],
//!                              last_seq2_off[d - 1][k    ]) + 1;
//!             seq2_index = MAX(seq2_index, last_seq2_off[d - 1][k - 1]);
//!             seq1_index = seq2_index + k - diag_origin;
//!             if (seq2_index < 0 || seq1_index + seq2_index < xdrop_score) {
//! ...
//!             index = s_FindFirstMismatch(seq1, seq2, len1, len2,
//! ...
//!             last_seq2_off[d][k] = seq2_index;
//! ```
//!
//! What changes. `x_greedy_row_impl` (the three passes) reads off, in its pass 2, the slide
//! length of every diagonal from one eight-byte comparison, as a vector loop over the block, and
//! stores it in an array that pass 3 reads in the order of the diagonals. Here pass 1 also stores
//! the two load offsets of each diagonal and whether they lie inside both sequences, and the
//! eight-byte comparison runs as scalar code inside the pass 3 loop, for the surviving diagonals
//! only. The vector form spent most of its time moving lanes between vector and general registers
//! (line F s5, probe P6: kernel 17.5 -> 12.3 cycles per cell of the band; the measurements are in
//! the `LOSAT_X_GREEDYROWFUSE` entry of `docs/agents/diagnostics.md`).
//!
//! Why the result is the same. Pass 3 is the pass 3 of `x_greedy_row_impl` statement by
//! statement: the same cells written with the same values in increasing `k`, the same
//! `x_greedy_cell_slow` calls in the same order, the same updates of `diag_lower`,
//! `curr_extent` and the cell that set it, and the same state left at the end. It differs only in
//! where the slide length of a surviving diagonal comes from, and that length is the same value:
//!
//! - Pass 2 of `x_greedy_row_impl` gives -1 ("`x_greedy_cell_slow`") when `unfit < 0`, i.e. when
//!   `seq1_index < 0`, `seq2_index < 0`, `len1 - 8 - seq1_index < 0` or `len2 - 8 - seq2_index <
//!   0`. Here `(seq1_index as u32) > (len1 - 8) as u32` holds exactly when `seq1_index < 0` or
//!   `seq1_index > len1 - 8` (`len1 >= 8` wherever the comparison is used, so `len1 - 8` is not
//!   negative), and the same for `seq2_index`; such a diagonal goes to `x_greedy_cell_slow`
//!   without a load. With `len1 < 8`, `len2 < 8` or the fence already hit (`inline_ok` false)
//!   both give -1 for every diagonal.
//! - Otherwise both read the same eight bytes of each sequence (`seq[at .. at + 8]`, or
//!   `seq[len - 8 - at .. len - at]` for `REVERSE`), and compute the same `diff`, `bits`,
//!   `stop_byte` and `index` with the same integer operations. The `REVERSE` offset `(len - 8) -
//!   at` cannot wrap: `at` is at most `len - 8` inside, and 0 outside. The bytes are read later
//!   than in pass 2, after the `x_greedy_cell_slow` calls of the diagonals before; those write only
//!   `current` and the cold state, never the sequences or the pass 1 arrays.
//! - `x_greedy_row_impl` gives -1 when `diff == 0`, when `stop_byte` is the fence sentry, or
//!   when `index > longest_match_run` (the value when the block started). Here it gives -1 when
//!   `index > min(longest_match_run, 7)` or `stop_byte` is the sentry. `diff == 0` is `bits ==
//!   64`, `index == 8`, which is above the minimum; an `index` of 0-7 is above the minimum exactly
//!   when it is above `longest_match_run`. Both read `longest_match_run` once per block, before
//!   pass 3 (which may raise it through `x_greedy_cell_slow`).
//!
//! So both forms hand pass 3 the same slide length (or -1) for every surviving diagonal, and pass
//! 3 does the same with it. A dropped diagonal (`start < 0`) never reads the slide length in
//! either form. The unit tests below compare this function with `x_greedy_row_avx2` and with the
//! plain build on single distances (every cell of the row, the state, the return value) and on
//! whole traceback alignments (rows, `max_score`, edit script, seed, lengths, return value).

#![allow(clippy::too_many_arguments)]

use super::{
    x_greedy_cell_slow, XGreedyRowCold, XGreedyRowState, FENCE_SENTRY, X_GREEDY_BLOCK, X_NO_CELL,
};
use crate::algorithm::blastn::constants::INVALID_OFFSET;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-553,571-614
/// ```c
///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
///             seq2_index = MAX(last_seq2_off[d - 1][k + 1],
///                              last_seq2_off[d - 1][k    ]) + 1;
/// ...
///             index = s_FindFirstMismatch(seq1, seq2, len1, len2,
/// ...
///             last_seq2_off[d][k] = seq2_index;
/// ```
/// One distance, as `x_greedy_row_impl` (same arguments, same result; module doc). Two passes
/// per block of the band:
///
/// 1. For every diagonal, the offset the reference loop starts its slide from, or -1 when the
///    diagonal fails the X-drop test; then the eight-byte load offsets of both sequences
///    (clamped to 0 when outside, reversed for `REVERSE`) and whether either was outside.
/// 2. The diagonals in the reference order, as pass 3 of `x_greedy_row_impl`; a surviving
///    diagonal whose offsets lie inside both sequences tests its eight bytes here.
///
/// SAFETY: `seq1` and `seq2` must be readable for `len1` and `len2` bytes.
#[inline(always)]
unsafe fn x_greedy_row_fused<const REVERSE: bool>(
    seq1: *const u8,
    seq2: *const u8,
    len1: i32,
    len2: i32,
    previous: &[i32],
    current: &mut [i32],
    tmp_diag_lower: i32,
    diag_origin: i32,
    xdrop_score: i32,
    st: &mut XGreedyRowState,
) -> bool {
    const NOT_A_BASE: u64 = 0xFCFC_FCFC_FCFC_FCFC;
    let band_len = current.len();
    debug_assert!(previous.len() >= band_len + 2);
    let previous = previous.as_ptr();
    let current = current.as_mut_ptr();
    // seq1_index = seq2_index + k - diag_origin = seq2_index + first + cell
    let first = tmp_diag_lower - diag_origin;
    // The eight-byte test needs eight bytes in both sequences. With the fence already hit, no
    // diagonal may skip `x_greedy_cell_slow`: the reference loop returns at its first surviving
    // diagonal.
    let inline_ok = len1 >= 8 && len2 >= 8 && !st.fence_hit;
    let mut cold = XGreedyRowCold {
        diag_lower: st.diag_lower,
        end1_reached: st.end1_reached,
        end2_reached: st.end2_reached,
        curr_extent: st.curr_extent,
        best_cell: X_NO_CELL,
        end1_cell: X_NO_CELL,
        longest_match_run: st.longest_match_run,
        seed_start_q: st.seed_start_q,
        seed_start_s: st.seed_start_s,
        fence_hit: st.fence_hit,
        stopped: false,
    };
    let mut last_valid = X_NO_CELL;

    // All four arrays are written for `0..block_len` before they are read.
    let mut start = [std::mem::MaybeUninit::<i32>::uninit(); X_GREEDY_BLOCK];
    let mut at1 = [std::mem::MaybeUninit::<u32>::uninit(); X_GREEDY_BLOCK];
    let mut at2 = [std::mem::MaybeUninit::<u32>::uninit(); X_GREEDY_BLOCK];
    let mut outside = [std::mem::MaybeUninit::<i32>::uninit(); X_GREEDY_BLOCK];
    let start = start.as_mut_ptr() as *mut i32;
    let at1 = at1.as_mut_ptr() as *mut u32;
    let at2 = at2.as_mut_ptr() as *mut u32;
    let outside = outside.as_mut_ptr() as *mut i32;
    // The last eight-byte load offsets; used only with `inline_ok` (both lengths >= 8).
    let last1 = (len1 - 8) as u32;
    let last2 = (len2 - 8) as u32;
    let mut block = 0usize;
    'row: while block < band_len {
        let block_len = (band_len - block).min(X_GREEDY_BLOCK);
        let block_first = first + block as i32;

        // pass 1
        for i in 0..block_len {
            let p = previous.add(block + i);
            let mut seq2_index = (*p.add(2)).max(*p.add(1)) + 1;
            seq2_index = seq2_index.max(*p);
            let extent = seq2_index + seq2_index + block_first + i as i32;
            *start.add(i) = if seq2_index < 0 || extent < xdrop_score {
                -1
            } else {
                seq2_index
            };
        }
        for i in 0..block_len {
            let seq2_index = *start.add(i);
            let seq1_index = seq2_index + block_first + i as i32;
            // unsigned compares catch both the negative and the too large offsets
            let out1 = (seq1_index as u32) > last1;
            let out2 = (seq2_index as u32) > last2;
            let a1 = if out1 { 0 } else { seq1_index as u32 };
            let a2 = if out2 { 0 } else { seq2_index as u32 };
            *at1.add(i) = if REVERSE { last1 - a1 } else { a1 };
            *at2.add(i) = if REVERSE { last2 - a2 } else { a2 };
            *outside.add(i) = -((out1 | out2) as i32);
        }

        // pass 2
        // 8 when all eight positions match, which is above `run7`
        let run7 = cold.longest_match_run.min(7);
        let mut diag_lower = cold.diag_lower;
        let mut curr_extent = cold.curr_extent;
        let mut best_cell = cold.best_cell;
        for i in 0..block_len {
            let cell = block + i;
            let seq2_index = *start.add(i);
            if seq2_index < 0 {
                if tmp_diag_lower + cell as i32 == diag_lower {
                    diag_lower += 1;
                } else {
                    *current.add(cell) = INVALID_OFFSET;
                }
                continue;
            }
            let index = if inline_ok && *outside.add(i) == 0 {
                let w1 =
                    u64::from_le((seq1.add(*at1.add(i) as usize) as *const u64).read_unaligned());
                let w2 =
                    u64::from_le((seq2.add(*at2.add(i) as usize) as *const u64).read_unaligned());
                let diff = (w1 ^ w2) | (w1 & NOT_A_BASE);
                let bits = if REVERSE {
                    diff.leading_zeros()
                } else {
                    diff.trailing_zeros()
                };
                let shift = bits & 0x38;
                let stop_byte = if REVERSE {
                    (w2 >> (56 - shift)) as u8
                } else {
                    (w2 >> shift) as u8
                };
                let index = (bits >> 3) as i32;
                let mut index_or_slow = index;
                if index > run7 {
                    index_or_slow = -1;
                }
                if stop_byte == FENCE_SENTRY {
                    index_or_slow = -1;
                }
                index_or_slow
            } else {
                -1
            };
            if index < 0 {
                cold.diag_lower = diag_lower;
                cold.curr_extent = curr_extent;
                cold.best_cell = best_cell;
                x_greedy_cell_slow::<REVERSE>(
                    seq1,
                    seq2,
                    len1,
                    len2,
                    seq2_index + block_first + i as i32,
                    seq2_index,
                    cell,
                    tmp_diag_lower + cell as i32,
                    current,
                    &mut cold,
                );
                if cold.stopped {
                    last_valid = cell;
                    break 'row;
                }
                diag_lower = cold.diag_lower;
                curr_extent = cold.curr_extent;
                best_cell = cold.best_cell;
                continue;
            }
            let stop2 = seq2_index + index;
            *current.add(cell) = stop2;
            // seq1_index + index + seq2_index + index
            let extent = stop2 + stop2 + block_first + i as i32;
            if extent > curr_extent {
                curr_extent = extent;
                best_cell = cell;
            }
        }
        cold.diag_lower = diag_lower;
        cold.curr_extent = curr_extent;
        cold.best_cell = best_cell;

        // last surviving diagonal of this block
        let mut j = block_len;
        while j > 0 {
            j -= 1;
            if *start.add(j) >= 0 {
                last_valid = block + j;
                break;
            }
        }
        block += block_len;
    }

    st.diag_lower = cold.diag_lower;
    if last_valid != X_NO_CELL {
        let k = tmp_diag_lower + last_valid as i32;
        st.diag_upper = if cold.end1_cell == last_valid {
            k - 1
        } else {
            k
        };
    }
    st.end1_reached = cold.end1_reached;
    st.end2_reached = cold.end2_reached;
    st.curr_extent = cold.curr_extent;
    if cold.best_cell != X_NO_CELL {
        st.curr_seq2_index = *current.add(cold.best_cell);
        st.curr_diag = tmp_diag_lower + cold.best_cell as i32;
    }
    st.longest_match_run = cold.longest_match_run;
    st.seed_start_q = cold.seed_start_q;
    st.seed_start_s = cold.seed_start_s;
    st.fence_hit = cold.fence_hit;
    cold.stopped
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537
/// ```c
///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
/// ```
/// `x_greedy_row_fused` for the C `reverse` argument, compiled for AVX2/BMI (`tzcnt` / `lzcnt`).
///
/// SAFETY: as `x_greedy_row_avx2` (the CPU must support the features).
#[target_feature(enable = "avx2,bmi1,bmi2,lzcnt")]
pub(super) unsafe fn x_greedy_row_fused_avx2(
    reverse: bool,
    seq1: *const u8,
    seq2: *const u8,
    len1: i32,
    len2: i32,
    previous: &[i32],
    current: &mut [i32],
    tmp_diag_lower: i32,
    diag_origin: i32,
    xdrop_score: i32,
    st: &mut XGreedyRowState,
) -> bool {
    // No NCBI counterpart (memory placement only): see `x_frame_align`.
    crate::utils::x_frame::x_frame_align();
    if reverse {
        x_greedy_row_fused::<true>(
            seq1,
            seq2,
            len1,
            len2,
            previous,
            current,
            tmp_diag_lower,
            diag_origin,
            xdrop_score,
            st,
        )
    } else {
        x_greedy_row_fused::<false>(
            seq1,
            seq2,
            len1,
            len2,
            previous,
            current,
            tmp_diag_lower,
            diag_origin,
            xdrop_score,
            st,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::super::{
        blast_greedy_align, format_gap_prelim_edit_block_for_trace, x_greedy_row_any,
        x_greedy_row_avx2, GapPrelimEditBlock, GreedyNonAffineMem, GreedySeed, XGreedyRowKernel,
        X_GREEDY_ROW_TEST_KERNEL,
    };
    use super::*;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
        fn pick<T: Copy>(&mut self, values: &[T]) -> T {
            values[self.below(values.len() as u64) as usize]
        }
        fn range(&mut self, low: i32, high: i32) -> i32 {
            low + self.below((high - low + 1) as u64) as i32
        }
    }

    fn avx2() -> bool {
        std::is_x86_feature_detected!("avx2")
            && std::is_x86_feature_detected!("bmi1")
            && std::is_x86_feature_detected!("bmi2")
            && std::is_x86_feature_detected!("lzcnt")
    }

    /// A pair of related sequences of `len` template bases (identity, substitutions, indels,
    /// ambiguity bytes in either, sometimes fence bytes in `seq2`), reversed for `reverse` (the
    /// left extension reads `seq[len - 1 - i]`), with matching bases and fence bytes after the
    /// lengths that a read past the end would turn into a different result.
    fn related(rng: &mut Rng, len: usize, reverse: bool) -> (Vec<u8>, i32, Vec<u8>, i32) {
        let identity = rng.pick(&[1000u64, 1000, 995, 980, 950, 900, 800, 700, 600]);
        let indel_share = rng.pick(&[0u64, 0, 100, 300, 600]); // of the differences, per mille
        let ambiguity1 = rng.pick(&[0u64, 0, 0, 5, 30]); // per mille
        let ambiguity2 = rng.pick(&[0u64, 0, 0, 5, 30]);
        let (mut seq1, mut seq2) = (Vec::new(), Vec::new());
        for _ in 0..len {
            let base = rng.below(4) as u8;
            if rng.below(1000) < identity {
                seq1.push(base);
                seq2.push(base);
                continue;
            }
            let roll = rng.below(1000);
            if roll < indel_share / 2 {
                seq1.push(base);
            } else if roll < indel_share {
                seq2.push(base);
            } else {
                seq1.push(base);
                seq2.push((base + 1 + rng.below(3) as u8) & 3);
            }
        }
        for byte in seq1.iter_mut() {
            if rng.below(1000) < ambiguity1 {
                *byte = 4 + rng.below(12) as u8;
            }
        }
        for byte in seq2.iter_mut() {
            if rng.below(1000) < ambiguity2 {
                *byte = 4 + rng.below(12) as u8;
            }
        }
        if !seq2.is_empty() {
            for _ in 0..rng.pick(&[0u64, 0, 0, 1, 3]) {
                let at = rng.below(seq2.len() as u64) as usize;
                seq2[at] = FENCE_SENTRY;
            }
        }
        let (len1, len2) = (seq1.len() as i32, seq2.len() as i32);
        if reverse {
            seq1.reverse();
            seq2.reverse();
        }
        for _ in 0..rng.below(16) {
            seq1.push(rng.pick(&[0u8, 1, 2, 3, FENCE_SENTRY]));
        }
        for _ in 0..rng.below(16) {
            seq2.push(rng.pick(&[0u8, 1, 2, 3, FENCE_SENTRY]));
        }
        (seq1, len1, seq2, len2)
    }

    /// One distance: the sequences, the row of distance `d - 1` and the arguments.
    struct Row {
        reverse: bool,
        seq1: Vec<u8>,
        len1: i32,
        seq2: Vec<u8>,
        len2: i32,
        previous: Vec<i32>,
        /// stale cells the call may leave as they are
        current: Vec<i32>,
        tmp_diag_lower: i32,
        diag_origin: i32,
        xdrop_score: i32,
        st: XGreedyRowState,
    }

    /// Bands of 1-600 diagonals (one to three blocks of `X_GREEDY_BLOCK`), lengths below eight
    /// and up to about 300, offsets of distance `d - 1` from rows that rise along the band and
    /// from noise (`INVALID_OFFSET`, -1, past the ends), the X-drop bound from always passing to
    /// dropping most diagonals, the fence hit on entry, a seed run from 0 to above eight.
    fn row(rng: &mut Rng) -> Row {
        let reverse = rng.below(2) == 1;
        let len = match rng.below(6) {
            0 => rng.below(12),
            1 | 2 => 12 + rng.below(60),
            _ => 72 + rng.below(240),
        } as usize;
        let (seq1, len1, seq2, len2) = related(rng, len, reverse);
        let band_len = match rng.below(6) {
            0 => 1 + rng.below(8),
            1 | 2 | 3 => 9 + rng.below(250),
            _ => 250 + rng.below(350),
        } as usize;
        // seq1_index = seq2_index + first + cell
        let first = rng.range(-(band_len as i32) - 10, 20);
        let tmp_diag_lower = rng.range(-30, 30);
        let diag_origin = tmp_diag_lower - first;
        let noise = rng.pick(&[0u64, 0, 20, 100, 400]); // per mille of the cells
        let base = rng.range(-3, len2.max(1));
        let mut previous: Vec<i32> = (0..band_len + 2)
            .map(|j| {
                if rng.below(1000) < noise {
                    let any = rng.range(-2, len2 + 3);
                    return rng.pick(&[INVALID_OFFSET, -1, len2, len2 + 1, any]);
                }
                // near the main diagonal of the band the rows reach furthest
                let o = first + j as i32 - 1;
                (base - o.abs() / rng.pick(&[1, 2, 4]) + rng.range(-2, 2)).max(INVALID_OFFSET)
            })
            .collect();
        // The reference loop never takes a diagonal with seq1_index < 0 (its slide would index
        // before the sequence): raise the offsets that would make one survive with it until none
        // does (raising an offset only raises the starts of its three diagonals).
        loop {
            let mut raised = false;
            for i in 0..band_len {
                let s2 = previous[i + 2]
                    .max(previous[i + 1])
                    .saturating_add(1)
                    .max(previous[i]);
                let o = first + i as i32;
                if s2 >= 0 && s2 + o < 0 {
                    previous[i + 1] = previous[i + 1].max(-o - 1);
                    raised = true;
                }
            }
            if !raised {
                break;
            }
        }
        let xdrop_score = match rng.below(4) {
            0 => i32::MIN / 4,
            1 => rng.range(-20, 20),
            _ => rng.range(0, 2 * (len1 + len2) + 10),
        };
        let diag_lower = match rng.below(4) {
            0 | 1 => tmp_diag_lower,
            2 => tmp_diag_lower + rng.range(0, 3),
            _ => rng.range(tmp_diag_lower - 5, tmp_diag_lower + band_len as i32 + 5),
        };
        let st = XGreedyRowState {
            diag_lower,
            diag_upper: tmp_diag_lower + band_len as i32 - 1,
            end1_reached: rng.below(4) == 0,
            end2_reached: rng.below(4) == 0,
            curr_extent: rng.range(-5, len1 + len2 + 5),
            curr_seq2_index: rng.range(-3, len2),
            curr_diag: rng.range(-3, 3),
            longest_match_run: rng.pick(&[0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 20, 1000]),
            seed_start_q: rng.range(-1, len1),
            seed_start_s: rng.range(-1, len2),
            fence_hit: rng.below(20) == 0,
        };
        let current = (0..band_len).map(|_| rng.range(-9, 999)).collect();
        Row {
            reverse,
            seq1,
            len1,
            seq2,
            len2,
            previous,
            current,
            tmp_diag_lower,
            diag_origin,
            xdrop_score,
            st,
        }
    }

    type RowFn = unsafe fn(
        bool,
        *const u8,
        *const u8,
        i32,
        i32,
        &[i32],
        &mut [i32],
        i32,
        i32,
        i32,
        &mut XGreedyRowState,
    ) -> bool;

    /// The row `run` leaves: (return value, every cell, the state).
    fn run(f: RowFn, r: &Row) -> (bool, Vec<i32>, XGreedyRowState) {
        let mut current = r.current.clone();
        let mut st = r.st;
        // SAFETY: the lengths lie within the vectors (bytes follow them) and the CPU has AVX2.
        let stopped = unsafe {
            f(
                r.reverse,
                r.seq1.as_ptr(),
                r.seq2.as_ptr(),
                r.len1,
                r.len2,
                &r.previous,
                &mut current,
                r.tmp_diag_lower,
                r.diag_origin,
                r.xdrop_score,
                &mut st,
            )
        };
        (stopped, current, st)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-553,571-614
    /// ```c
    ///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
    /// ...
    ///             index = s_FindFirstMismatch(seq1, seq2, len1, len2,
    /// ...
    ///             last_seq2_off[d][k] = seq2_index;
    /// ```
    /// Single distances on random rows: the fused build leaves the return value, every cell of
    /// the row (stale cells included) and every field of the state that `x_greedy_row_avx2`
    /// (`LOSAT_X_GREEDYROWFUSE=0`) and the plain build leave. The generator must reach the
    /// eight-byte test's outcomes (a slide length, a long slide, a stop at the fence byte, offsets
    /// outside a sequence), rows of several blocks, rows with fewer than eight bytes or the fence
    /// hit on entry, and fence stops.
    #[test]
    fn x_greedy_row_fused_matches_three_passes_on_random_rows() {
        if !avx2() {
            eprintln!("x_greedy_row_fused_matches_three_passes_on_random_rows: no AVX2, skipped");
            return;
        }
        let rows: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(20_000);
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        // slide length read off, longer than the run, fence byte, outside a sequence
        let mut cells = [0usize; 4];
        let (mut blocks, mut short, mut gated, mut stopped) = (0, 0, 0, 0);
        for n in 0..rows {
            let r = row(&mut rng);
            let what = format!(
                "row {n}: reverse={} len1={} len2={} band={} first={} xdrop_score={} st={:?}",
                r.reverse,
                r.len1,
                r.len2,
                r.current.len(),
                r.tmp_diag_lower - r.diag_origin,
                r.xdrop_score,
                r.st
            );
            let expected = run(x_greedy_row_avx2, &r);
            assert_eq!(run(x_greedy_row_fused_avx2, &r), expected, "{what}: fused");
            assert_eq!(run(x_greedy_row_any, &r), expected, "{what}: plain");

            blocks += usize::from(r.current.len() > X_GREEDY_BLOCK);
            short += usize::from(r.len1 < 8 || r.len2 < 8);
            gated += usize::from(r.st.fence_hit);
            stopped += usize::from(expected.0 && !r.st.fence_hit);
            if r.len1 < 8 || r.len2 < 8 || r.st.fence_hit {
                continue;
            }
            // what the eight-byte test sees on the surviving diagonals (with the run on entry)
            let first = r.tmp_diag_lower - r.diag_origin;
            for i in 0..r.current.len() {
                let p = &r.previous[i..i + 3];
                let s2 = p[2].max(p[1]).saturating_add(1).max(p[0]);
                let s1 = s2 + first + i as i32;
                if s2 < 0 || s1 + s2 < r.xdrop_score {
                    continue;
                }
                if s1 > r.len1 - 8 || s2 > r.len2 - 8 {
                    cells[3] += 1;
                    continue;
                }
                let (at1, at2) = if r.reverse {
                    ((r.len1 - 8 - s1) as usize, (r.len2 - 8 - s2) as usize)
                } else {
                    (s1 as usize, s2 as usize)
                };
                let w1 = u64::from_le_bytes(r.seq1[at1..at1 + 8].try_into().unwrap());
                let w2 = u64::from_le_bytes(r.seq2[at2..at2 + 8].try_into().unwrap());
                let diff = (w1 ^ w2) | (w1 & 0xFCFC_FCFC_FCFC_FCFC);
                let index = if r.reverse {
                    diff.leading_zeros()
                } else {
                    diff.trailing_zeros()
                } as i32
                    >> 3;
                let stop = if r.reverse { 7 - index } else { index };
                if index < 8 && r.seq2[at2 + stop as usize] == FENCE_SENTRY {
                    cells[2] += 1;
                } else if index > r.st.longest_match_run.min(7) {
                    cells[1] += 1;
                } else {
                    cells[0] += 1;
                }
            }
        }
        let summary = format!(
            "rows {rows}: several blocks {blocks}, fewer than eight bytes {short}, fence on \
             entry {gated}, fence stops {stopped}; diagonals: slide read off {}, longer than the \
             run {}, fence byte {}, outside {}",
            cells[0], cells[1], cells[2], cells[3]
        );
        eprintln!("x_greedy_row_fused_matches_three_passes_on_random_rows: {summary}");
        assert!(cells.iter().all(|&c| c * 10 > rows), "{summary}");
        assert!(blocks * 10 > rows && short * 20 > rows, "{summary}");
        assert!(gated * 40 > rows && stopped * 40 > rows, "{summary}");
    }

    /// What one `blast_greedy_align` call returns and leaves.
    #[derive(Debug, PartialEq, Eq)]
    struct Outcome {
        ret: i32,
        seq1_align_len: i32,
        seq2_align_len: i32,
        seed: (i32, i32, i32),
        fence_hit: bool,
        edit: String,
    }

    /// `blast_greedy_align` with an edit block (traceback: `x_greedy_row` runs every distance), or
    /// score only with the fence hit on entry (the fused loop declines it: `x_greedy_row` runs
    /// distance 1 and stops at its first surviving diagonal), with the build `kernel`.
    fn align(
        kernel: XGreedyRowKernel,
        mem: &mut GreedyNonAffineMem,
        (seq1, len1, seq2, len2): (&[u8], i32, &[u8], i32),
        (reverse, xdrop, match_cost, mismatch_cost, max_dist): (bool, i32, i32, i32, i32),
        traceback: bool,
    ) -> Outcome {
        X_GREEDY_ROW_TEST_KERNEL.with(|k| k.set(Some(kernel)));
        let mut edit = GapPrelimEditBlock::new();
        // values the call overwrites
        let mut seed = GreedySeed {
            start_q: -7,
            start_s: -7,
            match_length: -7,
        };
        let mut fence_hit = !traceback;
        let (mut seq1_align_len, mut seq2_align_len) = (-5, -5);
        let ret = blast_greedy_align(
            seq1,
            len1,
            seq2,
            len2,
            reverse,
            xdrop,
            match_cost,
            mismatch_cost,
            &mut seq1_align_len,
            &mut seq2_align_len,
            max_dist,
            mem,
            traceback.then_some(&mut edit),
            &mut fence_hit,
            &mut seed,
        );
        X_GREEDY_ROW_TEST_KERNEL.with(|k| k.set(None));
        Outcome {
            ret,
            seq1_align_len,
            seq2_align_len,
            seed: (seed.start_q, seed.start_s, seed.match_length),
            fence_hit,
            edit: format_gap_prelim_edit_block_for_trace(&edit),
        }
    }

    /// The scratch two calls leave, whole: both rows, `max_score` and the traceback rows.
    fn assert_same_memory(a: &GreedyNonAffineMem, b: &GreedyNonAffineMem, what: &str) {
        assert!(a.last_seq2_off == b.last_seq2_off, "{what}: rows differ");
        assert!(a.max_score == b.max_score, "{what}: max_score differs");
        assert!(
            a.traceback_used == b.traceback_used && a.traceback_pool == b.traceback_pool,
            "{what}: traceback rows differ"
        );
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:509-537,682-751
    /// ```c
    ///     for (d = 1; d <= max_dist; d++) {
    /// ...
    ///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
    /// ...
    ///     if (edit_block == NULL)
    ///         return best_dist;
    /// ...
    ///     s_GetNextNonAffineTback(last_seq2_off, d, diag, &new_seq2_index);
    /// ```
    /// Sequences of whole alignments on one scratch per build, so that the cells earlier calls
    /// left matter: traceback calls and score-only calls with the fence hit on entry, on related
    /// sequences of up to about 3,000 bases with megablast's and other costs. The fused build
    /// returns the values, the edit script and the seed of `x_greedy_row_avx2`
    /// (`LOSAT_X_GREEDYROWFUSE=0`) and of the plain build, and leaves the same scratch.
    #[test]
    fn x_greedy_row_fused_matches_three_passes_on_alignments() {
        if !avx2() {
            eprintln!("x_greedy_row_fused_matches_three_passes_on_alignments: no AVX2, skipped");
            return;
        }
        let histories: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .map_or(300, |n: usize| (n / 10).max(1));
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let (mut calls, mut converged, mut stopped, mut fences, mut long) = (0, 0, 0, 0, 0);
        for history in 0..histories {
            let mut three = GreedyNonAffineMem::new();
            let mut fused = GreedyNonAffineMem::new();
            let mut plain = GreedyNonAffineMem::new();
            for step in 0..1 + rng.below(8) {
                let len = match rng.below(4) {
                    0 => rng.below(40),
                    1 => 40 + rng.below(300),
                    _ => 340 + rng.below(2700),
                } as usize;
                let reverse = rng.below(2) == 1;
                let (seq1, len1, seq2, len2) = related(&mut rng, len, reverse);
                // megablast's costs (reward 1 doubled), and costs with match_cost / 2 != 1
                let (match_cost, mismatch_cost) =
                    rng.pick(&[(2, 4), (2, 4), (2, 4), (2, 6), (2, 2), (4, 6), (6, 8)]);
                let xdrop = rng.pick(&[0, 1, 5, 13, 20, 26, 27, 40, 54, 100, 300]);
                let max_dist = match rng.below(3) {
                    0 => 1 + rng.below(12) as i32,
                    _ => 1000.min(len2 / 2 + 1),
                };
                let traceback = rng.below(6) != 0;
                let input = (&seq1[..], len1, &seq2[..], len2);
                let costs = (reverse, xdrop, match_cost, mismatch_cost, max_dist);
                let what = format!(
                    "history {history} step {step}: len1={len1} len2={len2} reverse={reverse} \
                     xdrop={xdrop} costs={match_cost}/{mismatch_cost} max_dist={max_dist} \
                     traceback={traceback}"
                );
                let expected = align(XGreedyRowKernel::Avx2, &mut three, input, costs, traceback);
                let got = align(XGreedyRowKernel::Fused, &mut fused, input, costs, traceback);
                assert_eq!(got, expected, "{what}: fused");
                assert_same_memory(&fused, &three, &format!("{what}: fused"));
                let got = align(XGreedyRowKernel::Plain, &mut plain, input, costs, traceback);
                assert_eq!(got, expected, "{what}: plain");
                assert_same_memory(&plain, &three, &format!("{what}: plain"));

                calls += 1;
                long += usize::from(len >= 340);
                if expected.fence_hit && traceback {
                    fences += 1;
                } else if expected.ret < 0 {
                    stopped += 1;
                } else {
                    converged += 1;
                }
            }
        }
        let summary = format!(
            "calls {calls}: converged {converged}, past max_dist {stopped}, fence {fences}; \
             340+ template bases {long}"
        );
        eprintln!("x_greedy_row_fused_matches_three_passes_on_alignments: {summary}");
        assert!(converged * 3 > calls && long * 3 > calls, "{summary}");
        assert!(stopped * 40 > calls && fences * 40 > calls, "{summary}");
    }
}
