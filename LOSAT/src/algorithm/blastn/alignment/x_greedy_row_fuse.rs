//! No NCBI counterpart: one distance of the megablast greedy alignment (`x_greedy_row`) in two
//! passes per block of the band instead of three, for x86-64 CPUs with AVX2/BMI/LZCNT
//! (`LOSAT_X_GREEDYROWFUSE`; `=0` runs `x_greedy_row_avx2`, the three passes).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-553,571-614
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
//! `x_greedy_cell_slow` calls and the same `x_past_seq1_end` diagonals in the same order, the
//! same updates of `diag_lower`, `curr_extent` and the cell that set it, and the same state left
//! at the end. It differs only in
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
    x_greedy_cell_slow, x_past_seq1_end, XGreedyRowCold, XGreedyRowState, FENCE_SENTRY,
    X_GREEDY_BLOCK, X_NO_CELL,
};
use crate::algorithm::blastn::constants::INVALID_OFFSET;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-553,571-614
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
                // A diagonal past the end of `seq1` (`x_past_seq1_end`), unless the fence was
                // already hit on entry (C returns at the first kept diagonal then).
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:565-600
                let seq1_index = seq2_index + block_first + i as i32;
                if !cold.fence_hit
                    && x_past_seq1_end::<REVERSE>(seq2, len1, len2, seq1_index, seq2_index)
                {
                    *current.add(cell) = seq2_index;
                    let extent = seq1_index + seq2_index;
                    if extent > curr_extent {
                        curr_extent = extent;
                        best_cell = cell;
                    }
                    continue;
                }
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
#[path = "x_greedy_row_fuse_tests.rs"]
mod tests;
