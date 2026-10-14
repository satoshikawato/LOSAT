//! EXPERIMENT (LOSAT_X_GREEDYFUSE, shadow LOSAT_X_GREEDYFUSESHADOW): the loop over distances of
//! the score-only `BLAST_GreedyAlign` on an unpacked subject in one function, each distance in one
//! pass over its diagonals.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:509-537,680-686
//! ```c
//!     for (d = 1; d <= max_dist; d++) {
//! ...
//!         last_seq2_off[d - 1][diag_lower-1] = kInvalidOffset;
//! ...
//!         xdrop_score = (Int4)ceil((double)xdrop_score / (match_cost / 2));
//! ...
//!         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
//! ...
//!     }   /* end loop over distinct distances */
//!
//!     if (!converged)
//!         return -1;
//!
//!     if (edit_block == NULL)
//!         return best_dist;
//! ```
//!
//! Where it runs. `blast_greedy_align` (the port of `BLAST_GreedyAlign`) calls
//! `x_greedy_fused_dispatch` after its set-up: `max_score[0..xdrop_offset]` cleared, the first
//! `s_FindFirstMismatch`, the early return, `last_seq2_off[0][diag_origin]`, `max_score[0]` and the
//! first diagonal bounds (greedy_align.c:463-505). The call takes over the rest of the function
//! when all of these hold: the switch is on, there is no edit block (score only, the preliminary
//! stage), `rem == 4` (an unpacked subject: `LOSAT_X_GREEDYPRELIM`, on by default), the fence was
//! not hit by the first scan, `1 <= max_dist`, and `len1` / `len2` lie within the slices.
//! Otherwise the ported loop runs, unchanged. Traceback, packed subjects and the affine
//! alignment keep their paths.
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
//! - The same slide. `s_FindFirstMismatch` with `rem == 4` continues while
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
//!   fence test with its bound and the C updates in order.
//! - The same integers. `seq2_index`, `seq1_index`, the extent and the scores are the C
//!   expressions on `i32`. `xdrop_score` is the reference's
//!   `((raw as f64) / (match_cost as f64 / 2.0)).ceil() as i32`, except for `match_cost == 2`
//!   (megablast's reward 1, doubled by `BLAST_AffineGreedyAlign`): the divisor is then 1.0, and
//!   an `i32` divided by 1.0 and rounded up is itself (exact in `f64`), so `raw` is used as it is.
//! - The same exits, with the same values in the outputs: converged -> `best_dist`; past
//!   `max_dist` -> -1; a fence at the stop byte of a slide -> `fence_hit` and 0, at once (the cell
//!   of that diagonal and the seed are not written, as in C). The outputs `seq1_align_len`,
//!   `seq2_align_len` and the seed keep the values of their last update in every case.
//!
//! The variables are read once. `LOSAT_X_GREEDYFUSE` (off unless set to any value; not in the
//! strict set, `utils/xswitch.rs`) runs this module. `LOSAT_X_GREEDYFUSESHADOW` (implies it) runs
//! this module on copies of the rows, of `max_score` and of the outputs, then the reference
//! (`blast_greedy_align` again, from its start, with this module off for that call), and panics
//! unless the return value, both align lengths, the seed, `fence_hit`, both rows and `max_score`
//! are equal; `main` prints one `[X_SHADOW] GREEDYFUSE compared=N` line.

#![allow(clippy::too_many_arguments)]

use super::{blast_greedy_align, x_first_mismatch, GreedyNonAffineMem, GreedySeed, FENCE_SENTRY};
use crate::algorithm::blastn::constants::INVALID_OFFSET;
use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

#[cfg(test)]
thread_local! {
    /// Lets a test choose the mode per call (0 off, 1 this module, 2 shadow).
    static X_GREEDY_FUSE_TEST_MODE: Cell<Option<u8>> = const { Cell::new(None) };
    /// Lets a test run the plain build on a CPU with AVX2.
    static X_GREEDY_FUSE_TEST_PLAIN: Cell<bool> = const { Cell::new(false) };
    /// Diagonals that took `x_cell_cold` with fewer than eight positions left, and with a slide
    /// of eight or more.
    static X_GREEDY_FUSE_TEST_REST: Cell<(u64, u64)> = const { Cell::new((0, 0)) };
}

thread_local! {
    // No NCBI counterpart: set while the shadow runs the reference, so that the reference call
    // does not take this module again.
    static X_SHADOW_REFERENCE: Cell<bool> = const { Cell::new(false) };
}

// No NCBI counterpart: the number of calls the shadow compared (relaxed; it orders nothing).
static X_SHADOW_COMPARED: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: reads LOSAT_X_GREEDYFUSE / LOSAT_X_GREEDYFUSESHADOW once (off unless set
/// to any value; not in the strict set). 0 = the ported loop, 1 = this module, 2 = both and
/// compare; SHADOW wins. It chooses a LOSAT code path and changes no value NCBI computes.
fn x_greedy_fuse_mode() -> u8 {
    #[cfg(test)]
    if let Some(mode) = X_GREEDY_FUSE_TEST_MODE.with(Cell::get) {
        return mode;
    }
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_GREEDYFUSESHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_GREEDYFUSE").is_some() {
            1
        } else {
            0
        }
    })
}

/// The summary line `main` prints at exit with LOSAT_X_GREEDYFUSESHADOW.
// No NCBI counterpart: prints a counter; it does not change any value NCBI computes.
pub fn print_summary() {
    if x_greedy_fuse_mode() != 2 {
        return;
    }
    let compared = X_SHADOW_COMPARED.load(Ordering::Relaxed);
    eprintln!(
        "[X_SHADOW] GREEDYFUSE compared={compared} ({})",
        if compared == 0 {
            "nothing compared: no score-only greedy alignment on an unpacked subject ran"
        } else {
            "every score-only greedy alignment identical"
        }
    );
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:380-387
/// ```c
/// Int4 BLAST_GreedyAlign(const Uint1* seq1, Int4 len1,
///                        const Uint1* seq2, Int4 len2,
///                        Boolean reverse, Int4 xdrop_threshold,
///                        Int4 match_cost, Int4 mismatch_cost,
///                        Int4* seq1_align_len, Int4* seq2_align_len,
///                        SGreedyAlignMem* aux_data,
///                        GapPrelimEditBlock *edit_block, Uint1 rem,
///                        Boolean * fence_hit, SGreedySeed *seed)
/// ```
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
/// ```c
/// typedef struct SGreedyAlignMem {
///    Int4 max_dist;
///    Int4 xdrop;
///    Int4** last_seq2_off;
///    Int4* max_score;
/// ...
/// } SGreedyAlignMem;
/// ```
/// The two rolling rows of `last_seq2_off` (origin 0) and `aux_data->max_score`.
#[derive(Clone, Copy)]
struct XFuseMem {
    rows: [*mut i32; 2],
    max_score: *mut i32,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:463-470
/// ```c
///     *seq1_align_len = index;
///     *seq2_align_len = index;
///     seq1_index = index;
///
///     seed->start_q = 0;
///     seed->start_s = 0;
///     seed->match_length = longest_match_run = index;
/// ```
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

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:396-404,537-541
/// ```c
///     Int4 diag_lower, diag_upper;
/// ...
///     Int4 longest_match_run;
/// ...
///         Int4 curr_extent = 0;
///         Int4 curr_seq2_index = 0;
///         Int4 curr_diag = 0;
/// ```
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
/// ```c
///             if (index > longest_match_run) {
///                 seed->start_q = seq1_index;
///                 seed->start_s = seq2_index;
/// ...
///             if (seq2_index == len2) {
///                 diag_lower = k + 1;
///                 end2_reached = TRUE;
/// ```
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
/// ```c
///             diag_upper = k;
/// ...
///             index = s_FindFirstMismatch(seq1, seq2, len1, len2,
///                                         seq1_index, seq2_index,
///                                         fence_hit, reverse, rem);
///             if(fence_hit && *fence_hit){
///             	return 0;
///             }
///
///             if (index > longest_match_run) {
///                 seed->start_q = seq1_index;
///                 seed->start_s = seq2_index;
///                 seed->match_length = longest_match_run = index;
///             }
///             seq1_index += index;
///             seq2_index += index;
/// ...
///             last_seq2_off[d][k] = seq2_index;
/// ...
///             if (seq1_index + seq2_index > curr_extent) {
///                 curr_extent = seq1_index + seq2_index;
///                 curr_seq2_index = seq2_index;
///                 curr_diag = k;
///             }
/// ...
///             if (seq2_index == len2) {
///                 diag_lower = k + 1;
///                 end2_reached = TRUE;
///             }
///             if (seq1_index == len1) {
///                 diag_upper = k - 1;
///                 end1_reached = TRUE;
///             }
/// ```
/// A kept diagonal, as the C loop handles it after its X-drop test: the slide by
/// `s_FindFirstMismatch` for `rem == 4` (`x_first_mismatch`: the same predicate, eight bytes at a
/// time, then byte by byte), its fence test with its bound, and every update in the C order. The
/// diagonals whose slide one eight-byte comparison does not settle (fewer than eight positions
/// left, eight or more matches, the fence at the stop, a longer run than the seed's). With a
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
/// ```c
///     for (d = 1; d <= max_dist; d++) {
///         Int4 xdrop_score;
///         Int4 curr_score;
///         Int4 curr_extent = 0;
///         Int4 curr_seq2_index = 0;
///         Int4 curr_diag = 0;
///         Int4 tmp_diag_lower = diag_lower;
///         Int4 tmp_diag_upper = diag_upper;
/// ...
///     }   /* end loop over distinct distances */
///
///     if (!converged)
///         return -1;
///
///     if (edit_block == NULL)
///         return best_dist;
/// ```
/// The distance loop of `BLAST_GreedyAlign` without an edit block, for `rem == 4`, from the
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
    // ```c
    //     max_score = aux_data->max_score + xdrop_offset;
    // ```
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
        // ```c
        //         last_seq2_off[d - 1][diag_lower-1] = kInvalidOffset;
        //         last_seq2_off[d - 1][diag_lower] = kInvalidOffset;
        //         last_seq2_off[d - 1][diag_upper] = kInvalidOffset;
        //         last_seq2_off[d - 1][diag_upper+1] = kInvalidOffset;
        // ```
        *previous.add((hot.diag_lower - 1) as usize) = INVALID_OFFSET;
        *previous.add(hot.diag_lower as usize) = INVALID_OFFSET;
        *previous.add(hot.diag_upper as usize) = INVALID_OFFSET;
        *previous.add((hot.diag_upper + 1) as usize) = INVALID_OFFSET;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:531-533
        // ```c
        //         xdrop_score = max_score[d - xdrop_offset] +
        //                       (match_cost + mismatch_cost) * d - xdrop_threshold;
        //         xdrop_score = (Int4)ceil((double)xdrop_score / (match_cost / 2));
        // ```
        let raw = *max_score_raw.add(d as usize) + op_cost * d - xdrop_threshold;
        let xdrop_score = if unit_divisor {
            raw
        } else {
            ((raw as f64) / divisor).ceil() as i32
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:513-518
        // ```c
        //         Int4 curr_extent = 0;
        //         Int4 curr_seq2_index = 0;
        //         Int4 curr_diag = 0;
        //         Int4 tmp_diag_lower = diag_lower;
        //         Int4 tmp_diag_upper = diag_upper;
        // ```
        let tmp_diag_lower = hot.diag_lower;
        let tmp_diag_upper = hot.diag_upper;
        let mut diag_lower = hot.diag_lower;
        let mut diag_upper = hot.diag_upper;
        let mut curr_extent: i32 = 0;
        let mut curr_diag: i32 = 0;
        let mut longest_match_run = hot.longest_match_run;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-553
        // ```c
        //         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
        // ...
        //             seq2_index = MAX(last_seq2_off[d - 1][k + 1],
        //                              last_seq2_off[d - 1][k    ]) + 1;
        //             seq2_index = MAX(seq2_index, last_seq2_off[d - 1][k - 1]);
        //             seq1_index = seq2_index + k - diag_origin;
        // ```
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
            // ```c
            //             if (seq2_index < 0 || seq1_index + seq2_index < xdrop_score) {
            // ...
            //                 if (k == diag_lower)
            //                     diag_lower++;
            //                 else
            //                     last_seq2_off[d][k] = kInvalidOffset;
            //                 continue;
            //             }
            //             diag_upper = k;
            // ```
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
            // ```c
            //             diag_upper = k;
            // ...
            //             index = s_FindFirstMismatch(seq1, seq2, len1, len2,
            // ...
            //             last_seq2_off[d][k] = seq2_index;
            // ...
            //             if (seq1_index + seq2_index > curr_extent) {
            // ```
            // One eight-byte comparison (module doc) settles the slide when both sequences have
            // eight positions left and it stops within them. The C code then runs, for this
            // diagonal: `diag_upper = k`; no fence (the stop byte of `seq2` is not the sentry);
            // no seed (`index <= longest_match_run`); the cell; the extent; no end test fires
            // (the stop lies inside both sequences). Every other kept diagonal is `x_cell_cold`.
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
        // ```c
        //             if (seq1_index + seq2_index > curr_extent) {
        //                 curr_extent = seq1_index + seq2_index;
        //                 curr_seq2_index = seq2_index;
        //                 curr_diag = k;
        //             }
        // ...
        //         curr_score = curr_extent * (match_cost / 2) -
        //                         d * (match_cost + mismatch_cost);
        // ...
        //         if (curr_score >= max_score[d - 1]) {
        //             max_score[d] = curr_score;
        //             best_dist = d;
        //             best_diag = curr_diag;
        //             *seq2_align_len = curr_seq2_index;
        //             *seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        //         }
        //         else {
        //             max_score[d] = max_score[d - 1];
        //         }
        // ```
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
        // ```c
        //         if (diag_lower > diag_upper) {
        //             converged = TRUE;
        //             break;
        //         }
        // ...
        //         if (!end2_reached)
        //             diag_lower--;
        //         if (!end1_reached)
        //             diag_upper++;
        // ...
        //         if (edit_block == NULL) {
        // ...
        //             last_seq2_off[d + 1] = last_seq2_off[d - 1];
        // ```
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
    // ```c
    //             if(fence_hit && *fence_hit){
    //             	return 0;
    //             }
    // ...
    //     if (!converged)
    //         return -1;
    //
    //     if (edit_block == NULL)
    //         return best_dist;
    // ```
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
/// ```c
///     for (d = 1; d <= max_dist; d++) {
/// ```
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
/// ```c
///     for (d = 1; d <= max_dist; d++) {
/// ```
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
/// ```c
///     for (d = 1; d <= max_dist; d++) {
/// ```
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
        "LOSAT_X_GREEDYFUSE: greedy scratch smaller than the set-up allocated"
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
/// ```c
///     last_seq2_off[0][diag_origin] = seq1_index;
///     max_score[0] = seq1_index * match_cost;
///     diag_lower = diag_origin - 1;
///     diag_upper = diag_origin + 1;
///     end1_reached = end2_reached = FALSE;
///
///     /* for each distance */
///     for (d = 1; d <= max_dist; d++) {
/// ```
/// Dispatch point of LOSAT_X_GREEDYFUSE / LOSAT_X_GREEDYFUSESHADOW in `blast_greedy_align`, after
/// the lines above. `None`: the ported loop runs (switch off, or a case this module does not
/// take, module doc). `Some(value)`: the call's return value; the outputs, rows and `max_score`
/// hold what the ported loop leaves.
/// EXPERIMENT (LOSAT_X_GREEDYFUSE / LOSAT_X_GREEDYFUSESHADOW)
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
    rem: u8,
    fence_hit: &mut bool,
    seed: &mut GreedySeed,
    diag_lower: i32,
    diag_upper: i32,
    max_score_offset: usize,
) -> Option<i32> {
    let mode = x_greedy_fuse_mode();
    if mode == 0
        || !score_only
        || rem != 4
        || *fence_hit
        || max_dist < 1
        || len1 < 0
        || len2 < 0
        || len1 as usize > seq1.len()
        || len2 as usize > seq2.len()
    {
        return None;
    }
    if mode == 2 && X_SHADOW_REFERENCE.with(Cell::get) {
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
    if mode == 1 {
        let result = x_fused_on(
            reverse,
            &input,
            &mut non_affine_mem.last_seq2_off,
            &mut non_affine_mem.max_score,
            &mut out,
        );
        out.write(seq1_align_len, seq2_align_len, seed, fence_hit);
        return Some(result);
    }

    // LOSAT_X_GREEDYFUSESHADOW: this module on copies of the state the set-up left ...
    let mut rows = non_affine_mem.last_seq2_off.clone();
    let mut max_score = non_affine_mem.max_score.clone();
    let fused = x_fused_on(reverse, &input, &mut rows, &mut max_score, &mut out);
    // ... and the reference: the whole call again on the scratch itself. Its set-up writes the
    // values it wrote before this dispatch (the same inputs), so it reaches its loop from the
    // same state.
    X_SHADOW_REFERENCE.with(|flag| flag.set(true));
    let reference = blast_greedy_align(
        seq1,
        len1,
        seq2,
        len2,
        reverse,
        xdrop_threshold,
        match_cost,
        mismatch_cost,
        seq1_align_len,
        seq2_align_len,
        max_dist,
        non_affine_mem,
        None,
        rem,
        fence_hit,
        seed,
    );
    X_SHADOW_REFERENCE.with(|flag| flag.set(false));
    let reference_out = XFuseOut::read(*seq1_align_len, *seq2_align_len, seed, *fence_hit);
    assert!(
        (fused, out) == (reference, reference_out),
        "LOSAT_X_GREEDYFUSESHADOW: result differs (len1={len1} len2={len2} reverse={reverse} \
         xdrop={xdrop_threshold} costs={match_cost}/{mismatch_cost} max_dist={max_dist}): fused \
         {fused} {out:?} reference {reference} {reference_out:?}"
    );
    for (row, (fused_row, reference_row)) in rows
        .iter()
        .zip(non_affine_mem.last_seq2_off.iter())
        .enumerate()
    {
        assert!(
            fused_row == reference_row,
            "LOSAT_X_GREEDYFUSESHADOW: row {row} differs (len1={len1} len2={len2} \
             reverse={reverse} max_dist={max_dist})"
        );
    }
    assert!(
        max_score == non_affine_mem.max_score,
        "LOSAT_X_GREEDYFUSESHADOW: max_score differs (len1={len1} len2={len2} reverse={reverse} \
         max_dist={max_dist})"
    );
    X_SHADOW_COMPARED.fetch_add(1, Ordering::Relaxed);
    Some(reference)
}

#[cfg(test)]
mod tests {
    use super::super::X_GREEDY_TEST_MODE;
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
    }

    /// One call of `BLAST_GreedyAlign` without an edit block on unpacked sequences.
    struct Problem {
        seq1: Vec<u8>,
        len1: i32,
        seq2: Vec<u8>,
        len2: i32,
        reverse: bool,
        xdrop: i32,
        match_cost: i32,
        mismatch_cost: i32,
        max_dist: i32,
        /// `fence_hit` on entry (callers pass false; true must leave the ported loop running).
        fence_in: bool,
    }

    #[derive(Debug, PartialEq, Eq)]
    struct Outcome {
        ret: i32,
        seq1_align_len: i32,
        seq2_align_len: i32,
        seed: (i32, i32, i32),
        fence_hit: bool,
    }

    /// `blast_greedy_align` with the LOSAT_X_GREEDYFUSE mode `fuse_mode`, the LOSAT_X_GREEDYFAST
    /// mode `greedy_mode`, and the plain build of this module when `plain`.
    fn call(
        fuse_mode: u8,
        greedy_mode: u8,
        plain: bool,
        mem: &mut GreedyNonAffineMem,
        p: &Problem,
    ) -> Outcome {
        X_GREEDY_FUSE_TEST_MODE.with(|m| m.set(Some(fuse_mode)));
        X_GREEDY_TEST_MODE.with(|m| m.set(Some(greedy_mode)));
        X_GREEDY_FUSE_TEST_PLAIN.with(|m| m.set(plain));
        // values the call overwrites
        let mut seed = GreedySeed {
            start_q: -7,
            start_s: -7,
            match_length: -7,
        };
        let mut fence_hit = p.fence_in;
        let (mut seq1_align_len, mut seq2_align_len) = (-5, -5);
        let ret = blast_greedy_align(
            &p.seq1,
            p.len1,
            &p.seq2,
            p.len2,
            p.reverse,
            p.xdrop,
            p.match_cost,
            p.mismatch_cost,
            &mut seq1_align_len,
            &mut seq2_align_len,
            p.max_dist,
            mem,
            None,
            4,
            &mut fence_hit,
            &mut seed,
        );
        X_GREEDY_FUSE_TEST_MODE.with(|m| m.set(None));
        X_GREEDY_TEST_MODE.with(|m| m.set(None));
        X_GREEDY_FUSE_TEST_PLAIN.with(|m| m.set(false));
        Outcome {
            ret,
            seq1_align_len,
            seq2_align_len,
            seed: (seed.start_q, seed.start_s, seed.match_length),
            fence_hit,
        }
    }

    /// The memory a call leaves: both rows and `max_score`, whole.
    fn assert_same_memory(fused: &GreedyNonAffineMem, reference: &GreedyNonAffineMem, what: &str) {
        for row in 0..2 {
            assert!(
                fused.last_seq2_off[row] == reference.last_seq2_off[row],
                "{what}: row {row} differs"
            );
        }
        assert!(
            fused.max_score == reference.max_score,
            "{what}: max_score differs"
        );
    }

    /// Related sequences (60-100 % identity, substitutions and indels), ambiguity bytes in one or
    /// both, sometimes a fence byte in `seq2`, unrelated tails, bytes after `len1` / `len2` that
    /// must not be read; lengths 0-20 often, up to about 4,000.
    fn problem(rng: &mut Rng) -> Problem {
        let template_len = match rng.below(8) {
            0..=2 => rng.below(21),
            3 | 4 => 21 + rng.below(100),
            5 | 6 => 121 + rng.below(800),
            _ => 1000 + rng.below(3000),
        } as usize;
        let identity = rng.pick(&[
            1000u64, 1000, 995, 980, 950, 900, 850, 800, 750, 700, 650, 600,
        ]);
        let indel_share = rng.pick(&[0u64, 0, 100, 300, 600]); // of the differences, per mille
        let ambiguity1 = rng.pick(&[0u64, 0, 0, 5, 30]); // per mille
        let ambiguity2 = rng.pick(&[0u64, 0, 0, 5, 30]);
        let mut seq1 = Vec::new();
        let mut seq2 = Vec::new();
        for _ in 0..template_len {
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
        if rng.below(3) == 0 {
            for _ in 0..rng.below(200) {
                seq1.push(rng.below(4) as u8);
            }
            for _ in 0..rng.below(200) {
                seq2.push(rng.below(4) as u8);
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
        if !seq2.is_empty() && rng.below(5) == 0 {
            let at = rng.below(seq2.len() as u64) as usize;
            seq2[at] = FENCE_SENTRY;
        }
        let len1 = seq1.len() as i32;
        let len2 = seq2.len() as i32;
        let reverse = rng.below(2) == 1;
        if reverse {
            // the left extension reads seq[len - 1 - i]
            seq1.reverse();
            seq2.reverse();
        }
        // bytes after the lengths: matching bases and fence bytes, which a read past the end
        // would turn into a different result
        for _ in 0..rng.below(16) {
            seq1.push(rng.pick(&[0u8, 1, 2, 3, FENCE_SENTRY]));
        }
        for _ in 0..rng.below(16) {
            seq2.push(rng.pick(&[0u8, 1, 2, 3, FENCE_SENTRY]));
        }
        // megablast's costs (reward 1 doubled), and costs with match_cost / 2 != 1
        let (match_cost, mismatch_cost) = rng.pick(&[
            (2, 4),
            (2, 4),
            (2, 4),
            (2, 6),
            (2, 2),
            (4, 6),
            (6, 8),
            (4, 4),
        ]);
        let xdrop = if template_len <= 500 {
            rng.pick(&[
                0, 1, 2, 5, 10, 13, 20, 26, 27, 40, 54, 100, 300, 1000, 100_000,
            ])
        } else {
            rng.pick(&[0, 1, 5, 13, 20, 26, 27, 40, 54, 100])
        };
        // GREEDY_MAX_COST and GREEDY_MAX_COST_FRACTION as the caller uses them, and small values
        // that stop the loop at max_dist
        let initial = 1000.min(len2 / 2 + 1);
        let max_dist = match rng.below(6) {
            0 => 1 + rng.below(4) as i32,
            1 => 1 + rng.below(12) as i32,
            2 | 3 => initial,
            4 => 2 * initial,
            _ => initial.min(30) + rng.below(40) as i32,
        };
        Problem {
            seq1,
            len1,
            seq2,
            len2,
            reverse,
            xdrop,
            match_cost,
            mismatch_cost,
            max_dist,
            fence_in: rng.below(25) == 0,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:509-537,682-686
    /// ```c
    ///     for (d = 1; d <= max_dist; d++) {
    /// ...
    ///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
    /// ...
    ///     if (!converged)
    ///         return -1;
    ///
    ///     if (edit_block == NULL)
    ///         return best_dist;
    /// ```
    /// EXPERIMENT (LOSAT_X_GREEDYFUSE / LOSAT_X_GREEDYFUSESHADOW).
    ///
    /// Sequences of calls on one scratch per path, so that the cells earlier calls left matter:
    /// the ported loop (LOSAT_X_GREEDYFAST off), the fast rows (on), this module in the AVX2 build
    /// and in the plain build, and the shadow. After every call: the same return value, align
    /// lengths, seed and `fence_hit`, and the same rows and `max_score` as the ported loop.
    #[test]
    fn x_fused_greedy_matches_reference_on_call_sequences() {
        let histories: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(2000);
        let mut rng = Rng(0x7F4A_7C15_9E37_79B9);
        let (mut converged, mut stopped, mut fences, mut early, mut gated) = (0, 0, 0, 0, 0);
        let mut calls = 0usize;
        let mut unit_costs = 0usize;
        let mut other_costs = 0usize;
        X_GREEDY_FUSE_TEST_REST.with(|n| n.set((0, 0)));
        let compared_before = X_SHADOW_COMPARED.load(Ordering::Relaxed);
        for history in 0..histories {
            let mut reference = GreedyNonAffineMem::new();
            let mut fast_rows = GreedyNonAffineMem::new();
            let mut fused = GreedyNonAffineMem::new();
            let mut plain = GreedyNonAffineMem::new();
            let mut shadow = GreedyNonAffineMem::new();
            for step in 0..1 + rng.below(12) {
                let p = problem(&mut rng);
                let what = format!(
                    "history {history} step {step}: len1={} len2={} reverse={} xdrop={} costs={}/{} \
                     max_dist={} fence_in={}",
                    p.len1,
                    p.len2,
                    p.reverse,
                    p.xdrop,
                    p.match_cost,
                    p.mismatch_cost,
                    p.max_dist,
                    p.fence_in
                );
                let expected = call(0, 0, false, &mut reference, &p);
                let rows = call(0, 1, false, &mut fast_rows, &p);
                assert_eq!(rows, expected, "{what}: fast rows");
                assert_same_memory(&fast_rows, &reference, &format!("{what}: fast rows"));
                let got = call(1, 0, false, &mut fused, &p);
                assert_eq!(got, expected, "{what}: fused");
                assert_same_memory(&fused, &reference, &format!("{what}: fused"));
                let got = call(1, 0, true, &mut plain, &p);
                assert_eq!(got, expected, "{what}: fused, plain build");
                assert_same_memory(&plain, &reference, &format!("{what}: fused, plain build"));
                let got = call(2, 1, false, &mut shadow, &p);
                assert_eq!(got, expected, "{what}: shadow");
                assert_same_memory(&shadow, &reference, &format!("{what}: shadow"));

                calls += 1;
                if p.fence_in {
                    gated += 1;
                } else if expected.seq1_align_len == expected.seed.2
                    && expected.ret == 0
                    && (expected.seed.2 == p.len1 || expected.seed.2 == p.len2)
                {
                    early += 1;
                } else if expected.fence_hit {
                    fences += 1;
                } else if expected.ret < 0 {
                    stopped += 1;
                } else {
                    converged += 1;
                }
                if p.match_cost == 2 {
                    unit_costs += 1;
                } else {
                    other_costs += 1;
                }
            }
        }
        let (short, long) = X_GREEDY_FUSE_TEST_REST.with(Cell::get);
        let summary = format!(
            "calls {calls}: converged {converged}, past max_dist {stopped}, fence {fences}, early \
             return {early}, fence on entry {gated}; slides with < 8 positions left {short}, \
             slides of 8 or more {long}; match_cost 2: {unit_costs}, other: {other_costs}"
        );
        eprintln!("x_fused_greedy_matches_reference_on_call_sequences: {summary}");
        // the generator must reach every exit and every path
        assert!(converged * 3 > calls, "{summary}");
        assert!(stopped * 20 > calls, "{summary}");
        assert!(fences * 40 > calls, "{summary}");
        assert!(early > 0 && gated > 0, "{summary}");
        assert!(short > 1000 && long > 1000, "{summary}");
        assert!(other_costs * 5 > calls, "{summary}");
        assert!(
            X_SHADOW_COMPARED.load(Ordering::Relaxed) > compared_before,
            "shadow comparisons not counted"
        );
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:331-334,356-359,571-576
    /// ```c
    ///             if (seq2_index < len2 && seq2[seq2_index] == FENCE_SENTRY) {
    ///                 ASSERT(fence_hit);
    ///                 *fence_hit = TRUE;
    ///             }
    /// ...
    ///             if(fence_hit && *fence_hit){
    ///             	return 0;
    ///             }
    /// ```
    /// EXPERIMENT (LOSAT_X_GREEDYFUSE).
    ///
    /// A fence inside the loop, reached by a slide that stops within eight positions (the
    /// eight-byte test), by a slide of eight or more, and with fewer than eight positions left:
    /// 0 and `fence_hit` at once, with the rows, `max_score`, align lengths and seed the ported
    /// loop leaves, on scratch with stale cells.
    #[test]
    fn x_fused_greedy_fence_inside_the_loop() {
        // (length, substitution at, fence at, slides of the fence diagonal: (short, long))
        let cases = [
            (40usize, 5usize, 12usize, (0u64, 0u64)),
            (40, 5, 30, (0, 1)),
            (10, 3, 8, (1, 0)),
        ];
        for reverse in [false, true] {
            for &(len, substitution, fence, rest) in &cases {
                let mut seq1: Vec<u8> = (0..len).map(|i| (i * 7 % 4) as u8).collect();
                let mut seq2 = seq1.clone();
                seq2[substitution] = (seq2[substitution] + 1) & 3;
                seq2[fence] = FENCE_SENTRY;
                if reverse {
                    seq1.reverse();
                    seq2.reverse();
                }
                let p = Problem {
                    seq1,
                    len1: len as i32,
                    seq2,
                    len2: len as i32,
                    reverse,
                    xdrop: 26,
                    match_cost: 2,
                    mismatch_cost: 4,
                    max_dist: 21,
                    fence_in: false,
                };
                let mut reference = GreedyNonAffineMem::new();
                let mut fused = GreedyNonAffineMem::new();
                reference.ensure_base_capacity(48);
                fused.ensure_base_capacity(48);
                for row in 0..2 {
                    reference.last_seq2_off[row].fill(123);
                    fused.last_seq2_off[row].fill(123);
                }
                let expected = call(0, 0, false, &mut reference, &p);
                assert!(expected.fence_hit && expected.ret == 0, "{expected:?}");
                X_GREEDY_FUSE_TEST_REST.with(|n| n.set((0, 0)));
                let got = call(1, 0, false, &mut fused, &p);
                let taken = X_GREEDY_FUSE_TEST_REST.with(Cell::get);
                assert_eq!(got, expected, "len {len} fence {fence} reverse {reverse}");
                assert_same_memory(&fused, &reference, "fence");
                // the slides of distance 1 before the fence: the shifted diagonals stop early
                assert!(
                    taken.0 >= rest.0 && taken.1 == rest.1,
                    "len {len} fence {fence} reverse {reverse}: slides {taken:?}"
                );
            }
        }
    }
}
