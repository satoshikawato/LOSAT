//! The gapped kernels of a PSSM search (scalar, `Int4` scores): the start of a preliminary gapped
//! alignment, the start of a traceback, the score-only `Blast_SemiGappedAlign`, `ALIGN_EX` and
//! `BLAST_GappedAlignmentWithTraceback`, each reading the PSSM row of the absolute query position
//! (`gap_align->positionBased`). Restricted alignment is blastp's only (`blast_parameters.c:851-855`)
//! and has no PSSM version here.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:354

use super::rows::PssmRows;
use crate::algorithm::blastp::gapalign::{
    adjust_subject_range, gap_dp_reserve_band, gap_dp_reserve_initial, BlastGapDp,
};
use crate::common::GapEditOp;
// The traceback's result, which the Smith-Waterman round loop reads (one owner of the loop).
use crate::core::composition_adjustment::smith_waterman::SwXdropTraceback;

/// NCBI's `MININT` (c++/src/algo/blast/core/blast_gapalign.c:58).
const MININT: i32 = i32::MIN / 2;
/// NCBI's `HSP_MAX_WINDOW` (c++/include/algo/blast/core/blast_gapalign.h).
const HSP_MAX_WINDOW: i32 = 11;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:363-376
const SCRIPT_SUB: u8 = 3;
const SCRIPT_GAP_IN_A: u8 = 0;
const SCRIPT_GAP_IN_B: u8 = 6;
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;
/// NCBI's `FENCE_SENTRY` (c++/include/algo/blast/core/blast_gapalign.h): no protein letter has this value, so a PSSM
/// alignment never meets it; kept for the statement order of `ALIGN_EX`.
const FENCE_SENTRY: u8 = 201;

/// The memory of the PSSM gapped kernels (`gap_align->dp_mem` with its `dp_mem_alloc`, and the
/// traceback rows of `ALIGN_EX`), kept from one alignment to the next as NCBI keeps its
/// `BlastGapAlignStruct`. The DP cells grow with NCBI's window (`gap_dp_reserve_initial`,
/// `gap_dp_reserve_band`: the band that the X-drop test reaches), not with the subject's letters.
pub(crate) struct PssmGapScratch {
    dp: Vec<BlastGapDp>,
    /// NCBI's `dp_mem_alloc` (`BLAST_GapAlignStructNew` starts it at 1000).
    dp_mem_alloc: usize,
    /// The traceback rows of `ALIGN_EX`, one after the other.
    cells: Vec<u8>,
    /// For each row: its first cell in `cells` and its `edit_start_offset`.
    rows: Vec<(usize, usize)>,
    /// The traceback's operations, as `GapPrelimEditBlockAdd` merges them.
    ops: Vec<(u8, u32)>,
}

impl PssmGapScratch {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:313-319
    pub(crate) fn new() -> Self {
        Self {
            dp: Vec::new(),
            dp_mem_alloc: 1000,
            cells: Vec::new(),
            rows: Vec::new(),
            ops: Vec::new(),
        }
    }
}

/// NCBI's `num_extra_cells`: the cells that the X-drop test can reach past the band's last cell.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:797-800
#[inline(always)]
fn num_extra_cells(x_dropoff: i32, gap_extend: i32, n: usize) -> usize {
    if gap_extend > 0 {
        (x_dropoff / gap_extend + 3) as usize
    } else {
        n + 3
    }
}

impl Default for PssmGapScratch {
    fn default() -> Self {
        Self::new()
    }
}

/// The subject letter that the DP reads at `b_index`: `B[b_index + 1]` forward, `B[N - 1 -
/// b_index]` backward (`b_ptr` moves before it is read), 0 (the sentinel) outside the slice.
#[inline(always)]
fn subject_letter(subject: &[u8], base: usize, n: usize, b_index: usize, reverse: bool) -> u8 {
    let index = if reverse {
        if b_index >= n {
            return 0;
        }
        base + n - b_index - 1
    } else {
        base + b_index + 1
    };
    subject.get(index).copied().unwrap_or(0)
}

/// The PSSM row of DP row `a_index`: `pssm[M - a_index]` backward, `pssm[a_index +
/// query_offset]` forward.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:839-849
#[inline(always)]
fn pssm_row<R: PssmRows>(
    rows: &R,
    m: usize,
    a_index: usize,
    query_offset: usize,
    reverse: bool,
) -> &[i32; crate::core::blast_util::BLASTAA_SIZE] {
    if reverse {
        rows.row(m - a_index)
    } else {
        rows.row(a_index + query_offset)
    }
}

/// `BlastGetStartForGappedAlignment` with a PSSM: the middle of the best 11-letter window of an
/// ungapped HSP. `q_length` is NCBI's `Uint4` of the HSP's `Int4` length: a negative length
/// (never made by a protein search) is a huge `Uint4`, whose window loop reads only the first
/// window.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3393-3439
pub(crate) fn pssm_get_start_for_gapped_alignment<R: PssmRows>(
    rows: &R,
    subject: &[u8],
    q_start: i32,
    length: i32,
    s_start: i32,
) -> i32 {
    let letter = |index: i32| -> usize {
        usize::try_from(index)
            .ok()
            .and_then(|index| subject.get(index))
            .copied()
            .unwrap_or(0) as usize
    };
    let row = |pos: i32| rows.row(pos as usize);
    if (0..=HSP_MAX_WINDOW).contains(&length) {
        return q_start + length / 2;
    }
    let mut score = 0i32;
    for i in 0..HSP_MAX_WINDOW {
        score += row(q_start + i)[letter(s_start + i)];
    }
    let mut max_score = score;
    let mut max_offset = q_start + HSP_MAX_WINDOW - 1;
    if length > 0 {
        // The ungapped HSP's query and subject lengths are the same.
        let hsp_end = q_start + length;
        for index1 in q_start + HSP_MAX_WINDOW..hsp_end {
            let s_index = s_start + (index1 - q_start);
            score -= row(index1 - HSP_MAX_WINDOW)[letter(s_index - HSP_MAX_WINDOW)];
            score += row(index1)[letter(s_index)];
            if score > max_score {
                max_score = score;
                max_offset = index1;
            }
        }
    }
    if max_score > 0 {
        max_offset - HSP_MAX_WINDOW / 2
    } else {
        q_start
    }
}

/// `BlastGetOffsetsForGappedAlignment` with a PSSM: the traceback's start for an HSP whose gapped
/// start is not set, or `None` when no window scores above 0. The windows read 11 letters from
/// the HSP's starts and ends; a window outside the subject gives no start (as for BLASTP, this is
/// not reached: a preliminary gapped start is never both 0).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3248-3320
pub(crate) fn pssm_get_offsets_for_gapped_alignment<R: PssmRows>(
    rows: &R,
    subject: &[u8],
    query_start: i32,
    query_end: i32,
    subject_start: i32,
    subject_end: i32,
) -> Option<(i32, i32)> {
    let q_length = query_end - query_start;
    let s_length = subject_end - subject_start;
    if q_length <= HSP_MAX_WINDOW {
        return Some((query_start + q_length / 2, subject_start + q_length / 2));
    }
    let score_at = |q: i32, s: i32| -> Option<i32> {
        let row = (q >= 0 && (q as usize) < rows.len()).then(|| rows.row(q as usize))?;
        let s = subject.get(usize::try_from(s).ok()?)?;
        if rows.reads_by_letter(row) {
            rows.letters_read(row, std::iter::once(*s));
        }
        Some(row[*s as usize])
    };
    let mut score = 0;
    for i in 0..HSP_MAX_WINDOW {
        score += score_at(query_start + i, subject_start + i)?;
    }
    let mut max_score = score;
    let mut max_offset = query_start + HSP_MAX_WINDOW - 1;
    for index1 in query_start + HSP_MAX_WINDOW..query_start + q_length.min(s_length) {
        let s_index = index1 - query_start + subject_start;
        score -= score_at(index1 - HSP_MAX_WINDOW, s_index - HSP_MAX_WINDOW)?;
        score += score_at(index1, s_index)?;
        if score > max_score {
            max_score = score;
            max_offset = index1;
        }
    }
    if max_score > 0 {
        return Some((max_offset, max_offset - query_start + subject_start));
    }
    let mut end_score = 0;
    for i in 0..HSP_MAX_WINDOW {
        end_score += score_at(
            query_end - HSP_MAX_WINDOW + i,
            subject_start + s_length - HSP_MAX_WINDOW + i,
        )?;
    }
    (end_score > 0).then_some((
        query_end - HSP_MAX_WINDOW / 2,
        subject_end - HSP_MAX_WINDOW / 2,
    ))
}

/// The score-only `Blast_SemiGappedAlign` with a PSSM: (`*a_offset`, `*b_offset`, the best
/// score). Backward (`reverse_sequence`) the rows are `pssm[M - a_index]` and the subject is
/// read from `s_base + n - 1` down; forward they are `pssm[a_index + query_offset]` and the
/// subject is read from `s_base + 1` up.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:776-985
/// The direction is a constant (`REVERSE`, one function for each), and each row reads its columns
/// in three runs (series PR, item PS12s, D5): before, inside and after the columns whose subject
/// letter is in `subject` (`subject_columns`, once per call), the runs outside reading 0, the
/// value `subject_letter` gives there (the sentinel, NCBI's read past it). The statements, their
/// order and the values read are those above; the cells' reads need no bounds checks.
// Out of line, as before PS12s (one function for each direction): inlined into a caller's loop
// (RPS-BLAST's whole preliminary search) its registers spill.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn semi_gapped_align_score_only<R: PssmRows, const REVERSE: bool>(
    rows: &R,
    m: usize,
    query_offset: usize,
    subject: &[u8],
    s_base: usize,
    n: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut PssmGapScratch,
) -> (usize, usize, i32) {
    let gap_open_extend = gap_open + gap_extend;
    let x_dropoff = x_drop.max(gap_open_extend);
    if n == 0 || m == 0 {
        return (0, 0, 0);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:797-807
    let extra = num_extra_cells(x_dropoff, gap_extend, n);
    let PssmGapScratch {
        dp, dp_mem_alloc, ..
    } = scratch;
    gap_dp_reserve_initial(dp, dp_mem_alloc, extra, n);
    // The subject letters of the DP's columns, taken once (D5 of PS12s): column `b_index` reads
    // `subject_letter(subject, s_base, n, b_index, REVERSE)`, which is a letter of `subject` for
    // the columns `inside` and 0 (the sentinel's value, NCBI's read past it) for the others.
    let (inside, letters) = subject_columns::<REVERSE>(subject, s_base, n);
    let mut score = -gap_open_extend;
    dp[0].best = 0;
    dp[0].best_gap = -gap_open_extend;
    let mut i = 1usize;
    while i <= n {
        if score < -x_dropoff {
            break;
        }
        dp[i].best = score;
        dp[i].best_gap = score - gap_open_extend;
        score -= gap_extend;
        i += 1;
    }
    let mut b_size = i;
    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;
    for a_index in 1..=m {
        let matrix_row = pssm_row(rows, m, a_index, query_offset, REVERSE);
        let mut score = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;
        // The column of the row's last new best score (NCBI's `*a_offset = a_index; *b_offset =
        // b_index` at each one), `usize::MAX` while the row has none.
        let mut row_best = usize::MAX;
        // The row's columns `first_b_index..b_size` in order, in three runs: before the columns
        // whose letter is in `subject`, those, after them.
        let start = first_b_index;
        let end = b_size.max(start);
        let lo = inside.start.clamp(start, end);
        let hi = inside.end.clamp(lo, end);
        {
            // One cell of the row (the statements of NCBI's inner loop, in its order).
            let mut cell = |b_index: usize, dp_cell: &mut BlastGapDp, letter: u8| {
                let mut score_gap_col = dp_cell.best_gap;
                let next_score = dp_cell.best + matrix_row[letter as usize];
                if score < score_gap_col {
                    score = score_gap_col;
                }
                if score < score_gap_row {
                    score = score_gap_row;
                }
                if best_score - score > x_dropoff {
                    if b_index == first_b_index {
                        first_b_index += 1;
                    } else {
                        dp_cell.best = MININT;
                    }
                } else {
                    last_b_index = b_index;
                    if score > best_score {
                        best_score = score;
                        row_best = b_index;
                    }
                    score_gap_row -= gap_extend;
                    score_gap_col -= gap_extend;
                    dp_cell.best_gap = (score - gap_open_extend).max(score_gap_col);
                    score_gap_row = (score - gap_open_extend).max(score_gap_row);
                    dp_cell.best = score;
                }
                score = next_score;
            };
            for (b_index, dp_cell) in (start..lo).zip(&mut dp[start..lo]) {
                cell(b_index, dp_cell, 0);
            }
            if lo < hi {
                let width = hi - lo;
                let cells = &mut dp[lo..hi];
                if REVERSE {
                    // Column `b` reads `letters[n - 1 - b]`: the run reads its letters backward.
                    let run = &letters[n - hi..n - lo];
                    for k in 0..width {
                        cell(lo + k, &mut cells[k], run[width - 1 - k]);
                    }
                } else {
                    let run = &letters[lo..hi];
                    for k in 0..width {
                        cell(lo + k, &mut cells[k], run[k]);
                    }
                }
            }
            for (b_index, dp_cell) in (hi..end).zip(&mut dp[hi..end]) {
                cell(b_index, dp_cell, 0);
            }
        }
        if row_best != usize::MAX {
            a_offset = a_index;
            b_offset = row_best;
        }
        if first_b_index == b_size {
            break;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:923-931
        gap_dp_reserve_band(dp, dp_mem_alloc, last_b_index, extra, n);
        if last_b_index + 1 < b_size {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= (best_score - x_dropoff) && b_size <= n {
                dp[b_size].best = score_gap_row;
                dp[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size += 1;
            }
        }
        if b_size <= n {
            dp[b_size].best = MININT;
            dp[b_size].best_gap = MININT;
            b_size += 1;
        }
    }
    (a_offset, b_offset, best_score)
}

/// The DP columns of `semi_gapped_align_score_only` whose subject letter is in `subject`, and
/// the letters they read; every other column reads 0 (`subject_letter`). Forward, column `b` (0
/// to N) reads `subject[base + b + 1]`: the columns `0..hi` while that index is in `subject`,
/// column `b` reading `letters[b]` (`letters` from `base + 1`). Backward, column `b` reads
/// `subject[base + N - 1 - b]` for `b` below N and 0 at N (NCBI's `b_ptr` before the first
/// letter): the columns `lo..N` whose index is in `subject`, column `b` reading
/// `letters[N - 1 - b]` (`letters` from `base`, `N - lo` of them).
fn subject_columns<const REVERSE: bool>(
    subject: &[u8],
    base: usize,
    n: usize,
) -> (std::ops::Range<usize>, &[u8]) {
    if REVERSE {
        let lo = (base + n).saturating_sub(subject.len()).min(n);
        match subject.get(base..base + n - lo) {
            Some(letters) => (lo..n, letters),
            None => (n..n, &[]),
        }
    } else {
        let hi = (n + 1).min(subject.len().saturating_sub(base + 1));
        match subject.get(base + 1..base + 1 + hi) {
            Some(letters) => (0..hi, letters),
            None => (0..0, &[]),
        }
    }
}

/// A score-only gapped alignment: the alignment's ends and score.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PssmScoreOnly {
    pub(crate) query_start: i32,
    pub(crate) query_stop: i32,
    pub(crate) subject_start: i32,
    pub(crate) subject_stop: i32,
    pub(crate) score: i32,
}

/// `s_BlastProtGappedAlignment` with a PSSM, without restricted alignment: the score-only
/// extensions to the left (including the start) and to the right of `(q_off, s_off)`. `subject`
/// is `subject->sequence` (its letters, then its trailing sentinel); `subject_length` is
/// `subject->length`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4245-4347
#[allow(clippy::too_many_arguments)]
pub(crate) fn pssm_score_only_gapped_alignment<R: PssmRows>(
    rows: &R,
    subject: &[u8],
    subject_length: usize,
    q_off: usize,
    s_off: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut PssmGapScratch,
) -> PssmScoreOnly {
    let query_length = rows.len();
    let q_length = q_off + 1;
    let (s_length, adjusted_length, subject_shift) =
        adjust_subject_range(s_off + 1, subject_length, q_length, query_length);
    // `q_length` and `s_length` are at least 1 here.
    let (private_q_start, private_s_start, score_left) = semi_gapped_align_score_only::<R, true>(
        rows,
        q_length,
        q_off,
        subject,
        subject_shift,
        s_length,
        gap_open,
        gap_extend,
        x_drop,
        scratch,
    );
    let query_start = (q_length - private_q_start) as i32;
    let subject_start = (s_length - private_s_start + subject_shift) as i32;
    let (query_stop, subject_stop, score_right) =
        if q_length < query_length && s_length < adjusted_length {
            let (q_stop, s_stop, score_right) = semi_gapped_align_score_only::<R, false>(
                rows,
                query_length - q_length,
                q_off,
                subject,
                s_off,
                adjusted_length - s_length,
                gap_open,
                gap_extend,
                x_drop,
                scratch,
            );
            (
                (q_stop + q_off + 1) as i32,
                (s_stop + s_off + 1) as i32,
                score_right,
            )
        } else {
            (q_off as i32, s_off as i32, 0)
        };
    PssmScoreOnly {
        query_start,
        query_stop,
        subject_start,
        subject_stop,
        score: score_left + score_right,
    }
}

/// `ALIGN_EX` with a PSSM: (`*a_offset`, `*b_offset`, the best score, the traceback's operations
/// as `GapPrelimEditBlockAdd` writes them, from the best cell back to the start).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:420-735
#[allow(clippy::too_many_arguments)]
fn align_ex<R: PssmRows>(
    rows: &R,
    m: usize,
    query_offset: usize,
    subject: &[u8],
    s_base: usize,
    n: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    scratch: &mut PssmGapScratch,
    fence_hit: &mut bool,
) -> (usize, usize, i32, Vec<GapEditOp>) {
    let gap_open_extend = gap_open + gap_extend;
    let x_dropoff = x_drop.max(gap_open_extend);
    if n == 0 || m == 0 {
        return (0, 0, 0, Vec::new());
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:457-468
    let extra = num_extra_cells(x_dropoff, gap_extend, n);
    let PssmGapScratch {
        dp,
        dp_mem_alloc,
        cells,
        rows: trace_rows,
        ops,
    } = scratch;
    gap_dp_reserve_initial(dp, dp_mem_alloc, extra, n);
    cells.clear();
    trace_rows.clear();
    // Row 0 (`edit_script[0]`): its cells 1.. are gaps in A; cell 0 is never read.
    trace_rows.push((0, 0));
    cells.push(SCRIPT_GAP_IN_A);
    let mut score = -gap_open_extend;
    dp[0].best = 0;
    dp[0].best_gap = -gap_open_extend;
    let mut i = 1usize;
    while i <= n {
        if score < -x_dropoff {
            break;
        }
        dp[i].best = score;
        dp[i].best_gap = score - gap_open_extend;
        score -= gap_extend;
        cells.push(SCRIPT_GAP_IN_A);
        i += 1;
    }
    let mut b_size = i;
    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;
    for a_index in 1..=m {
        let row_start = cells.len();
        let orig_b_index = first_b_index;
        trace_rows.push((row_start, orig_b_index));
        let matrix_row = pssm_row(rows, m, a_index, query_offset, reverse);
        let mut score = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;
        let mut b_index = first_b_index;
        while b_index < b_size {
            let letter = subject_letter(subject, s_base, n, b_index, reverse);
            let mut score_gap_col = dp[b_index].best_gap;
            if letter == FENCE_SENTRY {
                *fence_hit = true;
                break;
            }
            let next_score = dp[b_index].best + matrix_row[letter as usize];
            let mut script = SCRIPT_SUB;
            if score < score_gap_col {
                script = SCRIPT_GAP_IN_B;
                score = score_gap_col;
            }
            if score < score_gap_row {
                script = SCRIPT_GAP_IN_A;
                score = score_gap_row;
            }
            if best_score - score > x_dropoff {
                if first_b_index == b_index {
                    first_b_index += 1;
                } else {
                    dp[b_index].best = MININT;
                }
            } else {
                last_b_index = b_index;
                if score > best_score {
                    best_score = score;
                    a_offset = a_index;
                    b_offset = b_index;
                }
                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                if score_gap_col < (score - gap_open_extend) {
                    dp[b_index].best_gap = score - gap_open_extend;
                } else {
                    dp[b_index].best_gap = score_gap_col;
                    script += SCRIPT_EXTEND_GAP_B;
                }
                if score_gap_row < (score - gap_open_extend) {
                    score_gap_row = score - gap_open_extend;
                } else {
                    script += SCRIPT_EXTEND_GAP_A;
                }
                dp[b_index].best = score;
            }
            score = next_score;
            cells.push(script);
            b_index += 1;
        }
        // The letters of the cells this row read (before the band's end or the fence), for a row
        // whose letters are read one by one (`PssmRows::reads_by_letter`).
        if rows.reads_by_letter(matrix_row) {
            rows.letters_read(
                matrix_row,
                (orig_b_index..b_index).map(|b| subject_letter(subject, s_base, n, b, reverse)),
            );
        }
        if first_b_index == b_size || *fence_hit {
            break;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:641-649
        gap_dp_reserve_band(dp, dp_mem_alloc, last_b_index, extra, n);
        if last_b_index + 1 < b_size {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= (best_score - x_dropoff) && b_size <= n {
                dp[b_size].best = score_gap_row;
                dp[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                // `edit_script_row[b_size]`: the row's cells run on from the main loop's.
                let at = row_start + (b_size - orig_b_index);
                if cells.len() <= at {
                    cells.resize(at + 1, SCRIPT_GAP_IN_A);
                }
                cells[at] = SCRIPT_GAP_IN_A;
                b_size += 1;
            }
        }
        if b_size <= n {
            dp[b_size].best = MININT;
            dp[b_size].best_gap = MININT;
            b_size += 1;
        }
    }
    ops.clear();
    if *fence_hit {
        return (a_offset, b_offset, best_score, Vec::new());
    }
    let mut a_index = a_offset;
    let mut b_index = b_offset;
    let mut script = SCRIPT_SUB;
    while a_index > 0 || b_index > 0 {
        let (start, first) = trace_rows[a_index];
        let next_script = cells[start + (b_index - first)];
        script = match script {
            SCRIPT_GAP_IN_A if next_script & SCRIPT_EXTEND_GAP_A != 0 => SCRIPT_GAP_IN_A,
            SCRIPT_GAP_IN_B if next_script & SCRIPT_EXTEND_GAP_B != 0 => SCRIPT_GAP_IN_B,
            _ => next_script & SCRIPT_OP_MASK,
        };
        if script == SCRIPT_GAP_IN_A {
            b_index -= 1;
        } else if script == SCRIPT_GAP_IN_B {
            a_index -= 1;
        } else {
            a_index -= 1;
            b_index -= 1;
        }
        // GapPrelimEditBlockAdd: one more of the last operation, or a new one.
        match ops.last_mut() {
            Some(last) if last.0 == script => last.1 += 1,
            _ => ops.push((script, 1)),
        }
    }
    let edit_ops = ops
        .drain(..)
        .map(|(op, count)| match op {
            SCRIPT_GAP_IN_A => GapEditOp::Del(count),
            SCRIPT_GAP_IN_B => GapEditOp::Ins(count),
            _ => GapEditOp::Sub(count),
        })
        .collect();
    (a_offset, b_offset, best_score, edit_ops)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:671-736
/// PHI-BLAST's score-only extension (`phi_search::align`): `semi_gapped_align_score_only` with the
/// rows of the query's score matrix, whose cells NCBI's matrix branch reads when the call's query
/// offset is its query pointer's (DESIGN_ps8 §2.1). Returns (`*a_offset`, `*b_offset`, score).
#[allow(clippy::too_many_arguments)]
pub(crate) fn phi_score_only(
    rows: &super::rows::QueryPssm,
    m: usize,
    query_offset: usize,
    subject: &[u8],
    s_base: usize,
    n: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    scratch: &mut PssmGapScratch,
) -> (usize, usize, i32) {
    if reverse {
        semi_gapped_align_score_only::<super::rows::QueryPssm, true>(
            rows,
            m,
            query_offset,
            subject,
            s_base,
            n,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
        )
    } else {
        semi_gapped_align_score_only::<super::rows::QueryPssm, false>(
            rows,
            m,
            query_offset,
            subject,
            s_base,
            n,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:831-906
/// PHI-BLAST's `ALIGN_EX` (`phi_search::align`), with the query's matrix rows; no fence.
#[allow(clippy::too_many_arguments)]
pub(crate) fn phi_align_ex(
    rows: &super::rows::QueryPssm,
    m: usize,
    query_offset: usize,
    subject: &[u8],
    s_base: usize,
    n: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    scratch: &mut PssmGapScratch,
) -> (usize, usize, i32, Vec<GapEditOp>) {
    let mut fence_hit = false;
    align_ex(
        rows,
        m,
        query_offset,
        subject,
        s_base,
        n,
        gap_open,
        gap_extend,
        x_drop,
        reverse,
        scratch,
        &mut fence_hit,
    )
}

/// PHI-BLAST's edit script of its left operations (with the pattern's appended) and right ones.
pub(crate) fn phi_edit_script(rev: Vec<GapEditOp>, fwd: Vec<GapEditOp>) -> Vec<GapEditOp> {
    prelim_edit_blocks_to_script(rev, fwd)
}

/// `Blast_PrelimEditBlockToGapEditScript`: the left extension's operations as they are, then the
/// right extension's reversed, merging the two last operations when they are of one kind.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2494-2534
fn prelim_edit_blocks_to_script(rev: Vec<GapEditOp>, fwd: Vec<GapEditOp>) -> Vec<GapEditOp> {
    let mut esp = rev;
    let merge = match (esp.last_mut(), fwd.last()) {
        (Some(GapEditOp::Sub(a)), Some(GapEditOp::Sub(b)))
        | (Some(GapEditOp::Del(a)), Some(GapEditOp::Del(b)))
        | (Some(GapEditOp::Ins(a)), Some(GapEditOp::Ins(b))) => {
            *a += *b;
            true
        }
        _ => false,
    };
    esp.extend(fwd.into_iter().rev().skip(usize::from(merge)));
    esp
}

/// The result of a traceback: the alignment's ends, score and edit script.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PssmTraceback {
    pub(crate) query_start: i32,
    pub(crate) query_stop: i32,
    pub(crate) subject_start: i32,
    pub(crate) subject_stop: i32,
    pub(crate) score: i32,
    pub(crate) edit_script: Vec<GapEditOp>,
    /// `*fence_hit`: the alignment read a `FENCE_SENTRY` of a partly translated subject (tblastn
    /// `-in_pssm`, whose traceback then translates the subject whole; a protein subject has
    /// none, and psiblast does not read it).
    pub(crate) fence_hit: bool,
}

/// `BLAST_GappedAlignmentWithTraceback` with a PSSM (not out-of-frame): `ALIGN_EX` to the left
/// (including the start) and to the right of `(q_start, s_start)`, the edit script, and the
/// pruning of leading and trailing gaps. `subject` is the adjusted subject (its letters from the
/// traceback's start shift; `subject_length` of them, then the bytes after them), `q_start` an
/// absolute query position.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4599-4718
#[allow(clippy::too_many_arguments)]
pub(crate) fn pssm_gapped_alignment_with_traceback<R: PssmRows>(
    rows: &R,
    subject: &[u8],
    subject_length: usize,
    q_start: usize,
    s_start: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut PssmGapScratch,
) -> PssmTraceback {
    let query_length = rows.len();
    let mut fence_hit = false;
    let (left_q, left_s, mut score_left, left_ops) = align_ex(
        rows,
        q_start + 1,
        q_start,
        subject,
        0,
        s_start + 1,
        gap_open,
        gap_extend,
        x_drop,
        true,
        scratch,
        &mut fence_hit,
    );
    let mut query_start = q_start as i32 - left_q as i32 + 1;
    let mut subject_start = s_start as i32 - left_s as i32 + 1;
    let (mut query_stop, mut subject_stop, mut score_right, right_ops) =
        if !fence_hit && q_start < query_length && s_start < subject_length {
            let (right_q, right_s, score_right, right_ops) = align_ex(
                rows,
                query_length - q_start - 1,
                q_start,
                subject,
                s_start,
                subject_length - s_start - 1,
                gap_open,
                gap_extend,
                x_drop,
                false,
                scratch,
                &mut fence_hit,
            );
            (
                q_start as i32 + right_q as i32 + 1,
                s_start as i32 + right_s as i32 + 1,
                score_right,
                right_ops,
            )
        } else {
            (q_start as i32 - 1, s_start as i32 - 1, 0, Vec::new())
        };
    let mut edit_script = prelim_edit_blocks_to_script(left_ops, right_ops);
    let mut leading = 0usize;
    while leading < edit_script.len() {
        match edit_script[leading] {
            GapEditOp::Sub(_) => break,
            GapEditOp::Del(count) => {
                score_left += gap_open + count as i32 * gap_extend;
                subject_start += count as i32;
            }
            GapEditOp::Ins(count) => {
                score_left += gap_open + count as i32 * gap_extend;
                query_start += count as i32;
            }
        }
        leading += 1;
    }
    edit_script.drain(..leading);
    while let Some(&last) = edit_script.last() {
        match last {
            GapEditOp::Sub(_) => break,
            GapEditOp::Del(count) => {
                score_right += gap_open + count as i32 * gap_extend;
                subject_stop -= count as i32;
            }
            GapEditOp::Ins(count) => {
                score_right += gap_open + count as i32 * gap_extend;
                query_stop -= count as i32;
            }
        }
        edit_script.pop();
    }
    PssmTraceback {
        query_start,
        query_stop,
        subject_start,
        subject_stop,
        score: score_left + score_right,
        edit_script,
        fence_hit,
    }
}

/// The rows of a PSSM as `ALIGN_EX` reads them from a `query_offset` of one less than
/// `queryStart` (`s_SWFindFinalEndsUsingXdrop`): its row `a_index` (from 1) is the PSSM's row
/// `start + a_index - 1`.
struct RowsFrom<'r, R> {
    rows: &'r R,
    start: usize,
}

impl<R: PssmRows> PssmRows for RowsFrom<'_, R> {
    #[inline(always)]
    fn row(&self, pos: usize) -> &[i32; crate::core::blast_util::BLASTAA_SIZE] {
        self.rows.row(self.start + pos - 1)
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.rows.len() + 1 - self.start
    }

    #[inline(always)]
    fn reads_by_letter(&self, row: &[i32; crate::core::blast_util::BLASTAA_SIZE]) -> bool {
        self.rows.reads_by_letter(row)
    }

    #[inline(always)]
    fn letters_read(
        &self,
        row: &[i32; crate::core::blast_util::BLASTAA_SIZE],
        letters: impl Iterator<Item = u8>,
    ) {
        self.rows.letters_read(row, letters);
    }
}

/// `s_SWFindFinalEndsUsingXdrop` with a PSSM, as `s_NewAlignmentUsingXdrop` calls it: the forward
/// `ALIGN_EX` from the Smith-Waterman alignment's start `(query_start, match_start)` over the
/// rectangle that ends at `(query_end, match_end)` (both included), its X-drop doubled up to two
/// times while the score stays below the Smith-Waterman `score`, and the edit script of the last
/// try (`Blast_PrelimEditBlockToGapEditScript` of an empty left part). `subject` is NCBI's
/// `subject->data` less one: the letters after one sentinel byte.
///
/// NCBI's doubling of `gap_x_dropoff` overflows `Int4` above 2^30 (undefined behavior); the doubled
/// value is kept here as its mathematical value, which saturates at `i32::MAX` without a change:
/// `ALIGN_EX` compares it with differences of `Int4` scores only.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:844-883
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1849-1857
#[allow(clippy::too_many_arguments)]
pub(crate) fn pssm_sw_xdrop_traceback<R: PssmRows>(
    rows: &R,
    subject: &[u8],
    (query_start, query_end): (usize, usize),
    (match_start, match_end): (usize, usize),
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    score: i32,
    scratch: &mut PssmGapScratch,
) -> SwXdropTraceback {
    let rows = RowsFrom {
        rows,
        start: query_start,
    };
    let mut x_drop = x_dropoff;
    let mut doubling_count = 0;
    let mut fence_hit = false;
    loop {
        let (a_offset, b_offset, xdrop_align_score, forward) = align_ex(
            &rows,
            query_end - query_start + 1,
            0,
            subject,
            match_start,
            match_end - match_start + 1,
            gap_open,
            gap_extend,
            x_drop,
            false,
            scratch,
            &mut fence_hit,
        );
        x_drop = x_drop.saturating_mul(2);
        doubling_count += 1;
        if !(xdrop_align_score < score && doubling_count < 3) {
            return SwXdropTraceback {
                query_extent: a_offset as i32,
                match_extent: b_offset as i32,
                score: xdrop_align_score,
                edit_script: prelim_edit_blocks_to_script(Vec::new(), forward),
            };
        }
    }
}

#[cfg(test)]
#[path = "gapalign_tests.rs"]
mod tests;
