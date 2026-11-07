//! The gapped kernels of a PSSM search (scalar, `Int4` scores): the start of a preliminary gapped
//! alignment, the start of a traceback, the score-only `Blast_SemiGappedAlign`, `ALIGN_EX` and
//! `BLAST_GappedAlignmentWithTraceback`, each reading the PSSM row of the absolute query position
//! (`gap_align->positionBased`). Restricted alignment is blastp's only (`blast_parameters.c:851-855`)
//! and has no PSSM version here.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:354
//! ```c
//!    gap_align->positionBased = (sbp->psi_matrix != NULL);
//! ```

use super::rows::PssmRows;
use crate::algorithm::blastp::gapalign::adjust_subject_range;
use crate::common::GapEditOp;

/// NCBI's `MININT` (blast_gapalign.c:58).
const MININT: i32 = i32::MIN / 2;
/// NCBI's `HSP_MAX_WINDOW` (blast_gapalign.h).
const HSP_MAX_WINDOW: i32 = 11;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:368-376
// ```c
// enum {
//     SCRIPT_SUB           = eGapAlignSub,     /**< Substitution */
//     SCRIPT_GAP_IN_A      = eGapAlignDel,     /**< Deletion */
//     SCRIPT_GAP_IN_B      = eGapAlignIns,     /**< Insertion */
//     SCRIPT_OP_MASK       = 0x07, /**< Mask for edit script operations */
//
//     SCRIPT_EXTEND_GAP_A  = 0x10, /**< continue a gap in A */
//     SCRIPT_EXTEND_GAP_B  = 0x40  /**< continue a gap in B */
// };
// ```
const SCRIPT_SUB: u8 = 3;
const SCRIPT_GAP_IN_A: u8 = 0;
const SCRIPT_GAP_IN_B: u8 = 6;
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;
/// NCBI's `FENCE_SENTRY` (blast_gapalign.h): no protein letter has this value, so a PSSM
/// alignment never meets it; kept for the statement order of `ALIGN_EX`.
const FENCE_SENTRY: u8 = 201;

/// The DP cell of NCBI's `BlastGapDP`.
#[derive(Clone, Copy)]
struct GapDp {
    best: i32,
    best_gap: i32,
}

/// The memory of the PSSM gapped kernels (`gap_align->dp_mem` and the traceback rows of
/// `ALIGN_EX`), kept from one alignment to the next as NCBI keeps its `BlastGapAlignStruct`.
/// NCBI's `dp_mem` grows with `dp_mem_alloc` and its cells are written before they are read;
/// here the vector holds `N + 2` cells of each extension (the most the DP addresses), so the
/// values are the same and only the allocation differs.
pub(crate) struct PssmGapScratch {
    dp: Vec<GapDp>,
    /// The traceback rows of `ALIGN_EX`, one after the other.
    cells: Vec<u8>,
    /// For each row: its first cell in `cells` and its `edit_start_offset`.
    rows: Vec<(usize, usize)>,
    /// The traceback's operations, as `GapPrelimEditBlockAdd` merges them.
    ops: Vec<(u8, u32)>,
}

impl PssmGapScratch {
    pub(crate) fn new() -> Self {
        Self {
            dp: Vec::new(),
            cells: Vec::new(),
            rows: Vec::new(),
            ops: Vec::new(),
        }
    }

    fn reserve_dp(&mut self, n: usize) {
        if self.dp.len() < n + 2 {
            self.dp.resize(
                n + 2,
                GapDp {
                    best: MININT,
                    best_gap: MININT,
                },
            );
        }
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
/// ```c
///         if (!(gap_align->positionBased)) {
///             ...
///         }
///         else {
///             if(reversed || reverse_sequence)
///                 matrix_row = pssm[M - a_index];
///             else
///                 matrix_row = pssm[a_index + query_offset];
///         }
/// ```
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
/// ```c
///     if (q_length <= HSP_MAX_WINDOW) {
///         max_offset = q_start + q_length/2;
///         return max_offset;
///     }
///
///     hsp_end = q_start + HSP_MAX_WINDOW;
///     query_var = query + q_start;
///     subject_var = subject + s_start;
///     score=0;
///     for (index1=q_start; index1<hsp_end; index1++) {
///         if (!(positionBased))
///             score += sbp->matrix->data[*query_var][*subject_var];
///         else
///             score += sbp->psi_matrix->pssm->data[index1][*subject_var];
///         query_var++; subject_var++;
///     }
///     max_score = score;
///     max_offset = hsp_end - 1;
///     hsp_end = q_start + MIN(q_length, s_length);
///     for (index1=q_start + HSP_MAX_WINDOW; index1<hsp_end; index1++) {
///         if (!(positionBased)) {
///             ...
///         } else {
///             score -= sbp->psi_matrix->pssm->data[index1-HSP_MAX_WINDOW][*(subject_var-HSP_MAX_WINDOW)];
///             score += sbp->psi_matrix->pssm->data[index1][*subject_var];
///         }
///         if (score > max_score) {
///             max_score = score;
///             max_offset = index1;
///         }
///         query_var++; subject_var++;
///     }
///     if (max_score > 0)
///        max_offset -= HSP_MAX_WINDOW/2;
///     else
///        max_offset = q_start;
///
///     return max_offset;
/// ```
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
/// ```c
///     if (q_length <= HSP_MAX_WINDOW) {
///         *q_retval = q_start + q_length/2;
///         *s_retval = s_start + q_length/2;
///         return TRUE;
///     }
///
///     hsp_end = q_start + HSP_MAX_WINDOW;
///     ...
///         else
///             score += sbp->psi_matrix->pssm->data[index1][*subject_var];
///     ...
///     hsp_end = q_start + MIN(q_length, s_length);
///     for (index1=q_start + HSP_MAX_WINDOW; index1<hsp_end; index1++) {
///         ...
///         } else {
///             score -= sbp->psi_matrix->pssm->data[index1-HSP_MAX_WINDOW][*(subject_var-HSP_MAX_WINDOW)];
///             score += sbp->psi_matrix->pssm->data[index1][*subject_var];
///         }
///     ...
///     if (max_score > 0)
///     {
///         *q_retval = max_offset;
///         *s_retval = (max_offset - q_start) + s_start;
///         return TRUE;
///     }
///     else  /* Test the window around the ends of the HSP. */
///     {
///         score=0;
///         query_var = query + q_start + q_length - HSP_MAX_WINDOW;
///         subject_var = subject + s_start + s_length - HSP_MAX_WINDOW;
///         for (index1=hsp->query.end-HSP_MAX_WINDOW; index1<hsp->query.end; index1++) {
///             if (!(positionBased))
///                 score += sbp->matrix->data[*query_var][*subject_var];
///             else
///                 score += sbp->psi_matrix->pssm->data[index1][*subject_var];
///             query_var++; subject_var++;
///         }
///         if (score > 0)
///         {
///             *q_retval = hsp->query.end - HSP_MAX_WINDOW/2;
///             *s_retval = hsp->subject.end - HSP_MAX_WINDOW/2;
///             return TRUE;
///         }
///     }
///     return FALSE;
/// ```
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
/// ```c
///     *a_offset = 0;
///     *b_offset = 0;
///     gap_open = score_params->gap_open;
///     gap_extend = score_params->gap_extend;
///     gap_open_extend = gap_open + gap_extend;
///     x_dropoff = gap_align->gap_x_dropoff;
///
///     if (x_dropoff < gap_open_extend)
///         x_dropoff = gap_open_extend;
///
///     if(N <= 0 || M <= 0)
///         return 0;
///     ...
///     if (gap_extend > 0)
///         num_extra_cells = x_dropoff / gap_extend + 3;
///     else
///         num_extra_cells = N + 3;
///     ...
///     score = -gap_open_extend;
///     score_array[0].best = 0;
///     score_array[0].best_gap = -gap_open_extend;
///
///     for (i = 1; i <= N; i++) {
///         if (score < -x_dropoff)
///             break;
///
///         score_array[i].best = score;
///         score_array[i].best_gap = score - gap_open_extend;
///         score -= gap_extend;
///     }
///
///     b_size = i;
///     best_score = 0;
///     first_b_index = 0;
///     ...
///     for (a_index = 1; a_index <= M; a_index++) {
///         ...
///         score = MININT;
///         score_gap_row = MININT;
///         last_b_index = first_b_index;
///
///         for (b_index = first_b_index; b_index < b_size; b_index++) {
///
///             b_ptr += b_increment;
///             score_gap_col = score_array[b_index].best_gap;
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
///
///             if (score < score_gap_col)
///                 score = score_gap_col;
///
///             if (score < score_gap_row)
///                 score = score_gap_row;
///
///             if (best_score - score > x_dropoff) {
///                 if (b_index == first_b_index)
///                     first_b_index++;
///                 else
///                     score_array[b_index].best = MININT;
///             }
///             else {
///                 last_b_index = b_index;
///                 if (score > best_score) {
///                     best_score = score;
///                     *a_offset = a_index;
///                     *b_offset = b_index;
///                 }
///                 score_gap_row -= gap_extend;
///                 score_gap_col -= gap_extend;
///                 score_array[b_index].best_gap = MAX(score - gap_open_extend,
///                                                     score_gap_col);
///                 score_gap_row = MAX(score - gap_open_extend, score_gap_row);
///                 score_array[b_index].best = score;
///             }
///
///             score = next_score;
///         }
///
///         if (first_b_index == b_size)
///             break;
///         ...
///         if (last_b_index < b_size - 1) {
///             b_size = last_b_index + 1;
///         }
///         else {
///             while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {
///                 score_array[b_size].best = score_gap_row;
///                 score_array[b_size].best_gap = score_gap_row - gap_open_extend;
///                 score_gap_row -= gap_extend;
///                 b_size++;
///             }
///         }
///
///         if (b_size <= N) {
///             score_array[b_size].best = MININT;
///             score_array[b_size].best_gap = MININT;
///             b_size++;
///         }
///     }
///
///     return best_score;
/// ```
#[allow(clippy::too_many_arguments)]
fn semi_gapped_align_score_only<R: PssmRows>(
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
) -> (usize, usize, i32) {
    let gap_open_extend = gap_open + gap_extend;
    let x_dropoff = x_drop.max(gap_open_extend);
    if n == 0 || m == 0 {
        return (0, 0, 0);
    }
    scratch.reserve_dp(n);
    let dp = &mut scratch.dp;
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
        let matrix_row = pssm_row(rows, m, a_index, query_offset, reverse);
        let mut score = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;
        let mut b_index = first_b_index;
        while b_index < b_size {
            let letter = subject_letter(subject, s_base, n, b_index, reverse);
            let mut score_gap_col = dp[b_index].best_gap;
            let next_score = dp[b_index].best + matrix_row[letter as usize];
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
                dp[b_index].best_gap = (score - gap_open_extend).max(score_gap_col);
                score_gap_row = (score - gap_open_extend).max(score_gap_row);
                dp[b_index].best = score;
            }
            score = next_score;
            b_index += 1;
        }
        if first_b_index == b_size {
            break;
        }
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
/// ```c
///       q_length = init_hsp->offsets.qs_offsets.q_off + 1;
///       s_length = init_hsp->offsets.qs_offsets.s_off + 1;
///       query = query_blk->sequence;
///       subject = subject_blk->sequence;
///    ...
///    AdjustSubjectRange(&s_length, &subject_length, q_length, query_length,
///                       &subject_shift);
///
///    found_start = FALSE;
///    found_end = FALSE;
///
///    /* Looking for "left" score */
///    score_left = 0;
///    if (q_length != 0 && s_length != 0) {
///       found_start = TRUE;
///       ...
///             score_left = Blast_SemiGappedAlign(query, subject+subject_shift,
///                                    q_length, s_length,
///                                    &private_q_start, &private_s_start, TRUE,
///                                    NULL, gap_align, score_params,
///                                    init_hsp->offsets.qs_offsets.q_off,
///                                    FALSE, TRUE,
///                                    fence_hit);
///       gap_align->query_start = q_length - private_q_start;
///       gap_align->subject_start = s_length - private_s_start + subject_shift;
///    }
///
///    score_right = 0;
///    if (q_length < query_length && s_length < subject_length) {
///       found_end = TRUE;
///       ...
///             score_right = Blast_SemiGappedAlign(
///                                    query+init_hsp->offsets.qs_offsets.q_off,
///                                    subject+init_hsp->offsets.qs_offsets.s_off,
///                                    query_length-q_length,
///                                    subject_length-s_length,
///                                    &(gap_align->query_stop),
///                                    &(gap_align->subject_stop),
///                                    TRUE, NULL, gap_align, score_params,
///                                    init_hsp->offsets.qs_offsets.q_off, FALSE,
///                                    FALSE,
///                                    fence_hit);
///          ...
///          gap_align->query_stop += init_hsp->offsets.qs_offsets.q_off + 1;
///          gap_align->subject_stop += init_hsp->offsets.qs_offsets.s_off + 1;
///       }
///    }
///
///    if (found_start == FALSE) {  /* impossible for in-frame */
///       gap_align->query_start = init_hsp->offsets.qs_offsets.q_off;
///       gap_align->subject_start = init_hsp->offsets.qs_offsets.s_off;
///    }
///    if (found_end == FALSE) {    /* impossible for out-of-frame */
///       gap_align->query_stop = init_hsp->offsets.qs_offsets.q_off;
///       gap_align->subject_stop = init_hsp->offsets.qs_offsets.s_off;
///    }
///
///    gap_align->score = score_right+score_left;
/// ```
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
    let (private_q_start, private_s_start, score_left) = semi_gapped_align_score_only(
        rows,
        q_length,
        q_off,
        subject,
        subject_shift,
        s_length,
        gap_open,
        gap_extend,
        x_drop,
        true,
        scratch,
    );
    let query_start = (q_length - private_q_start) as i32;
    let subject_start = (s_length - private_s_start + subject_shift) as i32;
    let (query_stop, subject_stop, score_right) =
        if q_length < query_length && s_length < adjusted_length {
            let (q_stop, s_stop, score_right) = semi_gapped_align_score_only(
                rows,
                query_length - q_length,
                q_off,
                subject,
                s_off,
                adjusted_length - s_length,
                gap_open,
                gap_extend,
                x_drop,
                false,
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
/// ```c
///     if(N <= 0 || M <= 0)
///         return 0;
///     ...
///     for (i = 1; i <= N; i++) {
///         if (score < -x_dropoff)
///             break;
///
///         score_array[i].best = score;
///         score_array[i].best_gap = score - gap_open_extend;
///         score -= gap_extend;
///         edit_script_row[i] = SCRIPT_GAP_IN_A;
///     }
///     state_struct->used = i + 1;
///     ...
///     for (a_index = 1; a_index <= M; a_index++) {
///         ...
///         edit_script[a_index] = state_struct->state_array +
///                                 state_struct->used + 1;
///         edit_start_offset[a_index] = first_b_index;
///
///         /* the inner loop assumes the current traceback
///            row begins at offset zero of B */
///
///         edit_script_row = edit_script[a_index] - first_b_index;
///         orig_b_index = first_b_index;
///         ...
///         for (b_index = first_b_index; b_index < b_size; b_index++) {
///             int matrix_index = 0;
///
///             b_ptr += b_increment;
///             score_gap_col = score_array[b_index].best_gap;
///
///             matrix_index = *b_ptr;
///
///             if (matrix_index == FENCE_SENTRY) {
///                 if (fence_hit) {
///                     *fence_hit = 1;
///                 }
///                 break;
///             }
///
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
///
///             script = SCRIPT_SUB;
///             script_col = SCRIPT_EXTEND_GAP_B;
///             script_row = SCRIPT_EXTEND_GAP_A;
///
///             if (score < score_gap_col) {
///                 script = SCRIPT_GAP_IN_B;
///                 score = score_gap_col;
///             }
///             if (score < score_gap_row) {
///                 script = SCRIPT_GAP_IN_A;
///                 score = score_gap_row;
///             }
///
///             if (best_score - score > x_dropoff) {
///
///                 if (first_b_index == b_index)
///                     first_b_index++;
///                 else
///                     score_array[b_index].best = MININT;
///             }
///             else {
///                 last_b_index = b_index;
///                 if (score > best_score) {
///                     best_score = score;
///                     *a_offset = a_index;
///                     *b_offset = b_index;
///                 }
///
///                 score_gap_row -= gap_extend;
///                 score_gap_col -= gap_extend;
///                 if (score_gap_col < (score - gap_open_extend)) {
///                     score_array[b_index].best_gap = score - gap_open_extend;
///                 }
///                 else {
///                     score_array[b_index].best_gap = score_gap_col;
///                     script += script_col;
///                 }
///
///                 if (score_gap_row < (score - gap_open_extend))
///                     score_gap_row = score - gap_open_extend;
///                 else
///                     script += script_row;
///
///                 score_array[b_index].best = score;
///             }
///
///             score = next_score;
///             edit_script_row[b_index] = script;
///         }
///
///         if (first_b_index == b_size || (fence_hit && *fence_hit))
///             break;
///         ...
///         if (last_b_index < b_size - 1) {
///             b_size = last_b_index + 1;
///         }
///         else {
///             while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {
///
///                 score_array[b_size].best = score_gap_row;
///                 score_array[b_size].best_gap = score_gap_row - gap_open_extend;
///                 score_gap_row -= gap_extend;
///                 edit_script_row[b_size] = SCRIPT_GAP_IN_A;
///                 b_size++;
///             }
///         }
///         ...
///         if (b_size <= N) {
///             score_array[b_size].best = MININT;
///             score_array[b_size].best_gap = MININT;
///             b_size++;
///         }
///     }
///
///     a_index = *a_offset;
///     b_index = *b_offset;
///     script = SCRIPT_SUB;
///
///     if (fence_hit && *fence_hit)
///         goto done;
///
///     while (a_index > 0 || b_index > 0) {
///         next_script =
///             edit_script[a_index][b_index - edit_start_offset[a_index]];
///
///         switch(script) {
///         case SCRIPT_GAP_IN_A:
///             script = next_script & SCRIPT_OP_MASK;
///             if (next_script & SCRIPT_EXTEND_GAP_A)
///                 script = SCRIPT_GAP_IN_A;
///             break;
///
///         case SCRIPT_GAP_IN_B:
///             script = next_script & SCRIPT_OP_MASK;
///             if (next_script & SCRIPT_EXTEND_GAP_B)
///                 script = SCRIPT_GAP_IN_B;
///             break;
///
///         default:
///             script = next_script & SCRIPT_OP_MASK;
///             break;
///         }
///
///         if (script == SCRIPT_GAP_IN_A) {
///             b_index--;
///         }
///         else if (script == SCRIPT_GAP_IN_B) {
///             a_index--;
///         }
///         else {
///             a_index--;
///             b_index--;
///         }
///         GapPrelimEditBlockAdd(edit_block, (EGapAlignOpType)script, 1);
///     }
/// ```
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
    scratch.reserve_dp(n);
    let PssmGapScratch {
        dp,
        cells,
        rows: trace_rows,
        ops,
    } = scratch;
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
        if first_b_index == b_size || *fence_hit {
            break;
        }
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

/// `Blast_PrelimEditBlockToGapEditScript`: the left extension's operations as they are, then the
/// right extension's reversed, merging the two last operations when they are of one kind.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2494-2534
/// ```c
///    if (fwd_prelim_tback->num_ops > 0 && rev_prelim_tback->num_ops > 0 &&
///        fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].op_type ==
///          rev_prelim_tback->edit_ops[(rev_prelim_tback->num_ops)-1].op_type)
///      merge_ops = TRUE;
///    ...
///    for (i=0; i < rev_prelim_tback->num_ops; i++) {
///       ...
///    }
///
///    if (fwd_prelim_tback->num_ops == 0)
///       return esp;
///
///    if (merge_ops)
///        esp->num[index-1] += fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].num;
///
///    /* If we merge, then we skip the first one. */
///    if (merge_ops)
///       i = fwd_prelim_tback->num_ops - 2;
///    else
///       i = fwd_prelim_tback->num_ops - 1;
///
///    for (; i >= 0; i--) {
/// ```
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
}

/// `BLAST_GappedAlignmentWithTraceback` with a PSSM (not out-of-frame): `ALIGN_EX` to the left
/// (including the start) and to the right of `(q_start, s_start)`, the edit script, and the
/// pruning of leading and trailing gaps. `subject` is the adjusted subject (its letters from the
/// traceback's start shift; `subject_length` of them, then the bytes after them), `q_start` an
/// absolute query position.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4599-4718
/// ```c
///         score_left = ALIGN_EX(query, subject, q_start+1, s_start+1,  &private_q_length, &private_s_length, rev_prelim_tback,
///                               gap_align,
///                               score_params, q_start, FALSE /*reversed*/, TRUE /*reverse_sequence*/,
///              fence_hit);
///
///        gap_align->query_start = q_start - private_q_length + 1;
///        gap_align->subject_start = s_start - private_s_length + 1;
///     }
///
///     score_right = 0;
///
///     if ((! (fence_hit && *fence_hit)) &&
///         (q_start < q_length) &&
///         (s_start < s_length)) {
///
///        found_end = TRUE;
///        ...
///             score_right =
///                ALIGN_EX(query+q_start, subject+s_start,
///                   q_length-q_start-1, s_length-s_start-1, &private_q_length,
///                   &private_s_length, fwd_prelim_tback, gap_align,
///                   score_params, q_start, FALSE, FALSE,
///                   fence_hit);
///         }
///
///         gap_align->query_stop = q_start + private_q_length + 1;
///         gap_align->subject_stop = s_start + private_s_length + 1;
///     }
///
///     if (found_end == FALSE) {
///         gap_align->query_stop = q_start - 1;
///         gap_align->subject_stop = s_start - 1;
///     }
///     ...
///         GapEditScript *esp = Blast_PrelimEditBlockToGapEditScript(
///                                                         rev_prelim_tback,
///                                                         fwd_prelim_tback);
///         ...
///         gap_align->edit_script = esp;
///         if (esp) {
///             while (esp->size && esp->op_type[0] != eGapAlignSub) {
///                 score_left += score_params->gap_open +
///                              esp->num[0] * score_params->gap_extend;
///
///                 if (esp->op_type[0] == eGapAlignDel)
///                     gap_align->subject_start += esp->num[0];
///                 else
///                     gap_align->query_start += esp->num[0];
///                 ...
///             }
///             i = esp->size;
///             while (esp->size && esp->op_type[i-1] != eGapAlignSub) {
///                 score_right += score_params->gap_open +
///                              esp->num[i-1] * score_params->gap_extend;
///
///                 if (esp->op_type[i-1] == eGapAlignDel)
///                     gap_align->subject_stop -= esp->num[i-1];
///                 else
///                     gap_align->query_stop -= esp->num[i-1];
///                 ...
///             }
///         }
///     }
///
///     gap_align->score = score_right + score_left;
/// ```
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
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::pssm_search::rows::QueryPssm;
    use crate::core::blast_util::BLASTAA_SIZE;

    /// A PSSM built from BLOSUM62's rows of `query` (NCBIstdaa): the PSSM kernels must then give
    /// the alignments of blastp's matrix kernels.
    fn blosum62_pssm(query: &[u8]) -> QueryPssm {
        QueryPssm {
            rows: query
                .iter()
                .map(|&q| {
                    let mut row = [0i32; BLASTAA_SIZE];
                    for (s, cell) in row.iter_mut().enumerate() {
                        *cell = crate::core::blast_util::blosum62_score(q, s as u8);
                    }
                    row
                })
                .collect(),
        }
    }

    fn encode(text: &str) -> Vec<u8> {
        text.bytes()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect()
    }

    // The PSSM kernels on BLOSUM62 rows agree with blastp's matrix kernels (gapalign.rs), which
    // the BLASTP fixtures check against NCBI: same ends, score and edit script, on a pair with a
    // gap in each sequence.
    #[test]
    fn blosum62_rows_give_blastps_alignments() {
        let query = encode("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ");
        let subject_letters =
            encode("MKTAYIAKQRQISFVKSHFSRQDLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVKVKVKALPDAQ");
        let mut subject = subject_letters.clone();
        subject.push(0);
        let rows = blosum62_pssm(&query);
        let mut scratch = PssmGapScratch::new();
        let (q_off, s_off) = (30usize, 31usize);
        let ours = pssm_score_only_gapped_alignment(
            &rows,
            &subject,
            subject_letters.len(),
            q_off,
            s_off,
            11,
            1,
            38,
            &mut scratch,
        );
        let theirs = crate::algorithm::blastp::gapalign::blastp_score_only_gapped_alignment(
            &query,
            &subject_letters,
            q_off,
            s_off,
            crate::config::ScoringMatrix::Blosum62,
            11,
            1,
            38,
            crate::algorithm::blastp::gapalign::BlastpGappedAlignmentMode::Exact,
        );
        assert_eq!(
            (
                ours.query_start,
                ours.query_stop,
                ours.subject_start,
                ours.subject_stop,
                ours.score
            ),
            (
                theirs.query_start,
                theirs.query_stop,
                theirs.subject_start,
                theirs.subject_stop,
                theirs.score
            )
        );
        let traced = pssm_gapped_alignment_with_traceback(
            &rows,
            &subject,
            subject_letters.len(),
            q_off,
            s_off,
            11,
            1,
            64,
            &mut scratch,
        );
        let reference = crate::algorithm::blastp::gapalign::blast_gapped_alignment_with_traceback(
            &query,
            &subject_letters,
            q_off,
            s_off,
            crate::config::ScoringMatrix::Blosum62,
            None,
            11,
            1,
            64,
        )
        .expect("an alignment");
        assert_eq!(traced.query_start, reference.query_start);
        assert_eq!(traced.query_stop, reference.query_stop);
        assert_eq!(traced.subject_start, reference.subject_start);
        assert_eq!(traced.subject_stop, reference.subject_stop);
        assert_eq!(traced.score, reference.score);
        assert_eq!(traced.edit_script, reference.edit_script);
        assert!(traced
            .edit_script
            .iter()
            .any(|op| matches!(op, GapEditOp::Del(_))));
    }

    #[test]
    fn gapped_start_is_the_middle_of_the_best_window() {
        let query = encode("AAAAAAAAAAWWWWWWWWWWWWAAAA");
        let rows = blosum62_pssm(&query);
        // The subject equals the query: the windows grow to the first all-W one (ending at 20).
        let start = pssm_get_start_for_gapped_alignment(&rows, &query, 0, 26, 0);
        assert_eq!(start, 20 - HSP_MAX_WINDOW / 2);
        // A short HSP starts at its middle.
        assert_eq!(
            pssm_get_start_for_gapped_alignment(&rows, &query, 4, 10, 4),
            9
        );
    }
}
