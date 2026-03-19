//! BLASTP gapped alignment with traceback.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c

use crate::common::{GapEditOp, Hit};
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::AdjustedProteinMatrix;
use crate::stats::karlin::{bit_score as calc_bit_score, evalue_from_raw_score};
use crate::stats::KarlinParams;
use crate::utils::matrix::protein_score;

use super::alignment::build_alignment_view_with_matrix;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:363-371
// ```c
// enum {
//     SCRIPT_SUB           = eGapAlignSub,
//     SCRIPT_GAP_IN_A      = eGapAlignDel,
//     SCRIPT_GAP_IN_B      = eGapAlignIns,
//     SCRIPT_OP_MASK       = 0x07,
//     SCRIPT_EXTEND_GAP_A  = 0x10,
//     SCRIPT_EXTEND_GAP_B  = 0x40
// };
// ```
const SCRIPT_SUB: u8 = 3;
const SCRIPT_GAP_IN_A: u8 = 0;
const SCRIPT_GAP_IN_B: u8 = 6;
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:59-60
// ```c
// /** Lower bound for scores. Divide by two to prevent underflows. */
// #define MININT INT4_MIN/2
// ```
const GAP_MININT: i32 = i32::MIN / 2;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:339-360
// ```c
// #define MAX_SUBJECT_OFFSET 90000
// #define MAX_TOTAL_GAPS 3000
// ```
const MAX_SUBJECT_OFFSET: usize = 90_000;
const MAX_TOTAL_GAPS: usize = 3_000;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:963-965
// ```c
// /** For restricted gapped alignment, gaps may only start once in
//     this many sequence offsets */
// #define RESTRICT_SIZE 10
// ```
const RESTRICT_SIZE: usize = 10;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:399-414
// ```c
// matrix = gap_align->sbp->matrix->data;
// if (gap_align->positionBased) {
//     pssm = gap_align->sbp->psi_matrix->pssm->data;
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1748-1770
// ```c
// s_NewAlignmentFromGapAlign(...);
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1902-1945
// ```c
// BLAST_GappedAlignmentWithTraceback(...);
// ```
#[derive(Clone, Copy)]
enum BlastpScoreMatrix<'a> {
    Standard(ScoringMatrix),
    Adjusted(&'a AdjustedProteinMatrix),
}

impl BlastpScoreMatrix<'_> {
    #[inline]
    fn score(self, query_residue: u8, subject_residue: u8) -> i32 {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:830-831
        // ```c
        // next_score = score_array[b_index].best + matrix_row[*b_ptr];
        // ```
        //
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3410-3426
        // ```c
        // score += sbp->matrix->data[*query_var][*subject_var];
        // ...
        // score += sbp->matrix->data[*query_var][*subject_var];
        // ```
        match self {
            Self::Standard(matrix) => protein_score(matrix, query_residue, subject_residue),
            Self::Adjusted(matrix) => matrix.score(query_residue, subject_residue),
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:523-550
// ```c
// score = 0;
// sum = 0;
// ...
// if (esp->op_type[index] == eGapAlignSub) {
//     sum += factor*matrix[*query & kResidueMask][*subject];
//     query++;
//     subject++;
//     op_index++;
// } else if (esp->op_type[index] == eGapAlignDel) {
//     sum -= gap_open + gap_extend * esp->num[index];
//     subject += esp->num[index];
//     op_index += esp->num[index];
// } else if (esp->op_type[index] == eGapAlignIns) {
//     sum -= gap_open + gap_extend * esp->num[index];
//     query += esp->num[index];
//     op_index += esp->num[index];
// }
// if (sum < 0) {
//     sum = 0;
// } else if (sum > score) {
//     score = sum;
// }
// ```
#[cfg(test)]
fn score_edit_script_local_alignment(
    query: &[u8],
    subject: &[u8],
    query_start: usize,
    subject_start: usize,
    edit_script: &[GapEditOp],
    score_matrix: BlastpScoreMatrix<'_>,
    gap_open: i32,
    gap_extend: i32,
) -> i32 {
    let mut query_pos = query_start;
    let mut subject_pos = subject_start;
    let mut score = 0i32;
    let mut sum = 0i32;

    for op in edit_script {
        match *op {
            GapEditOp::Sub(n) => {
                for _ in 0..n as usize {
                    if query_pos >= query.len() || subject_pos >= subject.len() {
                        return score;
                    }
                    sum += score_matrix.score(query[query_pos], subject[subject_pos]);
                    query_pos += 1;
                    subject_pos += 1;
                    if sum < 0 {
                        sum = 0;
                    } else if sum > score {
                        score = sum;
                    }
                }
            }
            GapEditOp::Del(n) => {
                let n = n as usize;
                sum -= gap_open + gap_extend * n as i32;
                subject_pos = subject_pos.saturating_add(n);
                if sum < 0 {
                    sum = 0;
                } else if sum > score {
                    score = sum;
                }
            }
            GapEditOp::Ins(n) => {
                let n = n as usize;
                sum -= gap_open + gap_extend * n as i32;
                query_pos = query_pos.saturating_add(n);
                if sum < 0 {
                    sum = 0;
                } else if sum > score {
                    score = sum;
                }
            }
        }
    }

    score
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:56-62
// ```c
// typedef struct BlastGapDP {
//    Int4 best;
//    Int4 best_gap;
// } BlastGapDP;
// ```
#[derive(Clone, Copy)]
struct BlastGapDp {
    best: i32,
    best_gap: i32,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80
// ```c
// typedef struct BlastGapAlignStruct {
//    BlastGapDP* dp_mem;
//    Int4 dp_mem_alloc;
//    Uint1** edit_script;
//    Int4* edit_start_offset;
// } BlastGapAlignStruct;
// ```
struct GapAlignScratch {
    dp_mem: Vec<BlastGapDp>,
    dp_mem_alloc: usize,
    trace_rows: Vec<Vec<u8>>,
    trace_offsets: Vec<usize>,
    trace_rows_used: usize,
}

impl GapAlignScratch {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:313-319
    // ```c
    // gap_align->dp_mem_alloc = 1000;
    // gap_align->dp_mem = (BlastGapDP *) calloc(gap_align->dp_mem_alloc,
    //                                           sizeof(BlastGapDP));
    // ```
    fn new() -> Self {
        let dp_mem_alloc = 1000usize;
        Self {
            dp_mem: vec![
                BlastGapDp {
                    best: GAP_MININT,
                    best_gap: GAP_MININT,
                };
                dp_mem_alloc
            ],
            dp_mem_alloc,
            trace_rows: Vec::with_capacity(100),
            trace_offsets: Vec::with_capacity(100),
            trace_rows_used: 0,
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4077-4081
// ```c
// Blast_HSPInit(gap_align->query_start, gap_align->query_stop,
//               gap_align->subject_start, gap_align->subject_stop,
//               ..., gap_align->score, &(gap_align->edit_script), &new_hsp);
// ```
#[derive(Debug, Clone)]
pub(crate) struct BlastpGapAlignResult {
    pub query_start: usize,
    pub query_stop: usize,
    pub subject_start: usize,
    pub subject_stop: usize,
    pub score: i32,
    pub edit_script: Vec<GapEditOp>,
    pub num_ident: usize,
    pub num_positives: usize,
    pub mismatches: usize,
    pub gap_opens: usize,
    pub gap_letters: usize,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4077-4091
// ```c
// Blast_HSPInit(gap_align->query_start, gap_align->query_stop,
//               gap_align->subject_start, gap_align->subject_stop,
//               init_hsp->offsets.qs_offsets.q_off,
//               init_hsp->offsets.qs_offsets.s_off, context,
//               query_frame, subject->frame, gap_align->score,
//               &(gap_align->edit_script), &new_hsp);
// ```
#[derive(Debug, Clone, Copy)]
pub(crate) struct BlastpPreliminaryHsp {
    pub query_start: i32,
    pub query_end: i32,
    pub subject_start: i32,
    pub subject_end: i32,
    pub gapped_query_start: i32,
    pub gapped_subject_start: i32,
    pub raw_score: i32,
    pub query_context: i32,
    pub query_frame: i32,
    pub subject_frame: i32,
    pub query_length: usize,
    pub q_idx: u32,
    pub s_idx: u32,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319
// ```c
// static Int2
// s_BlastProtGappedAlignment(..., Boolean restricted_alignment, ...)
// {
//    ...
//    score_left = Blast_SemiGappedAlign(...) or s_RestrictedGappedAlign(...);
//    ...
//    score_right = Blast_SemiGappedAlign(...) or s_RestrictedGappedAlign(...);
//    ...
//    gap_align->score = score_right+score_left;
// }
// ```
#[derive(Debug, Clone, Copy)]
pub(crate) struct BlastpScoreOnlyGappedAlignment {
    pub query_start: i32,
    pub query_stop: i32,
    pub subject_start: i32,
    pub subject_stop: i32,
    pub score: i32,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319
// ```c
// s_BlastProtGappedAlignment(..., Boolean restricted_alignment, ...)
// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BlastpGappedAlignmentMode {
    Exact,
    Restricted,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:797-808
// ```c
// if (num_extra_cells > gap_align->dp_mem_alloc) {
//     gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
//                                   2 * gap_align->dp_mem_alloc);
//     ...
// }
// ```
fn gap_dp_reserve_initial(
    dp_mem: &mut Vec<BlastGapDp>,
    dp_mem_alloc: &mut usize,
    num_extra_cells: usize,
) {
    if num_extra_cells > *dp_mem_alloc {
        *dp_mem_alloc = (num_extra_cells + 100).max(2 * *dp_mem_alloc);
        dp_mem.resize(
            *dp_mem_alloc,
            BlastGapDp {
                best: GAP_MININT,
                best_gap: GAP_MININT,
            },
        );
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:923-931
// ```c
// if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {
//     gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
//                                   2 * gap_align->dp_mem_alloc);
//     ...
// }
// ```
fn gap_dp_reserve_band(
    dp_mem: &mut Vec<BlastGapDp>,
    dp_mem_alloc: &mut usize,
    last_b_index: usize,
    num_extra_cells: usize,
) {
    if last_b_index + num_extra_cells + 3 >= *dp_mem_alloc {
        *dp_mem_alloc = (last_b_index + num_extra_cells + 100).max(2 * *dp_mem_alloc);
        dp_mem.resize(
            *dp_mem_alloc,
            BlastGapDp {
                best: GAP_MININT,
                best_gap: GAP_MININT,
            },
        );
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:128-139
// ```c
// static void s_GapPurgeState(BlastGapAlignStruct* gap_align)
// {
//     gap_align->state_struct->used = 0;
// }
// ```
fn gap_reset_traceback_state(scratch: &mut GapAlignScratch) {
    scratch.trace_rows_used = 0;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:70-115
// ```c
// static GapStateArrayStruct*
// s_GapGetState(...)
// {
//     Int4 chunksize = MAX(CHUNKSIZE, length + length/3);
//     ...
// }
// ```
fn gap_alloc_trace_row<'a>(
    trace_rows: &'a mut Vec<Vec<u8>>,
    trace_offsets: &mut Vec<usize>,
    trace_rows_used: &mut usize,
    row_capacity: usize,
    start_offset: usize,
) -> &'a mut Vec<u8> {
    let row_index = *trace_rows_used;
    *trace_rows_used += 1;

    let row_capacity_slack = row_capacity.saturating_add(row_capacity / 3);
    if row_index < trace_offsets.len() {
        trace_offsets[row_index] = start_offset;
    } else {
        trace_offsets.push(start_offset);
    }

    if row_index < trace_rows.len() {
        let row = &mut trace_rows[row_index];
        row.clear();
        if row.capacity() < row_capacity_slack {
            row.reserve(row_capacity_slack.saturating_sub(row.capacity()));
        }
        if row_capacity > 0 {
            row.resize(row_capacity, 0);
        }
        row
    } else {
        let mut row = Vec::with_capacity(row_capacity_slack);
        if row_capacity > 0 {
            row.resize(row_capacity, 0);
        }
        trace_rows.push(row);
        trace_rows.last_mut().unwrap()
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3410-3438
// ```c
// if (q_length <= HSP_MAX_WINDOW) {
//     max_offset = q_start + q_length/2;
//     return max_offset;
// }
// ...
// ```
pub(crate) fn blastp_get_start_for_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    q_length: usize,
    s_start: usize,
    s_length: usize,
    matrix: ScoringMatrix,
) -> usize {
    const HSP_MAX_WINDOW: usize = 11;
    if q_length <= HSP_MAX_WINDOW {
        return q_start + q_length / 2;
    }

    let mut score = 0i32;
    for offset in 0..HSP_MAX_WINDOW {
        score += protein_score(matrix, query[q_start + offset], subject[s_start + offset]);
    }
    let mut max_score = score;
    let mut max_offset = q_start + HSP_MAX_WINDOW - 1;
    let hsp_end = q_start + q_length.min(s_length);

    for index in (q_start + HSP_MAX_WINDOW)..hsp_end {
        let subject_index = s_start + (index - q_start);
        score -= protein_score(
            matrix,
            query[index - HSP_MAX_WINDOW],
            subject[subject_index - HSP_MAX_WINDOW],
        );
        score += protein_score(matrix, query[index], subject[subject_index]);
        if score > max_score {
            max_score = score;
            max_offset = index;
        }
    }

    if max_score > 0 {
        max_offset - HSP_MAX_WINDOW / 2
    } else {
        q_start
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:340-377
// ```c
// void AdjustSubjectRange(Int4* subject_offset_ptr, Int4* subject_length_ptr,
//                        Int4 query_offset, Int4 query_length, Int4* start_shift)
// {
//     ...
// }
// ```
pub(crate) fn adjust_subject_range(
    subject_offset: usize,
    subject_length: usize,
    query_offset: usize,
    query_length: usize,
) -> (usize, usize, usize) {
    if subject_length < MAX_SUBJECT_OFFSET {
        return (subject_offset, subject_length, 0);
    }

    let max_extension_left = query_offset + MAX_TOTAL_GAPS;
    let max_extension_right = query_length - query_offset + MAX_TOTAL_GAPS;

    let (adjusted_offset, start_shift) = if subject_offset <= max_extension_left {
        (subject_offset, 0)
    } else {
        (max_extension_left, subject_offset - max_extension_left)
    };

    let adjusted_length = subject_length.min(subject_offset + max_extension_right) - start_shift;
    (adjusted_offset, adjusted_length, start_shift)
}

fn push_trace_op(ops_reversed: &mut Vec<(u8, u32)>, op_type: u8) {
    if let Some(last) = ops_reversed.last_mut() {
        if last.0 == op_type {
            last.1 += 1;
            return;
        }
    }
    ops_reversed.push((op_type, 1));
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2481-2534
// ```c
// Blast_PrelimEditBlockToGapEditScript (GapPrelimEditBlock* rev_prelim_tback,
//                                       GapPrelimEditBlock* fwd_prelim_tback)
// {
//    ...
//    if (fwd_prelim_tback->num_ops > 0 && rev_prelim_tback->num_ops > 0 &&
//        fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].op_type ==
//          rev_prelim_tback->edit_ops[(rev_prelim_tback->num_ops)-1].op_type)
//      merge_ops = TRUE;
//    ...
// }
// ```
fn prelim_edit_ops_to_gap_edit_script(
    rev_prelim_tback: Vec<GapEditOp>,
    fwd_prelim_tback: Vec<GapEditOp>,
) -> Vec<GapEditOp> {
    let mut esp = rev_prelim_tback;
    let merge_ops = matches!(
        (esp.last(), fwd_prelim_tback.last()),
        (Some(GapEditOp::Sub(_)), Some(GapEditOp::Sub(_)))
            | (Some(GapEditOp::Del(_)), Some(GapEditOp::Del(_)))
            | (Some(GapEditOp::Ins(_)), Some(GapEditOp::Ins(_)))
    );

    if merge_ops {
        if let (Some(last), Some(first_fwd)) = (esp.last_mut(), fwd_prelim_tback.last()) {
            match (last, first_fwd) {
                (GapEditOp::Sub(lhs), GapEditOp::Sub(rhs))
                | (GapEditOp::Del(lhs), GapEditOp::Del(rhs))
                | (GapEditOp::Ins(lhs), GapEditOp::Ins(rhs)) => {
                    *lhs += *rhs;
                }
                _ => unreachable!("merge_ops only allows equal operation kinds"),
            }
        }
    }

    let skip = usize::from(merge_ops);
    for op in fwd_prelim_tback.into_iter().rev().skip(skip) {
        esp.push(op);
    }
    esp
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4681-4707
// ```c
// while (esp->size && esp->op_type[0] != eGapAlignSub) {
//     score_left += score_params->gap_open + esp->num[0] * score_params->gap_extend;
//     ...
// }
// ...
// while (esp->size && esp->op_type[i-1] != eGapAlignSub) {
//     score_right += score_params->gap_open + esp->num[i-1] * score_params->gap_extend;
//     ...
// }
// ```
fn prune_terminal_gap_ops(
    edit_script: &mut Vec<GapEditOp>,
    gap_open: i32,
    gap_extend: i32,
    score_left: &mut i32,
    score_right: &mut i32,
    query_start: &mut usize,
    query_stop: &mut usize,
    subject_start: &mut usize,
    subject_stop: &mut usize,
) {
    while let Some(first) = edit_script.first().copied() {
        match first {
            GapEditOp::Sub(_) => break,
            GapEditOp::Del(n) => {
                *score_left += gap_open + gap_extend * n as i32;
                *subject_start += n as usize;
            }
            GapEditOp::Ins(n) => {
                *score_left += gap_open + gap_extend * n as i32;
                *query_start += n as usize;
            }
        }
        edit_script.remove(0);
    }

    while let Some(last) = edit_script.last().copied() {
        match last {
            GapEditOp::Sub(_) => break,
            GapEditOp::Del(n) => {
                *score_right += gap_open + gap_extend * n as i32;
                *subject_stop = subject_stop.saturating_sub(n as usize);
            }
            GapEditOp::Ins(n) => {
                *score_right += gap_open + gap_extend * n as i32;
                *query_stop = query_stop.saturating_sub(n as usize);
            }
        }
        edit_script.pop();
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:745-818
// ```c
// if (*q == *s)
//    num_ident++;
// else if (matrix[*q][*s] > 0)
//    num_pos ++;
// ```
fn stats_from_edit_ops_protein(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    edit_ops: &[GapEditOp],
    matrix: ScoringMatrix,
) -> (usize, usize, usize, usize, usize) {
    let mut num_ident = 0usize;
    let mut num_positives = 0usize;
    let mut mismatches = 0usize;
    let mut gap_opens = 0usize;
    let mut gap_letters = 0usize;
    let mut qi = q_start;
    let mut si = s_start;

    for op in edit_ops {
        let num = op.num() as usize;
        match *op {
            GapEditOp::Sub(_) => {
                for _ in 0..num {
                    let q = q_seq[qi];
                    let s = s_seq[si];
                    if q == s {
                        num_ident += 1;
                        num_positives += 1;
                    } else {
                        mismatches += 1;
                        if protein_score(matrix, q, s) > 0 {
                            num_positives += 1;
                        }
                    }
                    qi += 1;
                    si += 1;
                }
            }
            GapEditOp::Del(_) => {
                gap_opens += 1;
                gap_letters += num;
                si += num;
            }
            GapEditOp::Ins(_) => {
                gap_opens += 1;
                gap_letters += num;
                qi += num;
            }
        }
    }

    (num_ident, num_positives, mismatches, gap_opens, gap_letters)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:839-849
// ```c
// if(reverse_sequence)
//     matrix_row = matrix[ A[ M - a_index ] ];
// else
//     matrix_row = matrix[ A[ a_index ] ];
// ```
#[inline]
fn blastp_query_residue(
    query: &[u8],
    query_base: usize,
    query_length: usize,
    a_index: usize,
    reverse_sequence: bool,
) -> u8 {
    let absolute_index = if reverse_sequence {
        query_base + query_length - a_index
    } else {
        query_base + a_index
    };
    query[absolute_index]
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:852-866
// ```c
// if(reverse_sequence)
//     b_ptr = &B[N - first_b_index];
// else
//     b_ptr = &B[first_b_index];
// ...
// b_ptr += b_increment;
// next_score = score_array[b_index].best + matrix_row[*b_ptr];
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:155-161
// ```c
// seq_blk->sequence_start = (Uint1*) sequence;
// seq_blk->sequence = (Uint1*) sequence + 1;
// seq_blk->length = seqlen;
// ```
#[inline]
fn blastp_subject_residue(
    subject: &[u8],
    subject_base: usize,
    subject_length: usize,
    b_index: usize,
    reverse_sequence: bool,
) -> u8 {
    let absolute_index = if reverse_sequence {
        subject_base as isize + subject_length as isize - b_index as isize - 1
    } else {
        subject_base as isize + b_index as isize + 1
    };
    if absolute_index < 0 {
        0
    } else {
        subject.get(absolute_index as usize).copied().unwrap_or(0)
    }
}

fn align_ex_protein_score_only(
    query: &[u8],
    query_base: usize,
    subject: &[u8],
    subject_base: usize,
    len1: usize,
    len2: usize,
    score_matrix: BlastpScoreMatrix<'_>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    scratch: &mut GapAlignScratch,
) -> (usize, usize, i32) {
    if len1 == 0 || len2 == 0 {
        return (0, 0, 0);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:737-957
    // ```c
    // Int4 Blast_SemiGappedAlign(..., Boolean score_only, ...) {
    //    if (!score_only) {
    //        return ALIGN_EX(...);
    //    }
    //    ...
    //    for (a_index = 1; a_index <= M; a_index++) {
    //        ...
    //        for (b_index = first_b_index; b_index < b_size; b_index++) {
    //            ...
    //            if (score < score_gap_col) score = score_gap_col;
    //            if (score < score_gap_row) score = score_gap_row;
    //            ...
    //            score_array[b_index].best_gap =
    //                MAX(score - gap_open_extend, score_gap_col - gap_extend);
    //            score_gap_row = MAX(score - gap_open_extend, score_gap_row - gap_extend);
    //        }
    //    }
    // }
    // ```
    let gap_open_extend = gap_open + gap_extend;
    let mut x_dropoff = x_drop;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend + 3) as usize
    } else {
        len2 + 3
    };
    gap_dp_reserve_initial(
        &mut scratch.dp_mem,
        &mut scratch.dp_mem_alloc,
        num_extra_cells,
    );

    let mut score = -gap_open_extend;
    scratch.dp_mem[0].best = 0;
    scratch.dp_mem[0].best_gap = -gap_open_extend;

    let mut b_size = 1usize;
    for index in 1..=len2.min(scratch.dp_mem_alloc.saturating_sub(1)) {
        if score < -x_dropoff {
            break;
        }
        scratch.dp_mem[index].best = score;
        scratch.dp_mem[index].best_gap = score - gap_open_extend;
        score -= gap_extend;
        b_size = index + 1;
    }

    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;

    for a_index in 1..=len1 {
        let qc = blastp_query_residue(query, query_base, len1, a_index, reverse);
        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;

        for b_index in first_b_index..b_size {
            let sc = blastp_subject_residue(subject, subject_base, len2, b_index, reverse);
            let score_gap_col = scratch.dp_mem[b_index].best_gap;
            let next_score = scratch.dp_mem[b_index].best + score_matrix.score(qc, sc);

            if score_val < score_gap_col {
                score_val = score_gap_col;
            }
            if score_val < score_gap_row {
                score_val = score_gap_row;
            }

            if best_score - score_val > x_dropoff {
                if b_index == first_b_index {
                    first_b_index += 1;
                } else {
                    scratch.dp_mem[b_index].best = GAP_MININT;
                }
            } else {
                last_b_index = b_index;
                if score_val > best_score {
                    best_score = score_val;
                    a_offset = a_index;
                    b_offset = b_index;
                }

                score_gap_row -= gap_extend;
                let score_gap_col_ext = score_gap_col - gap_extend;
                scratch.dp_mem[b_index].best_gap =
                    (score_val - gap_open_extend).max(score_gap_col_ext);
                score_gap_row = (score_val - gap_open_extend).max(score_gap_row);
                scratch.dp_mem[b_index].best = score_val;
            }

            score_val = next_score;
        }

        if first_b_index == b_size {
            break;
        }

        gap_dp_reserve_band(
            &mut scratch.dp_mem,
            &mut scratch.dp_mem_alloc,
            last_b_index,
            num_extra_cells,
        );

        if last_b_index < b_size.saturating_sub(1) {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - x_dropoff
                && b_size <= len2
                && b_size < scratch.dp_mem_alloc
            {
                scratch.dp_mem[b_size].best = score_gap_row;
                scratch.dp_mem[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size += 1;
            }
        }

        if b_size <= len2 && b_size < scratch.dp_mem_alloc {
            scratch.dp_mem[b_size].best = GAP_MININT;
            scratch.dp_mem[b_size].best_gap = GAP_MININT;
            b_size += 1;
        }
    }

    (a_offset, b_offset, best_score)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:991-1285
// ```c
// static Int4
// s_RestrictedGappedAlign(const Uint1* A, const Uint1* B, Int4 M, Int4 N,
//    Int4* a_offset, Int4* b_offset, BlastGapAlignStruct* gap_align,
//    const BlastScoringParameters* score_params, Int4 query_offset,
//    Boolean reverse_sequence)
// {
//     ...
// }
// ```
fn restricted_gapped_align_protein(
    query: &[u8],
    query_base: usize,
    subject: &[u8],
    subject_base: usize,
    len1: usize,
    len2: usize,
    score_matrix: BlastpScoreMatrix<'_>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    dp_mem: &mut Vec<BlastGapDp>,
    dp_mem_alloc: &mut usize,
) -> (usize, usize, i32) {
    if len1 == 0 || len2 == 0 {
        return (0, 0, 0);
    }

    let gap_open_extend = gap_open + gap_extend;
    let mut x_dropoff = x_drop;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend + 3) as usize
    } else {
        len2 + 3
    };
    gap_dp_reserve_initial(dp_mem, dp_mem_alloc, num_extra_cells);

    let mut score = -gap_open_extend;
    dp_mem[0].best = 0;
    dp_mem[0].best_gap = -gap_open_extend;

    let mut b_size = 1usize;
    for index in 1..=len2.min(dp_mem.len().saturating_sub(1)) {
        if score < -x_dropoff {
            break;
        }
        dp_mem[index].best = score;
        dp_mem[index].best_gap = score - gap_open_extend;
        score -= gap_extend;
        b_size = index + 1;
    }

    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;
    let mut b_gap = 0usize;

    for a_index in 1..=len1 {
        let matrix_row_char = blastp_query_residue(query, query_base, len1, a_index, reverse);
        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;

        if a_index % RESTRICT_SIZE != 0 {
            for b_index in first_b_index..b_size {
                let sc = blastp_subject_residue(subject, subject_base, len2, b_index, reverse);
                let next_score = dp_mem[b_index].best + score_matrix.score(matrix_row_char, sc);

                if b_index != b_gap {
                    if best_score - score_val > x_dropoff {
                        dp_mem[b_index].best = GAP_MININT;
                        if b_index == first_b_index {
                            first_b_index += 1;
                        }
                    } else {
                        last_b_index = b_index;
                        if score_val > best_score {
                            best_score = score_val;
                            a_offset = a_index;
                            b_offset = b_index;
                        }
                        dp_mem[b_index].best = score_val;
                    }
                } else {
                    b_gap += RESTRICT_SIZE;
                    let mut score_gap_col = dp_mem[b_index].best_gap;
                    if score_val < score_gap_col {
                        score_val = score_gap_col;
                    }

                    if best_score - score_val > x_dropoff {
                        dp_mem[b_index].best = GAP_MININT;
                        if b_index == first_b_index {
                            first_b_index += 1;
                        }
                    } else {
                        last_b_index = b_index;
                        if score_val > best_score {
                            best_score = score_val;
                            a_offset = a_index;
                            b_offset = b_index;
                        }

                        score_gap_col -= gap_extend;
                        dp_mem[b_index].best_gap = (score_val - gap_open_extend).max(score_gap_col);
                        dp_mem[b_index].best = score_val;
                    }
                }
                score_val = next_score;
            }
            score_gap_row = score_val;
        } else {
            for b_index in first_b_index..b_size {
                let sc = blastp_subject_residue(subject, subject_base, len2, b_index, reverse);
                let next_score = dp_mem[b_index].best + score_matrix.score(matrix_row_char, sc);

                if b_index != b_gap {
                    if score_val < score_gap_row {
                        score_val = score_gap_row;
                    }

                    if best_score - score_val > x_dropoff {
                        dp_mem[b_index].best = GAP_MININT;
                        if b_index == first_b_index {
                            first_b_index += 1;
                        }
                    } else {
                        last_b_index = b_index;
                        if score_val > best_score {
                            best_score = score_val;
                            a_offset = a_index;
                            b_offset = b_index;
                        }

                        score_gap_row -= gap_extend;
                        score_gap_row = (score_val - gap_open_extend).max(score_gap_row);
                        dp_mem[b_index].best = score_val;
                    }
                } else {
                    b_gap += RESTRICT_SIZE;
                    let mut score_gap_col = dp_mem[b_index].best_gap;

                    if score_val < score_gap_col {
                        score_val = score_gap_col;
                    }
                    if score_val < score_gap_row {
                        score_val = score_gap_row;
                    }

                    if best_score - score_val > x_dropoff {
                        dp_mem[b_index].best = GAP_MININT;
                        if b_index == first_b_index {
                            first_b_index += 1;
                        }
                    } else {
                        last_b_index = b_index;
                        if score_val > best_score {
                            best_score = score_val;
                            a_offset = a_index;
                            b_offset = b_index;
                        }

                        score_gap_row -= gap_extend;
                        score_gap_col -= gap_extend;
                        dp_mem[b_index].best_gap = (score_val - gap_open_extend).max(score_gap_col);
                        score_gap_row = (score_val - gap_open_extend).max(score_gap_row);
                        dp_mem[b_index].best = score_val;
                    }
                }
                score_val = next_score;
            }
        }

        if first_b_index == b_size {
            break;
        }

        let b_index_mod = first_b_index % RESTRICT_SIZE;
        b_gap = first_b_index;
        if b_index_mod > 0 {
            b_gap += RESTRICT_SIZE - b_index_mod;
        }

        gap_dp_reserve_band(dp_mem, dp_mem_alloc, last_b_index, num_extra_cells);
        if last_b_index < b_size.saturating_sub(1) {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - x_dropoff && b_size <= len2 {
                dp_mem[b_size].best = score_gap_row;
                dp_mem[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size += 1;
            }
        }

        if b_size <= len2 {
            dp_mem[b_size].best = GAP_MININT;
            dp_mem[b_size].best_gap = GAP_MININT;
            b_size += 1;
        }
    }

    (a_offset, b_offset, best_score)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:364-733
// ```c
// static Int4 ALIGN_EX(const Uint1* A, const Uint1* B, Int4 M, Int4 N,
//                      Int4* a_offset, Int4* b_offset, GapPrelimEditBlock* edit_block, ...)
// {
//     ...
// }
// ```
fn align_ex_protein(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    score_matrix: BlastpScoreMatrix<'_>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    scratch: &mut GapAlignScratch,
) -> (usize, usize, i32, Vec<GapEditOp>) {
    if len1 == 0 || len2 == 0 {
        return (0, 0, 0, Vec::new());
    }

    let gap_open_extend = gap_open + gap_extend;
    let mut x_dropoff = x_drop;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend + 3) as usize
    } else {
        len2 + 3
    };
    gap_dp_reserve_initial(
        &mut scratch.dp_mem,
        &mut scratch.dp_mem_alloc,
        num_extra_cells,
    );
    gap_reset_traceback_state(scratch);

    let mut edit_script_num_rows = 100usize;
    let mut score = -gap_open_extend;
    scratch.dp_mem[0].best = 0;
    scratch.dp_mem[0].best_gap = -gap_open_extend;

    let row0 = gap_alloc_trace_row(
        &mut scratch.trace_rows,
        &mut scratch.trace_offsets,
        &mut scratch.trace_rows_used,
        0,
        0,
    );
    row0.reserve(num_extra_cells);
    row0.push(SCRIPT_GAP_IN_A);

    let mut b_size = 1usize;
    for i in 1..=len2.min(scratch.dp_mem_alloc - 1) {
        if score < -x_dropoff {
            break;
        }
        scratch.dp_mem[i].best = score;
        scratch.dp_mem[i].best_gap = score - gap_open_extend;
        score -= gap_extend;
        row0.push(SCRIPT_GAP_IN_A);
        b_size = i + 1;
    }

    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:471-482
    // ```c
    // if(reverse_sequence)
    //     matrix_row = matrix[ A[ M - a_index ] ];
    // else
    //     matrix_row = matrix[ A[ a_index ] ];
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:484-487
    // ```c
    // if(reverse_sequence)
    //     b_ptr = &B[N - first_b_index];
    // else
    //     b_ptr = &B[first_b_index];
    // ```
    let get_q = |index: usize| -> u8 {
        if reverse {
            q_seq[len1 - index]
        } else {
            q_seq[index]
        }
    };
    let get_s = |index: usize| -> u8 {
        if reverse {
            s_seq[len2 - 1 - index]
        } else {
            s_seq[index + 1]
        }
    };

    for a_index in 1..=len1 {
        if a_index >= edit_script_num_rows {
            edit_script_num_rows *= 2;
            if scratch.trace_rows.len() < edit_script_num_rows {
                scratch
                    .trace_rows
                    .reserve(edit_script_num_rows - scratch.trace_rows.len());
            }
            if scratch.trace_offsets.len() < edit_script_num_rows {
                scratch
                    .trace_offsets
                    .reserve(edit_script_num_rows - scratch.trace_offsets.len());
            }
        }

        let orig_b_index = first_b_index;
        let row_capacity = if gap_extend > 0 {
            b_size.saturating_sub(first_b_index) + num_extra_cells
        } else {
            len2.saturating_add(3).saturating_sub(first_b_index)
        };
        let edit_script_row = gap_alloc_trace_row(
            &mut scratch.trace_rows,
            &mut scratch.trace_offsets,
            &mut scratch.trace_rows_used,
            row_capacity,
            orig_b_index,
        );
        let qc = get_q(a_index);

        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;

        for b_index in first_b_index..b_size {
            let sc = if b_index < len2 { get_s(b_index) } else { 0 };
            let score_gap_col = scratch.dp_mem[b_index].best_gap;
            let next_score = scratch.dp_mem[b_index].best + score_matrix.score(qc, sc);

            let mut script = SCRIPT_SUB;
            if score_val < score_gap_col {
                script = SCRIPT_GAP_IN_B;
                score_val = score_gap_col;
            }
            if score_val < score_gap_row {
                script = SCRIPT_GAP_IN_A;
                score_val = score_gap_row;
            }

            if best_score - score_val > x_dropoff {
                if b_index == first_b_index {
                    first_b_index += 1;
                } else {
                    scratch.dp_mem[b_index].best = GAP_MININT;
                }
            } else {
                last_b_index = b_index;
                if score_val > best_score {
                    best_score = score_val;
                    a_offset = a_index;
                    b_offset = b_index;
                }

                score_gap_row -= gap_extend;
                let score_gap_col_ext = score_gap_col - gap_extend;
                let open_gap_col = score_val - gap_open_extend;
                if score_gap_col_ext < open_gap_col {
                    scratch.dp_mem[b_index].best_gap = open_gap_col;
                } else {
                    scratch.dp_mem[b_index].best_gap = score_gap_col_ext;
                    script += SCRIPT_EXTEND_GAP_B;
                }

                let open_gap_row = score_val - gap_open_extend;
                if score_gap_row < open_gap_row {
                    score_gap_row = open_gap_row;
                } else {
                    script += SCRIPT_EXTEND_GAP_A;
                }
                scratch.dp_mem[b_index].best = score_val;
            }

            score_val = next_score;
            let row_idx = b_index.saturating_sub(orig_b_index);
            if row_idx < edit_script_row.len() {
                edit_script_row[row_idx] = script;
            }
        }

        if first_b_index == b_size {
            break;
        }

        gap_dp_reserve_band(
            &mut scratch.dp_mem,
            &mut scratch.dp_mem_alloc,
            last_b_index,
            num_extra_cells,
        );
        let old_b_size = b_size;
        if last_b_index < b_size.saturating_sub(1) {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - x_dropoff
                && b_size <= len2
                && b_size < scratch.dp_mem_alloc
            {
                scratch.dp_mem[b_size].best = score_gap_row;
                scratch.dp_mem[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                let row_idx = b_size.saturating_sub(orig_b_index);
                if row_idx < edit_script_row.len() {
                    edit_script_row[row_idx] = SCRIPT_GAP_IN_A;
                } else {
                    edit_script_row.push(SCRIPT_GAP_IN_A);
                }
                b_size += 1;
            }
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:658-665
        // ```c
        // /* update the memory allocator to reflect the exact number
        //    of traceback cells this row needed */
        // state_struct->used += MAX(b_index, b_size) - orig_b_index + 1;
        // ```
        //
        // The C code only exposes the traceback cells written for this row.
        // Keep the Rust row length to the same written span so zero-filled
        // tail capacity cannot be read back as `SCRIPT_GAP_IN_A`.
        let used_cells = old_b_size
            .max(b_size)
            .saturating_sub(orig_b_index)
            .saturating_add(1);
        edit_script_row.truncate(used_cells);

        if b_size <= len2 && b_size < scratch.dp_mem_alloc {
            scratch.dp_mem[b_size].best = GAP_MININT;
            scratch.dp_mem[b_size].best_gap = GAP_MININT;
            b_size += 1;
        }
    }

    if best_score <= 0 {
        return (0, 0, 0, Vec::new());
    }

    let mut a_index = a_offset;
    let mut b_index = b_offset;
    let mut script = SCRIPT_SUB;
    let mut ops_reversed: Vec<(u8, u32)> = Vec::new();

    while a_index > 0 || b_index > 0 {
        let row_idx = a_index;
        if row_idx >= scratch.trace_rows_used {
            break;
        }
        let start_offset = scratch.trace_offsets[row_idx];
        let Some(b_rel) = b_index.checked_sub(start_offset) else {
            break;
        };
        let row = &scratch.trace_rows[row_idx];
        if b_rel >= row.len() {
            break;
        }
        let next_script = row[b_rel];
        match script & SCRIPT_OP_MASK {
            x if x == SCRIPT_GAP_IN_A => {
                script = next_script & SCRIPT_OP_MASK;
                if next_script & SCRIPT_EXTEND_GAP_A != 0 {
                    script = SCRIPT_GAP_IN_A;
                }
            }
            x if x == SCRIPT_GAP_IN_B => {
                script = next_script & SCRIPT_OP_MASK;
                if next_script & SCRIPT_EXTEND_GAP_B != 0 {
                    script = SCRIPT_GAP_IN_B;
                }
            }
            _ => {
                script = next_script & SCRIPT_OP_MASK;
            }
        }

        let op_type = script & SCRIPT_OP_MASK;
        if op_type == SCRIPT_GAP_IN_A {
            if b_index > 0 {
                b_index -= 1;
            }
        } else if op_type == SCRIPT_GAP_IN_B {
            if a_index > 0 {
                a_index -= 1;
            }
        } else {
            if a_index > 0 {
                a_index -= 1;
            }
            if b_index > 0 {
                b_index -= 1;
            }
        }
        push_trace_op(&mut ops_reversed, op_type);
    }

    let mut edit_ops = Vec::new();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:689-726
    // ```c
    // while (a_index > 0 || b_index > 0) {
    //     ...
    //     GapPrelimEditBlockAdd(edit_block, (EGapAlignOpType)script, 1);
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2494-2534
    // ```c
    // /* The fwd_prelim_tback script will get reversed here as the traceback
    //  * started from the highest scoring point and worked backwards. The
    //  * rev_prelim_tback script does NOT get reversed. Since it was reversed
    //  * when the traceback was produced it's already "forward". */
    // ```
    //
    // `ALIGN_EX` must therefore return the preliminary traceback in the
    // `GapPrelimEditBlockAdd` insertion order for both directions. Only
    // `Blast_PrelimEditBlockToGapEditScript` reverses the forward half.
    for (op_type, count) in ops_reversed {
        match op_type {
            x if x == SCRIPT_SUB => edit_ops.push(GapEditOp::Sub(count)),
            x if x == SCRIPT_GAP_IN_A => edit_ops.push(GapEditOp::Del(count)),
            x if x == SCRIPT_GAP_IN_B => edit_ops.push(GapEditOp::Ins(count)),
            _ => edit_ops.push(GapEditOp::Sub(count)),
        }
    }

    (a_offset, b_offset, best_score, edit_ops)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4549-4712
// ```c
// Int2 BLAST_GappedAlignmentWithTraceback(...);
// ```
pub(crate) fn blast_gapped_alignment_with_traceback(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    matrix: ScoringMatrix,
    adjusted_matrix: Option<&AdjustedProteinMatrix>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> Option<BlastpGapAlignResult> {
    if q_start >= query.len() || s_start >= subject.len() {
        return None;
    }

    let score_matrix = adjusted_matrix.map_or(BlastpScoreMatrix::Standard(matrix), |matrix| {
        BlastpScoreMatrix::Adjusted(matrix)
    });

    let mut scratch = GapAlignScratch::new();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4599-4611
    // ```c
    // score_left = ALIGN_EX(query, subject, q_start+1, s_start+1,
    //                       &private_q_length, &private_s_length, rev_prelim_tback, ...,
    //                       FALSE, TRUE, fence_hit);
    // gap_align->query_start = q_start - private_q_length + 1;
    // gap_align->subject_start = s_start - private_s_length + 1;
    // ```
    let (left_q, left_s, mut score_left, left_edit_ops) = align_ex_protein(
        query,
        subject,
        q_start + 1,
        s_start + 1,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop,
        true,
        &mut scratch,
    );
    let mut query_start = q_start + 1 - left_q;
    let mut subject_start = s_start + 1 - left_s;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4629-4642
    // ```c
    // score_right = ALIGN_EX(query+q_start, subject+s_start,
    //                        q_length-q_start-1, s_length-s_start-1, ...,
    //                        FALSE, FALSE, fence_hit);
    // gap_align->query_stop = q_start + private_q_length + 1;
    // gap_align->subject_stop = s_start + private_s_length + 1;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4624-4642
    // ```c
    // if ((! (fence_hit && *fence_hit)) &&
    //     (q_start < q_length) &&
    //     (s_start < s_length)) {
    //    found_end = TRUE;
    //    score_right = ALIGN_EX(query+q_start, subject+s_start,
    //       q_length-q_start-1, s_length-s_start-1, ...);
    //    gap_align->query_stop = q_start + private_q_length + 1;
    //    gap_align->subject_stop = s_start + private_s_length + 1;
    // }
    // ```
    let (mut query_stop, mut subject_stop, mut score_right, right_edit_ops) =
        if q_start < query.len() && s_start < subject.len() {
            let (right_q, right_s, score_right, right_edit_ops) = align_ex_protein(
                &query[q_start..],
                &subject[s_start..],
                query.len().saturating_sub(q_start + 1),
                subject.len().saturating_sub(s_start + 1),
                score_matrix,
                gap_open,
                gap_extend,
                x_drop,
                false,
                &mut scratch,
            );
            (
                q_start + right_q + 1,
                s_start + right_s + 1,
                score_right,
                right_edit_ops,
            )
        } else {
            (
                q_start.saturating_sub(1),
                s_start.saturating_sub(1),
                0,
                Vec::new(),
            )
        };

    let mut edit_script = prelim_edit_ops_to_gap_edit_script(left_edit_ops, right_edit_ops);
    prune_terminal_gap_ops(
        &mut edit_script,
        gap_open,
        gap_extend,
        &mut score_left,
        &mut score_right,
        &mut query_start,
        &mut query_stop,
        &mut subject_start,
        &mut subject_stop,
    );

    if edit_script.is_empty() || query_stop <= query_start || subject_stop <= subject_start {
        return None;
    }

    let score = score_left + score_right;

    let (num_ident, num_positives, mismatches, gap_opens, gap_letters) =
        stats_from_edit_ops_protein(
            query,
            subject,
            query_start,
            subject_start,
            &edit_script,
            matrix,
        );

    Some(BlastpGapAlignResult {
        query_start,
        query_stop,
        subject_start,
        subject_stop,
        score,
        edit_script,
        num_ident,
        num_positives,
        mismatches,
        gap_opens,
        gap_letters,
    })
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319
// ```c
// static Int2
// s_BlastProtGappedAlignment(..., Boolean restricted_alignment, ...)
// {
//    ...
//    gap_align->query_start = q_length - private_q_start;
//    gap_align->subject_start = s_length - private_s_start + subject_shift;
//    ...
//    gap_align->query_stop += init_hsp->offsets.qs_offsets.q_off + 1;
//    gap_align->subject_stop += init_hsp->offsets.qs_offsets.s_off + 1;
//    ...
//    gap_align->score = score_right+score_left;
// }
// ```
pub(crate) fn blastp_score_only_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    mode: BlastpGappedAlignmentMode,
) -> BlastpScoreOnlyGappedAlignment {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4251-4258
    // ```c
    // q_length = init_hsp->offsets.qs_offsets.q_off + 1;
    // s_length = init_hsp->offsets.qs_offsets.s_off + 1;
    // AdjustSubjectRange(&s_length, &subject_length, q_length, query_length,
    //                    &subject_shift);
    // ```
    assert!(
        q_start < query.len() && s_start < subject.len(),
        "NCBI BLAST requires preliminary gapped starts to lie within query/subject bounds"
    );

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4247-4258
    // ```c
    // q_length = init_hsp->offsets.qs_offsets.q_off + 1;
    // s_length = init_hsp->offsets.qs_offsets.s_off + 1;
    // AdjustSubjectRange(&s_length, &subject_length, q_length, query_length,
    //                    &subject_shift);
    // ```
    let q_length = q_start + 1;
    let s_length = s_start + 1;
    let (subject_offset, subject_length, subject_shift) =
        adjust_subject_range(s_length, subject.len(), q_length, query.len());
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4257-4258
    // ```c
    // AdjustSubjectRange(&s_length, &subject_length, q_length, query_length,
    //                    &subject_shift);
    // ```
    assert!(
        subject_shift < subject.len() && subject_length > 0,
        "NCBI BLAST requires AdjustSubjectRange to yield a non-empty subject window"
    );
    let mut scratch = GapAlignScratch::new();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4266-4278
    // ```c
    // score_left = s_RestrictedGappedAlign(query, subject+subject_shift,
    //                                      q_length, s_length, ...);
    // ...
    // score_left = Blast_SemiGappedAlign(query, subject+subject_shift,
    //                                    q_length, s_length, ...);
    // ```
    let (private_q_start, private_s_start, score_left) = match mode {
        BlastpGappedAlignmentMode::Exact => align_ex_protein_score_only(
            query,
            0,
            subject,
            subject_shift,
            q_length,
            subject_offset,
            BlastpScoreMatrix::Standard(matrix),
            gap_open,
            gap_extend,
            x_drop,
            true,
            &mut scratch,
        ),
        BlastpGappedAlignmentMode::Restricted => restricted_gapped_align_protein(
            query,
            0,
            subject,
            subject_shift,
            q_length,
            subject_offset,
            BlastpScoreMatrix::Standard(matrix),
            gap_open,
            gap_extend,
            x_drop,
            true,
            &mut scratch.dp_mem,
            &mut scratch.dp_mem_alloc,
        ),
    };
    let query_hit_start = i32::try_from(q_length - private_q_start)
        .expect("NCBI BLAST preliminary query start must fit in Int4");
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4280-4281
    // ```c
    // gap_align->query_start = q_length - private_q_start;
    // gap_align->subject_start = s_length - private_s_start + subject_shift;
    // ```
    let subject_hit_start = i32::try_from(subject_offset - private_s_start + subject_shift)
        .expect("NCBI BLAST preliminary subject start must fit in Int4");

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4300-4313
    // ```c
    // if (q_length < query_length && s_length < subject_length) {
    //     score_right = ...(..., query_length-q_length, subject_length-s_length, ...);
    // }
    // gap_align->query_stop += init_hsp->offsets.qs_offsets.q_off + 1;
    // gap_align->subject_stop += init_hsp->offsets.qs_offsets.s_off + 1;
    // ```
    let (query_hit_stop, subject_hit_stop, score_right) =
        if q_length < query.len() && subject_offset < subject_length {
            let (private_q_stop, private_s_stop, score_right) = match mode {
                BlastpGappedAlignmentMode::Exact => align_ex_protein_score_only(
                    query,
                    q_start,
                    subject,
                    s_start,
                    query.len() - q_length,
                    subject_length - subject_offset,
                    BlastpScoreMatrix::Standard(matrix),
                    gap_open,
                    gap_extend,
                    x_drop,
                    false,
                    &mut scratch,
                ),
                BlastpGappedAlignmentMode::Restricted => restricted_gapped_align_protein(
                    query,
                    q_start,
                    subject,
                    s_start,
                    query.len() - q_length,
                    subject_length - subject_offset,
                    BlastpScoreMatrix::Standard(matrix),
                    gap_open,
                    gap_extend,
                    x_drop,
                    false,
                    &mut scratch.dp_mem,
                    &mut scratch.dp_mem_alloc,
                ),
            };
            (
                i32::try_from(q_start + private_q_stop + 1)
                    .expect("NCBI BLAST preliminary query stop must fit in Int4"),
                i32::try_from(s_start + private_s_stop + 1)
                    .expect("NCBI BLAST preliminary subject stop must fit in Int4"),
                score_right,
            )
        } else {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4342-4344
            // ```c
            // if (found_end == FALSE) {
            //     gap_align->query_stop = init_hsp->offsets.qs_offsets.q_off;
            //     gap_align->subject_stop = init_hsp->offsets.qs_offsets.s_off;
            // }
            // ```
            (
                i32::try_from(q_start)
                    .expect("NCBI BLAST preliminary no-right query stop must fit in Int4"),
                i32::try_from(s_start)
                    .expect("NCBI BLAST preliminary no-right subject stop must fit in Int4"),
                0,
            )
        };

    BlastpScoreOnlyGappedAlignment {
        query_start: query_hit_start,
        query_stop: query_hit_stop,
        subject_start: subject_hit_start,
        subject_stop: subject_hit_stop,
        score: score_left + score_right,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1935-1948
// ```c
// BLAST_GappedAlignmentWithTraceback(context->prog_number,
//                                    query_data->data,
//                                    subject_data->data, gapAlign,
//                                    context->scoringParams,
//                                    q_start, s_start,
//                                    query_data->length,
//                                    subject_data->length,
//                                    &fence_hit);
// ```
pub(crate) fn traceback_preliminary_blastp_hsp(
    preliminary_hsp: &BlastpPreliminaryHsp,
    query: &[u8],
    subject: &[u8],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    params: &KarlinParams,
    effective_space: f64,
) -> Hit {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1924-1946
    // ```c
    // q_start = hsp->query.gapped_start - query_range->begin;
    // s_start = hsp->subject.gapped_start - subject_range->begin;
    // status = BLAST_GappedAlignmentWithTraceback(..., q_start, s_start, ...);
    // ```
    let q_start = usize::try_from(preliminary_hsp.gapped_query_start)
        .expect("NCBI BLAST traceback requires non-negative preliminary gapped query start");
    let s_start = usize::try_from(preliminary_hsp.gapped_subject_start)
        .expect("NCBI BLAST traceback requires non-negative preliminary gapped subject start");
    let aligned = blast_gapped_alignment_with_traceback(
        query, subject, q_start, s_start, matrix, None, gap_open, gap_extend, x_drop,
    )
    .expect("NCBI BLAST traceback helper must produce a gapped alignment");

    let alignment_len = aligned
        .edit_script
        .iter()
        .map(|op| op.num() as usize)
        .sum::<usize>();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:151-182
    // ```c
    // Blast_HSPInit(Int4 query_start, Int4 query_end, Int4 subject_start,
    //               Int4 subject_end, ...)
    // {
    //    new_hsp->query.offset = query_start;
    //    new_hsp->query.end = query_end;
    //    new_hsp->subject.offset = subject_start;
    //    new_hsp->subject.end = subject_end;
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:151-178
    // ```c
    // new_hsp->query.offset = query_start;
    // new_hsp->subject.offset = subject_start;
    // new_hsp->query.end = query_end;
    // new_hsp->subject.end = subject_end;
    // ```
    assert!(
        alignment_len > 0,
        "NCBI BLAST traceback helper must not discard zero-length alignments as None"
    );

    let identity = 100.0 * (aligned.num_ident as f64) / (alignment_len as f64);
    Hit {
        identity,
        length: alignment_len,
        mismatch: aligned.mismatches,
        gapopen: aligned.gap_opens,
        q_start: aligned.query_start + 1,
        q_end: aligned.query_stop,
        s_start: aligned.subject_start + 1,
        s_end: aligned.subject_stop,
        e_value: evalue_from_raw_score(aligned.score, params, effective_space),
        bit_score: calc_bit_score(aligned.score, params),
        num_ident: aligned.num_ident,
        query_frame: preliminary_hsp.query_frame,
        query_length: preliminary_hsp.query_length,
        q_idx: preliminary_hsp.q_idx,
        s_idx: preliminary_hsp.s_idx,
        raw_score: aligned.score,
        gap_info: Some(aligned.edit_script),
        num_positives: aligned.num_positives,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:647-662
// ```c
// for (index=extra_start; index < hsp_list->hspcnt; index++) {
//     ...
//     delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(hsp, query, ...);
//     if (!delete_hsp)
//         delete_hsp = Blast_HSPTestIdentityAndLength(program_number, hsp,
//                                                     query_nomask, subject, ...);
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:234-250
// ```c
// Blast_HSPListGetEvalues(program_number, query_info, subject_length,
//                         hsp_list, kGapped, FALSE, sbp, 0, scale_factor);
// Blast_HSPListReapByEvalue(hsp_list, hit_params->options);
// ...
// Blast_HSPListGetBitScores(hsp_list, kGapped, sbp);
// ```
fn refresh_blastp_hit_statistics(
    hit: &mut Hit,
    query_nomask: &[u8],
    subject: &[u8],
    matrix: ScoringMatrix,
) -> bool {
    let Some(view) = build_alignment_view_with_matrix(query_nomask, subject, hit, matrix) else {
        return false;
    };
    hit.length = view.alignment_len;
    hit.num_ident = view.num_ident;
    hit.num_positives = view.num_positives;
    hit.gapopen = view.gap_opens;
    hit.mismatch = view
        .alignment_len
        .saturating_sub(view.num_ident.saturating_add(view.gap_letters));
    hit.identity = if view.alignment_len > 0 {
        (view.num_ident as f64 * 100.0) / view.alignment_len as f64
    } else {
        0.0
    };
    true
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:440-477
// ```c
// static Boolean
// s_UpdateReevaluatedHSP(BlastHSP* hsp, Boolean gapped,
//                        Int4 cutoff_score,
//                        Int4 score, const Uint1* query_start, const Uint1* subject_start,
//                        const Uint1* best_q_start, const Uint1* best_q_end,
//                        const Uint1* best_s_start, const Uint1* best_s_end,
//                        int best_start_esp_index,
//                        int best_end_esp_index,
//                        int best_end_esp_num)
// {
//     ...
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:440-642
// ```c
// static Boolean s_UpdateReevaluatedHSP(...);
// Boolean Blast_HSPReevaluateWithAmbiguitiesGapped(...);
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:399-414
// ```c
// matrix = gap_align->sbp->matrix->data;
// ```
pub(crate) fn reevaluate_blastp_hit_with_traceback(
    hit: &mut Hit,
    query: &[u8],
    query_nomask: &[u8],
    subject: &[u8],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    cutoff_score: i32,
) -> bool {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_options.c:787-804
    // ```c
    // if (!Blast_ProgramIsNucleotide(program_number)) {/*protein-protein options.*/
    //     (*options)->shift_pen = INT2_MAX;
    //     (*options)->is_ooframe = FALSE;
    //     (*options)->gap_open = BLAST_GAP_OPEN_PROT;
    //     (*options)->gap_extend = BLAST_GAP_EXTN_PROT;
    //     (*options)->matrix = strdup(BLAST_DEFAULT_MATRIX);
    // } else {
    //     (*options)->penalty = BLAST_PENALTY;
    //     (*options)->reward = BLAST_REWARD;
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:555-564
    // ```c
    // params->options = (BlastScoringOptions *)score_options;
    // ...
    // params->reward = score_options->reward;
    // params->penalty = score_options->penalty;
    // ```
    //
    // For blastp, `BlastScoringOptionsNew` leaves `reward` at zero, and
    // `BlastScoringParametersNew` copies that zero into `score_params->reward`.
    let blastp_reward = 0i32;
    let gap_ops = match hit.gap_info.as_mut() {
        Some(info) if !info.is_empty() => info,
        _ => return true,
    };

    let q_offset = hit.q_start.saturating_sub(1);
    let s_offset = hit.s_start.min(hit.s_end).saturating_sub(1);
    let mut query_pos = q_offset;
    let mut subject_pos = s_offset;

    let mut score = 0i32;
    let mut sum = 0i32;
    let mut best_q_start = query_pos;
    let mut best_q_end = query_pos;
    let mut best_s_start = subject_pos;
    let mut best_s_end = subject_pos;
    let mut current_q_start = query_pos;
    let mut current_s_start = subject_pos;
    let mut best_start_esp_index = 0usize;
    let mut best_end_esp_index = 0usize;
    let mut current_start_esp_index = 0usize;
    let mut best_end_esp_num = None;

    for index in 0..gap_ops.len() {
        let mut op_index = 0usize;
        while op_index < gap_ops[index].num() as usize {
            match gap_ops[index] {
                GapEditOp::Sub(_) => {
                    if query_pos >= query.len() || subject_pos >= subject.len() {
                        return true;
                    }
                    sum += protein_score(matrix, query[query_pos], subject[subject_pos]);
                    query_pos += 1;
                    subject_pos += 1;
                    op_index += 1;
                }
                GapEditOp::Del(n) => {
                    let n = n as usize;
                    sum -= gap_open + gap_extend * n as i32;
                    subject_pos += n;
                    op_index += n;
                }
                GapEditOp::Ins(n) => {
                    let n = n as usize;
                    sum -= gap_open + gap_extend * n as i32;
                    query_pos += n;
                    op_index += n;
                }
            }

            if sum < 0 {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:562-569
                // ```c
                // if (op_index < esp->num[index]) {
                //     esp->num[index] -= op_index;
                //     current_start_esp_index = index;
                //     op_index = 0;
                // } else {
                //     current_start_esp_index = index + 1;
                // }
                // ```
                if op_index < gap_ops[index].num() as usize {
                    if let GapEditOp::Sub(num) = &mut gap_ops[index] {
                        *num -= op_index as u32;
                    }
                    current_start_esp_index = index;
                    op_index = 0;
                } else {
                    current_start_esp_index = index + 1;
                }
                sum = 0;
                current_q_start = query_pos;
                current_s_start = subject_pos;
                if score < cutoff_score {
                    best_q_start = query_pos;
                    best_q_end = query_pos;
                    best_s_start = subject_pos;
                    best_s_end = subject_pos;
                    score = 0;
                    best_start_esp_index = current_start_esp_index;
                    best_end_esp_index = current_start_esp_index;
                    best_end_esp_num = None;
                }
            } else if sum > score {
                score = sum;
                best_q_start = current_q_start;
                best_s_start = current_s_start;
                best_q_end = query_pos;
                best_s_end = subject_pos;
                best_start_esp_index = current_start_esp_index;
                best_end_esp_index = index;
                best_end_esp_num = Some(op_index);
            }
        }
    }

    hit.raw_score = score;
    let Some(mut best_end_esp_num) = best_end_esp_num else {
        return true;
    };
    if hit.raw_score < cutoff_score
        || best_start_esp_index >= gap_ops.len()
        || best_end_esp_index >= gap_ops.len()
    {
        return true;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:618-637
    // ```c
    // if (best_start_esp_index < esp->size && best_end_esp_index < esp->size) {
    //     ASSERT(esp->op_type[best_start_esp_index] == eGapAlignSub);
    //     ASSERT(esp->op_type[best_end_esp_index] == eGapAlignSub);
    //
    //     qp = (Int4)(best_q_start - q);
    //     sp = (Int4)(best_s_start - s);
    //     ext = 0;
    //     while(qp > 0 && sp > 0 && (q[--qp] == s[--sp]) && q[qp]<4) ext++;
    //     best_q_start -= ext;
    //     best_s_start -= ext;
    //     esp->num[best_start_esp_index] += ext;
    //     if (best_end_esp_index == best_start_esp_index) best_end_esp_num += ext;
    //     score += ext * score_params->reward;
    //
    //     qp = (Int4)(best_q_end - q);
    //     sp = (Int4)(best_s_end - s);
    //     ext = 0;
    //     while(qp < qlen && sp < slen && q[qp]<4 && (q[qp++] == s[sp++])) ext++;
    //     best_q_end += ext;
    //     best_s_end += ext;
    //     esp->num[best_end_esp_index] += ext;
    //     best_end_esp_num += ext;
    //     score += ext * score_params->reward;
    // }
    // ```
    if matches!(gap_ops.get(best_start_esp_index), Some(GapEditOp::Sub(_)))
        && matches!(gap_ops.get(best_end_esp_index), Some(GapEditOp::Sub(_)))
    {
        let mut qp = best_q_start;
        let mut sp = best_s_start;
        let mut ext = 0usize;
        while qp > 0 && sp > 0 {
            let next_q = qp - 1;
            let next_s = sp - 1;
            if query[next_q] != subject[next_s] || query[next_q] >= 4 {
                break;
            }
            qp = next_q;
            sp = next_s;
            ext += 1;
        }
        best_q_start -= ext;
        best_s_start -= ext;
        if let Some(GapEditOp::Sub(num)) = gap_ops.get_mut(best_start_esp_index) {
            *num += ext as u32;
        }
        if best_end_esp_index == best_start_esp_index {
            best_end_esp_num += ext;
        }
        score += ext as i32 * blastp_reward;

        let mut qp = best_q_end;
        let mut sp = best_s_end;
        let mut ext = 0usize;
        while qp < query.len() && sp < subject.len() {
            if query[qp] >= 4 || query[qp] != subject[sp] {
                break;
            }
            qp += 1;
            sp += 1;
            ext += 1;
        }
        best_q_end += ext;
        best_s_end += ext;
        if let Some(GapEditOp::Sub(num)) = gap_ops.get_mut(best_end_esp_index) {
            *num += ext as u32;
        }
        best_end_esp_num += ext;
        score += ext as i32 * blastp_reward;
    }

    hit.q_start = best_q_start + 1;
    hit.q_end = best_q_end;
    hit.s_start = best_s_start + 1;
    hit.s_end = best_s_end;

    if best_end_esp_index != gap_ops.len().saturating_sub(1) || best_start_esp_index > 0 {
        let new_ops = gap_ops[best_start_esp_index..=best_end_esp_index].to_vec();
        *gap_ops = new_ops;
    }
    if let Some(last) = gap_ops.last_mut() {
        match last {
            GapEditOp::Sub(n) | GapEditOp::Del(n) | GapEditOp::Ins(n) => {
                *n = best_end_esp_num as u32;
            }
        }
    }
    if gap_ops.is_empty() || gap_ops.last().is_some_and(|op| op.num() == 0) {
        return true;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:649-660
    // ```c
    // query = query_blk->sequence + query_info->contexts[hsp->context].query_offset;
    // query_nomask = query_blk->sequence_nomask +
    //     query_info->contexts[hsp->context].query_offset;
    // delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(hsp, query, ...);
    // if (!delete_hsp)
    //     delete_hsp = Blast_HSPTestIdentityAndLength(program_number, hsp,
    //                                                 query_nomask, subject, ...);
    // ```
    !refresh_blastp_hit_statistics(hit, query_nomask, subject, matrix)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::composition_adjustment::adjust_scores::{
        blast_adjust_scores, build_matrix_info, read_aa_composition,
    };
    use crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode;
    use crate::utils::matrix::{aa_char_to_ncbistdaa, ncbistdaa};
    use crate::utils::seg::{SegMasker, SegParams};
    use std::fs;

    fn fasta_sequence_by_id(contents: &str, wanted_id: &str) -> Vec<u8> {
        let mut current_id = None::<&str>;
        let mut sequence = Vec::new();
        for line in contents.lines() {
            if let Some(header) = line.strip_prefix('>') {
                let id = header.split_whitespace().next().unwrap_or_default();
                if current_id == Some(wanted_id) {
                    break;
                }
                current_id = Some(id);
                continue;
            }
            if current_id == Some(wanted_id) {
                sequence.extend(line.trim().bytes().map(aa_char_to_ncbistdaa));
            }
        }
        sequence
    }

    fn mask_query_sequence(sequence: &[u8]) -> Vec<u8> {
        let mut masked = sequence.to_vec();
        let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));
        for interval in masker.mask_sequence(sequence) {
            for pos in interval.start..interval.end {
                masked[pos] = ncbistdaa::X;
            }
        }
        masked
    }

    fn mask_subject_sequence(sequence: &[u8]) -> Vec<u8> {
        let mut masked = sequence.to_vec();
        let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));
        for interval in masker.mask_sequence(sequence) {
            for pos in interval.start..interval.end {
                masked[pos] = ncbistdaa::X;
            }
        }
        masked
    }

    #[test]
    fn test_adjust_subject_range_keeps_short_subject() {
        let (offset, length, shift) = adjust_subject_range(100, 1000, 50, 400);
        assert_eq!(offset, 100);
        assert_eq!(length, 1000);
        assert_eq!(shift, 0);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4176-4190
    // ```c
    // s_offset = *subject_offset_ptr;
    // max_extension_left = query_offset + MAX_TOTAL_GAPS;
    // ...
    // *start_shift = s_offset - max_extension_left;
    // *subject_offset_ptr = max_extension_left;
    // *subject_length_ptr =
    //    MIN(subject_length, s_offset + max_extension_right) - *start_shift;
    // ```
    #[test]
    fn test_adjust_subject_range_shifts_long_subject_ncbi_style() {
        let subject_offset = 7_000usize;
        let subject_length = 100_000usize;
        let query_offset = 5usize;
        let query_length = 8usize;
        let (adjusted_offset, adjusted_length, start_shift) =
            adjust_subject_range(subject_offset, subject_length, query_offset, query_length);

        assert_eq!(adjusted_offset, query_offset + MAX_TOTAL_GAPS);
        assert_eq!(
            start_shift,
            subject_offset - (query_offset + MAX_TOTAL_GAPS)
        );
        assert_eq!(
            adjusted_length,
            subject_length.min(subject_offset + 3_003) - start_shift
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2488-2533
    // ```c
    // /* The fwd_prelim_tback script will get reversed here ... */
    // if (merge_ops)
    //     esp->num[index-1] += fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].num;
    // for (; i >= 0; i--) {
    //     op = fwd_prelim_tback->edit_ops + i;
    //     esp->op_type[index] = op->op_type;
    //     esp->num[index] = op->num;
    // }
    // ```
    #[test]
    fn test_prelim_edit_ops_to_gap_edit_script_reverses_forward_ops_ncbi_style() {
        let merged = prelim_edit_ops_to_gap_edit_script(
            vec![GapEditOp::Sub(2), GapEditOp::Del(1)],
            vec![GapEditOp::Ins(3), GapEditOp::Del(4), GapEditOp::Sub(5)],
        );
        assert_eq!(
            merged,
            vec![
                GapEditOp::Sub(2),
                GapEditOp::Del(1),
                GapEditOp::Sub(5),
                GapEditOp::Del(4),
                GapEditOp::Ins(3),
            ]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2496-2533
    // ```c
    // if (merge_ops)
    //     esp->num[index-1] += fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].num;
    // /* If we merge, then we skip the first one. */
    // if (merge_ops)
    //     i = fwd_prelim_tback->num_ops - 2;
    // ```
    #[test]
    fn test_prelim_edit_ops_to_gap_edit_script_merges_last_forward_op_ncbi_style() {
        let merged = prelim_edit_ops_to_gap_edit_script(
            vec![GapEditOp::Sub(2), GapEditOp::Del(1)],
            vec![GapEditOp::Ins(3), GapEditOp::Del(4), GapEditOp::Del(5)],
        );
        assert_eq!(
            merged,
            vec![
                GapEditOp::Sub(2),
                GapEditOp::Del(6),
                GapEditOp::Del(4),
                GapEditOp::Ins(3),
            ]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4681-4707
    // ```c
    // while (esp->size && esp->op_type[0] != eGapAlignSub) {
    //     score_left += score_params->gap_open + esp->num[0] * score_params->gap_extend;
    //     ...
    // }
    // while (esp->size && esp->op_type[i-1] != eGapAlignSub) {
    //     score_right += score_params->gap_open + esp->num[i-1] * score_params->gap_extend;
    //     ...
    // }
    // ```
    #[test]
    fn test_prune_terminal_gap_ops_updates_score_and_boundaries_ncbi_style() {
        let mut edit_script = vec![GapEditOp::Ins(2), GapEditOp::Sub(3), GapEditOp::Del(1)];
        let mut score_left = 30;
        let mut score_right = 40;
        let mut query_start = 10usize;
        let mut query_stop = 15usize;
        let mut subject_start = 20usize;
        let mut subject_stop = 24usize;

        prune_terminal_gap_ops(
            &mut edit_script,
            11,
            1,
            &mut score_left,
            &mut score_right,
            &mut query_start,
            &mut query_stop,
            &mut subject_start,
            &mut subject_stop,
        );

        assert_eq!(edit_script, vec![GapEditOp::Sub(3)]);
        assert_eq!(score_left, 43);
        assert_eq!(score_right, 52);
        assert_eq!(query_start, 12);
        assert_eq!(query_stop, 15);
        assert_eq!(subject_start, 20);
        assert_eq!(subject_stop, 23);
    }

    #[test]
    fn test_blastp_traceback_builds_gap_script() {
        let query = [1u8, 2, 3, 4, 5];
        let subject = [1u8, 2, 9, 4, 5];
        let result = blast_gapped_alignment_with_traceback(
            &query,
            &subject,
            1,
            1,
            ScoringMatrix::Blosum62,
            None,
            11,
            1,
            25,
        )
        .expect("traceback alignment");
        assert!(!result.edit_script.is_empty());
        assert!(result.score > 0);
        assert!(result.query_stop > result.query_start);
    }

    #[test]
    fn test_reevaluate_blastp_hit_keeps_best_subalignment() {
        let query = [1u8, 2, 3, 4];
        let subject = [1u8, 2, 21, 4];
        let mut hit = Hit {
            identity: 0.0,
            length: 4,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 4,
            s_start: 1,
            s_end: 4,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: 0,
            query_length: 4,
            q_idx: 0,
            s_idx: 0,
            raw_score: 0,
            gap_info: Some(vec![GapEditOp::Sub(4)]),
            num_positives: 0,
        };
        let delete = reevaluate_blastp_hit_with_traceback(
            &mut hit,
            &query,
            &query,
            &subject,
            ScoringMatrix::Blosum62,
            11,
            1,
            1,
        );
        assert!(!delete);
        assert!(hit.raw_score > 0);
        assert!(hit.q_end >= hit.q_start);
    }

    #[test]
    fn test_reevaluate_blastp_hit_trims_partial_leading_sub_run() {
        let matrix = ScoringMatrix::Blosum62;
        let mismatch = (1u8..=u8::MAX)
            .find(|&residue| protein_score(matrix, 1, residue) < 0)
            .expect("negative mismatch for residue 1");
        let match_score = protein_score(matrix, 1, 1);
        let query = [1u8, 1, 1, 1];
        let subject = [mismatch, 1, 1, 1];
        let mut hit = Hit {
            identity: 0.0,
            length: 4,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 4,
            s_start: 1,
            s_end: 4,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: 0,
            query_length: 4,
            q_idx: 0,
            s_idx: 0,
            raw_score: 0,
            gap_info: Some(vec![GapEditOp::Sub(4)]),
            num_positives: 0,
        };
        let delete = reevaluate_blastp_hit_with_traceback(
            &mut hit, &query, &query, &subject, matrix, 11, 1, 1,
        );

        assert!(!delete);
        assert_eq!(hit.q_start, 2);
        assert_eq!(hit.q_end, 4);
        assert_eq!(hit.s_start, 2);
        assert_eq!(hit.s_end, 4);
        assert_eq!(hit.raw_score, match_score * 3);
        assert_eq!(hit.gap_info, Some(vec![GapEditOp::Sub(3)]));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:618-637
    // ```c
    // while(qp > 0 && sp > 0 && (q[--qp] == s[--sp]) && q[qp]<4) ext++;
    // best_q_start -= ext;
    // best_s_start -= ext;
    // esp->num[best_start_esp_index] += ext;
    // ```
    #[test]
    fn test_reevaluate_blastp_hit_extends_left_over_exact_low_residue_match() {
        let matrix = ScoringMatrix::Blosum62;
        let mismatch = (1u8..=u8::MAX)
            .find(|&residue| protein_score(matrix, 5, residue) < 0)
            .expect("negative mismatch for residue 5");
        let match_score = protein_score(matrix, 1, 1);
        let query = [9u8, 9, 1, 1, 5];
        let subject = [1u8, 1, mismatch];
        let mut hit = Hit {
            identity: 0.0,
            length: 3,
            mismatch: 0,
            gapopen: 1,
            q_start: 3,
            q_end: 5,
            s_start: 2,
            s_end: 3,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: 0,
            query_length: query.len(),
            q_idx: 0,
            s_idx: 0,
            raw_score: 0,
            gap_info: Some(vec![GapEditOp::Ins(1), GapEditOp::Sub(2)]),
            num_positives: 0,
        };

        let delete = reevaluate_blastp_hit_with_traceback(
            &mut hit, &query, &query, &subject, matrix, 11, 1, 1,
        );

        assert!(!delete);
        assert_eq!(hit.q_start, 3);
        assert_eq!(hit.q_end, 4);
        assert_eq!(hit.s_start, 1);
        assert_eq!(hit.s_end, 2);
        assert_eq!(hit.raw_score, match_score);
        assert_eq!(hit.gap_info, Some(vec![GapEditOp::Sub(2)]));
    }

    #[test]
    fn test_traceback_keeps_seed_when_start_is_last_residue() {
        let query = [1u8, 2, 3];
        let subject = [8u8, 9, 3];
        let result = blast_gapped_alignment_with_traceback(
            &query,
            &subject,
            2,
            2,
            ScoringMatrix::Blosum62,
            None,
            11,
            1,
            25,
        )
        .expect("traceback keeps terminal seed");
        assert_eq!(result.query_start, 2);
        assert_eq!(result.query_stop, 3);
        assert_eq!(result.subject_start, 2);
        assert_eq!(result.subject_stop, 3);
    }

    #[test]
    fn test_score_only_alignment_keeps_seed_when_start_is_last_residue() {
        let query = [1u8, 2, 3];
        let subject = [8u8, 9, 3];
        let result = blastp_score_only_gapped_alignment(
            &query,
            &subject,
            2,
            2,
            ScoringMatrix::Blosum62,
            11,
            1,
            25,
            BlastpGappedAlignmentMode::Exact,
        );
        assert_eq!(result.query_start, 2);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4338-4343
        // ```c
        // if (found_end == FALSE) {
        //     gap_align->query_stop = init_hsp->offsets.qs_offsets.q_off;
        //     gap_align->subject_stop = init_hsp->offsets.qs_offsets.s_off;
        // }
        // ```
        assert_eq!(result.query_stop, 2);
        assert_eq!(result.subject_start, 2);
        assert_eq!(result.subject_stop, 2);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4257-4291
    // ```c
    // AdjustSubjectRange(&s_length, &subject_length, q_length, query_length,
    //                    &subject_shift);
    // ...
    // score_left = Blast_SemiGappedAlign(query, subject+subject_shift,
    //                                    q_length, s_length, ...);
    // gap_align->query_start = q_length - private_q_start;
    // gap_align->subject_start = s_length - private_s_start + subject_shift;
    // ```
    #[test]
    fn test_score_only_alignment_exact_preserves_original_subject_coordinates_with_shift() {
        let query = vec![1u8; 8];
        let subject = vec![1u8; 100_000];
        let q_start = 5usize;
        let s_start = 7_000usize;

        let result = blastp_score_only_gapped_alignment(
            &query,
            &subject,
            q_start,
            s_start,
            ScoringMatrix::Blosum62,
            11,
            1,
            25,
            BlastpGappedAlignmentMode::Exact,
        );

        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_stop, 8);
        assert_eq!(result.subject_start, 6_995);
        assert_eq!(result.subject_stop, 7_003);
        assert!(result.score > 0);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4273-4291
    // ```c
    // if (restricted_alignment) {
    //    score_left = s_RestrictedGappedAlign(query, subject+subject_shift,
    //                                         q_length, s_length, ...);
    // }
    // gap_align->subject_start = s_length - private_s_start + subject_shift;
    // ```
    #[test]
    fn test_score_only_alignment_restricted_preserves_original_subject_coordinates_with_shift() {
        let query = vec![1u8; 8];
        let subject = vec![1u8; 100_000];
        let q_start = 5usize;
        let s_start = 7_000usize;

        let result = blastp_score_only_gapped_alignment(
            &query,
            &subject,
            q_start,
            s_start,
            ScoringMatrix::Blosum62,
            11,
            1,
            25,
            BlastpGappedAlignmentMode::Restricted,
        );

        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_stop, 8);
        assert_eq!(result.subject_start, 6_995);
        assert_eq!(result.subject_stop, 7_003);
        assert!(result.score > 0);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4708-4711
    // ```c
    // gap_align->score = score_right + score_left;
    // return status;
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:523-550
    // ```c
    // if (esp->op_type[index] == eGapAlignSub) {
    //     sum += factor*matrix[*query & kResidueMask][*subject];
    // } else if (esp->op_type[index] == eGapAlignDel) {
    //     sum -= gap_open + gap_extend * esp->num[index];
    // } else if (esp->op_type[index] == eGapAlignIns) {
    //     sum -= gap_open + gap_extend * esp->num[index];
    // }
    // ```
    #[test]
    fn test_traceback_score_matches_rescored_edit_script_wssv_537() {
        let query_fasta = fs::read_to_string("tests/fasta/WSSV.faa").expect("query fasta");
        let subject_fasta = fs::read_to_string("tests/fasta/PajaWSV.faa").expect("subject fasta");
        let query = mask_query_sequence(&fasta_sequence_by_id(&query_fasta, "YP_009220537.1"));
        let subject = mask_subject_sequence(&fasta_sequence_by_id(&subject_fasta, "BDT63528.1"));

        let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176 / 32.0).unwrap();
        let query_composition = read_aa_composition(&query);
        let subject_composition = read_aa_composition(&subject);
        let adjusted = blast_adjust_scores(
            &matrix_info,
            &query_composition,
            query.len() as i32,
            &subject_composition,
            subject.len() as i32,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
            0,
        )
        .unwrap()
        .expect("adjusted matrix");

        let aligned = blast_gapped_alignment_with_traceback(
            &query,
            &subject,
            8,
            37,
            ScoringMatrix::Blosum62,
            Some(&adjusted.adjusted_matrix),
            11 * 32,
            32,
            100_000,
        )
        .expect("traceback alignment");
        let rescored = score_edit_script_local_alignment(
            &query,
            &subject,
            aligned.query_start,
            aligned.subject_start,
            &aligned.edit_script,
            BlastpScoreMatrix::Adjusted(&adjusted.adjusted_matrix),
            11 * 32,
            32,
        );

        assert_eq!(aligned.score, rescored);
    }
}
