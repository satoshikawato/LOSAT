//! BLASTP gapped alignment with traceback.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c

use crate::common::GapEditOp;
#[cfg(test)]
use crate::common::Hit;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::AdjustedProteinMatrix;
use crate::utils::matrix::{
    blosum62_ncbistdaa_score_row, blosum62_score_ncbistdaa_direct, protein_score, standard_score,
};

#[cfg(test)]
use super::alignment::build_alignment_view_with_matrix;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:363-371
const SCRIPT_SUB: u8 = 3;
const SCRIPT_GAP_IN_A: u8 = 0;
const SCRIPT_GAP_IN_B: u8 = 6;
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:57-60
const GAP_MININT: i32 = i32::MIN / 2;
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_util.h:364
const FENCE_SENTRY: u8 = 201;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4155-4160
const MAX_SUBJECT_OFFSET: usize = 90_000;
const MAX_TOTAL_GAPS: usize = 3_000;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:963-966
const RESTRICT_SIZE: usize = 10;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4273-4278
// EXPERIMENT (LOSAT_X_BPGAP): `s_RestrictedGappedAlign` on separate row
// arrays, eight cells at a time with AVX2. A child module, so that it calls the reference
// function and its helpers without changing them.
#[path = "x_restricted_gap.rs"]
mod x_restricted_gap;
// No NCBI counterpart: the summary line `main` prints at exit; it changes no value NCBI computes.
pub use x_restricted_gap::print_summary as x_restricted_gap_print_summary;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:399-414
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1748-1770
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1902-1945
#[derive(Clone, Copy)]
enum BlastpScoreMatrix<'a> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:737-957
    Blosum62,
    Standard(ScoringMatrix),
    Adjusted(&'a AdjustedProteinMatrix),
}

#[inline(always)]
fn blastp_standard_score_impl<const BLOSUM62: bool>(
    matrix: ScoringMatrix,
    query_residue: u8,
    subject_residue: u8,
) -> i32 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:547-549
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:779-790,802-807
    if BLOSUM62 {
        blosum62_score_ncbistdaa_direct(query_residue, subject_residue)
    } else {
        protein_score(matrix, query_residue, subject_residue)
    }
}

impl BlastpScoreMatrix<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:737-957
    #[inline]
    fn standard(matrix: ScoringMatrix) -> Self {
        if matrix == ScoringMatrix::Blosum62 {
            Self::Blosum62
        } else {
            Self::Standard(matrix)
        }
    }

    #[inline]
    fn score(self, query_residue: u8, subject_residue: u8) -> i32 {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:830-831
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3410-3426
        match self {
            Self::Blosum62 => blosum62_score_ncbistdaa_direct(query_residue, subject_residue),
            Self::Standard(matrix) => standard_score(matrix, query_residue, subject_residue),
            Self::Adjusted(matrix) => matrix.score(query_residue, subject_residue),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543-578
    #[inline(always)]
    fn is_blosum62(self) -> bool {
        matches!(self, Self::Blosum62)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:836-849
#[inline(always)]
fn blastp_score_row_for_query<const BLOSUM62: bool>(
    score_matrix: BlastpScoreMatrix<'_>,
    query_residue: u8,
) -> Option<&[i32]> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:836-849
    //
    // Kappa redo supplies an adjusted square matrix; it is still a normal
    // row-addressed protein score matrix in NCBI, so Rust keeps the same row
    // lookup outside the inner DP loop.
    if BLOSUM62 {
        Some(blosum62_ncbistdaa_score_row(query_residue as usize))
    } else {
        match score_matrix {
            BlastpScoreMatrix::Adjusted(matrix) => Some(&matrix.scores[query_residue as usize]),
            _ => None,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:859
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1099
#[inline(always)]
fn blastp_score_from_row<const BLOSUM62: bool>(
    score_matrix: BlastpScoreMatrix<'_>,
    matrix_row: Option<&[i32]>,
    query_residue: u8,
    subject_residue: u8,
) -> i32 {
    if BLOSUM62 {
        matrix_row.expect("BLOSUM62 score row is always present")[subject_residue as usize]
    } else if let Some(matrix_row) = matrix_row {
        matrix_row[subject_residue as usize]
    } else {
        score_matrix.score(query_residue, subject_residue)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:523-598
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

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:56-62
#[derive(Clone, Copy)]
struct BlastGapDp {
    best: i32,
    best_gap: i32,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:69-97
pub(crate) struct GapAlignScratch {
    // No NCBI counterpart: scratch rows for the SIMD kernel (NCBI keeps its DP rows in
    // gap_align->dp_mem); it does not change any value NCBI computes.
    // EXPERIMENT: SIMD X-drop kernel state.
    fast: crate::utils::xdrop_simd::XdropScratch,
    dp_mem: Vec<BlastGapDp>,
    dp_mem_alloc: usize,
    trace_rows: Vec<Vec<u8>>,
    trace_offsets: Vec<usize>,
    trace_rows_used: usize,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:689-726
    trace_ops_reversed: Vec<(u8, u32)>,
    // No NCBI counterpart: the row arrays of LOSAT_X_BPGAP (NCBI keeps its rows in
    // gap_align->dp_mem); empty unless the switch runs.
    x_restricted: x_restricted_gap::XRestrictedScratch,
}

impl GapAlignScratch {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:313-319
    pub(crate) fn new() -> Self {
        let dp_mem_alloc = 1000usize;
        Self {
            fast: crate::utils::xdrop_simd::XdropScratch::for_search(),
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:689-726
            trace_ops_reversed: Vec::new(),
            x_restricted: x_restricted_gap::XRestrictedScratch::new(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4077-4081
#[derive(Debug, Clone)]
pub(crate) struct BlastpGapAlignResult {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:85-88
    // Raw traceback state includes stop=-1 after a left fence at start0.
    pub query_start: i32,
    pub query_stop: i32,
    pub subject_start: i32,
    pub subject_stop: i32,
    pub score: i32,
    pub edit_script: Vec<GapEditOp>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4074-4091
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4209-4319
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BlastpScoreOnlyGappedAlignment {
    pub query_start: i32,
    pub query_stop: i32,
    pub subject_start: i32,
    pub subject_stop: i32,
    pub score: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4209-4319
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BlastpGappedAlignmentMode {
    Exact,
    Restricted,
}

#[inline(always)]
fn blastp_gapped_start_score<const BLOSUM62: bool>(
    matrix: ScoringMatrix,
    query_residue: u8,
    subject_residue: u8,
) -> i32 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3410-3426
    if BLOSUM62 {
        blosum62_score_ncbistdaa_direct(query_residue, subject_residue)
    } else {
        protein_score(matrix, query_residue, subject_residue)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:802-821
// `n` is NCBI's `N`, the subject letters of this extension (`len2`).
fn gap_dp_reserve_initial(
    dp_mem: &mut Vec<BlastGapDp>,
    dp_mem_alloc: &mut usize,
    num_extra_cells: usize,
    n: usize,
) {
    if num_extra_cells > *dp_mem_alloc {
        *dp_mem_alloc = (num_extra_cells + 100).max(2 * *dp_mem_alloc);
    }
    gap_dp_ensure_cells(dp_mem, *dp_mem_alloc, n);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:923-931
fn gap_dp_reserve_band(
    dp_mem: &mut Vec<BlastGapDp>,
    dp_mem_alloc: &mut usize,
    last_b_index: usize,
    num_extra_cells: usize,
    n: usize,
) {
    if last_b_index + num_extra_cells + 3 >= *dp_mem_alloc {
        *dp_mem_alloc = (last_b_index + num_extra_cells + 100).max(2 * *dp_mem_alloc);
        // `n` is fixed within one extension, so `dp_mem` already holds as many cells as the
        // smaller of the old `dp_mem_alloc` and one more than `n`, and grows only with the count.
        gap_dp_ensure_cells(dp_mem, *dp_mem_alloc, n);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:945-958
// The DP writes a cell before it reads it and never addresses a cell past `N`, so of NCBI's
// `malloc`/`realloc` block of `dp_mem_alloc` cells only the first (the smaller of `dp_mem_alloc`
// and one more than `N`) are ever touched; NCBI leaves the rest of the block untouched. `dp_mem`
// holds those cells and grows when a later extension has a longer `N`; `dp_mem_alloc` keeps NCBI's
// count, which decides the bounds below.
// Filling every cell of a huge `-xdrop_gap` block took seconds and gigabytes
// for an output that does not change (TN-4).
fn gap_dp_ensure_cells(dp_mem: &mut Vec<BlastGapDp>, dp_mem_alloc: usize, n: usize) {
    let cells = dp_mem_alloc.min(n.saturating_add(1));
    if dp_mem.len() < cells {
        dp_mem.resize(
            cells,
            BlastGapDp {
                best: GAP_MININT,
                best_gap: GAP_MININT,
            },
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:128-139
fn gap_reset_traceback_state(scratch: &mut GapAlignScratch) {
    scratch.trace_rows_used = 0;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:70-115
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
        row
    } else {
        let row = Vec::with_capacity(row_capacity_slack);
        trace_rows.push(row);
        trace_rows.last_mut().unwrap()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:560-665
#[inline(always)]
fn gap_trace_row_write(row: &mut Vec<u8>, row_idx: usize, script: u8) {
    if row_idx == row.len() {
        row.push(script);
    } else if row_idx < row.len() {
        row[row_idx] = script;
    } else {
        row.resize(row_idx, SCRIPT_GAP_IN_A);
        row.push(script);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3401-3403
pub(crate) fn blastp_get_start_for_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    q_length: usize,
    s_start: usize,
    s_length: usize,
    matrix: ScoringMatrix,
) -> usize {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3401-3438
    if matrix == ScoringMatrix::Blosum62 {
        blastp_get_start_for_gapped_alignment_impl::<true>(
            query, subject, q_start, q_length, s_start, s_length, matrix,
        )
    } else {
        blastp_get_start_for_gapped_alignment_impl::<false>(
            query, subject, q_start, q_length, s_start, s_length, matrix,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3401-3403
fn blastp_get_start_for_gapped_alignment_impl<const BLOSUM62: bool>(
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
        score += blastp_gapped_start_score::<BLOSUM62>(
            matrix,
            query[q_start + offset],
            subject[s_start + offset],
        );
    }
    let mut max_score = score;
    let mut max_offset = q_start + HSP_MAX_WINDOW - 1;
    let hsp_end = q_start + q_length.min(s_length);

    for index in (q_start + HSP_MAX_WINDOW)..hsp_end {
        let subject_index = s_start + (index - q_start);
        score -= blastp_gapped_start_score::<BLOSUM62>(
            matrix,
            query[index - HSP_MAX_WINDOW],
            subject[subject_index - HSP_MAX_WINDOW],
        );
        score +=
            blastp_gapped_start_score::<BLOSUM62>(matrix, query[index], subject[subject_index]);
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign_priv.h:120
/// `HSP_MAX_WINDOW`, the letters of a window of `BlastGetStartForGappedAlignment`.
pub(crate) const GAPPED_START_WINDOW: usize = 11;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3405-3413
/// The `HSP_MAX_WINDOW` letters of `sequence` from `start` that the first window of
/// `BlastGetStartForGappedAlignment` reads, a letter past the end of `sequence` read as the
/// sentinel byte (NULLB, 0). NCBI's buffers hold a NULLB right after each context of the query and
/// after the subject, whose score is `BLAST_SCORE_MIN` (blast_stat.c:1559-1579): once a window holds
/// it, its sum is below 0 whatever letters follow, as it is when they read as NULLBs.
#[inline]
pub(crate) fn gapped_start_window(sequence: &[u8], start: usize) -> [u8; GAPPED_START_WINDOW] {
    let mut window = [0u8; GAPPED_START_WINDOW];
    for (index, letter) in window.iter_mut().enumerate() {
        *letter = sequence
            .get(start.saturating_add(index))
            .copied()
            .unwrap_or(0);
    }
    window
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1054,1083 and
// c++/src/algo/blast/core/blast_gapalign.c:3394,3401,3405-3413,3418-3420,3434-3437
/// `BlastGetStartForGappedAlignment` of a one-hit HSP of negative `Int4` length (`-window_size 0`
/// with a low `-threshold`: the one-hit extension's best part ends before the last reset of its
/// running sum), from the bytes of its first window (`query_window`, `subject_window`). NCBI reads
/// the length as `Uint4`, above `HSP_MAX_WINDOW`, and the sum of `q_start` and the smaller of
/// `q_length` and `s_length` wraps back to `q_start` plus the length, so only the first window is
/// scored: the gapped start is `q_start + 5` when its sum is above 0, else `q_start`. TBLASTN and
/// BLASTX take it from here; BLASTP keeps its own copy inline in its gapped stage
/// (`blastp::blast_engine::blastp_get_start_for_gapped_alignment_int4_length`, out of line it cost
/// BLASTP's default search 0.6 %).
#[cold]
#[inline(never)]
pub(crate) fn gapped_start_int4(
    query_window: &[u8; GAPPED_START_WINDOW],
    subject_window: &[u8; GAPPED_START_WINDOW],
    q_start: usize,
    matrix: ScoringMatrix,
) -> usize {
    let score: i32 = query_window
        .iter()
        .zip(subject_window)
        .map(|(&q, &s)| standard_score(matrix, q, s))
        .sum();
    if score > 0 {
        q_start + GAPPED_START_WINDOW - 1 - GAPPED_START_WINDOW / 2
    } else {
        q_start
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:340-377
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2481-2534
fn prelim_edit_ops_to_gap_edit_script(
    rev_prelim_tback: Vec<GapEditOp>,
    fwd_prelim_tback: Vec<GapEditOp>,
) -> Vec<GapEditOp> {
    let mut esp = rev_prelim_tback;
    esp.reserve(fwd_prelim_tback.len());
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4681-4707
fn prune_terminal_gap_ops(
    edit_script: &mut Vec<GapEditOp>,
    gap_open: i32,
    gap_extend: i32,
    score_left: &mut i32,
    score_right: &mut i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4687-4690,4703-4706
    query_start: &mut i32,
    query_stop: &mut i32,
    subject_start: &mut i32,
    subject_stop: &mut i32,
) {
    let mut leading_gap_ops = 0usize;
    while leading_gap_ops < edit_script.len() {
        match edit_script[leading_gap_ops] {
            GapEditOp::Sub(_) => break,
            GapEditOp::Del(n) => {
                *score_left += gap_open + gap_extend * n as i32;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4687-4688
                *subject_start += n as i32;
            }
            GapEditOp::Ins(n) => {
                *score_left += gap_open + gap_extend * n as i32;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4689-4690
                *query_start += n as i32;
            }
        }
        leading_gap_ops += 1;
    }
    if leading_gap_ops > 0 {
        edit_script.drain(0..leading_gap_ops);
    }

    while let Some(last) = edit_script.last().copied() {
        match last {
            GapEditOp::Sub(_) => break,
            GapEditOp::Del(n) => {
                *score_right += gap_open + gap_extend * n as i32;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4703-4704
                *subject_stop -= n as i32;
            }
            GapEditOp::Ins(n) => {
                *score_right += gap_open + gap_extend * n as i32;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4705-4706
                *query_stop -= n as i32;
            }
        }
        edit_script.pop();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:745-818
pub(crate) fn stats_from_edit_ops_protein(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    edit_ops: &[GapEditOp],
    matrix: ScoringMatrix,
) -> (usize, usize, usize, usize, usize) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:779-790,802-807
    if matrix == ScoringMatrix::Blosum62 {
        stats_from_edit_ops_protein_impl::<true>(q_seq, s_seq, q_start, s_start, edit_ops, matrix)
    } else {
        stats_from_edit_ops_protein_impl::<false>(q_seq, s_seq, q_start, s_start, edit_ops, matrix)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-596:
//     adjusted_subject, hsp, score_options, &align_length, sbp);
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:767-811:
// Reuse the same edit-script walk after SEG-masked traceback.
pub(crate) fn protein_identities_from_edit_ops(
    query_nomask: &[u8],
    subject: &[u8],
    query_start: usize,
    subject_start: usize,
    edit_ops: &[GapEditOp],
    matrix: ScoringMatrix,
) -> usize {
    stats_from_edit_ops_protein(
        query_nomask,
        subject,
        query_start,
        subject_start,
        edit_ops,
        matrix,
    )
    .0
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:767-811
fn stats_from_edit_ops_protein_impl<const BLOSUM62: bool>(
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
                        if blastp_standard_score_impl::<BLOSUM62>(matrix, q, s) > 0 {
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:839-849
#[inline]
fn blastp_query_residue<const REVERSE: bool>(
    query: &[u8],
    query_base: usize,
    query_length: usize,
    a_index: usize,
) -> u8 {
    let absolute_index = if REVERSE {
        query_base + query_length - a_index
    } else {
        query_base + a_index
    };
    debug_assert!(absolute_index < query.len());
    // SAFETY: ALIGN_EX and Blast_SemiGappedAlign iterate a_index inside the
    // NCBI A sequence span; the debug assertion mirrors that precondition.
    unsafe { *query.get_unchecked(absolute_index) }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:852-866
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:155-161
#[inline]
fn blastp_subject_residue<const REVERSE: bool>(
    subject: &[u8],
    subject_base: usize,
    subject_length: usize,
    b_index: usize,
) -> u8 {
    let absolute_index = if REVERSE {
        if b_index >= subject_length {
            return 0;
        } else {
            subject_base + subject_length - b_index - 1
        }
    } else {
        subject_base + b_index + 1
    };
    if absolute_index >= subject.len() {
        return 0;
    }
    // SAFETY: the explicit absolute_index < subject.len() check preserves the
    // prior sentinel behavior while avoiding a second bounds check in the DP
    // cell hot loop.
    unsafe { *subject.get_unchecked(absolute_index) }
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:841-866
    // No NCBI counterpart: the SIMD kernel (`utils/xdrop_simd.rs`) runs this X-drop; the scalar
    // port below answers the calls it declines.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:862-912
    // Dispatch point: the SIMD kernel runs this score-only X-drop
    // recurrence of Blast_SemiGappedAlign
    // (same max/gap updates, same X-drop test, same best-cell choice) in vector lanes and returns
    // (a_offset, b_offset, score).
    // If the problem does not fit the 8/16/32-bit lanes (or the target has no vector kernel) it
    // returns None and the scalar port below runs.
    if crate::utils::xdrop_simd::KERNEL_COMPILED {
        use crate::utils::xdrop_simd::{xdrop_align, ColSeq, RowSeq};
        let q = if reverse {
            RowSeq::Bytes {
                data: query,
                base: (query_base + len1) as isize,
                step: -1,
            }
        } else {
            RowSeq::Bytes {
                data: query,
                base: query_base as isize,
                step: 1,
            }
        };
        let s = if reverse {
            ColSeq {
                data: subject,
                base: (subject_base + len2) as isize - 1,
                step: -1,
                zero_from: len2,
            }
        } else {
            ColSeq {
                data: subject,
                base: subject_base as isize + 1,
                step: 1,
                zero_from: usize::MAX,
            }
        };
        let fast = with_x_scores(score_matrix, |scores| {
            xdrop_align(
                &q,
                &s,
                scores,
                len1,
                len2,
                gap_open,
                gap_extend,
                x_drop,
                false,
                false,
                &mut scratch.fast,
            )
        });
        if let Some(r) = fast {
            crate::utils::xstats::add(&crate::utils::xstats::DP_FAST_CALLS, 1);
            return (r.a_offset, r.b_offset, r.score);
        }
        crate::utils::xstats::add(&crate::utils::xstats::DP_FALLBACK_CALLS, 1);
    }
    align_ex_protein_score_only_scalar(
        query,
        query_base,
        subject,
        subject_base,
        len1,
        len2,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop,
        reverse,
        scratch,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:841-843,866
// The tables hold the same BLOSUM62 scores that matrix[A[a]][B[b]] reads (built from
// blosum62_score_ncbistdaa_direct); they only change how the row is looked up.
// EXPERIMENT: lookup tables for the SIMD kernel.
fn x_blosum62_tables() -> &'static crate::utils::xdrop_simd::StaticTables {
    use std::sync::OnceLock;
    static T: OnceLock<Box<crate::utils::xdrop_simd::StaticTables>> = OnceLock::new();
    T.get_or_init(|| {
        let f = |q: u8, s: u8| blosum62_score_ncbistdaa_direct(q, s);
        crate::utils::xdrop_simd::build_tables(&f, 28).expect("BLOSUM62 fits the 16-bit tables")
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543,578
// The tables hold the scores that matrix[A[a]][B[b]] reads for a matrix other than BLOSUM62
// (`standard_score`, built as `x_blosum62_tables` is); they only change how the row is
// looked up, and give the kernel the matrix's largest score (its 8-bit lanes need it).
/// The lookup tables of the SIMD kernel for `matrix`, built at the first call for the matrix;
/// `None` if a score is outside the 16-bit tables (the kernel then builds each row as before).
fn x_standard_tables(
    matrix: ScoringMatrix,
) -> Option<&'static crate::utils::xdrop_simd::StaticTables> {
    use std::sync::OnceLock;
    static T: [OnceLock<Option<Box<crate::utils::xdrop_simd::StaticTables>>>;
        ScoringMatrix::ALL.len()] = [const { OnceLock::new() }; ScoringMatrix::ALL.len()];
    T[matrix as usize]
        .get_or_init(|| {
            let f = move |q: u8, s: u8| standard_score(matrix, q, s);
            crate::utils::xdrop_simd::build_tables(&f, 28)
        })
        .as_deref()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:839-850
// Selects the score source for the SIMD kernel: the same matrix rows (BLOSUM62, composition-
// adjusted, or standard) that the scalar port reads.
#[inline]
fn with_x_scores<R>(
    score_matrix: BlastpScoreMatrix<'_>,
    run: impl FnOnce(&crate::utils::xdrop_simd::Scores<'_>) -> R,
) -> R {
    use crate::utils::xdrop_simd::Scores;
    match score_matrix {
        BlastpScoreMatrix::Blosum62 => run(&Scores::Tables {
            t: x_blosum62_tables(),
            n: 28,
        }),
        BlastpScoreMatrix::Adjusted(m) => run(&Scores::Rows28(&m.scores)),
        BlastpScoreMatrix::Standard(matrix) => {
            let f = move |q: u8, s: u8| standard_score(matrix, q, s);
            run(&match x_standard_tables(matrix) {
                Some(t) => Scores::Tables { t, n: 28 },
                None => Scores::Func { f: &f, n: 28 },
            })
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:862-912
// The scalar port of the Blast_SemiGappedAlign score-only recurrence, in its own function: the
// dispatch above calls it for the calls the SIMD kernel does not take.
fn align_ex_protein_score_only_scalar(
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
    match (score_matrix.is_blosum62(), reverse) {
        (true, true) => align_ex_protein_score_only_impl::<true, true>(
            query,
            query_base,
            subject,
            subject_base,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
        ),
        (true, false) => align_ex_protein_score_only_impl::<true, false>(
            query,
            query_base,
            subject,
            subject_base,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
        ),
        (false, true) => align_ex_protein_score_only_impl::<false, true>(
            query,
            query_base,
            subject,
            subject_base,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
        ),
        (false, false) => align_ex_protein_score_only_impl::<false, false>(
            query,
            query_base,
            subject,
            subject_base,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
        ),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:841-866
fn align_ex_protein_score_only_impl<const BLOSUM62: bool, const REVERSE: bool>(
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
    scratch: &mut GapAlignScratch,
) -> (usize, usize, i32) {
    if len1 == 0 || len2 == 0 {
        return (0, 0, 0);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:737-957
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
        len2,
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
        let qc = blastp_query_residue::<REVERSE>(query, query_base, len1, a_index);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:841-849
        let matrix_row = blastp_score_row_for_query::<BLOSUM62>(score_matrix, qc);
        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;
        // No NCBI counterpart: row and cell counters (xstats feature); they do not change any value
        // NCBI computes.
        crate::utils::xstats::add(&crate::utils::xstats::DP_SO_ROWS, 1);
        crate::utils::xstats::add(
            &crate::utils::xstats::DP_SO_CELLS,
            (b_size - first_b_index) as u64,
        );

        #[expect(
            clippy::mut_range_bound,
            reason = "NCBI's loop reads its start once and moves first_b_index past a first column that fails X-drop (blast_gapalign.c:862,886-887)"
        )]
        for b_index in first_b_index..b_size {
            let sc = blastp_subject_residue::<REVERSE>(subject, subject_base, len2, b_index);
            let score_gap_col = scratch.dp_mem[b_index].best_gap;
            let next_score = scratch.dp_mem[b_index].best
                + blastp_score_from_row::<BLOSUM62>(score_matrix, matrix_row, qc, sc);

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
            len2,
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:990-1285
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1080-1115
    match (score_matrix.is_blosum62(), reverse) {
        (true, true) => restricted_gapped_align_protein_impl::<true, true>(
            query,
            query_base,
            subject,
            subject_base,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            dp_mem,
            dp_mem_alloc,
        ),
        (true, false) => restricted_gapped_align_protein_impl::<true, false>(
            query,
            query_base,
            subject,
            subject_base,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            dp_mem,
            dp_mem_alloc,
        ),
        (false, true) => restricted_gapped_align_protein_impl::<false, true>(
            query,
            query_base,
            subject,
            subject_base,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            dp_mem,
            dp_mem_alloc,
        ),
        (false, false) => restricted_gapped_align_protein_impl::<false, false>(
            query,
            query_base,
            subject,
            subject_base,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            dp_mem,
            dp_mem_alloc,
        ),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1080-1180
fn restricted_gapped_align_protein_impl<const BLOSUM62: bool, const REVERSE: bool>(
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
    gap_dp_reserve_initial(dp_mem, dp_mem_alloc, num_extra_cells, len2);

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
        let matrix_row_char = blastp_query_residue::<REVERSE>(query, query_base, len1, a_index);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1080-1087
        let matrix_row = blastp_score_row_for_query::<BLOSUM62>(score_matrix, matrix_row_char);
        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;

        if a_index % RESTRICT_SIZE != 0 {
            #[expect(
                clippy::mut_range_bound,
                reason = "NCBI's loop reads its start once and moves first_b_index past a first column that fails X-drop (blast_gapalign.c:1112,1124-1125)"
            )]
            for b_index in first_b_index..b_size {
                let sc = blastp_subject_residue::<REVERSE>(subject, subject_base, len2, b_index);
                let next_score = dp_mem[b_index].best
                    + blastp_score_from_row::<BLOSUM62>(
                        score_matrix,
                        matrix_row,
                        matrix_row_char,
                        sc,
                    );

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
            #[expect(
                clippy::mut_range_bound,
                reason = "NCBI's loop reads its start once and moves first_b_index past a first column that fails X-drop (blast_gapalign.c:1177,1192-1193)"
            )]
            for b_index in first_b_index..b_size {
                let sc = blastp_subject_residue::<REVERSE>(subject, subject_base, len2, b_index);
                let next_score = dp_mem[b_index].best
                    + blastp_score_from_row::<BLOSUM62>(
                        score_matrix,
                        matrix_row,
                        matrix_row_char,
                        sc,
                    );

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

        gap_dp_reserve_band(dp_mem, dp_mem_alloc, last_b_index, num_extra_cells, len2);
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:364-733
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
    fence_hit: &mut bool,
    read_end_sentinel: bool,
) -> (usize, usize, i32, Vec<GapEditOp>) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543-578
    // No NCBI counterpart: the SIMD kernel (`utils/xdrop_simd.rs`) runs this X-drop; the scalar
    // port below answers the calls it declines.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:563-635
    // Dispatch point: the SIMD kernel runs this traceback recurrence of
    // ALIGN_EX in vector lanes
    // (same per-cell script bits, same fence test, same X-drop test) and returns offsets, score,
    // edit script and fence_hit.
    // A call that does not fit the lanes, has fence_hit already set, or is on a target without a
    // vector kernel runs the scalar port below.
    if crate::utils::xdrop_simd::KERNEL_COMPILED && !*fence_hit {
        use crate::utils::xdrop_simd::{
            xdrop_align, ColSeq, RowSeq, SCRIPT_GAP_IN_A as X_GAP_A, SCRIPT_GAP_IN_B as X_GAP_B,
        };
        let q = if reverse {
            RowSeq::Bytes {
                data: q_seq,
                base: len1 as isize,
                step: -1,
            }
        } else {
            RowSeq::Bytes {
                data: q_seq,
                base: 0,
                step: 1,
            }
        };
        let s = if reverse {
            ColSeq {
                data: s_seq,
                base: len2 as isize - 1,
                step: -1,
                zero_from: len2,
            }
        } else {
            ColSeq {
                data: s_seq,
                base: 1,
                step: 1,
                zero_from: if read_end_sentinel { usize::MAX } else { len2 },
            }
        };
        let fast = with_x_scores(score_matrix, |scores| {
            xdrop_align(
                &q,
                &s,
                scores,
                len1,
                len2,
                gap_open,
                gap_extend,
                x_drop,
                true,
                true,
                &mut scratch.fast,
            )
        });
        // the chunks of the traceback rows past the scratch's own go back to the pool (`ops` holds
        // the script; `XTrace`)
        scratch.fast.x_trace_release();
        if let Some(r) = fast {
            let ops: Vec<GapEditOp> = scratch
                .fast
                .ops
                .iter()
                .map(|&(op, n)| match op {
                    x if x == X_GAP_A => GapEditOp::Del(n),
                    x if x == X_GAP_B => GapEditOp::Ins(n),
                    _ => GapEditOp::Sub(n),
                })
                .collect();
            if r.fence_hit {
                *fence_hit = true;
            }
            crate::utils::xstats::add(&crate::utils::xstats::DP_FAST_CALLS, 1);
            return (r.a_offset, r.b_offset, r.score, ops);
        }
        crate::utils::xstats::add(&crate::utils::xstats::DP_FALLBACK_CALLS, 1);
    }
    align_ex_protein_scalar(
        q_seq,
        s_seq,
        len1,
        len2,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop,
        reverse,
        scratch,
        fence_hit,
        read_end_sentinel,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:563-636
// The scalar port of ALIGN_EX, in its own function: the dispatch above calls it for the calls the
// SIMD kernel does not take.
fn align_ex_protein_scalar(
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
    fence_hit: &mut bool,
    read_end_sentinel: bool,
) -> (usize, usize, i32, Vec<GapEditOp>) {
    match (score_matrix.is_blosum62(), reverse) {
        (true, true) => align_ex_protein_impl::<true, true>(
            q_seq,
            s_seq,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
            fence_hit,
            read_end_sentinel,
        ),
        (true, false) => align_ex_protein_impl::<true, false>(
            q_seq,
            s_seq,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
            fence_hit,
            read_end_sentinel,
        ),
        (false, true) => align_ex_protein_impl::<false, true>(
            q_seq,
            s_seq,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
            fence_hit,
            read_end_sentinel,
        ),
        (false, false) => align_ex_protein_impl::<false, false>(
            q_seq,
            s_seq,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            scratch,
            fence_hit,
            read_end_sentinel,
        ),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543-578
fn align_ex_protein_impl<const BLOSUM62: bool, const REVERSE: bool>(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    score_matrix: BlastpScoreMatrix<'_>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut GapAlignScratch,
    fence_hit: &mut bool,
    read_end_sentinel: bool,
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
        len2,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:470,481-489
    // NCBI's `malloc`ed chunk is touched only up to `N`; the row's capacity is
    // capped there likewise (TN-4).
    row0.reserve(num_extra_cells.min(len2.saturating_add(2)));
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:514-519,669
        // A row holds at most `N + 2 - first_b_index` cells (`b_size <= N + 1`);
        // NCBI's chunk is not touched past them, so the capacity stops there (TN-4).
        let row_capacity = if gap_extend > 0 {
            b_size.saturating_sub(first_b_index) + num_extra_cells
        } else {
            len2.saturating_add(3).saturating_sub(first_b_index)
        }
        .min(len2.saturating_add(2).saturating_sub(first_b_index));
        let edit_script_row = gap_alloc_trace_row(
            &mut scratch.trace_rows,
            &mut scratch.trace_offsets,
            &mut scratch.trace_rows_used,
            row_capacity,
            orig_b_index,
        );
        let qc = blastp_query_residue::<REVERSE>(q_seq, 0, len1, a_index);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:542-551
        let matrix_row = blastp_score_row_for_query::<BLOSUM62>(score_matrix, qc);

        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:560-575
        //
        // Track the C loop variable `b_index` after the inner loop finishes.
        // `ALIGN_EX` uses that exact post-loop value in
        // the larger of `b_index` and `b_size`, minus `orig_b_index`, plus 1.
        let mut row_end_b_index = b_size;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:563-635
        // Capture the row's starting band independently of the X-drop cursor.
        // The slice borrow ends after this ascending scan, before DP reserve.
        let band_start = first_b_index;
        let dp_band = &mut scratch.dp_mem[band_start..b_size];
        // No NCBI counterpart: row and cell counters (xstats feature); they do not change any value
        // NCBI computes.
        crate::utils::xstats::add(&crate::utils::xstats::DP_TB_ROWS, 1);
        crate::utils::xstats::add(&crate::utils::xstats::DP_TB_CELLS, dp_band.len() as u64);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:513-518,531-540
        // The allocator has cleared this row. Reserve the actual main band
        // before borrowing spare capacity; reserve is relative to Vec length.
        debug_assert!(edit_script_row.is_empty());
        edit_script_row.reserve(dp_band.len());
        let trace_band = &mut edit_script_row.spare_capacity_mut()[..dp_band.len()];
        for (band_index, (cell, trace_cell)) in
            dp_band.iter_mut().zip(trace_band.iter_mut()).enumerate()
        {
            let b_index = band_start + band_index;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:563-569
            // Column `b_index` reads `B[b_index + 1]`, so the last column (`b_size == N + 1`)
            // reads the byte after the extension's letters. A translated subject window
            // keeps its sentinel (the right fence, or the end of a full translation) there
            // (`read_end_sentinel`); other callers read 0.
            let sc = if b_index < len2 || read_end_sentinel {
                blastp_subject_residue::<REVERSE>(s_seq, 0, len2, b_index)
            } else {
                0
            };
            let score_gap_col = cell.best_gap;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:569-575
            if sc == FENCE_SENTRY {
                row_end_b_index = b_index;
                *fence_hit = true;
                break;
            }
            let next_score =
                cell.best + blastp_score_from_row::<BLOSUM62>(score_matrix, matrix_row, qc, sc);

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
                    cell.best = GAP_MININT;
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
                    cell.best_gap = open_gap_col;
                } else {
                    cell.best_gap = score_gap_col_ext;
                    script += SCRIPT_EXTEND_GAP_B;
                }

                let open_gap_row = score_val - gap_open_extend;
                if score_gap_row < open_gap_row {
                    score_gap_row = open_gap_row;
                } else {
                    script += SCRIPT_EXTEND_GAP_A;
                }
                cell.best = score_val;
            }

            score_val = next_score;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:634-635
            trace_cell.write(script);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:563-575,634-639
        //     if (matrix_index == FENCE_SENTRY) { ... break; }
        // SAFETY: the row began empty and capacity covers the captured band.
        // Every visited non-fence cell, including X-drop failures, wrote one
        // byte. The only early break sets row_end_b_index before that write.
        // Thus exactly this prefix is initialized; the spare borrow has ended.
        unsafe { edit_script_row.set_len(row_end_b_index - band_start) };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:638-639
        if first_b_index == b_size || *fence_hit {
            break;
        }

        gap_dp_reserve_band(
            &mut scratch.dp_mem,
            &mut scratch.dp_mem_alloc,
            last_b_index,
            num_extra_cells,
            len2,
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
                let row_idx = b_size.saturating_sub(orig_b_index);
                gap_trace_row_write(edit_script_row, row_idx, SCRIPT_GAP_IN_A);
                b_size += 1;
            }
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:666-669
        //
        // The C code only exposes the traceback cells written for this row.
        // Keep the Rust row length to the same written span so zero-filled
        // tail capacity cannot be read back as `SCRIPT_GAP_IN_A`.
        let used_cells = row_end_b_index
            .max(b_size)
            .saturating_sub(orig_b_index)
            .saturating_add(1);
        if edit_script_row.len() < used_cells {
            edit_script_row.resize(used_cells, SCRIPT_GAP_IN_A);
        } else {
            edit_script_row.truncate(used_cells);
        }

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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:689-726
    scratch.trace_ops_reversed.clear();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:686-687
    if *fence_hit {
        return (a_offset, b_offset, best_score, Vec::new());
    }

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
            b_index = b_index.saturating_sub(1);
        } else if op_type == SCRIPT_GAP_IN_B {
            a_index = a_index.saturating_sub(1);
        } else {
            a_index = a_index.saturating_sub(1);
            b_index = b_index.saturating_sub(1);
        }
        push_trace_op(&mut scratch.trace_ops_reversed, op_type);
    }

    let mut edit_ops = Vec::with_capacity(scratch.trace_ops_reversed.len());
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:689-726
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2494-2534
    //
    // `ALIGN_EX` must therefore return the preliminary traceback in the
    // `GapPrelimEditBlockAdd` insertion order for both directions. Only
    // `Blast_PrelimEditBlockToGapEditScript` reverses the forward half.
    for (op_type, count) in scratch.trace_ops_reversed.drain(..) {
        match op_type {
            x if x == SCRIPT_SUB => edit_ops.push(GapEditOp::Sub(count)),
            x if x == SCRIPT_GAP_IN_A => edit_ops.push(GapEditOp::Del(count)),
            x if x == SCRIPT_GAP_IN_B => edit_ops.push(GapEditOp::Ins(count)),
            _ => edit_ops.push(GapEditOp::Sub(count)),
        }
    }

    (a_offset, b_offset, best_score, edit_ops)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4549-4712
#[cfg(test)]
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1914-1944
    let mut scratch = GapAlignScratch::new();
    blast_gapped_alignment_with_traceback_with_scratch(
        query,
        subject,
        q_start,
        s_start,
        matrix,
        adjusted_matrix,
        gap_open,
        gap_extend,
        x_drop,
        &mut scratch,
        None,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1914-1944
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4549-4554,4620-4643:
// Int2 BLAST_GappedAlignmentWithTraceback(..., Boolean * fence_hit);
// if ((! (fence_hit && *fence_hit)) && ...) { /* right extension */ }
pub(crate) fn blast_gapped_alignment_with_traceback_with_scratch(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    matrix: ScoringMatrix,
    adjusted_matrix: Option<&AdjustedProteinMatrix>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut GapAlignScratch,
    fence_hit_out: Option<&mut bool>,
) -> Option<BlastpGapAlignResult> {
    gapped_alignment_with_traceback(
        query,
        subject,
        None,
        q_start,
        s_start,
        matrix,
        adjusted_matrix,
        gap_open,
        gap_extend,
        x_drop,
        scratch,
        fence_hit_out,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:463-464,509-512
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1230-1235
/// The traceback of a translated subject window (TBLASTN): `subject` ends with the window's
/// sentinel, which is not one of the `subject_length` letters (NCBI's subject length is the
/// end of the translated window). The right extension stops before the sentinel, and its last
/// column reads it as NCBI's does (a right fence, or the end of a full translation).
pub(crate) fn blast_gapped_alignment_with_traceback_in_window(
    query: &[u8],
    subject: &[u8],
    subject_length: usize,
    q_start: usize,
    s_start: usize,
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut GapAlignScratch,
    fence_hit_out: Option<&mut bool>,
) -> Option<BlastpGapAlignResult> {
    gapped_alignment_with_traceback(
        query,
        subject,
        Some(subject_length),
        q_start,
        s_start,
        matrix,
        None,
        gap_open,
        gap_extend,
        x_drop,
        scratch,
        fence_hit_out,
    )
}

#[allow(clippy::too_many_arguments)]
fn gapped_alignment_with_traceback(
    query: &[u8],
    subject: &[u8],
    window_length: Option<usize>,
    q_start: usize,
    s_start: usize,
    matrix: ScoringMatrix,
    adjusted_matrix: Option<&AdjustedProteinMatrix>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut GapAlignScratch,
    fence_hit_out: Option<&mut bool>,
) -> Option<BlastpGapAlignResult> {
    let subject_length = window_length.unwrap_or(subject.len());
    if q_start >= query.len() || s_start >= subject_length {
        return None;
    }

    let score_matrix = adjusted_matrix.map_or(BlastpScoreMatrix::standard(matrix), |matrix| {
        BlastpScoreMatrix::Adjusted(matrix)
    });

    let mut fence_hit = false;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4599-4611
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
        scratch,
        &mut fence_hit,
        false,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4614-4615
    let raw_q_start = i32::try_from(q_start).expect("NCBI query offset fits Int4");
    let raw_s_start = i32::try_from(s_start).expect("NCBI subject offset fits Int4");
    let mut query_start =
        raw_q_start - i32::try_from(left_q).expect("NCBI query extent fits Int4") + 1;
    let mut subject_start =
        raw_s_start - i32::try_from(left_s).expect("NCBI subject extent fits Int4") + 1;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4629-4642
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4620-4647
    let (mut query_stop, mut subject_stop, mut score_right, right_edit_ops) =
        if !fence_hit && q_start < query.len() && s_start < subject_length {
            let (right_q, right_s, score_right, right_edit_ops) = align_ex_protein(
                &query[q_start..],
                &subject[s_start..],
                query.len().saturating_sub(q_start + 1),
                subject_length.saturating_sub(s_start + 1),
                score_matrix,
                gap_open,
                gap_extend,
                x_drop,
                false,
                scratch,
                &mut fence_hit,
                window_length.is_some(),
            );
            (
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4646-4647:
                raw_q_start + i32::try_from(right_q).expect("NCBI query extent fits Int4") + 1,
                raw_s_start + i32::try_from(right_s).expect("NCBI subject extent fits Int4") + 1,
                score_right,
                right_edit_ops,
            )
        } else {
            (
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4650-4653
                raw_q_start - 1,
                raw_s_start - 1,
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4620-4643:
    // if ((! (fence_hit && *fence_hit)) && ...) { /* right extension */ }
    // The caller must see the fence flag before deciding whether to retry.
    if let Some(out) = fence_hit_out {
        *out = fence_hit;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4714-4718
    // A successful NCBI traceback may have zero score, equal endpoints and
    // an empty script after a failed composition adjustment. The caller owns
    // cutoff/reap decisions; preserve that payload here.
    let score = score_left + score_right;

    Some(BlastpGapAlignResult {
        query_start,
        query_stop,
        subject_start,
        subject_stop,
        score,
        edit_script,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4209-4319,4333-4361
pub(crate) fn blastp_score_only_gapped_alignment_with_scratch(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    mode: BlastpGappedAlignmentMode,
    scratch: &mut GapAlignScratch,
) -> BlastpScoreOnlyGappedAlignment {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4251-4258
    assert!(
        q_start < query.len() && s_start < subject.len(),
        "NCBI BLAST requires preliminary gapped starts to lie within query/subject bounds"
    );

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4247-4258
    let q_length = q_start + 1;
    let s_length = s_start + 1;
    let (subject_offset, subject_length, subject_shift) =
        adjust_subject_range(s_length, subject.len(), q_length, query.len());
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4257-4258
    assert!(
        subject_shift < subject.len() && subject_length > 0,
        "NCBI BLAST requires AdjustSubjectRange to yield a non-empty subject window"
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4209-4319
    let score_matrix = BlastpScoreMatrix::standard(matrix);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4266-4281
    let (private_q_start, private_s_start, score_left) = match mode {
        BlastpGappedAlignmentMode::Exact => align_ex_protein_score_only(
            query,
            0,
            subject,
            subject_shift,
            q_length,
            subject_offset,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            true,
            scratch,
        ),
        // EXPERIMENT (LOSAT_X_BPGAP): the reference function unless the switch is on
        // (x_restricted_gap.rs).
        BlastpGappedAlignmentMode::Restricted => {
            x_restricted_gap::restricted_gapped_align_protein_x(
                query,
                0,
                subject,
                subject_shift,
                q_length,
                subject_offset,
                score_matrix,
                gap_open,
                gap_extend,
                x_drop,
                true,
                scratch,
            )
        }
    };
    let query_hit_start = i32::try_from(q_length - private_q_start)
        .expect("NCBI BLAST preliminary query start must fit in Int4");
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4290-4291
    let subject_hit_start = i32::try_from(subject_offset - private_s_start + subject_shift)
        .expect("NCBI BLAST preliminary subject start must fit in Int4");

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4295-4334
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
                    score_matrix,
                    gap_open,
                    gap_extend,
                    x_drop,
                    false,
                    scratch,
                ),
                // EXPERIMENT (LOSAT_X_BPGAP): as for the left extension.
                BlastpGappedAlignmentMode::Restricted => {
                    x_restricted_gap::restricted_gapped_align_protein_x(
                        query,
                        q_start,
                        subject,
                        s_start,
                        query.len() - q_length,
                        subject_length - subject_offset,
                        score_matrix,
                        gap_open,
                        gap_extend,
                        x_drop,
                        false,
                        scratch,
                    )
                }
            };
            (
                i32::try_from(q_start + private_q_stop + 1)
                    .expect("NCBI BLAST preliminary query stop must fit in Int4"),
                i32::try_from(s_start + private_s_stop + 1)
                    .expect("NCBI BLAST preliminary subject stop must fit in Int4"),
                score_right,
            )
        } else {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4342-4344
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

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:69-80
#[cfg(test)]
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
    let mut scratch = GapAlignScratch::new();
    blastp_score_only_gapped_alignment_with_scratch(
        query,
        subject,
        q_start,
        s_start,
        matrix,
        gap_open,
        gap_extend,
        x_drop,
        mode,
        &mut scratch,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:647-662
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:234-250
#[cfg(test)]
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
        .saturating_sub(view.num_ident.saturating_add(view.gap_letters)) as i64;
    hit.identity = if view.alignment_len > 0 {
        (view.num_ident as f64 * 100.0) / view.alignment_len as f64
    } else {
        0.0
    };
    true
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:439-477
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:440-642
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:399-414
#[cfg(test)]
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:547-549
    if matrix == ScoringMatrix::Blosum62 {
        reevaluate_blastp_hit_with_traceback_impl::<true>(
            hit,
            query,
            query_nomask,
            subject,
            matrix,
            gap_open,
            gap_extend,
            cutoff_score,
        )
    } else {
        reevaluate_blastp_hit_with_traceback_impl::<false>(
            hit,
            query,
            query_nomask,
            subject,
            matrix,
            gap_open,
            gap_extend,
            cutoff_score,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:523-551
#[cfg(test)]
fn reevaluate_blastp_hit_with_traceback_impl<const BLOSUM62: bool>(
    hit: &mut Hit,
    query: &[u8],
    query_nomask: &[u8],
    subject: &[u8],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    cutoff_score: i32,
) -> bool {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:509-510
    const BLAST_HITS_RESIDUE_MASK: u8 = 0x0f;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:787-846
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:555-564
    //
    // For blastp, `BlastScoringOptionsNew` leaves `reward` at zero, and `BlastScoringParametersNew`
    // copies that zero into `score_params->reward`, so the addition of `ext` times
    // `score_params->reward` to `score` in the end extension below adds 0 and the score stored
    // before it is NCBI's.
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
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:547-549
                    sum += blastp_standard_score_impl::<BLOSUM62>(
                        matrix,
                        query[query_pos] & BLAST_HITS_RESIDUE_MASK,
                        subject[subject_pos],
                    );
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
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:562-572
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:614-637
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
        // Adding `ext` times `score_params->reward` to `score` adds 0 for blastp (reward is 0, above).

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
        // Adding `ext` times `score_params->reward` to `score` adds 0 for blastp (reward is 0, above).
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:649-660
    !refresh_blastp_hit_statistics(hit, query_nomask, subject, matrix)
}

#[cfg(test)]
#[path = "gapalign_tests.rs"]
mod tests;
