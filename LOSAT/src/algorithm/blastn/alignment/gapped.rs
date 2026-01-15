//! Gapped extension algorithms for BLASTN
//!
//! This module implements gapped extension algorithms with affine gap penalties,
//! including both heuristic (bidirectional) and one-directional extensions.

use super::greedy::{greedy_align_one_direction, greedy_align_one_direction_ex};
use crate::common::GapEditOp;

/// Alignment statistics propagated alongside DP scores
#[derive(Clone, Copy, Default)]
pub struct AlnStats {
    matches: u32,
    mismatches: u32,
    gap_opens: u32,
    gap_letters: u32, // Total gap characters (for alignment length calculation)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign_priv.h:120
const HSP_MAX_WINDOW: usize = 11;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:82-83 (COMPRESSION_RATIO)
const COMPRESSION_RATIO: usize = 4;

/// BLASTNA alphabet size.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1060-1068
const BLASTNA_SIZE: usize = 16;

/// Map BLASTNA codes to NCBI4NA bitmasks (degeneracy).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:56-74
const BLASTNA_TO_NCBI4NA: [u8; BLASTNA_SIZE] = [
    1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0,
];

type BlastnaMatrix = [i32; BLASTNA_SIZE * BLASTNA_SIZE];

/// Round to nearest integer (half away from zero).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:437-441 (BLAST_Nint)
#[inline]
fn blast_nint(x: f64) -> i32 {
    if x >= 0.0 {
        (x + 0.5) as i32
    } else {
        (x - 0.5) as i32
    }
}

/// Build the BLASTNA scoring matrix from reward/penalty.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1052-1127 (BlastScoreBlkNuclMatrixCreate)
pub fn build_blastna_matrix(reward: i32, penalty: i32) -> BlastnaMatrix {
    let mut matrix = [0i32; BLASTNA_SIZE * BLASTNA_SIZE];
    let mut degeneracy = [0i32; BLASTNA_SIZE];

    // NCBI reference: blast_stat.c:1084-1099
    for index1 in 0..4 {
        degeneracy[index1] = 1;
    }
    for index1 in 4..BLASTNA_SIZE {
        let mut degen = 0;
        for index2 in 0..4 {
            if (BLASTNA_TO_NCBI4NA[index1] & BLASTNA_TO_NCBI4NA[index2]) != 0 {
                degen += 1;
            }
        }
        degeneracy[index1] = degen;
    }

    // NCBI reference: blast_stat.c:1102-1119
    for index1 in 0..BLASTNA_SIZE {
        for index2 in index1..BLASTNA_SIZE {
            let idx = index1 * BLASTNA_SIZE + index2;
            let val = if (BLASTNA_TO_NCBI4NA[index1] & BLASTNA_TO_NCBI4NA[index2]) != 0 {
                let degen = degeneracy[index2];
                let raw = ((degen - 1) * penalty + reward) as f64 / (degen as f64);
                blast_nint(raw)
            } else {
                penalty
            };
            matrix[idx] = val;
            if index1 != index2 {
                matrix[index2 * BLASTNA_SIZE + index1] = val;
            }
        }
    }

    // NCBI reference: blast_stat.c:1122-1127 (gap sentinel row/column)
    for index1 in 0..BLASTNA_SIZE {
        matrix[(BLASTNA_SIZE - 1) * BLASTNA_SIZE + index1] = i32::MIN / 2;
        matrix[index1 * BLASTNA_SIZE + (BLASTNA_SIZE - 1)] = i32::MIN / 2;
    }

    matrix
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:56-62 (BlastGapDP)
#[derive(Clone, Copy)]
struct BlastGapDP {
    best: i32,
    best_gap: i32,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:857-860 (MININT initialization)
const GAP_MININT: i32 = i32::MIN / 2;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_util.h:358-364
// ```c
// #define FENCE_SENTRY 201
// ```
const FENCE_SENTRY: u8 = 201;

/// Scratch memory mirroring NCBI's BlastGapAlignStruct.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80
pub struct GapAlignScratch {
    dp_mem: Vec<BlastGapDP>,
    dp_mem_alloc: usize,
    trace_rows: Vec<Vec<u8>>,
    trace_offsets: Vec<usize>,
    trace_rows_used: usize,
}

impl GapAlignScratch {
    pub fn new() -> Self {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:313-319 (BLAST_GapAlignStructNew)
        let dp_mem_alloc = 1000;
        let dp_mem = vec![
            BlastGapDP { best: GAP_MININT, best_gap: GAP_MININT };
            dp_mem_alloc
        ];

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:448-450 (edit_script/edit_start_offset arrays)
        let trace_rows = Vec::with_capacity(100);
        let trace_offsets = Vec::with_capacity(100);

        Self {
            dp_mem,
            dp_mem_alloc,
            trace_rows,
            trace_offsets,
            trace_rows_used: 0,
        }
    }
}

fn gap_dp_reserve_initial(
    score_array: &mut Vec<BlastGapDP>,
    dp_mem_alloc: &mut usize,
    num_extra_cells: usize,
) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:797-808 (Blast_SemiGappedAlign dp_mem realloc)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:462-468 (ALIGN_EX dp_mem realloc)
    if num_extra_cells + 3 >= *dp_mem_alloc {
        let new_alloc = (num_extra_cells + 100).max(*dp_mem_alloc * 2);
        *dp_mem_alloc = new_alloc;
        score_array.resize(new_alloc, BlastGapDP { best: GAP_MININT, best_gap: GAP_MININT });
    }
}

fn gap_dp_reserve_band(
    score_array: &mut Vec<BlastGapDP>,
    dp_mem_alloc: &mut usize,
    last_b_index: usize,
    num_extra_cells: usize,
) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:923-931 (Blast_SemiGappedAlign dynamic realloc)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:641-648 (ALIGN_EX dynamic realloc)
    if last_b_index + num_extra_cells + 3 >= *dp_mem_alloc {
        let new_alloc = (last_b_index + num_extra_cells + 100).max(*dp_mem_alloc * 2);
        *dp_mem_alloc = new_alloc;
        score_array.resize(new_alloc, BlastGapDP { best: GAP_MININT, best_gap: GAP_MININT });
    }
}

fn gap_reset_traceback_state(trace_rows_used: &mut usize) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:128-139 (s_GapPurgeState)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:438-439 (ALIGN_EX calls s_GapPurgeState)
    *trace_rows_used = 0;
}

fn gap_alloc_trace_row<'a>(
    trace_rows: &'a mut Vec<Vec<u8>>,
    trace_offsets: &mut Vec<usize>,
    trace_rows_used: &mut usize,
    row_capacity: usize,
    start_offset: usize,
) -> &'a mut Vec<u8> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:70-115 (s_GapGetState)
    let row_index = *trace_rows_used;
    *trace_rows_used += 1;

    if row_index < trace_offsets.len() {
        trace_offsets[row_index] = start_offset;
    } else {
        trace_offsets.push(start_offset);
    }

    if row_index < trace_rows.len() {
        let row = &mut trace_rows[row_index];
        row.clear();
        if row_capacity > 0 {
            row.resize(row_capacity, 0);
        }
        row
    } else {
        let mut row: Vec<u8> = Vec::new();
        if row_capacity > 0 {
            row.resize(row_capacity, 0);
        }
        trace_rows.push(row);
        trace_rows.last_mut().unwrap()
    }
}

/// Select a gapped-start seed within an HSP using a sliding window score.
///
/// NCBI reference: blast_gapalign.c:3248-3305 BlastGetOffsetsForGappedAlignment
pub fn blast_get_offsets_for_gapped_alignment(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    score_matrix: &BlastnaMatrix,
) -> Option<(usize, usize)> {
    // NCBI reference: blast_gapalign.c:3254-3257
    let q_length = q_end.saturating_sub(q_start);
    let s_length = s_end.saturating_sub(s_start);

    if q_length == 0 || s_length == 0 {
        return None;
    }

    if q_end > q_seq.len() || s_end > s_seq.len() {
        return None;
    }

    // NCBI reference: blast_gapalign.c:3259-3263
    if q_length <= HSP_MAX_WINDOW {
        let mid = q_start + q_length / 2;
        return Some((mid, s_start + q_length / 2));
    }

    // NCBI reference: blast_gapalign.c:3265-3278
    let mut score: i32 = 0;
    let mut q_idx = q_start;
    let mut s_idx = s_start;
    for _ in 0..HSP_MAX_WINDOW {
        let q_code = q_seq[q_idx] as usize;
        let s_code = s_seq[s_idx] as usize;
        // NCBI reference: blast_gapalign.c:3270-3273 (matrix-based scoring)
        score += score_matrix[q_code * BLASTNA_SIZE + s_code];
        q_idx += 1;
        s_idx += 1;
    }

    let mut max_score = score;
    let mut max_offset = q_start + HSP_MAX_WINDOW - 1;

    // NCBI reference: blast_gapalign.c:3278
    let hsp_end = q_start + q_length.min(s_length);

    // NCBI reference: blast_gapalign.c:3279-3292
    let mut index = q_start + HSP_MAX_WINDOW;
    while index < hsp_end {
        let q_leave = q_seq[q_idx - HSP_MAX_WINDOW] as usize;
        let s_leave = s_seq[s_idx - HSP_MAX_WINDOW] as usize;
        let q_enter = q_seq[q_idx] as usize;
        let s_enter = s_seq[s_idx] as usize;
        // NCBI reference: blast_gapalign.c:3280-3286 (matrix-based window update)
        score -= score_matrix[q_leave * BLASTNA_SIZE + s_leave];
        score += score_matrix[q_enter * BLASTNA_SIZE + s_enter];

        if score > max_score {
            max_score = score;
            max_offset = index;
        }

        q_idx += 1;
        s_idx += 1;
        index += 1;
    }

    // NCBI reference: blast_gapalign.c:3294-3299
    if max_score > 0 {
        let q_retval = max_offset;
        let s_retval = (max_offset - q_start) + s_start;
        return Some((q_retval, s_retval));
    }

    // NCBI reference: blast_gapalign.c:3300-3317
    score = 0;
    q_idx = q_end - HSP_MAX_WINDOW;
    s_idx = s_end - HSP_MAX_WINDOW;
    for _ in (q_end - HSP_MAX_WINDOW)..q_end {
        let q_code = q_seq[q_idx] as usize;
        let s_code = s_seq[s_idx] as usize;
        // NCBI reference: blast_gapalign.c:3304-3310 (matrix-based scoring)
        score += score_matrix[q_code * BLASTNA_SIZE + s_code];
        q_idx += 1;
        s_idx += 1;
    }

    if score > 0 {
        let q_retval = q_end - HSP_MAX_WINDOW / 2;
        let s_retval = s_end - HSP_MAX_WINDOW / 2;
        return Some((q_retval, s_retval));
    }

    None
}

/// Refine a gapped-start seed for nucleotide alignments based on identity runs.
///
/// NCBI reference: blast_gapalign.c:3323-3389 BlastGetStartForGappedAlignmentNucl
pub fn blast_get_start_for_gapped_alignment_nucl(
    q_seq: &[u8],
    s_seq: &[u8],
    q_offset: usize,
    q_end: usize,
    s_offset: usize,
    s_end: usize,
    q_gapped_start: usize,
    s_gapped_start: usize,
) -> (usize, usize) {
    // NCBI reference: blast_gapalign.c:3326-3332
    let mut hsp_max_ident_run: i32 = 10;
    let offset = (s_gapped_start - s_offset).min(q_gapped_start - q_offset);

    // NCBI reference: blast_gapalign.c:3334-3349
    let mut score: i32 = -1;
    let mut q_idx = q_gapped_start;
    let mut s_idx = s_gapped_start;
    let q_len_limit = q_end;

    while q_idx < q_len_limit && q_idx < q_seq.len() && s_idx < s_seq.len() && q_seq[q_idx] == s_seq[s_idx] {
        score += 1;
        if score > hsp_max_ident_run {
            return (q_gapped_start, s_gapped_start);
        }
        q_idx += 1;
        s_idx += 1;
    }

    q_idx = q_gapped_start;
    s_idx = s_gapped_start;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3344-3348
    // ```c
    // q = query + q_start;
    // s = subject + s_start;
    // while ((q-query >= 0) && (*q-- == *s--)) {
    //    score++;
    //    if (score > hspMaxIdentRun) return;
    // }
    // ```
    loop {
        if q_idx >= q_seq.len() || s_idx >= s_seq.len() {
            break;
        }
        if q_seq[q_idx] != s_seq[s_idx] {
            break;
        }
        score += 1;
        if score > hsp_max_ident_run {
            return (q_gapped_start, s_gapped_start);
        }
        if q_idx == 0 || s_idx == 0 {
            break;
        }
        q_idx -= 1;
        s_idx -= 1;
    }

    // NCBI reference: blast_gapalign.c:3350
    hsp_max_ident_run = (hsp_max_ident_run * 3) / 2;

    // NCBI reference: blast_gapalign.c:3352-3357
    let q_start = q_gapped_start - offset;
    let s_start = s_gapped_start - offset;
    let q_len = (s_end - s_start).min(q_end - q_start);

    if q_start + q_len > q_seq.len() || s_start + q_len > s_seq.len() {
        return (q_gapped_start, s_gapped_start);
    }

    // NCBI reference: blast_gapalign.c:3357-3389
    let mut max_score: i32 = 0;
    let mut max_offset = q_start;
    score = 0;
    let mut match_run = false;
    let mut prev_match = false;

    q_idx = q_start;
    s_idx = s_start;
    let mut index = q_start;
    let end = q_start + q_len;
    while index < end {
        match_run = q_seq[q_idx] == s_seq[s_idx];
        if match_run != prev_match {
            prev_match = match_run;
            if match_run {
                score = 1;
            } else if score > max_score {
                max_score = score;
                max_offset = index - (score as usize / 2);
            }
        } else if match_run {
            score += 1;
            if score > hsp_max_ident_run {
                let max_offset = index - (hsp_max_ident_run as usize / 2);
                let s_new = max_offset + s_start - q_start;
                return (max_offset, s_new);
            }
        }
        q_idx += 1;
        s_idx += 1;
        index += 1;
    }

    if match_run && score > max_score {
        max_score = score;
        max_offset = (q_start + q_len) - (score as usize / 2);
    }

    if max_score > 0 {
        let s_new = max_offset + s_start - q_start;
        return (max_offset, s_new);
    }

    (q_gapped_start, s_gapped_start)
}

pub fn extend_gapped_heuristic(
    q_seq: &[u8],
    s_seq: &[u8],
    s_len: usize,
    qs: usize,
    ss: usize,
    len: usize,
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    use_dp: bool, // If true, use DP-based extension (for blastn task); if false, use greedy (for megablast)
) -> (
    usize,
    usize,
    usize,
    usize,
    i32,
    usize,
    usize,
    usize,
    usize,
    usize,
) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:313-319 (BLAST_GapAlignStructNew)
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    let mut gap_scratch = GapAlignScratch::new();
    extend_gapped_heuristic_with_scratch(
        q_seq,
        s_seq,
        s_len,
        qs,
        ss,
        len,
        reward,
        penalty,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop,
        &mut gap_scratch,
        use_dp,
    )
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct reuse)
pub fn extend_gapped_heuristic_with_scratch(
    q_seq: &[u8],
    s_seq: &[u8],
    s_len: usize,
    qs: usize,
    ss: usize,
    len: usize,
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    gap_scratch: &mut GapAlignScratch,
    use_dp: bool, // If true, use DP-based extension (for blastn task); if false, use greedy (for megablast)
) -> (
    usize,
    usize,
    usize,
    usize,
    i32,
    usize,
    usize,
    usize,
    usize,
    usize,
) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2957-2960 (query/subject lengths)
    let subject_len = if use_dp { s_len } else { s_seq.len() };

    // Bounds validation: ensure seed coordinates are valid
    if qs >= q_seq.len() || ss >= subject_len {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0, 0);
    }

    if use_dp {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2953-3016 (s_BlastDynProgNtGappedAlignment)
        let offset_adjustment = COMPRESSION_RATIO - (ss % COMPRESSION_RATIO);
        let mut q_length = qs + offset_adjustment;
        let mut s_length = ss + offset_adjustment;

        // NCBI reference: blast_gapalign.c:2978-2984 (prevent start past end)
        if q_length > q_seq.len() || s_length > subject_len {
            q_length = q_length.saturating_sub(COMPRESSION_RATIO);
            s_length = s_length.saturating_sub(COMPRESSION_RATIO);
        }

        // NCBI reference: blast_gapalign.c:2986-2993 (left extension, packed subject)
        let (left_q_consumed, left_s_consumed, left_score, left_dp_cells) =
            blast_align_packed_nucl_with_scratch(
                q_seq,
                s_seq,
                q_length,
                s_length,
                0,
                score_matrix,
                gap_open,
                gap_extend,
                x_drop,
                true,
                0,
                gap_scratch,
            );

        // NCBI reference: blast_gapalign.c:2995-3007 (right extension)
        let (right_q_consumed, right_s_consumed, right_score, right_dp_cells) =
            if q_length < q_seq.len() && s_length < subject_len {
                // NCBI reference: blast_gapalign.c:2999-3003 (query+q_length-1, subject+(s_length+3)/4-1)
                let query_offset = q_length as isize - 1;
                let subject_byte_offset = (s_length + 3) / COMPRESSION_RATIO;
                blast_align_packed_nucl_with_scratch(
                    q_seq,
                    s_seq,
                    q_seq.len().saturating_sub(q_length),
                    subject_len.saturating_sub(s_length),
                    query_offset,
                    score_matrix,
                    gap_open,
                    gap_extend,
                    x_drop,
                    false,
                    subject_byte_offset,
                    gap_scratch,
                )
            } else {
                (0, 0, 0, 0)
            };

        let final_q_start = q_length.saturating_sub(left_q_consumed);
        let final_s_start = s_length.saturating_sub(left_s_consumed);
        let final_q_end = if q_length < q_seq.len() && s_length < subject_len {
            q_length + right_q_consumed
        } else {
            q_length
        };
        let final_s_end = if q_length < q_seq.len() && s_length < subject_len {
            s_length + right_s_consumed
        } else {
            s_length
        };

        if std::env::var("LOSAT_DEBUG_COORDS").is_ok() {
            eprintln!("[COORDS] dp_start: qs={}, ss={}", qs, ss);
            eprintln!(
                "[COORDS] dp_left: q_consumed={}, s_consumed={}",
                left_q_consumed, left_s_consumed
            );
            eprintln!(
                "[COORDS] dp_right: q_consumed={}, s_consumed={}",
                right_q_consumed, right_s_consumed
            );
            eprintln!(
                "[COORDS] dp_final: q={}-{}, s={}-{}",
                final_q_start, final_q_end, final_s_start, final_s_end
            );
        }

        return (
            final_q_start,
            final_q_end,
            final_s_start,
            final_s_end,
            left_score + right_score,
            0,
            0,
            0,
            0,
            left_dp_cells + right_dp_cells,
        );
    }

    // Clamp len to available sequence length
    let len = len.min(q_seq.len() - qs).min(subject_len - ss);
    if len == 0 {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0, 0);
    }

    // NCBI BLAST uses different extension algorithms based on task:
    // - megablast: greedy alignment (fast, good for high-identity sequences)
    // - blastn: DP-based alignment (slower, but handles divergent sequences better)
    // This follows NCBI BLAST's approach where blastn uses eDynProgScoreOnly
    // and megablast uses eGreedyScoreOnly.

    let seed_len = len;

    // First, extend to the right from the seed (start point excluded).
    let right_q_start = qs + seed_len;
    let right_s_start = ss + seed_len;
    let right_suffix_q = q_seq.get(right_q_start..).unwrap_or(&[]);
    let right_suffix_s = s_seq.get(right_s_start..).unwrap_or(&[]);

    // Greedy extension for megablast task (faster for high-identity sequences)
    let (
        right_q_consumed,
        right_s_consumed,
        right_score,
        right_matches,
        right_mismatches,
        right_gaps,
        right_gap_letters,
    ) = greedy_align_one_direction(
        right_suffix_q,
        right_suffix_s,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop,
    );

    // Then extend to the left from the seed
    let left_q_len = qs;
    let left_s_len = ss;
    let (
        left_q_consumed,
        left_s_consumed,
        left_score,
        left_matches,
        left_mismatches,
        left_gaps,
        left_gap_letters,
    ) = greedy_align_one_direction_ex(
        q_seq,
        s_seq,
        left_q_len,
        left_s_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop,
        true,  // reverse = true for left extension
    );

    let mut seed_score = 0;
    let mut seed_matches = 0;
    let mut seed_mismatches = 0;
    for k in 0..seed_len {
        if q_seq[qs + k] == s_seq[ss + k] {
            seed_score += reward;
            seed_matches += 1;
        } else {
            seed_score += penalty;
            seed_mismatches += 1;
        }
    }

    let total_score = left_score + seed_score + right_score;
    let total_matches = left_matches + seed_matches + right_matches;
    let total_mismatches = left_mismatches + seed_mismatches + right_mismatches;
    let total_gaps = left_gaps + right_gaps;
    let total_gap_letters = left_gap_letters + right_gap_letters;

    let final_q_start = qs.saturating_sub(left_q_consumed);
    let final_q_end = qs + seed_len + right_q_consumed;
    let final_s_start = ss.saturating_sub(left_s_consumed);
    let final_s_end = ss + seed_len + right_s_consumed;

    if std::env::var("LOSAT_DEBUG_COORDS").is_ok() {
        eprintln!("[COORDS] seed: qs={}, ss={}, len={}", qs, ss, seed_len);
        eprintln!(
            "[COORDS] left: q_consumed={}, s_consumed={}",
            left_q_consumed, left_s_consumed
        );
        eprintln!(
            "[COORDS] right: q_consumed={}, s_consumed={}",
            right_q_consumed, right_s_consumed
        );
        eprintln!(
            "[COORDS] final: q={}-{}, s={}-{}",
            final_q_start, final_q_end, final_s_start, final_s_end
        );
    }

    (
        final_q_start,
        final_q_end,
        final_s_start,
        final_s_end,
        total_score,
        total_matches,
        total_mismatches,
        total_gaps,
        total_gap_letters,
        0,
    )
}

/// Extend alignment in one direction using banded Smith-Waterman with affine gap penalties.
///
/// This implements a proper row-by-row DP algorithm with counts propagation for accurate statistics.
/// Instead of storing full traceback matrices, we propagate alignment statistics (matches, mismatches,
/// gap_opens) alongside the DP scores, keeping only 2 rows at a time for O(band_size) memory.
///
/// - Row i corresponds to query position i
/// - Band index k represents diagonal offset: j = i + (k - W) where W is half-bandwidth
/// - Three matrices track different alignment states:
///   - M[i,j]: best score ending with a match/mismatch at (i,j)
///   - Ix[i,j]: best score ending with a gap in subject (deletion from query)
///   - Iy[i,j]: best score ending with a gap in query (insertion in query)
///
/// Returns: (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters)
/// Extend alignment in one direction using NCBI-style semi-gapped DP with affine gap penalties.
///
/// This is a faithful transpilation of NCBI BLAST's Blast_SemiGappedAlign (score-only mode).
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:735-962
///
/// Key algorithm features:
/// - X-drop based dynamic window that expands/contracts based on score
/// - Tracks best score and position
/// - Minimal memory per cell (8 bytes vs 40 bytes with stats)
/// - No hard-coded extension limits (controlled by X-drop termination)
///
/// Returns: (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters, dp_cells)
pub fn extend_gapped_one_direction(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> (usize, usize, i32, usize, usize, usize, usize, usize) {
    // NCBI reference: blast_gapalign.c:735-962 Blast_SemiGappedAlign
    // This implements score-only mode (score_only = TRUE)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:313-319 (BLAST_GapAlignStructNew)
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    let mut gap_scratch = GapAlignScratch::new();
    extend_gapped_one_direction_with_scratch(
        q_seq,
        s_seq,
        reward,
        penalty,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop,
        &mut gap_scratch,
    )
}

fn extend_gapped_one_direction_with_scratch(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    gap_scratch: &mut GapAlignScratch,
) -> (usize, usize, i32, usize, usize, usize, usize, usize) {
    // NCBI reference: blast_gapalign.c:735-962 Blast_SemiGappedAlign
    // This implements score-only mode (score_only = TRUE)

    // NCBI reference: blast_gapalign.c:780-782
    // gap_open = score_params->gap_open;
    // gap_extend = score_params->gap_extend;
    // gap_open_extend = gap_open + gap_extend;
    let gap_open_extend = gap_open + gap_extend;

    // NCBI reference: blast_gapalign.c:783
    // x_dropoff = gap_align->gap_x_dropoff;
    let mut x_dropoff = x_drop;

    // NCBI reference: blast_gapalign.c:785-786
    // if (x_dropoff < gap_open_extend) x_dropoff = gap_open_extend;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    let m = q_seq.len();
    let n = s_seq.len();

    // NCBI reference: blast_gapalign.c:788-789
    // if(N <= 0 || M <= 0) return 0;
    if m == 0 || n == 0 {
        return (0, 0, 0, 0, 0, 0, 0, 0);
    }

    // NCBI reference: blast_gapalign.c:797-800
    // if (gap_extend > 0)
    //     num_extra_cells = x_dropoff / gap_extend + 3;
    // else
    //     num_extra_cells = N + 3;
    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend + 3) as usize
    } else {
        n + 3
    };

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (gap_align->dp_mem)
    let dp_mem_alloc = &mut gap_scratch.dp_mem_alloc;
    let score_array = &mut gap_scratch.dp_mem;

    // NCBI reference: blast_gapalign.c:797-808 (Blast_SemiGappedAlign dp_mem realloc)
    gap_dp_reserve_initial(score_array, dp_mem_alloc, num_extra_cells);

    // NCBI reference: blast_gapalign.c:811-822
    // Initialize row 0: leading gaps in subject
    // score = -gap_open_extend;
    // score_array[0].best = 0;
    // score_array[0].best_gap = -gap_open_extend;
    let mut score = -(gap_open_extend as i32);
    score_array[0].best = 0;
    score_array[0].best_gap = -(gap_open_extend as i32);

    // NCBI reference: blast_gapalign.c:815-822
    // for (i = 1; i <= N; i++) {
    //     if (score < -x_dropoff) break;
    //     score_array[i].best = score;
    //     score_array[i].best_gap = score - gap_open_extend;
    //     score -= gap_extend;
    // }
    let mut b_size = 1usize;
    for i in 1..=n.min(*dp_mem_alloc - 1) {
        if score < -x_dropoff {
            break;
        }
        score_array[i].best = score;
        score_array[i].best_gap = score - (gap_open_extend as i32);
        score -= gap_extend as i32;
        b_size = i + 1;
    }

    // NCBI reference: blast_gapalign.c:827-828
    // b_size = i;
    // best_score = 0;
    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;

    // NCBI reference: blast_gapalign.c:829
    // first_b_index = 0;
    let mut first_b_index = 0usize;

    // DP cell counter for diagnostics
    let mut dp_cells = 0usize;

    // NCBI reference: blast_gapalign.c:835-959
    // Main DP loop - for each query position (row)
    for a_index in 1..=m {
        // Debug: track last few rows
        if std::env::var("LOSAT_DEBUG_COORDS").is_ok() && m > 1000 && a_index >= m.saturating_sub(2) {
            eprintln!("[DP_ROW] a_index={}/{}, first_b={}, b_size={}", a_index, m, first_b_index, b_size);
        }

        // NCBI reference: blast_gapalign.c:839-843
        // matrix_row = matrix[ A[ a_index ] ];
        let qc = q_seq[a_index - 1];

        // NCBI reference: blast_gapalign.c:857-860
        // score = MININT;
        // score_gap_row = MININT;
        // last_b_index = first_b_index;
        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;

        // NCBI reference: blast_gapalign.c:862-912
        // Inner loop - for each subject position in the band
        for b_index in first_b_index..b_size {
            dp_cells += 1;

            // NCBI reference: blast_gapalign.c:862-871 (b_size can reach N+1; no b_index < N guard).
            // NCBI reference: blast_util.c:826 (NULLB sentinel at sequence ends).
            let sc = if b_index < n { s_seq[b_index] } else { 0 };
            let score_gap_col = score_array[b_index].best_gap;
            let row = (qc as usize) * BLASTNA_SIZE;
            // NCBI reference: blast_gapalign.c:864-866 (matrix_row lookup)
            let match_score = score_matrix[row + (sc as usize)];
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:866
            // ```c
            // next_score = score_array[b_index].best + matrix_row[*b_ptr];
            // ```
            let next_score = score_array[b_index].best + match_score;

            // NCBI reference: blast_gapalign.c:868-872
            // if (score < score_gap_col) score = score_gap_col;
            // if (score < score_gap_row) score = score_gap_row;
            if score_val < score_gap_col {
                score_val = score_gap_col;
            }
            if score_val < score_gap_row {
                score_val = score_gap_row;
            }

            // NCBI reference: blast_gapalign.c:874-890
            // X-drop check
            if best_score - score_val > x_dropoff {
                // Failed X-drop
                if b_index == first_b_index {
                    first_b_index += 1;
                } else {
                    score_array[b_index].best = GAP_MININT;
                }
            } else {
                // NCBI reference: blast_gapalign.c:892-908
                last_b_index = b_index;

                // Update best score and position
                if score_val > best_score {
                    best_score = score_val;
                    a_offset = a_index;
                    b_offset = b_index;
                }

                // NCBI reference: blast_gapalign.c:903-907
                // score_gap_row -= gap_extend;
                // score_gap_col -= gap_extend;
                // score_array[b_index].best_gap = MAX(score - gap_open_extend, score_gap_col);
                // score_gap_row = MAX(score - gap_open_extend, score_gap_row);
                score_gap_row -= gap_extend as i32;
                let score_gap_col_ext = score_gap_col - (gap_extend as i32);
                let open_gap_col = score_val - (gap_open_extend as i32);
                score_array[b_index].best_gap = open_gap_col.max(score_gap_col_ext);

                let open_gap_row = score_val - (gap_open_extend as i32);
                score_gap_row = open_gap_row.max(score_gap_row);

                // NCBI reference: blast_gapalign.c:908
                // score_array[b_index].best = score;
                score_array[b_index].best = score_val;
            }

            // NCBI reference: blast_gapalign.c:911
            // score = next_score;
            score_val = next_score;
        }

        // NCBI reference: blast_gapalign.c:918-919
        // if (first_b_index == b_size) break;
        if first_b_index >= b_size {
            break;
        }

        // NCBI reference: blast_gapalign.c:923-931
        // Enlarge window if necessary (dynamic reallocation)
        gap_dp_reserve_band(score_array, dp_mem_alloc, last_b_index, num_extra_cells);

        // NCBI reference: blast_gapalign.c:933-952
        // Band contraction: shorten loop bounds if X-dropoff failed earlier than last row
        if last_b_index < b_size.saturating_sub(1) {
            // NCBI: b_size = last_b_index + 1
            b_size = last_b_index + 1;
        } else {
            // NCBI reference: blast_gapalign.c:946-951
            // Extend window with gaps
            while score_gap_row >= best_score - x_dropoff && b_size <= n && b_size < *dp_mem_alloc {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - (gap_open_extend as i32);
                score_gap_row -= gap_extend as i32;
                b_size += 1;
            }
        }

        // NCBI reference: blast_gapalign.c:954-958
        // Sentinel
        if b_size <= n && b_size < *dp_mem_alloc {
            score_array[b_size].best = GAP_MININT;
            score_array[b_size].best_gap = GAP_MININT;
            b_size += 1;
        }
    }

    // NCBI reference: blast_gapalign.c:961
    // return best_score;
    if best_score <= 0 {
        return (0, 0, 0, 0, 0, 0, 0, dp_cells);
    }

    // Compute statistics from the alignment positions
    // For score-only mode, we estimate stats from score and positions
    // This matches NCBI's approach where stats are computed in traceback

    // Debug: check offset values
    if std::env::var("LOSAT_DEBUG_COORDS").is_ok() {
        eprintln!("[EXT_FWD] m={}, n={}, a_offset={}, b_offset={}, best_score={}", m, n, a_offset, b_offset, best_score);
    }

    // NCBI reference: blast_gapalign.c:893-896
    // a_offset and b_offset represent the loop index where best score was found.
    //
    // CRITICAL INSIGHT: In LOSAT, score_val at loop index b_index=k is the score from
    // the PREVIOUS iteration (next_score at b_index=k-1). This represents the diagonal
    // score to position k-1, not k. So when we set b_offset=k, it's actually 1 higher
    // than the consumed position.
    //
    // Therefore: s_consumed = b_offset (not b_offset + 1) to match the actual position.
    // For diagonal alignments: a_offset=k, b_offset=k, s_consumed=k (matches q_consumed=k)
    let q_consumed = a_offset;
    let s_consumed = b_offset;

    // Estimate matches/mismatches/gaps from score
    // For a gapless alignment: score = matches * reward + mismatches * penalty
    // alignment_len = q_consumed = s_consumed (for gapless)
    // For gapped: this is an approximation
    let alignment_len = q_consumed.max(s_consumed);
    let gap_letters = (q_consumed as i32 - s_consumed as i32).unsigned_abs() as usize;
    let gap_opens = if gap_letters > 0 { 1 } else { 0 };

    // Estimate matches from score (assuming no gaps for simplicity)
    // score = matches * reward + mismatches * penalty
    // matches + mismatches = alignment_len - gap_letters
    let non_gap_len = alignment_len.saturating_sub(gap_letters);
    let matches = if non_gap_len > 0 && reward > 0 {
        let estimated = (best_score + (non_gap_len as i32) * (-penalty)) / (reward - penalty);
        estimated.max(0) as usize
    } else {
        non_gap_len
    };
    let mismatches = non_gap_len.saturating_sub(matches);

    (
        q_consumed,
        s_consumed,
        best_score,
        matches,
        mismatches,
        gap_opens,
        gap_letters,
        dp_cells,
    )
}

/// Score-only gapped extension for packed nucleotide subjects.
///
/// This mirrors NCBI BLAST's s_BlastAlignPackedNucl implementation.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3033-3194
fn blast_align_packed_nucl_with_scratch(
    query: &[u8],
    subject_packed: &[u8],
    query_len: usize,
    subject_len: usize,
    query_offset: isize,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse_sequence: bool,
    subject_byte_offset: usize,
    gap_scratch: &mut GapAlignScratch,
) -> (usize, usize, i32, usize) {
    // NCBI reference: blast_gapalign.c:3063-3074
    let gap_open_extend = gap_open + gap_extend;
    let mut x_dropoff = x_drop;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    if query_len == 0 || subject_len == 0 {
        return (0, 0, 0, 0);
    }

    // NCBI reference: blast_gapalign.c:3082-3085
    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend + 3) as usize
    } else {
        query_len + 3
    };

    let dp_mem_alloc = &mut gap_scratch.dp_mem_alloc;
    let score_array = &mut gap_scratch.dp_mem;

    // NCBI reference: blast_gapalign.c:3087-3093 (dp_mem realloc)
    gap_dp_reserve_initial(score_array, dp_mem_alloc, num_extra_cells);

    // NCBI reference: blast_gapalign.c:3095-3105 (row 0 init)
    let mut score = -(gap_open_extend as i32);
    score_array[0].best = 0;
    score_array[0].best_gap = -(gap_open_extend as i32);

    let mut b_size = 1usize;
    for i in 1..=query_len {
        if score < -x_dropoff {
            break;
        }
        score_array[i].best = score;
        score_array[i].best_gap = score - (gap_open_extend as i32);
        score -= gap_extend as i32;
        b_size = i + 1;
    }

    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;
    let b_increment: isize = if reverse_sequence { -1 } else { 1 };

    let mut dp_cells = 0usize;

    // NCBI reference: blast_gapalign.c:3120-3185
    for a_index in 1..=subject_len {
        // NCBI reference: blast_gapalign.c:3125-3134 (packed subject base)
        let a_base_pair = if reverse_sequence {
            let byte_index = (subject_len - a_index) / COMPRESSION_RATIO;
            let shift = (a_index - 1) % COMPRESSION_RATIO;
            let byte = subject_packed.get(byte_index).copied().unwrap_or(0);
            (byte >> (2 * shift)) & 0x03
        } else {
            let byte_index = subject_byte_offset + (a_index - 1) / COMPRESSION_RATIO;
            let shift = 3 - ((a_index - 1) % COMPRESSION_RATIO);
            let byte = subject_packed.get(byte_index).copied().unwrap_or(0);
            (byte >> (2 * shift)) & 0x03
        };

        let row = (a_base_pair as usize) * BLASTNA_SIZE;

        // NCBI reference: blast_gapalign.c:3136-3140 (b_ptr init)
        let mut b_ptr_index: isize = if reverse_sequence {
            (query_len - first_b_index) as isize
        } else {
            first_b_index as isize
        };

        // NCBI reference: blast_gapalign.c:3141-3144
        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;

        for b_index in first_b_index..b_size {
            b_ptr_index += b_increment;
            let q_idx = query_offset + b_ptr_index;
            let q_code = if q_idx >= 0 {
                query.get(q_idx as usize).copied().unwrap_or(0) as usize
            } else {
                0usize
            };
            let score_gap_col = score_array[b_index].best_gap;
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3150
            // ```c
            // next_score = score_array[b_index].best + matrix_row[*b_ptr];
            // ```
            let next_score = score_array[b_index].best + score_matrix[row + q_code];

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
                    score_array[b_index].best = GAP_MININT;
                }
            } else {
                last_b_index = b_index;
                if score_val > best_score {
                    best_score = score_val;
                    a_offset = a_index;
                    b_offset = b_index;
                }

                score_gap_row -= gap_extend as i32;
                let score_gap_col_ext = score_gap_col - (gap_extend as i32);
                let open_gap_col = score_val - (gap_open_extend as i32);
                score_array[b_index].best_gap = open_gap_col.max(score_gap_col_ext);

                let open_gap_row = score_val - (gap_open_extend as i32);
                score_gap_row = open_gap_row.max(score_gap_row);

                score_array[b_index].best = score_val;
            }

            score_val = next_score;
            dp_cells += 1;
        }

        if first_b_index == b_size {
            break;
        }

        // NCBI reference: blast_gapalign.c:3176-3180
        gap_dp_reserve_band(score_array, dp_mem_alloc, last_b_index, num_extra_cells);

        if last_b_index < b_size.saturating_sub(1) {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - x_dropoff && b_size <= query_len {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - (gap_open_extend as i32);
                score_gap_row -= gap_extend as i32;
                b_size += 1;
            }
        }

        if b_size <= query_len {
            score_array[b_size].best = GAP_MININT;
            score_array[b_size].best_gap = GAP_MININT;
            b_size += 1;
        }
    }

    (b_offset, a_offset, best_score, dp_cells)
}

/// Extended version of extend_gapped_one_direction that supports reverse access.
/// When `reverse` is true, sequences are accessed from the end (for left extension),
/// avoiding the O(n) copy overhead of reversing the sequence.
///
/// This is a faithful transpilation of NCBI BLAST's Blast_SemiGappedAlign (score-only mode)
/// with reverse access support for left extension.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:735-962
///
/// Parameters:
/// - len1, len2: The actual lengths to use (can be less than q_seq.len(), s_seq.len())
/// - reverse: If true, access sequences in reverse order
pub fn extend_gapped_one_direction_ex(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
) -> (usize, usize, i32, usize, usize, usize, usize, usize) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:313-319 (BLAST_GapAlignStructNew)
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    let mut gap_scratch = GapAlignScratch::new();
    extend_gapped_one_direction_ex_with_scratch(
        q_seq,
        s_seq,
        len1,
        len2,
        reward,
        penalty,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop,
        reverse,
        &mut gap_scratch,
    )
}

fn extend_gapped_one_direction_ex_with_scratch(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    gap_scratch: &mut GapAlignScratch,
) -> (usize, usize, i32, usize, usize, usize, usize, usize) {
    // NCBI reference: blast_gapalign.c:780-782
    // gap_open_extend = gap_open + gap_extend;
    let gap_open_extend = gap_open + gap_extend;

    // NCBI reference: blast_gapalign.c:783-786
    // if (x_dropoff < gap_open_extend) x_dropoff = gap_open_extend;
    let mut x_dropoff = x_drop;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    let m = len1;
    let n = len2;

    // NCBI reference: blast_gapalign.c:788-789
    if m == 0 || n == 0 {
        return (0, 0, 0, 0, 0, 0, 0, 0);
    }

    // Helper to get sequence character with optional reverse access
    #[inline(always)]
    fn get_q(q_seq: &[u8], i: usize, len1: usize, reverse: bool) -> u8 {
        if reverse {
            q_seq[len1 - i]
        } else {
            q_seq[i - 1]
        }
    }

    #[inline(always)]
    fn get_s(s_seq: &[u8], j: usize, len2: usize, reverse: bool) -> u8 {
        if j >= len2 {
            // NCBI reference: blast_gapalign.c:862-871 (b_size can reach N+1; b_ptr walks onto NULLB).
            // NCBI reference: blast_util.c:826 (NULLB sentinel at sequence ends).
            return 0;
        }
        if reverse {
            s_seq[len2 - 1 - j]
        } else {
            s_seq[j]
        }
    }

    // NCBI reference: blast_gapalign.c:797-800
    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend + 3) as usize
    } else {
        n + 3
    };

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (gap_align->dp_mem)
    let dp_mem_alloc = &mut gap_scratch.dp_mem_alloc;
    let score_array = &mut gap_scratch.dp_mem;

    // NCBI reference: blast_gapalign.c:797-808 (Blast_SemiGappedAlign dp_mem realloc)
    gap_dp_reserve_initial(score_array, dp_mem_alloc, num_extra_cells);

    // NCBI reference: blast_gapalign.c:811-822
    let mut score = -(gap_open_extend as i32);
    score_array[0].best = 0;
    score_array[0].best_gap = -(gap_open_extend as i32);

    // NCBI reference: blast_gapalign.c:815-822
    let mut b_size = 1usize;
    for i in 1..=n.min(*dp_mem_alloc - 1) {
        if score < -x_dropoff {
            break;
        }
        score_array[i].best = score;
        score_array[i].best_gap = score - (gap_open_extend as i32);
        score -= gap_extend as i32;
        b_size = i + 1;
    }

    // NCBI reference: blast_gapalign.c:827-829
    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;
    let mut dp_cells = 0usize;

    // NCBI reference: blast_gapalign.c:835-959
    for a_index in 1..=m {
        // NCBI reference: blast_gapalign.c:839-843
        let qc = get_q(q_seq, a_index, len1, reverse);

        // NCBI reference: blast_gapalign.c:857-860
        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;

        // NCBI reference: blast_gapalign.c:862-912
        for b_index in first_b_index..b_size {
            dp_cells += 1;

            // NCBI reference: blast_gapalign.c:864-866
            let sc = get_s(s_seq, b_index, len2, reverse);
            let score_gap_col = score_array[b_index].best_gap;
            let row = (qc as usize) * BLASTNA_SIZE;
            // NCBI reference: blast_gapalign.c:864-866 (matrix_row lookup)
            let match_score = score_matrix[row + (sc as usize)];
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:866
            // ```c
            // next_score = score_array[b_index].best + matrix_row[*b_ptr];
            // ```
            let next_score = score_array[b_index].best + match_score;

            // NCBI reference: blast_gapalign.c:868-872
            if score_val < score_gap_col {
                score_val = score_gap_col;
            }
            if score_val < score_gap_row {
                score_val = score_gap_row;
            }

            // NCBI reference: blast_gapalign.c:874-890
            if best_score - score_val > x_dropoff {
                if b_index == first_b_index {
                    first_b_index += 1;
                } else {
                    score_array[b_index].best = GAP_MININT;
                }
            } else {
                // NCBI reference: blast_gapalign.c:892-908
                last_b_index = b_index;

                if score_val > best_score {
                    best_score = score_val;
                    a_offset = a_index;
                    b_offset = b_index;
                }

                // NCBI reference: blast_gapalign.c:903-907
                score_gap_row -= gap_extend as i32;
                let score_gap_col_ext = score_gap_col - (gap_extend as i32);
                let open_gap_col = score_val - (gap_open_extend as i32);
                score_array[b_index].best_gap = open_gap_col.max(score_gap_col_ext);

                let open_gap_row = score_val - (gap_open_extend as i32);
                score_gap_row = open_gap_row.max(score_gap_row);

                // NCBI reference: blast_gapalign.c:908
                score_array[b_index].best = score_val;
            }

            // NCBI reference: blast_gapalign.c:911
            score_val = next_score;
        }

        // NCBI reference: blast_gapalign.c:918-919
        if first_b_index >= b_size {
            break;
        }

        // NCBI reference: blast_gapalign.c:923-931
        gap_dp_reserve_band(score_array, dp_mem_alloc, last_b_index, num_extra_cells);

        // NCBI reference: blast_gapalign.c:933-952
        // Band contraction: shorten loop bounds if X-dropoff failed earlier than last row
        if last_b_index < b_size.saturating_sub(1) {
            // NCBI: b_size = last_b_index + 1
            b_size = last_b_index + 1;
        } else {
            // NCBI reference: blast_gapalign.c:946-951
            while score_gap_row >= best_score - x_dropoff && b_size <= n && b_size < *dp_mem_alloc {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - (gap_open_extend as i32);
                score_gap_row -= gap_extend as i32;
                b_size += 1;
            }
        }

        // NCBI reference: blast_gapalign.c:954-958
        if b_size <= n && b_size < *dp_mem_alloc {
            score_array[b_size].best = GAP_MININT;
            score_array[b_size].best_gap = GAP_MININT;
            b_size += 1;
        }
    }

    // NCBI reference: blast_gapalign.c:961
    if best_score <= 0 {
        return (0, 0, 0, 0, 0, 0, 0, dp_cells);
    }

    // Compute statistics from the alignment positions
    // For score-only mode, we estimate stats from score and positions

    // Debug: check offset values for extension
    if std::env::var("LOSAT_DEBUG_COORDS").is_ok() {
        eprintln!("[EXT_{}] m={}, n={}, a_offset={}, b_offset={}, best_score={}, expected_perfect_score={}",
            if reverse { "REV" } else { "FWD" },
            m, n, a_offset, b_offset, best_score, (m.min(n) * reward as usize) as i32);
    }

    // NCBI reference: blast_gapalign.c:893-896
    // a_offset and b_offset represent the loop index where best score was found.
    //
    // CRITICAL INSIGHT: In LOSAT, score_val at loop index b_index=k is the score from
    // the PREVIOUS iteration (next_score at b_index=k-1). This represents the diagonal
    // score to position k-1, not k. So when we set b_offset=k, it's actually 1 higher
    // than the consumed position.
    //
    // Therefore: s_consumed = b_offset (not b_offset + 1) to match the actual position.
    // For diagonal alignments: a_offset=k, b_offset=k, s_consumed=k (matches q_consumed=k)
    let q_consumed = a_offset;
    let s_consumed = b_offset;

    // Estimate matches/mismatches/gaps from score
    let alignment_len = q_consumed.max(s_consumed);
    let gap_letters = (q_consumed as i32 - s_consumed as i32).unsigned_abs() as usize;
    let gap_opens = if gap_letters > 0 { 1 } else { 0 };

    let non_gap_len = alignment_len.saturating_sub(gap_letters);
    let matches = if non_gap_len > 0 && reward > 0 {
        let estimated = (best_score + (non_gap_len as i32) * (-penalty)) / (reward - penalty);
        estimated.max(0) as usize
    } else {
        non_gap_len
    };
    let mismatches = non_gap_len.saturating_sub(matches);

    (
        q_consumed,
        s_consumed,
        best_score,
        matches,
        mismatches,
        gap_opens,
        gap_letters,
        dp_cells,
    )
}

// =============================================================================
// TRACEBACK-CAPTURING VERSION OF GAPPED EXTENSION
// NCBI Reference: blast_gapalign.c:364-733 (ALIGN_EX function)
// =============================================================================

/// Script operation codes for traceback
/// Reference: blast_gapalign.c:363-371
///
/// ```c
/// enum {
///     SCRIPT_SUB           = eGapAlignSub,     // Substitution
///     SCRIPT_GAP_IN_A      = eGapAlignDel,     // Deletion (gap in query)
///     SCRIPT_GAP_IN_B      = eGapAlignIns,     // Insertion (gap in subject)
///     SCRIPT_OP_MASK       = 0x07,             // Mask for opcode
///     SCRIPT_EXTEND_GAP_A  = 0x10,             // Continue a gap in A
///     SCRIPT_EXTEND_GAP_B  = 0x40              // Continue a gap in B
/// };
/// ```
const SCRIPT_SUB: u8 = 3;           // eGapAlignSub
const SCRIPT_GAP_IN_A: u8 = 0;      // eGapAlignDel (gap in query)
const SCRIPT_GAP_IN_B: u8 = 6;      // eGapAlignIns (gap in subject)
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;

/// Compute alignment statistics from an edit script.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:745-818 (identity)
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1055-1074 (gap counts)
fn stats_from_edit_ops(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    edit_ops: &[GapEditOp],
) -> (usize, usize, usize, usize) {
    let mut matches = 0usize;
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
                    if qi < q_seq.len() && si < s_seq.len() {
                        if q_seq[qi] == s_seq[si] {
                            matches += 1;
                        } else {
                            mismatches += 1;
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

    (matches, mismatches, gap_opens, gap_letters)
}

/// Extend alignment in one direction with TRACEBACK capture.
///
/// This is a faithful transpilation of NCBI BLAST's ALIGN_EX function.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:364-733
///
/// Unlike the score-only version, this function:
/// 1. Stores edit_script for each DP cell (1 byte per cell)
/// 2. After DP, reconstructs the traceback path from best position to origin
/// 3. Returns the edit script as a Vec<GapEditOp> with run-length encoding
///
/// Returns: (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters, edit_ops)
pub fn extend_gapped_one_direction_with_traceback(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> (usize, usize, i32, usize, usize, usize, usize, Vec<GapEditOp>) {
    // NCBI reference: blast_gapalign.c:364-733 ALIGN_EX function
    // This is the traceback-capturing version (score_only = FALSE case)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:313-319 (BLAST_GapAlignStructNew)
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    let mut gap_scratch = GapAlignScratch::new();
    extend_gapped_one_direction_with_traceback_with_scratch(
        q_seq,
        s_seq,
        reward,
        penalty,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop,
        &mut gap_scratch,
    )
}

fn extend_gapped_one_direction_with_traceback_with_scratch(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    gap_scratch: &mut GapAlignScratch,
) -> (usize, usize, i32, usize, usize, usize, usize, Vec<GapEditOp>) {
    // NCBI reference: blast_gapalign.c:364-733 ALIGN_EX function
    // This is the traceback-capturing version (score_only = FALSE case)

    // NCBI reference: blast_gapalign.c:423-426
    // gap_open_extend = gap_open + gap_extend;
    let gap_open_extend = gap_open + gap_extend;

    // NCBI reference: blast_gapalign.c:428-429
    // if (x_dropoff < gap_open_extend) x_dropoff = gap_open_extend;
    let mut x_dropoff = x_drop;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    let m = q_seq.len();
    let n = s_seq.len();

    // NCBI reference: blast_gapalign.c:431-432
    // if(N <= 0 || M <= 0) return 0;
    if m == 0 || n == 0 {
        return (0, 0, 0, 0, 0, 0, 0, Vec::new());
    }

    // NCBI reference: blast_gapalign.c:457-460
    // if (gap_extend > 0)
    //     num_extra_cells = x_dropoff / gap_extend + 3;
    // else
    //     num_extra_cells = N + 3;
    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend + 3) as usize
    } else {
        n + 3
    };

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (gap_align->dp_mem)
    let dp_mem_alloc = &mut gap_scratch.dp_mem_alloc;
    let score_array = &mut gap_scratch.dp_mem;
    let trace_rows = &mut gap_scratch.trace_rows;
    let trace_offsets = &mut gap_scratch.trace_offsets;
    let trace_rows_used = &mut gap_scratch.trace_rows_used;

    // NCBI reference: blast_gapalign.c:462-468 (ALIGN_EX dp_mem realloc)
    gap_dp_reserve_initial(score_array, dp_mem_alloc, num_extra_cells);

    // NCBI reference: blast_gapalign.c:448-450, 472-474 (edit_script/edit_start_offset arrays)
    // NCBI reference: blast_gapalign.c:438-439 (ALIGN_EX calls s_GapPurgeState)
    gap_reset_traceback_state(trace_rows_used);
    let mut edit_script_num_rows = 100usize;

    // NCBI reference: blast_gapalign.c:476-489
    // Initialize row 0 (gaps in subject at start)
    let mut score = -(gap_open_extend as i32);
    score_array[0].best = 0;
    score_array[0].best_gap = -(gap_open_extend as i32);

    // First edit script row for initial gap extension
    let row0 = gap_alloc_trace_row(trace_rows, trace_offsets, trace_rows_used, 0, 0);
    row0.reserve(num_extra_cells);
    row0.push(0); // Position 0

    // NCBI reference: blast_gapalign.c:481-490
    // for (i = 1; i <= N; i++) {
    //     if (score < -x_dropoff) break;
    //     score_array[i].best = score;
    //     score_array[i].best_gap = score - gap_open_extend;
    //     score -= gap_extend;
    //     edit_script_row[i] = SCRIPT_GAP_IN_A;
    // }
    let mut b_size = 1usize;
    for i in 1..=n.min(*dp_mem_alloc - 1) {
        if score < -x_dropoff {
            break;
        }
        score_array[i].best = score;
        score_array[i].best_gap = score - (gap_open_extend as i32);
        score -= gap_extend as i32;
        row0.push(SCRIPT_GAP_IN_A);
        b_size = i + 1;
    }

    // NCBI reference: blast_gapalign.c:492-494
    // b_size = i;
    // best_score = 0;
    // first_b_index = 0;
    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:373-379
    // ```c
    // ALIGN_EX(..., Boolean * fence_hit)
    // ```
    let mut fence_hit = false;

    // NCBI reference: blast_gapalign.c:500-676
    // Main DP loop
    for a_index in 1..=m {
        // NCBI reference: blast_gapalign.c:514-529
        // Allocate new row in edit_script
        if a_index >= edit_script_num_rows {
            edit_script_num_rows *= 2;
            if trace_rows.len() < edit_script_num_rows {
                trace_rows.reserve(edit_script_num_rows - trace_rows.len());
            }
            if trace_offsets.len() < edit_script_num_rows {
                trace_offsets.reserve(edit_script_num_rows - trace_offsets.len());
            }
        }

        // Create new row for this a_index
        let orig_b_index = first_b_index;
        let row_capacity = b_size.saturating_sub(first_b_index) + num_extra_cells + 10;
        let edit_script_row = gap_alloc_trace_row(
            trace_rows,
            trace_offsets,
            trace_rows_used,
            row_capacity,
            orig_b_index,
        );

        // NCBI reference: blast_gapalign.c:541-545
        // matrix_row = matrix[ A[ a_index ] ];
        let qc = q_seq[a_index - 1];

        // NCBI reference: blast_gapalign.c:559-561
        // score = MININT;
        // score_gap_row = MININT;
        // last_b_index = first_b_index;
        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;

        // NCBI reference: blast_gapalign.c:563-636
        // Inner loop for each subject position
        for b_index in first_b_index..b_size {
            // NCBI reference: blast_gapalign.c:563-578 (b_size can reach N+1; no b_index < N guard).
            // NCBI reference: blast_util.c:826 (NULLB sentinel at sequence ends).
            let sc = if b_index < n { s_seq[b_index] } else { 0 };
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:569-575
            // ```c
            // if (matrix_index == FENCE_SENTRY) {
            //     if (fence_hit) { *fence_hit = 1; }
            //     break;
            // }
            // ```
            if sc == FENCE_SENTRY {
                fence_hit = true;
                break;
            }
            let score_gap_col = score_array[b_index].best_gap;
            let row = (qc as usize) * BLASTNA_SIZE;
            // NCBI reference: blast_gapalign.c:563-567 (matrix_row lookup)
            let match_score = score_matrix[row + (sc as usize)];
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:578-579
            // ```c
            // next_score = score_array[b_index].best + matrix_row[*b_ptr];
            // ```
            let next_score = score_array[b_index].best + match_score;

            // NCBI reference: blast_gapalign.c:588-599
            // Determine script based on which path gives best score
            // script = SCRIPT_SUB;
            // script_col = SCRIPT_EXTEND_GAP_B;
            // script_row = SCRIPT_EXTEND_GAP_A;
            let mut script = SCRIPT_SUB;
            let mut script_col = SCRIPT_EXTEND_GAP_B;
            let mut script_row = SCRIPT_EXTEND_GAP_A;

            // NCBI reference: blast_gapalign.c:592-598
            // if (score < score_gap_col) { script = SCRIPT_GAP_IN_B; score = score_gap_col; }
            // if (score < score_gap_row) { script = SCRIPT_GAP_IN_A; score = score_gap_row; }
            if score_val < score_gap_col {
                script = SCRIPT_GAP_IN_B;
                score_val = score_gap_col;
            }
            if score_val < score_gap_row {
                script = SCRIPT_GAP_IN_A;
                score_val = score_gap_row;
            }

            // NCBI reference: blast_gapalign.c:601-632
            // X-drop check
            if best_score - score_val > x_dropoff {
                // NCBI reference: blast_gapalign.c:603-606
                if b_index == first_b_index {
                    first_b_index += 1;
                } else {
                    score_array[b_index].best = GAP_MININT;
                }
            } else {
                // NCBI reference: blast_gapalign.c:608-631
                last_b_index = b_index;

                // Update best score
                if score_val > best_score {
                    best_score = score_val;
                    a_offset = a_index;
                    b_offset = b_index;
                }

                // NCBI reference: blast_gapalign.c:616-628
                // Update gap scores with extend flags
                score_gap_row -= gap_extend as i32;
                let score_gap_col_ext = score_gap_col - (gap_extend as i32);
                let open_gap_col = score_val - (gap_open_extend as i32);

                if score_gap_col_ext < open_gap_col {
                    score_array[b_index].best_gap = open_gap_col;
                } else {
                    score_array[b_index].best_gap = score_gap_col_ext;
                    script += script_col;  // Mark as extending existing gap
                }

                let open_gap_row = score_val - (gap_open_extend as i32);
                if score_gap_row < open_gap_row {
                    score_gap_row = open_gap_row;
                } else {
                    script += script_row;  // Mark as extending existing gap
                }

                score_array[b_index].best = score_val;
            }

            // NCBI reference: blast_gapalign.c:634-635
            // score = next_score;
            // edit_script_row[b_index] = script;
            score_val = next_score;
            let row_idx = b_index.saturating_sub(orig_b_index);
            if row_idx < edit_script_row.len() {
                edit_script_row[row_idx] = script;
            }
        }

        // NCBI reference: blast_gapalign.c:638-639
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:638-639
        // ```c
        // if (first_b_index == b_size || (fence_hit && *fence_hit)) break;
        // ```
        if first_b_index >= b_size || fence_hit {
            break;
        }

        // NCBI reference: blast_gapalign.c:641-648
        // Reallocate DP memory if needed
        gap_dp_reserve_band(score_array, dp_mem_alloc, last_b_index, num_extra_cells);

        // NCBI reference: blast_gapalign.c:652-664
        // Band contraction or expansion
        if last_b_index < b_size.saturating_sub(1) {
            b_size = last_b_index + 1;
        } else {
            // NCBI reference: blast_gapalign.c:656-663
            // Extend band with gaps
            while score_gap_row >= best_score - x_dropoff && b_size <= n && b_size < *dp_mem_alloc {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - (gap_open_extend as i32);
                score_gap_row -= gap_extend as i32;
                // NCBI reference: blast_gapalign.c:661
                let row_idx = b_size.saturating_sub(orig_b_index);
                if row_idx < edit_script_row.len() {
                    edit_script_row[row_idx] = SCRIPT_GAP_IN_A;
                } else {
                    edit_script_row.push(SCRIPT_GAP_IN_A);
                }
                b_size += 1;
            }
        }

        // NCBI reference: blast_gapalign.c:671-675
        // Sentinel
        if b_size <= n && b_size < *dp_mem_alloc {
            score_array[b_size].best = GAP_MININT;
            score_array[b_size].best_gap = GAP_MININT;
            b_size += 1;
        }
    }

    // NCBI reference: blast_gapalign.c:678-727
    // Traceback: walk from best position back to origin
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:686-688
    // ```c
    // if (fence_hit && *fence_hit) goto done;
    // ```
    if fence_hit {
        return (0, 0, 0, 0, 0, 0, 0, Vec::new());
    }

    if best_score <= 0 {
        return (0, 0, 0, 0, 0, 0, 0, Vec::new());
    }

    // NCBI reference: blast_gapalign.c:682-684
    // a_index = *a_offset;
    // b_index = *b_offset;
    // script = SCRIPT_SUB;
    let mut a_index = a_offset;
    let mut b_index = b_offset;
    let mut script = SCRIPT_SUB;

    // Collect operations in reverse order, then reverse at the end
    let mut ops_reversed: Vec<(u8, u32)> = Vec::new(); // (op_type, count)

    // NCBI reference: blast_gapalign.c:689-726
    // while (a_index > 0 || b_index > 0) {
    //     next_script = edit_script[a_index][b_index - edit_start_offset[a_index]];
    //     switch(script) {
    //     case SCRIPT_GAP_IN_A: ...
    //     case SCRIPT_GAP_IN_B: ...
    //     default: ...
    //     }
    //     GapPrelimEditBlockAdd(edit_block, script, 1);
    // }
    // NCBI reference: blast_gapalign.c:689-726 (traceback uses edit_script rows/offsets)
    let used_rows = *trace_rows_used;
    while a_index > 0 || b_index > 0 {
        // Get next script from edit_script[a_index]
        let row_idx = if a_index < used_rows { a_index } else { used_rows.saturating_sub(1) };
        let start_offset = trace_offsets.get(row_idx).copied().unwrap_or(0);
        let b_rel = b_index.saturating_sub(start_offset);
        let next_script = trace_rows
            .get(row_idx)
            .and_then(|row| row.get(b_rel).copied())
            .unwrap_or(SCRIPT_SUB);

        // NCBI reference: blast_gapalign.c:698-714
        // Determine actual operation based on current script and extend flags
        match script & SCRIPT_OP_MASK {
            x if x == SCRIPT_GAP_IN_A => {
                // NCBI reference: blast_gapalign.c:699-703
                script = next_script & SCRIPT_OP_MASK;
                if next_script & SCRIPT_EXTEND_GAP_A != 0 {
                    script = SCRIPT_GAP_IN_A;
                }
            }
            x if x == SCRIPT_GAP_IN_B => {
                // NCBI reference: blast_gapalign.c:705-709
                script = next_script & SCRIPT_OP_MASK;
                if next_script & SCRIPT_EXTEND_GAP_B != 0 {
                    script = SCRIPT_GAP_IN_B;
                }
            }
            _ => {
                // NCBI reference: blast_gapalign.c:711-713
                script = next_script & SCRIPT_OP_MASK;
            }
        }

        // NCBI reference: blast_gapalign.c:716-725
        // Advance indices based on operation
        let op_type = script & SCRIPT_OP_MASK;
        if op_type == SCRIPT_GAP_IN_A {
            // Gap in query (deletion) - only subject advances
            if b_index > 0 {
                b_index -= 1;
            }
        } else if op_type == SCRIPT_GAP_IN_B {
            // Gap in subject (insertion) - only query advances
            if a_index > 0 {
                a_index -= 1;
            }
        } else {
            // Substitution - both advance
            if a_index > 0 {
                a_index -= 1;
            }
            if b_index > 0 {
                b_index -= 1;
            }
        }

        // NCBI reference: blast_gapalign.c:726
        // GapPrelimEditBlockAdd(edit_block, script, 1);
        // Add to ops (run-length encoded)
        if let Some(last) = ops_reversed.last_mut() {
            if last.0 == op_type {
                last.1 += 1;
            } else {
                ops_reversed.push((op_type, 1));
            }
        } else {
            ops_reversed.push((op_type, 1));
        }
    }

    // Reverse to get forward order
    ops_reversed.reverse();

    // Convert to GapEditOp
    let mut edit_ops: Vec<GapEditOp> = Vec::with_capacity(ops_reversed.len());
    for (op_type, count) in ops_reversed {
        match op_type {
            x if x == SCRIPT_SUB => edit_ops.push(GapEditOp::Sub(count)),
            x if x == SCRIPT_GAP_IN_A => edit_ops.push(GapEditOp::Del(count)), // Gap in query = Del
            x if x == SCRIPT_GAP_IN_B => edit_ops.push(GapEditOp::Ins(count)), // Gap in subject = Ins
            _ => edit_ops.push(GapEditOp::Sub(count)), // Default to Sub
        }
    }

    // Compute statistics from edit script
    let (matches, mismatches, gap_opens, gap_letters) =
        stats_from_edit_ops(q_seq, s_seq, 0, 0, &edit_ops);

    let q_consumed = a_offset;
    let s_consumed = b_offset;

    (
        q_consumed,
        s_consumed,
        best_score,
        matches,
        mismatches,
        gap_opens,
        gap_letters,
        edit_ops,
    )
}

/// Extended version of traceback-capturing extension that supports reverse access.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:364-733 (ALIGN_EX)
///
/// When `reverse` is true, sequences are accessed from the end (for left extension).
pub fn extend_gapped_one_direction_with_traceback_ex(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
) -> (usize, usize, i32, usize, usize, usize, usize, Vec<GapEditOp>) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:313-319 (BLAST_GapAlignStructNew)
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    let mut gap_scratch = GapAlignScratch::new();
    extend_gapped_one_direction_with_traceback_ex_with_scratch(
        q_seq,
        s_seq,
        len1,
        len2,
        reward,
        penalty,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop,
        reverse,
        &mut gap_scratch,
    )
}

fn extend_gapped_one_direction_with_traceback_ex_with_scratch(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    gap_scratch: &mut GapAlignScratch,
) -> (usize, usize, i32, usize, usize, usize, usize, Vec<GapEditOp>) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:364-733 (ALIGN_EX)
    if !reverse {
        // Forward extension - use subsequences
        let q_sub = &q_seq[..len1.min(q_seq.len())];
        let s_sub = &s_seq[..len2.min(s_seq.len())];

        return extend_gapped_one_direction_with_traceback_with_scratch(
            q_sub,
            s_sub,
            reward,
            penalty,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            gap_scratch,
        );
    }

    let m = len1;
    let n = len2;
    if m == 0 || n == 0 {
        return (0, 0, 0, 0, 0, 0, 0, Vec::new());
    }

    #[inline(always)]
    fn get_q(q_seq: &[u8], i: usize, len1: usize, reverse: bool) -> u8 {
        if reverse {
            q_seq[len1 - i]
        } else {
            q_seq[i - 1]
        }
    }

    #[inline(always)]
    fn get_s(s_seq: &[u8], j: usize, len2: usize, reverse: bool) -> u8 {
        if j >= len2 {
            // NCBI reference: blast_gapalign.c:563-578 (b_size can reach N+1; b_ptr walks onto NULLB).
            // NCBI reference: blast_util.c:826 (NULLB sentinel at sequence ends).
            return 0;
        }
        if reverse {
            s_seq[len2 - 1 - j]
        } else {
            s_seq[j]
        }
    }

    // NCBI reference: blast_gapalign.c:423-426
    let gap_open_extend = gap_open + gap_extend;

    // NCBI reference: blast_gapalign.c:428-429
    let mut x_dropoff = x_drop;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    // NCBI reference: blast_gapalign.c:457-460
    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend + 3) as usize
    } else {
        n + 3
    };

    let dp_mem_alloc = &mut gap_scratch.dp_mem_alloc;
    let score_array = &mut gap_scratch.dp_mem;
    let trace_rows = &mut gap_scratch.trace_rows;
    let trace_offsets = &mut gap_scratch.trace_offsets;
    let trace_rows_used = &mut gap_scratch.trace_rows_used;

    // NCBI reference: blast_gapalign.c:462-468 (ALIGN_EX dp_mem realloc)
    gap_dp_reserve_initial(score_array, dp_mem_alloc, num_extra_cells);

    // NCBI reference: blast_gapalign.c:438-439 (ALIGN_EX calls s_GapPurgeState)
    gap_reset_traceback_state(trace_rows_used);
    let mut edit_script_num_rows = 100usize;

    // NCBI reference: blast_gapalign.c:476-489
    let mut score = -(gap_open_extend as i32);
    score_array[0].best = 0;
    score_array[0].best_gap = -(gap_open_extend as i32);

    let row0 = gap_alloc_trace_row(trace_rows, trace_offsets, trace_rows_used, 0, 0);
    row0.reserve(num_extra_cells);
    row0.push(0);

    let mut b_size = 1usize;
    for i in 1..=n.min(*dp_mem_alloc - 1) {
        if score < -x_dropoff {
            break;
        }
        score_array[i].best = score;
        score_array[i].best_gap = score - (gap_open_extend as i32);
        score -= gap_extend as i32;
        row0.push(SCRIPT_GAP_IN_A);
        b_size = i + 1;
    }

    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:373-379
    // ```c
    // ALIGN_EX(..., Boolean * fence_hit)
    // ```
    let mut fence_hit = false;

    for a_index in 1..=m {
        if a_index >= edit_script_num_rows {
            edit_script_num_rows *= 2;
            if trace_rows.len() < edit_script_num_rows {
                trace_rows.reserve(edit_script_num_rows - trace_rows.len());
            }
            if trace_offsets.len() < edit_script_num_rows {
                trace_offsets.reserve(edit_script_num_rows - trace_offsets.len());
            }
        }

        let orig_b_index = first_b_index;
        let row_capacity = b_size.saturating_sub(first_b_index) + num_extra_cells + 10;
        let edit_script_row = gap_alloc_trace_row(
            trace_rows,
            trace_offsets,
            trace_rows_used,
            row_capacity,
            orig_b_index,
        );

        let qc = get_q(q_seq, a_index, len1, reverse);

        let mut score_val = GAP_MININT;
        let mut score_gap_row = GAP_MININT;
        let mut last_b_index = first_b_index;

        for b_index in first_b_index..b_size {
            let sc = get_s(s_seq, b_index, len2, reverse);
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:569-575
            // ```c
            // if (matrix_index == FENCE_SENTRY) {
            //     if (fence_hit) { *fence_hit = 1; }
            //     break;
            // }
            // ```
            if sc == FENCE_SENTRY {
                fence_hit = true;
                break;
            }
            let score_gap_col = score_array[b_index].best_gap;
            let row = (qc as usize) * BLASTNA_SIZE;
            let match_score = score_matrix[row + (sc as usize)];
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:578-579
            // ```c
            // next_score = score_array[b_index].best + matrix_row[*b_ptr];
            // ```
            let next_score = score_array[b_index].best + match_score;

            let mut script = SCRIPT_SUB;
            let mut script_col = SCRIPT_EXTEND_GAP_B;
            let mut script_row = SCRIPT_EXTEND_GAP_A;

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
                    score_array[b_index].best = GAP_MININT;
                }
            } else {
                last_b_index = b_index;

                if score_val > best_score {
                    best_score = score_val;
                    a_offset = a_index;
                    b_offset = b_index;
                }

                score_gap_row -= gap_extend as i32;
                let score_gap_col_ext = score_gap_col - (gap_extend as i32);
                let open_gap_col = score_val - (gap_open_extend as i32);

                if score_gap_col_ext < open_gap_col {
                    score_array[b_index].best_gap = open_gap_col;
                } else {
                    score_array[b_index].best_gap = score_gap_col_ext;
                    script += script_col;
                }

                let open_gap_row = score_val - (gap_open_extend as i32);
                if score_gap_row < open_gap_row {
                    score_gap_row = open_gap_row;
                } else {
                    script += script_row;
                }

                score_array[b_index].best = score_val;
            }

            score_val = next_score;
            let row_idx = b_index.saturating_sub(orig_b_index);
            if row_idx < edit_script_row.len() {
                edit_script_row[row_idx] = script;
            }
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:638-639
        // ```c
        // if (first_b_index == b_size || (fence_hit && *fence_hit)) break;
        // ```
        if first_b_index >= b_size || fence_hit {
            break;
        }

        gap_dp_reserve_band(score_array, dp_mem_alloc, last_b_index, num_extra_cells);

        if last_b_index < b_size.saturating_sub(1) {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - x_dropoff && b_size <= n && b_size < *dp_mem_alloc {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - (gap_open_extend as i32);
                score_gap_row -= gap_extend as i32;
                let row_idx = b_size.saturating_sub(orig_b_index);
                if row_idx < edit_script_row.len() {
                    edit_script_row[row_idx] = SCRIPT_GAP_IN_A;
                } else {
                    edit_script_row.push(SCRIPT_GAP_IN_A);
                }
                b_size += 1;
            }
        }

        if b_size <= n && b_size < *dp_mem_alloc {
            score_array[b_size].best = GAP_MININT;
            score_array[b_size].best_gap = GAP_MININT;
            b_size += 1;
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:686-688
    // ```c
    // if (fence_hit && *fence_hit) goto done;
    // ```
    if fence_hit {
        return (0, 0, 0, 0, 0, 0, 0, Vec::new());
    }

    if best_score <= 0 {
        return (0, 0, 0, 0, 0, 0, 0, Vec::new());
    }

    let mut a_index = a_offset;
    let mut b_index = b_offset;
    let mut script = SCRIPT_SUB;

    let mut ops_reversed: Vec<(u8, u32)> = Vec::new();
    let used_rows = *trace_rows_used;
    while a_index > 0 || b_index > 0 {
        let row_idx = if a_index < used_rows { a_index } else { used_rows.saturating_sub(1) };
        let start_offset = trace_offsets.get(row_idx).copied().unwrap_or(0);
        let b_rel = b_index.saturating_sub(start_offset);
        let next_script = trace_rows
            .get(row_idx)
            .and_then(|row| row.get(b_rel).copied())
            .unwrap_or(SCRIPT_SUB);

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

        if let Some(last) = ops_reversed.last_mut() {
            if last.0 == op_type {
                last.1 += 1;
            } else {
                ops_reversed.push((op_type, 1));
            }
        } else {
            ops_reversed.push((op_type, 1));
        }
    }

    ops_reversed.reverse();

    let mut edit_ops: Vec<GapEditOp> = Vec::with_capacity(ops_reversed.len());
    for (op_type, count) in ops_reversed {
        match op_type {
            x if x == SCRIPT_SUB => edit_ops.push(GapEditOp::Sub(count)),
            x if x == SCRIPT_GAP_IN_A => edit_ops.push(GapEditOp::Del(count)),
            x if x == SCRIPT_GAP_IN_B => edit_ops.push(GapEditOp::Ins(count)),
            _ => edit_ops.push(GapEditOp::Sub(count)),
        }
    }

    // Reverse edit ops to match left-to-right order on the original sequences.
    edit_ops.reverse();

    let q_sub = &q_seq[..len1.min(q_seq.len())];
    let s_sub = &s_seq[..len2.min(s_seq.len())];
    let (matches, mismatches, gap_opens, gap_letters) =
        stats_from_edit_ops(q_sub, s_sub, 0, 0, &edit_ops);

    let q_consumed = a_offset;
    let s_consumed = b_offset;

    (
        q_consumed,
        s_consumed,
        best_score,
        matches,
        mismatches,
        gap_opens,
        gap_letters,
        edit_ops,
    )
}

/// Bidirectional gapped extension with traceback capture.
/// Extends from a seed position in both directions.
///
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2778-2831
/// BLAST_GappedAlignmentWithTraceback
///
/// Returns: (qs, qe, ss, se, score, matches, mismatches, gaps, gap_letters, edit_ops)
pub fn extend_gapped_heuristic_with_traceback(
    q_seq: &[u8],
    s_seq: &[u8],
    qs: usize,
    ss: usize,
    len: usize,
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> (usize, usize, usize, usize, i32, usize, usize, usize, usize, Vec<GapEditOp>) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:313-319 (BLAST_GapAlignStructNew)
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    let mut gap_scratch = GapAlignScratch::new();
    extend_gapped_heuristic_with_traceback_with_scratch(
        q_seq,
        s_seq,
        qs,
        ss,
        len,
        reward,
        penalty,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop,
        &mut gap_scratch,
    )
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct reuse)
pub fn extend_gapped_heuristic_with_traceback_with_scratch(
    q_seq: &[u8],
    s_seq: &[u8],
    qs: usize,
    ss: usize,
    len: usize,
    reward: i32,
    penalty: i32,
    score_matrix: &BlastnaMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    gap_scratch: &mut GapAlignScratch,
) -> (usize, usize, usize, usize, i32, usize, usize, usize, usize, Vec<GapEditOp>) {
    // Bounds validation
    if qs >= q_seq.len() || ss >= s_seq.len() {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0, Vec::new());
    }

    let len = len.min(q_seq.len() - qs).min(s_seq.len() - ss);
    if len == 0 {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0, Vec::new());
    }

    // NCBI blastn uses a single gapped start point; left extension includes it.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4601-4647
    let seed_len = 1usize;

    // Left extension (reverse direction, start included)
    let left_q_len = qs + seed_len;
    let left_s_len = ss + seed_len;

    let (left_q_consumed, left_s_consumed, left_score, _left_matches, _left_mismatches, _left_gaps, _left_gap_letters, left_edit_ops) =
        if left_q_len > 0 && left_s_len > 0 {
            extend_gapped_one_direction_with_traceback_ex_with_scratch(
                q_seq,
                s_seq,
                left_q_len,
                left_s_len,
                reward,
                penalty,
                score_matrix,
                gap_open,
                gap_extend,
                x_drop,
                true, // reverse for left extension
                gap_scratch,
            )
        } else {
            (0, 0, 0, 0, 0, 0, 0, Vec::new())
        };

    // Right extension (forward direction, start excluded)
    let right_q_start = qs + seed_len;
    let right_s_start = ss + seed_len;
    let right_q_len = q_seq.len().saturating_sub(right_q_start);
    let right_s_len = s_seq.len().saturating_sub(right_s_start);

    let (right_q_consumed, right_s_consumed, right_score, _right_matches, _right_mismatches, _right_gaps, _right_gap_letters, right_edit_ops) =
        if right_q_len > 0 && right_s_len > 0 {
            extend_gapped_one_direction_with_traceback_with_scratch(
                &q_seq[right_q_start..],
                &s_seq[right_s_start..],
                reward,
                penalty,
                score_matrix,
                gap_open,
                gap_extend,
                x_drop,
                gap_scratch,
            )
        } else {
            (0, 0, 0, 0, 0, 0, 0, Vec::new())
        };

    // Combine coordinates (NCBI: start included on left, excluded on right)
    let mut final_qs = qs + seed_len - left_q_consumed;
    let mut final_qe = qs + seed_len + right_q_consumed;
    let mut final_ss = ss + seed_len - left_s_consumed;
    let mut final_se = ss + seed_len + right_s_consumed;

    let mut total_score = left_score + right_score;

    // Combine edit scripts: left_ops + right_ops (merge if needed)
    let mut combined_edit_ops: Vec<GapEditOp> = Vec::with_capacity(
        left_edit_ops.len() + right_edit_ops.len()
    );
    combined_edit_ops.extend(left_edit_ops);

    if !right_edit_ops.is_empty() {
        if let Some(last) = combined_edit_ops.last_mut() {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/gapinfo.c:174-185 (GapPrelimEditBlockAdd merge)
            let same_type = matches!(
                (&*last, &right_edit_ops[0]),
                (GapEditOp::Sub(_), GapEditOp::Sub(_))
                    | (GapEditOp::Del(_), GapEditOp::Del(_))
                    | (GapEditOp::Ins(_), GapEditOp::Ins(_))
            );
            if same_type {
                let add = right_edit_ops[0].num();
                match last {
                    GapEditOp::Sub(n) | GapEditOp::Del(n) | GapEditOp::Ins(n) => {
                        *n += add;
                    }
                }
                combined_edit_ops.extend_from_slice(&right_edit_ops[1..]);
            } else {
                combined_edit_ops.extend_from_slice(&right_edit_ops);
            }
        } else {
            combined_edit_ops.extend_from_slice(&right_edit_ops);
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4683-4712
    // ```c
    // while (esp->size && esp->op_type[0] != eGapAlignSub) {
    //     score_left += score_params->gap_open +
    //                  esp->num[0] * score_params->gap_extend;
    //     if (esp->op_type[0] == eGapAlignDel)
    //         gap_align->subject_start += esp->num[0];
    //     else
    //         gap_align->query_start += esp->num[0];
    //     ...
    // }
    // while (esp->size && esp->op_type[i-1] != eGapAlignSub) {
    //     score_right += score_params->gap_open +
    //                  esp->num[i-1] * score_params->gap_extend;
    //     if (esp->op_type[i-1] == eGapAlignDel)
    //         gap_align->subject_stop -= esp->num[i-1];
    //     else
    //         gap_align->query_stop -= esp->num[i-1];
    //     ...
    // }
    // ```
    if !combined_edit_ops.is_empty() {
        let mut start_idx = 0usize;
        while start_idx < combined_edit_ops.len() {
            match combined_edit_ops[start_idx] {
                GapEditOp::Sub(_) => break,
                GapEditOp::Del(n) => {
                    let n = n as usize;
                    total_score += gap_open + gap_extend * (n as i32);
                    final_ss = final_ss.saturating_add(n);
                }
                GapEditOp::Ins(n) => {
                    let n = n as usize;
                    total_score += gap_open + gap_extend * (n as i32);
                    final_qs = final_qs.saturating_add(n);
                }
            }
            start_idx += 1;
        }

        let mut end_idx = combined_edit_ops.len();
        while end_idx > start_idx {
            match combined_edit_ops[end_idx - 1] {
                GapEditOp::Sub(_) => break,
                GapEditOp::Del(n) => {
                    let n = n as usize;
                    total_score += gap_open + gap_extend * (n as i32);
                    final_se = final_se.saturating_sub(n);
                }
                GapEditOp::Ins(n) => {
                    let n = n as usize;
                    total_score += gap_open + gap_extend * (n as i32);
                    final_qe = final_qe.saturating_sub(n);
                }
            }
            end_idx -= 1;
        }

        if start_idx > 0 || end_idx < combined_edit_ops.len() {
            combined_edit_ops = combined_edit_ops[start_idx..end_idx].to_vec();
        }
    }

    let (total_matches, total_mismatches, total_gaps, total_gap_letters) =
        stats_from_edit_ops(q_seq, s_seq, final_qs, final_ss, &combined_edit_ops);

    (
        final_qs,
        final_qe,
        final_ss,
        final_se,
        total_score,
        total_matches,
        total_mismatches,
        total_gaps,
        total_gap_letters,
        combined_edit_ops,
    )
}

