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

pub fn extend_gapped_heuristic(
    q_seq: &[u8],
    s_seq: &[u8],
    qs: usize,
    ss: usize,
    len: usize,
    reward: i32,
    penalty: i32,
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
    // Bounds validation: ensure seed coordinates are valid
    if qs >= q_seq.len() || ss >= s_seq.len() {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0, 0);
    }

    // Clamp len to available sequence length
    let len = len.min(q_seq.len() - qs).min(s_seq.len() - ss);
    if len == 0 {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0, 0);
    }

    // NCBI BLAST uses different extension algorithms based on task:
    // - megablast: greedy alignment (fast, good for high-identity sequences)
    // - blastn: DP-based alignment (slower, but handles divergent sequences better)
    // This follows NCBI BLAST's approach where blastn uses eDynProgScoreOnly
    // and megablast uses eGreedyScoreOnly.

    // First, extend to the right from the seed
    let right_suffix_q = q_seq.get(qs + len..).unwrap_or(&[]);
    let right_suffix_s = s_seq.get(ss + len..).unwrap_or(&[]);
    
    let (
        right_q_consumed,
        right_s_consumed,
        right_score,
        right_matches,
        right_mismatches,
        right_gaps,
        right_gap_letters,
        right_dp_cells,
    ) = if use_dp {
        // DP-based extension for blastn task (handles divergent sequences better)
        extend_gapped_one_direction(
            right_suffix_q,
            right_suffix_s,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop,
        )
    } else {
        // Greedy extension for megablast task (faster for high-identity sequences)
        let (q_cons, s_cons, score, matches, mismatches, gaps, gap_letters) =
            greedy_align_one_direction(
                right_suffix_q,
                right_suffix_s,
                reward,
                penalty,
                gap_open,
                gap_extend,
                x_drop,
            );
        (q_cons, s_cons, score, matches, mismatches, gaps, gap_letters, 0)
    };

    // Then extend to the left from the seed
    // For greedy alignment, use reverse flag to avoid O(N) prefix copying
    // For DP-based alignment, still need to copy (DP doesn't support reverse flag yet)
    let (
        left_q_consumed,
        left_s_consumed,
        left_score,
        left_matches,
        left_mismatches,
        left_gaps,
        left_gap_letters,
        left_dp_cells,
    ) = if use_dp {
        // DP-based extension for blastn task - use reverse flag to avoid O(n) prefix copying
        extend_gapped_one_direction_ex(
            q_seq,
            s_seq,
            qs,  // len1 = prefix length (characters before seed)
            ss,  // len2 = prefix length (characters before seed)
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop,
            true,  // reverse = true for left extension
        )
    } else {
        // Greedy extension for megablast task - use reverse flag to avoid copying
        // Pass the full sequences with the prefix length and reverse=true
        let (q_cons, s_cons, score, matches, mismatches, gaps, gap_letters) =
            greedy_align_one_direction_ex(
                q_seq,
                s_seq,
                qs,  // len1 = prefix length (characters before seed)
                ss,  // len2 = prefix length (characters before seed)
                reward,
                penalty,
                gap_open,
                gap_extend,
                x_drop,
                true,  // reverse = true for left extension
            );
        (q_cons, s_cons, score, matches, mismatches, gaps, gap_letters, 0)
    };

    // Calculate seed score
    let mut seed_score = 0;
    let mut seed_matches = 0;
    let mut seed_mismatches = 0;
    for k in 0..len {
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
    let total_dp_cells = left_dp_cells + right_dp_cells;

    // Calculate final positions
    // NCBI reference: blast_gapalign.c:4614-4615, 4646-4647
    // NCBI uses the formula:
    //   Left:  gap_align->query_start = q_start - private_q_length + 1;
    //   Right: gap_align->query_stop = q_start + private_q_length + 1;
    //
    // The +1 compensates for the DP recording pattern where the best score
    // is recorded at index k but corresponds to position k-1 (for subject).
    // NCBI's a_offset is 1-indexed, b_offset is 0-indexed.
    //
    // In LOSAT's structure, seed is handled separately (not included in extensions),
    // so we need to adjust the formula:
    //   - left_q_consumed and left_s_consumed are raw offsets from extend_gapped_one_direction_ex
    //   - right_q_consumed and right_s_consumed are raw offsets from extend_gapped_one_direction
    //
    // For subject (b_offset is 0-indexed), we need +1 to match query (a_offset is 1-indexed).
    // The overall +1 in NCBI's formula for 1-indexed output is handled by the caller.
    let final_q_start = qs - left_q_consumed;
    let final_q_end = qs + len + right_q_consumed;
    let final_s_start = ss - left_s_consumed;
    let final_s_end = ss + len + right_s_consumed;

    // Debug output for coordinate tracking
    if std::env::var("LOSAT_DEBUG_COORDS").is_ok() {
        eprintln!("[COORDS] seed: qs={}, ss={}, len={}", qs, ss, len);
        eprintln!("[COORDS] left: q_consumed={}, s_consumed={}", left_q_consumed, left_s_consumed);
        eprintln!("[COORDS] right: q_consumed={}, s_consumed={}", right_q_consumed, right_s_consumed);
        eprintln!("[COORDS] final: q={}-{}, s={}-{}", final_q_start, final_q_end, final_s_start, final_s_end);
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
        total_dp_cells,
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
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
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

    // NCBI reference: blast_gapalign.c:746-749
    // BlastGapDP* score_array;
    // struct BlastGapDP { Int4 best; Int4 best_gap; };
    #[derive(Clone, Copy)]
    struct BlastGapDP {
        best: i32,
        best_gap: i32,
    }

    const MININT: i32 = i32::MIN / 2;

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

    // NCBI reference: blast_gapalign.c:802-808
    // Dynamic memory allocation - start with initial allocation, grow as needed
    let mut dp_mem_alloc = num_extra_cells.max(1000);
    let mut score_array: Vec<BlastGapDP> = vec![
        BlastGapDP { best: MININT, best_gap: MININT };
        dp_mem_alloc
    ];

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
    for i in 1..=n.min(dp_mem_alloc - 1) {
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
        let mut score_val = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;

        // NCBI reference: blast_gapalign.c:862-912
        // Inner loop - for each subject position in the band
        for b_index in first_b_index..b_size {
            if b_index >= n {
                break;
            }
            dp_cells += 1;

            // NCBI reference: blast_gapalign.c:864-866
            // b_ptr += b_increment;
            // score_gap_col = score_array[b_index].best_gap;
            // next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
            let sc = s_seq[b_index];
            let score_gap_col = score_array[b_index].best_gap;
            let match_score = if qc == sc { reward } else { penalty };
            let next_score = if score_array[b_index].best > MININT {
                score_array[b_index].best + match_score
            } else {
                MININT
            };

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
                    score_array[b_index].best = MININT;
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
        if last_b_index + num_extra_cells + 3 >= dp_mem_alloc {
            dp_mem_alloc = (last_b_index + num_extra_cells + 100).max(dp_mem_alloc * 2);
            score_array.resize(dp_mem_alloc, BlastGapDP { best: MININT, best_gap: MININT });
        }

        // NCBI reference: blast_gapalign.c:933-952
        // Band contraction: shorten loop bounds if X-dropoff failed earlier than last row
        if last_b_index < b_size.saturating_sub(1) {
            // NCBI: b_size = last_b_index + 1
            b_size = last_b_index + 1;
        } else {
            // NCBI reference: blast_gapalign.c:946-951
            // Extend window with gaps
            while score_gap_row >= best_score - x_dropoff && b_size <= n && b_size < dp_mem_alloc {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - (gap_open_extend as i32);
                score_gap_row -= gap_extend as i32;
                b_size += 1;
            }
        }

        // NCBI reference: blast_gapalign.c:954-958
        // Sentinel
        if b_size <= n && b_size < dp_mem_alloc {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
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
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
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

    // NCBI reference: blast_gapalign.c:746-749
    // struct BlastGapDP { Int4 best; Int4 best_gap; };
    #[derive(Clone, Copy)]
    struct BlastGapDP {
        best: i32,
        best_gap: i32,
    }

    const MININT: i32 = i32::MIN / 2;

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

    // NCBI reference: blast_gapalign.c:802-808
    let mut dp_mem_alloc = num_extra_cells.max(1000);
    let mut score_array: Vec<BlastGapDP> = vec![
        BlastGapDP { best: MININT, best_gap: MININT };
        dp_mem_alloc
    ];

    // NCBI reference: blast_gapalign.c:811-822
    let mut score = -(gap_open_extend as i32);
    score_array[0].best = 0;
    score_array[0].best_gap = -(gap_open_extend as i32);

    // NCBI reference: blast_gapalign.c:815-822
    let mut b_size = 1usize;
    for i in 1..=n.min(dp_mem_alloc - 1) {
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
        let mut score_val = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;

        // NCBI reference: blast_gapalign.c:862-912
        for b_index in first_b_index..b_size {
            if b_index >= n {
                break;
            }
            dp_cells += 1;

            // NCBI reference: blast_gapalign.c:864-866
            let sc = get_s(s_seq, b_index, len2, reverse);
            let score_gap_col = score_array[b_index].best_gap;
            let match_score = if qc == sc { reward } else { penalty };
            let next_score = if score_array[b_index].best > MININT {
                score_array[b_index].best + match_score
            } else {
                MININT
            };

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
                    score_array[b_index].best = MININT;
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
        if last_b_index + num_extra_cells + 3 >= dp_mem_alloc {
            dp_mem_alloc = (last_b_index + num_extra_cells + 100).max(dp_mem_alloc * 2);
            score_array.resize(dp_mem_alloc, BlastGapDP { best: MININT, best_gap: MININT });
        }

        // NCBI reference: blast_gapalign.c:933-952
        // Band contraction: shorten loop bounds if X-dropoff failed earlier than last row
        if last_b_index < b_size.saturating_sub(1) {
            // NCBI: b_size = last_b_index + 1
            b_size = last_b_index + 1;
        } else {
            // NCBI reference: blast_gapalign.c:946-951
            while score_gap_row >= best_score - x_dropoff && b_size <= n && b_size < dp_mem_alloc {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - (gap_open_extend as i32);
                score_gap_row -= gap_extend as i32;
                b_size += 1;
            }
        }

        // NCBI reference: blast_gapalign.c:954-958
        if b_size <= n && b_size < dp_mem_alloc {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
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

/// Re-extend an existing gapped alignment with a larger X-drop for final traceback.
///
/// NCBI reference: blast_traceback.c - BLAST_GappedAlignmentWithTraceback
/// Uses gap_x_dropoff_final (100) instead of preliminary gap_x_dropoff (25-30)
///
/// This allows the alignment to extend further through low-scoring regions,
/// producing longer alignments that match NCBI BLAST output.
///
/// The function picks a seed point within the preliminary alignment and
/// re-extends in both directions using the larger X-drop value.
///
/// Returns: (final_qs, final_qe, final_ss, final_se, score, matches, mismatches, gaps, gap_letters, dp_cells)
pub fn extend_final_traceback(
    q_seq: &[u8],
    s_seq: &[u8],
    prelim_qs: usize,  // Preliminary alignment query start
    prelim_qe: usize,  // Preliminary alignment query end
    prelim_ss: usize,  // Preliminary alignment subject start
    prelim_se: usize,  // Preliminary alignment subject end
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,  // Final X-drop (100 for nucleotide)
    use_dp: bool,
) -> (usize, usize, usize, usize, i32, usize, usize, usize, usize, usize) {
    // NCBI BLAST picks the "gapped_start" point (usually the highest-scoring position within the HSP)
    // For simplicity, we use the center of the preliminary alignment as the seed point
    // This ensures we have sequence on both sides to extend

    let prelim_q_len = prelim_qe.saturating_sub(prelim_qs);
    let prelim_s_len = prelim_se.saturating_sub(prelim_ss);

    if prelim_q_len == 0 || prelim_s_len == 0 {
        return (prelim_qs, prelim_qe, prelim_ss, prelim_se, 0, 0, 0, 0, 0, 0);
    }

    // Use the center of the preliminary alignment as the seed point
    // NCBI uses gapped_start which is typically near the best-scoring region
    let seed_q_pos = prelim_qs + prelim_q_len / 2;
    let seed_s_pos = prelim_ss + prelim_s_len / 2;

    // Ensure seed positions are within bounds
    if seed_q_pos >= q_seq.len() || seed_s_pos >= s_seq.len() {
        return (prelim_qs, prelim_qe, prelim_ss, prelim_se, 0, 0, 0, 0, 0, 0);
    }

    // Re-extend from the seed point with the larger X-drop
    // This is equivalent to calling extend_gapped_heuristic with x_drop_final
    // but starting from a known good position within the alignment

    let (final_qs, final_qe, final_ss, final_se, score, matches, mismatches, gaps, gap_letters, dp_cells) =
        extend_gapped_heuristic(
            q_seq,
            s_seq,
            seed_q_pos,
            seed_s_pos,
            1, // Use a minimal seed length - the extension will find the true boundaries
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop_final,
            use_dp,
        );

    // NCBI reference: blast_gapalign.c:2969-3055 BLAST_GappedAlignmentWithTraceback
    // NCBI simply returns the final traceback coordinates without comparing to preliminary
    // The final traceback with larger x_drop should always produce better or equal results
    (final_qs, final_qe, final_ss, final_se, score, matches, mismatches, gaps, gap_letters, dp_cells)
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
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
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

    // NCBI reference: blast_gapalign.c:389-394
    // struct BlastGapDP { Int4 best; Int4 best_gap; };
    #[derive(Clone, Copy)]
    struct BlastGapDP {
        best: i32,
        best_gap: i32,
    }

    const MININT: i32 = i32::MIN / 2;

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

    // NCBI reference: blast_gapalign.c:462-468
    // Allocate DP memory
    let mut dp_mem_alloc = num_extra_cells.max(1000);
    let mut score_array: Vec<BlastGapDP> = vec![
        BlastGapDP { best: MININT, best_gap: MININT };
        dp_mem_alloc
    ];

    // NCBI reference: blast_gapalign.c:448-450, 472-474
    // Allocate traceback storage
    // edit_script[a_index] points to the script row for query position a_index
    // edit_start_offset[a_index] gives the starting b_index for that row
    let mut edit_script_num_rows = 100usize;
    let mut edit_script: Vec<Vec<u8>> = Vec::with_capacity(edit_script_num_rows);
    let mut edit_start_offset: Vec<usize> = Vec::with_capacity(edit_script_num_rows);

    // NCBI reference: blast_gapalign.c:476-489
    // Initialize row 0 (gaps in subject at start)
    let mut score = -(gap_open_extend as i32);
    score_array[0].best = 0;
    score_array[0].best_gap = -(gap_open_extend as i32);

    // First edit script row for initial gap extension
    let mut row0: Vec<u8> = Vec::with_capacity(num_extra_cells);
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
    for i in 1..=n.min(dp_mem_alloc - 1) {
        if score < -x_dropoff {
            break;
        }
        score_array[i].best = score;
        score_array[i].best_gap = score - (gap_open_extend as i32);
        score -= gap_extend as i32;
        row0.push(SCRIPT_GAP_IN_A);
        b_size = i + 1;
    }

    edit_script.push(row0);
    edit_start_offset.push(0);

    // NCBI reference: blast_gapalign.c:492-494
    // b_size = i;
    // best_score = 0;
    // first_b_index = 0;
    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;

    // NCBI reference: blast_gapalign.c:500-676
    // Main DP loop
    for a_index in 1..=m {
        // NCBI reference: blast_gapalign.c:514-529
        // Allocate new row in edit_script
        if a_index >= edit_script_num_rows {
            edit_script_num_rows *= 2;
            edit_script.reserve(edit_script_num_rows - edit_script.len());
            edit_start_offset.reserve(edit_script_num_rows - edit_start_offset.len());
        }

        // Create new row for this a_index
        let row_capacity = b_size.saturating_sub(first_b_index) + num_extra_cells + 10;
        let mut edit_script_row: Vec<u8> = vec![0; row_capacity];
        let orig_b_index = first_b_index;

        // NCBI reference: blast_gapalign.c:541-545
        // matrix_row = matrix[ A[ a_index ] ];
        let qc = q_seq[a_index - 1];

        // NCBI reference: blast_gapalign.c:559-561
        // score = MININT;
        // score_gap_row = MININT;
        // last_b_index = first_b_index;
        let mut score_val = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;

        // NCBI reference: blast_gapalign.c:563-636
        // Inner loop for each subject position
        for b_index in first_b_index..b_size {
            if b_index >= n {
                break;
            }

            // NCBI reference: blast_gapalign.c:566-578
            let sc = s_seq[b_index];
            let score_gap_col = score_array[b_index].best_gap;
            let match_score = if qc == sc { reward } else { penalty };
            let next_score = if score_array[b_index].best > MININT {
                score_array[b_index].best + match_score
            } else {
                MININT
            };

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
                    score_array[b_index].best = MININT;
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
        if first_b_index >= b_size {
            break;
        }

        // NCBI reference: blast_gapalign.c:641-648
        // Reallocate DP memory if needed
        if last_b_index + num_extra_cells + 3 >= dp_mem_alloc {
            dp_mem_alloc = (last_b_index + num_extra_cells + 100).max(dp_mem_alloc * 2);
            score_array.resize(dp_mem_alloc, BlastGapDP { best: MININT, best_gap: MININT });
        }

        // NCBI reference: blast_gapalign.c:652-664
        // Band contraction or expansion
        if last_b_index < b_size.saturating_sub(1) {
            b_size = last_b_index + 1;
        } else {
            // NCBI reference: blast_gapalign.c:656-663
            // Extend band with gaps
            while score_gap_row >= best_score - x_dropoff && b_size <= n && b_size < dp_mem_alloc {
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
        if b_size <= n && b_size < dp_mem_alloc {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
            b_size += 1;
        }

        // Store this row's edit script
        edit_script.push(edit_script_row);
        edit_start_offset.push(orig_b_index);
    }

    // NCBI reference: blast_gapalign.c:678-727
    // Traceback: walk from best position back to origin
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
    while a_index > 0 || b_index > 0 {
        // Get next script from edit_script[a_index]
        let row_idx = if a_index < edit_script.len() { a_index } else { edit_script.len() - 1 };
        let start_offset = edit_start_offset.get(row_idx).copied().unwrap_or(0);
        let b_rel = b_index.saturating_sub(start_offset);
        let next_script = edit_script.get(row_idx)
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
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut gap_opens = 0usize;
    let mut gap_letters = 0usize;

    let mut qi = 0usize;
    let mut si = 0usize;
    for op in &edit_ops {
        match op {
            GapEditOp::Sub(n) => {
                // Count matches/mismatches
                for _ in 0..*n {
                    if qi < q_seq.len() && si < s_seq.len() {
                        if q_seq[qi] == s_seq[si] {
                            matches += 1;
                        } else {
                            mismatches += 1;
                        }
                        qi += 1;
                        si += 1;
                    }
                }
            }
            GapEditOp::Del(n) => {
                gap_opens += 1;
                gap_letters += *n as usize;
                si += *n as usize;
            }
            GapEditOp::Ins(n) => {
                gap_opens += 1;
                gap_letters += *n as usize;
                qi += *n as usize;
            }
        }
    }

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
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
) -> (usize, usize, i32, usize, usize, usize, usize, Vec<GapEditOp>) {
    // For reverse extension, we need to create reversed views of the sequences
    // and then reverse the edit script at the end

    if reverse {
        // Create reversed slices
        let q_rev: Vec<u8> = q_seq[..len1].iter().rev().copied().collect();
        let s_rev: Vec<u8> = s_seq[..len2].iter().rev().copied().collect();

        let (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters, mut edit_ops) =
            extend_gapped_one_direction_with_traceback(
                &q_rev,
                &s_rev,
                reward,
                penalty,
                gap_open,
                gap_extend,
                x_drop,
            );

        // Reverse the edit script for left extension
        edit_ops.reverse();

        (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters, edit_ops)
    } else {
        // Forward extension - use subsequences
        let q_sub = &q_seq[..len1.min(q_seq.len())];
        let s_sub = &s_seq[..len2.min(s_seq.len())];

        extend_gapped_one_direction_with_traceback(
            q_sub,
            s_sub,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop,
        )
    }
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
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> (usize, usize, usize, usize, i32, usize, usize, usize, usize, Vec<GapEditOp>) {
    // Bounds validation
    if qs >= q_seq.len() || ss >= s_seq.len() {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0, Vec::new());
    }

    let len = len.min(q_seq.len() - qs).min(s_seq.len() - ss);
    if len == 0 {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0, Vec::new());
    }

    // Left extension (reverse direction)
    let left_q_len = qs;
    let left_s_len = ss;

    let (left_q_consumed, left_s_consumed, left_score, left_matches, left_mismatches, left_gaps, left_gap_letters, left_edit_ops) =
        if left_q_len > 0 && left_s_len > 0 {
            extend_gapped_one_direction_with_traceback_ex(
                q_seq,
                s_seq,
                left_q_len,
                left_s_len,
                reward,
                penalty,
                gap_open,
                gap_extend,
                x_drop,
                true, // reverse for left extension
            )
        } else {
            (0, 0, 0, 0, 0, 0, 0, Vec::new())
        };

    // Right extension (forward direction)
    let right_q_start = qs + len;
    let right_s_start = ss + len;
    let right_q_len = q_seq.len().saturating_sub(right_q_start);
    let right_s_len = s_seq.len().saturating_sub(right_s_start);

    let (right_q_consumed, right_s_consumed, right_score, right_matches, right_mismatches, right_gaps, right_gap_letters, right_edit_ops) =
        if right_q_len > 0 && right_s_len > 0 {
            extend_gapped_one_direction_with_traceback(
                &q_seq[right_q_start..],
                &s_seq[right_s_start..],
                reward,
                penalty,
                gap_open,
                gap_extend,
                x_drop,
            )
        } else {
            (0, 0, 0, 0, 0, 0, 0, Vec::new())
        };

    // Seed score
    let mut seed_matches = 0usize;
    let mut seed_mismatches = 0usize;
    for i in 0..len {
        if q_seq[qs + i] == s_seq[ss + i] {
            seed_matches += 1;
        } else {
            seed_mismatches += 1;
        }
    }
    let seed_score = (seed_matches as i32) * reward + (seed_mismatches as i32) * penalty;

    // Combine coordinates
    let final_qs = qs.saturating_sub(left_q_consumed);
    let final_qe = qs + len + right_q_consumed;
    let final_ss = ss.saturating_sub(left_s_consumed);
    let final_se = ss + len + right_s_consumed;

    let total_score = left_score + seed_score + right_score;
    let total_matches = left_matches + seed_matches + right_matches;
    let total_mismatches = left_mismatches + seed_mismatches + right_mismatches;
    let total_gaps = left_gaps + right_gaps;
    let total_gap_letters = left_gap_letters + right_gap_letters;

    // Combine edit scripts: left_ops + seed_ops + right_ops
    let mut combined_edit_ops: Vec<GapEditOp> = Vec::with_capacity(
        left_edit_ops.len() + 1 + right_edit_ops.len()
    );

    // Add left extension ops
    combined_edit_ops.extend(left_edit_ops);

    // Add seed as substitutions
    if len > 0 {
        // Try to merge with last op if it's also a Sub
        if let Some(GapEditOp::Sub(n)) = combined_edit_ops.last_mut() {
            *n += len as u32;
        } else {
            combined_edit_ops.push(GapEditOp::Sub(len as u32));
        }
    }

    // Add right extension ops (try to merge first Sub with seed)
    for op in right_edit_ops {
        if let Some(GapEditOp::Sub(n)) = combined_edit_ops.last_mut() {
            if let GapEditOp::Sub(m) = op {
                *n += m;
                continue;
            }
        }
        combined_edit_ops.push(op);
    }

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

