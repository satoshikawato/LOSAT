//! Gapped extension algorithms for BLASTN
//!
//! This module implements gapped extension algorithms with affine gap penalties,
//! including both heuristic (bidirectional) and one-directional extensions.

use super::greedy::{greedy_align_one_direction, greedy_align_one_direction_ex};

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
    let final_q_start = qs - left_q_consumed;
    let final_q_end = qs + len + right_q_consumed;
    let final_s_start = ss - left_s_consumed;
    let final_s_end = ss + len + right_s_consumed;

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
        if last_b_index < b_size.saturating_sub(1) {
            // Shrink window
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
    // NCBI BLAST: gap penalties are specified as positive values (cost), but used as negative in calculations
    let gap_open = -gap_open;
    let gap_extend = -gap_extend;

    const NEG_INF: i32 = i32::MIN / 2;
    // NCBI BLAST dynamically reallocates memory as needed without a fixed upper limit
    // Reference: blast_gapalign.c - NCBI expands the window as needed
    const MAX_WINDOW_SIZE: usize = 1_000_000;

    let m = len1;
    let n = len2;

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

    let gap_extend_abs = gap_extend.abs().max(1);
    let initial_window = (x_drop / gap_extend_abs + 3) as usize;
    let mut alloc_size = initial_window.max(100).min(MAX_WINDOW_SIZE);

    #[derive(Clone, Default)]
    struct DpCell {
        best: i32,
        best_gap: i32,
        stats: AlnStats,
        gap_stats: AlnStats,
    }

    let mut score_array: Vec<DpCell> = vec![
        DpCell {
            best: NEG_INF,
            best_gap: NEG_INF,
            stats: AlnStats::default(),
            gap_stats: AlnStats::default()
        };
        alloc_size
    ];

    let gap_open_extend = gap_open + gap_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = gap_open_extend;

    let mut score = gap_open_extend;
    let mut b_size = 1usize;
    for j in 1..=n.min(alloc_size - 1) {
        if score < -x_drop {
            break;
        }
        score_array[j].best = score;
        score_array[j].best_gap = score + gap_open_extend;
        score_array[j].stats = AlnStats {
            matches: 0,
            mismatches: 0,
            gap_opens: 1,
            gap_letters: j as u32,
        };
        score += gap_extend;
        b_size = j + 1;
    }

    let mut best_score = 0;
    let mut best_i = 0;
    let mut best_j = 0;
    let mut best_stats = AlnStats::default();
    let mut first_b_index = 0usize;
    let mut dp_cells = 0usize;

    for i in 1..=m {
        let qc = get_q(q_seq, i, len1, reverse);

        let mut score_val = NEG_INF;
        let mut score_gap_row = NEG_INF;
        let mut score_gap_row_stats = AlnStats::default();
        let mut score_stats = AlnStats::default();
        let mut last_b_index = first_b_index;

        for j in first_b_index..b_size {
            if j >= n {
                break;
            }
            let sc = get_s(s_seq, j, len2, reverse);
            dp_cells += 1;

            let mut score_gap_col = score_array[j].best_gap;
            let score_gap_col_stats = score_array[j].gap_stats;

            let is_match = qc == sc;
            let match_score = if is_match { reward } else { penalty };
            let next_score = if score_array[j].best > NEG_INF {
                score_array[j].best + match_score
            } else {
                NEG_INF
            };
            let mut next_stats = score_array[j].stats;
            if is_match {
                next_stats.matches += 1;
            } else {
                next_stats.mismatches += 1;
            }

            if score_val < score_gap_col {
                score_val = score_gap_col;
                score_stats = score_gap_col_stats;
            }
            if score_val < score_gap_row {
                score_val = score_gap_row;
                score_stats = score_gap_row_stats;
            }

            if best_score - score_val > x_drop {
                if j == first_b_index {
                    first_b_index += 1;
                } else {
                    score_array[j].best = NEG_INF;
                }
            } else {
                last_b_index = j;

                if score_val > best_score {
                    best_score = score_val;
                    best_i = i;
                    best_j = j + 1;
                    best_stats = score_stats;
                }

                score_gap_row += gap_extend;
                let open_gap_row = score_val + gap_open_extend;
                if open_gap_row > score_gap_row {
                    score_gap_row = open_gap_row;
                    score_gap_row_stats = score_stats;
                    score_gap_row_stats.gap_opens += 1;
                    score_gap_row_stats.gap_letters += 1;
                } else {
                    score_gap_row_stats.gap_letters += 1;
                }

                score_gap_col += gap_extend;
                let open_gap_col = score_val + gap_open_extend;
                if open_gap_col > score_gap_col {
                    score_array[j].best_gap = open_gap_col;
                    score_array[j].gap_stats = score_stats;
                    score_array[j].gap_stats.gap_opens += 1;
                    score_array[j].gap_stats.gap_letters += 1;
                } else {
                    score_array[j].best_gap = score_gap_col;
                    score_array[j].gap_stats.gap_letters += 1;
                }

                score_array[j].best = score_val;
                score_array[j].stats = score_stats;
            }

            score_val = next_score;
            score_stats = next_stats;
        }

        if first_b_index >= b_size {
            break;
        }

        if last_b_index + initial_window + 3 >= alloc_size && alloc_size < MAX_WINDOW_SIZE {
            let new_alloc = (last_b_index + initial_window + 100)
                .max(alloc_size * 2)
                .min(MAX_WINDOW_SIZE);
            score_array.resize(
                new_alloc,
                DpCell {
                    best: NEG_INF,
                    best_gap: NEG_INF,
                    stats: AlnStats::default(),
                    gap_stats: AlnStats::default(),
                },
            );
            alloc_size = new_alloc;
        }

        if last_b_index < b_size.saturating_sub(1) {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - x_drop
                && b_size <= n
                && b_size < alloc_size
                && b_size < MAX_WINDOW_SIZE
            {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].stats = score_gap_row_stats;
                score_array[b_size].best_gap = score_gap_row + gap_open_extend;
                score_array[b_size].gap_stats = score_gap_row_stats;
                score_array[b_size].gap_stats.gap_opens += 1;
                score_array[b_size].gap_stats.gap_letters += 1;
                score_gap_row += gap_extend;
                score_gap_row_stats.gap_letters += 1;
                b_size += 1;
            }
        }

        if b_size <= n && b_size < alloc_size {
            score_array[b_size].best = NEG_INF;
            score_array[b_size].best_gap = NEG_INF;
            b_size += 1;
        }
    }

    if best_score <= 0 {
        return (0, 0, 0, 0, 0, 0, 0, dp_cells);
    }

    (
        best_i,
        best_j,
        best_score,
        best_stats.matches as usize,
        best_stats.mismatches as usize,
        best_stats.gap_opens as usize,
        best_stats.gap_letters as usize,
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
    )
}

