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
        // DP-based extension for blastn task - still needs prefix copying
        let q_prefix: Vec<u8> = q_seq[..qs].iter().rev().copied().collect();
        let s_prefix: Vec<u8> = s_seq[..ss].iter().rev().copied().collect();
        extend_gapped_one_direction(
            &q_prefix,
            &s_prefix,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop,
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
/// This implements NCBI BLAST's Blast_SemiGappedAlign approach:
/// - X-drop based dynamic window that expands/contracts based on score
/// - Tracks best score across ALL diagonals
/// - Propagates alignment statistics for accurate traceback-based calculation
/// - No hard-coded extension limits (controlled by X-drop termination)
///
/// Returns: (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters, dp_cells)
/// The dp_cells count is for diagnostic purposes.
pub fn extend_gapped_one_direction(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> (usize, usize, i32, usize, usize, usize, usize, usize) {
    // NCBI BLAST-style adaptive banding implementation
    // Key insight from Blast_SemiGappedAlign:
    // - Use dynamic window bounds (first_b_index to b_size) that expand/contract based on X-drop
    // - Window naturally expands as needed, but with a maximum limit to prevent explosion
    // - Reallocate memory when window needs to grow

    const NEG_INF: i32 = i32::MIN / 2;
    // Maximum window size to prevent computational explosion on very long high-identity alignments
    // Increased from 5000 to 50000 to allow alignments with more gaps (NCBI BLAST can produce
    // alignments with 800+ gap characters that require wider bands to traverse)
    // Note: NCBI BLAST uses additional mechanisms (greedy alignment, fences) that we don't have yet
    const MAX_WINDOW_SIZE: usize = 50000;

    let m = q_seq.len();
    let n = s_seq.len();

    if m == 0 || n == 0 {
        return (0, 0, 0, 0, 0, 0, 0, 0);
    }

    // Calculate initial window size based on X-drop (NCBI-style)
    // num_extra_cells = x_dropoff / gap_extend + 3
    let gap_extend_abs = gap_extend.abs().max(1);
    let initial_window = (x_drop / gap_extend_abs + 3) as usize;

    // Start with a reasonable initial allocation, will grow as needed (up to MAX_WINDOW_SIZE)
    let mut alloc_size = initial_window.max(100).min(MAX_WINDOW_SIZE);

    // DP score arrays - indexed by j (subject position)
    // We use a 1D array approach like NCBI BLAST
    #[derive(Clone, Default)]
    struct DpCell {
        best: i32,
        best_gap: i32, // best score ending in a gap
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

    // Initialize row 0
    let gap_open_extend = gap_open + gap_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = gap_open_extend; // Cost to open gap at position 0

    // Initialize leading gaps in subject (j > 0)
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

    // Track best score and position
    let mut best_score = 0;
    let mut best_i = 0;
    let mut best_j = 0;
    let mut best_stats = AlnStats::default();

    // Dynamic window bounds (NCBI-style)
    let mut first_b_index = 0usize;

    // DP cell counter for diagnostics
    let mut dp_cells = 0usize;

    // Process each row (query position)
    for i in 1..=m {
        let qc = q_seq[i - 1];

        // Running scores for this row
        let mut score_val = NEG_INF;
        let mut score_gap_row = NEG_INF; // Best score ending in gap in query (Ix)
        let mut score_gap_row_stats = AlnStats::default();
        let mut score_stats = AlnStats::default();
        let mut last_b_index = first_b_index;

        for j in first_b_index..b_size {
            let sc = if j < n { s_seq[j] } else { break };
            dp_cells += 1;

            // Get previous column's gap score
            let score_gap_col = score_array[j].best_gap;
            let score_gap_col_stats = score_array[j].gap_stats;

            // Compute match/mismatch score
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

            // Best of: continue from M, continue from Ix (gap in query), continue from Iy (gap in subject)
            if score_val < score_gap_col {
                score_val = score_gap_col;
                score_stats = score_gap_col_stats;
            }
            if score_val < score_gap_row {
                score_val = score_gap_row;
                score_stats = score_gap_row_stats;
            }

            // X-drop check
            if best_score - score_val > x_drop {
                // Failed X-drop - mark this cell as invalid
                if j == first_b_index {
                    first_b_index += 1;
                } else {
                    score_array[j].best = NEG_INF;
                }
            } else {
                last_b_index = j;

                // Update best score
                if score_val > best_score {
                    best_score = score_val;
                    best_i = i;
                    best_j = j + 1; // Convert to 1-based
                    best_stats = score_stats;
                }

                // Update gap scores for next iteration
                // Gap in query (Ix): extend existing gap or open new gap
                let extend_gap_row = score_gap_row + gap_extend;
                let open_gap_row = score_val + gap_open_extend;
                if open_gap_row > extend_gap_row {
                    score_gap_row = open_gap_row;
                    score_gap_row_stats = score_stats;
                    score_gap_row_stats.gap_opens += 1;
                    score_gap_row_stats.gap_letters += 1;
                } else {
                    score_gap_row = extend_gap_row;
                    score_gap_row_stats.gap_letters += 1;
                }

                // Gap in subject (Iy): extend existing gap or open new gap
                let extend_gap_col = score_gap_col + gap_extend;
                let open_gap_col = score_val + gap_open_extend;
                if open_gap_col > extend_gap_col {
                    score_array[j].best_gap = open_gap_col;
                    score_array[j].gap_stats = score_stats;
                    score_array[j].gap_stats.gap_opens += 1;
                    score_array[j].gap_stats.gap_letters += 1;
                } else {
                    score_array[j].best_gap = extend_gap_col;
                    score_array[j].gap_stats.gap_letters += 1;
                }

                // Store current score
                score_array[j].best = score_val;
                score_array[j].stats = score_stats;
            }

            // Move to next cell
            score_val = next_score;
            score_stats = next_stats;
        }

        // Check if all positions failed X-drop
        if first_b_index >= b_size {
            break;
        }

        // Expand window if needed (NCBI-style adaptive banding, with MAX_WINDOW_SIZE limit)
        if last_b_index + initial_window + 3 >= alloc_size && alloc_size < MAX_WINDOW_SIZE {
            // Need to grow the array (but respect MAX_WINDOW_SIZE)
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
            // This row ended earlier than last row - shrink window
            b_size = last_b_index + 1;
        } else {
            // Extend window if we can continue with gaps (respect MAX_WINDOW_SIZE)
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

        // Ensure we have a sentinel
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

