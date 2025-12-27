//! Greedy alignment algorithms for high-identity nucleotide sequences
//!
//! This module implements NCBI BLAST's greedy alignment algorithm (Zhang et al., 2000)
//! for fast alignment of high-identity sequences.
//!
//! ## NCBI BLAST Compatibility Notes
//!
//! This implementation follows NCBI BLAST's greedy_align.c and blast_gapalign.c:
//!
//! - **Non-affine greedy** (`greedy_align_one_direction_with_max_dist`):
//!   Implements BLAST_GreedyAlign with distance-based tracking and X-drop
//!
//! - **Affine greedy** (`affine_greedy_align_one_direction_with_max_dist`):
//!   Implements BLAST_AffineGreedyAlign with three-state tracking (insert/match/delete)
//!
//! Key NCBI BLAST behaviors preserved:
//! - Score normalization when reward is odd (doubles all scores)
//! - Dynamic max_dist doubling for long alignments
//! - X-drop scaled to match score units
//! - Convergence detection via diagonal bounds
//!
//! Statistics estimation:
//! - Since we don't perform full traceback, statistics (matches, mismatches,
//!   gap_opens, gap_letters) are estimated from the alignment endpoints
//! - NCBI BLAST's edit script traceback could be added for exact statistics if needed

use std::cell::RefCell;
use super::super::constants::{GREEDY_MAX_COST, GREEDY_MAX_COST_FRACTION, INVALID_OFFSET, INVALID_DIAG};
use super::super::sequence_compare::{find_first_mismatch_ex, find_first_mismatch};
use super::utilities::gdb3;

/// Thread-local memory pool for non-affine greedy alignment.
/// This avoids per-call allocation overhead by reusing memory across calls.
/// Similar to NCBI BLAST's SGreedyAlignMem structure.
pub struct GreedyAlignMem {
    /// Two rows for last_seq2_off (we swap between them)
    last_seq2_off_a: Vec<i32>,
    last_seq2_off_b: Vec<i32>,
    /// Array of maximum scores at each distance
    max_score: Vec<i32>,
    /// Current allocated size
    allocated_size: usize,
}

impl GreedyAlignMem {
    fn new() -> Self {
        Self {
            last_seq2_off_a: Vec::new(),
            last_seq2_off_b: Vec::new(),
            max_score: Vec::new(),
            allocated_size: 0,
        }
    }

    /// Ensure the memory pool has enough capacity for the given array_size and max_score_size
    fn ensure_capacity(&mut self, array_size: usize, max_score_size: usize) {
        if self.allocated_size < array_size {
            self.last_seq2_off_a.resize(array_size, -2);
            self.last_seq2_off_b.resize(array_size, -2);
            self.allocated_size = array_size;
        }
        if self.max_score.len() < max_score_size {
            self.max_score.resize(max_score_size, 0);
        }
    }

    /// Reset arrays to initial state (fill with sentinel values)
    /// Uses slice::fill for better performance than element-by-element loops
    fn reset(&mut self, array_size: usize, max_score_size: usize) {
        // Fill with sentinel values using slice::fill (more efficient than loops)
        let a_len = array_size.min(self.last_seq2_off_a.len());
        let b_len = array_size.min(self.last_seq2_off_b.len());
        let score_len = max_score_size.min(self.max_score.len());
        
        self.last_seq2_off_a[..a_len].fill(-2);
        self.last_seq2_off_b[..b_len].fill(-2);
        self.max_score[..score_len].fill(0);
    }
}

thread_local! {
    /// Thread-local memory pool for non-affine greedy alignment
    static GREEDY_MEM: RefCell<GreedyAlignMem> = RefCell::new(GreedyAlignMem::new());
}
/// Bookkeeping structure for affine greedy alignment (NCBI BLAST's SGreedyOffset).
/// When aligning two sequences, stores the largest offset into the second sequence
/// that leads to a high-scoring alignment for a given start point, tracking
/// different path endings separately for affine gap penalties.
#[derive(Clone, Copy, Default)]
pub struct GreedyOffset {
    insert_off: i32,  // Best offset for a path ending in an insertion (gap in seq2)
    match_off: i32,   // Best offset for a path ending in a match/mismatch
    delete_off: i32,  // Best offset for a path ending in a deletion (gap in seq1)
}

/// Signal that a diagonal/offset is invalid
// INVALID_OFFSET and INVALID_DIAG are imported from constants


/// Greedy alignment for high-identity nucleotide sequences.
///
/// This implements NCBI BLAST's greedy alignment algorithm (Zhang et al., 2000):
/// - Uses distance-based tracking instead of score-based
/// - Quickly scans for exact matches before doing any DP
/// - Only explores diagonals that can achieve a given distance
/// - Much faster than full DP for similar sequences
///
/// For sequences with >90% identity, this is typically 10-100x faster than full DP.
///
/// This is a wrapper that implements NCBI BLAST's dynamic max_dist doubling strategy:
/// - Start with initial max_dist based on sequence length
/// - If alignment doesn't converge, double max_dist and retry
/// - This allows finding arbitrarily long alignments
///
/// For affine gap penalties (gap_open != 0 || gap_extend != 0), uses the affine
/// greedy algorithm (NCBI BLAST's BLAST_AffineGreedyAlign).
///
/// Returns: (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters)
/// Greedy alignment for high-identity nucleotide sequences.
/// If `reverse` is true, sequences are accessed from the end (for left extension).
/// This avoids the need to copy and reverse sequences for left extension.
pub fn greedy_align_one_direction_ex(
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
) -> (usize, usize, i32, usize, usize, usize, usize) {
    // Calculate initial max_dist like NCBI BLAST:
    // max_dist = MIN(GREEDY_MAX_COST, max(len1, len2) / GREEDY_MAX_COST_FRACTION + 1)
    let max_len = len1.max(len2);
    let mut max_dist = GREEDY_MAX_COST.min(max_len / GREEDY_MAX_COST_FRACTION + 1);
    
    // Retry with doubled max_dist until convergence (NCBI BLAST approach)
    loop {
        // Choose between affine and non-affine greedy based on gap penalties
        // This matches NCBI BLAST's BLAST_AffineGreedyAlign which falls back to
        // BLAST_GreedyAlign when gap_open == 0 && gap_extend == 0
        let (result, converged) = if gap_open != 0 || gap_extend != 0 {
            affine_greedy_align_one_direction_with_max_dist(
                q_seq, s_seq, reward, penalty, gap_open, gap_extend, x_drop, max_dist,
            )
        } else {
            greedy_align_one_direction_with_max_dist(
                q_seq, s_seq, len1, len2, reward, penalty, gap_open, gap_extend, x_drop, max_dist, reverse,
            )
        };
        
        if converged {
            return result;
        }
        
        // Double max_dist and retry (NCBI BLAST's approach)
        // NCBI BLAST continues until convergence without an explicit upper limit
        max_dist *= 2;
        
        // Safety limit to prevent infinite loops on extremely divergent sequences
        // This is a practical limit that should never be reached for reasonable sequences
        // NCBI BLAST typically converges well before this limit
        if max_dist > 100000 {
            // Return best result found so far
            return result;
        }
    }
}

/// Wrapper for backward compatibility - uses slice lengths and forward direction
pub fn greedy_align_one_direction(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> (usize, usize, i32, usize, usize, usize, usize) {
    greedy_align_one_direction_ex(
        q_seq, s_seq, q_seq.len(), s_seq.len(),
        reward, penalty, gap_open, gap_extend, x_drop, false,
    )
}

/// Internal greedy alignment function with explicit max_dist parameter.
/// If `reverse` is true, sequences are accessed from the end (for left extension).
/// Returns: ((q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters), converged)
pub fn greedy_align_one_direction_with_max_dist(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    reward: i32,
    penalty: i32,
    _gap_open: i32,
    _gap_extend: i32,
    x_drop: i32,
    max_dist: usize,
    reverse: bool,
) -> ((usize, usize, i32, usize, usize, usize, usize), bool) {
    if len1 == 0 || len2 == 0 {
        return ((0, 0, 0, 0, 0, 0, 0), true);
    }

    // Calculate match and mismatch costs for distance-based tracking
    // NCBI BLAST doubles scores if reward is odd to avoid fractions
    // IMPORTANT: x_drop must also be doubled to maintain consistent units
    let (match_cost, mismatch_cost, scaled_xdrop) = if reward % 2 == 1 {
        (reward * 2, (-penalty) * 2, x_drop * 2)
    } else {
        (reward, -penalty, x_drop)
    };

    let op_cost = match_cost + mismatch_cost; // Cost of a mismatch in distance terms

    // X-drop offset for score comparison (using scaled_xdrop for consistent units)
    let xdrop_offset = ((scaled_xdrop + match_cost / 2) / op_cost.max(1) + 1) as usize;

    // Find initial run of matches
    let initial_matches = find_first_mismatch_ex(q_seq, s_seq, len1, len2, 0, 0, reverse);

    if initial_matches == len1 || initial_matches == len2 {
        // Perfect match - return immediately (converged)
        let score = (initial_matches as i32) * reward;
        return (
            (
                initial_matches,
                initial_matches,
                score,
                initial_matches,
                0,
                0,
                0,
            ),
            true,
        );
    }

    // Diagonal origin (center of diagonal space)
    // NCBI BLAST uses max_dist (not scaled_max_dist) for diag_origin
    let diag_origin = max_dist + 2;
    let array_size = 2 * diag_origin + 4;
    let max_score_size = max_dist + xdrop_offset + 2;

    // Use thread-local memory pool to avoid per-call allocation overhead
    // This is critical for performance - without pooling, non-affine greedy is 20x slower
    GREEDY_MEM.with(|mem_cell| {
        let mut mem = mem_cell.borrow_mut();
        mem.ensure_capacity(array_size, max_score_size);
        mem.reset(array_size, max_score_size);

        // Initialize distance 0
        mem.last_seq2_off_a[diag_origin] = initial_matches as i32;
        mem.max_score[xdrop_offset] = (initial_matches as i32) * match_cost;

        let mut best_dist = 0usize;
        let mut best_seq1_len = initial_matches;
        let mut best_seq2_len = initial_matches;

        let mut diag_lower = diag_origin as i32 - 1;
        let mut diag_upper = diag_origin as i32 + 1;
        let mut end1_reached = initial_matches == len1;
        let mut end2_reached = initial_matches == len2;

        // Track convergence
        let mut converged = false;

        // Use flag to track which array is "prev" (true = a is prev, false = b is prev)
        let mut use_a_as_prev = true;

        // For each distance (use max_dist for non-affine greedy)
        for d in 1..=max_dist {
            if diag_lower > diag_upper {
                converged = true;
                break; // Converged
            }

            // Set sentinel values on the "prev" array
            if use_a_as_prev {
                if diag_lower >= 1 {
                    mem.last_seq2_off_a[(diag_lower - 1) as usize] = -2;
                    mem.last_seq2_off_a[diag_lower as usize] = -2;
                }
                if (diag_upper as usize) < array_size - 1 {
                    mem.last_seq2_off_a[diag_upper as usize] = -2;
                    mem.last_seq2_off_a[(diag_upper + 1) as usize] = -2;
                }
            } else {
                if diag_lower >= 1 {
                    mem.last_seq2_off_b[(diag_lower - 1) as usize] = -2;
                    mem.last_seq2_off_b[diag_lower as usize] = -2;
                }
                if (diag_upper as usize) < array_size - 1 {
                    mem.last_seq2_off_b[diag_upper as usize] = -2;
                    mem.last_seq2_off_b[(diag_upper + 1) as usize] = -2;
                }
            }

            // X-drop score threshold (using scaled_xdrop for consistent units)
            let xdrop_idx = if d >= xdrop_offset {
                d - xdrop_offset
            } else {
                0
            };
            let xdrop_score = mem.max_score[xdrop_idx + xdrop_offset] + (op_cost * d as i32) - scaled_xdrop;
            let xdrop_score = (xdrop_score + match_cost / 2 - 1) / (match_cost / 2); // Ceiling division matching NCBI

            let mut curr_extent = 0i32;
            let mut curr_seq2_index = 0i32;
            let mut curr_diag = diag_origin as i32;
            let tmp_diag_lower = diag_lower;
            let tmp_diag_upper = diag_upper;

            // For each diagonal
            for k in tmp_diag_lower..=tmp_diag_upper {
                let ku = k as usize;
                if ku >= array_size || ku == 0 {
                    continue;
                }

                // Find largest seq2 offset that increases distance from d-1 to d
                // Access the "prev" array based on the flag
                let (prev_k_plus, prev_k, prev_k_minus) = if use_a_as_prev {
                    let p_plus = if ku + 1 < array_size { mem.last_seq2_off_a[ku + 1] } else { -2 };
                    let p = mem.last_seq2_off_a[ku];
                    let p_minus = if ku >= 1 { mem.last_seq2_off_a[ku - 1] } else { -2 };
                    (p_plus, p, p_minus)
                } else {
                    let p_plus = if ku + 1 < array_size { mem.last_seq2_off_b[ku + 1] } else { -2 };
                    let p = mem.last_seq2_off_b[ku];
                    let p_minus = if ku >= 1 { mem.last_seq2_off_b[ku - 1] } else { -2 };
                    (p_plus, p, p_minus)
                };

                let mut seq2_index = prev_k_plus.max(prev_k) + 1;
                seq2_index = seq2_index.max(prev_k_minus);

                let seq1_index = seq2_index + k - diag_origin as i32;

                if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                    // X-drop test failed or invalid diagonal
                    if k == diag_lower {
                        diag_lower += 1;
                    } else {
                        // Write to "curr" array
                        if use_a_as_prev {
                            mem.last_seq2_off_b[ku] = -2;
                        } else {
                            mem.last_seq2_off_a[ku] = -2;
                        }
                    }
                    continue;
                }

                diag_upper = k;

                // Slide down diagonal until mismatch
                let seq1_idx = seq1_index as usize;
                let seq2_idx = seq2_index as usize;

                if seq1_idx < len1 && seq2_idx < len2 {
                    let matches = find_first_mismatch_ex(q_seq, s_seq, len1, len2, seq1_idx, seq2_idx, reverse);
                    let new_seq1_index = seq1_index + matches as i32;
                    let new_seq2_index = seq2_index + matches as i32;

                    // Write to "curr" array
                    if use_a_as_prev {
                        mem.last_seq2_off_b[ku] = new_seq2_index;
                    } else {
                        mem.last_seq2_off_a[ku] = new_seq2_index;
                    }

                    // Track best extent
                    let extent = new_seq1_index + new_seq2_index;
                    if extent > curr_extent {
                        curr_extent = extent;
                        curr_seq2_index = new_seq2_index;
                        curr_diag = k;
                    }

                    // Clamp bounds
                    if new_seq2_index as usize >= len2 {
                        diag_lower = k + 1;
                        end2_reached = true;
                    }
                    if new_seq1_index as usize >= len1 {
                        diag_upper = k - 1;
                        end1_reached = true;
                    }
                } else {
                    // Write to "curr" array
                    if use_a_as_prev {
                        mem.last_seq2_off_b[ku] = seq2_index;
                    } else {
                        mem.last_seq2_off_a[ku] = seq2_index;
                    }
                }
            }

            // Compute max score for this distance
            let curr_score = (curr_extent * match_cost) / 2 - (d as i32) * op_cost;

            if curr_score >= mem.max_score[d - 1 + xdrop_offset] {
                mem.max_score[d + xdrop_offset] = curr_score;
                best_dist = d;
                best_seq2_len = curr_seq2_index as usize;
                best_seq1_len = (curr_seq2_index + curr_diag - diag_origin as i32) as usize;
            } else {
                mem.max_score[d + xdrop_offset] = mem.max_score[d - 1 + xdrop_offset];
            }

            // Check convergence
            if diag_lower > diag_upper {
                converged = true;
                break;
            }

            // Expand bounds for next distance
            if !end2_reached {
                diag_lower -= 1;
            }
            if !end1_reached {
                diag_upper += 1;
            }

            // Swap arrays by toggling the flag
            use_a_as_prev = !use_a_as_prev;
        }

        // Calculate final statistics
        // For greedy alignment, distance = mismatches + gaps
        // gap_letters is the difference in consumed lengths (indicates indels)
        let gap_letters = best_seq1_len.abs_diff(best_seq2_len);
        // If there are gap letters, there's at least one gap open
        // Without traceback, we estimate gap_opens = 1 if any gaps exist
        let gap_opens = if gap_letters > 0 { 1 } else { 0 };
        // Mismatches = distance - gap_letters (distance includes both mismatches and gaps)
        let mismatches = best_dist.saturating_sub(gap_letters);
        let matches = best_seq1_len.min(best_seq2_len).saturating_sub(mismatches);

        // Calculate score
        let score = (matches as i32) * reward + (mismatches as i32) * penalty;

        (
            (
                best_seq1_len,
                best_seq2_len,
                score,
                matches,
                mismatches,
                gap_opens,
                gap_letters,
            ),
            converged,
        )
    })
}

/// Affine greedy alignment function (NCBI BLAST's BLAST_AffineGreedyAlign).
/// This handles affine gap penalties properly by tracking three separate offsets
/// for each diagonal: paths ending in insertion, match, or deletion.
///
/// Returns: ((q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters), converged)
pub fn affine_greedy_align_one_direction_with_max_dist(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    in_gap_open: i32,
    in_gap_extend: i32,
    x_drop: i32,
    max_dist: usize,
) -> ((usize, usize, i32, usize, usize, usize, usize), bool) {
    let len1 = q_seq.len();
    let len2 = s_seq.len();

    if len1 == 0 || len2 == 0 {
        return ((0, 0, 0, 0, 0, 0, 0), true);
    }

    // Score normalization (NCBI BLAST approach):
    // Make sure bits of match_score don't disappear if divided by 2
    // IMPORTANT: In NCBI BLAST, mismatch_score is passed as -score_params->penalty,
    // where score_params->penalty is stored as a NEGATIVE value (e.g., -3).
    // So -(-3) = 3, meaning mismatch_score is the POSITIVE penalty magnitude.
    // In blemir, penalty is stored as a POSITIVE value (e.g., 3), so we use it directly.
    let (match_score, mismatch_penalty, xdrop_threshold, gap_open_in, gap_extend_in) =
        if reward % 2 == 1 {
            (reward * 2, penalty * 2, x_drop * 2, in_gap_open * 2, in_gap_extend * 2)
        } else {
            (reward, penalty, x_drop, in_gap_open, in_gap_extend)
        };

    // Fill in derived scores and penalties (NCBI BLAST approach)
    // op_cost = match_score + mismatch_penalty (both positive, so op_cost is positive)
    // This represents the "cost" of a mismatch in distance units
    let match_score_half = match_score / 2;
    let mut op_cost = match_score + mismatch_penalty;
    let mut gap_open = gap_open_in;
    let mut gap_extend = gap_extend_in + match_score_half;
    let score_common_factor = gdb3(&mut op_cost, &mut gap_open, &mut gap_extend);
    let gap_open_extend = gap_open + gap_extend;
    let max_penalty = op_cost.max(gap_open_extend);

    // Scaled max_dist for affine alignment
    // With the correct sign convention (op_cost positive), gap_extend should always be positive
    let scaled_max_dist = (max_dist as i32) * gap_extend;

    // Diagonal origin (center of diagonal space)
    let diag_origin = max_dist + 2;
    let array_size = 2 * diag_origin + 4;

    // For affine greedy, we need to track diag_lower and diag_upper for ALL distances
    // (not just current), because contributions can come from distances < d-1
    // With correct sign convention, max_penalty should always be positive
    
    // Debug assertions to catch sign issues
    debug_assert!(op_cost > 0, "op_cost must be positive, got {}", op_cost);
    debug_assert!(gap_extend > 0, "gap_extend must be positive, got {}", gap_extend);
    debug_assert!(max_penalty > 0, "max_penalty must be positive, got {}", max_penalty);
    debug_assert!(scaled_max_dist >= 0, "scaled_max_dist must be non-negative, got {}", scaled_max_dist);
    
    // Safe conversion from i32 to usize with bounds checking
    if max_penalty < 0 || scaled_max_dist < 0 {
        // Sign error - return early with non-convergence
        return ((0, 0, 0, 0, 0, 0, 0), false);
    }
    
    let max_penalty_usize = max_penalty as usize;
    let scaled_max_dist_usize = scaled_max_dist as usize;
    let bounds_size = scaled_max_dist_usize + 1 + max_penalty_usize + 1;
    
    let mut diag_lower_arr: Vec<i32> = vec![INVALID_DIAG; bounds_size];
    let mut diag_upper_arr: Vec<i32> = vec![-INVALID_DIAG; bounds_size];

    // Initialize negative distance bounds with empty ranges
    // (for distances < max_penalty, which map to indices 0..max_penalty_usize)
    for i in 0..max_penalty_usize {
        diag_lower_arr[i] = INVALID_DIAG;
        diag_upper_arr[i] = -INVALID_DIAG;
    }

    // last_seq2_off[d][k] stores GreedyOffset for diagonal k at distance d
    // We need to keep all rows for affine alignment (contributions from d - max_penalty)
    // Allocate enough rows: scaled_max_dist + max_penalty + 2
    let num_rows = scaled_max_dist_usize + max_penalty_usize + 2;
    
    // Check allocation sizes to prevent memory exhaustion
    let total_elements = num_rows.checked_mul(array_size);
    if total_elements.is_none() || total_elements.unwrap() > 100_000_000 {
        // Allocation would be too large, return early with non-convergence
        return ((0, 0, 0, 0, 0, 0, 0), false);
    }
    
    let mut last_seq2_off: Vec<Vec<GreedyOffset>> = vec![vec![GreedyOffset {
        insert_off: INVALID_OFFSET,
        match_off: INVALID_OFFSET,
        delete_off: INVALID_OFFSET,
    }; array_size]; num_rows];

    // Max score at each distance for X-drop
    let xdrop_offset = ((xdrop_threshold + match_score_half) / score_common_factor.max(1) + 1) as usize;
    let max_score_size = scaled_max_dist_usize + xdrop_offset + 2;
    let mut max_score: Vec<i32> = vec![0; max_score_size];

    // Find initial run of matches
    let initial_matches = find_first_mismatch(q_seq, s_seq, 0, 0);

    if initial_matches == len1 || initial_matches == len2 {
        // Perfect match - return immediately (converged)
        let score = (initial_matches as i32) * reward;
        return (
            (initial_matches, initial_matches, score, initial_matches, 0, 0, 0),
            true,
        );
    }

    // Initialize distance 0
    last_seq2_off[0][diag_origin].match_off = initial_matches as i32;
    last_seq2_off[0][diag_origin].insert_off = INVALID_OFFSET;
    last_seq2_off[0][diag_origin].delete_off = INVALID_OFFSET;
    max_score[xdrop_offset] = (initial_matches as i32) * match_score;
    diag_lower_arr[max_penalty_usize] = diag_origin as i32;
    diag_upper_arr[max_penalty_usize] = diag_origin as i32;

    let mut best_dist = 0i32;
    let mut best_diag = diag_origin as i32;
    let mut best_seq1_len = initial_matches;
    let mut best_seq2_len = initial_matches;
    let _ = best_diag; // Suppress unused warning for initial value

    // Set up for distance 1
    let mut curr_diag_lower = diag_origin as i32 - 1;
    let mut curr_diag_upper = diag_origin as i32 + 1;
    let mut end1_diag = 0i32;
    let mut end2_diag = 0i32;
    let mut num_nonempty_dist = 1i32;
    let mut d = 1i32;

    let mut converged = false;

    // Helper function to safely access diag bounds
    let get_diag_lower = |arr: &[i32], d: i32, max_pen: usize| -> i32 {
        let idx = (d + max_pen as i32) as usize;
        if idx < arr.len() { arr[idx] } else { INVALID_DIAG }
    };
    let get_diag_upper = |arr: &[i32], d: i32, max_pen: usize| -> i32 {
        let idx = (d + max_pen as i32) as usize;
        if idx < arr.len() { arr[idx] } else { -INVALID_DIAG }
    };

    // For each distance
    while d <= scaled_max_dist {
        // Compute X-dropoff score threshold
        let xdrop_idx = if d as usize >= xdrop_offset {
            (d as usize) - xdrop_offset
        } else {
            0
        };
        let xdrop_score_raw = max_score[xdrop_idx + xdrop_offset] 
            + score_common_factor * d - xdrop_threshold;
        let xdrop_score = ((xdrop_score_raw as f64) / (match_score_half as f64)).ceil() as i32;
        let xdrop_score = xdrop_score.max(0);

        let mut curr_extent = 0i32;
        let mut curr_seq2_index = 0i32;
        let mut curr_diag = 0i32;
        let tmp_diag_lower = curr_diag_lower;
        let tmp_diag_upper = curr_diag_upper;

        // For each valid diagonal
        for k in tmp_diag_lower..=tmp_diag_upper {
            let ku = k as usize;
            if ku >= array_size || ku == 0 {
                continue;
            }

            // Find best offset for DELETE (gap in seq1) - look at k+1 diagonal
            let mut seq2_index_del = INVALID_OFFSET;
            
            // From gap opening (match -> delete)
            let d_open = d - gap_open_extend;
            if d_open >= 0 {
                let dl = get_diag_lower(&diag_lower_arr, d_open, max_penalty_usize);
                let du = get_diag_upper(&diag_upper_arr, d_open, max_penalty_usize);
                if k + 1 >= dl && k + 1 <= du && (d_open as usize) < last_seq2_off.len() {
                    let ku1 = (k + 1) as usize;
                    if ku1 < array_size {
                        seq2_index_del = last_seq2_off[d_open as usize][ku1].match_off;
                    }
                }
            }
            
            // From gap extension (delete -> delete)
            let d_ext = d - gap_extend;
            if d_ext >= 0 {
                let dl = get_diag_lower(&diag_lower_arr, d_ext, max_penalty_usize);
                let du = get_diag_upper(&diag_upper_arr, d_ext, max_penalty_usize);
                if k + 1 >= dl && k + 1 <= du && (d_ext as usize) < last_seq2_off.len() {
                    let ku1 = (k + 1) as usize;
                    if ku1 < array_size {
                        let ext_off = last_seq2_off[d_ext as usize][ku1].delete_off;
                        if ext_off > seq2_index_del {
                            seq2_index_del = ext_off;
                        }
                    }
                }
            }
            
            // Save delete offset (deletion means seq2 offset slips by one)
            let du = d as usize;
            if du < last_seq2_off.len() {
                last_seq2_off[du][ku].delete_off = if seq2_index_del == INVALID_OFFSET {
                    INVALID_OFFSET
                } else {
                    seq2_index_del + 1
                };
            }

            // Find best offset for INSERT (gap in seq2) - look at k-1 diagonal
            let mut seq2_index_ins = INVALID_OFFSET;
            
            // From gap opening (match -> insert)
            if d_open >= 0 && k >= 1 {
                let dl = get_diag_lower(&diag_lower_arr, d_open, max_penalty_usize);
                let du_bound = get_diag_upper(&diag_upper_arr, d_open, max_penalty_usize);
                if k - 1 >= dl && k - 1 <= du_bound && (d_open as usize) < last_seq2_off.len() {
                    let km1 = (k - 1) as usize;
                    if km1 < array_size {
                        seq2_index_ins = last_seq2_off[d_open as usize][km1].match_off;
                    }
                }
            }
            
            // From gap extension (insert -> insert)
            if d_ext >= 0 && k >= 1 {
                let dl = get_diag_lower(&diag_lower_arr, d_ext, max_penalty_usize);
                let du_bound = get_diag_upper(&diag_upper_arr, d_ext, max_penalty_usize);
                if k - 1 >= dl && k - 1 <= du_bound && (d_ext as usize) < last_seq2_off.len() {
                    let km1 = (k - 1) as usize;
                    if km1 < array_size {
                        let ext_off = last_seq2_off[d_ext as usize][km1].insert_off;
                        if ext_off > seq2_index_ins {
                            seq2_index_ins = ext_off;
                        }
                    }
                }
            }
            
            // Save insert offset (insertion doesn't change seq2 offset)
            if du < last_seq2_off.len() {
                last_seq2_off[du][ku].insert_off = seq2_index_ins;
            }

            // Compare with mismatch path (from diagonal k at d - op_cost)
            let mut seq2_index = last_seq2_off[du][ku].insert_off.max(last_seq2_off[du][ku].delete_off);
            
            let d_mismatch = d - op_cost;
            if d_mismatch >= 0 {
                let dl = get_diag_lower(&diag_lower_arr, d_mismatch, max_penalty_usize);
                let du_bound = get_diag_upper(&diag_upper_arr, d_mismatch, max_penalty_usize);
                if k >= dl && k <= du_bound && (d_mismatch as usize) < last_seq2_off.len() {
                    let match_off = last_seq2_off[d_mismatch as usize][ku].match_off;
                    if match_off != INVALID_OFFSET {
                        seq2_index = seq2_index.max(match_off + 1);
                    }
                }
            }

            // Choose seq1 offset to remain on diagonal k
            let seq1_index = seq2_index + k - diag_origin as i32;

            // X-dropoff test
            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == curr_diag_lower {
                    curr_diag_lower += 1;
                } else if du < last_seq2_off.len() {
                    last_seq2_off[du][ku].match_off = INVALID_OFFSET;
                }
                continue;
            }
            curr_diag_upper = k;

            // Slide down diagonal until mismatch
            let seq1_idx = seq1_index as usize;
            let seq2_idx = seq2_index as usize;

            let matches = if seq1_idx < len1 && seq2_idx < len2 {
                find_first_mismatch(q_seq, s_seq, seq1_idx, seq2_idx)
            } else {
                0
            };

            let new_seq1_index = seq1_index + matches as i32;
            let new_seq2_index = seq2_index + matches as i32;

            // Save match offset
            if du < last_seq2_off.len() {
                last_seq2_off[du][ku].match_off = new_seq2_index;
            }

            // Track best extent
            let extent = new_seq1_index + new_seq2_index;
            if extent > curr_extent {
                curr_extent = extent;
                curr_seq2_index = new_seq2_index;
                curr_diag = k;
            }

            // Clamp bounds to avoid walking off sequences
            if new_seq1_index as usize >= len1 {
                curr_diag_upper = k;
                end1_diag = k - 1;
            }
            if new_seq2_index as usize >= len2 {
                curr_diag_lower = k;
                end2_diag = k + 1;
            }
        }

        // Compute maximum score for distance d
        let curr_score = curr_extent * match_score_half - d * score_common_factor;

        // Update best if this is better
        let prev_max = if d >= 1 && (d - 1) as usize + xdrop_offset < max_score.len() {
            max_score[(d - 1) as usize + xdrop_offset]
        } else {
            0
        };

        if curr_score > prev_max {
            if (d as usize) + xdrop_offset < max_score.len() {
                max_score[(d as usize) + xdrop_offset] = curr_score;
            }
            best_dist = d;
            best_diag = curr_diag;
            best_seq2_len = curr_seq2_index as usize;
            best_seq1_len = (curr_seq2_index + best_diag - diag_origin as i32) as usize;
        } else if (d as usize) + xdrop_offset < max_score.len() {
            max_score[(d as usize) + xdrop_offset] = prev_max;
        }

        // Save diagonal bounds for this distance
        let bounds_idx = (d + max_penalty) as usize;
        if bounds_idx < diag_lower_arr.len() {
            if curr_diag_lower <= curr_diag_upper {
                num_nonempty_dist += 1;
                diag_lower_arr[bounds_idx] = curr_diag_lower;
                diag_upper_arr[bounds_idx] = curr_diag_upper;
            } else {
                diag_lower_arr[bounds_idx] = INVALID_DIAG;
                diag_upper_arr[bounds_idx] = -INVALID_DIAG;
            }
        }

        // Check if we should decrement num_nonempty_dist
        let old_bounds_idx = (d - max_penalty + max_penalty) as usize;
        if old_bounds_idx < diag_lower_arr.len() {
            let old_lower = diag_lower_arr[old_bounds_idx];
            let old_upper = diag_upper_arr[old_bounds_idx];
            if old_lower <= old_upper {
                num_nonempty_dist -= 1;
            }
        }

        // Convergence check: max_penalty consecutive empty ranges
        if num_nonempty_dist == 0 {
            converged = true;
            break;
        }

        // Compute diagonal range for next distance
        d += 1;
        
        let d_goe = d - gap_open_extend;
        let d_ge = d - gap_extend;
        let d_op = d - op_cost;

        let lower_goe = get_diag_lower(&diag_lower_arr, d_goe, max_penalty_usize);
        let lower_ge = get_diag_lower(&diag_lower_arr, d_ge, max_penalty_usize);
        let lower_op = get_diag_lower(&diag_lower_arr, d_op, max_penalty_usize);
        
        curr_diag_lower = lower_goe.min(lower_ge) - 1;
        curr_diag_lower = curr_diag_lower.min(lower_op);
        
        if end2_diag > 0 {
            curr_diag_lower = curr_diag_lower.max(end2_diag);
        }

        let upper_goe = get_diag_upper(&diag_upper_arr, d_goe, max_penalty_usize);
        let upper_ge = get_diag_upper(&diag_upper_arr, d_ge, max_penalty_usize);
        let upper_op = get_diag_upper(&diag_upper_arr, d_op, max_penalty_usize);
        
        curr_diag_upper = upper_goe.max(upper_ge) + 1;
        curr_diag_upper = curr_diag_upper.max(upper_op);
        
        if end1_diag > 0 {
            curr_diag_upper = curr_diag_upper.min(end1_diag);
        }
    }

    if !converged {
        // Did not converge - return best result found so far
        // Calculate statistics from best alignment
        let alignment_len = best_seq1_len.max(best_seq2_len);
        let gap_letters = best_seq1_len.abs_diff(best_seq2_len);
        
        // Estimate statistics (without full traceback)
        // For affine greedy, we estimate based on distance and gap penalties
        let estimated_gaps = if gap_letters > 0 { 1 } else { 0 };
        let estimated_mismatches = (best_dist as usize).saturating_sub(estimated_gaps);
        let matches = alignment_len.saturating_sub(estimated_mismatches).saturating_sub(gap_letters);
        
        let score = (matches as i32) * reward 
            + (estimated_mismatches as i32) * penalty
            - (estimated_gaps as i32) * in_gap_open
            - (gap_letters as i32) * in_gap_extend;
        
        return (
            (best_seq1_len, best_seq2_len, score, matches, estimated_mismatches, estimated_gaps, gap_letters),
            false,
        );
    }

    // Calculate final statistics
    // For affine greedy, we need to estimate statistics from the alignment
    // Since we don't have full traceback, we estimate based on the distance and gap penalties
    let alignment_len = best_seq1_len.max(best_seq2_len);
    let gap_letters = best_seq1_len.abs_diff(best_seq2_len);
    
    // Estimate gap opens and mismatches from distance
    // In affine greedy: distance = sum of (gap_open_extend for each gap) + (gap_extend for each additional gap letter) + (op_cost for each mismatch)
    // We approximate: if there are gaps, assume one gap open
    let estimated_gap_opens = if gap_letters > 0 { 1 } else { 0 };
    
    // Remaining distance after accounting for gaps
    let gap_cost = if gap_letters > 0 {
        gap_open_extend + (gap_letters.saturating_sub(1) as i32) * gap_extend
    } else {
        0
    };
    let remaining_dist = (best_dist - gap_cost).max(0);
    let estimated_mismatches = if op_cost > 0 {
        (remaining_dist / op_cost) as usize
    } else {
        0
    };
    
    let matches = alignment_len.saturating_sub(estimated_mismatches).saturating_sub(gap_letters);
    
    // Calculate final score
    let score = (matches as i32) * reward 
        + (estimated_mismatches as i32) * penalty
        - (estimated_gap_opens as i32) * in_gap_open
        - (gap_letters as i32) * in_gap_extend;

    (
        (best_seq1_len, best_seq2_len, score, matches, estimated_mismatches, estimated_gap_opens, gap_letters),
        converged,
    )
}
