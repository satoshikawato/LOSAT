//! Alignment statistics and region alignment functions
//!
//! This module provides functions for calculating alignment statistics
//! and performing full region alignments (for cluster-then-extend chaining).

use crate::stats::KarlinParams;
use crate::algorithm::common::evalue::calculate_evalue_database_search;

fn calculate_evalue(
    score: i32,
    q_len: usize,
    db_len: usize,
    db_num_seqs: usize,
    params: &KarlinParams,
) -> (f64, f64) {
    calculate_evalue_database_search(score, q_len, db_len, db_num_seqs, params)
}

/// Merge overlapping HSPs and filter redundant hits.
///
/// This function:
/// 1. Groups hits by query-subject pair
/// 2. Sorts hits by bit score (descending)
/// 3. Removes hits that overlap significantly with higher-scoring hits
/// 4. Returns the filtered list of non-redundant hits
/// Align a region using banded Smith-Waterman (for cluster-then-extend chaining).
/// Unlike extend_gapped_heuristic which extends from a seed, this aligns the entire region.
/// Returns: (q_start, q_end, s_start, s_end, score, matches, mismatches, gap_opens, gap_letters)
/// where q_start/q_end and s_start/s_end are 0-based coordinates within the input slices.
///
/// TRUE LOCAL Smith-Waterman with affine gap penalties.
///
/// Key difference from previous version: includes max(0, ...) reset to allow alignments
/// to start fresh at any position. This prevents "bridging" through low-quality regions
/// and produces alignments more similar to NCBI BLAST.
///
/// OPTIMIZED VERSION: Uses rolling score rows (O(band_size) memory for scores) + compact traceback
/// instead of storing full DpCell structs for every cell. This reduces memory from ~300MB to ~15MB
/// for a 10kb region, and eliminates massive allocation overhead.
#[allow(dead_code)]
pub fn align_region(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> (usize, usize, usize, usize, i32, usize, usize, usize, usize) {
    const BAND_WIDTH: usize = 256; // Wider band for NCBI BLAST compatibility
    const NEG_INF: i32 = i32::MIN / 2;

    let m = q_seq.len();
    let n = s_seq.len();

    if m == 0 || n == 0 {
        return (0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    let band_size = 2 * BAND_WIDTH + 1;

    // Traceback state encoding
    // 0 = STOP (alignment starts here - from max(0,...) reset)
    // 1 = from M (diagonal), 2 = from Ix (up), 3 = from Iy (left)
    const TB_STOP: u8 = 0;
    const TB_DIAG: u8 = 1;
    const TB_UP: u8 = 2;
    const TB_LEFT: u8 = 3;

    // Rolling score arrays for affine gap DP (only need prev and current row)
    // Each array has band_size elements
    let mut prev_m: Vec<i32> = vec![NEG_INF; band_size];
    let mut prev_ix: Vec<i32> = vec![NEG_INF; band_size];
    let mut prev_iy: Vec<i32> = vec![NEG_INF; band_size];
    let mut cur_m: Vec<i32> = vec![NEG_INF; band_size];
    let mut cur_ix: Vec<i32> = vec![NEG_INF; band_size];
    let mut cur_iy: Vec<i32> = vec![NEG_INF; band_size];

    // Compact traceback storage: for each cell, store which state gave the best score
    // Layout: flat array indexed by i * band_size + k
    let tb_size = (m + 1) * band_size;
    let mut tb_state: Vec<u8> = vec![TB_STOP; tb_size];

    // Track best score position
    let mut best_score = 0;
    let mut best_i = 0;
    let mut best_j = 0;

    // Initialize row 0: M[0,0] = 0 (can start alignment here)
    // For true local alignment, we don't initialize leading gaps - they would be negative
    // and we'd rather start fresh (score 0) than pay gap penalties
    prev_m[BAND_WIDTH] = 0;
    // tb_state[BAND_WIDTH] = TB_STOP; // Already initialized to STOP

    // Process each row
    for i in 1..=m {
        // Reset current row
        for k in 0..band_size {
            cur_m[k] = NEG_INF;
            cur_ix[k] = NEG_INF;
            cur_iy[k] = NEG_INF;
        }

        // Determine valid k range
        let k_min = if i > BAND_WIDTH { 0 } else { BAND_WIDTH - i };
        let k_max = if n + BAND_WIDTH >= i {
            (n + BAND_WIDTH - i).min(band_size - 1)
        } else {
            0
        };

        if k_min > k_max {
            // Swap rows
            std::mem::swap(&mut prev_m, &mut cur_m);
            std::mem::swap(&mut prev_ix, &mut cur_ix);
            std::mem::swap(&mut prev_iy, &mut cur_iy);
            continue;
        }

        for k in k_min..=k_max {
            let j_offset = k as isize - BAND_WIDTH as isize;
            let j = (i as isize + j_offset) as usize;

            if j == 0 || j > n {
                continue;
            }

            let q_base = q_seq[i - 1];
            let s_base = s_seq[j - 1];
            let match_score = if q_base == s_base { reward } else { penalty };

            // M[i,j] = max(0, M[i-1,j-1] + match, Ix[i-1,j-1] + match, Iy[i-1,j-1] + match)
            // The max(0, ...) is the key for TRUE LOCAL alignment
            let diag_score = prev_m[k].max(prev_ix[k]).max(prev_iy[k]);
            let m_score = if diag_score > NEG_INF {
                (diag_score + match_score).max(0)
            } else {
                // Can always start fresh with a match (if positive) or 0
                match_score.max(0)
            };
            cur_m[k] = m_score;

            // Determine traceback for M state
            let tb_idx = i * band_size + k;
            if m_score == 0 {
                // Started fresh here (max(0,...) chose 0)
                tb_state[tb_idx] = TB_STOP;
            } else if diag_score > NEG_INF && m_score == diag_score + match_score {
                // Came from previous cell
                if prev_m[k] >= prev_ix[k] && prev_m[k] >= prev_iy[k] {
                    tb_state[tb_idx] = TB_DIAG;
                } else if prev_ix[k] >= prev_iy[k] {
                    tb_state[tb_idx] = TB_UP;
                } else {
                    tb_state[tb_idx] = TB_LEFT;
                }
            } else {
                // Started fresh with just match_score (no valid predecessor)
                tb_state[tb_idx] = TB_STOP;
            }

            // Ix[i,j] = max(M[i-1,j] + gap_open + gap_extend, Ix[i-1,j] + gap_extend)
            // Note: For local alignment, we don't apply max(0,...) to gap states
            // because gaps can only extend existing alignments, not start new ones
            if k + 1 < band_size {
                let open_score = if prev_m[k + 1] > 0 {
                    prev_m[k + 1] + gap_open + gap_extend
                } else {
                    NEG_INF
                };
                let extend_score = if prev_ix[k + 1] > NEG_INF {
                    prev_ix[k + 1] + gap_extend
                } else {
                    NEG_INF
                };
                cur_ix[k] = open_score.max(extend_score);
            }

            // Iy[i,j] = max(M[i,j-1] + gap_open + gap_extend, Iy[i,j-1] + gap_extend)
            if k > 0 {
                let open_score = if cur_m[k - 1] > 0 {
                    cur_m[k - 1] + gap_open + gap_extend
                } else {
                    NEG_INF
                };
                let extend_score = if cur_iy[k - 1] > NEG_INF {
                    cur_iy[k - 1] + gap_extend
                } else {
                    NEG_INF
                };
                cur_iy[k] = open_score.max(extend_score);
            }

            // Track best score (Smith-Waterman: best anywhere in matrix)
            // Only consider M state for best score (gaps can't end an alignment optimally)
            if cur_m[k] > best_score {
                best_score = cur_m[k];
                best_i = i;
                best_j = j;
            }
        }

        // Swap rows
        std::mem::swap(&mut prev_m, &mut cur_m);
        std::mem::swap(&mut prev_ix, &mut cur_ix);
        std::mem::swap(&mut prev_iy, &mut cur_iy);
    }

    if best_score <= 0 {
        return (0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    // Traceback to find alignment start and compute stats
    // Start from best_i, best_j and trace back until we hit TB_STOP
    let mut cur_i = best_i;
    let mut cur_j = best_j;
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut gap_opens = 0usize;
    let mut gap_letters = 0usize;
    let mut in_gap_ix = false;
    let mut in_gap_iy = false;

    // We need to track which state we're in during traceback
    // Start in M state (best score is always from M)
    #[derive(Clone, Copy, PartialEq)]
    enum TbState {
        M,
        Ix,
        Iy,
    }
    let mut current_state = TbState::M;

    while cur_i > 0 && cur_j > 0 {
        let k = (cur_j as isize - cur_i as isize + BAND_WIDTH as isize) as usize;
        if k >= band_size {
            break;
        }

        let tb_idx = cur_i * band_size + k;
        let state = tb_state[tb_idx];

        // Stop if we've reached a STOP state (alignment boundary from max(0,...))
        if state == TB_STOP {
            // Count the current position as part of alignment
            if q_seq[cur_i - 1] == s_seq[cur_j - 1] {
                matches += 1;
            } else {
                mismatches += 1;
            }
            cur_i -= 1;
            cur_j -= 1;
            break;
        }

        match current_state {
            TbState::M => {
                // We're in M state, count match/mismatch
                if q_seq[cur_i - 1] == s_seq[cur_j - 1] {
                    matches += 1;
                } else {
                    mismatches += 1;
                }

                // Determine where we came from
                match state {
                    TB_DIAG => {
                        // Came from M[i-1,j-1]
                        current_state = TbState::M;
                        cur_i -= 1;
                        cur_j -= 1;
                    }
                    TB_UP => {
                        // Came from Ix[i-1,j-1]
                        current_state = TbState::Ix;
                        cur_i -= 1;
                        cur_j -= 1;
                    }
                    TB_LEFT => {
                        // Came from Iy[i-1,j-1]
                        current_state = TbState::Iy;
                        cur_i -= 1;
                        cur_j -= 1;
                    }
                    _ => break,
                }
            }
            TbState::Ix => {
                // We're in Ix state (gap in subject)
                gap_letters += 1;
                if !in_gap_ix {
                    gap_opens += 1;
                    in_gap_ix = true;
                }
                in_gap_iy = false;

                // Ix can come from M (gap open) or Ix (gap extend)
                // We need to check the previous row's M vs Ix
                if k + 1 < band_size && cur_i > 0 {
                    // Check if we came from M or Ix in previous row
                    // This is a simplification - we assume gap extend if Ix was valid
                    cur_i -= 1;
                    // Stay in Ix or go to M based on scores (simplified: check if M was better)
                    // For now, assume we came from M (gap open already counted)
                    current_state = TbState::M;
                } else {
                    break;
                }
            }
            TbState::Iy => {
                // We're in Iy state (gap in query)
                gap_letters += 1;
                if !in_gap_iy {
                    gap_opens += 1;
                    in_gap_iy = true;
                }
                in_gap_ix = false;

                // Iy can come from M (gap open) or Iy (gap extend)
                if k > 0 && cur_j > 0 {
                    cur_j -= 1;
                    // Simplified: assume we came from M
                    current_state = TbState::M;
                } else {
                    break;
                }
            }
        }
    }

    // cur_i and cur_j are now at the start of the alignment (0-based in DP terms)
    let q_start = cur_i;
    let q_end = best_i;
    let s_start = cur_j;
    let s_end = best_j;

    (
        q_start,
        q_end,
        s_start,
        s_end,
        best_score,
        matches,
        mismatches,
        gap_opens,
        gap_letters,
    )
}
