//! Extension algorithms for TBLASTX
//!
//! This module implements ungapped and gapped extension algorithms for protein sequences,
//! following NCBI BLAST's approach with X-drop termination and affine gap penalties.
//!
//! Sequences are encoded in NCBISTDAA (0-27), and scoring uses BLOSUM62 via
//! the ncbistdaa_to_blosum62 conversion table.

use crate::stats::KarlinParams;
use crate::algorithm::common::evalue::calculate_evalue_alignment_length;
use crate::utils::matrix::{blosum62_score, ncbistdaa};
use super::constants::{
    GAP_EXTEND, GAP_OPEN, SENTINEL_BYTE, SENTINEL_PENALTY,
};

/// Get the substitution matrix score for two amino acids in NCBISTDAA encoding.
/// 
/// If either character is a sentinel byte (SENTINEL_BYTE = 0, NCBI NULLB), returns
/// SENTINEL_PENALTY (-4) to trigger X-drop termination at sequence boundaries.
/// 
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:120
///   const Uint1 kProtSentinel = NULLB;
/// Reference: ncbi-blast/c++/include/algo/blast/core/ncbi_std.h:181
///   #define NULLB '\0'
#[inline(always)]
pub fn get_score(a: u8, b: u8) -> i32 {
    // Check for sentinel bytes (NCBI BLAST style sequence boundary markers)
    if a == SENTINEL_BYTE || b == SENTINEL_BYTE {
        return SENTINEL_PENALTY;
    }
    blosum62_score(a, b)
}

/// Alignment statistics propagated alongside DP scores for traceback-based calculation
#[derive(Clone, Copy, Default)]
struct ProteinAlnStats {
    matches: u32,
    mismatches: u32,
    gap_opens: u32,
    gap_letters: u32,
}

/// One-hit ungapped extension: extends from a single seed position in both directions.
/// Used when a single high-scoring seed triggers extension without two-hit requirement.
/// Returns (q_start, q_end, s_start, s_end, score, s_last_off)
/// where s_last_off is the rightmost subject position scanned (for diagonal suppression)
///
/// This implementation follows NCBI BLAST's s_BlastAaExtendOneHit:
/// 1. First, find the best scoring position within the word
/// 2. Then extend left from that position
/// 3. Then extend right from that position
///
/// # Arguments
/// * `x_drop` - X-drop threshold (LOSAT default: 11, NCBI standard: 7)
pub fn extend_hit_ungapped(
    q_seq: &[u8],
    s_seq: &[u8],
    q_pos: usize,
    s_pos: usize,
    _seed_score: i32,
    x_drop: i32,
) -> (usize, usize, usize, usize, i32, usize) {
    let k_size = 3;

    // Step 1: Find the best scoring position within the word (NCBI BLAST style)
    // This handles cases where the score drops to 0 or below within the word
    let mut sum = 0i32;
    let mut score = 0i32;
    let mut q_left_off = q_pos;
    let mut q_best_left_off = q_pos;
    let mut q_right_off = q_pos;

    for i in 0..k_size {
        if q_pos + i >= q_seq.len() || s_pos + i >= s_seq.len() {
            break;
        }
        let q_char = unsafe { *q_seq.get_unchecked(q_pos + i) };
        let s_char = unsafe { *s_seq.get_unchecked(s_pos + i) };
        // NCBI BLAST does NOT break on stop codons in word scanning
        // Stop codons have very negative scores in BLOSUM62
        sum += get_score(q_char, s_char);

        if sum > score {
            score = sum;
            q_best_left_off = q_left_off;
            q_right_off = q_pos + i;
        } else if sum <= 0 {
            sum = 0;
            q_left_off = q_pos + i + 1;
        }
    }

    // Calculate the initial hit width and positions
    let init_hit_width = q_right_off - q_best_left_off + 1;
    let q_left_off = q_best_left_off;
    let s_left_off = q_left_off + (s_pos - q_pos);
    let s_right_off = q_right_off + (s_pos - q_pos);

    // Step 2: Left extension from the best position
    // NCBI BLAST does NOT break on stop codons in s_BlastAaExtendLeft
    // Reference: aa_ungapped.c:886-921
    let mut current_score = score;
    let mut max_score = score;
    let mut left_disp = 0usize;

    let max_left = q_left_off.min(s_left_off);
    let mut i = 0;

    while i < max_left {
        let q_char = unsafe { *q_seq.get_unchecked(q_left_off - 1 - i) };
        let s_char = unsafe { *s_seq.get_unchecked(s_left_off - 1 - i) };

        // Stop codons are handled by the matrix (*-* = -4, AA-* = -4)
        // X-drop will naturally terminate the extension
        current_score += get_score(q_char, s_char);

        if current_score > max_score {
            max_score = current_score;
            left_disp = i + 1;
        }

        if (max_score - current_score) >= x_drop {
            break;
        }
        i += 1;
    }

    // Step 3: Right extension from the best position
    let mut right_score = max_score;
    let mut max_score_total = max_score;
    let mut right_disp = 0usize;
    let mut j = 0;

    let q_start_r = q_right_off + 1;
    let s_start_r = s_right_off + 1;
    let q_limit = q_seq.len();
    let s_limit = s_seq.len();

    while q_start_r + j < q_limit && s_start_r + j < s_limit {
        let q_char = unsafe { *q_seq.get_unchecked(q_start_r + j) };
        let s_char = unsafe { *s_seq.get_unchecked(s_start_r + j) };

        // Stop codons are handled by the matrix (*-* = -4, AA-* = -4)
        right_score += get_score(q_char, s_char);

        if right_score > max_score_total {
            max_score_total = right_score;
            right_disp = j + 1;
        }

        if right_score <= 0 || (max_score_total - right_score) >= x_drop {
            break;
        }
        j += 1;
    }

    // s_last_off is the rightmost subject position scanned (NCBI BLAST style)
    let s_last_off = s_start_r + j;

    // Calculate final positions
    let final_q_start = q_left_off - left_disp;
    let final_q_end = q_right_off + 1 + right_disp;
    let final_s_start = s_left_off - left_disp;
    let final_s_end = s_right_off + 1 + right_disp;
    let final_len = left_disp + init_hit_width + right_disp;

    // Sanity check: ensure length matches
    debug_assert_eq!(final_q_end - final_q_start, final_len);

    (
        final_q_start,
        final_q_end,
        final_s_start,
        final_s_end,
        max_score_total,
        s_last_off,
    )
}

/// Two-hit ungapped extension: NCBI BLAST style.
/// Given two hits L (left/first) and R (right/second) on the same diagonal,
/// extend from R to the left first. Only if the left extension reaches L,
/// proceed with right extension. This prevents fragmentation of longer HSPs.
///
/// Unlike the one-hit extension, this always returns a result (the left extension
/// at minimum), matching NCBI BLAST's behavior.
/// Returns (q_start, q_end, s_start, s_end, score, right_extended, s_last_off)
/// where s_last_off is the rightmost subject position scanned (for diagonal suppression)
///
/// # Arguments
/// * `x_drop` - X-drop threshold (LOSAT default: 11, NCBI standard: 7)
pub fn extend_hit_two_hit(
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,  // Position of first hit (L) in subject
    s_right_off: usize, // Position of second hit (R) in subject
    q_right_off: usize, // Position of second hit (R) in query
    x_drop: i32,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    let k_size = 3;

    // First, find the best scoring position within the word at R
    let mut score = 0i32;
    let mut left_score = 0i32;
    let mut right_d = 0usize;

    for i in 0..k_size {
        if q_right_off + i >= q_seq.len() || s_right_off + i >= s_seq.len() {
            break;
        }
        let q_char = unsafe { *q_seq.get_unchecked(q_right_off + i) };
        let s_char = unsafe { *s_seq.get_unchecked(s_right_off + i) };
        score += get_score(q_char, s_char);

        if score > left_score {
            left_score = score;
            right_d = i + 1; // Position is one beyond the end of the word
        }
    }

    let q_right_off = q_right_off + right_d;
    let s_right_off = s_right_off + right_d;

    // Extend to the left from R, trying to reach L
    // Note: NCBI BLAST does NOT explicitly check for stop codons in s_BlastAaExtendLeft.
    // Stop codons have very negative scores in BLOSUM62 (-4 to -5), so X-drop handles them.
    // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:886-921
    let mut current_score = 0i32;
    let mut max_score = 0i32;
    let mut left_disp = 0usize;

    let max_left = q_right_off.min(s_right_off);
    let mut i = 0;

    while i < max_left {
        let q_char = unsafe { *q_seq.get_unchecked(q_right_off - 1 - i) };
        let s_char = unsafe { *s_seq.get_unchecked(s_right_off - 1 - i) };

        // Stop codons are handled by the matrix (*-* = -4, AA-* = -4)
        // X-drop will naturally terminate the extension
        current_score += get_score(q_char, s_char);

        if current_score > max_score {
            max_score = current_score;
            left_disp = i + 1;
        }

        if (max_score - current_score) >= x_drop {
            break;
        }
        i += 1;
    }

    // Check if left extension reached the first hit
    // The distance from R to L is (s_right_off - s_left_off)
    // We need left_disp >= this distance for the extension to reach L
    let distance_to_first_hit = s_right_off.saturating_sub(s_left_off);
    let reached_first_hit = left_disp >= distance_to_first_hit;

    // Only extend to the right if left extension reached the first hit (NCBI BLAST behavior)
    let mut right_disp = 0usize;
    let mut max_score_total = max_score;
    let mut right_extended = false;
    let mut s_last_off = s_right_off; // Default to current position if no right extension

    if reached_first_hit {
        // NCBI BLAST: right_extend is set to TRUE when we enter this block,
        // regardless of how many characters are actually extended.
        // Reference: aa_ungapped.c:1140
        right_extended = true;
        
        let mut right_score = max_score;
        let q_limit = q_seq.len();
        let s_limit = s_seq.len();
        let mut j = 0;

        while (q_right_off + j) < q_limit && (s_right_off + j) < s_limit {
            let q_char = unsafe { *q_seq.get_unchecked(q_right_off + j) };
            let s_char = unsafe { *s_seq.get_unchecked(s_right_off + j) };

            // Stop codons are handled by the matrix (*-* = -4, AA-* = -4)
            right_score += get_score(q_char, s_char);

            if right_score > max_score_total {
                max_score_total = right_score;
                right_disp = j + 1;
            }

            if right_score <= 0 || (max_score_total - right_score) >= x_drop {
                break;
            }
            j += 1;
        }
        // s_last_off is the rightmost subject position scanned (NCBI BLAST style)
        s_last_off = s_right_off + j;
    }

    let q_start = q_right_off - left_disp;
    let q_end = q_right_off + right_disp;
    let s_start = s_right_off - left_disp;
    let s_end = s_right_off + right_disp;

    (q_start, q_end, s_start, s_end, max_score_total, right_extended, s_last_off)
}

/// Banded Smith-Waterman gapped extension for protein sequences with affine gap penalties.
///
/// This implements a proper dynamic programming algorithm that:
/// - Uses affine gap penalties (gap_open + gap_extend * gap_length)
/// - Uses X-drop based dynamic window (NCBI Blast_SemiGappedAlign style)
/// - Uses the BLOSUM62 substitution matrix for scoring
/// - Handles stop codons appropriately
/// - Tracks best score across ALL diagonals (not just main diagonal)
/// - Propagates alignment statistics for accurate traceback-based calculation
///
/// Returns (q_start, q_end, s_start, s_end, score, matches, mismatches, gap_opens, gap_letters)
pub fn extend_gapped_protein(
    q_seq: &[u8],
    s_seq: &[u8],
    qs: usize,
    ss: usize,
    len: usize,
    x_drop: i32,
) -> (usize, usize, usize, usize, i32, usize, usize, usize, usize) {
    // Bounds validation: ensure seed coordinates are valid
    if qs >= q_seq.len() || ss >= s_seq.len() {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0);
    }

    // Clamp len to available sequence length
    let len = len.min(q_seq.len() - qs).min(s_seq.len() - ss);
    if len == 0 {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0);
    }

    // First, extend to the right from the seed using NCBI-style semi-gapped DP
    let (right_q_consumed, right_s_consumed, right_score, right_stats) =
        extend_gapped_protein_one_direction(
            q_seq,
            s_seq,
            (qs + len).min(q_seq.len()),
            (ss + len).min(s_seq.len()), // Start from end of seed (clamped)
            true,                        // Forward direction
            x_drop,
        );

    // Then extend to the left from the seed
    let (left_q_consumed, left_s_consumed, left_score, left_stats) =
        extend_gapped_protein_one_direction(
            q_seq, s_seq, qs, ss,    // Start from beginning of seed
            false, // Backward direction
            x_drop,
        );

    // Calculate seed score and actual matches/mismatches
    let mut seed_score = 0;
    let mut seed_matches = 0;
    let mut seed_mismatches = 0;
    for k in 0..len {
        let q_char = q_seq[qs + k];
        let s_char = s_seq[ss + k];
        if q_char == ncbistdaa::STOP || s_char == ncbistdaa::STOP {
            seed_score -= 100;
            seed_mismatches += 1;
        } else {
            seed_score += get_score(q_char, s_char);
            if q_char == s_char {
                seed_matches += 1;
            } else {
                seed_mismatches += 1;
            }
        }
    }

    let total_score = left_score + seed_score + right_score;
    let total_matches = left_stats.matches as usize + seed_matches + right_stats.matches as usize;
    let total_mismatches =
        left_stats.mismatches as usize + seed_mismatches + right_stats.mismatches as usize;
    let total_gap_opens = left_stats.gap_opens as usize + right_stats.gap_opens as usize;
    let total_gap_letters = left_stats.gap_letters as usize + right_stats.gap_letters as usize;

    // Calculate final positions
    let final_q_start = qs - left_q_consumed;
    let final_q_end = qs + len + right_q_consumed;
    let final_s_start = ss - left_s_consumed;
    let final_s_end = ss + len + right_s_consumed;

    // Return full alignment boundaries
    (
        final_q_start,
        final_q_end,
        final_s_start,
        final_s_end,
        total_score,
        total_matches,
        total_mismatches,
        total_gap_opens,
        total_gap_letters,
    )
}

/// Extend alignment in one direction using NCBI-style semi-gapped DP with affine gap penalties.
///
/// This implements NCBI BLAST's Blast_SemiGappedAlign approach:
/// - X-drop based dynamic window that expands/contracts based on score
/// - Tracks best score across ALL diagonals (not just main diagonal)
/// - Propagates alignment statistics for accurate traceback-based calculation
/// - No hard-coded extension limits (controlled by X-drop termination)
///
/// Returns: (q_consumed, s_consumed, score, stats)
fn extend_gapped_protein_one_direction(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    forward: bool,
    x_drop: i32,
) -> (usize, usize, i32, ProteinAlnStats) {
    const NEG_INF: i32 = i32::MIN / 2;
    const INITIAL_BAND_WIDTH: usize = 32;

    // Determine sequence bounds
    let (q_len, s_len) = if forward {
        (q_seq.len() - q_start, s_seq.len() - s_start)
    } else {
        (q_start, s_start)
    };

    if q_len == 0 || s_len == 0 {
        return (0, 0, 0, ProteinAlnStats::default());
    }

    // Dynamic band that can expand based on X-drop
    let max_band = q_len.max(s_len).min(1000);
    let band_size = (2 * INITIAL_BAND_WIDTH + 1).min(2 * max_band + 1);

    // DP score arrays (2 rows)
    let mut m_prev = vec![NEG_INF; band_size];
    let mut m_curr = vec![NEG_INF; band_size];
    let mut iq_prev = vec![NEG_INF; band_size];
    let mut iq_curr = vec![NEG_INF; band_size];
    let mut is_prev = vec![NEG_INF; band_size];
    let mut is_curr = vec![NEG_INF; band_size];

    // Statistics arrays (2 rows) - propagate counts alongside scores
    let mut m_stats_prev = vec![ProteinAlnStats::default(); band_size];
    let mut m_stats_curr = vec![ProteinAlnStats::default(); band_size];
    let mut iq_stats_prev = vec![ProteinAlnStats::default(); band_size];
    let mut iq_stats_curr = vec![ProteinAlnStats::default(); band_size];
    let mut is_stats_prev = vec![ProteinAlnStats::default(); band_size];
    let mut is_stats_curr = vec![ProteinAlnStats::default(); band_size];

    // Track best score, position, and stats
    let mut best_score = 0;
    let mut best_q_consumed = 0;
    let mut best_s_consumed = 0;
    let mut best_stats = ProteinAlnStats::default();

    // Initialize first position
    let center = band_size / 2;
    m_prev[center] = 0;

    // Initialize leading gaps in subject (Iy[0,j] for j > 0)
    for k in (center + 1)..band_size {
        let j = k - center;
        if j <= s_len {
            is_prev[k] = GAP_OPEN + GAP_EXTEND * (j as i32);
            is_stats_prev[k] = ProteinAlnStats {
                matches: 0,
                mismatches: 0,
                gap_opens: 1,
                gap_letters: j as u32,
            };
        }
    }

    // Track the active range of the band (NCBI-style dynamic window)
    let mut left_bound = center;
    let mut right_bound = center;

    // Process each row (query position)
    let max_rows = q_len.min(s_len * 2 + 100); // Allow some diagonal drift
    for i in 1..=max_rows {
        // Reset current row
        for k in 0..band_size {
            m_curr[k] = NEG_INF;
            iq_curr[k] = NEG_INF;
            is_curr[k] = NEG_INF;
            m_stats_curr[k] = ProteinAlnStats::default();
            iq_stats_curr[k] = ProteinAlnStats::default();
            is_stats_curr[k] = ProteinAlnStats::default();
        }

        // Initialize leading gap in query (Ix[i,0])
        if i <= center {
            let gap_k = center - i;
            if gap_k < band_size {
                iq_curr[gap_k] = GAP_OPEN + GAP_EXTEND * (i as i32);
                iq_stats_curr[gap_k] = ProteinAlnStats {
                    matches: 0,
                    mismatches: 0,
                    gap_opens: 1,
                    gap_letters: i as u32,
                };
            }
        }

        let mut row_max = NEG_INF;
        let mut new_left = band_size;
        let mut new_right = 0;

        // Expand bounds slightly to allow diagonal drift
        let k_min = left_bound.saturating_sub(1);
        let k_max = (right_bound + 1).min(band_size - 1);

        // First pass: compute M and Ix
        for k in k_min..=k_max {
            let j_offset = k as isize - center as isize;
            let j = (i as isize + j_offset) as usize;

            if j == 0 || j > s_len {
                continue;
            }

            // Get sequence characters
            let (qc, sc) = if forward {
                if i > q_len {
                    continue;
                }
                (q_seq[q_start + i - 1], s_seq[s_start + j - 1])
            } else {
                if i > q_start || j > s_start {
                    continue;
                }
                (q_seq[q_start - i], s_seq[s_start - j])
            };

            // Match/mismatch score using BLOSUM62 matrix
            let match_score = if qc == ncbistdaa::STOP || sc == ncbistdaa::STOP {
                -100
            } else {
                get_score(qc, sc)
            };
            let is_match = qc == sc && qc != ncbistdaa::STOP;

            // Compute M[i,j]
            let from_m = if m_prev[k] > NEG_INF {
                m_prev[k] + match_score
            } else {
                NEG_INF
            };
            let from_iq = if iq_prev[k] > NEG_INF {
                iq_prev[k] + match_score
            } else {
                NEG_INF
            };
            let from_is = if is_prev[k] > NEG_INF {
                is_prev[k] + match_score
            } else {
                NEG_INF
            };

            let best = from_m.max(from_iq).max(from_is);
            m_curr[k] = best;

            // Propagate stats from best predecessor
            let mut stats = if best == from_m && from_m > NEG_INF {
                m_stats_prev[k]
            } else if best == from_iq && from_iq > NEG_INF {
                iq_stats_prev[k]
            } else if best == from_is && from_is > NEG_INF {
                is_stats_prev[k]
            } else {
                ProteinAlnStats::default()
            };

            // Add current match/mismatch
            if is_match {
                stats.matches += 1;
            } else {
                stats.mismatches += 1;
            }
            m_stats_curr[k] = stats;

            // Compute Ix[i,j] (gap in subject, consume query)
            if k + 1 < band_size {
                let open_gap = if m_prev[k + 1] > NEG_INF {
                    m_prev[k + 1] + GAP_OPEN + GAP_EXTEND
                } else {
                    NEG_INF
                };
                let extend_gap = if iq_prev[k + 1] > NEG_INF {
                    iq_prev[k + 1] + GAP_EXTEND
                } else {
                    NEG_INF
                };

                let best_gap = open_gap.max(extend_gap);
                iq_curr[k] = best_gap;

                // Propagate stats
                if best_gap == open_gap && open_gap > NEG_INF {
                    let mut gap_stats = m_stats_prev[k + 1];
                    gap_stats.gap_opens += 1;
                    gap_stats.gap_letters += 1;
                    iq_stats_curr[k] = gap_stats;
                } else if best_gap == extend_gap && extend_gap > NEG_INF {
                    let mut gap_stats = iq_stats_prev[k + 1];
                    gap_stats.gap_letters += 1;
                    iq_stats_curr[k] = gap_stats;
                }
            }
        }

        // Second pass: compute Iy
        for k in k_min..=k_max {
            let j_offset = k as isize - center as isize;
            let j = (i as isize + j_offset) as usize;

            if j == 0 || j > s_len {
                continue;
            }

            if k > 0 {
                let open_gap = if m_curr[k - 1] > NEG_INF {
                    m_curr[k - 1] + GAP_OPEN + GAP_EXTEND
                } else {
                    NEG_INF
                };
                let extend_gap = if is_curr[k - 1] > NEG_INF {
                    is_curr[k - 1] + GAP_EXTEND
                } else {
                    NEG_INF
                };

                let best_gap = open_gap.max(extend_gap);
                is_curr[k] = best_gap;

                // Propagate stats
                if best_gap == open_gap && open_gap > NEG_INF {
                    let mut gap_stats = m_stats_curr[k - 1];
                    gap_stats.gap_opens += 1;
                    gap_stats.gap_letters += 1;
                    is_stats_curr[k] = gap_stats;
                } else if best_gap == extend_gap && extend_gap > NEG_INF {
                    let mut gap_stats = is_stats_curr[k - 1];
                    gap_stats.gap_letters += 1;
                    is_stats_curr[k] = gap_stats;
                }
            }
        }

        // Track best score and stats across ALL diagonals (NCBI-style)
        for k in k_min..=k_max {
            let j_offset = k as isize - center as isize;
            let j = (i as isize + j_offset) as usize;

            if j == 0 || j > s_len {
                continue;
            }

            let cell_max = m_curr[k].max(iq_curr[k]).max(is_curr[k]);

            // Update row max and active bounds
            if cell_max > NEG_INF && best_score - cell_max <= x_drop {
                if cell_max > row_max {
                    row_max = cell_max;
                }
                if k < new_left {
                    new_left = k;
                }
                if k > new_right {
                    new_right = k;
                }
            }

            // Update best score if this cell is better (across ALL diagonals)
            if cell_max > best_score {
                best_score = cell_max;
                best_q_consumed = i;
                best_s_consumed = j;
                // Get stats from the best state
                best_stats = if cell_max == m_curr[k] {
                    m_stats_curr[k]
                } else if cell_max == iq_curr[k] {
                    iq_stats_curr[k]
                } else {
                    is_stats_curr[k]
                };
            }
        }

        // X-drop termination
        if best_score - row_max > x_drop {
            break;
        }

        // Update active bounds for next row (NCBI-style dynamic window)
        if new_left < band_size && new_right > 0 {
            left_bound = new_left;
            right_bound = new_right;
        }

        // Swap rows
        std::mem::swap(&mut m_prev, &mut m_curr);
        std::mem::swap(&mut iq_prev, &mut iq_curr);
        std::mem::swap(&mut is_prev, &mut is_curr);
        std::mem::swap(&mut m_stats_prev, &mut m_stats_curr);
        std::mem::swap(&mut iq_stats_prev, &mut iq_stats_curr);
        std::mem::swap(&mut is_stats_prev, &mut is_stats_curr);
    }

    (best_q_consumed, best_s_consumed, best_score, best_stats)
}

/// Convert amino acid coordinates to DNA coordinates for a given frame
pub fn convert_coords(aa_start: usize, aa_end: usize, frame: i8, dna_len: usize) -> (usize, usize) {
    let f_abs = frame.abs() as usize;
    let shift = f_abs - 1;

    if frame > 0 {
        let start_bp = aa_start * 3 + shift + 1;
        let end_bp = (aa_end) * 3 + shift;
        (start_bp, end_bp)
    } else {
        let start_bp = dna_len - (aa_start * 3 + shift);
        let end_bp_calc = dna_len - (aa_end * 3 + shift - 1);
        (start_bp, end_bp_calc)
    }
}

/// Calculate bit score and E-value for a protein alignment
///
/// This function uses alignment length as the effective search space,
/// which is appropriate for protein alignments.
pub fn calculate_statistics(score: i32, aln_len: usize, params: &KarlinParams) -> (f64, f64) {
    calculate_evalue_alignment_length(score, aln_len, params)
}

