use super::result::{AlignmentResult, EditOp};
use super::traceback::{TracebackDir, TracebackMatrix};

/// Configuration for banded Smith-Waterman alignment
#[derive(Debug, Clone, Copy)]
pub struct BandedSwConfig {
    /// Band width (positions on each side of diagonal)
    pub band_width: usize,
    /// X-drop threshold for terminating extension
    pub x_drop: i32,
    /// Maximum extension length
    pub max_extension: usize,
    /// Whether to compute full traceback
    pub full_traceback: bool,
}

impl Default for BandedSwConfig {
    fn default() -> Self {
        Self {
            band_width: 32,
            x_drop: 30,
            max_extension: 10000,
            full_traceback: false,
        }
    }
}

/// Result of a gapped extension
#[derive(Debug, Clone)]
pub struct ExtensionResult {
    /// Final score
    pub score: i32,
    /// Query end position (0-based, exclusive)
    pub q_end: usize,
    /// Subject end position (0-based, exclusive)
    pub s_end: usize,
    /// Number of matches
    pub matches: usize,
    /// Number of mismatches
    pub mismatches: usize,
    /// Number of gap openings
    pub gap_opens: usize,
    /// Alignment length (columns)
    pub alignment_len: usize,
    /// Optional edit script
    pub edit_script: Option<Vec<EditOp>>,
}

/// Extend alignment in one direction using banded Smith-Waterman
///
/// This is a simplified version that tracks statistics without full traceback
pub fn extend_gapped_one_direction(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    config: &BandedSwConfig,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> ExtensionResult {
    let band = config.band_width;
    let x_drop = config.x_drop;
    let max_ext = config.max_extension;

    let q_len = query.len().saturating_sub(q_start).min(max_ext);
    let s_len = subject.len().saturating_sub(s_start).min(max_ext);

    if q_len == 0 || s_len == 0 {
        return ExtensionResult {
            score: 0,
            q_end: q_start,
            s_end: s_start,
            matches: 0,
            mismatches: 0,
            gap_opens: 0,
            alignment_len: 0,
            edit_script: None,
        };
    }

    // DP arrays for banded alignment
    let band_size = 2 * band + 1;
    let mut h_curr = vec![0i32; band_size];
    let mut h_prev = vec![0i32; band_size];
    let mut e_curr = vec![i32::MIN / 2; band_size];

    // Statistics tracking
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut gap_opens = 0usize;
    let mut alignment_len = 0usize;

    let mut best_score = 0i32;
    let mut best_q_end = q_start;
    let mut best_s_end = s_start;
    let mut best_matches = 0usize;
    let mut best_mismatches = 0usize;
    let mut best_gap_opens = 0usize;
    let mut best_alignment_len = 0usize;

    // Track gap state for statistics
    let mut in_gap = vec![false; band_size];

    for i in 0..q_len {
        let q_idx = q_start + i;
        if q_idx >= query.len() {
            break;
        }

        std::mem::swap(&mut h_curr, &mut h_prev);
        h_curr.fill(0);
        e_curr.fill(i32::MIN / 2);

        let mut row_max = i32::MIN;
        let mut row_matches = matches;
        let mut row_mismatches = mismatches;
        let mut row_gap_opens = gap_opens;
        let mut row_alignment_len = alignment_len;

        let diag_start = if i > band { i - band } else { 0 };
        let diag_end = (i + band + 1).min(s_len);

        for j in diag_start..diag_end {
            let s_idx = s_start + j;
            if s_idx >= subject.len() {
                break;
            }

            let band_idx = (j as isize - i as isize + band as isize) as usize;
            if band_idx >= band_size {
                continue;
            }

            // Match/mismatch score
            let match_score = if query[q_idx] == subject[s_idx] {
                reward
            } else {
                penalty
            };

            // Get previous diagonal value
            let prev_band_idx = if j > 0 {
                ((j - 1) as isize - (i as isize - 1) + band as isize) as usize
            } else {
                band_size // Invalid
            };

            let diag_score = if i > 0 && j > 0 && prev_band_idx < band_size {
                h_prev[prev_band_idx] + match_score
            } else if i == 0 && j == 0 {
                match_score
            } else {
                match_score
            };

            // Gap in subject (insertion)
            let up_band_idx = band_idx;
            let ins_open = if up_band_idx < band_size {
                h_prev[up_band_idx] + gap_open + gap_extend
            } else {
                i32::MIN / 2
            };
            let ins_extend = e_curr[band_idx] + gap_extend;
            e_curr[band_idx] = ins_open.max(ins_extend);

            // Gap in query (deletion)
            let left_band_idx = if band_idx > 0 {
                band_idx - 1
            } else {
                band_size
            };
            let del_open = if left_band_idx < band_size {
                h_curr[left_band_idx] + gap_open + gap_extend
            } else {
                i32::MIN / 2
            };

            // Best score for this cell
            let score = diag_score.max(e_curr[band_idx]).max(del_open).max(0);
            h_curr[band_idx] = score;

            // Track statistics
            if score > 0 {
                if score == diag_score {
                    if query[q_idx] == subject[s_idx] {
                        row_matches = matches + 1;
                    } else {
                        row_mismatches = mismatches + 1;
                    }
                    row_alignment_len = alignment_len + 1;
                    in_gap[band_idx] = false;
                } else if score == e_curr[band_idx] || score == del_open {
                    if !in_gap[band_idx] {
                        row_gap_opens = gap_opens + 1;
                    }
                    row_alignment_len = alignment_len + 1;
                    in_gap[band_idx] = true;
                }
            }

            if score > row_max {
                row_max = score;
            }

            // Update best if this is the highest score
            if score > best_score {
                best_score = score;
                best_q_end = q_idx + 1;
                best_s_end = s_idx + 1;
                best_matches = row_matches;
                best_mismatches = row_mismatches;
                best_gap_opens = row_gap_opens;
                best_alignment_len = row_alignment_len;
            }
        }

        // X-drop termination
        if best_score - row_max > x_drop {
            break;
        }

        matches = row_matches;
        mismatches = row_mismatches;
        gap_opens = row_gap_opens;
        alignment_len = row_alignment_len;
    }

    ExtensionResult {
        score: best_score,
        q_end: best_q_end,
        s_end: best_s_end,
        matches: best_matches,
        mismatches: best_mismatches,
        gap_opens: best_gap_opens,
        alignment_len: best_alignment_len,
        edit_script: None,
    }
}

/// Extend alignment bidirectionally from a seed position
pub fn extend_gapped_heuristic(
    query: &[u8],
    subject: &[u8],
    q_seed: usize,
    s_seed: usize,
    seed_len: usize,
    config: &BandedSwConfig,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> AlignmentResult {
    // Extend left (reverse sequences)
    let q_left: Vec<u8> = query[..q_seed].iter().rev().copied().collect();
    let s_left: Vec<u8> = subject[..s_seed].iter().rev().copied().collect();

    let left_result = extend_gapped_one_direction(
        &q_left, &s_left, 0, 0, config, reward, penalty, gap_open, gap_extend,
    );

    // Extend right
    let right_result = extend_gapped_one_direction(
        query,
        subject,
        q_seed + seed_len,
        s_seed + seed_len,
        config,
        reward,
        penalty,
        gap_open,
        gap_extend,
    );

    // Combine results
    let q_start = q_seed - left_result.q_end;
    let q_end = right_result.q_end;
    let s_start = s_seed - left_result.s_end;
    let s_end = right_result.s_end;

    // Calculate seed contribution
    let mut seed_matches = 0;
    let mut seed_mismatches = 0;
    for i in 0..seed_len {
        if q_seed + i < query.len() && s_seed + i < subject.len() {
            if query[q_seed + i] == subject[s_seed + i] {
                seed_matches += 1;
            } else {
                seed_mismatches += 1;
            }
        }
    }

    let total_score = left_result.score
        + right_result.score
        + (seed_matches as i32 * reward)
        + (seed_mismatches as i32 * penalty);
    let total_matches = left_result.matches + seed_matches + right_result.matches;
    let total_mismatches = left_result.mismatches + seed_mismatches + right_result.mismatches;
    let total_gap_opens = left_result.gap_opens + right_result.gap_opens;
    let total_alignment_len = left_result.alignment_len + seed_len + right_result.alignment_len;

    AlignmentResult::new(
        q_start + 1, // Convert to 1-based
        q_end,
        s_start + 1, // Convert to 1-based
        s_end,
        total_score,
        total_matches,
        total_mismatches,
        total_gap_opens,
        total_alignment_len,
    )
}

/// Extend alignment with full traceback for exact statistics
pub fn extend_with_traceback(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    config: &BandedSwConfig,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> ExtensionResult {
    let q_len = query
        .len()
        .saturating_sub(q_start)
        .min(config.max_extension);
    let s_len = subject
        .len()
        .saturating_sub(s_start)
        .min(config.max_extension);

    if q_len == 0 || s_len == 0 {
        return ExtensionResult {
            score: 0,
            q_end: q_start,
            s_end: s_start,
            matches: 0,
            mismatches: 0,
            gap_opens: 0,
            alignment_len: 0,
            edit_script: None,
        };
    }

    // Full DP with traceback
    let rows = q_len + 1;
    let cols = s_len + 1;

    let mut h = vec![vec![0i32; cols]; rows];
    let mut traceback = TracebackMatrix::new(rows, cols);

    let mut best_score = 0i32;
    let mut best_i = 0;
    let mut best_j = 0;

    for i in 1..rows {
        let q_idx = q_start + i - 1;
        if q_idx >= query.len() {
            break;
        }

        for j in 1..cols {
            let s_idx = s_start + j - 1;
            if s_idx >= subject.len() {
                break;
            }

            // Check band constraint
            let diag_diff = (i as isize - j as isize).unsigned_abs();
            if diag_diff > config.band_width {
                continue;
            }

            let match_score = if query[q_idx] == subject[s_idx] {
                reward
            } else {
                penalty
            };

            let diag = h[i - 1][j - 1] + match_score;
            let up = h[i - 1][j] + gap_open + gap_extend;
            let left = h[i][j - 1] + gap_open + gap_extend;

            let score = diag.max(up).max(left).max(0);
            h[i][j] = score;

            if score == 0 {
                traceback.set(i, j, TracebackDir::Stop);
            } else if score == diag {
                traceback.set(i, j, TracebackDir::Diag);
            } else if score == up {
                traceback.set(i, j, TracebackDir::Up);
            } else {
                traceback.set(i, j, TracebackDir::Left);
            }

            if score > best_score {
                best_score = score;
                best_i = i;
                best_j = j;
            }
        }

        // X-drop check
        let row_max = h[i].iter().max().copied().unwrap_or(0);
        if best_score - row_max > config.x_drop {
            break;
        }
    }

    // Traceback to get edit script
    let mut edit_script = Vec::new();
    let mut i = best_i;
    let mut j = best_j;

    while i > 0 && j > 0 {
        match traceback.get(i, j) {
            TracebackDir::Diag => {
                let q_idx = q_start + i - 1;
                let s_idx = s_start + j - 1;
                if query[q_idx] == subject[s_idx] {
                    edit_script.push(EditOp::Match);
                } else {
                    edit_script.push(EditOp::Mismatch);
                }
                i -= 1;
                j -= 1;
            }
            TracebackDir::Up => {
                edit_script.push(EditOp::Ins);
                i -= 1;
            }
            TracebackDir::Left => {
                edit_script.push(EditOp::Del);
                j -= 1;
            }
            TracebackDir::Stop => break,
        }
    }

    edit_script.reverse();

    // Compute statistics from edit script
    let mut matches = 0;
    let mut mismatches = 0;
    let mut gap_opens = 0;
    let mut prev_op: Option<EditOp> = None;

    for &op in &edit_script {
        match op {
            EditOp::Match => matches += 1,
            EditOp::Mismatch => mismatches += 1,
            EditOp::Ins => {
                if prev_op != Some(EditOp::Ins) {
                    gap_opens += 1;
                }
            }
            EditOp::Del => {
                if prev_op != Some(EditOp::Del) {
                    gap_opens += 1;
                }
            }
        }
        prev_op = Some(op);
    }

    ExtensionResult {
        score: best_score,
        q_end: q_start + best_i,
        s_end: s_start + best_j,
        matches,
        mismatches,
        gap_opens,
        alignment_len: edit_script.len(),
        edit_script: Some(edit_script),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extend_gapped_simple() {
        let query = b"ACGTACGTACGT";
        let subject = b"ACGTACGTACGT";

        let config = BandedSwConfig::default();
        let result = extend_gapped_one_direction(query, subject, 0, 0, &config, 1, -1, -2, -1);

        assert!(result.score > 0);
        assert!(result.matches > 0);
    }

    #[test]
    fn test_extend_with_traceback() {
        let query = b"ACGTACGT";
        let subject = b"ACGTACGT";

        let config = BandedSwConfig {
            full_traceback: true,
            ..Default::default()
        };

        let result = extend_with_traceback(query, subject, 0, 0, &config, 1, -1, -2, -1);

        assert!(result.edit_script.is_some());
        assert_eq!(result.matches, 8);
        assert_eq!(result.mismatches, 0);
    }

    #[test]
    fn test_extend_heuristic() {
        let query = b"NNNNACGTACGTNNNN";
        let subject = b"NNNNACGTACGTNNNN";

        let config = BandedSwConfig::default();
        let result = extend_gapped_heuristic(query, subject, 4, 4, 8, &config, 1, -1, -2, -1);

        assert!(result.score > 0);
        assert!(result.matches >= 8);
    }
}
