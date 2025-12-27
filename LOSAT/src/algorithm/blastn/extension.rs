use super::constants::X_DROP_UNGAPPED;

/// Extend an ungapped hit in both directions using X-drop termination.
/// 
/// # Index Semantics
/// - Left extension starts at i=0 (the seed position itself), extending leftward
/// - Right extension starts at j=1 (position after seed), avoiding double-counting
/// - The seed position (q_pos, s_pos) is scored once in the left extension loop
///
/// This follows NCBI BLAST's na_ungapped.c pattern where the extension includes
/// the starting position and extends outward in both directions.
///
/// # Arguments
/// * `q_seq` - Query sequence
/// * `s_seq` - Subject sequence
/// * `q_pos` - Seed position in query
/// * `s_pos` - Seed position in subject
/// * `reward` - Score for match (typically positive, e.g., 2)
/// * `penalty` - Score for mismatch (typically negative, e.g., -3)
/// * `x_drop` - X-drop threshold for termination (optional, uses default if None)
///
/// # Returns
/// (q_start, q_end, s_start, s_end, score) - Half-open interval [start, end)
pub fn extend_hit_ungapped(
    q_seq: &[u8],
    s_seq: &[u8],
    q_pos: usize,
    s_pos: usize,
    reward: i32,
    penalty: i32,
    x_drop: Option<i32>,
) -> (usize, usize, usize, usize, i32) {
    let x_drop_val = x_drop.unwrap_or(X_DROP_UNGAPPED);
    let mut current_score = 0;
    let mut max_score = 0;

    // Left extension: starts at i=0 (seed position) and extends leftward
    // When i=0: accesses q_seq[q_pos] and s_seq[s_pos] (the seed position)
    // When i=1: accesses q_seq[q_pos-1] and s_seq[s_pos-1] (one left of seed)
    let mut best_i = 0;
    let mut i = 0;
    let max_left = q_pos.min(s_pos);

    while i <= max_left {
        let q = unsafe { *q_seq.get_unchecked(q_pos - i) };
        let s = unsafe { *s_seq.get_unchecked(s_pos - i) };
        let sc = if q == s { reward } else { penalty };
        current_score += sc;

        if current_score > max_score {
            max_score = current_score;
            best_i = i;
        } else if (max_score - current_score) > x_drop_val {
            break;
        }
        i += 1;
    }

    // Right extension: starts at j=1 to avoid double-counting seed position
    // When j=1: accesses q_seq[q_pos+1] and s_seq[s_pos+1]
    let mut current_score_r = max_score;
    let mut max_score_total = max_score;
    let mut best_j = 0;
    let mut j = 1;

    while (q_pos + j) < q_seq.len() && (s_pos + j) < s_seq.len() {
        let q = unsafe { *q_seq.get_unchecked(q_pos + j) };
        let s = unsafe { *s_seq.get_unchecked(s_pos + j) };
        let sc = if q == s { reward } else { penalty };
        current_score_r += sc;

        if current_score_r > max_score_total {
            max_score_total = current_score_r;
            best_j = j;
        } else if (max_score_total - current_score_r) > x_drop_val {
            break;
        }
        j += 1;
    }

    (
        q_pos - best_i,
        q_pos + best_j + 1,
        s_pos - best_i,
        s_pos + best_j + 1,
        max_score_total,
    )
}

