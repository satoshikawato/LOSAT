use super::constants::X_DROP_UNGAPPED;

/// Extend an ungapped hit in both directions using X-drop termination.
/// Returns: (q_start, q_end, s_start, s_end, score)
pub fn extend_hit_ungapped(
    q_seq: &[u8],
    s_seq: &[u8],
    q_pos: usize,
    s_pos: usize,
    reward: i32,
    penalty: i32,
) -> (usize, usize, usize, usize, i32) {
    let mut current_score = 0;
    let mut max_score = 0;

    // Left
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
        } else if (max_score - current_score) > X_DROP_UNGAPPED {
            break;
        }
        i += 1;
    }

    // Right
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
        } else if (max_score_total - current_score_r) > X_DROP_UNGAPPED {
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

