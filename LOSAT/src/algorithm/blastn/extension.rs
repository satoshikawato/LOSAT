use super::constants::X_DROP_UNGAPPED;

/// Extend an ungapped hit in both directions using X-drop termination.
/// 
/// This implementation follows NCBI BLAST's s_NuclUngappedExtendExact:
/// Reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:148-243
/// 
/// Key differences from protein extension:
/// 1. Uses sum-based scoring: accumulate sum, reset to 0 when sum > 0
/// 2. X-drop check: sum < X (where X is negative) instead of max_score - current_score > X
/// 3. Right extension uses dynamic X_current: X_current = (-score > X) ? -score : X
/// 
/// Returns: (q_start, q_end, s_start, s_end, score)
pub fn extend_hit_ungapped(
    q_seq: &[u8],
    s_seq: &[u8],
    q_pos: usize,
    s_pos: usize,
    reward: i32,
    penalty: i32,
) -> (usize, usize, usize, usize, i32) {
    // X is passed as negative value in NCBI BLAST (line 744: -(cutoffs->x_dropoff))
    // We use positive value internally, so negate it for comparison
    let x_drop_neg = -X_DROP_UNGAPPED;
    
    let mut score = 0i32;
    let mut sum = 0i32;
    
    // Left extension
    // Reference: na_ungapped.c:188-204
    let mut q_beg = q_pos;
    let mut i = 0;
    let max_left = q_pos.min(s_pos);
    
    while i <= max_left {
        let q = unsafe { *q_seq.get_unchecked(q_pos - i) };
        let s = unsafe { *s_seq.get_unchecked(s_pos - i) };
        let sc = if q == s { reward } else { penalty };
        
        sum += sc;
        
        if sum > 0 {
            q_beg = q_pos - i;
            score += sum;
            sum = 0;
        } else if sum < x_drop_neg {
            // X-drop termination: sum becomes more negative than X
            // Reference: na_ungapped.c:201-202
            break;
        }
        i += 1;
    }
    
    let q_start = q_beg;
    let s_start = s_pos - (q_pos - q_beg);
    
    // Right extension
    // Reference: na_ungapped.c:217-239
    let mut q_end = q_pos;
    let mut j = 0;
    let max_right = (q_seq.len() - q_pos).min(s_seq.len() - s_pos);
    
    // X_current is used to break out of loop if score goes negative
    // Reference: na_ungapped.c:223-230
    let mut x_current = x_drop_neg;
    sum = 0;
    
    while j < max_right {
        let q = unsafe { *q_seq.get_unchecked(q_pos + j) };
        let s = unsafe { *s_seq.get_unchecked(s_pos + j) };
        let sc = if q == s { reward } else { penalty };
        
        sum += sc;
        
        if sum > 0 {
            q_end = q_pos + j + 1;
            score += sum;
            // Update X_current dynamically: X_current = (-score > X) ? -score : X
            // Reference: na_ungapped.c:230
            if -score > x_drop_neg {
                x_current = -score;
            } else {
                x_current = x_drop_neg;
            }
            sum = 0;
        } else if sum < x_current {
            // X-drop termination with dynamic threshold
            // Reference: na_ungapped.c:232-233
            break;
        }
        j += 1;
    }
    
    let q_end_pos = q_end;
    let s_end_pos = s_start + (q_end_pos - q_start);
    
    (
        q_start,
        q_end_pos,
        s_start,
        s_end_pos,
        score,
    )
}

