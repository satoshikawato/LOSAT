/// Find the first mismatch between two sequences, starting from given positions.
/// If `reverse` is true, sequences are accessed from the end (for left extension).
/// Returns the number of matching characters before the first mismatch.
/// Forward (`reverse == false`), `len1` and `len2` must be the slices' lengths; the crate's
/// callers pass them (not `pub`: outside code could pass longer ones).
#[inline(always)]
pub(crate) fn find_first_mismatch_ex(
    seq1: &[u8],
    seq2: &[u8],
    len1: usize,
    len2: usize,
    start1: usize,
    start2: usize,
    reverse: bool,
) -> usize {
    let max_len = (len1 - start1).min(len2 - start2);
    let mut count = 0;

    if reverse {
        // Access sequences from the end
        let s1 = &seq1[..len1 - start1];
        let s2 = &seq2[..len2 - start2];
        let mut i1 = s1.len();
        let mut i2 = s2.len();

        while count < max_len && i1 > 0 && i2 > 0 {
            i1 -= 1;
            i2 -= 1;
            // SAFETY: `i1 < s1.len()` and `i2 < s2.len()`: both start at their slice's length, and
            // the loop guard (`i1 > 0 && i2 > 0`) lets them drop by one just before this read.
            if unsafe { *s1.get_unchecked(i1) != *s2.get_unchecked(i2) } {
                break;
            }
            count += 1;
        }
    } else {
        let s1 = &seq1[start1..];
        let s2 = &seq2[start2..];
        while count < max_len {
            // SAFETY: `count < max_len <= len1 - start1 <= s1.len()` (and the same for `s2`): every
            // forward call in the crate passes the slices' own lengths as `len1` and `len2`
            // (`find_first_mismatch`, and `greedy_align_one_direction` through the greedy
            // walks); the slicing above checks `start1` and `start2`.
            if unsafe { *s1.get_unchecked(count) != *s2.get_unchecked(count) } {
                break;
            }
            count += 1;
        }
    }

    count
}

#[inline(always)]
pub fn find_first_mismatch(seq1: &[u8], seq2: &[u8], start1: usize, start2: usize) -> usize {
    find_first_mismatch_ex(seq1, seq2, seq1.len(), seq2.len(), start1, start2, false)
}
