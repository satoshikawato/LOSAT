#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:351-357
// while (seq1_index < len1 && seq2_index < len2 &&
//        seq1[seq1_index] < 4 &&
//        seq1[seq1_index] == seq2[seq2_index]) {
//     ++seq1_index;
//     ++seq2_index;
// }
// This Rust-only import makes the existing AArch64 runtime dispatch compile;
// it does not change NCBI's exact-match comparison condition.
#[cfg(target_arch = "aarch64")]
use std::arch::is_aarch64_feature_detected;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// SIMD-optimized sequence comparison for word_length match check
/// Compares two sequences of equal length using SIMD instructions when available
/// Falls back to scalar comparison if SIMD is not available
#[inline]
pub fn compare_sequences_simd(query: &[u8], subject: &[u8], len: usize) -> bool {
    if query.len() < len || subject.len() < len {
        return false;
    }

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("sse2") {
            // SAFETY: the CPU has SSE2 (checked on the line above); the 16-byte loads stay below
            // `len` (`i + 16 <= len`), and the check at the top keeps `len` within both slices.
            unsafe {
                return compare_sequences_sse2(query, subject, len);
            }
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if is_aarch64_feature_detected!("neon") {
            // SAFETY: the CPU has NEON (checked on the line above); the 16-byte loads stay below
            // `len` (`i + 16 <= len`), and the check at the top keeps `len` within both slices.
            unsafe {
                return compare_sequences_neon(query, subject, len);
            }
        }
    }

    // Fallback to scalar comparison
    compare_sequences_scalar(query, subject, len)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn compare_sequences_sse2(query: &[u8], subject: &[u8], len: usize) -> bool {
    let mut i = 0;

    // Process 16 bytes at a time using SSE2
    while i + 16 <= len {
        let q_chunk = _mm_loadu_si128(query.as_ptr().add(i) as *const __m128i);
        let s_chunk = _mm_loadu_si128(subject.as_ptr().add(i) as *const __m128i);
        let cmp = _mm_cmpeq_epi8(q_chunk, s_chunk);
        let mask = _mm_movemask_epi8(cmp);

        // If any bytes don't match, mask will have 0 bits
        if mask != 0xFFFF {
            return false;
        }
        i += 16;
    }

    // Handle remaining bytes with scalar comparison
    while i < len {
        if query[i] != subject[i] {
            return false;
        }
        i += 1;
    }

    true
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn compare_sequences_neon(query: &[u8], subject: &[u8], len: usize) -> bool {
    let mut i = 0;

    // Process 16 bytes at a time using NEON
    while i + 16 <= len {
        let q_chunk = vld1q_u8(query.as_ptr().add(i));
        let s_chunk = vld1q_u8(subject.as_ptr().add(i));
        let cmp = vceqq_u8(q_chunk, s_chunk);
        let all_eq = vminvq_u8(cmp);

        if all_eq == 0 {
            return false;
        }
        i += 16;
    }

    // Handle remaining bytes with scalar comparison
    while i < len {
        if query[i] != subject[i] {
            return false;
        }
        i += 1;
    }

    true
}

/// Scalar fallback for sequence comparison
#[inline]
pub fn compare_sequences_scalar(query: &[u8], subject: &[u8], len: usize) -> bool {
    if query.len() < len || subject.len() < len {
        return false;
    }

    // Unrolled loop for better performance
    let mut i = 0;
    while i + 8 <= len {
        if query[i] != subject[i]
            || query[i + 1] != subject[i + 1]
            || query[i + 2] != subject[i + 2]
            || query[i + 3] != subject[i + 3]
            || query[i + 4] != subject[i + 4]
            || query[i + 5] != subject[i + 5]
            || query[i + 6] != subject[i + 6]
            || query[i + 7] != subject[i + 7]
        {
            return false;
        }
        i += 8;
    }

    // Handle remaining bytes
    while i < len {
        if query[i] != subject[i] {
            return false;
        }
        i += 1;
    }

    true
}

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
