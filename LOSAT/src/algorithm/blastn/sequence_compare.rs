#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

/// SIMD-optimized 28-mer comparison
/// Optimized specifically for 28-byte comparisons used in megablast
/// Uses 16 bytes + 12 bytes with SSE2, or 32 bytes with AVX2
#[inline]
pub fn compare_28mer_simd(query: &[u8], subject: &[u8]) -> bool {
    const MER_LEN: usize = 28;
    
    if query.len() < MER_LEN || subject.len() < MER_LEN {
        return false;
    }

    #[cfg(target_arch = "x86_64")]
    {
        // Try AVX2 first (32 bytes at once, covers 28 bytes)
        if is_x86_feature_detected!("avx2") {
            unsafe {
                return compare_28mer_avx2(query, subject);
            }
        }
        
        // Fall back to SSE2 (16 bytes + 12 bytes)
        if is_x86_feature_detected!("sse2") {
            unsafe {
                return compare_28mer_sse2(query, subject);
            }
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if is_aarch64_feature_detected!("neon") {
            unsafe {
                return compare_28mer_neon(query, subject);
            }
        }
    }

    // Fallback to scalar comparison
    compare_28mer_scalar(query, subject)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn compare_28mer_avx2(query: &[u8], subject: &[u8]) -> bool {
    // Load 32 bytes (covers 28 bytes, with 4 extra bytes that we ignore)
    let q_chunk = _mm256_loadu_si256(query.as_ptr() as *const __m256i);
    let s_chunk = _mm256_loadu_si256(subject.as_ptr() as *const __m256i);
    let cmp = _mm256_cmpeq_epi8(q_chunk, s_chunk);
    
    // Get the movemask for all 32 bytes (returns 32-bit mask)
    let mask = _mm256_movemask_epi8(cmp);
    
    // Check if first 28 bytes match (bits 0-27 should all be set)
    // 0x0FFFFFFF = 28 bits set (0xFFFFFFFF >> 4)
    if (mask & 0x0FFFFFFF) != 0x0FFFFFFF {
        return false;
    }
    
    true
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn compare_28mer_sse2(query: &[u8], subject: &[u8]) -> bool {
    // First 16 bytes
    let q_chunk1 = _mm_loadu_si128(query.as_ptr() as *const __m128i);
    let s_chunk1 = _mm_loadu_si128(subject.as_ptr() as *const __m128i);
    let cmp1 = _mm_cmpeq_epi8(q_chunk1, s_chunk1);
    let mask1 = _mm_movemask_epi8(cmp1);
    
    if mask1 != 0xFFFF {
        return false;
    }
    
    // Remaining 12 bytes (bytes 16-27)
    // Load 16 bytes starting at offset 16, but only check first 12
    let q_chunk2 = _mm_loadu_si128(query.as_ptr().add(16) as *const __m128i);
    let s_chunk2 = _mm_loadu_si128(subject.as_ptr().add(16) as *const __m128i);
    let cmp2 = _mm_cmpeq_epi8(q_chunk2, s_chunk2);
    let mask2 = _mm_movemask_epi8(cmp2);
    
    // Check only the lower 12 bits (0x0FFF)
    if (mask2 & 0x0FFF) != 0x0FFF {
        return false;
    }
    
    true
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn compare_28mer_neon(query: &[u8], subject: &[u8]) -> bool {
    // First 16 bytes
    let q_chunk1 = vld1q_u8(query.as_ptr());
    let s_chunk1 = vld1q_u8(subject.as_ptr());
    let cmp1 = vceqq_u8(q_chunk1, s_chunk1);
    let all_eq1 = vminvq_u8(cmp1);
    
    if all_eq1 == 0 {
        return false;
    }
    
    // Remaining 12 bytes (bytes 16-27)
    // Load 16 bytes but only check first 12
    let q_chunk2 = vld1q_u8(query.as_ptr().add(16));
    let s_chunk2 = vld1q_u8(subject.as_ptr().add(16));
    let cmp2 = vceqq_u8(q_chunk2, s_chunk2);
    
    // Check first 12 bytes by extracting and comparing
    // Use vget_lane to check individual bytes, or use scalar for last 12 bytes
    // For simplicity, check the minimum value of the first 12 bytes
    let cmp2_low = vget_low_u8(cmp2);
    let all_eq2_low = vminv_u8(cmp2_low);
    
    if all_eq2_low == 0 {
        return false;
    }
    
    // Check bytes 24-27 (4 bytes from high half)
    let cmp2_high = vget_high_u8(cmp2);
    let all_eq2_high = vminv_u8(cmp2_high);
    
    if all_eq2_high == 0 {
        return false;
    }
    
    true
}

/// Scalar fallback for 28-mer comparison
#[inline]
pub fn compare_28mer_scalar(query: &[u8], subject: &[u8]) -> bool {
    const MER_LEN: usize = 28;
    
    if query.len() < MER_LEN || subject.len() < MER_LEN {
        return false;
    }
    
    // Unrolled comparison for 28 bytes
    // Compare in chunks of 8 bytes (3 chunks) + 4 bytes
    if query[0..8] != subject[0..8] {
        return false;
    }
    if query[8..16] != subject[8..16] {
        return false;
    }
    if query[16..24] != subject[16..24] {
        return false;
    }
    if query[24..28] != subject[24..28] {
        return false;
    }
    
    true
}

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
            unsafe {
                return compare_sequences_sse2(query, subject, len);
            }
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if is_aarch64_feature_detected!("neon") {
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
#[inline(always)]
pub fn find_first_mismatch_ex(
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
            if unsafe { *s1.get_unchecked(i1) != *s2.get_unchecked(i2) } {
                break;
            }
            count += 1;
        }
    } else {
        let s1 = &seq1[start1..];
        let s2 = &seq2[start2..];
        while count < max_len {
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


