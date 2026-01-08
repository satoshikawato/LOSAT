//! Offset pair storage and SIMD copy functions
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:99-115

/// Query/subject offset pair (NCBI `BlastOffsetPair` equivalent for the hot path).
#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct OffsetPair {
    pub q_off: i32,
    pub s_off: i32,
}

// ---------------------------------------------------------------------------
// SIMD copy helpers (k-mer scan hot path)
// ---------------------------------------------------------------------------

// NCBI BLAST reference (c++/src/algo/blast/core/blast_aascan.c:99-115):
//   /* copy the hits. */
//   {
//       Int4 i;
//       Int4 s_off = (Int4)(s - subject->sequence);
//       for (i = 0; i < numhits; i++) {
//           offset_pairs[i + totalhits].qs_offsets.q_off = src[i];
//           offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
//       }
//   }
//
// LOSAT performs the same copy, but can accelerate the (overflow) case where
// `numhits > AA_HITS_PER_CELL` by packing each (q_off, s_off) pair into u64 and
// storing 4 pairs at a time with AVX2. Output order is unchanged.

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn copy_offset_pairs_overflow_avx2(
    src: *const i32,
    dest: *mut OffsetPair,
    numhits: usize,
    s_off: i32,
) {
    use std::arch::x86_64::*;

    // Pack as: [low32=q_off | high32=s_off] (little-endian)
    let s_off_u64 = (s_off as u32 as u64) << 32;
    let s_vec = _mm256_set1_epi64x(s_off_u64 as i64);
    let low32_mask = _mm256_set1_epi64x(0xFFFF_FFFFu64 as i64);

    let mut i = 0usize;

    // 8 pairs per iteration (2x 4-pair stores)
    while i + 8 <= numhits {
        let q0 = _mm_loadu_si128(src.add(i) as *const __m128i);
        let q1 = _mm_loadu_si128(src.add(i + 4) as *const __m128i);

        // cvtepi32_epi64 sign-extends; mask to keep only low 32 bits.
        let q0_64 = _mm256_and_si256(_mm256_cvtepi32_epi64(q0), low32_mask);
        let q1_64 = _mm256_and_si256(_mm256_cvtepi32_epi64(q1), low32_mask);

        let p0 = _mm256_or_si256(q0_64, s_vec);
        let p1 = _mm256_or_si256(q1_64, s_vec);

        _mm256_storeu_si256(dest.add(i) as *mut __m256i, p0);
        _mm256_storeu_si256(dest.add(i + 4) as *mut __m256i, p1);

        i += 8;
    }

    // 4 pairs remainder
    if i + 4 <= numhits {
        let q = _mm_loadu_si128(src.add(i) as *const __m128i);
        let q64 = _mm256_and_si256(_mm256_cvtepi32_epi64(q), low32_mask);
        let p = _mm256_or_si256(q64, s_vec);
        _mm256_storeu_si256(dest.add(i) as *mut __m256i, p);
        i += 4;
    }

    // Tail
    while i < numhits {
        *dest.add(i) = OffsetPair { q_off: *src.add(i), s_off };
        i += 1;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
pub unsafe fn copy_offset_pairs_overflow_sse2(
    src: *const i32,
    dest: *mut OffsetPair,
    numhits: usize,
    s_off: i32,
) {
    use std::arch::x86_64::*;

    // Interleave q_off with constant s_off in 32-bit lanes:
    // q = [q3 q2 q1 q0]
    // s = [ s  s  s  s ]
    // lo = unpacklo(q,s) => [q1 s q0 s] (two OffsetPair)
    // hi = unpackhi(q,s) => [q3 s q2 s] (two OffsetPair)
    let s_vec = _mm_set1_epi32(s_off);

    let mut i = 0usize;
    while i + 4 <= numhits {
        let q = _mm_loadu_si128(src.add(i) as *const __m128i);
        let lo = _mm_unpacklo_epi32(q, s_vec);
        let hi = _mm_unpackhi_epi32(q, s_vec);

        _mm_storeu_si128(dest.add(i) as *mut __m128i, lo);
        _mm_storeu_si128(dest.add(i + 2) as *mut __m128i, hi);

        i += 4;
    }

    while i < numhits {
        *dest.add(i) = OffsetPair { q_off: *src.add(i), s_off };
        i += 1;
    }
}
