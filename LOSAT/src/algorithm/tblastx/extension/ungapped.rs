//! Ungapped extension algorithms with SIMD optimization
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c

use std::sync::OnceLock;
use super::get_score;

// ---------------------------------------------------------------------------
// Ungapped extension SIMD helpers (strict NCBI parity)
// ---------------------------------------------------------------------------

// NCBI BLAST reference (c++/src/algo/blast/core/aa_ungapped.c:846-866):
//   for (i = 0; i < n; i++) {
//       score += matrix[q[i]][s[i]];
//       if (score > maxscore) {
//           maxscore = score;
//           best_i = i;
//       }
//       if (score <= 0 || (maxscore - score) >= dropoff)
//           break;
//   }
//   *length = best_i + 1;
//   *s_last_off = s_off + i;
//
// We keep the exact control flow/termination conditions, but batch the score
// lookup for consecutive residue pairs using AVX2 gather. Accumulation and
// X-drop checks remain strictly sequential to preserve 1-bit parity.

static UNGAPPED_SCORE_TABLE_32: OnceLock<[i32; 32 * 32]> = OnceLock::new();

#[inline(always)]
fn ungapped_score_table_32() -> &'static [i32; 32 * 32] {
    UNGAPPED_SCORE_TABLE_32.get_or_init(|| {
        let mut t = [0i32; 32 * 32];
        for a in 0..32u8 {
            for b in 0..32u8 {
                t[(a as usize) * 32 + (b as usize)] = get_score(a, b);
            }
        }
        t
    })
}

// NCBI BLAST reference (c++/src/algo/blast/core/aa_ungapped.c:886-921):
//   n = MIN(s_off, q_off);
//   best_i = n + 1;
//   s = subject->sequence + s_off - n;
//   q = query->sequence + q_off - n;
//   for (i = n; i >= 0; i--) {
//       score += matrix[q[i]][s[i]];
//       if (score > maxscore) {
//           maxscore = score;
//           best_i = i;
//       }
//       if ((maxscore - score) >= dropoff)
//           break;
//   }
//   *length = n - best_i + 1;
//   return maxscore;
//
// LOSAT form uses (q_off-1-i, s_off-1-i) order; we keep that exact control flow,
// but batch the score lookup with AVX2 gather. Accumulation and X-drop checks
// remain strictly sequential to preserve 1-bit parity.

#[inline(always)]
fn extend_left_ungapped_scalar(
    q_seq: &[u8],
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize) {
    let max_left = q_off.min(s_off);

    let mut score = initial_score;
    let mut max_score = initial_score;
    let mut left_disp = 0usize;
    let mut i = 0usize;

    while i < max_left {
        let q_char = unsafe { *q_seq.get_unchecked(q_off - 1 - i) };
        let s_char = unsafe { *s_seq.get_unchecked(s_off - 1 - i) };

        score += get_score(q_char, s_char);

        if score > max_score {
            max_score = score;
            left_disp = i + 1;
        }

        if (max_score - score) >= x_drop {
            break;
        }
        i += 1;
    }

    (max_score, left_disp)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn extend_left_ungapped_avx2(
    q_seq: &[u8],
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize) {
    let table = ungapped_score_table_32().as_ptr();
    let max_left = q_off.min(s_off);

    let mut score = initial_score;
    let mut max_score = initial_score;
    let mut left_disp = 0usize;
    let mut i = 0usize;

    // Process blocks of 8 residues. We must apply scores in the same order as the
    // scalar loop: nearest first (q_off-1), then further left.
    while i + 8 <= max_left {
        let q_ptr = q_seq.as_ptr().add(q_off - i - 8);
        let s_ptr = s_seq.as_ptr().add(s_off - i - 8);
        let scores = gather_scores8_avx2(q_ptr, s_ptr, table); // [far..near]

        for k in 0..8usize {
            // Apply near..far: scores[7], scores[6], ...
            score += scores[7 - k];
            let ext_len = i + k + 1;

            if score > max_score {
                max_score = score;
                left_disp = ext_len;
            }

            if (max_score - score) >= x_drop {
                return (max_score, left_disp);
            }
        }

        i += 8;
    }

    while i < max_left {
        let q_char = *q_seq.get_unchecked(q_off - 1 - i);
        let s_char = *s_seq.get_unchecked(s_off - 1 - i);

        score += get_score(q_char, s_char);
        let ext_len = i + 1;

        if score > max_score {
            max_score = score;
            left_disp = ext_len;
        }

        if (max_score - score) >= x_drop {
            break;
        }
        i += 1;
    }

    (max_score, left_disp)
}

#[inline(always)]
fn extend_left_ungapped(
    q_seq: &[u8],
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize) {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                return extend_left_ungapped_avx2(q_seq, s_seq, q_off, s_off, initial_score, x_drop);
            }
        }
    }
    extend_left_ungapped_scalar(q_seq, s_seq, q_off, s_off, initial_score, x_drop)
}

#[inline(always)]
fn extend_right_ungapped_scalar(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize, usize) {
    let q_limit = q_seq.len();
    let s_limit = s_seq.len();

    let mut score = initial_score;
    let mut max_score = initial_score;
    let mut right_disp = 0usize;
    let mut j = 0usize;

    while q_start + j < q_limit && s_start + j < s_limit {
        let q_char = unsafe { *q_seq.get_unchecked(q_start + j) };
        let s_char = unsafe { *s_seq.get_unchecked(s_start + j) };

        score += get_score(q_char, s_char);

        if score > max_score {
            max_score = score;
            right_disp = j + 1;
        }

        if score <= 0 || (max_score - score) >= x_drop {
            break;
        }
        j += 1;
    }

    (max_score, right_disp, j)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn gather_scores8_avx2(q_ptr: *const u8, s_ptr: *const u8, table: *const i32) -> [i32; 8] {
    use std::arch::x86_64::*;

    // Load 8 bytes from each sequence.
    let q8 = _mm_loadl_epi64(q_ptr as *const __m128i);
    let s8 = _mm_loadl_epi64(s_ptr as *const __m128i);

    // Expand to 8x u32 lanes.
    let q32 = _mm256_cvtepu8_epi32(q8);
    let s32 = _mm256_cvtepu8_epi32(s8);

    // Index into 32x32 score table: (q<<5)|s.
    let idx = _mm256_or_si256(_mm256_slli_epi32(q32, 5), s32);

    // Gather 8 i32 scores.
    let v = _mm256_i32gather_epi32(table, idx, 4);

    let mut out = [0i32; 8];
    _mm256_storeu_si256(out.as_mut_ptr() as *mut __m256i, v);
    out
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn extend_right_ungapped_avx2(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize, usize) {
    let table = ungapped_score_table_32().as_ptr();
    let q_limit = q_seq.len();
    let s_limit = s_seq.len();
    let n = (q_limit - q_start).min(s_limit - s_start);

    let mut score = initial_score;
    let mut max_score = initial_score;
    let mut right_disp = 0usize;
    let mut j = 0usize;

    while j + 16 <= n {
        let base = j;
        let q_ptr = q_seq.as_ptr().add(q_start + base);
        let s_ptr = s_seq.as_ptr().add(s_start + base);

        let scores0 = gather_scores8_avx2(q_ptr, s_ptr, table);
        let scores1 = gather_scores8_avx2(q_ptr.add(8), s_ptr.add(8), table);

        for k in 0..8usize {
            let pos = base + k;
            score += scores0[k];
            if score > max_score {
                max_score = score;
                right_disp = pos + 1;
            }
            if score <= 0 || (max_score - score) >= x_drop {
                return (max_score, right_disp, pos);
            }
        }
        for k in 0..8usize {
            let pos = base + 8 + k;
            score += scores1[k];
            if score > max_score {
                max_score = score;
                right_disp = pos + 1;
            }
            if score <= 0 || (max_score - score) >= x_drop {
                return (max_score, right_disp, pos);
            }
        }

        j += 16;
    }

    while j + 8 <= n {
        let base = j;
        let q_ptr = q_seq.as_ptr().add(q_start + base);
        let s_ptr = s_seq.as_ptr().add(s_start + base);
        let scores = gather_scores8_avx2(q_ptr, s_ptr, table);

        for k in 0..8usize {
            let pos = base + k;
            score += scores[k];
            if score > max_score {
                max_score = score;
                right_disp = pos + 1;
            }
            if score <= 0 || (max_score - score) >= x_drop {
                return (max_score, right_disp, pos);
            }
        }

        j += 8;
    }

    while j < n {
        let q_char = *q_seq.get_unchecked(q_start + j);
        let s_char = *s_seq.get_unchecked(s_start + j);

        score += get_score(q_char, s_char);

        if score > max_score {
            max_score = score;
            right_disp = j + 1;
        }

        if score <= 0 || (max_score - score) >= x_drop {
            break;
        }
        j += 1;
    }

    (max_score, right_disp, j)
}

#[inline(always)]
fn extend_right_ungapped(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize, usize) {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                return extend_right_ungapped_avx2(q_seq, s_seq, q_start, s_start, initial_score, x_drop);
            }
        }
    }
    extend_right_ungapped_scalar(q_seq, s_seq, q_start, s_start, initial_score, x_drop)
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

    // NCBI BLAST assumes caller provides valid coordinates within bounds.
    // The lookup table and offset generation guarantee valid positions.
    // Reference: aa_ungapped.c does NOT have bounds checks in word scanning loop.
    for i in 0..k_size {
        // SAFETY: Lookup table guarantees q_pos and s_pos have at least k_size valid positions
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
    let (max_score, left_disp) =
        extend_left_ungapped(q_seq, s_seq, q_left_off, s_left_off, score, x_drop);

    // Step 3: Right extension from the best position
    let q_start_r = q_right_off + 1;
    let s_start_r = s_right_off + 1;
    let (max_score_total, right_disp, j) =
        extend_right_ungapped(q_seq, s_seq, q_start_r, s_start_r, max_score, x_drop);

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
