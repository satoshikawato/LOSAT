//! Two-hit ungapped extension algorithm
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c

use super::get_score;
use super::ungapped::extend_hit_ungapped;
// EXPERIMENT (LOSAT_X_TBXUNGAP): x_two_hit_fast.rs.
use super::x_two_hit_fast;

// Import the internal extension functions from ungapped module
// These are used for left/right extension
#[inline(always)]
fn extend_left_ungapped_scalar(
    q_seq: &[u8],
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize) {
    // Use the public extend_hit_ungapped to get left extension behavior
    // We need to replicate the left extension logic here since it's internal
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
#[cfg(all(
    target_arch = "wasm32",
    target_feature = "simd128",
    not(feature = "tblastx-wasm-scalar")
))]
#[target_feature(enable = "simd128")]
unsafe fn extend_left_ungapped_simd128(
    q_seq: &[u8],
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize) {
    use std::arch::wasm32::*;

    let max_left = q_off.min(s_off);

    let mut score = initial_score;
    let mut max_score = initial_score;
    let mut left_disp = 0usize;
    let mut i = 0usize;

    while i + 16 <= max_left {
        let q_ptr = q_seq.as_ptr().add(q_off - i - 16);
        let s_ptr = s_seq.as_ptr().add(s_off - i - 16);
        // Safety: q_ptr/s_ptr are in-bounds by loop condition.
        let qv = v128_load(q_ptr as *const v128);
        let sv = v128_load(s_ptr as *const v128);
        let mut pos = i;
        macro_rules! step_lane_rev {
            ($lane:expr, $qv:ident, $sv:ident, $pos:ident) => {{
                let q = i8x16_extract_lane::<$lane>($qv) as u8;
                let s = i8x16_extract_lane::<$lane>($sv) as u8;
                score += get_score(q, s);
                let ext_len = $pos + 1;
                if score > max_score {
                    max_score = score;
                    left_disp = ext_len;
                }
                if (max_score - score) >= x_drop {
                    return (max_score, left_disp);
                }
                $pos += 1;
            }};
        }
        step_lane_rev!(15, qv, sv, pos);
        step_lane_rev!(14, qv, sv, pos);
        step_lane_rev!(13, qv, sv, pos);
        step_lane_rev!(12, qv, sv, pos);
        step_lane_rev!(11, qv, sv, pos);
        step_lane_rev!(10, qv, sv, pos);
        step_lane_rev!(9, qv, sv, pos);
        step_lane_rev!(8, qv, sv, pos);
        step_lane_rev!(7, qv, sv, pos);
        step_lane_rev!(6, qv, sv, pos);
        step_lane_rev!(5, qv, sv, pos);
        step_lane_rev!(4, qv, sv, pos);
        step_lane_rev!(3, qv, sv, pos);
        step_lane_rev!(2, qv, sv, pos);
        step_lane_rev!(1, qv, sv, pos);
        step_lane_rev!(0, qv, sv, pos);
        i += 16;
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
    #[cfg(all(
        target_arch = "wasm32",
        target_feature = "simd128",
        not(feature = "tblastx-wasm-scalar")
    ))]
    {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:886-921
        // ```c
        // for (i = n; i >= 0; i--) {
        //    score += matrix[q[i]][s[i]];
        //    if (score > maxscore) { ... }
        //    if ((maxscore - score) >= dropoff) break;
        // }
        // ```
        // The `tblastx-wasm-scalar` feature is diagnostic-only: it keeps this
        // exact scalar NCBI control flow available for native-vs-Wasm parity
        // isolation while leaving production SIMD builds unchanged.
        unsafe {
            return extend_left_ungapped_simd128(q_seq, s_seq, q_off, s_off, initial_score, x_drop);
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
#[cfg(all(
    target_arch = "wasm32",
    target_feature = "simd128",
    not(feature = "tblastx-wasm-scalar")
))]
#[target_feature(enable = "simd128")]
unsafe fn extend_right_ungapped_simd128(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    initial_score: i32,
    x_drop: i32,
) -> (i32, usize, usize) {
    use std::arch::wasm32::*;

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
        // Safety: q_ptr/s_ptr are in-bounds by loop condition.
        let qv = v128_load(q_ptr as *const v128);
        let sv = v128_load(s_ptr as *const v128);
        let mut pos = base;
        macro_rules! step_lane {
            ($lane:expr, $qv:ident, $sv:ident, $pos:ident) => {{
                let q = i8x16_extract_lane::<$lane>($qv) as u8;
                let s = i8x16_extract_lane::<$lane>($sv) as u8;
                score += get_score(q, s);
                if score > max_score {
                    max_score = score;
                    right_disp = $pos + 1;
                }
                if score <= 0 || (max_score - score) >= x_drop {
                    return (max_score, right_disp, $pos);
                }
                $pos += 1;
            }};
        }
        step_lane!(0, qv, sv, pos);
        step_lane!(1, qv, sv, pos);
        step_lane!(2, qv, sv, pos);
        step_lane!(3, qv, sv, pos);
        step_lane!(4, qv, sv, pos);
        step_lane!(5, qv, sv, pos);
        step_lane!(6, qv, sv, pos);
        step_lane!(7, qv, sv, pos);
        step_lane!(8, qv, sv, pos);
        step_lane!(9, qv, sv, pos);
        step_lane!(10, qv, sv, pos);
        step_lane!(11, qv, sv, pos);
        step_lane!(12, qv, sv, pos);
        step_lane!(13, qv, sv, pos);
        step_lane!(14, qv, sv, pos);
        step_lane!(15, qv, sv, pos);
        j += 16;
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
    #[cfg(all(
        target_arch = "wasm32",
        target_feature = "simd128",
        not(feature = "tblastx-wasm-scalar")
    ))]
    {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:846-866
        // ```c
        // for (i = 0; i < n; i++) {
        //    score += matrix[q[i]][s[i]];
        //    if (score > maxscore) { ... }
        //    if (score <= 0 || (maxscore - score) >= dropoff) break;
        // }
        // ```
        // The `tblastx-wasm-scalar` feature is diagnostic-only: it keeps this
        // exact scalar NCBI control flow available for native-vs-Wasm parity
        // isolation while leaving production SIMD builds unchanged.
        unsafe {
            return extend_right_ungapped_simd128(
                q_seq,
                s_seq,
                q_start,
                s_start,
                initial_score,
                x_drop,
            );
        }
    }
    extend_right_ungapped_scalar(q_seq, s_seq, q_start, s_start, initial_score, x_drop)
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
// EXPERIMENT (LOSAT_X_TBXUNGAP / LOSAT_X_TBXUNGAPSHADOW): the same values from the table walks of
// x_two_hit_fast.rs, which hand back to `extend_hit_two_hit_reference` what they do not cover; the
// debug output (`debug_ext`) stays with the reference.
#[inline]
pub fn extend_hit_two_hit(
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
    debug_ext: bool,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    if !debug_ext {
        match x_two_hit_fast::mode() {
            0 => {}
            1 => {
                if let Some(r) = x_two_hit_fast::two_hit(
                    q_seq,
                    s_seq,
                    s_left_off,
                    s_right_off,
                    q_right_off,
                    x_drop,
                ) {
                    return r;
                }
            }
            _ => {
                return x_two_hit_fast::shadow(
                    q_seq,
                    s_seq,
                    s_left_off,
                    s_right_off,
                    q_right_off,
                    x_drop,
                )
            }
        }
    }
    extend_hit_two_hit_reference(
        q_seq,
        s_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
        debug_ext,
    )
}

/// EXPERIMENT (LOSAT_X_TBXUNGAPSENT): `extend_hit_two_hit(&q_full[1..q_full.len() - 1],
/// &s_full[1..s_full.len() - 1], s_left_off, s_right_off, q_right_off, x_drop, false)` with the
/// walks that the sentinels of the two buffers stop (x_two_hit_fast.rs).
///
/// SAFETY: `q_full` is the buffer of a context of a batch `x_two_hit_fast::sent_contexts`
/// accepted, and `s_full` a buffer `x_two_hit_fast::sent_subject` returned (not empty).
// NCBI reference: aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit`.
#[inline]
pub unsafe fn extend_hit_two_hit_sent(
    q_full: &[u8],
    s_full: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    // SAFETY: the caller's.
    match unsafe {
        x_two_hit_fast::two_hit_sent(q_full, s_full, s_left_off, s_right_off, q_right_off, x_drop)
    } {
        Some(r) => r,
        None => extend_hit_two_hit_sent_refused(
            q_full,
            s_full,
            s_left_off,
            s_right_off,
            q_right_off,
            x_drop,
        ),
    }
}

/// The calls `extend_hit_two_hit_sent` leaves to `extend_hit_two_hit` (an X-drop outside
/// 1..=127, a word that does not fit), out of the callers' code.
// NCBI reference: aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit`.
#[cold]
#[inline(never)]
fn extend_hit_two_hit_sent_refused(
    q_full: &[u8],
    s_full: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    extend_hit_two_hit(
        &q_full[1..q_full.len() - 1],
        &s_full[1..s_full.len() - 1],
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
        false,
    )
}

/// `extend_hit_two_hit` as written before LOSAT_X_TBXUNGAP (the path with the switch off, and
/// the debug output).
pub(super) fn extend_hit_two_hit_reference(
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,  // Position of first hit (L) in subject
    s_right_off: usize, // Position of second hit (R) in subject
    q_right_off: usize, // Position of second hit (R) in query
    x_drop: i32,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1089-1096
    // ```c
    // s_BlastAaExtendTwoHit(..., Int4 q_right_off, Int4 dropoff, ...)
    // ```
    // Existing LOSAT diagnostics only; supplied by the current search owner.
    debug_ext: bool,
) -> (usize, usize, usize, usize, i32, bool, usize) {
    // DEBUG: Print sequence segment for tracing
    static DEBUG_PRINTED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

    let k_size = 3;

    // First, find the best scoring position within the word at R
    // NCBI BLAST assumes caller provides valid coordinates within bounds.
    // The lookup table and offset generation guarantee valid positions.
    // Reference: aa_ungapped.c does NOT have bounds checks in word scanning loop.
    let mut score = 0i32;
    let mut seed_score = 0i32;
    let mut right_d = 0usize;

    for i in 0..k_size {
        // SAFETY: Lookup table guarantees q_right_off and s_right_off have at least k_size valid positions
        let q_char = unsafe { *q_seq.get_unchecked(q_right_off + i) };
        let s_char = unsafe { *s_seq.get_unchecked(s_right_off + i) };
        score += get_score(q_char, s_char);

        if score > seed_score {
            seed_score = score;
            right_d = i + 1; // Position is one beyond the end of the word
        }
    }

    let q_right_off = q_right_off + right_d;
    let s_right_off = s_right_off + right_d;

    // Extend to the left from R, trying to reach L
    // Note: NCBI BLAST does NOT explicitly check for stop codons in s_BlastAaExtendLeft.
    // Stop codons have very negative scores in BLOSUM62 (-4 to -5), so X-drop handles them.
    // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:886-921
    let (left_score, left_disp) =
        extend_left_ungapped(q_seq, s_seq, q_right_off, s_right_off, 0, x_drop);

    // Check if left extension reached the first hit
    // The distance from R to L is (s_right_off - s_left_off)
    // We need left_disp >= this distance for the extension to reach L
    let distance_to_first_hit = s_right_off.saturating_sub(s_left_off);
    let reached_first_hit = left_disp >= distance_to_first_hit;

    // Only extend to the right if left extension reached the first hit (NCBI BLAST behavior)
    let mut right_disp = 0usize;
    let mut right_score = 0i32;
    let mut right_extended = false;
    let mut s_last_off = s_right_off; // Default to current position if no right extension

    if reached_first_hit {
        // NCBI BLAST: right_extend is set to TRUE when we enter this block,
        // regardless of how many characters are actually extended.
        // Reference: aa_ungapped.c:1140
        right_extended = true;

        // NCBI BLAST extends right starting from left_score as initial score
        let (new_max_score, new_right_disp, j) =
            extend_right_ungapped(q_seq, s_seq, q_right_off, s_right_off, left_score, x_drop);
        right_score = new_max_score;
        right_disp = new_right_disp;
        // s_last_off is the rightmost subject position scanned (NCBI BLAST style)
        s_last_off = s_right_off + j;
    }

    // NCBI BLAST returns MAX(left_score, right_score)
    // Reference: aa_ungapped.c:1157
    let max_score_total = left_score.max(right_score);

    let q_start = q_right_off - left_disp;
    let q_end = q_right_off + right_disp;
    let s_start = s_right_off - left_disp;
    let s_end = s_right_off + right_disp;

    // DEBUG: Print alignment details for first few high-scoring extensions
    if debug_ext
        && max_score_total >= 100
        && !DEBUG_PRINTED.swap(true, std::sync::atomic::Ordering::Relaxed)
    {
        eprintln!(
            "[DEBUG_EXT] left_score={} right_score={} max_score_total={}",
            left_score, right_score, max_score_total
        );
        eprintln!(
            "[DEBUG_EXT] q_start={} q_end={} s_start={} s_end={} score={}",
            q_start, q_end, s_start, s_end, max_score_total
        );
        eprintln!(
            "[DEBUG_EXT] left_disp={} right_disp={} x_drop={}",
            left_disp, right_disp, x_drop
        );
        // Print position-by-position scores
        let align_len = q_end - q_start;
        eprintln!("[DEBUG_EXT] alignment_len={}", align_len);
        let mut total = 0i32;
        for i in 0..align_len.min(100) {
            let q_char = q_seq[q_start + i];
            let s_char = s_seq[s_start + i];
            let sc = get_score(q_char, s_char);
            total += sc;
            if i < 10 || i >= align_len - 5 {
                eprintln!(
                    "[DEBUG_EXT] pos={:2} q={:2} s={:2} score={:3} cumul={:4}",
                    i, q_char, s_char, sc, total
                );
            }
        }
        eprintln!(
            "[DEBUG_EXT] manual_total={} reported_score={}",
            total, max_score_total
        );
    }

    (
        q_start,
        q_end,
        s_start,
        s_end,
        max_score_total,
        right_extended,
        s_last_off,
    )
}

// No NCBI counterpart (LOSAT_X_STATS): the cells one call of `extend_hit_two_hit` scores, for
// the statistics of the unit path (tblastx/x_unit_probe.rs). It repeats the call's scalar loops
// (the same sums and the same stops, aa_ungapped.c:846-866 and 886-921) and counts their steps;
// the search never calls it. The test below checks that it returns what the call returns.
/// What one call of `extend_hit_two_hit` scored, and its result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct XTwoHitCells {
    /// Cells of the left extension (from the word's best end towards the first hit).
    pub(crate) left: usize,
    /// Cells of the right extension; 0 when it did not run.
    pub(crate) right: usize,
    /// First subject position scored (the word's cells included).
    pub(crate) s_lo: usize,
    /// One past the last subject position scored.
    pub(crate) s_hi: usize,
    /// The seven values `extend_hit_two_hit` returns for the same arguments.
    pub(crate) result: (usize, usize, usize, usize, i32, bool, usize),
}

/// The cells of `extend_hit_two_hit(q_seq, s_seq, s_left_off, s_right_off, q_right_off, x_drop, false)`.
pub(crate) fn x_two_hit_cells(
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> XTwoHitCells {
    let k_size = 3;
    let (mut score, mut seed_score, mut right_d) = (0i32, 0i32, 0usize);
    for i in 0..k_size {
        score += get_score(q_seq[q_right_off + i], s_seq[s_right_off + i]);
        if score > seed_score {
            seed_score = score;
            right_d = i + 1;
        }
    }
    let (q_r, s_r) = (q_right_off + right_d, s_right_off + right_d);

    let max_left = q_r.min(s_r);
    let (mut score, mut left_score, mut left_disp, mut left) = (0i32, 0i32, 0usize, 0usize);
    while left < max_left {
        score += get_score(q_seq[q_r - 1 - left], s_seq[s_r - 1 - left]);
        left += 1;
        if score > left_score {
            left_score = score;
            left_disp = left;
        }
        if (left_score - score) >= x_drop {
            break;
        }
    }

    let reached_first_hit = left_disp >= s_r.saturating_sub(s_left_off);
    let (mut right, mut right_disp, mut right_score, mut s_last_off) = (0usize, 0usize, 0i32, s_r);
    if reached_first_hit {
        let (mut score, mut max_score, mut j) = (left_score, left_score, 0usize);
        while q_r + j < q_seq.len() && s_r + j < s_seq.len() {
            score += get_score(q_seq[q_r + j], s_seq[s_r + j]);
            right += 1;
            if score > max_score {
                max_score = score;
                right_disp = j + 1;
            }
            if score <= 0 || (max_score - score) >= x_drop {
                break;
            }
            j += 1;
        }
        right_score = max_score;
        s_last_off = s_r + j;
    }

    XTwoHitCells {
        left,
        right,
        s_lo: s_right_off.min(s_r - left),
        s_hi: (s_right_off + k_size).max(s_r + right),
        result: (
            q_r - left_disp,
            q_r + right_disp,
            s_r - left_disp,
            s_r + right_disp,
            left_score.max(right_score),
            reached_first_hit,
            s_last_off,
        ),
    }
}

#[cfg(test)]
mod x_cells_tests {
    use super::*;

    // No NCBI counterpart (a statistics helper). The counting copy returns what the extension
    // returns, on random NCBIstdaa sequences, offsets and X-drops.
    #[test]
    fn the_counting_copy_returns_what_the_extension_returns() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move |n: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % n
        };
        for _ in 0..20_000 {
            let q_len = 4 + next(120) as usize;
            let s_len = 4 + next(120) as usize;
            // Mostly a few residues, so that runs of positive scores are common.
            let alphabet = 1 + next(24);
            let q: Vec<u8> = (0..q_len).map(|_| 1 + next(alphabet) as u8).collect();
            let s: Vec<u8> = (0..s_len).map(|_| 1 + next(alphabet) as u8).collect();
            let q_right = next((q_len - 3) as u64 + 1) as usize;
            let s_right = next((s_len - 3) as u64 + 1) as usize;
            let s_left = s_right.saturating_sub(next(50) as usize);
            let x_drop = 1 + next(40) as i32;
            let cells = x_two_hit_cells(&q, &s, s_left, s_right, q_right, x_drop);
            let result = extend_hit_two_hit(&q, &s, s_left, s_right, q_right, x_drop, false);
            assert_eq!(cells.result, result);
            if !result.5 {
                assert_eq!(cells.right, 0);
            }
            assert!(cells.s_lo <= s_right && cells.s_hi >= s_right + 3 && cells.s_hi <= s_len);
        }
    }
}
