//! Line FN s9 probe (`LOSAT_X_FNUNG`, a bit mask; not a switch of the strict set): the loops of
//! the ungapped extensions of blastn in the shape that line FR counted (one induction variable,
//! the bounds checked once per call, the subject read through a 64-bit window).
//!
//! - bit 0: `extend_hit_ungapped_approx_ncbi` runs `approx_blocks`: each 4-letter query block is
//!   one 32-bit load folded into NCBI's overlapping BLASTNA byte by two shifts and two ORs
//!   (`four_base_byte_of_word`), and each direction is a counted loop over a slice checked once.
//! - bit 2: `extend_hit_ungapped_exact_counted` runs `exact_window`: the subject bases come out of
//!   a 64-bit big-endian window of the packed bytes with a constant 2-bit shift per letter
//!   (one variable shift per 29 to 32 letters), instead of a byte load and a variable shift each.
//!
//! The steps, the C tests and the C updates on the same `i32` values are those of the loops they
//! replace, so the four results are the same (unit tests below against the loop-for-loop port and
//! against today's approximate loops).

use super::extension::UngappedData;
use crate::core::blast_encoding::COMPRESSION_RATIO;

/// BLASTNA alphabet size (`blast_encoding.h:91-92`).
const BLASTNA_SIZE: usize = 16;

/// The probe mask, read once.
#[inline(always)]
pub(crate) fn probe() -> u8 {
    use std::sync::OnceLock;
    static MASK: OnceLock<u8> = OnceLock::new();
    *MASK.get_or_init(|| {
        std::env::var("LOSAT_X_FNUNG")
            .ok()
            .and_then(|v| v.trim().parse::<u8>().ok())
            .unwrap_or(0)
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:294,324
/// The overlapping BLASTNA byte of `query_four_base_byte` from the four letters read as one
/// little-endian word `w` (letter `start` in its low byte). With `v = (w << 6) | (w >> 4)`, the low
/// byte of `v` is `(q0 << 6) | (q1 << 4)` (bits 4..8 of `q0` are 0) and bits 20..28 of `v` are
/// `(q2 << 2) | q3` (bits 6..8 of `q1` land on bits 20..22 and are 0): it holds for letters below
/// 16, and BLASTNA letters and its sentinel (15) are.
#[inline(always)]
fn four_base_byte_of_word(w: u32) -> u8 {
    debug_assert!(w & 0xF0F0 == 0, "BLASTNA letters below 16");
    let v = (w << 6) | (w >> 4);
    (v | (v >> 20)) as u8
}

#[inline(always)]
fn word_at(q: &[u8], at: usize) -> u32 {
    debug_assert!(at + 4 <= q.len());
    // SAFETY: the callers' loops keep `at + 4` within the slice (see each loop).
    u32::from_le_bytes(unsafe { *(q.as_ptr().add(at) as *const [u8; 4]) })
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:261-349
/// `extend_hit_ungapped_approx_ncbi` with the arguments it takes (bit 0 of the probe).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn approx_blocks(
    q_seq: &[u8],
    s_seq_packed: &[u8],
    q_off: usize,
    s_off: usize,
    s_match_end: usize,
    subject_len: usize,
    x_dropoff: i32,
    score_table: &[i32; 256],
    reduced_cutoff: i32,
    matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
) -> UngappedData {
    let x = x_dropoff.wrapping_neg();
    let q_len = q_seq.len();
    let len = (COMPRESSION_RATIO - (s_off % COMPRESSION_RATIO)) % COMPRESSION_RATIO;
    let q_ext = q_off + len;
    let s_ext = s_off + len;
    let s_ext_byte = s_ext / COMPRESSION_RATIO;

    let mut score = 0i32;
    let mut sum = 0i32;

    // Left: today's loop takes block k (from 0) at query `q_ext - 4 - 4k` and subject byte
    // `s_ext / 4 - 1 - k` for k below `left_len`; its tests `q_idx < 4` and `s_idx == 0` never fire
    // there (`4 * left_len <= q_ext`, `left_len <= s_ext / 4`). Here block `m = left_len - 1 - k`
    // of the slices that start at `q_lo` and `s_ext / 4 - left_len`; `best` is the block of `new_q`
    // (`left_len` = `q_ext`, no positive block).
    let left_len = q_ext.min(s_ext) / COMPRESSION_RATIO;
    let q_lo = q_ext - COMPRESSION_RATIO * left_len;
    let mut best = left_len;
    if left_len > 0 {
        let qw = &q_seq[q_lo..q_ext];
        let sw = &s_seq_packed[s_ext_byte - left_len..s_ext_byte];
        let mut m = left_len;
        while m > 0 {
            m -= 1;
            // SAFETY: m < left_len = sw.len(), and 4 * m + 4 <= qw.len().
            let sb = unsafe { *sw.get_unchecked(m) };
            let qb = four_base_byte_of_word(word_at(qw, COMPRESSION_RATIO * m));
            sum = sum.wrapping_add(score_table[(qb ^ sb) as usize]);
            if sum > 0 {
                best = m;
                score = score.wrapping_add(sum);
                sum = 0;
            }
            if sum < x {
                break;
            }
        }
    }
    let q_start = q_lo + COMPRESSION_RATIO * best;

    // Right: block k at query `q_ext + 4k` and subject byte `s_ext / 4 + k`; `end` counts the
    // blocks up to the last positive one (`new_q` = `q_ext + 4 * end - 1`, or `q_ext` for none).
    let right_len = q_len
        .saturating_sub(q_ext)
        .min(subject_len.saturating_sub(s_ext))
        / COMPRESSION_RATIO;
    sum = 0;
    let mut end = 0usize;
    if right_len > 0 {
        let qw = &q_seq[q_ext..q_ext + COMPRESSION_RATIO * right_len];
        let sw = &s_seq_packed[s_ext_byte..s_ext_byte + right_len];
        let mut m = 0usize;
        while m < right_len {
            // SAFETY: m < right_len = sw.len(), and 4 * m + 4 <= qw.len().
            let sb = unsafe { *sw.get_unchecked(m) };
            let qb = four_base_byte_of_word(word_at(qw, COMPRESSION_RATIO * m));
            m += 1;
            sum = sum.wrapping_add(score_table[(qb ^ sb) as usize]);
            if sum > 0 {
                end = m;
                score = score.wrapping_add(sum);
                sum = 0;
            }
            if sum < x {
                break;
            }
        }
    }

    if score >= reduced_cutoff {
        return super::extension::extend_hit_ungapped_exact_ncbi(
            q_seq,
            s_seq_packed,
            q_off,
            s_off,
            subject_len,
            x_dropoff,
            matrix,
        );
    }

    let new_q = if end == 0 {
        q_ext
    } else {
        q_ext + COMPRESSION_RATIO * end - 1
    };
    let s_start = s_ext - (q_ext - q_start);
    let s_match_len = s_match_end.saturating_sub(s_start);
    let right_len = if new_q >= q_start {
        new_q - q_start + 1
    } else {
        0
    };
    UngappedData {
        q_start,
        s_start,
        length: s_match_len.max(right_len),
        score,
    }
}

#[inline(always)]
fn cell(matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE], code: u8, base: usize) -> i32 {
    debug_assert!((code as usize) < BLASTNA_SIZE && base < 4);
    // SAFETY: `code` is a BLASTNA letter (`< 16`) and `base` an ncbi2na base (`< 4`), so the
    // index is below 256.
    unsafe { *matrix.get_unchecked(((code as usize) << 4) | base) }
}

/// The subject bases from position `p` down: base `p` in bits 0..2, base `p - 1` in bits 2..4, ...;
/// and how many positions it holds (`p` down to `4 * (p / 4 - 7)`, or down to 0).
#[inline(always)]
fn window_down(s: &[u8], p: usize) -> (u64, usize) {
    let b = p >> 2;
    let shift = ((!p & 3) << 1) as u32;
    if b >= 7 {
        let raw = u64::from_be_bytes(s[b - 7..=b].try_into().expect("eight bytes"));
        (raw >> shift, 29 + (p & 3))
    } else {
        let raw = s[..=b].iter().fold(0u64, |acc, &byte| (acc << 8) | byte as u64);
        (raw >> shift, p + 1)
    }
}

/// The subject bases from position `p` up: base `p` in bits 62..64, base `p + 1` in bits 60..62,
/// ...; and how many positions it holds (to the end of byte `p / 4 + 7`, or of the slice).
#[inline(always)]
fn window_up(s: &[u8], p: usize) -> (u64, usize) {
    let b = p >> 2;
    let shift = ((p & 3) << 1) as u32;
    if b + 8 <= s.len() {
        let raw = u64::from_be_bytes(s[b..b + 8].try_into().expect("eight bytes"));
        (raw << shift, 32 - (p & 3))
    } else {
        let tail = &s[b..];
        let raw = tail
            .iter()
            .fold(0u64, |acc, &byte| (acc << 8) | byte as u64)
            << (8 * (8 - tail.len()));
        (raw << shift, COMPRESSION_RATIO * tail.len() - (p & 3))
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:148-153
/// `extend_hit_ungapped_exact_counted` with the arguments it takes (bit 2 of the probe): the same
/// steps (left: step `i` from 1 reads subject `s_off - i` and query `q_off - i`; right: step `i`
/// from 0 reads `s_off + i` and `q_off + i`), the bases taken from a window.
#[inline(never)]
pub(crate) fn exact_window(
    q_seq: &[u8],
    s_seq_packed: &[u8],
    q_off: usize,
    s_off: usize,
    subject_len: usize,
    x_dropoff: i32,
    matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
) -> UngappedData {
    let x = x_dropoff.wrapping_neg();
    debug_assert!(q_off < q_seq.len() && s_off < subject_len);
    debug_assert!(subject_len.div_ceil(4) <= s_seq_packed.len());
    let q = q_seq.as_ptr();
    let mut score: i32 = 0;
    let mut sum: i32 = 0;

    let left = q_off.min(s_off);
    let mut q_beg = q_off;
    let mut i = 0;
    'left: while i < left {
        let (mut w, n) = window_down(s_seq_packed, s_off - 1 - i);
        let stop = left.min(i + n);
        while i < stop {
            i += 1;
            let base = (w & 3) as usize;
            w >>= 2;
            // SAFETY: q_off - i < q_off < q_seq.len() (i <= left <= q_off).
            let code = unsafe { *q.add(q_off - i) };
            sum = sum.wrapping_add(cell(matrix, code, base));
            if sum > 0 {
                q_beg = q_off - i;
                score = score.wrapping_add(sum);
                sum = 0;
            } else if sum < x {
                break 'left;
            }
        }
    }
    let q_start = q_beg;
    let s_start = s_off - (q_off - q_start);

    let right = (q_seq.len() - q_off).min(subject_len - s_off);
    let mut q_end = q_off;
    sum = 0;
    let mut x_current = x;
    let mut i = 0;
    'right: while i < right {
        let (mut w, n) = window_up(s_seq_packed, s_off + i);
        let stop = right.min(i + n);
        while i < stop {
            let base = (w >> 62) as usize;
            w <<= 2;
            // SAFETY: q_off + i < q_off + right <= q_seq.len().
            let code = unsafe { *q.add(q_off + i) };
            i += 1;
            sum = sum.wrapping_add(cell(matrix, code, base));
            if sum > 0 {
                q_end = q_off + i;
                score = score.wrapping_add(sum);
                x_current = if score.wrapping_neg() > x {
                    score.wrapping_neg()
                } else {
                    x
                };
                sum = 0;
            } else if sum < x_current {
                break 'right;
            }
        }
    }

    UngappedData {
        q_start,
        s_start,
        length: q_end - q_beg,
        score,
    }
}

#[cfg(test)]
#[path = "x_ungap_block_tests.rs"]
mod tests;
