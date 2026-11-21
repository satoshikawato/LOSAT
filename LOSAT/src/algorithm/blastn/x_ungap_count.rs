//! LOSAT_X_BNUNGAPCOUNT (graduated 2026-10-27, always on): the exact ungapped
//! extension (`s_NuclUngappedExtendExact`) with the number of steps of each direction worked out
//! before its loop.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:164-241
//!
//! Where it runs. `extension::extend_hit_ungapped_exact_ncbi` (every caller of the exact
//! extension: blastn-short and other word sizes below 11, and the exact pass that follows
//! `extend_hit_ungapped_approx_ncbi` when its score reaches the reduced cutoff) calls
//! `extend_hit_ungapped_exact_counted`; the loop-for-loop port (`extend_hit_ungapped_exact_ported`,
//! `cfg(test)`) is the reference of the unit test below.
//!
//! Why the result is the same. The pair `(s, base)` is subject position `p = 4 * s + 3 - base`
//! (`base = 3 - p % 4`, and `NCBI2NA_UNPACK_BASE(ch, base)` is the base at `p`):
//!
//! - Left. `start` / `remainder` stand for position `max(s_off - q_off, 0)`, and
//!   `(s > start) || (s == start && base < remainder)` is "the position is above it". Each step
//!   moves one position down before it reads, from `s_off`, so the loop takes at most
//!   `s_off - max(s_off - q_off, 0) = min(q_off, s_off)` steps; step `i` (from 1) reads subject
//!   position `s_off - i` and query offset `q_off - i`.
//! - Right. `sf` / `remainder` stand for position `s_off + min(q_avail, s_avail)` (`q_avail <
//!   s_avail` picks `s_off + q_avail`, else `subject->length = s_off + s_avail`), and
//!   `s < sf || (s == sf && base > remainder)` is "the position is below it". Each step reads and
//!   then moves one position up, from `s_off`, so the loop takes at most `min(q_avail, s_avail)`
//!   steps; step `i` (from 0) reads subject position `s_off + i` and query offset `q_off + i`.
//!
//! This module runs those counts as plain counted loops and reads the base at `p` as
//! `(byte[p / 4] >> (2 * (3 - p % 4))) & 3`, with the C tests (`sum > 0`, then `sum < X` /
//! `sum < X_current`) and updates in the C order on the same `i32` values, so every step, every
//! exit and the four results are the C ones. The query offsets stay in `0..q_len` within the
//! counts, so the ported loops' `q_idx == 0` / `q_idx >= q_len` exits (which never fire) are not
//! needed, and the query letters are BLASTNA codes (`< 16`), so the matrix cell
//! `code * 16 + base` is in the 256-cell matrix.

use super::extension::UngappedData;

/// BLASTNA alphabet size (`blast_encoding.h:91-92`).
const BLASTNA_SIZE: usize = 16;

/// `matrix[code][base]` of the BLASTNA matrix (`code < 16`, `base < 4`).
#[inline(always)]
fn cell(matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE], code: u8, base: usize) -> i32 {
    debug_assert!((code as usize) < BLASTNA_SIZE && base < 4);
    // SAFETY: `code` is a BLASTNA letter (`< 16`) and `base` an ncbi2na base (`< 4`), so the
    // index is below 256.
    unsafe { *matrix.get_unchecked(((code as usize) << 4) | base) }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:148-153
/// The arguments are those of `extend_hit_ungapped_exact_ncbi`: the query buffer without its outer
/// sentinels (BLASTNA, one byte a letter), the packed subject (ncbi2na, four bases a byte, at least
/// `subject_len.div_ceil(4)` bytes), the offsets of the hit, the subject length, the X-drop as a
/// positive number and the BLASTNA matrix.
#[inline(never)]
pub(crate) fn extend_hit_ungapped_exact_counted(
    q_seq: &[u8],
    s_seq_packed: &[u8],
    q_off: usize,
    s_off: usize,
    subject_len: usize,
    x_dropoff: i32,
    matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
) -> UngappedData {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:743-747 (the same at
    // 913-917 and 1755-1759)
    // A raw X-drop past `Int4` (`INT_MIN`, `core::x_drop`) negates to itself in NCBI's x86-64
    // build, so the extension never stops at the X-drop; LOSAT negates it the same way.
    // Such a walk can pass two of the query's sentinels (BLASTNA 15, `INT4_MIN / 2` against
    // every base, blast_stat.c:1124-1127): NCBI's `Int4` sums then wrap in its x86-64 build, and
    // LOSAT's sums wrap the same way (the line's decision 22: x86-64's two's complement), also in
    // a build with overflow checks.
    let x = x_dropoff.wrapping_neg();
    debug_assert!(q_off < q_seq.len() && s_off < subject_len);
    debug_assert!(subject_len.div_ceil(4) <= s_seq_packed.len());
    let q = q_seq.as_ptr();
    let s = s_seq_packed.as_ptr();
    let mut score: i32 = 0;
    let mut sum: i32 = 0;

    // Extend to the left: min(q_off, s_off) steps at most.
    let left = q_off.min(s_off);
    let mut q_beg = q_off;
    let mut i = 0;
    while i < left {
        i += 1;
        let p = s_off - i;
        // SAFETY: p < s_off < subject_len, so byte p / 4 is in the packed subject, and
        // q_off - i < q_off < q_seq.len() (the callers' offsets of a hit; i <= left).
        let ch = unsafe { *s.add(p >> 2) };
        let code = unsafe { *q.add(q_off - i) };
        sum = sum.wrapping_add(cell(matrix, code, ((ch >> ((!p & 3) << 1)) & 3) as usize));
        if sum > 0 {
            q_beg = q_off - i;
            score = score.wrapping_add(sum);
            sum = 0;
        } else if sum < x {
            break;
        }
    }
    let q_start = q_beg;
    let s_start = s_off - (q_off - q_start);

    // Extend to the right: min(q_avail, s_avail) steps at most.
    let right = (q_seq.len() - q_off).min(subject_len - s_off);
    let mut q_end = q_off;
    sum = 0;
    let mut x_current = x;
    let mut i = 0;
    while i < right {
        let p = s_off + i;
        // SAFETY: p < s_off + right <= subject_len and q_off + i < q_off + right <= q_seq.len().
        let ch = unsafe { *s.add(p >> 2) };
        let code = unsafe { *q.add(q_off + i) };
        i += 1;
        sum = sum.wrapping_add(cell(matrix, code, ((ch >> ((!p & 3) << 1)) & 3) as usize));
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
            break;
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
#[path = "x_ungap_count_tests.rs"]
mod tests;
