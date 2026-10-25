//! EXPERIMENT (LOSAT_X_BNUNGAPCOUNT): the exact ungapped
//! extension (`s_NuclUngappedExtendExact`) with the number of steps of each direction worked out
//! before its loop.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:164-241
//! ```c
//!     q_avail = query->length - q_off;
//!     s_avail = subject->length - s_off;
//! ...
//!     if (q_off < s_off) {
//!         start = subject0 + (s_off - q_off) / COMPRESSION_RATIO;
//!         remainder = 3 - ((s_off - q_off) % COMPRESSION_RATIO);
//!     } else {
//!         start = subject0;
//!         remainder = 3;
//!     }
//! ...
//!     while ((s > start) || (s == start && base < remainder)) {
//!         if (base == 3) {
//!             s--;
//!             base = 0;
//!         } else {
//!             ++base;
//!         }
//!         ch = *s;
//!         if ((sum += matrix[*--q][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
//! ...
//!     if (q_avail < s_avail) {
//!         sf = subject0 + (s_off + q_avail) / COMPRESSION_RATIO;
//!         remainder = 3 - ((s_off + q_avail) % COMPRESSION_RATIO);
//!     } else {
//!         sf = subject0 + (subject->length) / COMPRESSION_RATIO;
//!         remainder = 3 - ((subject->length) % COMPRESSION_RATIO);
//!     }
//! ...
//!     while (s < sf || (s == sf && base > remainder)) {
//!         ch = *s;
//!         if ((sum += matrix[*q++][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
//! ...
//!         if (base == 0) {
//!             base = 3;
//!             s++;
//!         } else
//!             base--;
//!     }
//!
//!     ungapped_data->length = (Int4)(q_end - q_beg);
//! ```
//!
//! Where it runs. `extension::extend_hit_ungapped_exact_ncbi` (every caller of the exact
//! extension: blastn-short and other word sizes below 11, and the exact pass that follows
//! `extend_hit_ungapped_approx_ncbi` when its score reaches the reduced cutoff) calls
//! `extend_hit_ungapped_exact_counted` when the switch is on, and the ported loops otherwise.
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
//!
//! The variable is read once. `LOSAT_X_BNUNGAPCOUNT` (a switch of the strict set,
//! `utils/xswitch.rs`: on unless it is `0`) runs this module.

use super::extension::UngappedData;
use std::sync::OnceLock;

/// BLASTNA alphabet size (`blast_encoding.h:91-92`).
const BLASTNA_SIZE: usize = 16;

/// No NCBI counterpart: reads LOSAT_X_BNUNGAPCOUNT once. It is a strict-set switch: on unless it
/// is `0`. It chooses a LOSAT code path and changes no value NCBI computes.
#[inline(always)]
pub(crate) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNUNGAPCOUNT"))
}

/// `matrix[code][base]` of the BLASTNA matrix (`code < 16`, `base < 4`).
#[inline(always)]
fn cell(matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE], code: u8, base: usize) -> i32 {
    debug_assert!((code as usize) < BLASTNA_SIZE && base < 4);
    // SAFETY: `code` is a BLASTNA letter (`< 16`) and `base` an ncbi2na base (`< 4`), so the
    // index is below 256.
    unsafe { *matrix.get_unchecked(((code as usize) << 4) | base) }
}

/// NCBI reference: c++/src/algo/blast/core/na_ungapped.c:149-153
/// ```c
/// static void
/// s_NuclUngappedExtendExact(BLAST_SequenceBlk * query,
///                           BLAST_SequenceBlk * subject, Int4 ** matrix,
///                           Int4 q_off, Int4 s_off, Int4 X,
///                           BlastUngappedData * ungapped_data)
/// ```
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
    let x = -x_dropoff;
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
        sum += cell(matrix, code, ((ch >> ((!p & 3) << 1)) & 3) as usize);
        if sum > 0 {
            q_beg = q_off - i;
            score += sum;
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
        sum += cell(matrix, code, ((ch >> ((!p & 3) << 1)) & 3) as usize);
        if sum > 0 {
            q_end = q_off + i;
            score += sum;
            x_current = if -score > x { -score } else { x };
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
mod tests {
    use super::super::extension::extend_hit_ungapped_exact_ported;
    use super::*;

    fn next(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn key(u: UngappedData) -> (usize, usize, usize, i32) {
        (u.q_start, u.s_start, u.length, u.score)
    }

    /// Every hit offset of random query blocks (contexts of 1 to 40 letters between sentinels,
    /// a few ambiguity codes) against random and periodic subjects of 1 to 140 bases packed into
    /// exactly `len.div_ceil(4)` bytes, at three X-drops: the counted loops give the ported
    /// loops' four results.
    #[test]
    fn every_offset_matches_the_ported_loops() {
        let matrix = crate::algorithm::blastn::alignment::build_blastna_matrix(1, -3);
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut calls = 0u64;
        for round in 0..400 {
            let mut query = Vec::new();
            for _ in 0..1 + next(&mut state) % 4 {
                for _ in 0..1 + next(&mut state) % 40 {
                    let r = next(&mut state);
                    query.push(if r % 50 == 0 {
                        (4 + r % 11) as u8
                    } else {
                        (r % 4) as u8
                    });
                }
                query.push(15);
            }
            query.pop();
            let subject_len = 1 + (next(&mut state) % 140) as usize;
            let period = 1 + (next(&mut state) % 6) as usize;
            let motif: Vec<u8> = (0..period).map(|_| (next(&mut state) % 4) as u8).collect();
            let bases: Vec<u8> = (0..subject_len)
                .map(|i| {
                    if round % 2 == 0 {
                        motif[i % period]
                    } else {
                        (next(&mut state) % 4) as u8
                    }
                })
                .collect();
            let mut packed = vec![0u8; subject_len.div_ceil(4)];
            for (i, b) in bases.iter().enumerate() {
                packed[i / 4] |= b << (2 * (3 - i % 4));
            }
            for q_off in 0..query.len() {
                for s_off in 0..subject_len {
                    for x_dropoff in [5, 16, 40] {
                        let args = (q_off, s_off, subject_len, x_dropoff);
                        assert_eq!(
                            key(extend_hit_ungapped_exact_counted(
                                &query,
                                &packed,
                                q_off,
                                s_off,
                                subject_len,
                                x_dropoff,
                                &matrix
                            )),
                            key(extend_hit_ungapped_exact_ported(
                                &query,
                                &packed,
                                q_off,
                                s_off,
                                subject_len,
                                x_dropoff,
                                &matrix
                            )),
                            "round {round} {args:?}"
                        );
                        calls += 1;
                    }
                }
            }
        }
        assert!(calls > 100_000, "{calls}");
    }
}
