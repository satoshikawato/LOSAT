//! LOSAT_X_BNUNGAPCOUNT (graduated 2026-10-27, always on): the exact ungapped
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:743-747 (the same at
    // 913-917 and 1755-1759)
    // ```c
    //                s_NuclUngappedExtendExact(query, subject, matrix, q_off,
    //                                   s_off, -(cutoffs->x_dropoff), ungapped_data);
    //             }else {
    //                s_NuclUngappedExtend(query, subject, matrix, q_off, s_end, s_off,
    //                                  -(cutoffs->x_dropoff), ungapped_data,
    // ```
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

/// EXPERIMENT (LOSAT_X_FNEXK, line FN s1 probe): `extend_hit_ungapped_exact_counted` with the
/// subject read four bases at a time: each group of four positions comes from one two-byte window
/// of the packed subject, shifted once by the call's phase (`2 * (3 - first % 4)` on the left,
/// `2 * (s_off % 4)` on the right), and each step takes its base with a constant shift. The query
/// is walked by pointer. The steps, their order, the tests and the `i32` updates are those of the
/// counted loops, so every exit and the four results are theirs.
#[inline(never)]
pub(crate) fn extend_hit_ungapped_exact_unrolled(
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
    let m = matrix.as_ptr();
    let s = s_seq_packed.as_ptr();
    let s_bytes = s_seq_packed.len();
    // SAFETY (every `cell` read below): `code` is a BLASTNA letter (< 16) and the base < 4, so
    // `(code << 4) | base` < 256.
    let cell = |code: u8, base: u32| -> i32 { unsafe { *m.add(((code as usize) << 4) | base as usize) } };
    let mut score: i32 = 0;
    let mut sum: i32 = 0;

    // Left: min(q_off, s_off) steps; step j (from 0) reads subject position s_off - 1 - j and query
    // offset q_off - 1 - j. Group g holds positions first - 4 g .. first - 4 g - 3 (first =
    // s_off - 1): in the big-endian window (byte b - g - 1) << 8 | (byte b - g), b = first / 4, the
    // position first - 4 g - t sits at bits 2 (3 - k) + 2 t (k = first % 4), so the window shifted
    // right by 2 (3 - k) holds the group's t-th base at bits 2 t.
    let left = q_off.min(s_off);
    let q0 = q_seq.as_ptr();
    // SAFETY: q_off < q_seq.len(), and the walk stays in q_off - left .. q_off.
    let mut qp = unsafe { q0.add(q_off) };
    let q_stop = unsafe { q0.add(q_off - left) };
    let mut q_beg_p = qp;
    'left: {
        if left == 0 {
            break 'left;
        }
        let first = s_off - 1;
        let sh = (2 * (3 - (first & 3))) as u32;
        let mut b = first >> 2;
        loop {
            if qp == q_stop {
                break 'left;
            }
            // SAFETY: b < s_bytes (b <= first / 4 < subject_len.div_ceil(4)); byte b - 1 is read
            // only when b >= 1 (a group that reaches below position 0 stops by the step count).
            let lo = unsafe { *s.add(b) } as u32;
            let hi = if b > 0 { (unsafe { *s.add(b - 1) }) as u32 } else { 0 };
            let v = ((hi << 8) | lo) >> sh;
            macro_rules! step {
                ($t:expr) => {{
                    if qp == q_stop {
                        break 'left;
                    }
                    // SAFETY: q_stop < qp <= q0 + q_off.
                    qp = unsafe { qp.sub(1) };
                    let code = unsafe { *qp };
                    sum = sum.wrapping_add(cell(code, (v >> (2 * $t)) & 3));
                    if sum > 0 {
                        q_beg_p = qp;
                        score = score.wrapping_add(sum);
                        sum = 0;
                    } else if sum < x {
                        break 'left;
                    }
                }};
            }
            step!(0);
            step!(1);
            step!(2);
            step!(3);
            if b == 0 {
                break 'left;
            }
            b -= 1;
        }
    }
    // SAFETY: q_beg_p is in q0 ..= q0 + q_off.
    let q_start = unsafe { q_beg_p.offset_from(q0) } as usize;
    let s_start = s_off - (q_off - q_start);

    // Right: min(q_avail, s_avail) steps; step j reads subject position s_off + j and query offset
    // q_off + j. Group g holds positions s_off + 4 g .. + 3: in the window (byte b + g) << 8 |
    // (byte b + g + 1), b = s_off / 4, shifted left by 2 k (k = s_off % 4) and taken as its top
    // byte, the group's t-th base sits at bits 6 - 2 t.
    let right = (q_seq.len() - q_off).min(subject_len - s_off);
    // SAFETY: the walk reads q_off .. q_off + right <= q_seq.len().
    let mut qp = unsafe { q0.add(q_off) };
    let q_stop = unsafe { q0.add(q_off + right) };
    let mut q_end_p = qp;
    sum = 0;
    let mut x_current = x;
    'right: {
        if right == 0 {
            break 'right;
        }
        let sh = (2 * (s_off & 3)) as u32;
        let mut b = s_off >> 2;
        loop {
            if qp == q_stop {
                break 'right;
            }
            // SAFETY: a step is left, so its position s_off + j < subject_len lies in byte b
            // < s_bytes; byte b + 1 is read only when it exists.
            let hi = unsafe { *s.add(b) } as u32;
            let lo = if b + 1 < s_bytes { (unsafe { *s.add(b + 1) }) as u32 } else { 0 };
            let v = (((hi << 8) | lo) << sh) >> 8;
            macro_rules! step {
                ($t:expr) => {{
                    if qp == q_stop {
                        break 'right;
                    }
                    // SAFETY: q0 + q_off <= qp < q_stop.
                    let code = unsafe { *qp };
                    qp = unsafe { qp.add(1) };
                    sum = sum.wrapping_add(cell(code, (v >> (6 - 2 * $t)) & 3));
                    if sum > 0 {
                        q_end_p = qp;
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
                }};
            }
            step!(0);
            step!(1);
            step!(2);
            step!(3);
            b += 1;
        }
    }
    // SAFETY: q_end_p is in q0 + q_off ..= q_stop.
    let q_end = unsafe { q_end_p.offset_from(q0) } as usize;

    UngappedData {
        q_start,
        s_start,
        length: q_end - q_start,
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:197 and
    // blast_stat.c:1124-1127
    // ```c
    //         if ((sum += matrix[*--q][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
    // ...
    //         matrix[BLASTNA_SIZE-1][index1] = INT4_MIN / 2;
    // ...
    //         matrix[index1][BLASTNA_SIZE-1] = INT4_MIN / 2;
    // ```
    /// `-xdrop_ungap 1e10` (X-drop `INT_MIN`, which negates to itself): the left walk from the hit
    /// at 3 passes the two sentinels (`INT4_MIN / 2` each: `INT_MIN`) and the mismatch A/C (-3),
    /// where the `Int4` sum wraps to 2^31 - 3, above 0, as in NCBI's x86-64 build; the right walk
    /// adds the two matches: `i32::MAX`. A build with overflow checks gives the same.
    #[test]
    fn a_walk_past_two_sentinels_wraps_as_ncbis_int4() {
        let matrix = crate::algorithm::blastn::alignment::build_blastna_matrix(1, -3);
        let query = [0u8, 15, 15, 0, 0];
        // C A A A A, packed (ncbi2na: A 0, C 1).
        let packed = [0b0100_0000u8, 0];
        let expected = (0, 0, 5, i32::MAX);
        let counted =
            extend_hit_ungapped_exact_counted(&query, &packed, 3, 3, 5, i32::MIN, &matrix);
        assert_eq!(key(counted), expected);
        let ported = extend_hit_ungapped_exact_ported(&query, &packed, 3, 3, 5, i32::MIN, &matrix);
        assert_eq!(key(ported), expected);
        let unrolled =
            extend_hit_ungapped_exact_unrolled(&query, &packed, 3, 3, 5, i32::MIN, &matrix);
        assert_eq!(key(unrolled), expected);
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
                    query.push(if r.is_multiple_of(50) {
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
                    for x_dropoff in [5, 16, 40, i32::MIN] {
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
                        assert_eq!(
                            key(extend_hit_ungapped_exact_unrolled(
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
                            "unrolled round {round} {args:?}"
                        );
                        calls += 1;
                    }
                }
            }
        }
        assert!(calls > 100_000, "{calls}");
    }
}
