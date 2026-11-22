use crate::core::blast_encoding::COMPRESSION_RATIO;

/// BLASTNA alphabet size.
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_encoding.h:91-92
const BLASTNA_SIZE: usize = 16;

/// Ungapped extension results (0-based coordinates, length in bases).
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:141-147
#[derive(Clone, Copy, Debug)]
pub struct UngappedData {
    pub q_start: usize,
    pub s_start: usize,
    pub length: usize,
    pub score: i32,
}

/// Build the nucleotide ungapped score table for 4-base blocks.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:236-259
pub fn build_nucl_score_table(reward: i32, penalty: i32) -> [i32; 256] {
    let mut table = [0i32; 256];
    for i in 0..256 {
        let mut score = 0;
        if (i & 3) != 0 {
            score += penalty;
        } else {
            score += reward;
        }
        if ((i >> 2) & 3) != 0 {
            score += penalty;
        } else {
            score += reward;
        }
        if ((i >> 4) & 3) != 0 {
            score += penalty;
        } else {
            score += reward;
        }
        if (i >> 6) != 0 {
            score += penalty;
        } else {
            score += reward;
        }
        table[i] = score;
    }
    table
}

/// Extract base N from packed ncbi2na byte.
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_util.h:51-55
#[cfg(test)]
#[inline]
fn ncbi2na_unpack_base(byte: u8, n: u8) -> u8 {
    (byte >> (2 * n)) & 0x03
}

/// Exact ungapped extension (per-base) for nucleotide sequences.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:148-242
///
/// LOSAT_X_BNUNGAPCOUNT (graduated 2026-10-27, always on): the same loops run with their step
/// counts worked out first (`x_ungap_count.rs`). The loop-for-loop port below stays as the
/// reference of that module's unit test.
#[inline(always)]
pub(crate) fn extend_hit_ungapped_exact_ncbi(
    q_seq: &[u8],        // BLASTNA (1 byte/base)
    s_seq_packed: &[u8], // ncbi2na packed (4 bases/byte)
    q_off: usize,
    s_off: usize,
    subject_len: usize,
    x_dropoff: i32,
    matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
) -> UngappedData {
    super::x_ungap_count::extend_hit_ungapped_exact_counted(
        q_seq,
        s_seq_packed,
        q_off,
        s_off,
        subject_len,
        x_dropoff,
        matrix,
    )
}

/// The port of `s_NuclUngappedExtendExact` above, loop for loop: the reference of the unit test of
/// `x_ungap_count.rs` (it was the path of `LOSAT_X_BNUNGAPCOUNT=0`).
#[cfg(test)]
#[inline(never)]
pub(crate) fn extend_hit_ungapped_exact_ported(
    q_seq: &[u8],        // BLASTNA (1 byte/base)
    s_seq_packed: &[u8], // ncbi2na packed (4 bases/byte)
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
    let x = x_dropoff.wrapping_neg();
    let q_len = q_seq.len();

    let mut base: isize = 3 - (s_off % COMPRESSION_RATIO) as isize;
    let mut s_idx: isize = (s_off / COMPRESSION_RATIO) as isize;
    let (start_idx, remainder) = if q_off < s_off {
        let diff = s_off - q_off;
        (
            (diff / COMPRESSION_RATIO) as isize,
            3 - (diff % COMPRESSION_RATIO) as isize,
        )
    } else {
        (0isize, 3isize)
    };

    let mut q_idx: isize = q_off as isize;
    let mut q_beg: isize = q_idx;
    let mut q_end: isize = q_idx;
    let mut score: i32 = 0;
    let mut sum: i32 = 0;

    // Left extension (per-base).
    while (s_idx > start_idx) || (s_idx == start_idx && base < remainder) {
        if base == 3 {
            s_idx -= 1;
            base = 0;
        } else {
            base += 1;
        }

        if q_idx == 0 {
            break;
        }
        q_idx -= 1;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:188-200
        // SAFETY: q_idx and s_idx are kept within bounds by the loop guards
        // and mirror NCBI's pointer arithmetic on query/subject buffers.
        let ch = unsafe { *s_seq_packed.get_unchecked(s_idx as usize) };
        let s_base = ncbi2na_unpack_base(ch, base as u8) as usize;
        let q_code = unsafe { *q_seq.get_unchecked(q_idx as usize) } as usize;
        sum = sum.wrapping_add(matrix[q_code * BLASTNA_SIZE + s_base]);
        if sum > 0 {
            q_beg = q_idx;
            score = score.wrapping_add(sum);
            sum = 0;
        } else if sum < x {
            break;
        }
    }

    let q_start = q_beg as usize;
    let s_start = s_off - (q_off - q_start);

    let q_avail = q_len - q_off;
    let s_avail = subject_len - s_off;
    let (sf_idx, remainder) = if q_avail < s_avail {
        let pos = s_off + q_avail;
        (
            (pos / COMPRESSION_RATIO) as isize,
            3 - (pos % COMPRESSION_RATIO) as isize,
        )
    } else {
        (
            (subject_len / COMPRESSION_RATIO) as isize,
            3 - (subject_len % COMPRESSION_RATIO) as isize,
        )
    };

    // Right extension (per-base).
    let mut q_idx: isize = q_off as isize;
    let mut s_idx: isize = (s_off / COMPRESSION_RATIO) as isize;
    let mut base: isize = 3 - (s_off % COMPRESSION_RATIO) as isize;
    sum = 0;
    let mut x_current = x;
    while (s_idx < sf_idx) || (s_idx == sf_idx && base > remainder) {
        if q_idx as usize >= q_len {
            break;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:217-229
        // SAFETY: q_idx and s_idx stay within bounds under the same
        // conditions as the NCBI pointer-based loop.
        let ch = unsafe { *s_seq_packed.get_unchecked(s_idx as usize) };
        let s_base = ncbi2na_unpack_base(ch, base as u8) as usize;
        let q_code = unsafe { *q_seq.get_unchecked(q_idx as usize) } as usize;
        q_idx += 1;

        sum = sum.wrapping_add(matrix[q_code * BLASTNA_SIZE + s_base]);
        if sum > 0 {
            q_end = q_idx;
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

        if base == 0 {
            base = 3;
            s_idx += 1;
        } else {
            base -= 1;
        }
    }

    let length = (q_end - q_beg).max(0) as usize;
    UngappedData {
        q_start,
        s_start,
        length,
        score,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:292-294,322-324
// NCBI's overlapping BLASTNA bit expression of the query letters `start..start + 4`,
// ambiguities included, computed where it is read as NCBI does (LOSAT's memory: no array
// of these bytes at every offset of the query).
#[inline(always)]
pub fn query_four_base_byte(query: &[u8], start: usize) -> u8 {
    let q: &[u8; 4] = query[start..start + 4]
        .try_into()
        .expect("four query letters");
    (q[0] << 6) | (q[1] << 4) | (q[2] << 2) | q[3]
}

// The byte of `query_four_base_byte` at every starting offset (the tests' reference).
#[cfg(test)]
pub fn build_query_four_base_bytes(query: &[u8]) -> Vec<u8> {
    (0..query.len().saturating_sub(3))
        .map(|start| query_four_base_byte(query, start))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:294,324
// C promotes each BLASTNA byte before truncating the expression to Uint1.
#[cfg(test)]
#[path = "extension_four_base_tests.rs"]
mod four_base_tests;

/// NCBI's compressed query of the small-query lookup table: byte `j` holds the query
/// letters `j` to `j + 3` as two-bit bases (`blastna & 3`, so an ambiguity code or a
/// sentinel reads as a base), with 0 before the first and after the last letter. The
/// returned vector holds NCBI's `new_seq[j]` at `j + 3` (`j` from -3).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:473-499
/// (`new_seq[i - max_start]` is `Uint1`; with `max_start` = 3 each byte is the last
/// four letters read, and the prologue writes -3..-1 for queries of three or more letters.)
pub fn build_compressed_query(query: &[u8]) -> Vec<u8> {
    let base = |index: isize| -> u8 {
        if index >= 0 && (index as usize) < query.len() {
            query[index as usize] & 3
        } else {
            0
        }
    };
    (-3..query.len() as isize)
        .map(|j| base(j) << 6 | base(j + 1) << 4 | base(j + 2) << 2 | base(j + 3))
        .collect()
}

/// The number of matching two-bit bases at the right end (`left`) or the left end
/// (`right`) of the XOR of a query byte and a subject byte, at most 4.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1326-1330
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1349-1354
fn exact_match_extend_left(xor: u8) -> i64 {
    i64::from(xor.trailing_zeros().min(8) / 2)
}

fn exact_match_extend_right(xor: u8) -> i64 {
    i64::from(xor.leading_zeros().min(8) / 2)
}

/// The start of a word hit of a small-query lookup table, from the start of its lookup
/// word: NCBI's exact-match extension over the compressed query (`build_compressed_query`)
/// and the packed subject. `None` where NCBI drops the hit.
pub struct SmallNaWord<'a> {
    pub compressed_query: &'a [u8],
    pub subject_packed: &'a [u8],
    /// NCBI's `query->length` (the concatenated query).
    pub query_length: usize,
    pub word_length: usize,
    pub lut_word_length: usize,
}

impl SmallNaWord<'_> {
    /// `q[offset]` of NCBI's compressed query (offset from -3).
    fn query_byte(&self, offset: i64) -> u8 {
        self.compressed_query[(offset + 3) as usize]
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1514-1566
    /// Returns the query and subject offsets of the word.
    pub fn extend(
        &self,
        q_offset: usize,
        s_offset: usize,
        q_start: usize,
        q_range: usize,
        s_range: usize,
    ) -> Option<(usize, usize)> {
        let ratio = COMPRESSION_RATIO as i64;
        let (q_offset, s_offset) = (q_offset as i64, s_offset as i64);
        let (q_start, q_range, s_range) = (q_start as i64, q_range as i64, s_range as i64);
        let word_length = self.word_length as i64;
        let mut ext_max = (word_length - self.lut_word_length as i64)
            .min(s_offset)
            .min(q_offset - q_start);
        let rsdl = ratio - s_offset % ratio;
        let (s_offset, q_offset) = (s_offset + rsdl, q_offset + rsdl);
        ext_max += rsdl;

        let (mut s_off, mut q_off) = (s_offset, q_offset);
        let mut ext_left = 0;
        while ext_left < ext_max {
            let q_byte = self.query_byte(q_off - 4);
            let s_byte = self.subject_packed[(s_off / ratio - 1) as usize];
            let bases = exact_match_extend_left(q_byte ^ s_byte);
            ext_left += bases;
            if bases < 4 {
                break;
            }
            q_off -= 4;
            s_off -= 4;
        }
        ext_left = ext_left.min(ext_max);

        let (mut s_off, mut q_off) = (s_offset, q_offset);
        let ext_max = (word_length - ext_left)
            .min(s_range - s_off)
            .min(q_range - q_off);
        let mut ext_right = 0;
        while ext_right < ext_max {
            let q_byte = self.query_byte(q_off);
            let s_byte = self.subject_packed[(s_off / ratio) as usize];
            let bases = exact_match_extend_right(q_byte ^ s_byte);
            ext_right += bases;
            if bases < 4 {
                break;
            }
            q_off += 4;
            s_off += 4;
        }
        ext_right = ext_right.min(ext_max);

        if ext_left + ext_right < word_length {
            return None;
        }
        Some((
            (q_offset - ext_left) as usize,
            (s_offset - ext_left) as usize,
        ))
    }

    /// The extension for a lookup width that is a multiple of 4, a scan step that is a
    /// multiple of 4 and a word at most 4 letters longer than the lookup word.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1415-1435
    pub fn extend_aligned_one_byte(
        &self,
        q_offset: usize,
        s_offset: usize,
        q_start: usize,
        q_range: usize,
        s_range: usize,
    ) -> Option<(usize, usize)> {
        let ratio = COMPRESSION_RATIO as i64;
        let (q_offset, s_offset) = (q_offset as i64, s_offset as i64);
        let (q_start, q_range, s_range) = (q_start as i64, q_range as i64, s_range as i64);
        let lut_word_length = self.lut_word_length as i64;
        let ext_to = self.word_length as i64 - lut_word_length;
        let mut ext_left = 0;
        if s_offset > 0 && q_offset > 0 {
            let q_byte = self.query_byte(q_offset - 4);
            let s_byte = self.subject_packed[(s_offset / ratio - 1) as usize];
            ext_left = exact_match_extend_left(q_byte ^ s_byte)
                .min(ext_to)
                .min(q_offset - q_start);
        }
        if ext_left < ext_to && q_offset + lut_word_length < self.query_length as i64 {
            let q_byte = self.query_byte(q_offset + lut_word_length);
            let s_byte = self.subject_packed[((s_offset + lut_word_length) / ratio) as usize];
            let ext_right = exact_match_extend_right(q_byte ^ s_byte)
                .min(s_range - (s_offset + lut_word_length))
                .min(q_range - (q_offset + lut_word_length));
            if ext_left + ext_right < ext_to {
                return None;
            }
        }
        Some((
            (q_offset - ext_left) as usize,
            (s_offset - ext_left) as usize,
        ))
    }
}

#[cfg(test)]
#[path = "extension_small_na_tests.rs"]
mod small_na_tests;

/// Approximate ungapped extension (4-base blocks) with optional exact recomputation.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:261-349
pub(crate) fn extend_hit_ungapped_approx_ncbi(
    q_seq: &[u8],        // BLASTNA (1 byte/base)
    s_seq_packed: &[u8], // ncbi2na packed (4 bases/byte)
    q_off: usize,
    s_off: usize,
    s_match_end: usize,
    subject_len: usize,
    x_dropoff: i32,
    score_table: &[i32; 256],
    reduced_cutoff: i32,
    matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
) -> UngappedData {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:743-747 (the same at
    // 913-917 and 1755-1759)
    // A raw X-drop past `Int4` (`INT_MIN`, `core::x_drop`) negates to itself in NCBI's x86-64
    // build, so the extension never stops at the X-drop; LOSAT negates it the same way.
    let x = x_dropoff.wrapping_neg();
    let q_len = q_seq.len();

    let len = (COMPRESSION_RATIO - (s_off % COMPRESSION_RATIO)) % COMPRESSION_RATIO;
    let q_ext = q_off + len;
    let s_ext = s_off + len;

    let mut q_idx = q_ext as isize;
    let mut s_idx = (s_ext / COMPRESSION_RATIO) as isize;
    let left_len = q_ext.min(s_ext) / COMPRESSION_RATIO;

    let mut score = 0i32;
    let mut sum = 0i32;
    let mut new_q = q_idx;

    // Left extension in 4-base blocks.
    // The sums wrap as NCBI's `Int4` do in its x86-64 build (with the X-drop `INT_MIN` of a raw X-drop
    // past `Int4`; `x_ungap_count::extend_hit_ungapped_exact_counted`).
    for _ in 0..left_len {
        if q_idx < 4 || s_idx == 0 {
            break;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:292-304
        // SAFETY: q_idx >= 4 and s_idx > 0 are enforced by the loop guard.
        let s_byte = unsafe { *s_seq_packed.get_unchecked((s_idx - 1) as usize) };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:294
        let q_byte = query_four_base_byte(q_seq, (q_idx - 4) as usize);
        sum = sum.wrapping_add(score_table[(q_byte ^ s_byte) as usize]);
        if sum > 0 {
            new_q = q_idx - 4;
            score = score.wrapping_add(sum);
            sum = 0;
        }
        if sum < x {
            break;
        }
        q_idx -= 4;
        s_idx -= 1;
    }

    let q_start = new_q.max(0) as usize;
    let s_start = s_ext - (q_ext - q_start);

    // Right extension in 4-base blocks.
    let mut q_idx = q_ext as isize;
    let mut s_idx = (s_ext / COMPRESSION_RATIO) as isize;
    // NCBI's `Int4` lengths: one below 0 gives no block. A seed of a MegaBLAST database index can
    // start in the subject's last 3 letters (its seed runs on into the next sequence of the
    // index's local id), so that `s_ext`, the next multiple of 4, lies past the last letter.
    let right_len = q_len
        .saturating_sub(q_ext)
        .min(subject_len.saturating_sub(s_ext))
        / COMPRESSION_RATIO;
    sum = 0;
    new_q = q_idx;

    for _ in 0..right_len {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:322-334
        // SAFETY: q_idx + 3 and s_idx stay within bounds by the loop length
        // derived from query/subject lengths (NCBI uses the same limits).
        let s_byte = unsafe { *s_seq_packed.get_unchecked(s_idx as usize) };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:324
        let q_byte = query_four_base_byte(q_seq, q_idx as usize);
        sum = sum.wrapping_add(score_table[(q_byte ^ s_byte) as usize]);
        if sum > 0 {
            new_q = q_idx + 3;
            score = score.wrapping_add(sum);
            sum = 0;
        }
        if sum < x {
            break;
        }
        q_idx += 4;
        s_idx += 1;
    }

    if score >= reduced_cutoff {
        return extend_hit_ungapped_exact_ncbi(
            q_seq,
            s_seq_packed,
            q_off,
            s_off,
            subject_len,
            x_dropoff,
            matrix,
        );
    }

    let s_match_len = s_match_end.saturating_sub(s_start);
    let right_len = if new_q >= q_start as isize {
        (new_q - q_start as isize + 1) as usize
    } else {
        0
    };
    let length = s_match_len.max(right_len);

    UngappedData {
        q_start,
        s_start,
        length,
        score,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:286,292-305,336-341
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1122-1127
// The approximate extension works on the concatenated query: its left loop is bounded by
// the offset in the whole query, and the truncated `q_byte` reads the sentinel (0x0F) as
// bases, so it runs into the previous context. The exact extension stops at the sentinel.
#[cfg(test)]
#[path = "extension_approx_context_tests.rs"]
mod approx_context_tests;

/// Determine word type for two-stage lookup
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:508-607 (s_TypeOfWord function)
///
/// This is a FAITHFUL TRANSPILE of NCBI BLAST's s_TypeOfWord function.
/// It does NOT perform actual sequence extension - only calculates ext_to
/// and checks masked regions.
///
/// # Arguments
/// * `q_off` - Query offset (position of lut_word_length match start) - WILL BE MODIFIED if masked regions found
/// * `s_off` - Subject offset (position of lut_word_length match start) - WILL BE MODIFIED if masked regions found
/// * `has_locations` - Whether NCBI's `locations` (the lookup table's `masked_locations`)
///   is set: the checks of the lookup words run only then
/// * `q_range` - Exclusive end offset of the current query context
/// * `s_range` - Exclusive end offset of the current subject scan range
/// * `word_length` - Full word length (e.g., 28 for megablast)
/// * `lut_word_length` - Lookup table word length (e.g., 8)
/// * `check_double` - Whether to check for double word (default: true)
/// * `is_seed_masked` - Lookup-based seed check (s_IsSeedMasked)
///
/// # Returns
/// `(word_type, extended, q_off, s_off)` where:
/// - `word_type`: 0 = non-word, 1 = single word, 2 = double word
/// - `extended`: Number of bases extended (set to ext_to, NOT actual sequence checking)
/// - `q_off`, `s_off`: possibly adjusted offsets after masking checks
///
/// Note: reward and penalty parameters are NOT used (NCBI BLAST does not check actual matches)
pub fn type_of_word<F>(
    q_off: usize,
    s_off: usize,
    has_locations: bool,
    q_range: usize,
    s_range: usize,
    word_length: usize,
    lut_word_length: usize,
    check_double: bool,
    is_seed_masked: &mut F,
) -> (u8, usize, usize, usize)
where
    F: FnMut(usize, usize) -> bool,
{
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:508-607
    // FAITHFUL TRANSPILE - DO NOT ASSUME

    // NCBI: *extended = 0;
    let mut extended = 0;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:531
    if word_length == lut_word_length {
        return (1, 0, q_off, s_off);
    }

    // NCBI: Int4 q_end = *q_off + word_length;
    // NCBI: Int4 s_end = *s_off + word_length;
    let mut q_end = q_off + word_length;
    let mut s_end = s_off + word_length;

    // NCBI: if (locations) { ... }
    // Check masked regions and adjust q_off/s_off if needed
    let mut q_off_adjusted = q_off;
    let mut s_off_adjusted = s_off;
    if has_locations {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:452-490

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:540-544
        if is_seed_masked(s_end - lut_word_length, q_end - lut_word_length) {
            return (0, 0, q_off_adjusted, s_off_adjusted);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:546-549
        while q_off_adjusted < q_range && s_off_adjusted < s_range {
            if !is_seed_masked(s_off_adjusted, q_off_adjusted) {
                break;
            }
            q_off_adjusted += 1;
            s_off_adjusted += 1;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:546-549
        // Rust must stop if the absolute context/subject boundary is reached
        // before finding an unmasked seed.
        if q_off_adjusted >= q_range || s_off_adjusted >= s_range {
            return (0, 0, q_off_adjusted, s_off_adjusted);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:553-554
    let ext_to = word_length - (q_end - q_off_adjusted);
    // NCBI's `Int4 ext_max` is negative for a word that runs past the context or the scan
    // range (the aligned small-table extension does not check the right end at the end of
    // the query); then `ext_to > ext_max`, and no extension or double word is possible.
    let ext_max = (q_range as i64 - q_end as i64).min(s_range as i64 - s_end as i64);
    if ext_max < 0 {
        return if ext_to > 0 || has_locations {
            (0, 0, q_off_adjusted, s_off_adjusted)
        } else {
            (1, 0, q_off_adjusted, s_off_adjusted)
        };
    }
    let ext_max = ext_max as usize;

    // NCBI: if (ext_to || locations) {
    if ext_to > 0 || has_locations {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:559
        if ext_to > ext_max {
            return (0, 0, q_off_adjusted, s_off_adjusted);
        }

        // NCBI: q_end += ext_to;
        // NCBI: s_end += ext_to;
        q_end += ext_to;
        s_end += ext_to;

        // NCBI: for (s_pos = s_end - lut_word_length, q_pos = q_end - lut_word_length;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:563-569
        let mut s_pos = s_end - lut_word_length;
        let mut q_pos = q_end - lut_word_length;

        // NCBI BLAST: Loop condition is s_pos > *s_off (no underflow check)
        // We use saturating_sub to prevent underflow while maintaining NCBI behavior
        while s_pos > s_off_adjusted {
            if is_seed_masked(s_pos, q_pos) {
                return (0, 0, q_off_adjusted, s_off_adjusted);
            }
            // NCBI BLAST: s_pos -= lut_word_length (no bounds check)
            // Use saturating_sub to prevent underflow while maintaining loop termination
            s_pos = s_pos.saturating_sub(lut_word_length);
            q_pos = q_pos.saturating_sub(lut_word_length);
        }

        // NCBI: (*extended) = ext_to;
        extended = ext_to;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:575-576
    if !check_double {
        return (1, extended, q_off_adjusted, s_off_adjusted);
    }

    // NCBI: ext_to += word_length;
    let double_ext_to = ext_to + word_length;
    // NCBI: ext_max = MIN(ext_max, ext_to);
    let double_ext_max = ext_max.min(double_ext_to);

    // NCBI: try seed by seed
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:585-586
    let mut s_pos = s_end;
    let mut q_pos = q_end;

    // NCBI BLAST: Loop runs regardless of locations status
    // If locations is NULL, s_IsSeedMasked returns FALSE, so loop continues to ext_max
    while extended + lut_word_length <= double_ext_max {
        // NCBI BLAST: s_IsSeedMasked returns FALSE if locations is NULL
        if is_seed_masked(s_pos, q_pos) {
            break;
        }
        extended += lut_word_length;
        s_pos += lut_word_length;
        q_pos += lut_word_length;
    }

    // NCBI: try base by base
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:595-596
    // NCBI BLAST: Position adjustment happens UNCONDITIONALLY (before the while loop)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:594-603
    // NCBI BLAST: s_pos -= (lut_word_length - 1) (unconditional, no bounds check)
    // Use saturating_sub to prevent underflow while maintaining NCBI behavior
    s_pos = s_pos.saturating_sub(lut_word_length - 1);
    q_pos = q_pos.saturating_sub(lut_word_length - 1);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:600-602
    // NCBI BLAST: This loop runs regardless of locations status
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:597-603
    while extended < double_ext_max {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:597-599
        if is_seed_masked(s_pos, q_pos) {
            return (1, extended, q_off_adjusted, s_off_adjusted);
        }
        extended += 1;
        s_pos += 1;
        q_pos += 1;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:605
    if double_ext_max == double_ext_to {
        (2, extended, q_off_adjusted, s_off_adjusted)
    } else {
        (1, extended, q_off_adjusted, s_off_adjusted)
    }
}
