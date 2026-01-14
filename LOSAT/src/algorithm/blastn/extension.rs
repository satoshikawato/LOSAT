use super::constants::X_DROP_UNGAPPED;
use crate::utils::dust::MaskedInterval;
use crate::core::blast_encoding::COMPRESSION_RATIO;

/// BLASTNA alphabet size.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_encoding.h:91-92
/// ```c
/// #define BLASTNA_SIZE 16     /**< Size of nucleic acid alphabet */
/// ```
const BLASTNA_SIZE: usize = 16;

/// Ungapped extension results (0-based coordinates, length in bases).
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:141-147
/// ```c
/// typedef struct BlastUngappedData {
///    Int4 q_start;
///    Int4 s_start;
///    Int4 length;
///    Int4 score;
/// } BlastUngappedData;
/// ```
#[derive(Clone, Copy, Debug)]
pub struct UngappedData {
    pub q_start: usize,
    pub s_start: usize,
    pub length: usize,
    pub score: i32,
}

/// Build the nucleotide ungapped score table for 4-base blocks.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:236-259
/// ```c
/// for (i = 0; i < 256; i++) {
///    Int4 score = 0;
///    if (i & 3) score += penalty; else score += reward;
///    if ((i >> 2) & 3) score += penalty; else score += reward;
///    if ((i >> 4) & 3) score += penalty; else score += reward;
///    if (i >> 6) score += penalty; else score += reward;
///    table[i] = score;
/// }
/// ```
pub fn build_nucl_score_table(reward: i32, penalty: i32) -> [i32; 256] {
    let mut table = [0i32; 256];
    for i in 0..256 {
        let mut score = 0;
        if (i & 3) != 0 { score += penalty; } else { score += reward; }
        if ((i >> 2) & 3) != 0 { score += penalty; } else { score += reward; }
        if ((i >> 4) & 3) != 0 { score += penalty; } else { score += reward; }
        if (i >> 6) != 0 { score += penalty; } else { score += reward; }
        table[i] = score;
    }
    table
}

/// Extract base N from packed ncbi2na byte.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_util.h:51-55
/// ```c
/// #define NCBI2NA_UNPACK_BASE(x, N) (((x)>>(2*(N))) & NCBI2NA_MASK)
/// ```
#[inline]
fn ncbi2na_unpack_base(byte: u8, n: u8) -> u8 {
    (byte >> (2 * n)) & 0x03
}

/// Exact ungapped extension (per-base) for nucleotide sequences.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:148-242
/// ```c
/// s_NuclUngappedExtendExact(..., Int4 q_off, Int4 s_off, Int4 X, ...)
/// {
///    base = 3 - (s_off % 4);
///    ...
///    while ((s > start) || (s == start && base < remainder)) {
///        if (base == 3) { s--; base = 0; } else { ++base; }
///        if ((sum += matrix[*--q][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
///            q_beg = q; score += sum; sum = 0;
///        } else if (sum < X) break;
///    }
///    ...
///    while (s < sf || (s == sf && base > remainder)) {
///        if ((sum += matrix[*q++][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
///            q_end = q; score += sum; X_current = (-score > X) ? -score : X; sum = 0;
///        } else if (sum < X_current) break;
///        if (base == 0) { base = 3; s++; } else base--;
///    }
///    ungapped_data->length = (Int4)(q_end - q_beg);
/// }
/// ```
pub fn extend_hit_ungapped_exact_ncbi(
    q_seq: &[u8],             // BLASTNA (1 byte/base)
    s_seq_packed: &[u8],      // ncbi2na packed (4 bases/byte)
    q_off: usize,
    s_off: usize,
    subject_len: usize,
    x_dropoff: i32,
    matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
) -> UngappedData {
    let x = -x_dropoff;
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

        let ch = s_seq_packed[s_idx as usize];
        let s_base = ncbi2na_unpack_base(ch, base as u8) as usize;
        let q_code = q_seq[q_idx as usize] as usize;
        sum += matrix[q_code * BLASTNA_SIZE + s_base];
        if sum > 0 {
            q_beg = q_idx;
            score += sum;
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
        let ch = s_seq_packed[s_idx as usize];
        let s_base = ncbi2na_unpack_base(ch, base as u8) as usize;
        let q_code = q_seq[q_idx as usize] as usize;
        q_idx += 1;

        sum += matrix[q_code * BLASTNA_SIZE + s_base];
        if sum > 0 {
            q_end = q_idx;
            score += sum;
            x_current = if -score > x { -score } else { x };
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

/// Approximate ungapped extension (4-base blocks) with optional exact recomputation.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:261-349
/// ```c
/// s_NuclUngappedExtend(..., Int4 q_off, Int4 s_match_end, Int4 s_off, Int4 X, ...)
/// {
///    len = (4 - (s_off % 4)) % 4;
///    q_ext = q_off + len; s_ext = s_off + len;
///    ...
///    for (i = 0; i < len; s--, q -= 4, i++) {
///        Uint1 s_byte = s[-1];
///        Uint1 q_byte = (q[-4] << 6) | (q[-3] << 4) | (q[-2] << 2) | q[-1];
///        sum += score_table[q_byte ^ s_byte];
///        ...
///    }
///    ...
///    if (score >= reduced_cutoff) {
///        s_NuclUngappedExtendExact(...);
///    } else {
///        ungapped_data->length = MAX(s_match_end - s_start, ...);
///    }
/// }
/// ```
pub fn extend_hit_ungapped_approx_ncbi(
    q_seq: &[u8],             // BLASTNA (1 byte/base)
    s_seq_packed: &[u8],      // ncbi2na packed (4 bases/byte)
    q_off: usize,
    s_off: usize,
    s_match_end: usize,
    subject_len: usize,
    x_dropoff: i32,
    score_table: &[i32; 256],
    reduced_cutoff: i32,
    matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
) -> UngappedData {
    let x = -x_dropoff;
    let q_len = q_seq.len();

    let len = (COMPRESSION_RATIO - (s_off % COMPRESSION_RATIO)) % COMPRESSION_RATIO;
    let q_ext = q_off + len;
    let s_ext = s_off + len;

    let mut q_idx = q_ext as isize;
    let mut s_idx = (s_ext / COMPRESSION_RATIO) as isize;
    let left_len = (q_ext.min(s_ext) / COMPRESSION_RATIO) as usize;

    let mut score = 0i32;
    let mut sum = 0i32;
    let mut new_q = q_idx;

    // Left extension in 4-base blocks.
    for _ in 0..left_len {
        if q_idx < 4 || s_idx == 0 {
            break;
        }
        let s_byte = s_seq_packed[(s_idx - 1) as usize];
        let q_byte = (q_seq[(q_idx - 4) as usize] << 6)
            | (q_seq[(q_idx - 3) as usize] << 4)
            | (q_seq[(q_idx - 2) as usize] << 2)
            | q_seq[(q_idx - 1) as usize];
        sum += score_table[(q_byte ^ s_byte) as usize];
        if sum > 0 {
            new_q = q_idx - 4;
            score += sum;
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
    let right_len = ((q_len - q_ext).min(subject_len - s_ext) / COMPRESSION_RATIO) as usize;
    sum = 0;
    new_q = q_idx;

    for _ in 0..right_len {
        let s_byte = s_seq_packed[s_idx as usize];
        let q_byte = (q_seq[q_idx as usize] << 6)
            | (q_seq[(q_idx + 1) as usize] << 4)
            | (q_seq[(q_idx + 2) as usize] << 2)
            | q_seq[(q_idx + 3) as usize];
        sum += score_table[(q_byte ^ s_byte) as usize];
        if sum > 0 {
            new_q = q_idx + 3;
            score += sum;
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

/// Extend an ungapped hit in both directions using X-drop termination.
/// 
/// # Index Semantics
/// - Left extension starts at i=0 (the seed position itself), extending leftward
/// - Right extension starts at j=1 (position after seed), avoiding double-counting
/// - The seed position (q_pos, s_pos) is scored once in the left extension loop
///
/// This follows NCBI BLAST's na_ungapped.c pattern where the extension includes
/// the starting position and extends outward in both directions.
///
/// # Arguments
/// * `q_seq` - Query sequence
/// * `s_seq` - Subject sequence
/// * `q_pos` - Seed position in query
/// * `s_pos` - Seed position in subject
/// * `reward` - Score for match (typically positive, e.g., 2)
/// * `penalty` - Score for mismatch (typically negative, e.g., -3)
/// * `x_drop` - X-drop threshold for termination (optional, uses default if None)
///
/// # Returns
/// (q_start, q_end, s_start, s_end, score) - Half-open interval [start, end)
pub fn extend_hit_ungapped(
    q_seq: &[u8],
    s_seq: &[u8],
    q_pos: usize,
    s_pos: usize,
    reward: i32,
    penalty: i32,
    x_drop: Option<i32>,
) -> (usize, usize, usize, usize, i32) {
    let x_drop_val = x_drop.unwrap_or(X_DROP_UNGAPPED);
    let mut current_score = 0;
    let mut max_score = 0;

    // Left extension: starts at i=0 (seed position) and extends leftward
    // When i=0: accesses q_seq[q_pos] and s_seq[s_pos] (the seed position)
    // When i=1: accesses q_seq[q_pos-1] and s_seq[s_pos-1] (one left of seed)
    let mut best_i = 0;
    let mut i = 0;
    let max_left = q_pos.min(s_pos);

    while i <= max_left {
        let q = unsafe { *q_seq.get_unchecked(q_pos - i) };
        let s = unsafe { *s_seq.get_unchecked(s_pos - i) };
        let sc = if q == s { reward } else { penalty };
        current_score += sc;

        if current_score > max_score {
            max_score = current_score;
            best_i = i;
        } else if (max_score - current_score) > x_drop_val {
            break;
        }
        i += 1;
    }

    // Right extension: starts at j=1 to avoid double-counting seed position
    // When j=1: accesses q_seq[q_pos+1] and s_seq[s_pos+1]
    let mut current_score_r = max_score;
    let mut max_score_total = max_score;
    let mut best_j = 0;
    let mut j = 1;

    while (q_pos + j) < q_seq.len() && (s_pos + j) < s_seq.len() {
        let q = unsafe { *q_seq.get_unchecked(q_pos + j) };
        let s = unsafe { *s_seq.get_unchecked(s_pos + j) };
        let sc = if q == s { reward } else { penalty };
        current_score_r += sc;

        if current_score_r >= max_score_total {
            max_score_total = current_score_r;
            best_j = j;
        } else if (max_score_total - current_score_r) > x_drop_val {
            // NCBI BLAST uses only X-drop termination, NOT score <= 0 check
            // Reference: na_ungapped.c:220-233
            break;
        }
        j += 1;
    }

    let result = (
        q_pos - best_i,
        q_pos + best_j + 1,
        s_pos - best_i,
        s_pos + best_j + 1,
        max_score_total,
    );

    // Debug coordinate tracking
    if std::env::var("LOSAT_DEBUG_COORDS").is_ok() {
        eprintln!(
            "[UNGAPPED] q_pos={}, s_pos={}, best_i={}, best_j={} -> q=[{}, {}), s=[{}, {}), score={}",
            q_pos, s_pos, best_i, best_j, result.0, result.1, result.2, result.3, result.4
        );
    }

    result
}

/// Determine word type for two-stage lookup
/// 
/// NCBI reference: na_ungapped.c:508-607 (s_TypeOfWord function)
/// 
/// This is a FAITHFUL TRANSPILE of NCBI BLAST's s_TypeOfWord function.
/// It does NOT perform actual sequence extension - only calculates ext_to
/// and checks masked regions.
/// 
/// # Arguments
/// * `q_seq` - Query sequence
/// * `s_seq` - Subject sequence
/// * `q_off` - Query offset (position of lut_word_length match start) - WILL BE MODIFIED if masked regions found
/// * `s_off` - Subject offset (position of lut_word_length match start) - WILL BE MODIFIED if masked regions found
/// * `query_mask` - Masked intervals for query sequence
/// * `word_length` - Full word length (e.g., 28 for megablast)
/// * `lut_word_length` - Lookup table word length (e.g., 8)
/// * `check_double` - Whether to check for double word (default: true)
/// 
/// # Returns
/// `(word_type, extended, q_off, s_off)` where:
/// - `word_type`: 0 = non-word, 1 = single word, 2 = double word
/// - `extended`: Number of bases extended (set to ext_to, NOT actual sequence checking)
/// - `q_off`, `s_off`: possibly adjusted offsets after masking checks
/// 
/// Note: reward and penalty parameters are NOT used (NCBI BLAST does not check actual matches)
pub fn type_of_word(
    q_seq: &[u8],
    s_seq: &[u8],
    q_off: usize,
    s_off: usize,
    query_mask: &[MaskedInterval],
    word_length: usize,
    lut_word_length: usize,
    check_double: bool,
) -> (u8, usize, usize, usize) {
    // NCBI reference: na_ungapped.c:508-607
    // FAITHFUL TRANSPILE - DO NOT ASSUME
    
    // NCBI: *extended = 0;
    let mut extended = 0;
    
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:531
    // ```c
    // if (word_length == lut_word_length) return 1;
    // ```
    if word_length == lut_word_length {
        return (1, 0, q_off, s_off);
    }
    
    // NCBI: Int4 q_end = *q_off + word_length;
    // NCBI: Int4 s_end = *s_off + word_length;
    let mut q_end = q_off + word_length;
    let mut s_end = s_off + word_length;
    
    // NCBI: context = BSearchContextInfo(q_end, query_info);
    // NCBI: q_range = query_info->contexts[context].query_offset + query_info->contexts[context].query_length;
    // For blastn, we typically have a single context, so q_range = q_seq.len()
    let q_range = q_seq.len();
    let s_range = s_seq.len();
    
    // NCBI: if (locations) { ... }
    // Check masked regions and adjust q_off/s_off if needed
    let mut q_off_adjusted = q_off;
    let mut s_off_adjusted = s_off;
    
    if !query_mask.is_empty() {
        use crate::algorithm::blastn::lookup::is_kmer_masked;
        
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:540-544
        // ```c
        // if (s_IsSeedMasked(lookup_wrap, subject,
        //                    s_end - lut_word_length,
        //                    lut_word_length,
        //                    q_end - lut_word_length)) return 0;
        // ```
        if is_kmer_masked(query_mask, q_end - lut_word_length, lut_word_length) {
            return (0, 0, q_off_adjusted, s_off_adjusted);
        }
        
        // NCBI: search for valid left end and reposition q_off
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:546-549
        // ```c
        // for (; TRUE; ++(*s_off), ++(*q_off)) {
        //     if (!s_IsSeedMasked(lookup_wrap, subject,
        //                         *s_off, lut_word_length, *q_off)) break;
        // }
        // ```
        while q_off_adjusted < q_seq.len() && s_off_adjusted < s_seq.len() {
            if !is_kmer_masked(query_mask, q_off_adjusted, lut_word_length) {
                break;
            }
            q_off_adjusted += 1;
            s_off_adjusted += 1;
        }
    }
    
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:553-554
    // ```c
    // ext_to = word_length - (q_end - (*q_off));
    // ext_max = MIN(q_range - q_end, s_range - s_end);
    // ```
    let ext_to = word_length - (q_end - q_off_adjusted);
    let ext_max = (q_range - q_end).min(s_range - s_end);
    
    // NCBI: if (ext_to || locations) {
    if ext_to > 0 || !query_mask.is_empty() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:559
        // ```c
        // if (ext_to > ext_max) return 0;
        // ```
        if ext_to > ext_max {
            return (0, 0, q_off_adjusted, s_off_adjusted);
        }
        
        // NCBI: q_end += ext_to;
        // NCBI: s_end += ext_to;
        q_end += ext_to;
        s_end += ext_to;
        
        // NCBI: for (s_pos = s_end - lut_word_length, q_pos = q_end - lut_word_length;
        //      s_pos > *s_off; 
        //      s_pos -= lut_word_length, q_pos -= lut_word_length) {
        //     if (s_IsSeedMasked(lookup_wrap, subject, s_pos, lut_word_length, q_pos)) return 0;
        // }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:563-569
        // ```c
        // for (s_pos = s_end - lut_word_length,
        //      q_pos = q_end - lut_word_length;
        //      s_pos > *s_off;
        //      s_pos -= lut_word_length,
        //      q_pos -= lut_word_length) {
        //     if (s_IsSeedMasked(lookup_wrap, subject,
        //                       s_pos, lut_word_length, q_pos)) return 0;
        // }
        // ```
        if !query_mask.is_empty() {
            use crate::algorithm::blastn::lookup::is_kmer_masked;
            let mut s_pos = s_end - lut_word_length;
            let mut q_pos = q_end - lut_word_length;
            
            // NCBI BLAST: Loop condition is s_pos > *s_off (no underflow check)
            // We use saturating_sub to prevent underflow while maintaining NCBI behavior
            while s_pos > s_off_adjusted {
                if is_kmer_masked(query_mask, q_pos, lut_word_length) {
                    return (0, 0, q_off_adjusted, s_off_adjusted);
                }
                // NCBI BLAST: s_pos -= lut_word_length (no bounds check)
                // Use saturating_sub to prevent underflow while maintaining loop termination
                s_pos = s_pos.saturating_sub(lut_word_length);
                q_pos = q_pos.saturating_sub(lut_word_length);
            }
        }
        
        // NCBI: (*extended) = ext_to;
        extended = ext_to;
    }
    
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:575-576
    // ```c
    // if (!check_double) return 1;
    // ```
    if !check_double {
        return (1, extended, q_off_adjusted, s_off_adjusted);
    }
    
    // NCBI: ext_to += word_length;
    let double_ext_to = ext_to + word_length;
    // NCBI: ext_max = MIN(ext_max, ext_to);
    let double_ext_max = ext_max.min(double_ext_to);
    
    // NCBI: try seed by seed
    // for (s_pos = s_end, q_pos = q_end; 
    //      *extended + lut_word_length <= ext_max; 
    //      s_pos += lut_word_length, q_pos += lut_word_length, (*extended) += lut_word_length) {
    //     if (s_IsSeedMasked(lookup_wrap, subject, s_pos, lut_word_length, q_pos)) break;
    // }
    // NCBI BLAST: This loop runs regardless of locations (query_mask) status
    // Reference: na_ungapped.c:584-592
    // When locations is NULL, s_IsSeedMasked always returns FALSE, so the loop continues
    use crate::algorithm::blastn::lookup::is_kmer_masked;
    let mut s_pos = s_end;
    let mut q_pos = q_end;
    
    // NCBI BLAST: Loop runs regardless of locations status
    // If locations is NULL, s_IsSeedMasked returns FALSE, so loop continues to ext_max
    while extended + lut_word_length <= double_ext_max {
        // NCBI BLAST: s_IsSeedMasked returns FALSE if locations is NULL
        if !query_mask.is_empty() && is_kmer_masked(query_mask, q_pos, lut_word_length) {
            break;
        }
        extended += lut_word_length;
        s_pos += lut_word_length;
        q_pos += lut_word_length;
    }
    
    // NCBI: try base by base
    // s_pos -= (lut_word_length - 1);
    // q_pos -= (lut_word_length - 1);
    // NCBI BLAST: Position adjustment happens UNCONDITIONALLY (before the while loop)
    // Reference: na_ungapped.c:594-603
    // NCBI BLAST: s_pos -= (lut_word_length - 1) (unconditional, no bounds check)
    // Use saturating_sub to prevent underflow while maintaining NCBI behavior
    s_pos = s_pos.saturating_sub(lut_word_length - 1);
    q_pos = q_pos.saturating_sub(lut_word_length - 1);
    
    // NCBI: while (*extended < ext_max) {
    //     if (s_IsSeedMasked(lookup_wrap, subject, s_pos, lut_word_length, q_pos)) return 1;
    //     (*extended)++;
    //     ++s_pos;
    //     ++q_pos;
    // }
    // NCBI BLAST: This loop runs regardless of locations status
    // Reference: na_ungapped.c:597-603
    while extended < double_ext_max {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:597-599
        // ```c
        // if (s_IsSeedMasked(lookup_wrap, subject, s_pos,
        //                        lut_word_length, q_pos)) return 1;
        // ```
        if !query_mask.is_empty() && is_kmer_masked(query_mask, q_pos, lut_word_length) {
            return (1, extended, q_off_adjusted, s_off_adjusted);
        }
        extended += 1;
        s_pos += 1;
        q_pos += 1;
    }
    
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:605
    // ```c
    // return ((ext_max == ext_to) ? 2 : 1);
    // ```
    if double_ext_max == double_ext_to {
        return (2, extended, q_off_adjusted, s_off_adjusted);
    } else {
        return (1, extended, q_off_adjusted, s_off_adjusted);
    }
}
