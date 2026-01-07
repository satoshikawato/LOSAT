use super::constants::X_DROP_UNGAPPED;
use crate::utils::dust::MaskedInterval;

/// Extend an ungapped hit in both directions using X-drop termination.
/// 
/// NCBI reference: na_ungapped.c:148-243 (s_NuclUngappedExtendExact)
/// 
/// FAITHFUL TRANSPILE of NCBI BLAST's ungapped extension algorithm.
/// 
/// # NCBI Implementation Details:
/// - X is negative and when sum becomes more negative than X we break out of the loop
/// - Left extension: sum += matrix[*--q][...], if (sum > 0) { score += sum; sum = 0; } else if (sum < X) break;
/// - Right extension: sum += matrix[*q++][...], if (sum > 0) { score += sum; X_current = (-score > X) ? -score : X; sum = 0; } else if (sum < X_current) break;
/// 
/// # Index Semantics
/// - Left extension starts at i=0 (the seed position itself), extending leftward
/// - Right extension starts at j=1 (position after seed), avoiding double-counting
/// - The seed position (q_pos, s_pos) is scored once in the left extension loop
///
/// # Arguments
/// * `q_seq` - Query sequence (uncompressed, unlike NCBI which uses compressed)
/// * `s_seq` - Subject sequence (uncompressed, unlike NCBI which uses compressed)
/// * `q_pos` - Seed position in query
/// * `s_pos` - Seed position in subject
/// * `reward` - Score for match (typically positive, e.g., 2)
/// * `penalty` - Score for mismatch (typically negative, e.g., -3)
/// * `x_drop` - X-drop threshold for termination (optional, uses default if None)
///              NOTE: NCBI passes negative value: -(cutoffs->x_dropoff)
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
    // NCBI reference: na_ungapped.c:744
    // s_NuclUngappedExtendExact(..., -(cutoffs->x_dropoff), ungapped_data);
    // CRITICAL: X is negative in NCBI implementation
    // LOSAT receives positive value, so we negate it to match NCBI
    let x_drop_positive = x_drop.unwrap_or(X_DROP_UNGAPPED);
    let x = -x_drop_positive; // NCBI: X is negative

    // NCBI reference: na_ungapped.c:177-178
    // score = 0;
    // sum = 0;
    let mut score = 0i32;
    let mut sum = 0i32;

    // NCBI reference: na_ungapped.c:188-204 (Left extension)
    // /* extend to the left */
    // while ((s > start) || (s == start && base < remainder)) {
    //     if (base == 3) { s--; base = 0; } else { ++base; }
    //     ch = *s;
    //     if ((sum += matrix[*--q][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
    //         q_beg = q;
    //         score += sum;
    //         sum = 0;
    //     } else if (sum < X) {
    //         break;
    //     }
    // }
    let mut q_beg = q_pos;
    let mut i = 0;
    let max_left = q_pos.min(s_pos);

    while i <= max_left {
        let q = unsafe { *q_seq.get_unchecked(q_pos - i) };
        let s = unsafe { *s_seq.get_unchecked(s_pos - i) };
        let sc = if q == s { reward } else { penalty };
        
        // NCBI: sum += matrix[*--q][NCBI2NA_UNPACK_BASE(ch, base)]
        sum += sc;
        
        // NCBI: if ((sum += ...) > 0) {
        //     q_beg = q;
        //     score += sum;
        //     sum = 0;
        // } else if (sum < X) {
        //     break;
        // }
        if sum > 0 {
            q_beg = q_pos - i;
            score += sum;
            sum = 0;
        } else if sum < x {
            // NCBI reference: na_ungapped.c:201-202
            // else if (sum < X) {
            //     break;
            // }
            break;
        }
        i += 1;
    }

    // NCBI reference: na_ungapped.c:206-207
    // ungapped_data->q_start = (Int4)(q_beg - query->sequence);
    // ungapped_data->s_start = s_off - (q_off - ungapped_data->q_start);
    let q_start = q_beg;
    let s_start = s_pos - (q_pos - q_start);

    // NCBI reference: na_ungapped.c:217-239 (Right extension)
    // /* extend to the right */
    // q = query->sequence + q_off;
    // s = subject0 + s_off / COMPRESSION_RATIO;
    // sum = 0;
    // base = 3 - (s_off % COMPRESSION_RATIO);
    //
    // /* X_current is used to break out of loop if score goes negative */
    // Int4 X_current=X;
    // while (s < sf || (s == sf && base > remainder)) {
    //     ch = *s;
    //     if ((sum += matrix[*q++][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
    //         q_end = q;
    //         score += sum;
    //         X_current = (-score > X) ? -score : X;
    //         sum = 0;
    //     } else if (sum < X_current)
    //         break;
    //     ...
    // }
    let mut q_end = q_pos;
    sum = 0;
    
    // NCBI reference: na_ungapped.c:223-224
    // /* X_current is used to break out of loop if score goes negative */
    // Int4 X_current=X;
    let mut x_current = x;
    
    // NCBI reference: na_ungapped.c:218-219
    // /* extend to the right */
    // q = query->sequence + q_off;
    // s = subject0 + s_off / COMPRESSION_RATIO;
    // CRITICAL: Right extension starts at q_off (seed position), NOT q_off+1
    // NCBI uses *q++ which means it scores q_off first, then moves right
    let mut j = 0; // NCBI: Right extension starts at seed position (q_off)
    while (q_pos + j) < q_seq.len() && (s_pos + j) < s_seq.len() {
        let q = unsafe { *q_seq.get_unchecked(q_pos + j) };
        let s = unsafe { *s_seq.get_unchecked(s_pos + j) };
        let sc = if q == s { reward } else { penalty };
        
        // NCBI: sum += matrix[*q++][NCBI2NA_UNPACK_BASE(ch, base)]
        // *q++ means: use current q, then increment
        // So seed position (j=0) IS scored in right extension
        sum += sc;
        
        // NCBI: if ((sum += ...) > 0) {
        //     q_end = q;
        //     score += sum;
        //     X_current = (-score > X) ? -score : X;
        //     sum = 0;
        // } else if (sum < X_current)
        //     break;
        if sum > 0 {
            // NCBI reference: na_ungapped.c:227-228
            // if ((sum += matrix[*q++][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
            //     q_end = q;
            // CRITICAL: *q++ means use current q, then increment
            // So q_end = q points to position AFTER the last scored position
            // We need to use q_pos + j + 1 to match NCBI behavior
            q_end = q_pos + j + 1;
            score += sum;
            // NCBI reference: na_ungapped.c:230
            // X_current = (-score > X) ? -score : X;
            x_current = if -score > x { -score } else { x };
            sum = 0;
        } else if sum < x_current {
            // NCBI reference: na_ungapped.c:232-233
            // else if (sum < X_current)
            //     break;
            break;
        }
        j += 1;
    }

    // NCBI reference: na_ungapped.c:241-242
    // ungapped_data->length = (Int4)(q_end - q_beg);
    // ungapped_data->score = score;
    // CRITICAL: q_end already points to position AFTER last scored position (due to *q++)
    // So q_end - q_beg is the length (inclusive of both endpoints)
    // For half-open interval [start, end), we use q_end directly
    let q_end_pos = q_end; // q_end already points to position after last scored
    let s_end_pos = s_start + (q_end_pos - q_start);

    (
        q_start,
        q_end_pos,
        s_start,
        s_end_pos,
        score,
    )
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
/// `(word_type, extended)` where:
/// - `word_type`: 0 = non-word, 1 = single word, 2 = double word
/// - `extended`: Number of bases extended (set to ext_to, NOT actual sequence checking)
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
) -> (u8, usize) {
    // NCBI reference: na_ungapped.c:508-607
    // FAITHFUL TRANSPILE - DO NOT ASSUME
    
    // NCBI: *extended = 0;
    let mut extended = 0;
    
    // NCBI: if (word_length == lut_word_length) return 1;
    if word_length == lut_word_length {
        return (1, 0);
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
        
        // NCBI: check for right end first
        // if (s_IsSeedMasked(lookup_wrap, subject, s_end - lut_word_length, lut_word_length, q_end - lut_word_length)) return 0;
        if is_kmer_masked(query_mask, q_end - lut_word_length, lut_word_length) {
            return (0, 0);
        }
        
        // NCBI: search for valid left end and reposition q_off
        // for (; TRUE; ++(*s_off), ++(*q_off)) {
        //     if (!s_IsSeedMasked(lookup_wrap, subject, *s_off, lut_word_length, *q_off)) break;
        // }
        while q_off_adjusted < q_seq.len() && s_off_adjusted < s_seq.len() {
            if !is_kmer_masked(query_mask, q_off_adjusted, lut_word_length) {
                break;
            }
            q_off_adjusted += 1;
            s_off_adjusted += 1;
        }
        
        // Recalculate q_end and s_end after adjustment
        q_end = q_off_adjusted + word_length;
        s_end = s_off_adjusted + word_length;
    }
    
    // NCBI: ext_to = word_length - (q_end - (*q_off));
    let ext_to = word_length - (q_end - q_off_adjusted);
    
    // NCBI: ext_max = MIN(q_range - q_end, s_range - s_end);
    let ext_max = (q_range - q_end).min(s_range - s_end);
    
    // NCBI: if (ext_to || locations) {
    if ext_to > 0 || !query_mask.is_empty() {
        // NCBI: if (ext_to > ext_max) return 0;
        if ext_to > ext_max {
            return (0, 0);
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
        // NCBI BLAST: This loop runs ONLY if locations (query_mask) is not NULL
        // Reference: na_ungapped.c:563-570
        if !query_mask.is_empty() {
            use crate::algorithm::blastn::lookup::is_kmer_masked;
            let mut s_pos = s_end - lut_word_length;
            let mut q_pos = q_end - lut_word_length;
            
            // NCBI BLAST: Loop condition is s_pos > *s_off (no underflow check)
            // We use saturating_sub to prevent underflow while maintaining NCBI behavior
            while s_pos > s_off_adjusted {
                if is_kmer_masked(query_mask, q_pos, lut_word_length) {
                    return (0, 0);
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
    
    // NCBI: if (!check_double) return 1;
    if !check_double {
        return (1, extended);
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
        if !query_mask.is_empty() && is_kmer_masked(query_mask, q_pos, lut_word_length) {
            return (1, extended);
        }
        extended += 1;
        s_pos += 1;
        q_pos += 1;
    }
    
    // NCBI: return ((ext_max == ext_to) ? 2 : 1);
    if double_ext_max == double_ext_to {
        return (2, extended);
    } else {
        return (1, extended);
    }
}
