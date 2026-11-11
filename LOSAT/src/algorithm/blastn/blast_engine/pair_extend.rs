//! Steps of one offset pair of a blastn subject chunk's scan, as `search_query_batch`'s
//! `collect_prelim_hits_for_chunk` runs them for each pair, in the order of the pairs, on one
//! thread: the word extension of a lookup table hit of a two-stage table (`TwoStageLookup`) to
//! `word_length` letters (`word_extend_pair`, in the closure `process_offset_pairs`), and the
//! off-diagonal search of a single word in two-hit mode (`off_diagonal_found`, called by that
//! closure and by the scan of the other tables, `na_lookup` and `pv_direct_lookup`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1143-1150,1160-1161 (the end of
//! the loop of `s_BlastNaExtend`, :1052-1174; `s_BlastSmallNaExtend`, :1486-1593, ends the same
//! way for the small-query table), which hands each extended pair to
//! `s_BlastnDiagHashExtendInitialHit` (:799-944, the diagonal hash) or
//! `s_BlastnDiagTableExtendInitialHit` (:631-775, the diagonal array; both search the off
//! diagonals of a single word, :685-717 and :852-886)
//! ```c
//!         q_offset -= ext_left;
//!         s_offset -= ext_left;
//!         /* check the diagonal on which the hit lies. The boundaries
//!            extend from the first match of the hit to one beyond the last
//!            match */
//!
//!         if (word_params->container_type == eDiagHash) {
//!             hits_extended += s_BlastnDiagHashExtendInitialHit(query, subject,
//! ...
//!         } else {
//!             hits_extended += s_BlastnDiagTableExtendInitialHit(query, subject,
//! ```
//!
//! This module is a child of `run` for the types it keeps private (`DiagHashTable`). Its
//! functions are statements of the scan's closures that were there, in their order; the strips
//! folds of `x_pair_seed.rs` and `x_pair_seed_disc.rs` copy those statements one by one, so their
//! order is kept here. What the word extension reads for every pair of a chunk is
//! `PairExtendInputs`, built once per chunk; what changes from pair to pair is `LutHit` and
//! `OffDiagonal`. Nothing here runs on the pool. Both functions are `#[inline(always)]`, so the
//! scan's loops keep the shape they had (callgrind, n1). The ungapped extension and the save of
//! the initial hit stay in `process_offset_pairs`: as a function of their own, the 16 values of
//! the hit they read crossed a struct on every hit and cost `-task blastn -ungapped` 0.3 % more
//! instructions even inlined.

use super::{DiagHashTable, DiagStruct, SmallNaWord, COMPRESSION_RATIO};

/// What the word extension of one offset pair reads that is the same for every pair of one
/// subject chunk: the small-query table's form (`small_na_word`, `None` for the other tables),
/// whether its aligned one-byte form applies, whether the chunk is masked, the word and lookup word
/// lengths, the packed subject, the query with its sentinels, and the switch of the debug counters.
/// Built once per chunk by `collect_prelim_hits_for_chunk` before `process_offset_pairs`, with the
/// names there; every field is a shared reference or a `Copy` value. The other steps of
/// `process_offset_pairs` read the closure's captures.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1052-1060,1064,1081
/// ```c
/// s_BlastNaExtend(const BlastOffsetPair * offset_pairs, Int4 num_hits,
///                 const BlastInitialWordParameters * word_params,
///                 LookupTableWrap * lookup_wrap,
///                 BLAST_SequenceBlk * query,
///                 BLAST_SequenceBlk * subject, Int4 ** matrix,
///                 BlastQueryInfo * query_info,
///                 Blast_ExtendWord * ewp,
///                 BlastInitHitList * init_hitlist,
///                 Uint4 s_range)
/// ...
///     Int4 word_length, lut_word_length, ext_to;
/// ...
///     ext_to = word_length - lut_word_length;
/// ```
/// The debug switch has no NCBI counterpart (diagnostics only).
pub(super) struct PairExtendInputs<'a> {
    pub(super) small_na_word: &'a Option<SmallNaWord<'a>>,
    pub(super) small_na_aligned_one_byte: bool,
    pub(super) subject_masked: bool,
    pub(super) word_length: usize,
    pub(super) lut_word_length: usize,
    pub(super) search_seq_packed: &'a [u8],
    pub(super) encoded_query_concat_blastna_with_sentinels: &'a [u8],
    pub(super) debug_enabled: bool,
}

/// One lookup table hit as the word extension reads it: its query offset in the concatenated
/// query (`q_off0`) and its subject offset in the chunk (`kmer_start`), the range of its context
/// in the concatenated query, and the end of the subject range it was scanned in (`s_range`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1094-1095,1511-1513
/// ```c
///         Int4 s_offset = offset_pairs[index].qs_offsets.s_off;
///         Int4 q_offset = offset_pairs[index].qs_offsets.q_off;
/// ...
///         Int4 context = BSearchContextInfo(q_offset, query_info);
///         Int4 q_start = query_info->contexts[context].query_offset;
///         Int4 q_range = q_start + query_info->contexts[context].query_length;
/// ```
pub(super) struct LutHit {
    pub(super) q_off0: usize,
    pub(super) kmer_start: usize,
    pub(super) q_context_start: usize,
    pub(super) q_context_end: usize,
    pub(super) s_range: usize,
}

/// The word extension of one lookup table hit of a two-stage table: with the small-query table
/// the four-letters-at-a-time extension (`SmallNaWord`); otherwise, when the word is longer than
/// the lookup word, the hit is extended letter by letter on the packed subject to the left and
/// then to the right until `word_length` letters match. `None` when they do not (the pair is
/// dropped); otherwise the query and subject offsets of the word's start and the left and right
/// extensions (the trace prints these two).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1093-1095,1120-1124,1138-1144
/// (`s_BlastNaExtend`; `s_BlastSmallNaExtend` at :1504-1566 for the small-query table)
/// ```c
///     for (; index < num_hits; ++index) {
///         Int4 s_offset = offset_pairs[index].qs_offsets.s_off;
///         Int4 q_offset = offset_pairs[index].qs_offsets.q_off;
/// ...
///         if (ext_left < ext_to) {
///             Int4 ext_right = 0;
///             s_off = s_offset + lut_word_length;
///             if (s_off + ext_to - ext_left > s_range)
///                 continue;
/// ...
///             /* check if enough extra matches were found */
///             if (ext_left + ext_right < ext_to)
///                 continue;
///         }
///
///         q_offset -= ext_left;
///         s_offset -= ext_left;
/// ```
/// `#[inline(always)]`: the scan calls it for every lookup table hit; out of line it took 35
/// instructions a call (callgrind, n1: dc-megablast +1.4 %, with the off-diagonal search's call),
/// so it keeps the shape it had inside `process_offset_pairs`.
#[inline(always)]
pub(super) fn word_extend_pair(
    x: &PairExtendInputs<'_>,
    hit: LutHit,
    dbg_left_ext_iters: &mut usize,
    dbg_right_ext_iters: &mut usize,
) -> Option<(usize, usize, usize, usize)> {
    let &PairExtendInputs {
        small_na_word,
        small_na_aligned_one_byte,
        subject_masked,
        word_length,
        lut_word_length,
        search_seq_packed,
        encoded_query_concat_blastna_with_sentinels,
        debug_enabled,
        ..
    } = x;
    let LutHit {
        q_off0,
        kmer_start,
        q_context_start,
        q_context_end,
        s_range,
    } = hit;
    let (q_ext_start, s_ext_start, word_ext_left, word_ext_right) = if let Some(small_na_word) =
        &small_na_word
    {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1651-1663
        // ```c
        //     /* if sequence is masked, fall back to generic scanner and extender */
        //     if (subject->mask_type != eNoSubjMasking) {
        //         ...
        //             if (extend != (TNaExtendFunction)s_BlastNaExtendDirect) {
        //                  extend = (lookup_wrap->lut_type == eSmallNaLookupTable)
        //                     ? (TNaExtendFunction)s_BlastSmallNaExtend
        //                     : (TNaExtendFunction)s_BlastNaExtend;
        //             }
        // ```
        let word = if small_na_aligned_one_byte && !subject_masked {
            small_na_word.extend_aligned_one_byte(
                q_off0,
                kmer_start,
                q_context_start,
                q_context_end,
                s_range,
            )
        } else {
            small_na_word.extend(q_off0, kmer_start, q_context_start, q_context_end, s_range)
        };
        let (q_word, s_word) = word?;
        (q_word, s_word, q_off0.saturating_sub(q_word), 0)
    } else if word_length > lut_word_length {
        // NCBI BLAST: s_BlastNaExtend left/right extension uses packed subject
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1093-1144
        // ```c
        // Int4 ext_left = 0;
        // Int4 s_off = s_offset;
        // Uint1 *q = query->sequence + q_offset;
        // Uint1 *s = subject->sequence + s_off / COMPRESSION_RATIO;
        // for (; ext_left < MIN(ext_to, s_offset); ++ext_left) {
        //     s_off--;
        //     q--;
        //     if (s_off % COMPRESSION_RATIO == 3)
        //         s--;
        //     if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >> 6) != *q)
        //         break;
        // }
        // if (ext_left < ext_to) {
        //     Int4 ext_right = 0;
        //     s_off = s_offset + lut_word_length;
        //     if (s_off + ext_to - ext_left > s_range) continue;
        //     q = query->sequence + q_offset + lut_word_length;
        //     s = subject->sequence + s_off / COMPRESSION_RATIO;
        //     for (; ext_right < ext_to - ext_left; ++ext_right) {
        //         if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >> 6) != *q)
        //             break;
        //         s_off++;
        //         q++;
        //         if (s_off % COMPRESSION_RATIO == 0)
        //             s++;
        //     }
        //     if (ext_left + ext_right < ext_to) continue;
        // }
        let ext_to = word_length - lut_word_length;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1101-1114
        // ```c
        // Int4 ext_left = 0;
        // Int4 s_off = s_offset;
        // Uint1 *q = query->sequence + q_offset;
        // Uint1 *s = subject->sequence + s_off / COMPRESSION_RATIO;
        // for (; ext_left < MIN(ext_to, s_offset); ++ext_left) {
        //     s_off--;
        //     q--;
        //     if (s_off % COMPRESSION_RATIO == 3)
        //         s--;
        //     if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >> 6) != *q)
        //         break;
        // }
        // ```
        let max_ext_left = ext_to.min(kmer_start).min(q_off0.saturating_add(1));
        let mut ext_left = 0usize;
        let mut q_left = q_off0 + 1;
        let mut s_off = kmer_start;
        let mut s_idx = s_off / COMPRESSION_RATIO;
        while ext_left < max_ext_left {
            if debug_enabled {
                *dbg_left_ext_iters += 1;
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1101-1114
            // ```c
            // Uint1 *q = query->sequence + q_offset;
            // ...
            // s_off--;
            // q--;
            // ```
            // LOSAT's sentinel query buffer stores logical query offset
            // `n` at `n + 1`, so after NCBI's pre-decrement the first
            // left comparison at q_offset - 1 is buffer index q_offset.
            q_left -= 1;
            s_off -= 1;
            if s_off % COMPRESSION_RATIO == COMPRESSION_RATIO - 1 {
                s_idx -= 1;
            }
            // SAFETY: s_idx tracks s_off / COMPRESSION_RATIO within bounds by max_ext_left.
            let s_byte = unsafe { *search_seq_packed.get_unchecked(s_idx) };
            let s_base = (s_byte << (2 * (s_off % COMPRESSION_RATIO))) >> 6;
            // SAFETY: q_left is bounded by the outer-sentinel query buffer.
            let q_base =
                unsafe { *encoded_query_concat_blastna_with_sentinels.get_unchecked(q_left) };
            if s_base != q_base {
                break;
            }
            ext_left += 1;
        }

        // RIGHT extension (forwards from lut_word_length end)
        // NCBI BLAST: na_ungapped.c:1120-1136
        // if (ext_left < ext_to) {
        //     Int4 ext_right = 0;
        //     s_off = s_offset + lut_word_length;
        //     if (s_off + ext_to - ext_left > s_range) continue;
        //     q = query->sequence + q_offset + lut_word_length;
        //     s = subject->sequence + s_off / COMPRESSION_RATIO;
        //     for (; ext_right < ext_to - ext_left; ++ext_right) {
        //         if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >> 6) != *q)
        //             break;
        //         s_off++;
        //         q++;
        //         if (s_off % COMPRESSION_RATIO == 0)
        //             s++;
        //     }
        // }
        // Only do right extension if left didn't get all bases
        // NCBI reference: na_ungapped.c:1120-1136
        // if (ext_left < ext_to) {
        //     Int4 ext_right = 0;
        //     s_off = s_offset + lut_word_length;
        //     if (s_off + ext_to - ext_left > s_range) continue;
        //     q = query->sequence + q_offset + lut_word_length;
        //     s = subject->sequence + s_off / COMPRESSION_RATIO;
        //     for (; ext_right < ext_to - ext_left; ++ext_right) {
        //         if (mismatch) break;
        //         s_off++;
        //         q++;
        //     }
        // }
        let mut ext_right = 0usize;
        if ext_left < ext_to {
            // NCBI BLAST: s_off = s_offset + lut_word_length;
            // NCBI BLAST: if (s_off + ext_to - ext_left > s_range) continue;
            // Reference: na_ungapped.c:1122-1124
            let mut s_off = kmer_start + lut_word_length;
            if s_off + (ext_to - ext_left) > s_range {
                // Not enough room in subject, skip this seed (NCBI behavior)
                return None;
            }

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1120-1136
            // ```c
            // q = query->sequence + q_offset + lut_word_length;
            // ...
            // if (packed_subject_base != *q)
            //     break;
            // s_off++;
            // q++;
            // ```
            // The right comparison checks q_offset + lut_word_length
            // before incrementing q.
            let mut q_right = q_off0 + lut_word_length + 1;
            if q_right + (ext_to - ext_left) > encoded_query_concat_blastna_with_sentinels.len() {
                return None;
            }
            let mut s_idx = s_off / COMPRESSION_RATIO;

            // NCBI BLAST: for (; ext_right < ext_to - ext_left; ++ext_right)
            // Reference: na_ungapped.c:1128-1136
            while ext_right < (ext_to - ext_left) {
                if debug_enabled {
                    *dbg_right_ext_iters += 1;
                }
                // NCBI BLAST: if (base mismatch) break;
                // Reference: na_ungapped.c:1129-1131
                // SAFETY: s_idx tracks s_off / COMPRESSION_RATIO within bounds by s_len check above.
                let s_byte = unsafe { *search_seq_packed.get_unchecked(s_idx) };
                let s_base = (s_byte << (2 * (s_off % COMPRESSION_RATIO))) >> 6;
                // SAFETY: q_right is bounded by the outer-sentinel query buffer.
                let q_base =
                    unsafe { *encoded_query_concat_blastna_with_sentinels.get_unchecked(q_right) };
                if s_base != q_base {
                    break;
                }
                ext_right += 1;
                q_right += 1;
                s_off += 1;
                if s_off.is_multiple_of(COMPRESSION_RATIO) {
                    s_idx += 1;
                }
            }

            // NCBI BLAST: if (ext_left + ext_right < ext_to) continue;
            // Reference: na_ungapped.c:1125
            if ext_left + ext_right < ext_to {
                // Word_length match failed, skip this seed
                return None;
            }
        }

        // Adjust positions for type_of_word and ungapped extension
        // NCBI BLAST: q_offset -= ext_left; s_offset -= ext_left;
        // Reference: na_ungapped.c:1143-1144
        (
            q_off0 - ext_left,
            kmer_start - ext_left,
            ext_left,
            ext_right,
        )
    } else {
        // word_length == lut_word_length, no extension needed
        (q_off0, kmer_start, 0, 0)
    };
    Some((q_ext_start, s_ext_start, word_ext_left, word_ext_right))
}

/// What the off-diagonal search of one seed reads: the diagonal table's kind (array or hash), the
/// window, the off-diagonal range (`scan_range`, NCBI's `scan_range` of the word options) and the
/// word length, the array's length, mask and size, the array's levels and lengths and the hash
/// (read only), and the seed's diagonal (masked for the array), start and end in the table's
/// coordinates. Built at each call, from the names of the caller (the na branch passes `safe_k`
/// as the word length).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:658,663-669
/// ```c
///     Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
/// ...
///     diag = s_off + diag_table->diag_array_length - q_off;
///     real_diag = diag & diag_table->diag_mask;
///     last_hit = hit_level_array[real_diag].last_hit;
///     hit_saved = hit_level_array[real_diag].flag;
///     s_end = s_off + word_length;
///     s_off_pos = s_off + diag_table->offset;
///     s_end_pos = s_end + diag_table->offset;
/// ```
pub(super) struct OffDiagonal<'d> {
    pub(super) use_array_indexing: bool,
    pub(super) window_size: usize,
    pub(super) scan_range: usize,
    pub(super) word_length: usize,
    pub(super) diag_array_length: isize,
    pub(super) diag_mask: isize,
    pub(super) diag_array_size: usize,
    pub(super) hit_level_array: &'d [DiagStruct],
    pub(super) hit_len_array: &'d [u8],
    pub(super) diag_hash: &'d DiagHashTable,
    pub(super) diag: isize,
    pub(super) s_off_pos: usize,
    pub(super) s_end_pos: usize,
}

/// The off-diagonal search of a single word (`word_type == 1`) in two-hit mode: whether a
/// diagonal at most `Delta` away from the seed's (`Delta` is `scan_range`, at most `window_size -
/// word_length`) holds a first hit (a hit length that is not 0) whose end and start pass NCBI's
/// window tests (`s_a`, `s_b`), in the diagonal array or in the hash. The caller sets `hit_ready`
/// to false when none does (a new hit).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:685-712
/// (`s_BlastnDiagTableExtendInitialHit`, the diagonal array)
/// ```c
///         /* for single word, also try off diagonals */
///         if (word_type == 1) {
///             /* try off-diagonals */
///             Int4 orig_diag = real_diag + diag_table->diag_array_length;
///             Int4 s_a = s_off_pos + word_length - window_size;
///             Int4 s_b = s_end_pos - 2 * word_length;
///             Int4 delta;
///             if (Delta < 0) Delta = 0;
///             for (delta = 1; delta <= Delta ; ++delta) {
///                 Int4 off_diag  = (orig_diag + delta) & diag_table->diag_mask;
///                 Int4 off_s_end = hit_level_array[off_diag].last_hit;
///                 Int4 off_s_l   = diag_table->hit_len_array[off_diag];
///                 if ( off_s_l
///                  && off_s_end - delta >= s_a
///                  && off_s_end - off_s_l <= s_b) {
///                     off_found = TRUE;
///                     break;
///                 }
///                 off_diag  = (orig_diag - delta) & diag_table->diag_mask;
///                 off_s_end = hit_level_array[off_diag].last_hit;
///                 off_s_l   = diag_table->hit_len_array[off_diag];
///                 if ( off_s_l
///                  && off_s_end >= s_a
///                  && off_s_end - off_s_l + delta <= s_b) {
///                     off_found = TRUE;
///                     break;
///                 }
///             }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:852-881
/// (`s_BlastnDiagHashExtendInitialHit`, the diagonal hash)
/// ```c
///         /* for single word, also try off diagonals */
///         if (word_type == 1) {
///             /* try off-diagonal */
///             Int4 s_a = s_off_pos + word_length - window_size;
///             Int4 s_b = s_end_pos - 2 * word_length;
///             Int4 delta;
///             if (Delta < 0) Delta = 0;
///             for (delta = 1; delta <= Delta; ++delta) {
///                 Int4 off_s_end = 0;
///                 Int4 off_s_l = 0;
///                 Int4 off_hit_saved = 0;
///                 Int4 off_rc = s_BlastDiagHashRetrieve(hash_table, diag + delta,
///                               &off_s_end, &off_s_l, &off_hit_saved);
///                 if ( off_rc
///                   && off_s_l
///                   && off_s_end - delta >= s_a
///                   && off_s_end - off_s_l <= s_b) {
///                      off_found = TRUE;
///                      break;
///                 }
///                 off_rc = s_BlastDiagHashRetrieve(hash_table, diag - delta,
///                               &off_s_end, &off_s_l, &off_hit_saved);
///                 if ( off_rc
///                   && off_s_l
///                   && off_s_end >= s_a
///                   && off_s_end - off_s_l + delta <= s_b) {
///                      off_found = TRUE;
///                      break;
///                 }
///             }
/// ```
/// `#[inline(always)]`: the scan calls it for every single word in two-hit mode; out of line it
/// took 38 instructions a call with `scan_range` 0 (callgrind, n1 dc-megablast), so it keeps the
/// shape it had inside the scan's closures.
#[inline(always)]
pub(super) fn off_diagonal_found(od: OffDiagonal<'_>) -> bool {
    let OffDiagonal {
        use_array_indexing,
        window_size,
        scan_range,
        word_length,
        diag_array_length,
        diag_mask,
        diag_array_size,
        hit_level_array,
        hit_len_array,
        diag_hash,
        diag,
        s_off_pos,
        s_end_pos,
    } = od;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:657
    // ```c
    //     Boolean off_found = FALSE;
    // ```
    let mut off_found = false;
    // NCBI reference: na_ungapped.c:658, 692
    // Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
    // if (Delta < 0) Delta = 0;
    let delta_calc = window_size as isize - word_length as isize;
    let delta_max = if delta_calc < 0 {
        0
    } else {
        scan_range.min(delta_calc as usize) as isize
    };

    // NCBI reference: na_ungapped.c:689-690
    // Int4 s_a = s_off_pos + word_length - window_size;
    // Int4 s_b = s_end_pos - 2 * word_length;
    // CRITICAL: NCBI uses signed arithmetic (Int4), so s_a and s_b can be negative
    // LOSAT must use signed arithmetic to match NCBI behavior
    let s_a = s_off_pos as isize + word_length as isize - window_size as isize;
    let s_b = s_end_pos as isize - 2 * word_length as isize;

    if use_array_indexing {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:688-694
        // ```c
        // orig_diag = real_diag + diag_table->diag_array_length;
        // off_diag  = (orig_diag + delta) & diag_table->diag_mask;
        // ```
        let orig_diag = diag + diag_array_length;
        // NCBI reference: na_ungapped.c:693
        // for (delta = 1; delta <= Delta ; ++delta) {
        for delta in 1..=delta_max {
            // NCBI reference: na_ungapped.c:694-702
            // Int4 off_diag  = (orig_diag + delta) & diag_table->diag_mask;
            // Int4 off_s_end = hit_level_array[off_diag].last_hit;
            // Int4 off_s_l   = diag_table->hit_len_array[off_diag];
            // if ( off_s_l
            //  && off_s_end - delta >= s_a
            //  && off_s_end - off_s_l <= s_b) {
            //     off_found = TRUE;
            //     break;
            // }
            let off_diag = (orig_diag + delta) & diag_mask;
            let (off_s_end, off_s_l) = {
                let off_diag_idx = off_diag as usize;
                if off_diag_idx < diag_array_size {
                    let off_entry = &hit_level_array[off_diag_idx];
                    (off_entry.last_hit, hit_len_array[off_diag_idx])
                } else {
                    (0, 0)
                }
            };
            // NCBI: off_s_end - delta >= s_a (signed comparison)
            // Convert to signed for comparison to match NCBI behavior
            if off_s_l > 0
                && (off_s_end as isize - delta) >= s_a
                && (off_s_end as isize - off_s_l as isize) <= s_b
            {
                off_found = true;
                break;
            }

            // NCBI reference: na_ungapped.c:703-711
            // off_diag  = (orig_diag - delta) & diag_table->diag_mask;
            // off_s_end = hit_level_array[off_diag].last_hit;
            // off_s_l   = diag_table->hit_len_array[off_diag];
            // if ( off_s_l
            //  && off_s_end >= s_a
            //  && off_s_end - off_s_l + delta <= s_b) {
            //     off_found = TRUE;
            //     break;
            // }
            let off_diag = (orig_diag - delta) & diag_mask;
            let (off_s_end, off_s_l) = {
                let off_diag_idx = off_diag as usize;
                if off_diag_idx < diag_array_size {
                    let off_entry = &hit_level_array[off_diag_idx];
                    (off_entry.last_hit, hit_len_array[off_diag_idx])
                } else {
                    (0, 0)
                }
            };
            // NCBI: off_s_end >= s_a (signed comparison)
            // Convert to signed for comparison to match NCBI behavior
            if off_s_l > 0
                && (off_s_end as isize) >= s_a
                && (off_s_end as isize - off_s_l as isize + delta) <= s_b
            {
                off_found = true;
                break;
            }
        }
    } else {
        let diag_i32 = diag as i32;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:859-874
        // ```c
        // off_rc = s_BlastDiagHashRetrieve(hash_table, diag + delta,
        //           &off_s_end, &off_s_l, &off_hit_saved);
        // ...
        // off_rc = s_BlastDiagHashRetrieve(hash_table, diag - delta,
        //           &off_s_end, &off_s_l, &off_hit_saved);
        // ```
        for delta in 1..=delta_max {
            if let Some((off_s_end, off_s_l, _)) = diag_hash.retrieve(diag_i32 + delta as i32) {
                if off_s_l > 0
                    && (off_s_end as isize - delta) >= s_a
                    && (off_s_end as isize - off_s_l as isize) <= s_b
                {
                    off_found = true;
                    break;
                }
            }
            if let Some((off_s_end, off_s_l, _)) = diag_hash.retrieve(diag_i32 - delta as i32) {
                if off_s_l > 0
                    && (off_s_end as isize) >= s_a
                    && (off_s_end as isize - off_s_l as isize + delta) <= s_b
                {
                    off_found = true;
                    break;
                }
            }
        }
    }
    off_found
}
