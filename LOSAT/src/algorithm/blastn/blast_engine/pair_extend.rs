//! The steps of the offset pairs of a blastn subject chunk's two-stage scan, as
//! `two_stage_scan::scan_two_stage_lookup` (the two-stage branch of
//! `prelim_subject::prelim_chunk_hits`) runs them for each batch of pairs its scan collects, in the
//! order of the pairs, on one thread (`extend_offset_pairs`): for each pair the word extension of
//! the lookup table hit of a two-stage table (`TwoStageLookup`) to
//! `word_length` letters (`word_extend_pair`), the diagonal's last hit, the word type, the
//! off-diagonal search of a single word in two-hit mode (`off_diagonal_found`, which the scan of
//! the other tables, `na_lookup` and `pv_direct_lookup`, calls too: `na_scan::scan_na_lookup`),
//! the ungapped extension and the save of the initial hit.
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
//! This module is a child of `run` for the items it keeps private (the types `DiagHashTable`,
//! `OffsetPair`, `SubjectChunk`, `QueryContext`, `QueryContextIndex`, `UngappedHit` and
//! `BlastnTiming`; `x_push_ungapped`, `mask_lookup_index`, `packed_kmer_at_seed_mask` and
//! `X_HIT_AHEAD`). Its functions are the statements of the scan's closures that were there, in
//! their order; the strips folds of `x_pair_seed.rs` and `x_pair_seed_disc.rs` copy those
//! statements one by one, so their order is kept here. What the steps read that is the same for
//! every pair of a chunk is `PairExtendInputs`, and what they write is `PairExtendState`, both
//! built once per chunk; what changes from pair to pair is `LutHit` and `OffDiagonal`. Nothing
//! here runs on the pool. `word_extend_pair` and `off_diagonal_found` are `#[inline(always)]`,
//! so the loop of `extend_offset_pairs` keeps the shape it had as a closure (callgrind, n1). The
//! ungapped extension and the save of the initial hit are statements of that loop, whose values
//! they read as locals: as a function of their own, the 16 values of the hit they read crossed a
//! struct on every hit and cost `-task blastn -ungapped` 0.3 % more instructions even inlined.

use super::x_diag_lv::LvSeeds;
use super::{
    blastn_trace, extend_hit_ungapped_approx_ncbi, extend_hit_ungapped_exact_ncbi,
    mask_lookup_index, packed_kmer_at_seed_mask, type_of_word, x_dc_csr, x_push_ungapped,
    BlastnTiming, DiagHashTable, DiagStruct, OffsetPair, QueryContext, QueryContextIndex,
    SmallNaWord, SubjectChunk, UngappedHit, COMPRESSION_RATIO, X_HIT_AHEAD,
};
use crate::algorithm::blastn::alignment::gapped::BlastnaMatrix;
use crate::algorithm::blastn::lookup::TwoStageLookup;

/// What the steps of the offset pairs of one subject chunk read that is the same for every pair:
/// for the word extension the small-query table's form (`small_na_word`, `None` for the other
/// tables), whether its aligned one-byte form applies, whether the chunk is masked, the word and
/// lookup word lengths, the packed subject, the query with its sentinels and the switch of the
/// debug counters; then what the other steps read, in groups below. Built once per chunk by
/// `two_stage_scan::scan_two_stage_lookup` before the scan, with the names there; every field is a
/// shared reference or a `Copy` value.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1052-1060,1064,1081
/// (`s_BlastNaExtend`) and :651,656,658 (`s_BlastnDiagTableExtendInitialHit`; the hash's
/// function has the same lines at :819,823,825)
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
/// ...
///     Int4 window_size = word_params->options->window_size;
/// ...
///     Boolean two_hits = (window_size > 0);
/// ...
///     Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
/// ```
/// The prefetch, the level bound and the diagnostics have no NCBI counterpart.
pub(super) struct PairExtendInputs<'a> {
    /// The word extension (`word_extend_pair`).
    pub(super) small_na_word: &'a Option<SmallNaWord<'a>>,
    pub(super) small_na_aligned_one_byte: bool,
    pub(super) subject_masked: bool,
    pub(super) word_length: usize,
    pub(super) lut_word_length: usize,
    pub(super) search_seq_packed: &'a [u8],
    pub(super) encoded_query_concat_blastna_with_sentinels: &'a [u8],
    pub(super) debug_enabled: bool,
    /// The prefetch of the query letters of a later pair (LOSAT_X_MBHITPF), the context of a
    /// pair's query offset, the contexts and their cutoffs, the chunk's length and the chunk (its
    /// offset in the subject), and the debug window.
    pub(super) x_hit_prefetch: bool,
    pub(super) query_context_index: &'a QueryContextIndex,
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) cutoff_scores: &'a [i32],
    pub(super) s_len: usize,
    pub(super) chunk: &'a SubjectChunk,
    pub(super) debug_window: Option<(usize, usize, usize, usize)>,
    /// The diagonal table: whether it is the array (else the hash), the array's length, mask and
    /// size, the table's offset for this chunk, and the level bound of LOSAT_X_DIAGLV (whether
    /// the retrieve takes it, and the constants of the chunk's seeds).
    pub(super) use_array_indexing: bool,
    pub(super) diag_array_length: isize,
    pub(super) diag_mask: isize,
    pub(super) diag_array_size: usize,
    pub(super) diag_offset: isize,
    pub(super) x_lv_ref: bool,
    pub(super) x_lv_seeds: LvSeeds,
    /// The word type and the window: the table (the lookups of `s_TypeOfWord`), whether a hit
    /// takes two words, the window, the off-diagonal range, and the window of a hash insert.
    pub(super) two_stage: &'a TwoStageLookup,
    pub(super) two_hits: bool,
    pub(super) window_size: usize,
    pub(super) scan_range: usize,
    pub(super) diag_hash_window: i32,
    /// The ungapped extension: the concatenated query, the X-dropoff and the reduced cutoff of
    /// each context, the score matrix, and the score table of the approximate extension.
    pub(super) encoded_query_concat_blastna: &'a [u8],
    pub(super) x_dropoff_scores: &'a [i32],
    pub(super) reduced_cutoff_scores: &'a [i32],
    pub(super) score_matrix: &'a BlastnaMatrix,
    pub(super) nucl_score_table: &'a [i32; 256],
    /// The diagnostics: the trace lines (`blastn_trace`) with the subject's index and trace id,
    /// the debug window's lines with the minimum ungapped score they print, and LOSAT_TIMING.
    pub(super) blastn_trace_enabled: bool,
    pub(super) s_idx: usize,
    pub(super) s_id: &'a str,
    pub(super) debug_mode: bool,
    pub(super) min_ungapped_score: i32,
    pub(super) timing_enabled: bool,
    pub(super) timing_ref: Option<&'a BlastnTiming>,
}

/// What the steps write: the chunk's initial hits, the diagonal table (the array's levels and
/// lengths, or the hash), and the debug counters (the seeds in the debug window, the word
/// extension's left and right steps, the ungapped extensions, the hits under the cutoff and the
/// hits saved). Built once per chunk by `two_stage_scan::scan_two_stage_lookup` after
/// `PairExtendInputs`, from its scratch and its counters (references only, the types the steps had
/// as a closure), and held by the scan's callbacks for the whole scan.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1058-1059 (`s_BlastNaExtend`),
/// :641-642 (`s_BlastnDiagTableExtendInitialHit`) and :809 (`s_BlastnDiagHashExtendInitialHit`)
/// ```c
///                 Blast_ExtendWord * ewp,
///                 BlastInitHitList * init_hitlist,
/// ...
///                              BLAST_DiagTable * diag_table,
///                              BlastInitHitList * init_hitlist,
/// ...
///                                BLAST_DiagHash * hash_table,
/// ```
/// The counters have no NCBI counterpart (diagnostics only).
pub(super) struct PairExtendState<'s> {
    /// The chunk's initial hits (NCBI's `init_hitlist`).
    pub(super) ungapped_hits: &'s mut Vec<UngappedHit>,
    /// The diagonal array (a query block of at most 8000 letters) or the diagonal hash.
    pub(super) hit_level_array: &'s mut Vec<DiagStruct>,
    pub(super) hit_len_array: &'s mut Vec<u8>,
    pub(super) diag_hash: &'s mut DiagHashTable,
    /// The debug counters (written only when `debug_enabled`).
    pub(super) dbg_window_seeds: &'s mut usize,
    pub(super) dbg_left_ext_iters: &'s mut usize,
    pub(super) dbg_right_ext_iters: &'s mut usize,
    pub(super) dbg_ungapped_ext_calls: &'s mut usize,
    pub(super) dbg_ungapped_low: &'s mut usize,
    pub(super) dbg_gapped_attempted: &'s mut usize,
}

/// The steps of the offset pairs that the scan of one chunk with a two-stage table
/// (`TwoStageLookup`) collected, in their order: for each pair its context, the bounds test, the
/// word extension (`word_extend_pair`), the diagonal's last hit and the `last_hit` test, the word
/// type (`type_of_word`: in two-hit mode with `check_double`, and for a single word the
/// off-diagonal search, `off_diagonal_found`), the diagonal update of a hit that is not ready,
/// the ungapped extension, the cutoff test, and the save of the initial hit with its diagonal
/// update. Then the batch is empty (`*offset_pairs_len = 0`). The scan calls it when its batch
/// is full and once after the scan; the time of the ungapped extensions is taken off the scan's
/// time (`scan_pause_ns`) while the scan runs (`track_scan_pause`), not at the last call, which
/// comes after the scan's time is taken.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1093-1095,1149-1173
/// (`s_BlastNaExtend`, :1052-1174; the functions it calls are `s_BlastnDiagHashExtendInitialHit`,
/// :799-944, and `s_BlastnDiagTableExtendInitialHit`, :631-775)
/// ```c
///     for (; index < num_hits; ++index) {
///         Int4 s_offset = offset_pairs[index].qs_offsets.s_off;
///         Int4 q_offset = offset_pairs[index].qs_offsets.q_off;
/// ...
///         if (word_params->container_type == eDiagHash) {
///             hits_extended += s_BlastnDiagHashExtendInitialHit(query, subject,
/// ...
///         } else {
///             hits_extended += s_BlastnDiagTableExtendInitialHit(query, subject,
/// ...
///         }
///     }
///     return hits_extended;
/// ```
/// No `#[inline]`: as a closure of `collect_prelim_hits_for_chunk` (now
/// `prelim_subject::prelim_chunk_hits`) the steps were an out-of-line function already (the scan
/// calls them from 11 places).
pub(super) fn extend_offset_pairs(
    x: &PairExtendInputs<'_>,
    state: &mut PairExtendState<'_>,
    offset_pairs: &[OffsetPair],
    offset_pairs_len: &mut usize,
    scan_pause_ns: &mut u64,
    track_scan_pause: bool,
) {
    // `small_na_aligned_one_byte` and `subject_masked` are read by `word_extend_pair` only.
    let &PairExtendInputs {
        small_na_word,
        word_length,
        lut_word_length,
        search_seq_packed,
        encoded_query_concat_blastna_with_sentinels,
        debug_enabled,
        x_hit_prefetch,
        query_context_index,
        query_contexts,
        cutoff_scores,
        s_len,
        chunk,
        debug_window,
        use_array_indexing,
        diag_array_length,
        diag_mask,
        diag_array_size,
        diag_offset,
        x_lv_ref,
        x_lv_seeds,
        two_stage,
        two_hits,
        window_size,
        scan_range,
        diag_hash_window,
        encoded_query_concat_blastna,
        x_dropoff_scores,
        reduced_cutoff_scores,
        score_matrix,
        nucl_score_table,
        blastn_trace_enabled,
        s_idx,
        s_id,
        debug_mode,
        min_ungapped_score,
        timing_enabled,
        timing_ref,
        ..
    } = x;
    let ungapped_hits = &mut *state.ungapped_hits;
    let hit_level_array = &mut *state.hit_level_array;
    let hit_len_array = &mut *state.hit_len_array;
    let diag_hash = &mut *state.diag_hash;
    let dbg_window_seeds = &mut *state.dbg_window_seeds;
    let dbg_left_ext_iters = &mut *state.dbg_left_ext_iters;
    let dbg_right_ext_iters = &mut *state.dbg_right_ext_iters;
    let dbg_ungapped_ext_calls = &mut *state.dbg_ungapped_ext_calls;
    let dbg_ungapped_low = &mut *state.dbg_ungapped_low;
    let dbg_gapped_attempted = &mut *state.dbg_gapped_attempted;
    let hitsfound = *offset_pairs_len;
    for idx in 0..hitsfound {
        // EXPERIMENT (LOSAT_X_MBHITPF): the query's letters of a hit
        // `X_HIT_AHEAD` places on.
        if x_hit_prefetch && idx + X_HIT_AHEAD < hitsfound {
            x_dc_csr::prefetch(
                encoded_query_concat_blastna_with_sentinels
                    .as_ptr()
                    .wrapping_add(offset_pairs[idx + X_HIT_AHEAD].q_off),
            );
        }
        let pair = offset_pairs[idx];
        let q_off0 = pair.q_off;
        let kmer_start = pair.s_off;
        let s_range = pair.s_range;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:730-733
        // ```c
        // Int4 context = BSearchContextInfo(q_off, query_info);
        // ```
        let context_idx = query_context_index.context_for_offset(q_off0);
        let ctx = &query_contexts[context_idx];
        let q_idx = context_idx as u32;
        let q_context_start = ctx.query_offset as usize;
        let q_context_end = q_context_start + ctx.len;
        let q_pos_usize = q_off0 - q_context_start;

        // Use pre-computed cutoff score (computed once per query-subject pair)
        let cutoff_score = cutoff_scores[context_idx];

        // NCBI reference: na_ungapped.c:1081-1140
        // CRITICAL: In two-stage lookup, seed finding phase ONLY finds lut_word_length matches.
        // The word_length matching happens LATER in s_BlastnExtendInitialHit (extension phase).
        //
        // NCBI BLAST flow:
        // 1. Seed finding: Find all lut_word_length matches (this phase)
        // 2. Extension: For each seed, call s_BlastnExtendInitialHit which:
        //    - Does LEFT extension (backwards from lut_word_length start)
        //    - Does RIGHT extension (forwards from lut_word_length end)
        //    - Verifies word_length match (ext_left + ext_right >= ext_to)
        //
        // We should NOT filter seeds during seed finding - pass ALL seeds to extension.
        // The extension phase will handle word_length matching.

        // Skip only if sequences are too short for lut_word_length match
        // (word_length check happens in extension phase)
        if q_off0 + lut_word_length > q_context_end || kmer_start + lut_word_length > s_len {
            continue;
        }

        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:663-664
        // ```c
        // diag = s_off + diag_table->diag_array_length - q_off;
        // real_diag = diag & diag_table->diag_mask;
        // ```
        // This value is used only for seed trace output before the
        // word-length extension below adjusts q_off/s_off.
        let lookup_diag = if use_array_indexing {
            (kmer_start as isize + diag_array_length - q_off0 as isize) & diag_mask
        } else {
            kmer_start as isize - q_off0 as isize
        };

        // Check if this seed is in the debug window
        let s_pos = kmer_start + chunk.offset;
        let in_window = if let Some((q_start, q_end, s_start, s_end)) = debug_window {
            q_pos_usize >= q_start && q_pos_usize <= q_end && s_pos >= s_start && s_pos <= s_end
        } else {
            false
        };

        if debug_enabled && in_window {
            *dbg_window_seeds += 1;
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:193-207
        // ```c
        // num_hits = s_BlastLookupGetNumHits(lookup, index);
        // if (num_hits == 0)
        //     continue;
        // s_BlastLookupRetrieve(lookup, index,
        //                       offset_pairs + total_hits, ...);
        // ```
        if blastn_trace_enabled
            && blastn_trace::should_trace_range(
                "seed",
                q_idx,
                s_idx,
                s_id,
                q_pos_usize,
                q_pos_usize.saturating_add(lut_word_length),
                kmer_start.saturating_add(chunk.offset),
                kmer_start
                    .saturating_add(lut_word_length)
                    .saturating_add(chunk.offset),
                ctx.len,
                ctx.frame,
            )
        {
            blastn_trace::log(
                "seed",
                format!(
                    "subject={}({}) context={} q_off={} s_off={} lut_word_length={} lookup_index={} diag={} last_hit_pending",
                    s_id,
                    s_idx,
                    q_idx,
                    q_pos_usize,
                    kmer_start.saturating_add(chunk.offset),
                    lut_word_length,
                    pair.q_off,
                    lookup_diag
                ),
            );
        }

        // NCBI BLAST does NOT check mask_array/mask_hash at seed level
        // NCBI reference: na_ungapped.c:671-672
        // NCBI only checks: if (s_off_pos < last_hit) return 0;
        // This check happens below (line 753) using hit_level_array[real_diag].last_hit
        // The mask_array/mask_hash check above was LOSAT-specific and caused excessive filtering
        // REMOVED to match NCBI BLAST behavior

        // NCBI reference: na_ungapped.c:1081-1140
        // For two-stage lookup, verify word_length match BEFORE ungapped extension
        // This is done by left+right extension from the lut_word_length match position
        let Some((q_ext_start, s_ext_start, word_ext_left, word_ext_right)) = word_extend_pair(
            x,
            LutHit {
                q_off0,
                kmer_start,
                q_context_start,
                q_context_end,
                s_range,
            },
            dbg_left_ext_iters,
            dbg_right_ext_iters,
        ) else {
            continue;
        };

        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1093-1155
        // ```c
        // q_offset -= ext_left;
        // s_offset -= ext_left;
        // ...
        // hits_extended += s_BlastnDiagHashExtendInitialHit(query, subject,
        //                                  q_offset, s_offset,
        //                                  masked_locations,
        //                                  query_info, s_range,
        //                                  word_length, lut_word_length,
        //                                  lookup_wrap,
        // ```
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:828-839
        // ```c
        // diag = s_off - q_off;
        // s_end = s_off + word_length;
        // s_off_pos = s_off + hash_table->offset;
        // rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
        // if(!rc)  last_hit = 0;
        // if (s_off_pos < last_hit) return 0;
        // ```
        // For two-stage lookup NCBI computes the diagonal after the
        // word-length extension has adjusted q_offset/s_offset.
        let diag = if use_array_indexing {
            (s_ext_start as isize + diag_array_length - q_ext_start as isize) & diag_mask
        } else {
            s_ext_start as isize - q_ext_start as isize
        };
        let diag_idx = if use_array_indexing { diag as usize } else { 0 };
        let (last_hit, hit_saved) = if use_array_indexing && diag_idx < diag_array_size {
            let diag_entry = &hit_level_array[diag_idx];
            (diag_entry.last_hit, diag_entry.flag != 0)
        } else if !use_array_indexing {
            // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the same
            // decisions without the walk of a bucket below the bound.
            let (level, _hit_len, hit_saved) = if x_lv_ref {
                diag_hash.retrieve_below::<false>(
                    diag as i32,
                    &x_lv_seeds,
                    s_ext_start + diag_offset as usize,
                )
            } else {
                diag_hash.retrieve(diag as i32)
            }
            .unwrap_or((0, 0, false));
            (level, hit_saved)
        } else {
            (0, false)
        };
        let s_off_pos = s_ext_start + diag_offset as usize;
        let s_off_pos_i32 = s_off_pos as i32;

        // Hit within explored area should be rejected
        if s_off_pos_i32 < last_hit {
            continue;
        }

        // NCBI reference: na_ungapped.c:674-683
        // After word_length verification, call type_of_word for two-hit mode
        // if (two_hits && (hit_saved || s_end_pos > last_hit + window_size)) {
        //     word_type = s_TypeOfWord(...);
        //     if (!word_type) return 0;
        //     s_end += extended;
        //     s_end_pos += extended;
        // }
        let mut q_off = q_ext_start;
        let mut s_off = s_ext_start;
        let mut s_end = s_ext_start + word_length;
        let mut s_end_pos = s_end + diag_offset as usize;
        let mut word_type = 1u8; // Default: single word (when word_length == lut_word_length)
        let mut extended = 0usize;
        let mut off_found = false;
        let mut hit_ready = true;
        let query_mask = ctx.masks.as_slice();

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1093-1155
        // ```c
        // q_offset -= ext_left;
        // s_offset -= ext_left;
        // hits_extended += s_BlastnDiagHashExtendInitialHit(query, subject,
        //                                  q_offset, s_offset,
        //                                  masked_locations, query_info,
        //                                  s_range, word_length, lut_word_length,
        //                                  lookup_wrap, word_params, matrix,
        //                                  ewp->hash_table, init_hitlist,
        //                                  check_masks);
        // ```
        if blastn_trace_enabled
            && blastn_trace::should_trace_range(
                "seed",
                q_idx,
                s_idx,
                s_id,
                q_ext_start.saturating_sub(q_context_start),
                q_ext_start
                    .saturating_add(word_length)
                    .saturating_sub(q_context_start),
                s_ext_start.saturating_add(chunk.offset),
                s_ext_start
                    .saturating_add(word_length)
                    .saturating_add(chunk.offset),
                ctx.len,
                ctx.frame,
            )
        {
            blastn_trace::log(
                "seed",
                format!(
                    "subject={}({}) context={} word_extend raw_seed=({}, {}) adjusted_seed=({}, {}) ext_left={} ext_right={} diag={} last_hit={} hit_saved={} s_off_pos={} s_end_pos={}",
                    s_id,
                    s_idx,
                    q_idx,
                    q_off0.saturating_sub(q_context_start),
                    kmer_start.saturating_add(chunk.offset),
                    q_ext_start.saturating_sub(q_context_start),
                    s_ext_start.saturating_add(chunk.offset),
                    word_ext_left,
                    word_ext_right,
                    diag,
                    last_hit,
                    hit_saved,
                    s_off_pos,
                    s_end_pos
                ),
            );
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:41-69 (s_MBLookup)
        // ```c
        // if (! PV_TEST(pv, index, mb_lt->pv_array_bts)) {
        //     return FALSE;
        // }
        // q_off = mb_lt->hashtable[index];
        // while (q_off) {
        //     if (q_off == q_pos) return TRUE;
        //     q_off = mb_lt->next_pos[q_off];
        // }
        // ```
        let mut is_seed_masked = |s_pos: usize, q_pos: usize| -> bool {
            if s_pos + lut_word_length > s_len {
                return true;
            }
            let kmer = mask_lookup_index(
                packed_kmer_at_seed_mask(search_seq_packed, s_pos, lut_word_length),
                lut_word_length,
            );
            if !two_stage.has_hits(kmer) {
                return true;
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1027-1034
            // ```c
            // /* Also add 1 to all indices, because lookup table indices count
            //    from 1. */
            // mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
            // mb_lt->hashtable[ecode] = index;
            // ```
            let q_off_1 = (q_pos + 1) as u32;
            !two_stage.contains_hit(kmer, q_off_1)
        };

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:674-676
        // ```c
        // if (two_hits && (hit_saved || s_end_pos > last_hit + window_size)) {
        //     word_type = s_TypeOfWord(...);
        // ```
        let s_end_pos_i32 = s_end_pos as i32;
        let window_end = last_hit + window_size as i32;
        if two_hits && (hit_saved || s_end_pos_i32 > window_end) {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:674-680
            // ```c
            // word_type = s_TypeOfWord(query, subject, &q_off, &s_off,
            //                          query_mask, query_info, s_range,
            //                          word_length, lut_word_length, lut, TRUE, &extended);
            // ```
            // s_TypeOfWord uses query_mask to skip masked seeds
            let (wt, ext, q_off_adj, s_off_adj) = type_of_word(
                q_off,
                s_off,
                small_na_word.is_some() || !query_mask.is_empty(),
                q_context_end,
                s_range,
                word_length,
                lut_word_length,
                true, // check_double = TRUE
                &mut is_seed_masked,
            );
            word_type = wt;
            extended = ext;
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:674-680
            // ```c
            // word_type = s_TypeOfWord(query, subject, &q_off, &s_off,
            //                          query_mask, query_info, s_range,
            //                          word_length, lut_word_length, lut, TRUE, &extended);
            // ```
            q_off = q_off_adj;
            s_off = s_off_adj;

            // NCBI: if (!word_type) return 0;
            if word_type == 0 {
                // Non-word, skip this hit
                continue;
            }

            // NCBI: s_end += extended;
            // NCBI: s_end_pos += extended;
            s_end += extended;
            s_end_pos += extended;

            // NCBI reference: na_ungapped.c:685-717
            // for single word, also try off diagonals
            // if (word_type == 1) {
            //     /* try off-diagonals */
            //     Int4 orig_diag = real_diag + diag_table->diag_array_length;
            //     Int4 s_a = s_off_pos + word_length - window_size;
            //     Int4 s_b = s_end_pos - 2 * word_length;
            //     Int4 delta;
            //     if (Delta < 0) Delta = 0;
            //     for (delta = 1; delta <= Delta ; ++delta) {
            //         // Check diag + delta and diag - delta
            //         ...
            //     }
            //     if (!off_found) {
            //         hit_ready = 0;
            //     }
            // }
            if word_type == 1 {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:686-687
                // ```c
                //         if (word_type == 1) {
                //             /* try off-diagonals */
                // ```
                off_found = off_diagonal_found(OffDiagonal {
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
                });

                // NCBI: if (!off_found) {
                //     hit_ready = 0;
                // }
                if !off_found {
                    // NCBI reference: na_ungapped.c:713-716
                    // if (!off_found) {
                    //     /* This is a new hit */
                    //     hit_ready = 0;
                    // }
                    hit_ready = false;
                }
            }
        } else {
            // NCBI reference: na_ungapped.c:718-726
            // else if (check_masks) {
            //     /* check the masks for the word */
            //     if(!s_TypeOfWord(query, subject, &q_off, &s_off,
            //                     query_mask, query_info, s_range,
            //                     word_length, lut_word_length, lut, FALSE, &extended)) return 0;
            //     /* update the right end*/
            //     s_end += extended;
            //     s_end_pos += extended;
            // }
            // In NCBI, check_masks is TRUE by default (only FALSE when lut->stride is true)
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:718-725
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:718-725
            // ```c
            // if(!s_TypeOfWord(query, subject, &q_off, &s_off,
            //                 query_mask, query_info, s_range,
            //                 word_length, lut_word_length, lut, FALSE, &extended)) return 0;
            // ```
            let (wt, ext, q_off_adj, s_off_adj) = type_of_word(
                q_off,
                s_off,
                small_na_word.is_some() || !query_mask.is_empty(),
                q_context_end,
                s_range,
                word_length,
                lut_word_length,
                false, // check_double = FALSE (not in two-hit block)
                &mut is_seed_masked,
            );
            if wt == 0 {
                // Non-word, skip this hit
                continue;
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:718-725
            // ```c
            // if(!s_TypeOfWord(query, subject, &q_off, &s_off,
            //                 query_mask, query_info, s_range,
            //                 word_length, lut_word_length, lut, FALSE, &extended)) return 0;
            // ```
            q_off = q_off_adj;
            s_off = s_off_adj;
            // NCBI: s_end += extended;
            // NCBI: s_end_pos += extended;
            s_end += ext;
            s_end_pos += ext;
            // hit_ready remains true (default) - extension will proceed
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:841-921
        // ```c
        // if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
        //     word_type = s_TypeOfWord(..., TRUE, &extended);
        //     ...
        //     if (word_type == 1) {
        //         ...
        //         if (!off_found) {
        //             hit_ready = 0;
        //         }
        //     }
        // } else if (check_masks) {
        //     if (!s_TypeOfWord(..., FALSE, &extended)) return 0;
        // }
        // ```
        if blastn_trace_enabled
            && blastn_trace::should_trace_range(
                "seed",
                q_idx,
                s_idx,
                s_id,
                q_off.saturating_sub(q_context_start),
                q_off
                    .saturating_add(word_length)
                    .saturating_add(extended)
                    .saturating_sub(q_context_start),
                s_off.saturating_add(chunk.offset),
                s_end_pos
                    .saturating_sub(diag_offset as usize)
                    .saturating_add(chunk.offset),
                ctx.len,
                ctx.frame,
            )
        {
            blastn_trace::log(
                "seed",
                format!(
                    "subject={}({}) context={} type_of_word seed=({}, {}) word_type={} extended={} two_hits={} hit_ready={} off_found={} window_end={} s_end_pos={}",
                    s_id,
                    s_idx,
                    q_idx,
                    q_off.saturating_sub(q_context_start),
                    s_off.saturating_add(chunk.offset),
                    word_type,
                    extended,
                    two_hits,
                    hit_ready,
                    off_found,
                    window_end,
                    s_end_pos
                ),
            );
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:728-766
        // ```c
        // if (hit_ready) {
        //     if (word_params->ungapped_extension) {
        //         ...
        //     }
        // }
        // ```
        if !hit_ready {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:768-772
            // ```c
            // hit_level_array[real_diag].last_hit = s_end_pos;
            // hit_level_array[real_diag].flag = hit_ready;
            // if (two_hits) {
            //     diag_table->hit_len_array[real_diag] =
            //         (hit_ready) ? 0 : s_end_pos - s_off_pos;
            // }
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:939-941
            // ```c
            // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
            //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
            //                       hit_ready, s_off_pos, window_size + Delta + 1);
            // ```
            if use_array_indexing && diag_idx < diag_array_size {
                hit_level_array[diag_idx].last_hit = s_end_pos as i32;
                hit_level_array[diag_idx].flag = 0;
                if two_hits {
                    hit_len_array[diag_idx] = (s_end_pos - s_off_pos) as u8;
                }
            } else if !use_array_indexing {
                diag_hash.insert_plain(
                    diag as i32,
                    s_end_pos as i32,
                    (s_end_pos - s_off_pos) as i32,
                    false,
                    s_off_pos as i32,
                    diag_hash_window,
                );
            }
            continue;
        }

        // Now do ungapped extension from the adjusted position
        // NCBI BLAST: s_BlastnExtendInitialHit calls ungapped extension after word_length verification
        // Use q_off and s_off (adjusted by type_of_word if called)
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:740-749
        // ```c
        // if (word_params->matrix_only_scoring || word_length < 11)
        //    s_NuclUngappedExtendExact(..., -(cutoffs->x_dropoff), ...);
        // else
        //    s_NuclUngappedExtend(..., s_end, s_off, -(cutoffs->x_dropoff),
        //                         word_params->nucl_score_table,
        //                         cutoffs->reduced_nucl_cutoff_score);
        // ```
        if debug_enabled {
            *dbg_ungapped_ext_calls += 1;
        }
        let x_dropoff = x_dropoff_scores[context_idx];
        let reduced_cutoff = reduced_cutoff_scores[context_idx];
        let ungapped_start = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let ungapped = if word_length < 11 {
            extend_hit_ungapped_exact_ncbi(
                encoded_query_concat_blastna,
                search_seq_packed,
                q_off,
                s_off,
                s_len,
                x_dropoff,
                score_matrix,
            )
        } else {
            extend_hit_ungapped_approx_ncbi(
                encoded_query_concat_blastna,
                search_seq_packed,
                q_off,
                s_off,
                s_end,
                s_len,
                x_dropoff,
                nucl_score_table,
                reduced_cutoff,
                score_matrix,
            )
        };
        if let Some(ungapped_start) = ungapped_start {
            let elapsed_ns = ungapped_start.elapsed().as_nanos() as u64;
            if let Some(timing) = timing_ref {
                timing
                    .ungapped_ns
                    .fetch_add(elapsed_ns, std::sync::atomic::Ordering::Relaxed);
                timing
                    .ungapped_calls
                    .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            if track_scan_pause {
                *scan_pause_ns = scan_pause_ns.saturating_add(elapsed_ns);
            }
        }
        let ss = ungapped.s_start;
        let ungapped_se = ungapped.s_start + ungapped.length;
        let ungapped_score = ungapped.score;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:740-758
        // ```c
        // s_NuclUngappedExtendExact(..., -(cutoffs->x_dropoff), ungapped_data);
        // ...
        // if (off_found || ungapped_data->score >= cutoffs->cutoff_score) {
        //     BLAST_SaveInitialHit(init_hitlist, q_off, s_off, final_data);
        // }
        // ```
        if blastn_trace_enabled {
            // No NCBI counterpart (LOSAT_TRACE_BLASTN_* diagnostics): a rejected approximate
            // result may start left of its context (see the NCBI reference at the assertion
            // below), so the trace uses signed context offsets.
            let trace_qs = ungapped.q_start as isize - q_context_start as isize;
            let trace_qe = trace_qs + ungapped.length as isize;
            if blastn_trace::should_trace_range(
                "ungapped",
                q_idx,
                s_idx,
                s_id,
                trace_qs.max(0) as usize,
                trace_qe.max(0) as usize,
                ss.saturating_add(chunk.offset),
                ungapped_se.saturating_add(chunk.offset),
                ctx.len,
                ctx.frame,
            ) {
                blastn_trace::log(
                    "ungapped",
                    format!(
                        "subject={}({}) context={} seed=({}, {}) ungapped=q{}..{} s{}..{} raw_score={} cutoff={} x_dropoff={} off_found={} accepted={}",
                        s_id,
                        s_idx,
                        q_idx,
                        q_off.saturating_sub(q_context_start),
                        s_off.saturating_add(chunk.offset),
                        trace_qs,
                        trace_qe,
                        ss.saturating_add(chunk.offset),
                        ungapped_se.saturating_add(chunk.offset),
                        ungapped_score,
                        cutoff_score,
                        x_dropoff,
                        off_found,
                        off_found || ungapped_score >= cutoff_score
                    ),
                );
            }
        }

        // NCBI reference: na_ungapped.c:757-758
        // s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;
        // This is the END of the UNGAPPED extension, used for last_hit update
        let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

        // Skip if ungapped score is too low
        // NCBI reference: na_ungapped.c:752
        // if (off_found || ungapped_data->score >= cutoffs->cutoff_score)
        // Use dynamically calculated cutoff_score instead of fixed threshold
        // off_found is set by off-diagonal search (Step 3)
        // NCBI: if (off_found || ungapped_data->score >= cutoffs->cutoff_score)
        if !(off_found || ungapped_score >= cutoff_score) {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:752-760
            // ```c
            // if (off_found || ungapped_data->score >= cutoffs->cutoff_score) {
            //     ...
            // } else {
            //     hit_ready = 0;
            // }
            // ```
            hit_ready = false;
            if debug_enabled {
                *dbg_ungapped_low += 1;
            }
            if in_window && debug_mode {
                eprintln!(
                    "[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}",
                    q_pos_usize, s_pos, ungapped_score, min_ungapped_score
                );
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:768-772
            // ```c
            // hit_level_array[real_diag].last_hit = s_end_pos;
            // hit_level_array[real_diag].flag = hit_ready;
            // if (two_hits) {
            //     diag_table->hit_len_array[real_diag] =
            //         (hit_ready) ? 0 : s_end_pos - s_off_pos;
            // }
            // ```
            if use_array_indexing && diag_idx < diag_array_size {
                hit_level_array[diag_idx].last_hit = s_end_pos as i32;
                hit_level_array[diag_idx].flag = if hit_ready { 1 } else { 0 };
                if two_hits {
                    hit_len_array[diag_idx] = if hit_ready {
                        0
                    } else {
                        (s_end_pos - s_off_pos) as u8
                    };
                }
            } else if !use_array_indexing {
                diag_hash.insert_plain(
                    diag as i32,
                    s_end_pos as i32,
                    if hit_ready {
                        0
                    } else {
                        (s_end_pos - s_off_pos) as i32
                    },
                    hit_ready,
                    s_off_pos as i32,
                    diag_hash_window,
                );
            }
            continue;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2384-2392
        // ```c
        //     init_hsp->offsets.qs_offsets.q_off -= query_start;
        //     if (init_hsp->ungapped_data) {
        //         init_hsp->ungapped_data->q_start -= query_start;
        //     }
        //     ASSERT(init_hsp->ungapped_data == NULL ||
        //            init_hsp->ungapped_data->q_start >= 0);
        // ```
        // NCBI checks the start against the context only for saved hits. A rejected
        // approximate result may start left of its context: the left loop of
        // s_NuclUngappedExtend runs over the concatenated query and reads the sentinel as
        // bases (na_ungapped.c:284-305, `len = MIN(q_ext, s_ext) / COMPRESSION_RATIO;`),
        // and that result is kept only below reduced_nucl_cutoff_score = 0.8 * cutoff_score
        // (blast_parameters.c:412), so it is rejected above and never reaches this point.
        debug_assert!(ungapped.q_start >= q_context_start);
        let qs = ungapped.q_start - q_context_start;
        let qe = qs + ungapped.length;

        if debug_enabled {
            *dbg_gapped_attempted += 1;
        }

        if in_window && debug_mode {
            eprintln!(
                "[DEBUG WINDOW] Seed at q={}, s={} -> COLLECT UNGAPPED (score={}, len={})",
                q_pos_usize,
                s_pos,
                ungapped_score,
                qe - qs
            );
        }

        // NCBI architecture: Collect ungapped hit for batch processing
        // NCBI reference: blast_gapalign.c:3824 - init_hsp_array is sorted by score DESCENDING
        // Gapped extension will be done later in score order with containment check
        x_push_ungapped(
            ungapped_hits,
            UngappedHit {
                context_idx: q_idx,
                query_idx: ctx.query_idx,
                query_frame: ctx.frame,
                query_context_offset: ctx.query_offset,
                // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:135-149
                // ```c
                // struct { Uint4 q_off; Uint4 s_off; } qs_offsets;
                // ```
                seed_q_off: (q_off - q_context_start) as u32,
                seed_s_off: s_off as u32,
                qs: qs as u32,
                qe: qe as u32,
                ss: ss as u32,
                se: ungapped_se as u32,
                score: ungapped_score,
            },
        );

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:768-772
        // ```c
        // hit_level_array[real_diag].last_hit = s_end_pos;
        // hit_level_array[real_diag].flag = hit_ready;
        // if (two_hits) {
        //     diag_table->hit_len_array[real_diag] =
        //         (hit_ready) ? 0 : s_end_pos - s_off_pos;
        // }
        // ```
        // Update hit_level_array after collecting ungapped hit
        // (hit_ready is 1 here: flag 1, hit length 0)
        if use_array_indexing && diag_idx < diag_array_size {
            hit_level_array[diag_idx].last_hit = ungapped_s_end_pos as i32;
            hit_level_array[diag_idx].flag = 1;
            if two_hits {
                hit_len_array[diag_idx] = 0;
            }
        } else if !use_array_indexing {
            diag_hash.insert_plain(
                diag as i32,
                ungapped_s_end_pos as i32,
                0,
                true,
                s_off_pos as i32,
                diag_hash_window,
            );
        }
        // Gapped extension deferred to batch processing phase
    }
    *offset_pairs_len = 0;
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
/// Preconditions, which the caller (`extend_offset_pairs`) checks before the call: the hit is
/// dropped unless `q_off0 + lut_word_length <= q_context_end` and `kmer_start + lut_word_length <=
/// s_len`, and `s_range` is at most the subject's length. The unchecked reads of the packed
/// subject and of the sentinel query buffer below rely on them.
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
/// `#[inline(always)]`: `extend_offset_pairs` calls it for every lookup table hit; out of line it
/// took 35 instructions a call (callgrind, n1: dc-megablast with `LOSAT_X_BNDCONE=0` +1.4 %, with
/// the off-diagonal search's call), so it keeps the shape it had inside the scan's closure.
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
                // SAFETY: s_idx tracks s_off / COMPRESSION_RATIO within bounds: `s_off + (ext_to -
                // ext_left) <= s_range` (tested above) and `s_range` is at most the subject's length
                // (the fn's preconditions).
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
            // Reference: na_ungapped.c:1139-1140
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
/// coordinates. Built at each call, from the names of the caller (`na_scan::scan_na_lookup`
/// passes `safe_k` as the word length).
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
/// took 38 instructions a call with `scan_range` 0 (callgrind, n1 dc-megablast with
/// `LOSAT_X_BNDCONE=0`), so it keeps the shape it had inside the scan's closures.
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:657 (the diagonal table's
    // `s_BlastnDiagTableExtendInitialHit`), :824 (the hash's `s_BlastnDiagHashExtendInitialHit`)
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
