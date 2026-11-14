//! The scan of a blastn subject chunk with a lookup table that is not two-stage (the direct
//! table `PvDirectLookup`, or the na table `NaLookupTable`), whose lookup width is the word size
//! (word sizes 4 to 8, `coordination.rs` `choose_na_lookup_table`), as `search_query_batch`'s
//! `collect_prelim_hits_for_chunk` runs it when LOSAT_X_PAIRPARNA (`x_pair_seed_na.rs`) does not
//! take the chunk, on one thread (`scan_na_lookup`): the subject's words in scan order
//! (`scan_subject_kmers_with_ranges`), and for each lookup hit of a word, in the table's order,
//! the seed's context, its diagonal and the `last_hit` test, the word type, in two-hit mode the
//! off-diagonal search of a single word (`pair_extend::off_diagonal_found`, which the steps of
//! the offset pairs of a two-stage table call too), the ungapped extension, the cutoff test and
//! the save of the initial hit. A seed goes through these steps as the scan finds it (no batch
//! of offset pairs) and has no word extension (the word is the lookup word).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684 (the loop of
//! `BlastNaWordFinder`, :1597-1695), where `scansub` is a scanner of
//! c++/src/algo/blast/core/blast_nascan.c and `extend` is `s_BlastNaExtendDirect` (:963-1030),
//! which `BlastChooseNaExtend` picks when the lookup width is the word length (:1811-1812 for the
//! small table, :1828-1829 for the na table)
//! ```c
//!     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
//!
//!         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
//!
//!         if (hitsfound == 0)
//!             continue;
//!
//!         total_hits += hitsfound;
//!         hits_extended += extend(offset_pairs, hitsfound, word_params,
//!                                 lookup_wrap, query, subject, matrix,
//!                                 query_info, ewp, init_hitlist, scan_range[2] + lut_word_length);
//!     }
//! ```
//!
//! This module is a child of `run` for the items it keeps private (the types `DiagHashTable`,
//! `SubjectChunk`, `QueryContext`, `QueryContextIndex`, `UngappedHit` and `BlastnTiming`;
//! `scan_subject_kmers_with_ranges`, `select_mb_scan_kind`, `diag_hash_insert_window`,
//! `x_push_ungapped`, `mask_lookup_index` and `packed_kmer_at_seed_mask`). Its function is the
//! statements of the branch of `collect_prelim_hits_for_chunk` that were there, in their order;
//! the fold of LOSAT_X_PAIRPARNA (`x_pair_seed_na.rs` `fold_na_t`) copies those of one-hit mode
//! with the diagonal hash one by one, so their order is kept here. What the scan reads that is
//! the same for every seed of a chunk is `NaScanInputs`, and what it writes is `NaScanState`,
//! both built once per chunk. Nothing here runs on the pool.

use super::pair_extend;
use super::x_diag_lv;
use super::{
    blastn_trace, diag_hash_insert_window, extend_hit_ungapped_approx_ncbi,
    extend_hit_ungapped_exact_ncbi, mask_lookup_index, packed_kmer_at_seed_mask,
    scan_subject_kmers_with_ranges, select_mb_scan_kind, type_of_word, x_push_ungapped,
    BlastnTiming, DiagHashTable, DiagStruct, QueryContext, QueryContextIndex, SubjectChunk,
    UngappedHit,
};
use crate::algorithm::blastn::alignment::gapped::BlastnaMatrix;
use crate::algorithm::blastn::lookup::{NaLookupTable, PvDirectLookup};

/// What the scan of one subject chunk with a table that is not two-stage reads that is the same
/// for every seed, in groups below. Built once per chunk by `collect_prelim_hits_for_chunk`, after
/// LOSAT_X_PAIRPARNA declines the chunk, with the names there; every field is a shared reference
/// or a `Copy` value.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:964-969,985-992
/// (`s_BlastNaExtendDirect`) and :651,656,658 (`s_BlastnDiagTableExtendInitialHit`; the hash's
/// function has the same lines at :819,823,825)
/// ```c
/// s_BlastNaExtendDirect(const BlastOffsetPair * offset_pairs, Int4 num_hits,
///                       const BlastInitialWordParameters * word_params,
///                       LookupTableWrap * lookup_wrap,
///                       BLAST_SequenceBlk * query,
///                       BLAST_SequenceBlk * subject, Int4 ** matrix,
///                       BlastQueryInfo * query_info,
/// ...
///     else if (lookup_wrap->lut_type == eSmallNaLookupTable) {
///         BlastSmallNaLookupTable *lut =
///                         (BlastSmallNaLookupTable *) lookup_wrap->lut;
///         word_length = lut->word_length;
///     }
///     else {
///         BlastNaLookupTable *lut = (BlastNaLookupTable *) lookup_wrap->lut;
///         word_length = lut->word_length;
/// ...
///     Int4 window_size = word_params->options->window_size;
/// ...
///     Boolean two_hits = (window_size > 0);
/// ...
///     Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
/// ```
/// The level bound and the diagnostics have no NCBI counterpart.
pub(super) struct NaScanInputs<'a> {
    /// The scan: the packed subject, the chunk's length, the word length (the lookup word length
    /// too: these tables have no word extension), the scan step, the chunk's ranges, and whether
    /// the subject is masked.
    pub(super) search_seq_packed: &'a [u8],
    pub(super) s_len: usize,
    pub(super) safe_k: usize,
    pub(super) scan_step: usize,
    pub(super) subject_seq_ranges: &'a [(i32, i32)],
    pub(super) subject_masked: bool,
    /// The lookup table: whether it is the direct table (word sizes up to 13), else the na table.
    pub(super) use_direct_lookup: bool,
    pub(super) pv_direct_lookup_ref: Option<&'a PvDirectLookup>,
    pub(super) na_lookup_ref: Option<&'a NaLookupTable>,
    /// The seed's context: the context of a query offset, the contexts, each context's query, the
    /// contexts' cutoffs, and the chunk (its offset in the subject).
    pub(super) query_context_index: &'a QueryContextIndex,
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) encoded_queries_blastna: &'a [&'a [u8]],
    pub(super) cutoff_scores: &'a [i32],
    pub(super) chunk: &'a SubjectChunk,
    /// The diagonal table: whether it is the array (else the hash), the array's length, mask and
    /// size, the table's offset for this chunk, and whether the hash's retrieve takes the level
    /// bound of LOSAT_X_DIAGLV.
    pub(super) use_array_indexing: bool,
    pub(super) diag_table_length: usize,
    pub(super) diag_table_mask: isize,
    pub(super) diag_array_size: usize,
    pub(super) diag_offset: isize,
    pub(super) x_lv_ref: bool,
    /// The window (0: one hit) and the off-diagonal range.
    pub(super) window_size: usize,
    pub(super) scan_range: usize,
    /// The ungapped extension: the X-dropoff and the reduced cutoff of each context, the score
    /// matrix, and the score table of the approximate extension.
    pub(super) x_dropoff_scores: &'a [i32],
    pub(super) reduced_cutoff_scores: &'a [i32],
    pub(super) score_matrix: &'a BlastnaMatrix,
    pub(super) nucl_score_table: &'a [i32; 256],
    /// The diagnostics: the switch of the debug counters, the debug window, the trace lines
    /// (`blastn_trace`) with the subject's index and trace id, the debug window's lines with the
    /// minimum ungapped score they print, and LOSAT_TIMING.
    pub(super) debug_enabled: bool,
    pub(super) debug_window: Option<(usize, usize, usize, usize)>,
    pub(super) blastn_trace_enabled: bool,
    pub(super) s_idx: usize,
    pub(super) s_id: &'a str,
    pub(super) debug_mode: bool,
    pub(super) min_ungapped_score: i32,
    pub(super) timing_enabled: bool,
    pub(super) timing_ref: Option<&'a BlastnTiming>,
}

/// What the scan writes: the chunk's initial hits, the diagonal table (the array's levels and
/// lengths, or the hash), and the debug counters (the subject positions scanned, the lookup hits,
/// the seeds in the debug window, the hits under the cutoff and the hits saved). Built once per
/// chunk by `collect_prelim_hits_for_chunk` with `NaScanInputs`, from its scratch and its counters
/// (references only, the types the scan had as a branch of that closure); its borrows end with
/// the call.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:970-971
/// (`s_BlastNaExtendDirect`), :641-642 (`s_BlastnDiagTableExtendInitialHit`) and :809
/// (`s_BlastnDiagHashExtendInitialHit`)
/// ```c
///                       Blast_ExtendWord * ewp,
///                       BlastInitHitList * init_hitlist,
/// ...
///                              BLAST_DiagTable * diag_table,
///                              BlastInitHitList * init_hitlist,
/// ...
///                                BLAST_DiagHash * hash_table,
/// ```
/// The counters have no NCBI counterpart (diagnostics only).
pub(super) struct NaScanState<'s> {
    /// The chunk's initial hits (NCBI's `init_hitlist`).
    pub(super) ungapped_hits: &'s mut Vec<UngappedHit>,
    /// The diagonal array (a query block of at most 8000 letters) or the diagonal hash.
    pub(super) hit_level_array: &'s mut Vec<DiagStruct>,
    pub(super) hit_len_array: &'s mut Vec<u8>,
    pub(super) diag_hash: &'s mut DiagHashTable,
    /// The debug counters (written only when `debug_enabled`).
    pub(super) dbg_total_s_positions: &'s mut usize,
    pub(super) dbg_seeds_found: &'s mut usize,
    pub(super) dbg_window_seeds: &'s mut usize,
    pub(super) dbg_ungapped_low: &'s mut usize,
    pub(super) dbg_gapped_attempted: &'s mut usize,
}

/// The scan of one subject chunk with a table that is not two-stage and the steps of each seed
/// it finds, in the scan's order: the subject's words (`scan_subject_kmers_with_ranges` with the
/// word length as the lookup word length), the lookup hits of each word (the direct table, else
/// the na table; `q_off_1 == 0` is skipped), and for each hit its context and cutoff, its
/// diagonal, the diagonal's last hit and the `last_hit` test, the word type (`type_of_word`: in
/// two-hit mode with `check_double`, and for a single word the off-diagonal search,
/// `pair_extend::off_diagonal_found`), the diagonal update of a hit that is not ready, the
/// ungapped extension, the cutoff test, and the save of the initial hit with its diagonal update.
/// The time of the ungapped extensions is taken off the scan's time (`scan_pause_ns`), which
/// LOSAT_TIMING adds up at the end.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:995-1011
/// (`s_BlastNaExtendDirect`, :963-1030; the diagonal array's loop, :1012-1028, calls
/// `s_BlastnDiagTableExtendInitialHit`, :631-775, the same way), which hands each pair, with the
/// word length as the lookup word length, to `s_BlastnDiagHashExtendInitialHit` (:799-944)
/// ```c
///     if (word_params->container_type == eDiagHash) {
///         for (; index < num_hits; ++index) {
///             Int4 s_offset = offset_pairs[index].qs_offsets.s_off;
///             Int4 q_offset = offset_pairs[index].qs_offsets.q_off;
///
///             hits_extended += s_BlastnDiagHashExtendInitialHit(query, subject,
///                                                 q_offset, s_offset,
///                                                 NULL,
///                                                 query_info, s_range,
///                                                 word_length, word_length,
///                                                 lookup_wrap,
///                                                 word_params, matrix,
///                                                 ewp->hash_table,
///                                                 init_hitlist,
///                                                 check_masks);
///         }
///     }
/// ```
/// No `#[inline]`: one call site, the branch of `collect_prelim_hits_for_chunk` it was. Called out
/// of line, it moves the instructions of no measured case by 0.1 % (callgrind, n1; SW7 bn-10, whose
/// move of the two-stage scan ended the register spill this was `#[inline(always)]` for). The
/// scan's callback stays a function of its own, as it was.
pub(super) fn scan_na_lookup(x: &NaScanInputs<'_>, state: &mut NaScanState<'_>) {
    let &NaScanInputs {
        search_seq_packed,
        s_len,
        safe_k,
        scan_step,
        subject_seq_ranges,
        subject_masked,
        use_direct_lookup,
        pv_direct_lookup_ref,
        na_lookup_ref,
        query_context_index,
        query_contexts,
        encoded_queries_blastna,
        cutoff_scores,
        chunk,
        use_array_indexing,
        diag_table_length,
        diag_table_mask,
        diag_array_size,
        diag_offset,
        x_lv_ref,
        window_size,
        scan_range,
        x_dropoff_scores,
        reduced_cutoff_scores,
        score_matrix,
        nucl_score_table,
        debug_enabled,
        debug_window,
        blastn_trace_enabled,
        s_idx,
        s_id,
        debug_mode,
        min_ungapped_score,
        timing_enabled,
        timing_ref,
    } = x;
    let ungapped_hits = &mut *state.ungapped_hits;
    let hit_level_array = &mut *state.hit_level_array;
    let hit_len_array = &mut *state.hit_len_array;
    let diag_hash = &mut *state.diag_hash;
    let dbg_total_s_positions = &mut *state.dbg_total_s_positions;
    let dbg_seeds_found = &mut *state.dbg_seeds_found;
    let dbg_window_seeds = &mut *state.dbg_window_seeds;
    let dbg_ungapped_low = &mut *state.dbg_ungapped_low;
    let dbg_gapped_attempted = &mut *state.dbg_gapped_attempted;
    let scan_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    let mut scan_pause_ns: u64 = 0;
    let mb_scan_kind = select_mb_scan_kind(safe_k, scan_step, subject_masked);
    // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the level bound's constants
    // of this chunk's seeds (`s_end` adds `safe_k`).
    let x_lv_seeds = x_diag_lv::LvSeeds::new(window_size > 0, safe_k, window_size);
    scan_subject_kmers_with_ranges(
        search_seq_packed,
        s_len,
        safe_k,
        safe_k,
        scan_step,
        subject_seq_ranges,
        subject_masked,
        mb_scan_kind,
        |kmer_start, s_range, current_kmer| {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:193-207
            // ```c
            // for (; s <= s_end; s += scan_step) {
            //     num_hits = s_BlastLookupGetNumHits(lookup, index);
            //     if (num_hits == 0)
            //         continue;
            //     s_BlastLookupRetrieve(lookup,
            //                           index,
            //                           offset_pairs + total_hits,
            //                           ...);
            //     total_hits += num_hits;
            // }
            // ```
            if debug_enabled {
                *dbg_total_s_positions += 1;
            }

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:41-79
            // ```c
            // num_hits = s_BlastLookupGetNumHits(lookup, index);
            // if (num_hits == 0)
            //     continue;
            // s_BlastLookupRetrieve(lookup,
            //                       index,
            //                       offset_pairs + total_hits,
            //                       s_off);
            // ```
            // Phase 2: Use PV-based direct lookup (O(1) with fast PV filtering) for word_size <= 13
            // For word_size > 13, use NCBI NaLookupTable array-backed lookup
            let matches_slice: &[u32] = if use_direct_lookup {
                // Use PV for fast filtering before accessing the lookup table
                pv_direct_lookup_ref
                    .map(|pv_dl| pv_dl.get_hits_checked(current_kmer))
                    .unwrap_or(&[])
            } else {
                // Use NaLookupTable for larger word sizes
                na_lookup_ref
                    .map(|na_lt| na_lt.get_hits_checked(current_kmer))
                    .unwrap_or(&[])
            };

            for &q_off_1 in matches_slice {
                if debug_enabled {
                    *dbg_seeds_found += 1;
                }
                if q_off_1 == 0 {
                    continue;
                }
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1412
                // ```c
                // offset_pairs[i].qs_offsets.q_off = q_off - 1;
                // ```
                let q_off0 = q_off_1 as usize - 1;
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:900-903
                // ```c
                // Int4 context = BSearchContextInfo(q_off, query_info);
                // ```
                let context_idx = query_context_index.context_for_offset(q_off0);
                let ctx = &query_contexts[context_idx];
                let q_idx = context_idx as u32;

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:268-269
                // ```c
                // Uint1 *q_start = query->sequence;
                // ```
                let q_seq_blastna = encoded_queries_blastna[context_idx];
                let q_pos_usize = q_off0 - ctx.query_offset as usize;

                // Use pre-computed cutoff score (computed once per query-subject pair)
                let cutoff_score = cutoff_scores[context_idx];

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:663-664
                // ```c
                // diag = s_off + diag_table->diag_array_length - q_off;
                // real_diag = diag & diag_table->diag_mask;
                // ```
                let (diag_array_length, diag_mask) = if use_array_indexing {
                    (diag_table_length as isize, diag_table_mask)
                } else {
                    (0isize, 0isize)
                };
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:663-664
                // ```c
                // diag = s_off + diag_table->diag_array_length - q_off;
                // real_diag = diag & diag_table->diag_mask;
                // ```
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:828-829
                // ```c
                // diag = s_off - q_off;
                // s_end = s_off + word_length;
                // ```
                // `q_off` is the concatenated BLAST query offset returned by the
                // lookup table, not the per-context query position used for output
                // and sequence indexing. Keeping the diagonal in the global query
                // coordinate space preserves NCBI's two-hit last_hit history across
                // plus/minus blastn query contexts.
                let diag = if use_array_indexing {
                    (kmer_start as isize + diag_array_length - q_off0 as isize) & diag_mask
                } else {
                    kmer_start as isize - q_off0 as isize
                };

                // Check if this seed is in the debug window
                let s_pos = kmer_start + chunk.offset;
                let in_window = if let Some((q_start, q_end, s_start, s_end)) = debug_window {
                    q_pos_usize >= q_start
                        && q_pos_usize <= q_end
                        && s_pos >= s_start
                        && s_pos <= s_end
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
                        q_pos_usize.saturating_add(safe_k),
                        kmer_start.saturating_add(chunk.offset),
                        kmer_start
                            .saturating_add(safe_k)
                            .saturating_add(chunk.offset),
                        ctx.len,
                        ctx.frame,
                    )
                {
                    blastn_trace::log(
                        "seed",
                        format!(
                            "subject={}({}) context={} q_off={} s_off={} word_length={} lookup_index={} diag={} hits_for_lookup={}",
                            s_id,
                            s_idx,
                            q_idx,
                            q_pos_usize,
                            s_pos,
                            safe_k,
                            current_kmer,
                            diag,
                            matches_slice.len()
                        ),
                    );
                }

                // NCBI BLAST does NOT check mask_array/mask_hash at seed level
                // NCBI reference: na_ungapped.c:671-672
                // NCBI only checks: if (s_off_pos < last_hit) return 0;
                // This check happens below using hit_level_array[real_diag].last_hit
                // The mask_array/mask_hash check above was LOSAT-specific and caused excessive filtering
                // REMOVED to match NCBI BLAST behavior

                // NCBI reference: na_ungapped.c:656-672
                // Two-hit filter: check if there's a previous hit within window_size
                // NCBI: Boolean two_hits = (window_size > 0);
                // NCBI: last_hit = hit_level_array[real_diag].last_hit;
                // NCBI: hit_saved = hit_level_array[real_diag].flag;
                // NCBI: if (s_off_pos < last_hit) return 0;  // hit within explored area
                // NCBI: if (two_hits && (hit_saved || s_end_pos > last_hit + window_size)) { ... }
                let diag_idx = if use_array_indexing {
                    diag as usize
                } else {
                    0 // Not used for hash indexing
                };

                // Get current diagonal state
                let (last_hit, hit_saved) = if use_array_indexing && diag_idx < diag_array_size {
                    let diag_entry = &hit_level_array[diag_idx];
                    (diag_entry.last_hit, diag_entry.flag != 0)
                } else if !use_array_indexing {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:833-836
                    // ```c
                    // rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
                    // if(!rc)  last_hit = 0;
                    // ```
                    // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the same
                    // decisions without the walk of a bucket below the bound
                    // (`s_off_pos` below).
                    let (level, _hit_len, hit_saved) = if x_lv_ref {
                        diag_hash.retrieve_below::<false>(
                            diag as i32,
                            &x_lv_seeds,
                            kmer_start + diag_offset as usize,
                        )
                    } else {
                        diag_hash.retrieve(diag as i32)
                    }
                    .unwrap_or((0, 0, false));
                    (level, hit_saved)
                } else {
                    (0, false)
                };

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:656-672
                // ```c
                // last_hit = hit_level_array[real_diag].last_hit;
                // s_off_pos = s_off + diag_table->offset;
                // if (s_off_pos < last_hit) return 0;
                // ```
                // In LOSAT, we use kmer_start directly (0-based), but need to add offset for comparison
                let s_off_pos = kmer_start + diag_offset as usize;
                let s_off_pos_i32 = s_off_pos as i32;

                // Hit within explored area should be rejected
                if s_off_pos_i32 < last_hit {
                    continue;
                }

                // NCBI reference: na_ungapped.c:674-683
                // After two-hit check, call type_of_word for two-hit mode
                // if (two_hits && (hit_saved || s_end_pos > last_hit + window_size)) {
                //     word_type = s_TypeOfWord(...);
                //     if (!word_type) return 0;
                //     s_end += extended;
                //     s_end_pos += extended;
                // }
                let two_hits = window_size > 0;
                let mut q_off = q_pos_usize;
                let mut s_off = kmer_start;
                let mut s_end = kmer_start + safe_k;
                let mut s_end_pos = s_end + diag_offset as usize;
                let mut off_found = false;
                let mut hit_ready = true;
                let query_mask = ctx.masks.as_slice();
                let diag_hash_window = if !use_array_indexing {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:825-826, 939-941
                    // ```c
                    // Delta = MIN(word_params->options->scan_range, window_size - word_length);
                    // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
                    //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
                    //                       hit_ready, s_off_pos, window_size + Delta + 1);
                    // ```
                    diag_hash_insert_window(window_size, scan_range, safe_k)
                } else {
                    0
                };

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:459-489
                // ```c
                // static NCBI_INLINE Boolean s_IsSeedMasked(...)
                // {
                //     ...
                //     return !(((T_Lookup_Callback)(lookup_wrap->lookup_callback))
                //                                      (lookup_wrap, index, q_pos));
                // }
                // ```
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:121-133
                // ```c
                // if (! PV_TEST(pv, index, PV_ARRAY_BTS)) {
                //     return FALSE;
                // }
                // num_hits = lookup->thick_backbone[index].num_used;
                // lookup_pos = (num_hits <= NA_HITS_PER_CELL) ?
                //              lookup->thick_backbone[index].payload.entries :
                //              lookup->overflow + lookup->thick_backbone[index].payload.overflow_cursor;
                // for (i=0; i<num_hits; ++i) {
                //     if (lookup_pos[i] == q_pos) return TRUE;
                // }
                // ```
                #[expect(
                    clippy::manual_contains,
                    reason = "a linear scan of the cell's hits with an early exit, as NCBI's loop (na_ungapped.c:131-133); slice::contains would switch this per-seed check to std's chunked search"
                )]
                let mut is_seed_masked = |s_pos: usize, q_pos: usize| -> bool {
                    if s_pos + safe_k > s_len {
                        return true;
                    }
                    let kmer = mask_lookup_index(
                        packed_kmer_at_seed_mask(search_seq_packed, s_pos, safe_k),
                        safe_k,
                    );
                    let hits: &[u32] = if use_direct_lookup {
                        pv_direct_lookup_ref
                            .map(|pv_dl| pv_dl.get_hits_checked(kmer))
                            .unwrap_or(&[])
                    } else {
                        na_lookup_ref
                            .map(|na_lt| na_lt.get_hits_checked(kmer))
                            .unwrap_or(&[])
                    };
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1027-1034
                    // ```c
                    // /* Also add 1 to all indices, because lookup table indices count
                    //    from 1. */
                    // mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
                    // mb_lt->hashtable[ecode] = index;
                    // ```
                    let q_off_1 = (ctx.query_offset as usize + q_pos + 1) as u32;
                    !hits.iter().any(|&hit_q_off| hit_q_off == q_off_1)
                };

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:674-676
                // ```c
                // if (two_hits && (hit_saved || s_end_pos > last_hit + window_size)) {
                //     word_type = s_TypeOfWord(...);
                // ```
                let s_end_pos_i32 = s_end_pos as i32;
                let window_end = last_hit + window_size as i32;
                if two_hits && (hit_saved || s_end_pos_i32 > window_end) {
                    // NCBI reference: na_ungapped.c:677-680
                    // word_type = s_TypeOfWord(query, subject, &q_off, &s_off,
                    //                          query_mask, query_info, s_range,
                    //                          word_length, lut_word_length, lut, TRUE, &extended);
                    // For non-two-stage lookup, word_length == lut_word_length, so type_of_word returns (1, 0)
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:674-680
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:674-680
                    // ```c
                    // word_type = s_TypeOfWord(query, subject, &q_off, &s_off,
                    //                          query_mask, query_info, s_range,
                    //                          word_length, lut_word_length, lut, TRUE, &extended);
                    // ```
                    let (wt, ext, q_off_adj, s_off_adj) = type_of_word(
                        q_off,
                        s_off,
                        !query_mask.is_empty(),
                        ctx.len,
                        s_range,
                        safe_k, // word_length
                        safe_k, // lut_word_length (same for non-two-stage)
                        true,   // check_double = TRUE
                        &mut is_seed_masked,
                    );
                    let word_type = wt;
                    let extended = ext;
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

                    // NCBI reference: na_ungapped.c:685-717 (the diagonal array's,
                    // `s_BlastnDiagTableExtendInitialHit`) and :852-886 (the hash's,
                    // `s_BlastnDiagHashExtendInitialHit`) - off-diagonal search
                    if word_type == 1 {
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:686-687,853-854
                        // ```c
                        //         if (word_type == 1) {
                        //             /* try off-diagonals */
                        // ...
                        //         if (word_type == 1) {
                        //             /* try off-diagonal */
                        // ```
                        off_found = pair_extend::off_diagonal_found(pair_extend::OffDiagonal {
                            use_array_indexing,
                            window_size,
                            scan_range,
                            word_length: safe_k,
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
                        !query_mask.is_empty(),
                        ctx.len,
                        s_range,
                        safe_k, // word_length
                        safe_k, // lut_word_length (same for non-two-stage)
                        false,  // check_double = FALSE (not in two-hit block)
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
                let x_dropoff = x_dropoff_scores[context_idx];
                let reduced_cutoff = reduced_cutoff_scores[context_idx];
                let ungapped_start = if timing_enabled {
                    Some(std::time::Instant::now())
                } else {
                    None
                };
                let ungapped = if safe_k < 11 {
                    extend_hit_ungapped_exact_ncbi(
                        q_seq_blastna,
                        search_seq_packed,
                        q_off,
                        s_off,
                        s_len,
                        x_dropoff,
                        score_matrix,
                    )
                } else {
                    extend_hit_ungapped_approx_ncbi(
                        q_seq_blastna,
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
                    scan_pause_ns = scan_pause_ns.saturating_add(elapsed_ns);
                }
                let qs = ungapped.q_start;
                let qe = ungapped.q_start + ungapped.length;
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
                if blastn_trace_enabled
                    && blastn_trace::should_trace_range(
                        "ungapped",
                        q_idx,
                        s_idx,
                        s_id,
                        qs,
                        qe,
                        ss.saturating_add(chunk.offset),
                        ungapped_se.saturating_add(chunk.offset),
                        ctx.len,
                        ctx.frame,
                    )
                {
                    blastn_trace::log(
                        "ungapped",
                        format!(
                            "subject={}({}) context={} seed=({}, {}) ungapped=q{}..{} s{}..{} raw_score={} cutoff={} x_dropoff={} off_found={} accepted={}",
                            s_id,
                            s_idx,
                            q_idx,
                            q_off,
                            s_off.saturating_add(chunk.offset),
                            qs,
                            qe,
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

                // NCBI reference: na_ungapped.c:757-758
                // s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;
                // This is the END of the UNGAPPED extension, used for last_hit update
                let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

                // Skip if ungapped score is too low
                // NCBI reference: na_ungapped.c:752
                // if (off_found || ungapped_data->score >= cutoffs->cutoff_score)
                // Use dynamically calculated cutoff_score instead of fixed threshold
                // off_found is set by off-diagonal search
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
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:939-941
                        // ```c
                        // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
                        //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
                        //                       hit_ready, s_off_pos, window_size + Delta + 1);
                        // ```
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
                    if debug_enabled {
                        *dbg_ungapped_low += 1;
                    }
                    if in_window && debug_mode {
                        eprintln!(
                            "[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}",
                            q_pos_usize, s_pos, ungapped_score, min_ungapped_score
                        );
                    }
                    continue;
                }

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
                        seed_q_off: q_off as u32,
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
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:939-941
                    // ```c
                    // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
                    //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
                    //                       hit_ready, s_off_pos, window_size + Delta + 1);
                    // ```
                    diag_hash.insert_plain(
                        diag as i32,
                        ungapped_s_end_pos as i32,
                        0,
                        true,
                        s_off_pos as i32,
                        diag_hash_window,
                    );
                }
            }
        },
    );
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:193-207
    // ```c
    // for (; s <= s_end; s += scan_step) {
    //     num_hits = s_BlastLookupGetNumHits(lookup, index);
    //     if (num_hits == 0)
    //         continue;
    //     s_BlastLookupRetrieve(lookup,
    //                           index,
    //                           offset_pairs + total_hits,
    //                           ...);
    //     total_hits += num_hits;
    // }
    // ```
    if let Some(scan_start) = scan_start {
        let elapsed_ns = scan_start.elapsed().as_nanos() as u64;
        let scan_ns = elapsed_ns.saturating_sub(scan_pause_ns);
        if let Some(timing) = timing_ref {
            timing
                .scan_ns
                .fetch_add(scan_ns, std::sync::atomic::Ordering::Relaxed);
            timing
                .scan_calls
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }
}
