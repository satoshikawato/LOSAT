//! The scan of a blastn subject chunk with a lookup table that is not two-stage (the direct
//! table `PvDirectLookup`, or the na table `NaLookupTable`), whose lookup width is the word size
//! (word sizes 4 to 8, `coordination.rs` `choose_na_lookup_table`), as
//! `prelim_subject::prelim_chunk_hits` runs it when LOSAT_X_PAIRPARNA (`x_pair_seed_na.rs`) does
//! not take the chunk, on one thread (`scan_na_lookup`): the subject's words in scan order
//! (`scan_subject_kmers_with_ranges`), and for each lookup hit of a word, in the table's order,
//! the seed's context, its diagonal and the `last_hit` test, in two-hit mode the
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
//!
//! This module is a child of `run` for the items it keeps private (the types `DiagHashTable`,
//! `SubjectChunk` and `UngappedHit`; `scan_subject_kmers_with_ranges`, `select_mb_scan_kind`,
//! `diag_hash_insert_window`, `x_push_ungapped`, `mask_lookup_index` and
//! `packed_kmer_at_seed_mask`). Its function is the statements of the branch of `prelim_chunk_hits`
//! that were there, in their order; the fold of LOSAT_X_PAIRPARNA (`x_pair_seed_na.rs` `fold_na_t`)
//! copies those of one-hit mode with the diagonal hash one by one, so their order is kept here.
//! What the scan reads is `NaScanInputs` (the batch's word finder, `prelim_subject::ScanBatch` or
//! `special_subject`'s copy, and the chunk's own values), and what it writes is `NaScanState`, both
//! built once per chunk. Nothing here runs on the pool.

use super::pair_extend;
use super::prelim_subject::ScanBatch;
use super::{
    blastn_trace, diag_hash_insert_window, extend_hit_ungapped_approx_ncbi,
    extend_hit_ungapped_exact_ncbi, scan_subject_kmers_with_ranges, select_mb_scan_kind,
    x_push_ungapped, DiagHashTable, DiagStruct, SubjectChunk, UngappedHit,
};

/// What the scan of one subject chunk with a table that is not two-stage reads: the batch's word
/// finder (`ScanBatch` or `special_subject`'s copy: the tables, the word length `safe_k`,
/// the scan step, the query block, the cutoffs, the scoring, the diagonal table's form and sizes,
/// the window and the off-diagonal range, the level bound's constants, the diagnostics), and the
/// chunk's own values: the packed subject from the chunk's start, the chunk's length, ranges and
/// whether the subject is masked, the chunk (its offset in the subject), the diagonal table's
/// offset for this chunk, whether the hash's retrieve takes the level bound of LOSAT_X_DIAGLV, and
/// the subject's index and trace id. Built once per chunk by `prelim_chunk_hits`, after
/// LOSAT_X_PAIRPARNA declines the chunk, with the names there.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:964-969,985-992
/// (`s_BlastNaExtendDirect`) and :651,656,658 (`s_BlastnDiagTableExtendInitialHit`; the hash's
/// function has the same lines at :819,823,825)
/// The level bound and the diagnostics have no NCBI counterpart.
pub(super) struct NaScanInputs<'a> {
    pub(super) scan: &'a ScanBatch<'a>,
    pub(super) search_seq_packed: &'a [u8],
    pub(super) s_len: usize,
    pub(super) subject_seq_ranges: &'a [(i32, i32)],
    pub(super) subject_masked: bool,
    pub(super) chunk: &'a SubjectChunk,
    pub(super) diag_offset: isize,
    pub(super) x_lv_ref: bool,
    pub(super) s_idx: usize,
    pub(super) s_id: &'a str,
}

/// What the scan writes: the chunk's initial hits, the diagonal table (the array's levels and
/// lengths, or the hash), and the debug counters (the subject positions scanned, the lookup hits,
/// the seeds in the debug window, the hits under the cutoff and the hits saved). Built once per
/// chunk by `prelim_chunk_hits` with `NaScanInputs`, from its scratch and its counters
/// (references only, the types the scan had as a branch of that fn); its borrows end with
/// the call.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:970-971
/// (`s_BlastNaExtendDirect`), :641-642 (`s_BlastnDiagTableExtendInitialHit`) and :809
/// (`s_BlastnDiagHashExtendInitialHit`)
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
/// diagonal, the diagonal's last hit and the `last_hit` test, in two-hit mode the off-diagonal
/// search of the single word (`pair_extend::off_diagonal_found`; NCBI's `s_TypeOfWord` returns a
/// single word with nothing extended here, na_ungapped.c:531), the diagonal update of a hit that
/// is not ready, the
/// ungapped extension, the cutoff test, and the save of the initial hit with its diagonal update.
/// The time of the ungapped extensions is taken off the scan's time (`scan_pause_ns`), which
/// LOSAT_TIMING adds up at the end.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:995-1011
/// (`s_BlastNaExtendDirect`, :963-1030; the diagonal array's loop, :1012-1028, calls
/// `s_BlastnDiagTableExtendInitialHit`, :631-775, the same way), which hands each pair, with the
/// word length as the lookup word length, to `s_BlastnDiagHashExtendInitialHit` (:799-944)
/// No `#[inline]`: one call site, the branch of `prelim_chunk_hits` it was. Called out
/// of line, it moves the instructions of no measured case by 0.1 % (callgrind, n1; SW7 bn-10, whose
/// move of the two-stage scan ended the register spill this was `#[inline(always)]` for). The
/// scan's callback stays a function of its own, as it was.
pub(super) fn scan_na_lookup(x: &NaScanInputs<'_>, state: &mut NaScanState<'_>) {
    let &NaScanInputs {
        scan,
        search_seq_packed,
        s_len,
        subject_seq_ranges,
        subject_masked,
        chunk,
        diag_offset,
        x_lv_ref,
        s_idx,
        s_id,
    } = x;
    // The batch's values, read once per chunk (the scan's callback reads them as locals, as it read
    // the fields of `x` before they moved to the batch's word finder).
    let &ScanBatch {
        safe_k,
        scan_step,
        use_direct_lookup,
        pv_direct_lookup_ref,
        na_lookup_ref,
        query_context_index,
        query_contexts,
        encoded_query_concat_blastna,
        cutoff_scores,
        use_array_indexing,
        diag_array_length,
        diag_mask,
        diag_array_size,
        window_size,
        scan_range,
        x_dropoff_scores,
        reduced_cutoff_scores,
        score_matrix,
        nucl_score_table,
        exact_ungapped,
        debug_enabled,
        debug_window,
        blastn_trace_enabled,
        debug_mode,
        min_ungapped_score,
        timing_enabled,
        timing_ref,
        x_lv_seeds,
        ..
    } = scan;
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
    // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the level bound's constants of the seeds
    // (`s_end` adds `safe_k`) are the batch's (`scan.x_lv_seeds`: the word length is `safe_k`
    // without a two-stage table).
    scan_subject_kmers_with_ranges(
        search_seq_packed,
        s_len,
        safe_k,
        safe_k,
        scan_step,
        subject_seq_ranges,
        subject_masked,
        mb_scan_kind,
        |kmer_start, _s_range, current_kmer| {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
            if debug_enabled {
                *dbg_total_s_positions += 1;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:41-79
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
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1412
                let q_off0 = q_off_1 as usize - 1;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:900-903
                let context_idx = query_context_index.context_for_offset(q_off0);
                let ctx = &query_contexts[context_idx];
                let q_idx = context_idx as u32;

                let q_pos_usize = q_off0 - ctx.query_offset as usize;

                // Use pre-computed cutoff score (computed once per query-subject pair)
                let cutoff_score = cutoff_scores[context_idx];

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:663-664
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-829
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

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
                if blastn_trace_enabled
                    && blastn_trace::should_trace_range(
                        "seed",
                        q_idx,
                        s_idx,
                        s_id,
                        q_pos_usize as i32,
                        q_pos_usize.saturating_add(safe_k) as i32,
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
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:671-672
                // NCBI only checks whether the hit is before the diagonal's last hit (`s_off_pos < last_hit`: no extension).
                // This check happens below using hit_level_array[real_diag].last_hit
                // The mask_array/mask_hash check above was LOSAT-specific and caused excessive filtering
                // REMOVED to match NCBI BLAST behavior

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:656-672
                // Two-hit filter: check if there's a previous hit within window_size
                // NCBI: if (s_off_pos < last_hit) return 0;  // hit within explored area
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
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:833-836
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

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:656-672
                // In LOSAT, we use kmer_start directly (0-based), but need to add offset for comparison
                let s_off_pos = kmer_start + diag_offset as usize;
                let s_off_pos_i32 = s_off_pos as i32;

                // Hit within explored area should be rejected
                if s_off_pos_i32 < last_hit {
                    continue;
                }

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:674-726
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:531
                // `s_TypeOfWord` returns a single word with nothing extended when the word length
                // is the lookup word length, as here (both `safe_k`: the table is not two-stage),
                // whether or not the double word is checked: the word keeps the hit's offsets,
                // and the masks and the length of the hit's context are not read. In two-hit mode
                // a single word then gets the off-diagonal search.
                let two_hits = window_size > 0;
                let q_off = q_pos_usize;
                let s_off = kmer_start;
                let s_end = kmer_start + safe_k;
                let s_end_pos = s_end + diag_offset as usize;
                let mut off_found = false;
                let mut hit_ready = true;
                let diag_hash_window = if !use_array_indexing {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825-826, 939-941
                    diag_hash_insert_window(window_size, scan_range, safe_k)
                } else {
                    0
                };

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:674-676
                let s_end_pos_i32 = s_end_pos as i32;
                let window_end = last_hit + window_size as i32;
                if two_hits && (hit_saved || s_end_pos_i32 > window_end) {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:685-717 (the diagonal array's,
                    // `s_BlastnDiagTableExtendInitialHit`) and :852-886 (the hash's,
                    // `s_BlastnDiagHashExtendInitialHit`) - off-diagonal search
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:686-687,853-854
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
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:713-716
                        hit_ready = false;
                    }
                }

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:728-766
                if !hit_ready {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:768-772
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
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

                // Now do ungapped extension from the word's offsets
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:740-749
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:268-269
                // NCBI's `query` is the concatenated query and `q_off` an offset in it
                // (na_ungapped.c:164-206 too): a walk that does not stop at the context's
                // sentinel (a raw X-drop of 2^30 or more) reads the contexts before it. The word
                // keeps the lookup hit's offsets, so `q_off` is the hit's `q_off0`.
                let q_off_concat = q_off0;
                let x_dropoff = x_dropoff_scores[context_idx];
                let reduced_cutoff = reduced_cutoff_scores[context_idx];
                let ungapped_start = if timing_enabled {
                    Some(std::time::Instant::now())
                } else {
                    None
                };
                let ungapped = if exact_ungapped {
                    extend_hit_ungapped_exact_ncbi(
                        encoded_query_concat_blastna,
                        search_seq_packed,
                        q_off_concat,
                        s_off,
                        s_len,
                        x_dropoff,
                        score_matrix,
                    )
                } else {
                    extend_hit_ungapped_approx_ncbi(
                        encoded_query_concat_blastna,
                        search_seq_packed,
                        q_off_concat,
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
                let ss = ungapped.s_start;
                let ungapped_se = ungapped.s_start + ungapped.length;
                let ungapped_score = ungapped.score;

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:740-758
                if blastn_trace_enabled {
                    // No NCBI counterpart (LOSAT_TRACE_BLASTN_* diagnostics): the hit's
                    // context-relative range, NCBI's `Int4` one (below 0 for a hit that starts
                    // left of its context), computed only for the trace.
                    let qs = (ungapped.q_start as isize - ctx.query_offset as isize) as i32;
                    let qe = qs + ungapped.length as i32;
                    if blastn_trace::should_trace_range(
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
                    ) {
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
                }

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:757-758
                // This is the END of the UNGAPPED extension, used for last_hit update
                let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

                // Skip if ungapped score is too low
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:752
                // Use dynamically calculated cutoff_score instead of fixed threshold
                // off_found is set by off-diagonal search
                if !(off_found || ungapped_score >= cutoff_score) {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:752-760
                    hit_ready = false;
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:768-772
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
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
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
                        q_pos_usize, s_pos, ungapped_score, ungapped.length
                    );
                }

                // NCBI architecture: Collect ungapped hit for batch processing
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3824 - init_hsp_array is sorted by score DESCENDING
                // `UngappedHit::saved` holds NCBI's saving of the hit (`BLAST_SaveInitialHit`) and
                // its offsets in the context (`s_AdjustInitialHSPOffsets`, a negative `Int4` start
                // for a hit that starts left of its context).
                x_push_ungapped(
                    ungapped_hits,
                    UngappedHit::saved(q_idx, ctx, q_off_concat, s_off, &ungapped),
                );

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:768-772
                // Update hit_level_array after collecting ungapped hit
                // (hit_ready is 1 here: flag 1, hit length 0)
                if use_array_indexing && diag_idx < diag_array_size {
                    hit_level_array[diag_idx].last_hit = ungapped_s_end_pos as i32;
                    hit_level_array[diag_idx].flag = 1;
                    if two_hits {
                        hit_len_array[diag_idx] = 0;
                    }
                } else if !use_array_indexing {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
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
