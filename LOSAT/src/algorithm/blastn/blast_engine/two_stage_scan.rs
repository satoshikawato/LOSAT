//! The scan of a blastn subject chunk with a two-stage lookup table (`TwoStageLookup`: every table
//! of megablast and dc-megablast, and the tables of `-task blastn` whose lookup width is below the
//! word size), as `search_query_batch`'s `collect_prelim_hits_for_chunk` runs it when PAIRPAR
//! (`x_pair_seed.rs`; for a discontiguous table `x_pair_seed_disc.rs`) does not take the chunk, on
//! one thread (`scan_two_stage_lookup`): the subject's words in scan order, through the
//! discontiguous scan (`scan_subject_disc_words_with_ranges`; with LOSAT_X_DCCSR, a table laid out
//! as CSR through `x_dc_csr::Ring`), the scan of a contiguous table laid out as CSR (LOSAT_X_BNCSR,
//! one scan per form of its head, `x_scan_csr!`) or the contiguous scan
//! (`scan_subject_kmers_with_ranges`). The lookup hits of each word, in the table's order, become
//! offset pairs of a batch, which goes to the steps of the pairs (`pair_extend::extend_offset_pairs`:
//! the word extension, the diagonal, the word type, the ungapped extension and the save of the
//! initial hit) when it is full and once after the scan.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684 (the loop of
//! `BlastNaWordFinder`, :1597-1695, with the word lengths of a megablast table, :1625-1637), where
//! `scansub` is the megablast scanner that `s_MBChooseScanSubject`
//! (c++/src/algo/blast/core/blast_nascan.c:2602-2680) picks and `extend` is the variant
//! `BlastChooseNaExtend` (c++/src/algo/blast/core/na_ungapped.c:1791-1836) picks:
//! `s_BlastNaExtendDirect` for a discontiguous table or one whose width is the word size
//! (:1798-1799), else `s_BlastNaExtendAligned` or `s_BlastNaExtend` (:1052-1174). For the small-query table the
//! scanner is `s_SmallNaChooseScanSubject`'s (blast_nascan.c:1298) and the extension the one
//! `BlastChooseNaExtend` picks at na_ungapped.c:1806-1818 (`s_BlastSmallNaExtendAlignedOneByte` or
//! `s_BlastSmallNaExtend`); its word lengths are at na_ungapped.c:1617-1624
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
//! NCBI's scanner returns when its batch is full and the loop calls it again where it stopped;
//! this scan hands the batch to the steps of the pairs at that word and goes on (the same pairs in
//! the same batches, in the same order).
//!
//! This module is a child of `run` for the items it keeps private (the types `DiagHashTable`,
//! `OffsetPair`, `SubjectChunk`, `QueryContext`, `QueryContextIndex`, `UngappedHit` and
//! `BlastnTiming`; `scan_subject_kmers_with_ranges`, `scan_subject_disc_words_with_ranges`,
//! `select_mb_scan_kind`, `diag_hash_insert_window`, `x_mb_hit_prefetch`, `x_bn_pv_once` and
//! `OFFSET_ARRAY_SIZE`). Its function is the statements of the branch of
//! `collect_prelim_hits_for_chunk` that were there, in their order; the strips scans of PAIRPAR
//! (`x_pair_seed.rs` `scan_strip`, `x_pair_seed_disc.rs` `scan_strip_disc`) walk the lookup hits
//! in the order `on_word` does. What the scan reads that is the same for every word of a chunk is
//! `TwoStageScanInputs`, and what it writes is `TwoStageScanState`, both built once per chunk.
//! Nothing here runs on the pool.

use super::pair_extend;
use super::x_diag_lv;
use super::{
    diag_hash_insert_window, scan_subject_disc_words_with_ranges, scan_subject_kmers_with_ranges,
    select_mb_scan_kind, x_bn_pv_once, x_dc_csr, x_mb_hit_prefetch, BlastnTiming, DiagHashTable,
    DiagStruct, OffsetPair, QueryContext, QueryContextIndex, SmallNaWord, SubjectChunk,
    UngappedHit, OFFSET_ARRAY_SIZE,
};
use crate::algorithm::blastn::alignment::gapped::BlastnaMatrix;
use crate::algorithm::blastn::lookup::TwoStageLookup;

/// What the scan of one subject chunk with a two-stage table reads that is the same for every
/// word, in groups below. Built once per chunk by `collect_prelim_hits_for_chunk`, after PAIRPAR
/// declines the chunk, with the names there; every field is a shared reference or a `Copy` value.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1597-1602,1625-1634
/// (`BlastNaWordFinder`) and c++/src/algo/blast/core/blast_engine.c:451 (the batch's size)
/// ```c
/// Int2 BlastNaWordFinder(BLAST_SequenceBlk * subject,
///                        BLAST_SequenceBlk * query,
///                        BlastQueryInfo * query_info,
///                        LookupTableWrap * lookup_wrap,
///                        Int4 ** matrix,
///                        const BlastInitialWordParameters * word_params,
/// ...
///     else if (lookup_wrap->lut_type == eMBLookupTable) {
///         BlastMBLookupTable *lookup =
///                                 (BlastMBLookupTable *) lookup_wrap->lut;
///         if (lookup->discontiguous) {
///             word_length = lookup->template_length;
///             lut_word_length = lookup->template_length;
///         } else {
///             word_length = lookup->word_length;
///             lut_word_length = lookup->lut_word_length;
/// ...
///     const int kScanSubjectOffsetArraySize = GetOffsetArraySize(lookup);
/// ```
/// The level bound and the diagnostics have no NCBI counterpart.
pub(super) struct TwoStageScanInputs<'a> {
    /// The table and its scan: the two-stage table, the packed subject, the chunk's length, the
    /// scan step, the chunk's ranges, whether the subject is masked, and the size at which a batch
    /// of offset pairs goes to the steps of the pairs (`OFFSET_ARRAY_SIZE` and the table's longest
    /// chain).
    pub(super) two_stage: &'a TwoStageLookup,
    pub(super) search_seq_packed: &'a [u8],
    pub(super) s_len: usize,
    pub(super) scan_step: usize,
    pub(super) subject_seq_ranges: &'a [(i32, i32)],
    pub(super) subject_masked: bool,
    pub(super) offset_array_size: usize,
    /// The word extension of the small-query table: its compressed query (empty for the other
    /// tables), whether its aligned one-byte form applies, and the query with its sentinels.
    pub(super) small_na_compressed_query: &'a [u8],
    pub(super) small_na_aligned_one_byte: bool,
    pub(super) encoded_query_concat_blastna_with_sentinels: &'a [u8],
    /// The pair's context: the context of a query offset, the contexts, their cutoffs, and the
    /// chunk (its offset in the subject).
    pub(super) query_context_index: &'a QueryContextIndex,
    pub(super) query_contexts: &'a [QueryContext],
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
    /// The ungapped extension: the concatenated query, the X-dropoff and the reduced cutoff of each
    /// context, the score matrix, and the score table of the approximate extension.
    pub(super) encoded_query_concat_blastna: &'a [u8],
    pub(super) x_dropoff_scores: &'a [i32],
    pub(super) reduced_cutoff_scores: &'a [i32],
    pub(super) score_matrix: &'a BlastnaMatrix,
    pub(super) nucl_score_table: &'a [i32; 256],
    /// The diagnostics: the switch of the debug counters, BLEMIR_DEBUG and LOSAT_DEBUG_BLASTN (the
    /// warning on a long chain and the `[PERF]` line), the debug window, the trace lines
    /// (`blastn_trace`) with the subject's index and trace id, the debug window's lines with the
    /// minimum ungapped score they print, and LOSAT_TIMING.
    pub(super) debug_enabled: bool,
    pub(super) debug_mode: bool,
    pub(super) blastn_debug: bool,
    pub(super) debug_window: Option<(usize, usize, usize, usize)>,
    pub(super) blastn_trace_enabled: bool,
    pub(super) s_idx: usize,
    pub(super) s_id: &'a str,
    pub(super) min_ungapped_score: i32,
    pub(super) timing_enabled: bool,
    pub(super) timing_ref: Option<&'a BlastnTiming>,
}

/// What the scan writes: the chunk's initial hits, the diagonal table (the array's levels and
/// lengths, or the hash), the batch of offset pairs, and the debug counters (the subject positions
/// scanned, the lookup hits, the seeds in the debug window, the hits under the cutoff and the hits
/// saved). Built once per chunk by `collect_prelim_hits_for_chunk` with `TwoStageScanInputs`, from
/// its scratch and its counters (references only, the types the scan had as a branch of that
/// closure); its borrows end with the call.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1603-1606
/// (`BlastNaWordFinder`)
/// ```c
///                        Blast_ExtendWord * ewp,
///                        BlastOffsetPair * offset_pairs,
///                        Int4 max_hits,
///                        BlastInitHitList * init_hitlist,
/// ```
/// The counters have no NCBI counterpart (diagnostics only).
pub(super) struct TwoStageScanState<'s> {
    /// The chunk's initial hits (NCBI's `init_hitlist`).
    pub(super) ungapped_hits: &'s mut Vec<UngappedHit>,
    /// The diagonal array (a query block of at most 8000 letters) or the diagonal hash.
    pub(super) hit_level_array: &'s mut Vec<DiagStruct>,
    pub(super) hit_len_array: &'s mut Vec<u8>,
    pub(super) diag_hash: &'s mut DiagHashTable,
    /// The batch of offset pairs (NCBI's `offset_pairs`, allocated for `offset_array_size` pairs).
    pub(super) offset_pairs: &'s mut Vec<OffsetPair>,
    /// The debug counters (written only when `debug_enabled`).
    pub(super) dbg_total_s_positions: &'s mut usize,
    pub(super) dbg_seeds_found: &'s mut usize,
    pub(super) dbg_window_seeds: &'s mut usize,
    pub(super) dbg_ungapped_low: &'s mut usize,
    pub(super) dbg_gapped_attempted: &'s mut usize,
}

/// The scan of one subject chunk with a two-stage table and the steps of the offset pairs it
/// collects, in the scan's order: the chunk's constants (the word and lookup word lengths, the
/// small-query table's word, the diagonal mask, the hash's insert window, two hits or one, the
/// level bound's constants), the inputs and the state of the steps of the pairs
/// (`pair_extend::PairExtendInputs`, `PairExtendState`), then the subject's words (`on_word`: the
/// lookup hits of the word, and of the second template's word, appended to the batch; the batch
/// goes to `pair_extend::extend_offset_pairs` before a word when it holds `OFFSET_ARRAY_SIZE` pairs
/// and after it when it holds `offset_array_size`) through the discontiguous, CSR or contiguous
/// scan, the scan's time (LOSAT_TIMING, without the time of the ungapped extensions,
/// `scan_pause_ns`), the last batch, and the debug `[PERF]` line.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1504-1506,1523-1537
/// (`s_MBScanSubject_Any`, :1482-1626, one of the scanners that `s_MBChooseScanSubject`,
/// :2602-2680, picks for the loop of `BlastNaWordFinder`, na_ungapped.c:1673-1684), whose batch
/// ends before a word with hits once it holds `max_hits - longest_chain` pairs; the loop then hands
/// the batch to `extend`, `s_BlastNaExtend` (na_ungapped.c:1052-1174)
/// ```c
///    /* Since the test for number of hits here is done after adding them,
///       subtract the longest chain length from the allowed offset array size. */
///    max_hits -= mb_lt->longest_chain;
/// ...
///       for ( ; s <= s_end; s += scan_step) {
///
///          index = s[0] << 16 | s[1] << 8 | s[2];
///          index = index >> shift;
///
///          if (s_BlastMBLookupHasHits(mb_lt, index)) {
///              if (total_hits >= max_hits)
///                  break;
///              s_off = (Int4)((s - abs_start)*COMPRESSION_RATIO);
///              total_hits += s_BlastMBLookupRetrieve(mb_lt,
///                  index,
///                  offset_pairs + total_hits,
///                  s_off);
///          }
///       }
/// ```
/// `#[inline(always)]`: the branch of `collect_prelim_hits_for_chunk` that this was stays in that
/// closure's code (one copy in each instance of `search_query_batch`). Out of line, its
/// dc-megablast scan loop (LOSAT_X_BNDCONE=0, then inlined in this fn) took another register
/// allocation and block layout: one more instruction and five executed instructions of
/// loop-alignment padding per subject byte (+0.84 % for a 3,000-letter query, callgrind, n1).
/// `on_word` and the other callbacks stay functions of their own, as they were. Their closures no
/// longer sit in the generic over `UNGAPPED`, so they, and the generic scans called with them, have
/// one instance for both instances of `search_query_batch` and run out of line: the discontiguous
/// scan's loop (13.46 M -> 13.39 M instructions for that query) and the LOSAT_X_DCCSR ring's scan
/// (26.67 M -> 27.04 M for a 600 kb x 300 kb pair).
#[inline(always)]
pub(super) fn scan_two_stage_lookup(x: &TwoStageScanInputs<'_>, state: &mut TwoStageScanState<'_>) {
    let &TwoStageScanInputs {
        two_stage,
        search_seq_packed,
        s_len,
        scan_step,
        subject_seq_ranges,
        subject_masked,
        offset_array_size,
        small_na_compressed_query,
        small_na_aligned_one_byte,
        encoded_query_concat_blastna_with_sentinels,
        query_context_index,
        query_contexts,
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
        encoded_query_concat_blastna,
        x_dropoff_scores,
        reduced_cutoff_scores,
        score_matrix,
        nucl_score_table,
        debug_enabled,
        debug_mode,
        blastn_debug,
        debug_window,
        blastn_trace_enabled,
        s_idx,
        s_id,
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
    // For two-stage lookup, use lut_word_length (8) for scanning
    let lut_word_length = two_stage.lut_word_length();
    let word_length = two_stage.word_length();
    // NCBI's small-query table always keeps `masked_locations`, so
    // `s_TypeOfWord` checks the lookup words of every word hit (the
    // compressed query reads ambiguity codes as bases, and the lookup
    // table has no words with ambiguity codes).
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:407-411
    // ```c
    //     if (locations &&
    //         lookup->word_length > lookup->lut_word_length ) {
    //         /* because we use compressed query, we must always check masked location*/
    //         lookup->masked_locations = s_SeqLocListInvert(locations, query->length);
    //     }
    // ```
    #[expect(
        clippy::unnecessary_lazy_evaluations,
        reason = "the word is built for the small-query table only, as before the move (its query was a `Vec` there, whose borrow the lint does not count as cheap)"
    )]
    let small_na_word = (!small_na_compressed_query.is_empty()).then(|| SmallNaWord {
        compressed_query: small_na_compressed_query,
        subject_packed: search_seq_packed,
        query_length: encoded_query_concat_blastna.len(),
        word_length,
        lut_word_length,
    });

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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:825-826, 939-941
    // ```c
    // Delta = MIN(word_params->options->scan_range, window_size - word_length);
    // if (Delta < 0) Delta = 0;
    // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
    //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
    //                       hit_ready, s_off_pos, window_size + Delta + 1);
    // ```
    let diag_hash_window = if !use_array_indexing {
        diag_hash_insert_window(window_size, scan_range, word_length)
    } else {
        0
    };

    let two_hits = window_size > 0;
    // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the level bound's constants
    // of this chunk's seeds.
    let x_lv_seeds = x_diag_lv::LvSeeds::new(two_hits, word_length, window_size);

    // DEBUG: Count loop iterations
    let mut dbg_left_ext_iters = 0usize;
    let mut dbg_right_ext_iters = 0usize;
    let mut dbg_ungapped_ext_calls = 0usize;
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
    let scan_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    let mut scan_pause_ns: u64 = 0;
    let dbg_start_time = if debug_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:451-497
    // ```c
    // const int kScanSubjectOffsetArraySize = GetOffsetArraySize(lookup);
    // hitsfound = scansub(lookup_wrap, subject, offset_pairs,
    //                     kScanSubjectOffsetArraySize, &scan_range[1]);
    // if (hitsfound == 0) continue;
    // hits_extended += extend(offset_pairs, hitsfound, ...);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1482-1534
    // ```c
    // Int4 total_hits = 0;
    // ...
    // if (total_hits >= max_hits)
    //     break;
    // ...
    // total_hits += s_BlastMBLookupRetrieve(mb_lt,
    //     index, offset_pairs + total_hits, s_off);
    // ```
    let offset_pairs = state.offset_pairs.as_mut_slice();
    let mut offset_pairs_len = 0usize;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:451-497
    // ```c
    // const int kScanSubjectOffsetArraySize = GetOffsetArraySize(lookup);
    // hitsfound = scansub(lookup_wrap, subject, offset_pairs,
    //                     kScanSubjectOffsetArraySize, &scan_range[1]);
    // if (hitsfound == 0) continue;
    // hits_extended += extend(offset_pairs, hitsfound, ...);
    // ```
    let x_hit_prefetch = x_mb_hit_prefetch();
    // What the steps of the offset pairs of this chunk
    // (`pair_extend::extend_offset_pairs`) read, built once per chunk.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1052-1060
    // ```c
    // s_BlastNaExtend(const BlastOffsetPair * offset_pairs, Int4 num_hits,
    //                 const BlastInitialWordParameters * word_params,
    //                 LookupTableWrap * lookup_wrap,
    //                 BLAST_SequenceBlk * query,
    //                 BLAST_SequenceBlk * subject, Int4 ** matrix,
    //                 BlastQueryInfo * query_info,
    //                 Blast_ExtendWord * ewp,
    //                 BlastInitHitList * init_hitlist,
    //                 Uint4 s_range)
    // ```
    let pair_inputs = pair_extend::PairExtendInputs {
        small_na_word: &small_na_word,
        small_na_aligned_one_byte,
        subject_masked,
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
    };
    // What the steps of the offset pairs of this chunk write, held by the scan's
    // callbacks for the whole scan.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1058-1059
    // ```c
    //                 Blast_ExtendWord * ewp,
    //                 BlastInitHitList * init_hitlist,
    // ```
    let mut pair_state = pair_extend::PairExtendState {
        ungapped_hits: &mut *ungapped_hits,
        hit_level_array: &mut *hit_level_array,
        hit_len_array: &mut *hit_len_array,
        diag_hash: &mut *diag_hash,
        dbg_window_seeds: &mut *dbg_window_seeds,
        dbg_left_ext_iters: &mut dbg_left_ext_iters,
        dbg_right_ext_iters: &mut dbg_right_ext_iters,
        dbg_ungapped_ext_calls: &mut dbg_ungapped_ext_calls,
        dbg_ungapped_low: &mut *dbg_ungapped_low,
        dbg_gapped_attempted: &mut *dbg_gapped_attempted,
    };

    let mut on_word = |kmer_start: usize,
                       s_range: usize,
                       current_lut_kmer: u64,
                       second_template_index: Option<u64>| {
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
        if debug_mode || blastn_debug {
            // DEBUG: Track processing time
            let _dbg_start = std::time::Instant::now();
        }

        // Lookup using lut_word_length k-mer
        let matches_len = if debug_mode || blastn_debug {
            two_stage.count_hits(current_lut_kmer)
        } else {
            0
        };

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
        // DEBUG: Log max matches
        if (debug_mode || blastn_debug) && matches_len > 1000 {
            eprintln!(
                "[WARN] Large matches_slice: len={} for kmer at position {}",
                matches_len, kmer_start
            );
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1504-1533
        // ```c
        // max_hits -= mb_lt->longest_chain;
        // if (total_hits >= max_hits)
        //     break;
        // total_hits += s_BlastMBLookupRetrieve(mb_lt,
        //     index, offset_pairs + total_hits, s_off);
        // ```
        if offset_pairs_len >= OFFSET_ARRAY_SIZE {
            pair_extend::extend_offset_pairs(
                &pair_inputs,
                &mut pair_state,
                offset_pairs,
                &mut offset_pairs_len,
                &mut scan_pause_ns,
                true,
            );
        }

        // For each match, add to the offset_pairs buffer.
        two_stage.for_each_hit(current_lut_kmer, |q_off_1| {
            if debug_enabled {
                *dbg_seeds_found += 1;
            }
            if q_off_1 != 0 {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1418
                // ```c
                // offset_pairs[i].qs_offsets.q_off   = q_off - 1;
                // offset_pairs[i++].qs_offsets.s_off = s_off;
                // ```
                offset_pairs[offset_pairs_len] = OffsetPair {
                    q_off: q_off_1 as usize - 1,
                    s_off: kmer_start,
                    s_range,
                };
                offset_pairs_len += 1;
            }
        });
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1456-1469
        // ```c
        // /** access the megablast lookup table with two templates */
        // #define MB_ACCESS_HITS2()                               \
        //    if (total_hits >= max_hits)                          \
        //        break;                                           \
        //    if (s_BlastMBLookupHasHits(mb_lt, index)) {          \
        //        total_hits += s_BlastMBLookupRetrieve(mb_lt,     \
        //                   index, offset_pairs + total_hits,     \
        //                   scan_range[0]);                               \
        //    }                                                    \
        //    if (s_BlastMBLookupHasHits(mb_lt, index2)) {         \
        //        total_hits += s_BlastMBLookupRetrieve2(mb_lt,    \
        //                   index2, offset_pairs + total_hits,    \
        //                   scan_range[0]);                               \
        //    }
        // ```
        // The second template's hits of the same subject offset follow
        // the first template's (`longest_chain` counts both chains).
        if let Some(second_template_index) = second_template_index {
            two_stage.for_each_hit2(second_template_index, |q_off_1| {
                if debug_enabled {
                    *dbg_seeds_found += 1;
                }
                offset_pairs[offset_pairs_len] = OffsetPair {
                    q_off: q_off_1 as usize - 1,
                    s_off: kmer_start,
                    s_range,
                };
                offset_pairs_len += 1;
            });
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
        // ```c
        // offset_array_size = OFFSET_ARRAY_SIZE +
        //     ((BlastMBLookupTable*)lookup->lut)->longest_chain;
        // ```
        if offset_pairs_len >= offset_array_size {
            pair_extend::extend_offset_pairs(
                &pair_inputs,
                &mut pair_state,
                offset_pairs,
                &mut offset_pairs_len,
                &mut scan_pause_ns,
                true,
            );
        }
    };
    if let Some(disc) = two_stage.disc() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1647-1684
        // ```c
        //     scan_range[0] = 0;  /* subject seq mask index */
        //     scan_range[1] = 0;	/* start pos of scan */
        //     scan_range[2] = subject->length - lut_word_length; /*end pos (inclusive) of scan*/
        //
        //     /* if sequence is masked, fall back to generic scanner and extender */
        //     if (subject->mask_type != eNoSubjMasking) {
        //         if (lookup_wrap->lut_type == eMBLookupTable &&
        //             ((BlastMBLookupTable *) lookup_wrap->lut)->discontiguous) {
        //             /* discontiguous scan subs assumes any (non-aligned starting offset */
        //         } else {
        //             ...
        //         }
        //         /* generic scanner permits any (non-aligned) starting offset */
        //         scan_range[1] = subject->seq_ranges[0].left + word_length - lut_word_length;
        //         scan_range[2] = subject->seq_ranges[0].right - lut_word_length;
        //     }
        //     ...
        //     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
        //
        //         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
        //         ...
        //         hits_extended += extend(offset_pairs, hitsfound, word_params,
        //                                 lookup_wrap, query, subject, matrix,
        //                                 query_info, ewp, init_hitlist, scan_range[2] + lut_word_length);
        //     }
        // ```
        // EXPERIMENT (LOSAT_X_DCCSR): a table laid out as CSR has no chains; its
        // words go through `x_on_word_csr`, `on_word`'s statements for one template
        // with the cell's run (the chain's offsets in the chain's order) instead of
        // the chain walk, in the same order, each `RING` words late through the ring
        // that prefetches their cells' offsets and runs (`lookup::x_dc_csr`). `on_word`
        // keeps its two callers, so its code is the one without the switch.
        if let Some(csr) = disc.x_csr() {
            let mut x_on_word_csr = |(kmer_start, s_range, index): (usize, usize, u64)| {
                if debug_enabled {
                    *dbg_total_s_positions += 1;
                }
                let run: &[u32] = if two_stage.has_hits(index) {
                    csr.run(index as usize)
                } else {
                    &[]
                };
                if (debug_mode || blastn_debug) && run.len() > 1000 {
                    eprintln!(
                        "[WARN] Large matches_slice: len={} for kmer at position {}",
                        run.len(),
                        kmer_start
                    );
                }
                if offset_pairs_len >= OFFSET_ARRAY_SIZE {
                    pair_extend::extend_offset_pairs(
                        &pair_inputs,
                        &mut pair_state,
                        offset_pairs,
                        &mut offset_pairs_len,
                        &mut scan_pause_ns,
                        true,
                    );
                }
                for &q_off_1 in run {
                    if debug_enabled {
                        *dbg_seeds_found += 1;
                    }
                    offset_pairs[offset_pairs_len] = OffsetPair {
                        q_off: q_off_1 as usize - 1,
                        s_off: kmer_start,
                        s_range,
                    };
                    offset_pairs_len += 1;
                }
                if offset_pairs_len >= offset_array_size {
                    pair_extend::extend_offset_pairs(
                        &pair_inputs,
                        &mut pair_state,
                        offset_pairs,
                        &mut offset_pairs_len,
                        &mut scan_pause_ns,
                        true,
                    );
                }
            };
            let mut ring: x_dc_csr::Ring<'_, (usize, usize, u64)> = x_dc_csr::Ring::new(csr);
            scan_subject_disc_words_with_ranges(
                search_seq_packed,
                s_len,
                disc,
                subject_seq_ranges,
                subject_masked,
                |s_off, s_range, index, index2| {
                    debug_assert!(index2.is_none());
                    ring.push(
                        (s_off, s_range, index),
                        |w| w.2 as usize,
                        &mut x_on_word_csr,
                    );
                },
            );
            ring.flush(&mut x_on_word_csr);
        } else {
            #[expect(
                clippy::redundant_closure,
                reason = "a callback closure of its own, as before the move (whose braces rustfmt drops now): `&mut on_word` would make the generic scan's callback another type"
            )]
            scan_subject_disc_words_with_ranges(
                search_seq_packed,
                s_len,
                disc,
                subject_seq_ranges,
                subject_masked,
                |s_off, s_range, index, index2| on_word(s_off, s_range, index, index2),
            );
        }
    } else if let Some(csr) = two_stage.x_mb_csr() {
        // EXPERIMENT (LOSAT_X_BNCSR): a contiguous table laid out as CSR has no
        // chains; its words with hits (the presence test of `for_each_hit` first)
        // go through `x_on_word_csr`, `on_word`'s statements with the cell's run
        // (the chain's offsets in the chain's order) instead of the chain walk, in
        // the same order, each `RING` words late through the ring that prefetches
        // their cells' offsets and runs (`lookup::x_dc_csr`). Every word still
        // counts for the debug summary, as in `on_word`. A batch of pairs that
        // reached `OFFSET_ARRAY_SIZE` goes before the next word with hits (`on_word`:
        // before the next word), with the same pairs, and before that word's warning,
        // so the debug output keeps its order. `on_word` keeps its two callers, so its
        // code is the one without the switch.
        // EXPERIMENT (LOSAT_X_BNCSRPACK): the form of the table's head is chosen
        // once (`x_dc_csr::XHead`), and `x_scan_csr!` expands to each form's scan,
        // so the scan of a head of one element per cell is the one without the
        // switch.
        // EXPERIMENT (LOSAT_X_BNPVONCE): a packed head's presence bits hold one bit
        // per cell, so its scan (`$pv` false) tests the word there only, not first
        // in the coarser presence vector, whose bit stands for several cells; for
        // a head whose bits fit the cache (`x_bn_pv_once`).
        // The trait's methods are called in two of the three scans below only, so
        // it is imported here and not in the macro's body.
        use x_dc_csr::XRing as _;
        macro_rules! x_scan_csr {
            ($head:expr, $pv:expr) => {{
                use x_dc_csr::XHead;
                let head = $head;
                let mut x_on_word_csr = |(kmer_start, s_range, index): (usize, usize, u64)| {
                    if offset_pairs_len >= OFFSET_ARRAY_SIZE {
                        pair_extend::extend_offset_pairs(
                            &pair_inputs,
                            &mut pair_state,
                            offset_pairs,
                            &mut offset_pairs_len,
                            &mut scan_pause_ns,
                            true,
                        );
                    }
                    let run: &[u32] = head.run(index as usize);
                    if (debug_mode || blastn_debug) && run.len() > 1000 {
                        eprintln!(
                            "[WARN] Large matches_slice: len={} for kmer at position {}",
                            run.len(),
                            kmer_start
                        );
                    }
                    for &q_off_1 in run {
                        if debug_enabled {
                            *dbg_seeds_found += 1;
                        }
                        offset_pairs[offset_pairs_len] = OffsetPair {
                            q_off: q_off_1 as usize - 1,
                            s_off: kmer_start,
                            s_range,
                        };
                        offset_pairs_len += 1;
                    }
                    if offset_pairs_len >= offset_array_size {
                        pair_extend::extend_offset_pairs(
                            &pair_inputs,
                            &mut pair_state,
                            offset_pairs,
                            &mut offset_pairs_len,
                            &mut scan_pause_ns,
                            true,
                        );
                    }
                };
                let mut ring = head.ring::<(usize, usize, u64)>();
                let mb_scan_kind = select_mb_scan_kind(lut_word_length, scan_step, subject_masked);
                scan_subject_kmers_with_ranges(
                    search_seq_packed,
                    s_len,
                    word_length,
                    lut_word_length,
                    scan_step,
                    &subject_seq_ranges,
                    subject_masked,
                    mb_scan_kind,
                    |kmer_start, s_range, current_lut_kmer| {
                        if debug_enabled {
                            *dbg_total_s_positions += 1;
                        }
                        if !$pv || two_stage.has_hits(current_lut_kmer) {
                            ring.push(
                                (kmer_start, s_range, current_lut_kmer),
                                |w| w.2 as usize,
                                &mut x_on_word_csr,
                            );
                        }
                    },
                );
                ring.flush(&mut x_on_word_csr);
            }};
        }
        match csr.x_packed() {
            Some(packed) if x_bn_pv_once(&packed) => x_scan_csr!(packed, false),
            Some(packed) => x_scan_csr!(packed, true),
            None => x_scan_csr!(csr, true),
        }
    } else {
        let mb_scan_kind = select_mb_scan_kind(lut_word_length, scan_step, subject_masked);
        scan_subject_kmers_with_ranges(
            search_seq_packed,
            s_len,
            word_length,
            lut_word_length,
            scan_step,
            subject_seq_ranges,
            subject_masked,
            mb_scan_kind,
            |kmer_start, s_range, current_lut_kmer| {
                on_word(kmer_start, s_range, current_lut_kmer, None)
            },
        );
    }

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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:451-497
    // ```c
    // hitsfound = scansub(lookup_wrap, subject, offset_pairs,
    //                     kScanSubjectOffsetArraySize, &scan_range[1]);
    // if (hitsfound == 0) continue;
    // hits_extended += extend(offset_pairs, hitsfound, ...);
    // ```
    if offset_pairs_len > 0 {
        pair_extend::extend_offset_pairs(
            &pair_inputs,
            &mut pair_state,
            offset_pairs,
            &mut offset_pairs_len,
            &mut scan_pause_ns,
            false,
        );
    }

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
    // DEBUG: Log stats for this subject
    if debug_mode || blastn_debug {
        if let Some(dbg_start_time) = dbg_start_time {
            let elapsed = dbg_start_time.elapsed();
            eprintln!("[PERF] Subject scan took {:?}: left_ext={}, right_ext={}, ungapped={}, seeds={}, valid_pos={}",
                    elapsed, dbg_left_ext_iters, dbg_right_ext_iters, dbg_ungapped_ext_calls,
                    dbg_seeds_found, dbg_total_s_positions);
        }
    }
}
