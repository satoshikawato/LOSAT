//! The preliminary search of one blastn subject's chunks, as `search_query_batch`'s
//! `process_subject` runs it for the preliminary stage (`collect_subject_prelim`): the subject in
//! chunks of at most `MAX_DBSEQ_LEN` letters that overlap by `DBSEQ_CHUNK_OVERLAP`
//! (`SubjectSplitState::next_chunk`), each chunk searched (`prelim_chunk_hits`: its scan and
//! ungapped extensions, by `two_stage_scan::scan_two_stage_lookup` or `na_scan::scan_na_lookup`
//! unless a strips dispatch of `x_pair_seed.rs` takes the chunk, then its HSP list,
//! `prelim_hsp_list.rs`) on this thread, or, with a pool of two threads or more, in batches of
//! chunks on the pool (a scratch per worker), and each chunk's list merged into the subject's in
//! chunk order (`merge_prelim_hit_lists`), then the subject best-hit filter.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-599 (the chunk loop of
//! `s_BlastSearchEngineOneContext`, :422-602)
//! ```c
//!     while (TRUE) {
//!         status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
//!                                        dbseq_chunk_overlap);
//!
//!         if (status == SUBJECT_SPLIT_DONE) break;
//!         if (status == SUBJECT_SPLIT_NO_RANGE) continue;
//! ...
//!         if (aux_struct->WordFinder) {
//!             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
//! ...
//!         status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
//!                      kHspNumMax, &(backup.offset), INT4_MIN,
//!                      overlap, score_options->gapped_calculation,
//!                      Blast_ProgramIsMapping(program_number));
//! ...
//!     } /* End loop on chunks of subject sequence */
//! ...
//!     *hsp_list_out_ptr = combined_hsp_list;
//! ```
//! NCBI searches a subject's chunks one after another (its threads take whole subjects); the pool's
//! chunk batches are LOSAT's, and their lists are merged in chunk order as that loop merges them.
//!
//! This module is a child of `run` for the items it keeps private (the types `SubjectChunk`,
//! `SubjectSplitState`, `SubjectScratch`, `PrelimHit`, `QueryContext`, `QueryContextIndex`,
//! `DiagStruct` and `BlastnTiming`; `advance_diag_table_offset`, `diag_hash_insert_window`,
//! `merge_prelim_hit_lists`, `prelim_subject_best_hit` and the child modules of the scans and of
//! the HSP list). Its functions are the statements of what were the closure
//! `collect_prelim_hits_for_chunk` of `process_subject` and the chunk loop after it, in their
//! order. What they read that is the same for every subject of a batch is `PrelimBatch`, built once
//! per batch, and what is the same for every chunk of a subject is `PrelimSubject`, built once per
//! subject; both are `Sync`, for the chunks on the pool.

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;
use std::sync::atomic::AtomicU64;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use super::x_few_subj;

use super::{
    advance_diag_table_offset, diag_hash_insert_window, merge_prelim_hit_lists, na_scan,
    prelim_hsp_list, prelim_subject_best_hit, two_stage_scan, x_diag_lv, x_pair_seed, BlastnTiming,
    DiagStruct, GapAlignScratch, HspListSplit, PrelimHit, QueryContext, QueryContextIndex,
    SubjectChunk, SubjectChunkStatus, SubjectScratch, SubjectSplitState, COMPRESSION_RATIO,
};
use crate::algorithm::blastn::alignment::gapped::BlastnaMatrix;
use crate::algorithm::blastn::alignment::x_greedy_prelim::XPrelimSubject;
use crate::algorithm::blastn::lookup::{NaLookupTable, PvDirectLookup, TwoStageLookup};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::utils::threading::SearchPool;

/// What the preliminary search of a subject's chunks reads that is the same for every subject of
/// the batch, in groups below. Built once per batch by `search_query_batch`, before
/// `process_subject`, with the names there; every field is a shared reference or a `Copy` value, so
/// the struct is `Sync` for the chunks on the pool.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:422-432
/// (`s_BlastSearchEngineOneContext`'s parameters)
/// ```c
/// s_BlastSearchEngineOneContext(EBlastProgramType program_number,
///         BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
///         BLAST_SequenceBlk* subject, Int4 orig_length, LookupTableWrap* lookup,
///         BlastGapAlignStruct* gap_align,
///         const BlastScoringParameters* score_params,
///         const BlastInitialWordParameters* word_params,
///         const BlastExtensionParameters* ext_params,
///         const BlastHitSavingParameters* hit_params,
///         BlastDiagnostics* diagnostics,
///         BlastCoreAuxStruct* aux_struct,
/// ```
/// The parallel shape and the diagnostics have no NCBI counterpart.
pub(super) struct PrelimBatch<'a> {
    /// The query block: the context of a query offset, the length of the concatenated query (every
    /// query, both strands) and of the block that chooses the diagonal table (array or hash), the
    /// concatenated query with and without its sentinels, each context's letters, and each query's
    /// length (the subject best-hit filter).
    pub(super) query_context_index: &'a QueryContextIndex,
    pub(super) query_concat_length: usize,
    pub(super) diag_container_query_length: usize,
    pub(super) encoded_query_concat_blastna: &'a [u8],
    pub(super) encoded_query_concat_blastna_with_sentinels: &'a [u8],
    pub(super) encoded_queries_blastna: &'a [&'a [u8]],
    pub(super) query_lengths: &'a [usize],
    /// The word finder: the word size, the scan step, the off-diagonal range, the window (0: one
    /// hit), the size at which a batch of offset pairs goes to the steps of the pairs, whether a
    /// table whose width is the word size is the direct one, the small-query table's compressed
    /// query (empty for the other tables) and whether its aligned one-byte form applies, and the
    /// minimum ungapped score that the debug window's lines print.
    pub(super) effective_word_size: usize,
    pub(super) scan_step: usize,
    pub(super) scan_range: usize,
    pub(super) window_size: usize,
    pub(super) offset_array_size: usize,
    pub(super) use_direct_lookup: bool,
    pub(super) small_na_compressed_query: &'a [u8],
    pub(super) small_na_aligned_one_byte: bool,
    pub(super) min_ungapped_score: i32,
    /// The diagonal array's length and mask (NCBI's `diag_array_length`, `diag_mask`).
    pub(super) diag_array_length: usize,
    pub(super) diag_mask: usize,
    /// The scoring and the gapped aligner: the score matrix, the score table of the approximate
    /// extension, the reward and penalty, the gap costs, the aligner (dynamic programming or
    /// greedy) and its X-dropoff, and the interval tree's diagonal separation.
    pub(super) score_matrix: &'a BlastnaMatrix,
    pub(super) nucl_score_table: &'a [i32; 256],
    pub(super) reward: i32,
    pub(super) penalty: i32,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    pub(super) use_dp: bool,
    pub(super) x_drop_gapped: i32,
    pub(super) min_diag_separation: i32,
    /// The initial hits that the word finder saved in every chunk of the batch (NCBI's
    /// `good_init_extends`; they size the next batch).
    pub(super) good_init_extends: &'a AtomicU64,
    /// The parallel shape: whether the pool has two threads or more, whether it searches the
    /// subjects in parallel, its threads and the pool, the query group tasks (LOSAT_X_QGROUPPAR),
    /// the look-ahead's window (`None` when the helpers do not run) and the rule of
    /// LOSAT_X_GREEDYSPEC (kept as its `x_reach!` key), and the strips of LOSAT_X_BNFEWSUBJ and
    /// their kind (read by the pool's chunk batches only, so a build without the parallel feature
    /// has neither).
    pub(super) requested_parallel: bool,
    pub(super) use_parallel: bool,
    pub(super) num_threads: usize,
    pub(super) parallel_pool: &'a SearchPool<'a>,
    pub(super) x_qgroup_batch: bool,
    pub(super) x_ahead_window: Option<usize>,
    pub(super) x_greedy_spec: bool,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) x_strips: bool,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) x_kind: x_few_subj::Kind,
    /// The subject best-hit filter (`-subject_besthit`).
    pub(super) subject_besthit: bool,
    /// The diagnostics: BLEMIR_DEBUG, LOSAT_DEBUG_BLASTN, the trace lines (`blastn_trace`), the
    /// verbose lines and the debug window.
    pub(super) debug_mode: bool,
    pub(super) blastn_debug: bool,
    pub(super) blastn_trace_enabled: bool,
    pub(super) verbose: bool,
    pub(super) debug_window: Option<(usize, usize, usize, usize)>,
}

/// What the preliminary search of a subject's chunks reads that is the same for every chunk of the
/// subject, in groups below. Built once per subject by `process_subject`'s preliminary call, with
/// the names there (references and `Copy` values only, `Sync` as `PrelimBatch`); its borrows end
/// with the call.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:424,460-469,475
/// (`s_BlastSearchEngineOneContext`: the subject, its chunks' overlap and the backup of the subject
/// that its chunks are cut from)
/// ```c
///         BLAST_SequenceBlk* subject, Int4 orig_length, LookupTableWrap* lookup,
/// ...
///     if (Blast_ProgramIsMapping(program_number) &&
///         (Int4)query_info->max_length < 110) {
///
///         dbseq_chunk_overlap = (Int4)query_info->max_length
///             + (query_info->max_length / 2);
///     }
///     else {
///         dbseq_chunk_overlap = DBSEQ_CHUNK_OVERLAP;
///     }
/// ...
///     s_BackupSubject(subject, &backup);
/// ```
/// The trace id, the one-byte letters and the timing counters have no NCBI counterpart.
pub(super) struct PrelimSubject<'a> {
    /// The subject: its index in the batch, its trace id (made only for the debug and trace lines),
    /// its packed letters, its letters one per byte for the greedy aligner (`XPrelimSubject`),
    /// whether it has lowercase masks, and the overlap of its chunks.
    pub(super) s_idx: usize,
    pub(super) s_id: &'a str,
    pub(super) s_seq_packed_full: &'a [u8],
    pub(super) x_prelim_subject: &'a XPrelimSubject<'a>,
    pub(super) subject_masked: bool,
    pub(super) dbseq_chunk_overlap: usize,
    /// The batch's lookup tables, which the preliminary stage's calls of `process_subject` pass:
    /// the two-stage table, or the direct table or the na table.
    pub(super) two_stage_lookup_ref: Option<&'a TwoStageLookup>,
    pub(super) pv_direct_lookup_ref: Option<&'a PvDirectLookup>,
    pub(super) na_lookup_ref: Option<&'a NaLookupTable>,
    /// The batch's values as `process_subject` reads them: the queries, their contexts, each
    /// context's cutoff, hit-saving cutoff, X-dropoff and reduced cutoff, and the LOSAT_TIMING
    /// counters (`None` when timing is off).
    pub(super) queries: &'a [FastaRecord],
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) cutoff_scores: &'a [i32],
    pub(super) hit_saving_cutoff_scores: &'a [i32],
    pub(super) x_dropoff_scores: &'a [i32],
    pub(super) reduced_cutoff_scores: &'a [i32],
    pub(super) timing_ref: Option<&'a BlastnTiming>,
    pub(super) timing_enabled: bool,
}

/// The search of one chunk of a subject (one turn of the chunk loop): the chunk's ranges and packed
/// letters, the debug counters, the interval tree, the diagonal table (array or hash) set up for
/// the chunk, then the word finder: for a two-stage table the strips dispatch of PAIRPAR
/// (`x_pair_seed::run`) or the scan (`two_stage_scan::scan_two_stage_lookup`), for the other tables
/// the dispatch of LOSAT_X_PAIRPARNA (`x_pair_seed::run_na`) or the scan
/// (`na_scan::scan_na_lookup`); then the chunk's HSP list, of an ungapped search from its initial
/// hits (`prelim_hsp_list::ungapped_hsp_list`), else from the initial hits sorted
/// (`prelim_hsp_list::sort_ungapped_hits_for_gapped`) by the gapped loop with its purge and sort
/// (`prelim_hsp_list::gapped_hsp_list`); and the diagonal table's offset moved on. The list is in
/// the subject's coordinates. `reuse_prelim_hits` is set for the owner's serial calls (not a chunk
/// of a batch on the pool); the scratches are the subject's, or the pool worker's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:487-500,522-526,555,560-563
/// (one turn of the chunk loop of `s_BlastSearchEngineOneContext`)
/// ```c
///         /* Delete if not done in last loop iteration to prevent memory leak. */
///         hsp_list = Blast_HSPListFree(hsp_list);
///
///         BlastInitHitListReset(init_hitlist);
///
///         if (aux_struct->WordFinder) {
///             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
///                                    word_params, aux_struct->ewp,
///                                    aux_struct->offset_pairs,
///                                    kScanSubjectOffsetArraySize,
///                                    init_hitlist, ungapped_stats);
///
///             if (init_hitlist->total == 0) continue;
///         }
/// ...
///         if (aux_struct->GetGappedScore) {
///             status = aux_struct->GetGappedScore(program_number, query,
///                     query_info,
///                     subject, gap_align, score_params, ext_params, hit_params,
///                     word_params, init_hitlist, &hsp_list, gapped_stats, NULL);
///         }
/// ...
///         Blast_HSPListSortByScore(hsp_list);
/// ...
///         } else {
///             BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
///                     hit_params->options, &hsp_list);
///         }
/// ```
pub(super) fn prelim_chunk_hits<const UNGAPPED: bool>(
    batch: &PrelimBatch<'_>,
    subject: &PrelimSubject<'_>,
    chunk: &SubjectChunk,
    soft_ranges: &[(i32, i32)],
    gap_scratch: &mut GapAlignScratch,
    subject_scratch: &mut SubjectScratch,
    reuse_prelim_hits: bool,
) -> Vec<PrelimHit> {
    let &PrelimBatch {
        query_context_index,
        query_concat_length,
        diag_container_query_length,
        encoded_query_concat_blastna,
        encoded_query_concat_blastna_with_sentinels,
        encoded_queries_blastna,
        effective_word_size,
        scan_step,
        scan_range,
        window_size,
        offset_array_size,
        use_direct_lookup,
        small_na_compressed_query,
        small_na_aligned_one_byte,
        min_ungapped_score,
        diag_array_length,
        diag_mask,
        score_matrix,
        nucl_score_table,
        reward,
        penalty,
        gap_open,
        gap_extend,
        use_dp,
        x_drop_gapped,
        min_diag_separation,
        good_init_extends,
        use_parallel,
        num_threads,
        parallel_pool,
        x_qgroup_batch,
        x_ahead_window,
        x_greedy_spec,
        debug_mode,
        blastn_debug,
        blastn_trace_enabled,
        verbose,
        debug_window,
        ..
    } = batch;
    let &PrelimSubject {
        s_idx,
        s_id,
        s_seq_packed_full,
        x_prelim_subject,
        two_stage_lookup_ref,
        pv_direct_lookup_ref,
        na_lookup_ref,
        queries,
        query_contexts,
        cutoff_scores,
        hit_saving_cutoff_scores,
        x_dropoff_scores,
        reduced_cutoff_scores,
        timing_ref,
        timing_enabled,
        ..
    } = subject;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
    // ```c
    // /* Delete if not done in last loop iteration to prevent memory leak. */
    // hsp_list = Blast_HSPListFree(hsp_list);
    //
    // BlastInitHitListReset(init_hitlist);
    // ```
    let hits = &mut subject_scratch.hits;
    hits.clear();
    let prelim_hits = &mut subject_scratch.prelim_hits;
    prelim_hits.clear();
    let greedy_align_scratch = &mut subject_scratch.greedy_align_scratch;

    let s_len = chunk.length;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:264-307
    // ```c
    // if (backup->offset == 0 && residual == 0 && backup->next == backup->full_range.right) {
    //     subject->seq_ranges = backup->soft_ranges;
    //     subject->num_seq_ranges = backup->num_soft_ranges;
    //     return SUBJECT_SPLIT_OK;
    // }
    // ...
    // for (i=0; i<len; i++) {
    //     subject->seq_ranges[i].left = backup->soft_ranges[i+start].left - backup->offset;
    //     subject->seq_ranges[i].right = backup->soft_ranges[i+start].right - backup->offset;
    // }
    // ```
    let seq_ranges_scratch = &mut subject_scratch.seq_ranges_scratch;
    let subject_seq_ranges = chunk.seq_ranges.as_slice(
        soft_ranges,
        seq_ranges_scratch,
        chunk.offset as i32,
        s_len as i32,
    );
    let subject_masked = chunk.masked;

    let packed_offset = chunk.offset / COMPRESSION_RATIO;
    let s_seq_packed = &s_seq_packed_full[packed_offset..];

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
    let debug_enabled = debug_mode || blastn_debug;

    // Debug counters for this subject (only updated when debug_enabled)
    let mut dbg_total_s_positions = 0usize;
    let dbg_ambiguous_skipped = 0usize;
    let dbg_no_lookup_match = 0usize;
    let mut dbg_seeds_found = 0usize;
    let mut dbg_ungapped_low = 0usize;
    let dbg_two_hit_failed = 0usize;
    let mut dbg_gapped_attempted = 0usize;
    let mut dbg_window_seeds = 0usize;

    let safe_k = effective_word_size.min(31);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:147-259
    // ```c
    // /* Scan the compressed subject sequence */
    // index = s[0] << 16 | s[1] << 8 | s[2];
    // index = (index >> shift) & mask;
    // ```

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_stat.h:866-869 (query blastna, subject ncbi2na)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847 (subject blastna + ncbi2na)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2949-3016 (packed subject for score-only DP)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:503-507 (traceback uses uncompressed subject)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:806-833 (GetReverseNuclSequence uses ncbi4na/blastna)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85-93 (IUPACNA_TO_BLASTNA)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:148-349 (packed ncbi2na used for ungapped extension)

    // NCBI reference: blast_gapalign.c:3826-3831 Blast_IntervalTreeInit
    // Initialize interval tree for HSP containment checking
    // Tree is indexed by query offsets (primary) and subject offsets (midpoint subtrees)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3811-3831
    // ```c
    // tree = Blast_IntervalTreeInit(0, query->length+1,
    //                               0, subject->length+1);
    // ```
    // query->length is the concatenated length for both strands (2*len + 1).
    // The tree is the scratch's (`x_tree`), given this chunk's ranges.
    let interval_tree = &mut subject_scratch.x_tree;
    interval_tree.reset_with_bounds(
        0,                                // q_min
        (query_concat_length + 1) as i32, // q_max
        0,                                // s_min
        (s_len + 1) as i32,               // s_max
    );

    // NCBI architecture: Collect all ungapped hits first, then process in score-descending order
    // NCBI reference: blast_gapalign.c:3824 - ASSERT(Blast_InitHitListIsSortedByScore(init_hitlist))
    // This is critical for correct containment checking: high-score HSPs must be processed first
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
    // ```c
    // hsp_list = Blast_HSPListFree(hsp_list);
    // BlastInitHitListReset(init_hitlist);
    // ```
    let ungapped_hits = &mut subject_scratch.ungapped_hits;
    ungapped_hits.clear();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-64
    // ```c
    // diag_array_length = 1;
    // while (diag_array_length < (qlen+window_size))
    //     diag_array_length = diag_array_length << 1;
    // diag_table->diag_array_length = diag_array_length;
    // diag_table->diag_mask = diag_array_length-1;
    // diag_table->offset = window_size;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1002-1003
    // ```c
    // if ((status = BlastExtendWordNew(query->length, word_params,
    //                                 &aux_struct->ewp)) != 0)
    // ```
    // One diagonal table for the whole query block (every query and both
    // strands): `query->length` is the last context's offset plus its
    // length, `query_concat_length` here, and the diagonals use offsets in
    // the concatenated query.
    let (diag_table_length, diag_table_mask) = (diag_array_length, diag_mask as isize);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:166-233
    // ```c
    // const int kQueryLenForHashTable = 8000; /* For blastn, use hash table rather
    //                                         than diag array for any query longer
    //                                         than this */
    // if (Blast_ProgramIsNucleotide(program_number) &&
    //     !Blast_QueryIsPattern(program_number) &&
    //     (query_info->contexts[query_info->last_context].query_offset +
    //      query_info->contexts[query_info->last_context].query_length) > kQueryLenForHashTable)
    //     p->container_type = eDiagHash;
    // else
    //     p->container_type = eDiagArray;
    // ```
    let use_diag_hash = diag_container_query_length > 8000;
    let use_array_indexing = !use_diag_hash;
    let diag_window = window_size as i32;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:786-833
    // ```c
    // first_context = 1;
    // last_context = 1;
    // ...
    // subject->frame = context;
    // ```
    // BLASTN scans subject plus strand only; query contexts cover both strands.
    for _subject_strand in [false] {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
        // ```c
        // BlastSeqBlkSetSequence(subj, sequence.data.release(),
        //    ((sentinels == eSentinels) ? sequence.length - 2 :
        //     sequence.length));
        // ...
        // BlastSeqBlkSetCompressedSequence(subj,
        //                                  compressed_seq.data.release());
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:148-349 (packed subject used by ungapped extension)
        let search_seq_packed: &[u8] = s_seq_packed;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:656-662
        // ```c
        // s_off_pos = s_off + diag_table->offset;
        // ```
        let diag_offset = if use_array_indexing {
            subject_scratch.diag_table_offset as isize
        } else {
            subject_scratch.diag_hash.offset as isize
        };
        let diag_array_size = if use_array_indexing {
            diag_table_length
        } else {
            0
        };

        // NCBI BLAST does NOT use a separate mask array for diagonal suppression
        // NCBI only uses last_hit in hit_level_array (checked at line 753)
        // This mask_array/mask_hash was LOSAT-specific and caused excessive filtering
        // REMOVED to match NCBI BLAST behavior

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:84-105
        // ```c
        // n = diag->diag_array_length;
        // diag->offset = diag->window;
        // diag_struct_array = diag->hit_level_array;
        // for (i = 0; i < n; i++) {
        //     diag_struct_array[i].flag = 0;
        //     diag_struct_array[i].last_hit = -diag->window;
        //     if (diag->hit_len_array) diag->hit_len_array[i] = 0;
        // }
        // ```
        // NCBI reference: blast_extend.h:77-80, na_ungapped.c:660-666
        // Two-hit tracking: DiagStruct array for tracking last_hit and flag per diagonal
        // hit_level_array[diag] contains {last_hit, flag} (equivalent to NCBI's DiagStruct)
        // hit_len_array[diag] contains hit length (0 = no hit, >0 = hit length)
        // For single query: use Vec for O(1) access, otherwise use HashMap
        let hit_level_array = &mut subject_scratch.hit_level_array;
        let hit_len_array = &mut subject_scratch.hit_len_array;
        let diag_hash = &mut subject_scratch.diag_hash;
        // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the two retrieves of the
        // diagonal hash in the scans called below (`two_stage_scan::scan_two_stage_lookup`,
        // through `pair_extend::extend_offset_pairs`, for a two-stage table,
        // `na_scan::scan_na_lookup` for the other tables) skip the walk of
        // a bucket whose cells all lie below the level bound, unless diagnostics print what
        // they read.
        let x_lv_ref = diag_hash.x_lv_on()
            && !(debug_enabled || blastn_trace_enabled || debug_window.is_some());
        if diag_hash.x_lv_on() {
            crate::utils::xreach::x_reach!("LOSAT_X_DIAGLV", use_diag_hash && x_lv_ref);
            // EXPERIMENT (LOSAT_X_DIAGLV2): the reference folds called below
            // (`pair_extend::extend_offset_pairs`, `na_scan::scan_na_lookup`) never read
            // the memo (the strips stage of dc-megablast records its own side).
            if x_diag_lv::lv2_on() {
                crate::utils::xreach::x_reach!("LOSAT_X_DIAGLV2", false);
            }
        }
        if use_array_indexing {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:145-149
            // ```c
            // diag_table->hit_level_array = (DiagStruct *)
            //     calloc(diag_table->diag_array_length, sizeof(DiagStruct));
            // if (word_params->options->window_size) {
            //     diag_table->hit_len_array = (Uint1 *)
            //          calloc(diag_table->diag_array_length, sizeof(Uint1));
            // }
            // ```
            if hit_level_array.len() != diag_array_size {
                let had_entries = !hit_level_array.is_empty();
                let diag_default = DiagStruct::default();
                hit_level_array.resize(diag_array_size, diag_default);
                if had_entries {
                    hit_level_array.fill(diag_default);
                }
            }
            if window_size > 0 {
                if hit_len_array.len() != diag_array_size {
                    let had_entries = !hit_len_array.is_empty();
                    hit_len_array.resize(diag_array_size, 0);
                    if had_entries {
                        hit_len_array.fill(0);
                    }
                }
            } else if !hit_len_array.is_empty() {
                hit_len_array.clear();
            }
        } else {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:159-182
            // ```c
            // if (ewp->hash_table->offset >= INT4_MAX / 4) {
            //     ewp->hash_table->occupancy = 1;
            //     ewp->hash_table->offset = ewp->hash_table->window;
            //     memset(ewp->hash_table->backbone, 0,
            //            ewp->hash_table->num_buckets * sizeof(Int4));
            // } else {
            //     ewp->hash_table->offset += subject_length + ewp->hash_table->window;
            // }
            // ```
            // Hash-based diagonals persist across chunks; only clear on offset overflow.
            hit_level_array.clear();
            hit_len_array.clear();
        }

        if s_len < effective_word_size {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:159-176
            // ```c
            // if (ewp->diag_table->offset >= INT4_MAX / 4) {
            //     ewp->diag_table->offset = ewp->diag_table->window;
            //     s_BlastDiagClear(ewp->diag_table);
            // } else {
            //     ewp->diag_table->offset += subject_length + ewp->diag_table->window;
            // }
            // ```
            advance_diag_table_offset(
                &mut subject_scratch.diag_table_offset,
                diag_window,
                s_len,
                use_array_indexing,
                hit_level_array,
                hit_len_array,
                diag_hash,
            );
            return Vec::new();
        }

        // TWO-STAGE LOOKUP: Use separate rolling k-mer for lut_word_length
        if let Some(two_stage) = two_stage_lookup_ref {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
            // ```c
            //     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
            //
            //         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
            // ...
            //         hits_extended += extend(offset_pairs, hitsfound, word_params,
            //                                 lookup_wrap, query, subject, matrix,
            //                                 query_info, ewp, init_hitlist, scan_range[2] + lut_word_length);
            //     }
            // ```
            // Dispatch point of PAIRPAR: `x_pair_seed::run` does what
            // `two_stage_scan::scan_two_stage_lookup` does (scan, word extension, diagonal
            // table, ungapped extension; in the scopes of `x_pair_seed.rs` and, for a
            // discontiguous table, `x_pair_seed_disc.rs`) with the subject scan in strips
            // and the pairs folded by this thread in the scan order; `continue` ends the
            // one-iteration strand loop. Otherwise `two_stage_scan::scan_two_stage_lookup`
            // runs.
            if x_pair_seed::run(
                &x_pair_seed::Inputs {
                    two_stage,
                    search_seq_packed,
                    s_len,
                    subject_seq_ranges,
                    subject_masked,
                    scan_step,
                    query_context_index,
                    query_contexts,
                    encoded_query_concat_blastna,
                    encoded_query_concat_blastna_with_sentinels,
                    cutoff_scores,
                    x_dropoff_scores,
                    reduced_cutoff_scores,
                    score_matrix,
                    nucl_score_table,
                    diag_offset,
                    diag_hash_window: diag_hash_insert_window(
                        window_size,
                        scan_range,
                        two_stage.word_length(),
                    ),
                    use_array_indexing,
                    window_size,
                    scan_range,
                    small_na_word: !small_na_compressed_query.is_empty(),
                    diagnostics: debug_enabled || blastn_trace_enabled || debug_window.is_some(),
                    one_subject_owner: !use_parallel && reuse_prelim_hits,
                    timing: timing_ref,
                },
                diag_hash,
                x_pair_seed::DiagArray {
                    levels: &mut *hit_level_array,
                    lens: &mut *hit_len_array,
                    length: diag_table_length as isize,
                    mask: diag_table_mask,
                    size: diag_array_size,
                },
                ungapped_hits,
                &mut subject_scratch.x_pair_fold,
            ) {
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-498
            // ```c
            //         if (aux_struct->WordFinder) {
            //             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
            //                                    word_params, aux_struct->ewp,
            //                                    aux_struct->offset_pairs,
            //                                    kScanSubjectOffsetArraySize,
            //                                    init_hitlist, ungapped_stats);
            // ```
            // The reference scan of this chunk and the steps of the offset pairs it collects
            // (`two_stage_scan::scan_two_stage_lookup`), with what they read and write, built once per
            // chunk.
            two_stage_scan::scan_two_stage_lookup(
                &two_stage_scan::TwoStageScanInputs {
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
                },
                &mut two_stage_scan::TwoStageScanState {
                    ungapped_hits: &mut *ungapped_hits,
                    hit_level_array: &mut *hit_level_array,
                    hit_len_array: &mut *hit_len_array,
                    diag_hash: &mut *diag_hash,
                    offset_pairs: &mut subject_scratch.offset_pairs,
                    dbg_total_s_positions: &mut dbg_total_s_positions,
                    dbg_seeds_found: &mut dbg_seeds_found,
                    dbg_window_seeds: &mut dbg_window_seeds,
                    dbg_ungapped_low: &mut dbg_ungapped_low,
                    dbg_gapped_attempted: &mut dbg_gapped_attempted,
                },
            );
        }
        if two_stage_lookup_ref.is_none() {
            // Original lookup method (for non-two-stage lookup)
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
            // Dispatch point of LOSAT_X_PAIRPARNA: `x_pair_seed::run_na` does what
            // `na_scan::scan_na_lookup` does (scan, diagonal hash, ungapped extension;
            // one-hit diagonal-hash scope only) with the subject scan in strips on the pool
            // threads and the pairs folded in the scan order; `continue` ends the
            // one-iteration strand loop. Otherwise `na_scan::scan_na_lookup` runs.
            if x_pair_seed::run_na(
                &x_pair_seed::NaInputs {
                    pv_direct_lookup: pv_direct_lookup_ref,
                    na_lookup: na_lookup_ref,
                    use_direct_lookup,
                    word_length: safe_k,
                    search_seq_packed,
                    s_len,
                    subject_seq_ranges,
                    subject_masked,
                    scan_step,
                    query_context_index,
                    query_contexts,
                    encoded_queries_blastna,
                    cutoff_scores,
                    x_dropoff_scores,
                    reduced_cutoff_scores,
                    score_matrix,
                    nucl_score_table,
                    diag_offset,
                    diag_hash_window: diag_hash_insert_window(window_size, scan_range, safe_k),
                    use_array_indexing,
                    window_size,
                    scan_range,
                    diagnostics: debug_enabled || blastn_trace_enabled || debug_window.is_some(),
                    one_subject_owner: !use_parallel && reuse_prelim_hits,
                    timing: timing_ref,
                },
                diag_hash,
                ungapped_hits,
                &mut subject_scratch.x_pair_fold,
            ) {
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
            // ```c
            //     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
            //
            //         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
            // ...
            //         hits_extended += extend(offset_pairs, hitsfound, word_params,
            // ...
            //     }
            // ```
            // The reference scan of this chunk and the steps of the seeds it finds
            // (`na_scan::scan_na_lookup`), with what they read and write, built once per
            // chunk.
            na_scan::scan_na_lookup(
                &na_scan::NaScanInputs {
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
                },
                &mut na_scan::NaScanState {
                    ungapped_hits: &mut *ungapped_hits,
                    hit_level_array: &mut *hit_level_array,
                    hit_len_array: &mut *hit_len_array,
                    diag_hash: &mut *diag_hash,
                    dbg_total_s_positions: &mut dbg_total_s_positions,
                    dbg_seeds_found: &mut dbg_seeds_found,
                    dbg_window_seeds: &mut dbg_window_seeds,
                    dbg_ungapped_low: &mut dbg_ungapped_low,
                    dbg_gapped_attempted: &mut dbg_gapped_attempted,
                },
            );
        } // end of if/else for two-stage vs original lookup
    } // end of strand loop

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:560-563
    // ```c
    //         } else {
    //             BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
    //                     hit_params->options, &hsp_list);
    //         }
    // ```
    // An ungapped (`-ungapped`) search makes the chunk's HSP list from its initial hits
    // (`prelim_hsp_list::ungapped_hsp_list`); the diagonal table then moves on, as after
    // a gapped search's chunk.
    if UNGAPPED {
        let chunk_prelim_hits =
            prelim_hsp_list::ungapped_hsp_list(ungapped_hits, chunk, good_init_extends);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:159-176
        // (the diagonal table's offset moves on after each chunk, as below)
        advance_diag_table_offset(
            &mut subject_scratch.diag_table_offset,
            diag_window,
            s_len,
            use_array_indexing,
            &mut subject_scratch.hit_level_array,
            &mut subject_scratch.hit_len_array,
            &mut subject_scratch.diag_hash,
        );
        return chunk_prelim_hits;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1688-1692
    // ```c
    //     Blast_UngappedStatsUpdate(ungapped_stats, total_hits, hits_extended,
    //                               init_hitlist->total);
    //
    //     if (word_params->ungapped_extension)
    //         Blast_InitHitListSortByScore(init_hitlist);
    // ```
    // The initial hits sorted by score for the gapped loop, or grouped by query for its
    // query group tasks (LOSAT_X_QGROUPPAR), and counted
    // (`prelim_hsp_list::sort_ungapped_hits_for_gapped`).
    let (x_counts, x_order, total_ungapped) = prelim_hsp_list::sort_ungapped_hits_for_gapped(
        ungapped_hits,
        x_qgroup_batch,
        reuse_prelim_hits,
        queries,
        good_init_extends,
    );
    let gapped_start_time = std::time::Instant::now();

    // The values the gapped preliminary loop of this chunk
    // (`prelim_hsp_list::gapped_score_list`) and one hit's extension
    // (`prelim_hsp_list::x_prelim_extend`) read, for the query group tasks and the
    // look-ahead helpers on the pool as for this thread.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3664-3674
    // ```c
    // Int2 BLAST_GetGappedScore (EBlastProgramType program_number,
    //         BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
    //         BLAST_SequenceBlk* subject,
    //         BlastGapAlignStruct* gap_align,
    //         const BlastScoringParameters* score_params,
    //         const BlastExtensionParameters* ext_params,
    //         const BlastHitSavingParameters* hit_params,
    //         const BlastInitialWordParameters* word_params,
    //         BlastInitHitList* init_hitlist,
    //         BlastHSPList** hsp_list_ptr, BlastGappedStats* gapped_stats,
    //         Boolean * fence_hit)
    // ```
    // No NCBI counterpart: the look-ahead's schedule (scheduling only).
    let gapped = prelim_hsp_list::GappedPrelim {
        blastn_trace_enabled,
        verbose,
        s_idx,
        s_id,
        x_prelim_subject,
        chunk,
        s_seq_packed,
        s_len,
        encoded_queries_blastna,
        query_contexts,
        hit_saving_cutoff_scores,
        reward,
        penalty,
        gap_open,
        gap_extend,
        score_matrix,
        use_dp,
        x_drop_gapped,
        min_diag_separation,
        reuse_prelim_hits,
        total_ungapped,
        x_ahead_window,
        x_greedy_spec,
        num_threads,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:522-526,545,555
    // ```c
    //         if (aux_struct->GetGappedScore) {
    //             status = aux_struct->GetGappedScore(program_number, query,
    //                     query_info,
    //                     subject, gap_align, score_params, ext_params, hit_params,
    //                     word_params, init_hitlist, &hsp_list, gapped_stats, NULL);
    // ...
    //             Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    // ...
    //         Blast_HSPListSortByScore(hsp_list);
    // ```
    // The chunk's gapped HSP list (`prelim_hsp_list::gapped_hsp_list`): the gapped loop per
    // query on the pool (LOSAT_X_QGROUPPAR) or over the chunk's list on this thread with
    // the subject's scratch, then the purge, the sort and the chunk's offset.
    let chunk_prelim_hits = prelim_hsp_list::gapped_hsp_list(
        &gapped,
        &prelim_hsp_list::GappedDrive {
            x_counts: &x_counts,
            x_order: &x_order,
            queries,
            query_concat_length,
            parallel_pool,
            timing_ref,
            gapped_start_time,
        },
        ungapped_hits,
        prelim_hsp_list::GappedListState {
            interval_tree,
            gap_scratch,
            greedy_align_scratch,
            x_helper_scratch: Some(&mut subject_scratch.x_ahead_scratch),
            prelim_hits,
        },
    );
    // Suppress unused variable warnings when not in debug mode
    let _ = (
        dbg_total_s_positions,
        dbg_ambiguous_skipped,
        dbg_no_lookup_match,
        dbg_seeds_found,
        dbg_ungapped_low,
        dbg_two_hit_failed,
        dbg_gapped_attempted,
        dbg_window_seeds,
        total_ungapped,
    );
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:159-176
    // ```c
    // if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //     ewp->diag_table->offset = ewp->diag_table->window;
    //     s_BlastDiagClear(ewp->diag_table);
    // } else {
    //     ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    // }
    // ```
    let hit_level_array = &mut subject_scratch.hit_level_array;
    let hit_len_array = &mut subject_scratch.hit_len_array;
    let diag_hash = &mut subject_scratch.diag_hash;
    advance_diag_table_offset(
        &mut subject_scratch.diag_table_offset,
        diag_window,
        s_len,
        use_array_indexing,
        hit_level_array,
        hit_len_array,
        diag_hash,
    );
    chunk_prelim_hits
}

/// The preliminary search of one subject's chunks and the merge of their lists: the chunks that
/// `split_state` cuts (`SubjectSplitState::next_chunk`, skipping a chunk without a range), each
/// searched by `prelim_chunk_hits` with the subject's scratches on this thread, or, with a pool of
/// two threads or more (`requested_parallel`), taken in batches of twice the pool's threads: one
/// chunk, or the serial chunks of LOSAT_X_BNFEWSUBJ, on this thread; else the batch's chunks on the
/// pool, a scratch per worker, their lists collected in chunk order. Each chunk's list is merged
/// into the subject's (`merge_prelim_hit_lists`, with the chunk's offset and overlap), then the
/// subject best-hit filter (`-subject_besthit`). Returns the subject's list.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-483,578-593
/// (the chunk loop of `s_BlastSearchEngineOneContext`)
/// ```c
///     while (TRUE) {
///         status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
///                                        dbseq_chunk_overlap);
///
///         if (status == SUBJECT_SPLIT_DONE) break;
///         if (status == SUBJECT_SPLIT_NO_RANGE) continue;
/// ...
///         Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
///         overlap = (backup.offset == backup.hard_ranges[backup.hm_index].left) ?
///                   0 : dbseq_chunk_overlap;
///         status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
///                      kHspNumMax, &(backup.offset), INT4_MIN,
///                      overlap, score_options->gapped_calculation,
///                      Blast_ProgramIsMapping(program_number));
///
///         if((hit_params->options->hsp_filt_opt != NULL) &&
///            (hit_params->options->hsp_filt_opt->subject_besthit_opts != NULL)) {
/// ...
///     } /* End loop on chunks of subject sequence */
/// ```
/// The chunk batches on the pool have no NCBI counterpart (NCBI's threads take whole subjects).
pub(super) fn collect_subject_prelim<const UNGAPPED: bool>(
    batch: &PrelimBatch<'_>,
    subject: &PrelimSubject<'_>,
    split_state: &mut SubjectSplitState,
    gap_scratch: &mut GapAlignScratch,
    subject_scratch: &mut SubjectScratch,
) -> Vec<PrelimHit> {
    let &PrelimBatch {
        requested_parallel,
        subject_besthit,
        query_lengths,
        ..
    } = batch;
    let &PrelimSubject {
        subject_masked,
        dbseq_chunk_overlap,
        ..
    } = subject;
    let mut combined_prelim_hits: Vec<PrelimHit> = Vec::new();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:478-536
    // status = s_GetNextSubjectChunk(subject, &backup, kNucleotide, ...);
    // A single subject can contain independent chunks. Pool availability
    // is separate from the outer subject traversal's number of jobs.
    if requested_parallel {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        {
            // What only the pool's chunk batches read (a build without the parallel feature
            // has none of them).
            let &PrelimBatch {
                num_threads,
                use_parallel,
                x_strips,
                x_kind,
                offset_array_size,
                window_size,
                ..
            } = batch;
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:478-536
            // ```c
            // while (TRUE) {
            //     status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
            //                                    dbseq_chunk_overlap);
            //     if (status == SUBJECT_SPLIT_DONE) break;
            //     if (status == SUBJECT_SPLIT_NO_RANGE) continue;
            //     ...
            // }
            // ```
            let max_in_flight = num_threads.max(1) * 2;
            let mut chunk_batch: Vec<SubjectChunk> = Vec::with_capacity(max_in_flight);
            let mut done = false;
            loop {
                chunk_batch.clear();
                while chunk_batch.len() < max_in_flight && !done {
                    match split_state.next_chunk(subject_masked, dbseq_chunk_overlap) {
                        SubjectChunkStatus::Done => done = true,
                        SubjectChunkStatus::NoRange => continue,
                        SubjectChunkStatus::Ok(chunk) => chunk_batch.push(chunk),
                    }
                }
                if chunk_batch.is_empty() {
                    break;
                }

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:264-268
                // ```c
                // subject->seq_ranges = backup->soft_ranges;
                // subject->num_seq_ranges = backup->num_soft_ranges;
                // return SUBJECT_SPLIT_OK;
                // ```
                // NCBI reference: c++/src/algo/blast/core/blast_engine.c:478-500
                // status = s_GetNextSubjectChunk(subject, &backup, ...);
                // EXPERIMENT (LOSAT_X_BNFEWSUBJ): a batch whose longest chunk is a large
                // part of its letters (`x_few_subj::serial_chunks`), which this subject's
                // owner searches by itself, searches its chunks in order as one-chunk
                // batches do (strips for each chunk), with the one-thread loop's
                // statements below.
                let x_serial_chunks = !use_parallel
                    && !subject_masked
                    && x_strips
                    && x_few_subj::serial_chunks(
                        chunk_batch.iter().map(|chunk| chunk.length),
                        num_threads,
                        x_kind,
                    );
                if !use_parallel && !subject_masked && x_strips && chunk_batch.len() > 1 {
                    crate::utils::xreach::x_reach!("LOSAT_X_BNFEWSUBJ", x_serial_chunks);
                }
                crate::utils::threading::report_stage(
                    "blastn",
                    "subject_chunks",
                    chunk_batch.len(),
                    chunk_batch.len() > 1 && !x_serial_chunks,
                );
                let soft_ranges = split_state.soft_ranges.as_slice();
                if chunk_batch.len() == 1 || x_serial_chunks {
                    for chunk in chunk_batch.iter() {
                        let hits = prelim_chunk_hits::<UNGAPPED>(
                            batch,
                            subject,
                            chunk,
                            soft_ranges,
                            gap_scratch,
                            subject_scratch,
                            true,
                        );
                        let hits = merge_prelim_hit_lists(
                            &mut combined_prelim_hits,
                            hits,
                            HspListSplit::Subject {
                                offset: chunk.offset,
                            },
                            chunk.overlap,
                            !UNGAPPED,
                        );
                        if subject_besthit {
                            prelim_subject_best_hit(&mut combined_prelim_hits, query_lengths);
                        }
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
                        // ```c
                        // /* Delete if not done in last loop iteration to prevent memory leak. */
                        // hsp_list = Blast_HSPListFree(hsp_list);
                        //
                        // BlastInitHitListReset(init_hitlist);
                        // ```
                        subject_scratch.prelim_hits = hits;
                    }
                } else {
                    let mut chunk_results: Vec<(usize, usize, Vec<PrelimHit>)> =
                        Vec::with_capacity(chunk_batch.len());
                    chunk_batch
                        .par_iter()
                        .map_init(
                            || {
                                (
                                    GapAlignScratch::new(),
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:991-1041
                                    // ```c
                                    // Int4 offset_array_size = GetOffsetArraySize(lookup_wrap);
                                    // ...
                                    // aux_struct->offset_pairs =
                                    //   (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
                                    // ```
                                    SubjectScratch::new(offset_array_size, window_size),
                                )
                            },
                            |state, chunk| {
                                let (gap_scratch, subject_scratch) = state;
                                let hits = prelim_chunk_hits::<UNGAPPED>(
                                    batch,
                                    subject,
                                    chunk,
                                    soft_ranges,
                                    gap_scratch,
                                    subject_scratch,
                                    false,
                                );
                                (chunk.offset, chunk.overlap, hits)
                            },
                        )
                        .collect_into_vec(&mut chunk_results);
                    for (offset, overlap, hits) in chunk_results {
                        let _ = merge_prelim_hit_lists(
                            &mut combined_prelim_hits,
                            hits,
                            HspListSplit::Subject { offset },
                            overlap,
                            !UNGAPPED,
                        );
                        if subject_besthit {
                            prelim_subject_best_hit(&mut combined_prelim_hits, query_lengths);
                        }
                    }
                }

                if done {
                    break;
                }
            }
        }
        #[cfg(any(
            not(feature = "parallel"),
            all(target_arch = "wasm32", not(feature = "wasm-threads"))
        ))]
        {
            unreachable!("use_parallel is false when parallel threads are disabled");
        }
    } else {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:478-536
        // ```c
        // while (TRUE) {
        //     status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
        //                                    dbseq_chunk_overlap);
        //     if (status == SUBJECT_SPLIT_DONE) break;
        //     if (status == SUBJECT_SPLIT_NO_RANGE) continue;
        //     ...
        // }
        // ```
        loop {
            match split_state.next_chunk(subject_masked, dbseq_chunk_overlap) {
                SubjectChunkStatus::Done => break,
                SubjectChunkStatus::NoRange => continue,
                SubjectChunkStatus::Ok(chunk) => {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:264-268
                    // ```c
                    // subject->seq_ranges = backup->soft_ranges;
                    // subject->num_seq_ranges = backup->num_soft_ranges;
                    // return SUBJECT_SPLIT_OK;
                    // ```
                    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:478-500
                    // status = s_GetNextSubjectChunk(subject, &backup, ...);
                    crate::utils::threading::report_stage("blastn", "subject_chunks", 1, false);
                    let soft_ranges = split_state.soft_ranges.as_slice();
                    let hits = prelim_chunk_hits::<UNGAPPED>(
                        batch,
                        subject,
                        &chunk,
                        soft_ranges,
                        gap_scratch,
                        subject_scratch,
                        true,
                    );
                    let hits = merge_prelim_hit_lists(
                        &mut combined_prelim_hits,
                        hits,
                        HspListSplit::Subject {
                            offset: chunk.offset,
                        },
                        chunk.overlap,
                        !UNGAPPED,
                    );
                    if subject_besthit {
                        prelim_subject_best_hit(&mut combined_prelim_hits, query_lengths);
                    }
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
                    // ```c
                    // /* Delete if not done in last loop iteration to prevent memory leak. */
                    // hsp_list = Blast_HSPListFree(hsp_list);
                    //
                    // BlastInitHitListReset(init_hitlist);
                    // ```
                    subject_scratch.prelim_hits = hits;
                }
            }
        }
    }
    combined_prelim_hits
}
