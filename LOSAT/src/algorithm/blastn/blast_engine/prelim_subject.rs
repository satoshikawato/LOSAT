//! The preliminary search of one blastn subject, as `SubjectPasses::search_subject` calls
//! it for the preliminary stage (`prelim_subject`): the subject's lowercase masks or database masks
//! (`-db_soft_mask`, `-db_hard_mask`) and the set-up of its chunks, then the chunk loop (`collect_subject_prelim`):
//! the subject in chunks of at most `MAX_DBSEQ_LEN` letters that overlap by `DBSEQ_CHUNK_OVERLAP`
//! (`SubjectSplitState::next_chunk`), each chunk searched (`prelim_chunk_hits`: its scan and
//! ungapped extensions, by `two_stage_scan::scan_two_stage_lookup` or `na_scan::scan_na_lookup`
//! unless a strips dispatch of `x_pair_seed.rs` takes the chunk, then its HSP list,
//! `prelim_hsp_list.rs`) on this thread, or, with a pool of two threads or more, in batches of
//! jobs, runs of consecutive chunks of at most `MAX_DBSEQ_LEN` letters (a batch of one job, or
//! LOSAT_X_BNFEWSUBJ's serial chunks, on this thread; else the batch's jobs on the pool, a scratch
//! per worker, each job's chunks in order), and each chunk's list merged into the subject's in
//! chunk order (`merge_prelim_hit_lists`), then the subject best-hit filter; then the subject's
//! list reaped by the preliminary e-value or, under matrix-only scoring, by raw score
//! (`reap_prelim`), or, in an ungapped search, finished (`finish_ungapped_prelim`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-599 (the chunk loop of
//! `s_BlastSearchEngineOneContext`, :422-602)
//! NCBI searches a subject's chunks one after another (its threads take whole subjects); the pool's
//! chunk batches and their jobs are LOSAT's, and their lists are merged in chunk order as that loop
//! merges them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-900 (the e-values and the
//! reap of the subject's list at the end of `s_BlastSearchEngineCore`) and :1480-1542 (the finish
//! of an ungapped search's list in `BLAST_PreliminarySearchEngine`).
//!
//! This module is a child of `run` for the items it keeps private (the types `SubjectChunk`,
//! `SubjectSplitState`, `SubjectScratch`, `PrelimHit`, `QueryContext`, `QueryContextIndex`,
//! `DiagStruct`, `BlastnTiming`, `BlastnSubjectRef` and `SubjectCtx`; `advance_diag_table_offset`,
//! `diag_hash_insert_window`, `merge_prelim_hit_lists`, `prelim_subject_best_hit`,
//! `build_subject_seq_ranges_from_masks`, `read_db_soft_ranges`, `ungapped_subject_stage`, the
//! e-value helpers, `database` and the child modules of the scans and of the HSP list). Its
//! functions are the statements of `search_subject`'s preliminary side (the subject's set-up, the
//! reap and the `-ungapped` finish), of what was its closure `collect_prelim_hits_for_chunk` and of
//! the chunk loop, in their order. What they read of a batch is `PrelimBatch` and, with the batch's
//! lookup tables, `ScanBatch` (the batch's, or `special_subject`'s copies with the statistics of
//! each subject); what `search_subject` read of the call is `SubjectCtx`; what is the
//! same for every chunk of a subject is `PrelimSubject`, built once per subject. `PrelimBatch`,
//! `ScanBatch` and `PrelimSubject` are `Sync`, for the chunks on the pool.

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
use super::{x_few_subj, MAX_DBSEQ_LEN};

use super::{
    advance_diag_table_offset, build_subject_seq_ranges_from_masks,
    calculate_blastn_context_statistics, database, diag_hash_insert_window,
    encode_iupac_to_blastna, merge_prelim_hit_lists, na_scan, next_subject_chunk, prelim_hsp_list,
    prelim_subject_best_hit, read_db_seq_ranges, two_stage_scan, ungapped_subject_stage, x_diag_lv,
    x_pair_seed, BlastnSubjectRef, BlastnTiming, ContextKarlin, DiagStruct, GapAlignScratch,
    HspListSplit, PrelimHit, PrelimReap, QueryContext, QueryContextIndex, SubjectChunk,
    SubjectChunkStatus, SubjectCtx, SubjectScratch, SubjectSplitState, COMPRESSION_RATIO,
    DBSEQ_CHUNK_OVERLAP,
};
use crate::algorithm::blastn::alignment::gapped::BlastnaMatrix;
use crate::algorithm::blastn::alignment::x_greedy_prelim::XPrelimSubject;
use crate::algorithm::blastn::coordination::{collect_lowercase_masks, LookupTables};
use crate::algorithm::blastn::lookup::{NaLookupTable, PvDirectLookup, TwoStageLookup};
use crate::algorithm::blastn::x_subject_store::x_blastna_into;
use crate::blastinput::db_mask::MaskKind;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::utils::threading::SearchPool;

/// What the preliminary search of a subject reads of its batch, in groups below. Built once per
/// batch by `search_query_batch`, before `search_subject`, with the names there; with the
/// statistics of each subject also once per subject, by `SubjectPasses::special_subject`: a copy
/// with the subject's own cutoffs and search spaces (every field derived from the statistics must
/// be set in that copy, and in its `ScanBatch`). Every field is a shared reference or a `Copy`
/// value, so the struct is `Sync` for the chunks on the pool.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:422-432
/// (`s_BlastSearchEngineOneContext`'s parameters)
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
    /// The queries and their contexts, each context's effective search space, cutoff, hit-saving
    /// cutoff, X-dropoff and reduced cutoff of the preliminary stage (`batch_setup::batch_cutoffs`:
    /// the setup's; the stage's own where the statistics' total is 0 or less; the subject's own in
    /// `special_subject`'s copy), and the LOSAT_TIMING counters (`None` when timing is off).
    pub(super) queries: &'a [FastaRecord],
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) query_eff_searchsp: &'a [i64],
    pub(super) cutoff_scores: &'a [i32],
    pub(super) hit_saving_cutoff_scores: &'a [i32],
    pub(super) x_dropoff_scores: &'a [i32],
    pub(super) reduced_cutoff_scores: &'a [i32],
    pub(super) timing_ref: Option<&'a BlastnTiming>,
    /// The word finder: the word size, the scan step, the off-diagonal range, the window (0: one
    /// hit), the size at which a batch of offset pairs goes to the steps of the pairs, whether a
    /// table whose width is the word size is the direct one, the small-query table's compressed
    /// query (empty for the other tables) and whether its aligned one-byte form applies, and the
    /// minimum ungapped score that the debug window's lines print.
    pub(super) effective_word_size: usize,
    pub(super) scan_step: usize,
    pub(super) scan_range: usize,
    pub(super) window_size: usize,
    /// The window when it is below 0 (an imported search strategy's), else 0: the batch's word
    /// finder is then `neg_window.rs`'s.
    pub(super) neg_window: i32,
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
    /// NCBI's matrix-only scoring (`TaskConfig::matrix_only`): every ungapped extension is the
    /// exact one, and the preliminary lists are reaped by raw score (`PrelimReap`).
    pub(super) matrix_only: bool,
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
    /// The subject's masks: its lowercase masks (`-lcase_masking`); with `-db_soft_mask` or
    /// `-db_hard_mask`, the database, the algorithm whose masks a database sequence's search reads
    /// and the kind, and the first error of those reads by ordinal id (the batch returns it after
    /// the preliminary stage).
    pub(super) lcase_masking: bool,
    pub(super) db_mask: Option<(&'a database::BlastnDatabase, i32, MaskKind)>,
    pub(super) db_mask_error: &'a std::sync::Mutex<Option<(u32, anyhow::Error)>>,
    /// The seeds of the MegaBLAST database index that the batch is searched with (`None`
    /// without one, and for a batch split into query chunks: each chunk has its own), which
    /// replace the scan of every chunk (`mb_index_seeds::index_chunk_hits`).
    pub(super) mb_index: Option<&'a super::mb_index_seeds::IndexSeeds<'a>>,
    /// The subject's list after its chunks (the reap, the `-ungapped` finish): the concatenated
    /// query's letters without its masks (NCBI's `sequence_nomask`, as
    /// `encoded_query_concat_blastna` is laid out), the Karlin blocks by context, whether an
    /// e-value is rounded down from its score, the reap (by the e-value threshold, or by raw
    /// score under matrix-only scoring) and the minimum percent identity.
    pub(super) encoded_query_concat_nomask: &'a [u8],
    pub(super) search_karlin_ref: &'a [ContextKarlin],
    pub(super) round_down_evalue_score: bool,
    pub(super) prelim_reap: PrelimReap,
    pub(super) percent_identity: f64,
    /// The diagnostics: BLEMIR_DEBUG, LOSAT_DEBUG_BLASTN, the trace lines (`blastn_trace`), the
    /// verbose lines and the debug window.
    pub(super) debug_mode: bool,
    pub(super) blastn_debug: bool,
    pub(super) blastn_trace_enabled: bool,
    pub(super) verbose: bool,
    pub(super) debug_window: Option<(usize, usize, usize, usize)>,
}

/// What the word finder of every chunk of a batch reads (the scan, the steps of its offset pairs
/// and the ungapped extension; `two_stage_scan.rs`, `pair_extend.rs`, `na_scan.rs`), with the
/// batch's lookup tables and `PrelimBatch` for the rest of the chunk's search: the values of
/// `PrelimBatch` that the word finder reads, one load from the pointer that each chunk's word
/// finder takes, and what follows from them for every chunk (the diagonal table's form and sizes,
/// the word lengths, the level bound's constants). Built once per batch by the preliminary pass
/// (`SubjectPasses::prelim_pass`, in `prelim_into`; with the statistics of each subject also once
/// per subject, `special_subject`'s copy pointing at its `PrelimBatch` copy with the subject's
/// cutoff, X-dropoff and reduced cutoff), which holds the tables (the traceback drops
/// them, `x_lut_drop.rs`; `ScanBatch::new`); shared references and `Copy` values only (`Sync`, for
/// the subjects and the chunks on the pool). Each chunk's inputs are a reference to it and the
/// chunk's own values (its letters, length, ranges, masks, diagonal offset and ids).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-498 (the word finder of a
/// chunk takes pointers to the search's shared structs: the query, its contexts, the lookup table,
/// the matrix and the word parameters)
/// and c++/src/algo/blast/core/na_ungapped.c:1597-1608 (`BlastNaWordFinder`'s parameters, read
/// through those pointers where they are used). The level bound, the diagnostics and the timing
/// counters have no NCBI counterpart.
pub(super) struct ScanBatch<'a> {
    /// The rest of the batch's values (the gapped stage, the merges, the parallel shape).
    pub(super) batch: &'a PrelimBatch<'a>,
    /// The batch's lookup tables: the two-stage table, or the direct table or the na table.
    pub(super) two_stage_lookup_ref: Option<&'a TwoStageLookup>,
    pub(super) pv_direct_lookup_ref: Option<&'a PvDirectLookup>,
    pub(super) na_lookup_ref: Option<&'a NaLookupTable>,
    /// The query block (`PrelimBatch`'s): the context of a query offset, the contexts, the
    /// concatenated query with and without its sentinels (the scans read the concatenated query,
    /// as NCBI's), each context's cutoff, X-dropoff and reduced cutoff, the score matrix and the
    /// score table of the approximate extension.
    pub(super) query_context_index: &'a QueryContextIndex,
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) encoded_query_concat_blastna: &'a [u8],
    pub(super) encoded_query_concat_blastna_with_sentinels: &'a [u8],
    pub(super) cutoff_scores: &'a [i32],
    pub(super) x_dropoff_scores: &'a [i32],
    pub(super) reduced_cutoff_scores: &'a [i32],
    pub(super) score_matrix: &'a BlastnaMatrix,
    pub(super) nucl_score_table: &'a [i32; 256],
    /// Whether the ungapped extension is the exact one (`s_NuclUngappedExtendExact`) rather than
    /// the approximate one with its exact recomputation: under matrix-only scoring or for a word
    /// shorter than 11 letters. One value of the batch, read by every word finder.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:740-742 (the same at
    /// 910-912)
    pub(super) exact_ungapped: bool,
    /// The words: the word size capped at 31 letters (`safe_k`, the na tables' word), the word and
    /// lookup word lengths (the two-stage table's, else `safe_k`), the scan step, the off-diagonal
    /// range, the window and whether it is a two-hit search, the size at which a batch of offset
    /// pairs goes to the steps of the pairs, whether a table whose width is the word size is the
    /// direct one, the small-query table's compressed query (empty for the other tables) and
    /// whether its aligned one-byte form applies, and the level bound's constants of the seeds
    /// (LOSAT_X_DIAGLV, `x_diag_lv.rs`).
    pub(super) safe_k: usize,
    pub(super) word_length: usize,
    pub(super) lut_word_length: usize,
    pub(super) scan_step: usize,
    pub(super) scan_range: usize,
    pub(super) window_size: usize,
    /// `PrelimBatch::neg_window`.
    pub(super) neg_window: i32,
    /// The word finders that replace the scan of a chunk (`prelim_chunk_hits`), `None` for a
    /// batch with neither.
    pub(super) special_finders: Option<SpecialFinders<'a>>,
    pub(super) two_hits: bool,
    pub(super) offset_array_size: usize,
    pub(super) use_direct_lookup: bool,
    pub(super) small_na_compressed_query: &'a [u8],
    pub(super) small_na_aligned_one_byte: bool,
    pub(super) x_lv_seeds: x_diag_lv::LvSeeds,
    /// The diagonal table: whether it is the array (a query block of at most 8000 letters) or the
    /// hash; the array's length and mask (`PrelimBatch`'s), the same for the steps of the pairs
    /// (0 with the hash) and the array's size (0 with the hash); the hash's insert window of the
    /// two-stage table's word (0 with the array) and the window as the offset's step.
    pub(super) use_array_indexing: bool,
    pub(super) diag_table_length: usize,
    pub(super) diag_table_mask: isize,
    pub(super) diag_array_length: isize,
    pub(super) diag_mask: isize,
    pub(super) diag_array_size: usize,
    pub(super) diag_hash_window: i32,
    pub(super) diag_window: i32,
    /// The diagnostics: BLEMIR_DEBUG, LOSAT_DEBUG_BLASTN and either of them (the debug counters),
    /// the trace lines, the debug window, any of these (`diagnostics`), the minimum ungapped score
    /// that the debug window's lines print, and the LOSAT_TIMING counters (`None` when timing is
    /// off).
    pub(super) debug_mode: bool,
    pub(super) blastn_debug: bool,
    pub(super) debug_enabled: bool,
    pub(super) blastn_trace_enabled: bool,
    pub(super) debug_window: Option<(usize, usize, usize, usize)>,
    pub(super) diagnostics: bool,
    pub(super) min_ungapped_score: i32,
    pub(super) timing_ref: Option<&'a BlastnTiming>,
    pub(super) timing_enabled: bool,
}

/// The word finders of a batch that replace the scan of a chunk (`prelim_chunk_hits`): the
/// MegaBLAST database index's and that of a window below 0. A batch has them as one value
/// (`ScanBatch::special_finders`), so that the chunks of a batch with neither test one value, as
/// they tested the index alone.
#[derive(Clone, Copy)]
pub(super) struct SpecialFinders<'a> {
    /// The index's seeds (`PrelimBatch::mb_index`).
    pub(super) index: Option<&'a super::mb_index_seeds::IndexSeeds<'a>>,
    /// A window below 0 (`PrelimBatch::neg_window`): the word finder of `neg_window.rs`.
    pub(super) neg_window: bool,
}

impl<'a> ScanBatch<'a> {
    /// The word finder's values of `batch` with its lookup tables `tables`.
    pub(super) fn new(batch: &'a PrelimBatch<'a>, tables: &'a LookupTables) -> Self {
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:131-156
        let two_stage_lookup_ref = tables.two_stage_lookup.as_ref();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:166-233
        let use_array_indexing = batch.diag_container_query_length <= 8000;
        let (diag_table_length, diag_table_mask) =
            (batch.diag_array_length, batch.diag_mask as isize);
        let safe_k = batch.effective_word_size.min(31);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1625-1637
        let (word_length, lut_word_length) = two_stage_lookup_ref
            .map_or((safe_k, safe_k), |table| {
                (table.word_length(), table.lut_word_length())
            });
        let window_size = batch.window_size;
        let two_hits = window_size > 0;
        let debug_enabled = batch.debug_mode || batch.blastn_debug;
        let exact_ungapped = batch.matrix_only || word_length < 11;
        ScanBatch {
            batch,
            two_stage_lookup_ref,
            pv_direct_lookup_ref: tables.pv_direct_lookup.as_ref(),
            na_lookup_ref: tables.na_lookup.as_ref(),
            query_context_index: batch.query_context_index,
            query_contexts: batch.query_contexts,
            encoded_query_concat_blastna: batch.encoded_query_concat_blastna,
            encoded_query_concat_blastna_with_sentinels: batch
                .encoded_query_concat_blastna_with_sentinels,
            cutoff_scores: batch.cutoff_scores,
            x_dropoff_scores: batch.x_dropoff_scores,
            reduced_cutoff_scores: batch.reduced_cutoff_scores,
            score_matrix: batch.score_matrix,
            nucl_score_table: batch.nucl_score_table,
            exact_ungapped,
            safe_k,
            word_length,
            lut_word_length,
            scan_step: batch.scan_step,
            scan_range: batch.scan_range,
            window_size,
            neg_window: batch.neg_window,
            special_finders: (batch.mb_index.is_some() || batch.neg_window != 0).then_some(
                SpecialFinders {
                    index: batch.mb_index,
                    neg_window: batch.neg_window != 0,
                },
            ),
            two_hits,
            offset_array_size: batch.offset_array_size,
            use_direct_lookup: batch.use_direct_lookup,
            small_na_compressed_query: batch.small_na_compressed_query,
            small_na_aligned_one_byte: batch.small_na_aligned_one_byte,
            // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the level bound's constants of the
            // seeds (`s_end` adds the word length).
            x_lv_seeds: x_diag_lv::LvSeeds::new(two_hits, word_length, window_size),
            use_array_indexing,
            diag_table_length,
            diag_table_mask,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:663-664
            diag_array_length: if use_array_indexing {
                diag_table_length as isize
            } else {
                0
            },
            diag_mask: if use_array_indexing {
                diag_table_mask
            } else {
                0
            },
            diag_array_size: if use_array_indexing {
                diag_table_length
            } else {
                0
            },
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825-826, 939-941
            diag_hash_window: if use_array_indexing {
                0
            } else {
                diag_hash_insert_window(window_size, batch.scan_range, word_length)
            },
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:668
            diag_window: window_size as i32,
            debug_mode: batch.debug_mode,
            blastn_debug: batch.blastn_debug,
            debug_enabled,
            blastn_trace_enabled: batch.blastn_trace_enabled,
            debug_window: batch.debug_window,
            diagnostics: debug_enabled
                || batch.blastn_trace_enabled
                || batch.debug_window.is_some(),
            min_ungapped_score: batch.min_ungapped_score,
            timing_ref: batch.timing_ref,
            timing_enabled: batch.timing_ref.is_some(),
        }
    }
}

/// What the preliminary search of a subject's chunks reads that is the same for every chunk of the
/// subject, in groups below. Built once per subject by `prelim_subject`, with the names there
/// (references and `Copy` values only, `Sync` as `PrelimBatch`); its borrows end with the call.
/// What it reads of the batch comes from `ScanBatch` and `PrelimBatch` (or their copies).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:425,461-469,476
/// (`s_BlastSearchEngineOneContext`: the subject, its chunks' overlap and the backup of the subject
/// that its chunks are cut from)
/// The trace id and the one-byte letters have no NCBI counterpart.
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
    /// A database sequence of a search with a MegaBLAST database index and `-db_hard_mask`: the
    /// index's seeds placed into its hard chunks (`mb_index_seeds::HardSeeds`).
    pub(super) hard_index: Option<&'a super::mb_index_seeds::HardSeeds>,
}

/// The search of one chunk of a subject (one turn of the chunk loop): the chunk's ranges and packed
/// letters, the debug counters, the diagonal table (array or hash) set up for the chunk, then the
/// word finder: for a two-stage table the strips dispatch of PAIRPAR (`x_pair_seed::run`) or the
/// scan (`two_stage_scan::scan_two_stage_lookup`), for the other tables the dispatch of
/// LOSAT_X_PAIRPARNA (`x_pair_seed::run_na`) or the scan (`na_scan::scan_na_lookup`). A chunk
/// without an initial hit ends there with an empty list, after its diagonal table's offset moved
/// on (NCBI's `continue` after the word finder, below). Otherwise the chunk's HSP list, of an
/// ungapped search from its initial hits (`prelim_hsp_list::ungapped_hsp_list`), else from the
/// initial hits sorted (`prelim_hsp_list::sort_ungapped_hits_for_gapped`) by the gapped loop with
/// its interval tree, purge and sort (`prelim_hsp_list::gapped_hsp_list`); and the diagonal
/// table's offset moved on. The list is in the subject's coordinates. `reuse_prelim_hits` is set for the owner's serial calls (not a chunk
/// of a batch on the pool); the scratches are the subject's, or the pool worker's.
///
/// The word finder of a batch with an index, or with a window below 0 (`PrelimBatch::neg_window`),
/// replaces the scan (`ScanBatch::special_finders`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:487-500,522-527,555,560-563
/// (one turn of the chunk loop of `s_BlastSearchEngineOneContext`)
pub(super) fn prelim_chunk_hits<const UNGAPPED: bool>(
    scan: &ScanBatch<'_>,
    subject: &PrelimSubject<'_>,
    chunk: &SubjectChunk,
    soft_ranges: &[(i32, i32)],
    gap_scratch: &mut GapAlignScratch,
    subject_scratch: &mut SubjectScratch,
    reuse_prelim_hits: bool,
) -> Vec<PrelimHit> {
    // What every chunk of the batch reads is `scan`'s (its word finder) and `scan.batch`'s (the
    // rest), read where it is used; what the chunk's subject gives is `subject`'s.
    let batch = scan.batch;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1006-1012
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1727-1734
    // A batch searched with a MegaBLAST database index has the indexed word finder, which scans
    // a sequence of a database volume without an index as without the index.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1597-1695
    // A window below 0 has the word finder of `neg_window.rs`, after the index's.
    if let Some(special) = scan.special_finders {
        if let Some(index) = special.index {
            if index.is_indexed(subject.s_idx as u32) {
                return super::mb_index_seeds::index_chunk_hits::<UNGAPPED>(
                    scan,
                    index,
                    subject,
                    chunk,
                    gap_scratch,
                    subject_scratch,
                    reuse_prelim_hits,
                );
            }
        }
        if special.neg_window {
            return super::neg_window::neg_chunk_hits::<UNGAPPED>(
                scan,
                subject,
                chunk,
                soft_ranges,
                gap_scratch,
                subject_scratch,
                reuse_prelim_hits,
            );
        }
    }
    let &PrelimSubject {
        s_idx,
        s_id,
        s_seq_packed_full,
        ..
    } = subject;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:488-491
    let hits = &mut subject_scratch.hits;
    hits.clear();
    subject_scratch.prelim_hits.clear();

    let s_len = chunk.length;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:264-307
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
    // (`scan.debug_enabled`: BLEMIR_DEBUG or LOSAT_DEBUG_BLASTN.)

    // Debug counters for this subject (only updated when debug_enabled)
    let mut dbg_total_s_positions = 0usize;
    let dbg_ambiguous_skipped = 0usize;
    let dbg_no_lookup_match = 0usize;
    let mut dbg_seeds_found = 0usize;
    let mut dbg_ungapped_low = 0usize;
    let dbg_two_hit_failed = 0usize;
    let mut dbg_gapped_attempted = 0usize;
    let mut dbg_window_seeds = 0usize;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:147-259

    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_stat.h:866-869 (query blastna, subject ncbi2na)
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847 (subject blastna + ncbi2na)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2949-3016 (packed subject for score-only DP)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-507 (traceback uses uncompressed subject)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:806-833 (GetReverseNuclSequence uses ncbi4na/blastna)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:85-93 (IUPACNA_TO_BLASTNA)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:148-349 (packed ncbi2na used for ungapped extension)
    // The interval tree of the gapped stage is set up there, for a chunk with initial hits only
    // (`prelim_hsp_list::gapped_hsp_list`, as `BLAST_GetGappedScore` makes its tree after its
    // return for an empty list, blast_gapalign.c:3706-3707,3826-3831).

    // NCBI architecture: Collect all ungapped hits first, then process in score-descending order
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3824
    // Asserts that `init_hitlist` is sorted by score (`Blast_InitHitListIsSortedByScore`).
    // This is critical for correct containment checking: high-score HSPs must be processed first
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:488-491
    let ungapped_hits = &mut subject_scratch.ungapped_hits;
    ungapped_hits.clear();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-64
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1003
    // One diagonal table for the whole query block (every query and both
    // strands): `query->length` is the last context's offset plus its
    // length, `query_concat_length` here, and the diagonals use offsets in
    // the concatenated query (`scan.diag_table_length`, `scan.diag_table_mask`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:166-233
    // (`scan.use_array_indexing`, set once per batch, `ScanBatch::new`.)
    let use_array_indexing = scan.use_array_indexing;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:786-833
    // BLASTN scans subject plus strand only; query contexts cover both strands.
    for _subject_strand in [false] {
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:148-349 (packed subject used by ungapped extension)
        let search_seq_packed: &[u8] = s_seq_packed;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:668
        let diag_offset = if use_array_indexing {
            subject_scratch.diag_table_offset as isize
        } else {
            subject_scratch.diag_hash.offset as isize
        };
        let diag_array_size = scan.diag_array_size;

        // NCBI BLAST does NOT use a separate mask array for diagonal suppression
        // NCBI only uses last_hit in hit_level_array (checked at line 753)
        // This mask_array/mask_hash was LOSAT-specific and caused excessive filtering
        // REMOVED to match NCBI BLAST behavior

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:84-105
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:77-80, c++/src/algo/blast/core/na_ungapped.c:660-666
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
        let x_lv_ref = diag_hash.x_lv_on() && !scan.diagnostics;
        // The sides are recorded only when LOSAT_X_REACH is set (`xreach::recording`, one test for
        // the two marks of every chunk).
        if crate::utils::xreach::recording() && diag_hash.x_lv_on() {
            crate::utils::xreach::x_reach!("LOSAT_X_DIAGLV", !use_array_indexing && x_lv_ref);
            // EXPERIMENT (LOSAT_X_DIAGLV2): the reference folds called below
            // (`pair_extend::extend_offset_pairs`, `na_scan::scan_na_lookup`) never read
            // the memo (the strips stage of dc-megablast records its own side).
            if x_diag_lv::lv2_on() {
                crate::utils::xreach::x_reach!("LOSAT_X_DIAGLV2", false);
            }
        }
        if use_array_indexing {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-149
            if hit_level_array.len() != diag_array_size {
                let had_entries = !hit_level_array.is_empty();
                let diag_default = DiagStruct::default();
                hit_level_array.resize(diag_array_size, diag_default);
                if had_entries {
                    hit_level_array.fill(diag_default);
                }
            }
            if scan.window_size > 0 {
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:159-182
            // Hash-based diagonals persist across chunks; only clear on offset overflow.
            hit_level_array.clear();
            hit_len_array.clear();
        }

        if s_len < batch.effective_word_size {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:159-176
            advance_diag_table_offset(
                &mut subject_scratch.diag_table_offset,
                scan.diag_window,
                s_len,
                use_array_indexing,
                hit_level_array,
                hit_len_array,
                diag_hash,
            );
            return Vec::new();
        }

        // TWO-STAGE LOOKUP: Use separate rolling k-mer for lut_word_length
        if let Some(two_stage) = scan.two_stage_lookup_ref {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
            // Dispatch point of PAIRPAR: `x_pair_seed::run` does what
            // `two_stage_scan::scan_two_stage_lookup` does (scan, word extension, diagonal
            // table, ungapped extension; in the scopes of `x_pair_seed.rs` and, for a
            // discontiguous table, `x_pair_seed_disc.rs`) with the subject scan in strips
            // and the pairs folded by this thread in the scan order; `continue` ends the
            // one-iteration strand loop. Otherwise `two_stage_scan::scan_two_stage_lookup`
            // runs. `x_pair_seed::wants` asks `run`'s tests of a contiguous table first, so a
            // chunk that `run` would send back (every chunk of a database search) builds no
            // inputs and makes no call.
            if x_pair_seed::wants(
                two_stage.disc().is_some(),
                use_array_indexing,
                scan.window_size,
                scan.scan_range,
                !scan.small_na_compressed_query.is_empty(),
                scan.diagnostics,
                !batch.use_parallel && reuse_prelim_hits,
            ) && x_pair_seed::run(
                &x_pair_seed::Inputs {
                    two_stage,
                    search_seq_packed,
                    s_len,
                    subject_seq_ranges,
                    subject_masked,
                    scan_step: scan.scan_step,
                    query_context_index: scan.query_context_index,
                    query_contexts: scan.query_contexts,
                    encoded_query_concat_blastna: scan.encoded_query_concat_blastna,
                    encoded_query_concat_blastna_with_sentinels: scan
                        .encoded_query_concat_blastna_with_sentinels,
                    cutoff_scores: scan.cutoff_scores,
                    x_dropoff_scores: scan.x_dropoff_scores,
                    reduced_cutoff_scores: scan.reduced_cutoff_scores,
                    score_matrix: scan.score_matrix,
                    nucl_score_table: scan.nucl_score_table,
                    exact_ungapped: scan.exact_ungapped,
                    diag_offset,
                    diag_hash_window: diag_hash_insert_window(
                        scan.window_size,
                        scan.scan_range,
                        two_stage.word_length(),
                    ),
                    use_array_indexing,
                    window_size: scan.window_size,
                    scan_range: scan.scan_range,
                    small_na_word: !scan.small_na_compressed_query.is_empty(),
                    diagnostics: scan.diagnostics,
                    one_subject_owner: !batch.use_parallel && reuse_prelim_hits,
                    timing: scan.timing_ref,
                },
                diag_hash,
                x_pair_seed::DiagArray {
                    levels: &mut *hit_level_array,
                    lens: &mut *hit_len_array,
                    length: scan.diag_table_length as isize,
                    mask: scan.diag_table_mask,
                    size: diag_array_size,
                },
                ungapped_hits,
                &mut subject_scratch.x_pair_fold,
            ) {
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-498
            // The reference scan of this chunk and the steps of the offset pairs it collects
            // (`two_stage_scan::scan_two_stage_lookup`), with what they read and write, built once per
            // chunk.
            two_stage_scan::scan_two_stage_lookup(
                &two_stage_scan::TwoStageScanInputs {
                    scan,
                    two_stage,
                    search_seq_packed,
                    s_len,
                    subject_seq_ranges,
                    subject_masked,
                    chunk,
                    diag_offset,
                    x_lv_ref,
                    s_idx,
                    s_id,
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
                (),
            );
        }
        if scan.two_stage_lookup_ref.is_none() {
            // Original lookup method (for non-two-stage lookup)
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
            // Dispatch point of LOSAT_X_PAIRPARNA: `x_pair_seed::run_na` does what
            // `na_scan::scan_na_lookup` does (scan, diagonal hash, ungapped extension;
            // one-hit diagonal-hash scope only) with the subject scan in strips on the pool
            // threads and the pairs folded in the scan order; `continue` ends the
            // one-iteration strand loop. Otherwise `na_scan::scan_na_lookup` runs.
            // `x_pair_seed::wants_na` asks `run_na`'s tests first, so a chunk that `run_na`
            // would send back builds no inputs and makes no call.
            if x_pair_seed::wants_na(
                use_array_indexing,
                scan.window_size,
                scan.scan_range,
                scan.diagnostics,
                !batch.use_parallel && reuse_prelim_hits,
            ) && x_pair_seed::run_na(
                &x_pair_seed::NaInputs {
                    pv_direct_lookup: scan.pv_direct_lookup_ref,
                    na_lookup: scan.na_lookup_ref,
                    use_direct_lookup: scan.use_direct_lookup,
                    word_length: scan.safe_k,
                    search_seq_packed,
                    s_len,
                    subject_seq_ranges,
                    subject_masked,
                    scan_step: scan.scan_step,
                    query_context_index: scan.query_context_index,
                    query_contexts: scan.query_contexts,
                    encoded_query_concat_blastna: scan.encoded_query_concat_blastna,
                    cutoff_scores: scan.cutoff_scores,
                    x_dropoff_scores: scan.x_dropoff_scores,
                    reduced_cutoff_scores: scan.reduced_cutoff_scores,
                    score_matrix: scan.score_matrix,
                    nucl_score_table: scan.nucl_score_table,
                    exact_ungapped: scan.exact_ungapped,
                    diag_offset,
                    diag_hash_window: diag_hash_insert_window(
                        scan.window_size,
                        scan.scan_range,
                        scan.safe_k,
                    ),
                    use_array_indexing,
                    window_size: scan.window_size,
                    scan_range: scan.scan_range,
                    diagnostics: scan.diagnostics,
                    one_subject_owner: !batch.use_parallel && reuse_prelim_hits,
                    timing: scan.timing_ref,
                },
                diag_hash,
                ungapped_hits,
                &mut subject_scratch.x_pair_fold,
            ) {
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
            // The reference scan of this chunk and the steps of the seeds it finds
            // (`na_scan::scan_na_lookup`), with what they read and write, built once per
            // chunk.
            na_scan::scan_na_lookup(
                &na_scan::NaScanInputs {
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-500,560-563
    // The chunk's HSP list from its initial hits (`chunk_hsp_list`), whose count sizes the next
    // batch; then the diagonal table's offset moves on after the chunk, also after a chunk
    // without an initial hit (`Blast_ExtendWordExit` runs inside the word finder, before
    // NCBI's `continue`).
    let chunk_prelim_hits = chunk_hsp_list::<UNGAPPED>(
        batch,
        subject,
        chunk,
        s_seq_packed,
        gap_scratch,
        subject_scratch,
        reuse_prelim_hits,
        batch.good_init_extends,
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
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:159-176
    advance_diag_table_offset(
        &mut subject_scratch.diag_table_offset,
        scan.diag_window,
        s_len,
        use_array_indexing,
        &mut subject_scratch.hit_level_array,
        &mut subject_scratch.hit_len_array,
        &mut subject_scratch.diag_hash,
    );
    chunk_prelim_hits
}

/// The HSP list of a chunk from the initial hits that its word finder saved
/// (`subject_scratch.ungapped_hits`), counted into `good_init_extends` (they size the next
/// batch). A chunk without an initial hit ends with an empty list, as NCBI's chunk loop goes on to
/// the next chunk after its word finder: the saved hits counted are none, and nothing else runs
/// for the chunk: no sort, no gapped stage and its interval tree, no HSP list, and no merge into
/// the subject's list or best-hit filter (`collect_subject_prelim` skips both for an empty list,
/// as NCBI's `continue` does). The empty list returned is the scratch's own buffer when the owner
/// reuses it (`collect_subject_prelim` puts it back), as the gapped stage's empty list was.
/// LOSAT's `-verbose` line of the gapped stage (no NCBI counterpart) is still written for the
/// chunk, with no extension. An ungapped (`-ungapped`) search makes the chunk's HSP list from its
/// initial hits (`prelim_hsp_list::ungapped_hsp_list`); a gapped one sorts them
/// (`prelim_hsp_list::sort_ungapped_hits_for_gapped`) for the gapped loop with its interval tree,
/// purge and sort (`prelim_hsp_list::gapped_hsp_list`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-500,522-527,555,560-563
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) fn chunk_hsp_list<const UNGAPPED: bool>(
    batch: &PrelimBatch<'_>,
    subject: &PrelimSubject<'_>,
    chunk: &SubjectChunk,
    s_seq_packed: &[u8],
    gap_scratch: &mut GapAlignScratch,
    subject_scratch: &mut SubjectScratch,
    reuse_prelim_hits: bool,
    good_init_extends: &AtomicU64,
) -> Vec<PrelimHit> {
    let &PrelimSubject {
        s_idx,
        s_id,
        x_prelim_subject,
        ..
    } = subject;
    let s_len = chunk.length;
    let ungapped_hits = &mut subject_scratch.ungapped_hits;
    let prelim_hits = &mut subject_scratch.prelim_hits;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1686-1689
    if ungapped_hits.is_empty() {
        if !UNGAPPED && batch.verbose {
            prelim_hsp_list::verbose_gapped_line(std::time::Duration::ZERO, 0, 0, 0);
        }
        return if reuse_prelim_hits {
            std::mem::take(prelim_hits)
        } else {
            Vec::new()
        };
    }
    if UNGAPPED {
        return prelim_hsp_list::ungapped_hsp_list(ungapped_hits, chunk, good_init_extends);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1688-1692
    // The initial hits sorted by score for the gapped loop, or grouped by query for its
    // query group tasks (LOSAT_X_QGROUPPAR), and counted
    // (`prelim_hsp_list::sort_ungapped_hits_for_gapped`).
    let (x_counts, x_order, total_ungapped) = prelim_hsp_list::sort_ungapped_hits_for_gapped(
        ungapped_hits,
        batch.x_qgroup_batch,
        reuse_prelim_hits,
        batch.queries,
        good_init_extends,
    );
    // No NCBI counterpart: the clock of the gapped stage, read only by LOSAT_TIMING's
    // `gapped_ns` and the `-verbose` line (`prelim_hsp_list::gapped_hsp_list`), so it is read
    // only when one of them is on (NCBI's engine reads no clock; the scan's clock is kept the same
    // way, `two_stage_scan::scan_two_stage_lookup`).
    let gapped_start_time =
        (batch.timing_ref.is_some() || batch.verbose).then(std::time::Instant::now);

    // The values the gapped preliminary loop of this chunk
    // (`prelim_hsp_list::gapped_score_list`) and one hit's extension
    // (`prelim_hsp_list::x_prelim_extend`) read, for the query group tasks and the
    // look-ahead helpers on the pool as for this thread.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3664-3674
    // No NCBI counterpart: the look-ahead's schedule (scheduling only).
    let gapped = prelim_hsp_list::GappedPrelim {
        batch,
        s_idx,
        s_id,
        x_prelim_subject,
        chunk,
        s_seq_packed,
        s_len,
        reuse_prelim_hits,
        total_ungapped,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:522-526,545,555
    // The chunk's gapped HSP list (`prelim_hsp_list::gapped_hsp_list`): the gapped loop per
    // query on the pool (LOSAT_X_QGROUPPAR) or over the chunk's list on this thread with
    // the subject's scratch, then the purge, the sort and the chunk's offset.
    prelim_hsp_list::gapped_hsp_list(
        &gapped,
        &prelim_hsp_list::GappedDrive {
            x_counts: &x_counts,
            x_order: &x_order,
            gapped_start_time,
        },
        ungapped_hits,
        prelim_hsp_list::GappedListState {
            interval_tree: &mut subject_scratch.x_tree,
            gap_scratch,
            greedy_align_scratch: &mut subject_scratch.greedy_align_scratch,
            x_helper_scratch: Some(&mut subject_scratch.x_ahead_scratch),
            prelim_hits,
        },
    )
}

/// Whether a job of the pool's chunk batches that holds `letters` letters takes `chunk` as well:
/// a job is a run of consecutive chunks of a subject of at most `MAX_DBSEQ_LEN` letters in all
/// (at least one chunk), searched in order with one scratch (`collect_subject_prelim`). The job
/// size is NCBI's chunk length: every chunk of an unmasked or soft-masked subject but its last
/// has `MAX_DBSEQ_LEN` letters, so each is a job of its own, as before; the short chunks of a
/// hard-masked subject (one per hard range, series NP, M4a-2b) go together, so that a subject of
/// a few thousand letters is one job on its owner's scratch instead of a batch on the pool,
/// whose fresh scratches per task cost more than the chunks (review F4 of M4a-2b: 78 % of the
/// instructions of a 4-thread run of 4,000 hard-masked subjects of 1,500 letters went to
/// `SubjectScratch::new`). No NCBI counterpart (NCBI's threads take whole subjects).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn job_takes(letters: usize, chunk: &SubjectChunk) -> bool {
    letters + chunk.length <= MAX_DBSEQ_LEN
}

/// The next batch of a subject's chunks for the pool (`collect_subject_prelim`, `HARD`
/// instance) and its partition into jobs: up to `max_jobs` jobs (`job_takes`), from `kept` (the
/// chunk that would have started one more job of the batch before) and then the chunks that
/// `next_chunk` cuts (`None`: the subject's chunks are done), into `chunks` in order, with the end
/// of each job in `chunks` into `job_ends` (both replaced). The chunk that would start job
/// `max_jobs + 1` goes to `kept` for the next batch. Returns whether `next_chunk` ended (no batch
/// follows). A batch of one job is searched serially on the subject's scratch, its jobs on the
/// pool otherwise ([`batch_jobs`]).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn next_chunk_batch(
    chunks: &mut Vec<SubjectChunk>,
    job_ends: &mut Vec<usize>,
    kept: &mut Option<SubjectChunk>,
    max_jobs: usize,
    mut next_chunk: impl FnMut() -> Option<SubjectChunk>,
) -> bool {
    chunks.clear();
    job_ends.clear();
    let mut job_letters = 0usize;
    loop {
        let Some(chunk) = kept.take().or_else(&mut next_chunk) else {
            if !chunks.is_empty() {
                job_ends.push(chunks.len());
            }
            return true;
        };
        if chunks.is_empty() || !job_takes(job_letters, &chunk) {
            if !chunks.is_empty() {
                job_ends.push(chunks.len());
            }
            if job_ends.len() == max_jobs {
                *kept = Some(chunk);
                return false;
            }
            job_letters = 0;
        }
        job_letters += chunk.length;
        chunks.push(chunk);
    }
}

/// The jobs of a batch of `next_chunk_batch`: `chunks` cut at `job_ends`, in order.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn batch_jobs<'a>(chunks: &'a [SubjectChunk], job_ends: &[usize]) -> Vec<&'a [SubjectChunk]> {
    let mut start = 0;
    job_ends
        .iter()
        .map(|&end| {
            let job = &chunks[start..end];
            start = end;
            job
        })
        .collect()
}

/// The preliminary search of one subject's chunks and the merge of their lists: the chunks that
/// `split_state` cuts (`SubjectSplitState::next_chunk`, skipping a chunk without a range), each
/// searched by `prelim_chunk_hits` with the subject's scratches on this thread, or, with a pool of
/// two threads or more (`requested_parallel`), taken in batches of up to twice the pool's threads
/// jobs (`job_takes`: runs of consecutive chunks of at most `MAX_DBSEQ_LEN` letters in all;
/// `next_chunk_batch` for a hard-masked subject's): a batch of one job, or the serial chunks of
/// LOSAT_X_BNFEWSUBJ, on this thread; else the batch's
/// jobs on the pool, a scratch per worker, each job's chunks in order, their lists collected in
/// chunk order. Each chunk's list is merged into the subject's (`merge_prelim_hit_lists`, with the
/// chunk's offset and overlap), then the subject best-hit filter (`-subject_besthit`). Returns the
/// subject's list.
/// `HARD`: the chunks come from the subject's hard ranges (a database sequence of a search with
/// `-db_hard_mask`, `HardSplitState` in the worker's scratch; series NP, M4a-2b decision H1), else
/// from `split_state` (`next_subject_chunk`); an instance of its own, so the other subjects' loops
/// test nothing.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-483,578-593
/// (the chunk loop of `s_BlastSearchEngineOneContext`)
/// The chunk batches on the pool have no NCBI counterpart (NCBI's threads take whole subjects).
pub(super) fn collect_subject_prelim<const UNGAPPED: bool, const HARD: bool>(
    scan: &ScanBatch<'_>,
    subject: &PrelimSubject<'_>,
    split_state: &mut SubjectSplitState,
    gap_scratch: &mut GapAlignScratch,
    subject_scratch: &mut SubjectScratch,
) -> Vec<PrelimHit> {
    let batch = scan.batch;
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-536
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-536
            // A batch holds up to `max_in_flight` jobs. The chunks of an unmasked or soft-masked
            // subject are jobs of one chunk each (`job_takes`: every chunk but the last has
            // `MAX_DBSEQ_LEN` letters), so its instance takes `max_in_flight` chunks as they
            // come; a hard-masked subject's (`HARD`) are partitioned by `next_chunk_batch`.
            let max_in_flight = num_threads.max(1) * 2;
            let mut chunk_batch: Vec<SubjectChunk> = Vec::with_capacity(max_in_flight);
            let mut job_ends: Vec<usize> = Vec::new();
            let mut kept: Option<SubjectChunk> = None;
            let mut done = false;
            loop {
                let jobs = if HARD {
                    done = next_chunk_batch(
                        &mut chunk_batch,
                        &mut job_ends,
                        &mut kept,
                        max_in_flight,
                        || loop {
                            match next_subject_chunk::<HARD>(
                                split_state,
                                &mut subject_scratch.db_hard_split,
                                subject_masked,
                                dbseq_chunk_overlap,
                            ) {
                                SubjectChunkStatus::Done => return None,
                                SubjectChunkStatus::NoRange => continue,
                                SubjectChunkStatus::Ok(chunk) => return Some(chunk),
                            }
                        },
                    );
                    job_ends.len()
                } else {
                    chunk_batch.clear();
                    while chunk_batch.len() < max_in_flight && !done {
                        match next_subject_chunk::<HARD>(
                            split_state,
                            &mut subject_scratch.db_hard_split,
                            subject_masked,
                            dbseq_chunk_overlap,
                        ) {
                            SubjectChunkStatus::Done => done = true,
                            SubjectChunkStatus::NoRange => continue,
                            SubjectChunkStatus::Ok(chunk) => chunk_batch.push(chunk),
                        }
                    }
                    chunk_batch.len()
                };
                if chunk_batch.is_empty() {
                    break;
                }

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:264-268
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-500
                // EXPERIMENT (LOSAT_X_BNFEWSUBJ): a batch whose longest chunk is a large
                // part of its letters (`x_few_subj::serial_chunks`), which this subject's
                // owner searches by itself, searches its chunks in order as one-chunk
                // batches do (strips for each chunk), with the one-thread loop's
                // statements below.
                let x_serial_chunks = !use_parallel
                    && x_strips
                    && x_few_subj::serial_chunks(
                        chunk_batch.iter().map(|chunk| chunk.length),
                        num_threads,
                        x_kind,
                    );
                if !use_parallel && x_strips && chunk_batch.len() > 1 {
                    crate::utils::xreach::x_reach!("LOSAT_X_BNFEWSUBJ", x_serial_chunks);
                }
                crate::utils::threading::report_stage(
                    "blastn",
                    "subject_chunks",
                    jobs,
                    jobs > 1 && !x_serial_chunks,
                );
                let soft_ranges = split_state.soft_ranges.as_slice();
                if jobs == 1 || x_serial_chunks {
                    for chunk in chunk_batch.iter() {
                        let hits = prelim_chunk_hits::<UNGAPPED>(
                            scan,
                            subject,
                            chunk,
                            soft_ranges,
                            gap_scratch,
                            subject_scratch,
                            true,
                        );
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:500,565
                        // An empty list is neither merged nor followed by the best-hit filter
                        // (the merge of an empty list returns at once, and the filter of an
                        // unchanged list removes nothing).
                        let hits = if hits.is_empty() {
                            hits
                        } else {
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
                            hits
                        };
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:488-491
                        subject_scratch.prelim_hits = hits;
                    }
                } else {
                    let new_scratches = || {
                        (
                            GapAlignScratch::new(),
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:991-1041
                            SubjectScratch::new(offset_array_size, window_size),
                        )
                    };
                    let search_chunk =
                        |state: &mut (GapAlignScratch, SubjectScratch), chunk: &SubjectChunk| {
                            let (gap_scratch, subject_scratch) = state;
                            let hits = prelim_chunk_hits::<UNGAPPED>(
                                scan,
                                subject,
                                chunk,
                                soft_ranges,
                                gap_scratch,
                                subject_scratch,
                                false,
                            );
                            (chunk.offset, chunk.overlap, hits)
                        };
                    let mut chunk_results: Vec<(usize, usize, Vec<PrelimHit>)> =
                        Vec::with_capacity(chunk_batch.len());
                    if HARD {
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-500
                        // Each job's chunks one after another with one scratch, as the subject's
                        // own loop searches them.
                        let mut job_results: Vec<Vec<(usize, usize, Vec<PrelimHit>)>> =
                            Vec::with_capacity(jobs);
                        batch_jobs(&chunk_batch, &job_ends)
                            .par_iter()
                            .map_init(new_scratches, |state, job| {
                                job.iter()
                                    .map(|chunk| search_chunk(state, chunk))
                                    .collect::<Vec<_>>()
                            })
                            .collect_into_vec(&mut job_results);
                        chunk_results.extend(job_results.into_iter().flatten());
                    } else {
                        chunk_batch
                            .par_iter()
                            .map_init(new_scratches, search_chunk)
                            .collect_into_vec(&mut chunk_results);
                    }
                    for (offset, overlap, hits) in chunk_results {
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:500,565
                        // An empty list is neither merged nor followed by the best-hit filter.
                        if hits.is_empty() {
                            continue;
                        }
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-536
        loop {
            match next_subject_chunk::<HARD>(
                split_state,
                &mut subject_scratch.db_hard_split,
                subject_masked,
                dbseq_chunk_overlap,
            ) {
                SubjectChunkStatus::Done => break,
                SubjectChunkStatus::NoRange => continue,
                SubjectChunkStatus::Ok(chunk) => {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:264-268
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-500
                    crate::utils::threading::report_stage("blastn", "subject_chunks", 1, false);
                    let soft_ranges = split_state.soft_ranges.as_slice();
                    let hits = prelim_chunk_hits::<UNGAPPED>(
                        scan,
                        subject,
                        &chunk,
                        soft_ranges,
                        gap_scratch,
                        subject_scratch,
                        true,
                    );
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:500,565
                    // An empty list is neither merged nor followed by the best-hit filter
                    // (the merge of an empty list returns at once, and the filter of an
                    // unchanged list removes nothing).
                    let hits = if hits.is_empty() {
                        hits
                    } else {
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
                        hits
                    };
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:488-491
                    subject_scratch.prelim_hits = hits;
                }
            }
        }
    }
    combined_prelim_hits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1354,1429-1432
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1677
/// The shortest subject that the preliminary stage searches (the traceback skips the shorter
/// ones too). NCBI skips only an empty subject (`min_subj_seq_length` is 1); one shorter than a
/// word has no word to scan (`BlastNaWordFinder`'s scan range ends below its start), so a
/// non-negative window's scan can skip it, but:
/// - the seeds of a MegaBLAST database index (`indexed`) can run past the end of a short sequence
///   into the next one of the index's local id, so an indexed search searches it;
/// - with a negative window (`neg_window < 0`), `BlastNaWordFinder` still calls
///   `Blast_ExtendWordExit` with the subject length, moving the carried diagonal table's offset by
///   subject length plus window, which blocks later subjects' words, so any non-empty subject must
///   reach `collect_subject_prelim` -> `neg_chunk_hits` (which exits on a short chunk while
///   honoring masks and `-use_index`).
pub(super) fn shortest_subject(effective_word_size: usize, indexed_or_neg: bool) -> usize {
    if indexed_or_neg {
        1
    } else {
        effective_word_size
    }
}

/// The preliminary search of one subject, `search_subject`'s preliminary stage (`scan`: the batch's
/// values and lookup tables, which the preliminary pass passes): the subject's lowercase masks
/// (`-lcase_masking`) and the soft ranges between them, or a database sequence's soft ranges from
/// its stored masks (`-db_soft_mask`, `read_db_soft_ranges`: a read error ends the sequence's
/// search); the early return of a subject shorter than a word; its packed letters and their copy
/// one per byte for the greedy aligner, and the split state of its chunks; the debug line of its
/// cutoffs; then its chunks searched and their lists merged (`collect_subject_prelim`), the soft
/// ranges' buffer given back to the scratch, and the subject's list reaped by the preliminary
/// e-value or by raw score (`reap_prelim`) or, in an ungapped search, finished
/// (`finish_ungapped_prelim`). Returns the list that the preliminary stage gives to the collector
/// (empty after an early return). `subject_ctx` is what `search_subject` read of the call after
/// its early return without a valid context.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1469-1475,1480-1482,1541-1542
/// (one subject of the loop of `BLAST_PreliminarySearchEngine`)
/// `s_BlastSearchEngineCore` (:701-918) searches the subject's chunks
/// (`s_BlastSearchEngineOneContext`), then sets the e-values of the subject's list and reaps it
/// (:870-900).
pub(super) fn prelim_subject<const UNGAPPED: bool>(
    scan: &ScanBatch<'_>,
    subject_ctx: &SubjectCtx<'_>,
    subject: BlastnSubjectRef<'_>,
    gap_scratch: &mut GapAlignScratch,
    subject_scratch: &mut SubjectScratch,
) -> Vec<PrelimHit> {
    let batch = scan.batch;
    let &PrelimBatch {
        effective_word_size,
        use_dp,
        lcase_masking,
        db_mask,
        db_mask_error,
        blastn_debug,
        cutoff_scores,
        ..
    } = batch;
    let &SubjectCtx {
        s_idx,
        s_id,
        s_len_full,
    } = subject_ctx;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1647-1667
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:856-874
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1079-1089
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:812-826
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2550-2556
    // The lowercase masks of the subjects are those of the FASTA records of
    // `-subject`; a database sequence has none.
    let subject_masks = match (subject.letters, subject.store) {
        (Some(letters), _) if lcase_masking => collect_lowercase_masks(letters),
        (None, Some(store)) if lcase_masking => store.lowercase().to_vec(),
        _ => Vec::new(),
    };
    // Unmasked, every chunk carries its range inline (`SeqRanges::Inline`) and no one reads
    // the soft ranges, so none are made (NCBI points them at `full_range`).
    let mut subject_masked = !subject_masks.is_empty();
    let mut soft_ranges = if subject_masked {
        build_subject_seq_ranges_from_masks(subject_masks.as_slice(), s_len_full)
    } else {
        Vec::new()
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:399-405 (cited at
    // `SubjectScratch::db_masks`)
    // With `-db_soft_mask` or `-db_hard_mask` every database sequence is a masked subject
    // (its mask type is set whether or not it has masks; decision D3 of M4a-2): the
    // complement of its stored masks (`BlastnDatabase::seq_ranges_into`), read by the
    // preliminary stage only (the traceback scans nothing), is its soft ranges or, with
    // `-db_hard_mask`, its hard ranges (series NP, M4a-2b decision H1).
    // The read is out of line (`read_db_seq_ranges`, cold): a search without the options
    // pays one test of `db_mask` per sequence and keeps its registers.
    let db_masked = db_mask.is_some();
    if let (true, Some((database, algo_id, kind))) = (db_masked, db_mask) {
        match read_db_seq_ranges(
            database,
            algo_id,
            kind,
            subject.block,
            s_idx,
            s_len_full,
            subject_scratch,
            db_mask_error,
        ) {
            Some(ranges) => {
                soft_ranges = ranges;
                subject_masked = true;
            }
            None => return Vec::new(),
        }
    }

    if s_len_full
        < shortest_subject(
            effective_word_size,
            scan.neg_window < 0
                || batch
                    .mb_index
                    .is_some_and(|index| index.is_indexed(s_idx as u32)),
        )
    {
        if db_masked {
            subject_scratch.db_soft_ranges = soft_ranges;
            subject_scratch.db_hard_split.end();
        }
        return Vec::new();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
    // The packed letters of the subject set's record (`PreparedSubjects::packed`) or the
    // stored data of the database sequence (`BlastnSubjectRef::packed`).
    let s_seq_packed_full = subject.packed;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4019-4025
    // No NCBI counterpart: this subject's packed sequence widened to one base per byte for
    // the preliminary greedy alignment instead: a window around each call's ungapped hit, and
    // the whole subject once a call reads past its window (made then, shared by every chunk
    // and worker of this subject). Off unless the search uses the greedy aligner.
    let x_prelim_subject = XPrelimSubject::new(!use_dp, s_seq_packed_full, s_len_full);

    let dbseq_chunk_overlap = DBSEQ_CHUNK_OVERLAP;
    let mut split_state = SubjectSplitState::new(s_len_full, soft_ranges);

    // BLASTN debug: Log cutoff scores for this query-subject pair
    if blastn_debug {
        subject_ctx.debug_cutoffs(cutoff_scores);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-483,593-599
    // The subject's chunks searched and their lists merged (`collect_subject_prelim`), with what
    // the search of its chunks reads of the subject, built once per subject (what it reads of the
    // batch is `scan`'s, or `special_subject`'s copy's). A database sequence of a search with
    // `-db_hard_mask` takes its chunks from its hard ranges (`HardSplitState`, in the worker's
    // scratch; series NP, M4a-2b decision H1) in an instance of its own, so the chunk loops of
    // every other subject are the base's and the subject pays one test of the mask's kind (`HARD`).
    let hard = matches!(db_mask, Some((_, _, MaskKind::Hard)));
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1727-1734 (cited at
    // `IndexSeeds::is_indexed`)
    // The seeds of a sequence of an index volume with `-db_hard_mask`, placed into its hard
    // chunks before they are searched (the Owner's answer C to Q-HARDIDX).
    let hard_index = match batch.mb_index {
        Some(index) if hard && index.is_indexed(s_idx as u32) => Some(index.hard_seeds(
            s_idx as u32,
            s_len_full,
            &subject_scratch.db_hard_split.ranges,
            dbseq_chunk_overlap,
        )),
        _ => None,
    };
    let prelim_subject = PrelimSubject {
        s_idx,
        s_id,
        s_seq_packed_full,
        x_prelim_subject: &x_prelim_subject,
        subject_masked,
        dbseq_chunk_overlap,
        hard_index: hard_index.as_ref(),
    };
    let mut combined_prelim_hits = if hard {
        collect_subject_prelim::<UNGAPPED, true>(
            scan,
            &prelim_subject,
            &mut split_state,
            gap_scratch,
            subject_scratch,
        )
    } else {
        collect_subject_prelim::<UNGAPPED, false>(
            scan,
            &prelim_subject,
            &mut split_state,
            gap_scratch,
            subject_scratch,
        )
    };
    // The ranges' buffers go back to the scratch for the next database sequence.
    if db_masked {
        subject_scratch.db_soft_ranges = std::mem::take(&mut split_state.soft_ranges);
        subject_scratch.db_hard_split.end();
    }

    if !UNGAPPED && !combined_prelim_hits.is_empty() {
        reap_prelim(batch, &mut combined_prelim_hits);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1482-1541
    // An ungapped (`-ungapped`) search finishes each subject's HSP list here, in the
    // preliminary stage: the traceback only keeps what this leaves.
    if UNGAPPED && !combined_prelim_hits.is_empty() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:891-895
        // (`s_BlastSearchEngineCore`, for gapped and ungapped searches)
        // Under matrix-only scoring the list is reaped by raw score before its finish (which
        // reaps it again after the re-evaluation, blast_engine.c:1527-1535); the e-values that
        // NCBI computes before this reap are computed again in the finish. The e-value reap of
        // the core is not run for an ungapped search here.
        if let PrelimReap::RawScore(cutoff) = batch.prelim_reap {
            combined_prelim_hits.retain(|hit| hit.prelim_score >= cutoff);
        }
        if !combined_prelim_hits.is_empty() {
            finish_ungapped_prelim(batch, subject, &mut combined_prelim_hits);
        }
    }
    combined_prelim_hits
}

/// The reap of a gapped search's subject list, after its chunks are merged: each HSP's e-value
/// from its preliminary score, kept in the HSP, and the HSPs over the e-value threshold (under
/// matrix-only scoring: below the raw cutoff, `PrelimReap`) dropped, in place and in order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:885-888,897-899
/// (`s_BlastSearchEngineCore`) and :654,662-664 (`s_Blast_HSPListReapByPrelimEvalue`, :642-674)
pub(super) fn reap_prelim(batch: &PrelimBatch<'_>, combined_prelim_hits: &mut Vec<PrelimHit>) {
    let &PrelimBatch {
        query_eff_searchsp,
        search_karlin_ref,
        round_down_evalue_score,
        prelim_reap,
        ..
    } = batch;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:885-899
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:642-671
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:868,933
    // This reap happens after preliminary gapped scoring and before the
    // traceback stream is consumed, so low preliminary scores cannot survive
    // solely because traceback later finds a high-scoring alignment.
    prelim_reap.reap(combined_prelim_hits, |prelim| {
        let eff_searchsp = query_eff_searchsp
            .get(prelim.context_idx as usize)
            .copied()
            .unwrap_or(0);
        calculate_blastn_context_statistics(
            prelim.prelim_score,
            &search_karlin_ref[prelim.context_idx as usize].gapped,
            eff_searchsp,
            round_down_evalue_score,
        )
        .1
    });
}

/// The finish of an ungapped (`-ungapped`) search's subject list in the preliminary stage (the
/// traceback only keeps what it leaves): the subject again in BLASTNA without sentinel bytes (a
/// subject set's record; a database sequence has none here, and the preliminary search's caller
/// finishes its list, `x_ungapped_database_stage`), then each HSP re-evaluated on it with the
/// identity test, the score sort, and the e-values with the reap (`ungapped_subject_stage`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1492-1497,1521-1526,1534-1536
/// (`BLAST_PreliminarySearchEngine`)
/// The tab that indents line 1535 is written as spaces.
pub(super) fn finish_ungapped_prelim(
    batch: &PrelimBatch<'_>,
    subject: BlastnSubjectRef<'_>,
    combined_prelim_hits: &mut Vec<PrelimHit>,
) {
    let &PrelimBatch {
        encoded_query_concat_blastna,
        encoded_query_concat_nomask,
        score_matrix,
        percent_identity,
        search_karlin_ref,
        round_down_evalue_score,
        prelim_reap,
        query_eff_searchsp,
        cutoff_scores,
        ..
    } = batch;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2645-2668
    // The subject again, unpacked (BLASTNA, without sentinel bytes). A database
    // sequence has none here: the preliminary search's caller reads it again, for
    // the sequences that keep an HSP, and finishes the list itself
    // (`x_ungapped_database_stage`, lazily as NCBI does).
    let mut x_subject_vec: Vec<u8>;
    let s_blastna: Option<&[u8]> = match (subject.blastna, subject.store, subject.letters) {
        (Some(blastna), _, _) => Some(blastna),
        (None, Some(store), _) => {
            x_subject_vec = Vec::new();
            x_blastna_into(store, &mut x_subject_vec);
            Some(x_subject_vec.as_slice())
        }
        (None, None, Some(letters)) => {
            x_subject_vec = encode_iupac_to_blastna(letters);
            Some(x_subject_vec.as_slice())
        }
        (None, None, None) => None,
    };
    if let Some(s_blastna) = s_blastna {
        // (`percent_identity` and the cutoffs are the values the traceback
        // (`traceback_subject::traceback_subject`) reads too; the query coverage is reaped by the
        // collector's caller.)
        ungapped_subject_stage(
            combined_prelim_hits,
            s_blastna,
            encoded_query_concat_blastna,
            encoded_query_concat_nomask,
            score_matrix,
            cutoff_scores,
            percent_identity,
            |prelim| {
                let eff_searchsp = query_eff_searchsp
                    .get(prelim.context_idx as usize)
                    .copied()
                    .unwrap_or(0);
                calculate_blastn_context_statistics(
                    prelim.prelim_score,
                    &search_karlin_ref[prelim.context_idx as usize].gapped,
                    eff_searchsp,
                    round_down_evalue_score,
                )
                .1
            },
            prelim_reap,
        );
    }
}

#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[path = "prelim_subject_tests.rs"]
mod tests;
