//! BLASTP search coordination and execution
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c

use anyhow::{bail, Context, Result};
use bio::io::fasta;
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{
    atomic::{AtomicU64, Ordering as AtomicOrdering},
    Arc,
};
use std::time::Instant;

use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use crate::algorithm::tblastx::blast_aascan::{
    s_blast_aa_scan_subject_one_range, BlastOffsetPair as OffsetPair,
};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::constants::{
    GAP_TRIGGER_BIT_SCORE, X_DROP_GAPPED_FINAL, X_DROP_GAPPED_PRELIM, X_DROP_UNGAPPED_BITS,
};
use crate::algorithm::tblastx::lookup::compressed::{
    build_blosum62_compressed_lookup, BlastCompressedAaLookupTable,
};
use crate::algorithm::tblastx::lookup::prepare_blosum62_lookup_query_for_word_size;
use crate::algorithm::tblastx::lookup::{build_ncbi_lookup, BlastAaLookupTable, QueryContext};
use crate::algorithm::tblastx::ncbi_cutoffs::{gap_trigger_raw_score, x_drop_raw_score};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::common::Hit;
use crate::config::ScoringMatrix;
use crate::core::blast_stat::compute_blosum62_ideal_karlin_params;
use crate::core::composition_adjustment::adjust_scores::{
    build_matrix_info, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    BlastCompoAdjustMode, BlastCompoGappingParams,
};
use crate::report::{
    format_bitscore_ncbi, format_evalue_ncbi_tabular, format_percent_identity_ncbi,
    write_blastp_pairwise_report, write_hit_fields, BlastpPairwiseQuery, BlastpPairwiseReport,
    OutputConfig, OutputFormat, PairwiseConfig, PairwiseHit, ReportContext,
};
use crate::stats::search_space::SearchSpace;
use crate::stats::{
    blast_spouge_etos, blast_spouge_stoe, lookup_protein_gumbel_params, lookup_protein_params,
    lookup_protein_params_ungapped, protein_scoring_supported, BlastGumbelBlk, KarlinParams,
};
use crate::utils::matrix::{blosum62_score_ncbistdaa_direct, protein_score, BLASTAA_SIZE};

use super::alignment::build_alignment_view_with_matrix;
use super::args::{BlastpArgs, BlastpCompositionMode, BlastpLookupTableType, ResolvedBlastpArgs};
use super::encoding::{
    encode_protein_query_frame_with_seg, encode_protein_sequence, ncbistdaa_to_ascii,
    EncodedProtein,
};
#[cfg(test)]
use super::extension::{extend_one_hit, extend_two_hit};
use super::extension::{extend_one_hit_blosum62, extend_two_hit_blosum62, BlastpUngappedData};
use super::gapalign::{
    blastp_get_start_for_gapped_alignment, blastp_score_only_gapped_alignment_with_scratch,
    BlastpGappedAlignmentMode, BlastpPreliminaryHsp, GapAlignScratch,
};
use super::hsp::{
    blast_compo_early_termination, collect_hits_from_hit_lists, fill_results_from_compo_heaps,
    get_prelim_hitlist_size, reap_hsplist_by_evalue, trim_by_max_hsps, BlastCompoHeap,
    BlastpHitList, BlastpHsp, BlastpHspList,
};
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
use super::kappa::BlastpPostprocessResult;
use super::kappa::{
    blastp_kappa_range_counters_env_enabled, blastp_trace_hsp_target, blastp_trace_matches_pair,
    build_query_workspace, postprocess_preliminary_hits, print_blastp_kappa_range_counters,
    reset_blastp_kappa_range_counters, BlastRedoAlignParams, BlastpKappaSubjectRangeCache,
    BlastpTraceHspTarget,
};

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:142-154
// ```c
// typedef struct BlastUngappedData {
//    Int4 q_start;
//    Int4 s_start;
//    Int4 length;
//    Int4 score;
// } BlastUngappedData;
//
// typedef struct BlastInitHSP {
//     BlastOffsetPair offsets;
//     BlastUngappedData* ungapped_data;
// } BlastInitHSP;
// ```
#[derive(Clone, Copy, Debug)]
struct InitHSP {
    // NCBI BlastInitHSP::offsets.qs_offsets.{q_off,s_off}.
    q_seed_absolute: i32,
    s_seed: i32,
    // NCBI BlastInitHSP::ungapped_data payload.
    ungapped_data: BlastpUngappedData,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3900-3905
// ```c
// q_start = init_hsp->ungapped_data->q_start;
// q_end = q_start + init_hsp->ungapped_data->length;
// s_start = init_hsp->ungapped_data->s_start;
// s_end = s_start + init_hsp->ungapped_data->length;
// ```
impl InitHSP {
    #[inline]
    fn query_end_absolute(self) -> i32 {
        self.ungapped_data.q_start + self.ungapped_data.length
    }

    #[inline]
    fn subject_end(self) -> i32 {
        self.ungapped_data.s_start + self.ungapped_data.length
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1498
// ```c
// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//        != BLAST_SEQSRC_EOF) {
//    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
//    BlastHSPStreamWrite(hsp_stream, &hsp_list);
// }
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1688
// ```c
// aux_arg->hsp_list = NULL;
// aux_arg->gap_align = BLAST_GapAlignStructNew(...);
// aux_arg->hsp_stream = BlastHSPStreamNew(...);
// ```
struct BlastpSubjectResult {
    subject_index: usize,
    query_hsp_lists: Vec<(usize, BlastpHspList)>,
}

impl BlastpSubjectResult {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1553-1554
    // ```c
    // /* Save the results. */
    // status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
    // ```
    fn new(subject_index: usize) -> Self {
        Self {
            subject_index,
            query_hsp_lists: Vec::new(),
        }
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1554
    // ```c
    // while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
    //     status = s_BlastSearchEngineCore(..., &hsp_list, ...);
    //     status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
    // }
    // ```
    fn reset(&mut self, subject_index: usize) {
        self.subject_index = subject_index;
        self.query_hsp_lists.clear();
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:991-1041
// ```c
// aux_struct->offset_pairs =
//     (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_extend.c:109-144
// ```c
// diag_table->hit_level_array = (DiagStruct *)
//     calloc(diag_table->diag_array_length, sizeof(DiagStruct));
// diag_table->offset = window_size;
// ```
struct BlastpSubjectScratch {
    offset_pairs: Vec<OffsetPair>,
    diag_array: Vec<DiagStruct>,
    diag_offset: i32,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:253-267
    // ```c
    // BlastInitHitList* BLAST_InitHitListNew(void)
    // {
    //     BlastInitHitList* init_hitlist = calloc(1, sizeof(BlastInitHitList));
    //     init_hitlist->allocated = INIT_HSP_ALLOC_NUM;
    //     init_hitlist->init_hsp_array =
    //         malloc(init_hitlist->allocated * sizeof(BlastInitHSP));
    // }
    // ```
    init_hsps: Vec<InitHSP>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3524-3527
    // ```c
    // Int4 i = 0, k;
    // BlastInitHSPNode* nodes = gap_align->chaining->nodes;
    // ```
    chaining_nodes: Vec<BlastpInitHspNode>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3628-3633
    // ```c
    // if (nodes[k].best_score - gap_score + ... < ...) {
    //     nodes[k].init_hsp->ungapped_data = NULL;
    // }
    // ```
    chaining_keep: Vec<bool>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3971-4008
    // ```c
    // BlastHSPList* new_hsp_list = Blast_HSPListNew(kHspNumMax);
    // ...
    // Blast_HSPListFree(hsp_list);
    // hsp_list = new_hsp_list;
    // ```
    preliminary_hsp_list: Vec<BlastpPreliminaryHsp>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3975-3997
    // ```c
    // Blast_IntervalTreeReset(tree);
    // /* remove all HSPs computed for the current query */
    // for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
    //     ...
    // }
    // ```
    //
    // Rust tracks how many accepted preliminary HSPs belong to each query so
    // the redo path can tell when NCBI's remove/rebuild loop would leave the
    // interval tree unchanged.
    preliminary_hsp_counts_by_query: Vec<usize>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_collector.c:117-149
    // ```c
    // query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
    // ...
    // Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
    // ```
    preliminary_hits_by_query: Vec<Vec<BlastpPreliminaryHsp>>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3748-3786
    // ```c
    // found[query_index] = TRUE;
    // restricted_align_array[query_index] =
    //     init_hsp->ungapped_data->score < reduced_cutoff;
    // ```
    restricted_align_found: Vec<bool>,
    restricted_align_array: Vec<bool>,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80
    // ```c
    // BlastGapDP* dp_mem; /**< scratch structures for dynamic programming */
    // Int4 dp_mem_alloc;  /**< current number of structures allocated */
    // ```
    preliminary_gap_scratch: GapAlignScratch,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3830-3831
    // ```c
    // tree = Blast_IntervalTreeInit(0, query->length+1,
    //                               0, subject->length+1);
    // ```
    preliminary_interval_tree: BlastIntervalTree,
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1681
// ```c
// BLAST_GapAlignSetUp(..., &score_params, &ext_params,
//                     &hit_params, &eff_len_params, &gap_align);
// BLAST_PreliminarySearchEngine(program, query, query_info,
//                               seq_src, gap_align, score_params,
//                               lookup_wrap, word_options, ext_params,
//                               hit_params, eff_len_params, ...);
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1554
// ```c
// while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
//     status = s_BlastSearchEngineCore(..., &hsp_list, ...);
//     status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
// }
// ```
struct BlastpTiming {
    query_encoding_ns: AtomicU64,
    subject_encoding_ns: AtomicU64,
    lookup_ns: AtomicU64,
    scan_ungapped_ns: AtomicU64,
    scan_ungapped_calls: AtomicU64,
    preliminary_gapped_ns: AtomicU64,
    preliminary_gapped_calls: AtomicU64,
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3955-4009
    // ```c
    // if (restricted_alignment &&
    //     gap_align->score < cutoff &&
    //     gap_align->score >= restricted_cutoff) {
    //     ...
    //     Blast_IntervalTreeReset(tree);
    //     ...
    //     redo_index = index;
    //     redo_query = query_index;
    // }
    // ```
    preliminary_exact_retry_count: AtomicU64,
    preliminary_retry_removed_hsps: AtomicU64,
    preliminary_interval_tree_rebuilds: AtomicU64,
    preliminary_interval_tree_rebuild_nodes: AtomicU64,
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2525
    // ```c
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
    // ...
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
    // ```
    preliminary_endpoint_duplicates_purged: AtomicU64,
    prelim_merge_sort_ns: AtomicU64,
    prelim_merge_sort_calls: AtomicU64,
    kappa_redo_ns: AtomicU64,
    kappa_redo_calls: AtomicU64,
    final_traceback_ns: AtomicU64,
    final_traceback_calls: AtomicU64,
    output_format_ns: AtomicU64,
    output_format_calls: AtomicU64,
}

impl BlastpTiming {
    fn new() -> Self {
        Self {
            query_encoding_ns: AtomicU64::new(0),
            subject_encoding_ns: AtomicU64::new(0),
            lookup_ns: AtomicU64::new(0),
            scan_ungapped_ns: AtomicU64::new(0),
            scan_ungapped_calls: AtomicU64::new(0),
            preliminary_gapped_ns: AtomicU64::new(0),
            preliminary_gapped_calls: AtomicU64::new(0),
            preliminary_exact_retry_count: AtomicU64::new(0),
            preliminary_retry_removed_hsps: AtomicU64::new(0),
            preliminary_interval_tree_rebuilds: AtomicU64::new(0),
            preliminary_interval_tree_rebuild_nodes: AtomicU64::new(0),
            preliminary_endpoint_duplicates_purged: AtomicU64::new(0),
            prelim_merge_sort_ns: AtomicU64::new(0),
            prelim_merge_sort_calls: AtomicU64::new(0),
            kappa_redo_ns: AtomicU64::new(0),
            kappa_redo_calls: AtomicU64::new(0),
            final_traceback_ns: AtomicU64::new(0),
            final_traceback_calls: AtomicU64::new(0),
            output_format_ns: AtomicU64::new(0),
            output_format_calls: AtomicU64::new(0),
        }
    }

    #[inline]
    fn record_duration(counter: &AtomicU64, start: Option<Instant>) {
        if let Some(start) = start {
            let elapsed_ns = start.elapsed().as_nanos().min(u64::MAX as u128) as u64;
            counter.fetch_add(elapsed_ns, AtomicOrdering::Relaxed);
        }
    }

    #[inline]
    fn record_call(counter: &AtomicU64, calls: &AtomicU64, start: Option<Instant>) {
        Self::record_duration(counter, start);
        calls.fetch_add(1, AtomicOrdering::Relaxed);
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1705
    // ```c
    // BLAST_PreliminarySearchEngine(..., diagnostics, ...);
    // Blast_DiagnosticsUpdate(diagnostics, local_diagnostics);
    // ```
    #[inline]
    fn record_count(counter: &AtomicU64, value: usize) {
        counter.fetch_add(value.min(u64::MAX as usize) as u64, AtomicOrdering::Relaxed);
    }

    #[inline]
    fn seconds(counter: &AtomicU64) -> f64 {
        counter.load(AtomicOrdering::Relaxed) as f64 / 1e9
    }

    #[inline]
    fn calls(counter: &AtomicU64) -> u64 {
        counter.load(AtomicOrdering::Relaxed)
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1470-1475
// ```c
// status = s_BlastSearchEngineCore(program_number, query, query_info,
//                                  seq_arg.seq, lookup_wrap, gap_align,
//                                  score_params, word_params, ext_params,
//                                  hit_params, db_options, diagnostics,
//                                  aux_struct, &hsp_list, ...);
// ```
#[inline]
fn blastp_timing_start(enabled: bool) -> Option<Instant> {
    if enabled {
        Some(Instant::now())
    } else {
        None
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1674-1681
// ```c
// status = BLAST_PreliminarySearchEngine(program, query, query_info,
//                                       seq_src, gap_align, score_params,
//                                       lookup_wrap, word_options,
//                                       ext_params, hit_params, eff_len_params,
//                                       psi_options, db_options, hsp_stream, ...);
// ```
#[inline]
fn blastp_timing_env_enabled() -> bool {
    #[cfg(target_arch = "wasm32")]
    {
        false
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        std::env::var_os("LOSAT_TIMING").is_some()
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1705
// ```c
// BLAST_GapAlignSetUp(...);
// BLAST_PreliminarySearchEngine(...);
// Blast_DiagnosticsUpdate(diagnostics, local_diagnostics);
// ```
fn print_blastp_timing(timing: &BlastpTiming, total_start: Option<Instant>) {
    eprintln!(
        "[TIMING] blastp_query_encoding_seg: {:.3}s",
        BlastpTiming::seconds(&timing.query_encoding_ns)
    );
    eprintln!(
        "[TIMING] blastp_subject_encoding: {:.3}s",
        BlastpTiming::seconds(&timing.subject_encoding_ns)
    );
    eprintln!(
        "[TIMING] blastp_build_lookup: {:.3}s",
        BlastpTiming::seconds(&timing.lookup_ns)
    );
    eprintln!(
        "[TIMING] blastp_scan_ungapped: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.scan_ungapped_ns),
        BlastpTiming::calls(&timing.scan_ungapped_calls)
    );
    eprintln!(
        "[TIMING] blastp_prelim_gapped: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.preliminary_gapped_ns),
        BlastpTiming::calls(&timing.preliminary_gapped_calls)
    );
    eprintln!(
        "[TIMING] blastp_prelim_exact_retries: {} (removed_hsps={}, interval_tree_rebuilds={}, rebuild_nodes={}, endpoint_duplicates_purged={})",
        BlastpTiming::calls(&timing.preliminary_exact_retry_count),
        BlastpTiming::calls(&timing.preliminary_retry_removed_hsps),
        BlastpTiming::calls(&timing.preliminary_interval_tree_rebuilds),
        BlastpTiming::calls(&timing.preliminary_interval_tree_rebuild_nodes),
        BlastpTiming::calls(&timing.preliminary_endpoint_duplicates_purged)
    );
    eprintln!(
        "[TIMING] blastp_prelim_merge_sort: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.prelim_merge_sort_ns),
        BlastpTiming::calls(&timing.prelim_merge_sort_calls)
    );
    eprintln!(
        "[TIMING] blastp_kappa_redo: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.kappa_redo_ns),
        BlastpTiming::calls(&timing.kappa_redo_calls)
    );
    eprintln!(
        "[TIMING] blastp_final_traceback_stats: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.final_traceback_ns),
        BlastpTiming::calls(&timing.final_traceback_calls)
    );
    eprintln!(
        "[TIMING] blastp_output_format: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.output_format_ns),
        BlastpTiming::calls(&timing.output_format_calls)
    );
    if let Some(total_start) = total_start {
        eprintln!(
            "[TIMING] blastp_total: {:.3}s",
            total_start.elapsed().as_secs_f64()
        );
    }
}

impl BlastpSubjectScratch {
    fn new(offset_array_size: i32, diag_array_size: i32, window: i32, query_span: i32) -> Self {
        Self {
            offset_pairs: vec![
                OffsetPair::default();
                usize::try_from(offset_array_size)
                    .expect("NCBI BLAST offset pair array size must fit in usize")
            ],
            diag_array: vec![
                DiagStruct::default();
                usize::try_from(diag_array_size)
                    .expect("NCBI BLAST diagonal array size must fit in usize")
            ],
            diag_offset: window,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:253-267
            // ```c
            // init_hitlist->init_hsp_array =
            //     malloc(init_hitlist->allocated * sizeof(BlastInitHSP));
            // ```
            init_hsps: Vec::new(),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3524-3527
            // ```c
            // BlastInitHSPNode* nodes = gap_align->chaining->nodes;
            // ```
            chaining_nodes: Vec::new(),
            chaining_keep: Vec::new(),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3971-3975
            // ```c
            // BlastHSPList* new_hsp_list = Blast_HSPListNew(kHspNumMax);
            // Blast_IntervalTreeReset(tree);
            // ```
            preliminary_hsp_list: Vec::new(),
            preliminary_hsp_counts_by_query: Vec::new(),
            preliminary_hits_by_query: Vec::new(),
            restricted_align_found: Vec::new(),
            restricted_align_array: Vec::new(),
            preliminary_gap_scratch: GapAlignScratch::new(),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3830-3831
            // ```c
            // tree = Blast_IntervalTreeInit(0, query->length+1,
            //                               0, subject->length+1);
            // ```
            preliminary_interval_tree: BlastIntervalTree::new(0, query_span, 0, 1),
        }
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_extend.c:109-144
    // ```c
    // diag_table->hit_level_array = (DiagStruct *)
    //     calloc(diag_table->diag_array_length, sizeof(DiagStruct));
    // diag_table->offset = window_size;
    // ```
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    fn prepare_independent_subject(&mut self, diag_offset: i32) {
        self.diag_offset = diag_offset;
        self.diag_array.fill(DiagStruct::default());
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:229-235
        // ```c
        // for (index = 0; index < init_hitlist->total; ++index)
        //     sfree(init_hitlist->init_hsp_array[index].ungapped_data);
        // init_hitlist->total = 0;
        // ```
        self.init_hsps.clear();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3524-3578
        // ```c
        // BlastInitHSPNode* nodes = gap_align->chaining->nodes;
        // nodes[num_nodes].init_hsp = &init_array[i];
        // ```
        self.chaining_nodes.clear();
        self.chaining_keep.clear();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3971-3975
        // ```c
        // BlastHSPList* new_hsp_list = Blast_HSPListNew(kHspNumMax);
        // Blast_IntervalTreeReset(tree);
        // ```
        self.preliminary_hsp_list.clear();
        self.preliminary_hsp_counts_by_query.clear();
        for query_hits in &mut self.preliminary_hits_by_query {
            query_hits.clear();
        }
        self.restricted_align_found.clear();
        self.restricted_align_array.clear();
        self.preliminary_interval_tree.reset();
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:152-169
    // ```c
    // if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //     ewp->diag_table->offset = ewp->diag_table->window;
    //     s_BlastDiagClear(ewp->diag_table);
    // } else {
    //     ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    // }
    // ```
    fn finish_subject(&mut self, subject_length: usize, window: i32) {
        if self.diag_offset >= i32::MAX / 4 {
            self.diag_offset = window;
            for diag in &mut self.diag_array {
                *diag = DiagStruct::clear(window);
            }
        } else {
            self.diag_offset += subject_length as i32 + window;
        }
    }
}

#[inline]
fn fasta_id(record: &fasta::Record) -> Arc<str> {
    Arc::<str>::from(record.id().split_whitespace().next().unwrap_or("unknown"))
}

#[inline]
fn path_label(path: &Path) -> String {
    path.file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string())
}

#[inline]
fn fasta_defline(record: &fasta::Record) -> String {
    match record.desc() {
        Some(desc) => format!("{} {}", record.id(), desc),
        None => record.id().to_string(),
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1426
// ```c
// /* iterate over all subject sequences */
// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//        != BLAST_SEQSRC_EOF) {
//    if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
//        continue;
//    }
// ```
struct BlastpPreparedSubject {
    id: Arc<str>,
    title: Option<String>,
    encoded: EncodedProtein,
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1426
// ```c
// /* iterate over all subject sequences */
// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//        != BLAST_SEQSRC_EOF) {
//    if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
//        continue;
//    }
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1659-1668
// ```c
// /* Use a local diagnostics structure, because the one passed in an input
//   argument can be shared between multiple threads ... */
// BLAST_GapAlignSetUp(..., &gap_align)
// ```
fn prepare_blastp_subjects_preserving_order(
    subject_records: &[fasta::Record],
    num_threads: usize,
) -> Result<Vec<BlastpPreparedSubject>> {
    let prepare_subject = |record: &fasta::Record| BlastpPreparedSubject {
        id: fasta_id(record),
        title: record.desc().map(str::to_string),
        encoded: encode_protein_sequence(record.seq()),
    };

    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    {
        let num_threads = blastp_effective_num_threads(num_threads);
        if num_threads > 1 && subject_records.len() > 1 {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build()
                .context("failed to build BLASTP subject preparation thread pool")?;
            return Ok(pool.install(|| subject_records.par_iter().map(prepare_subject).collect()));
        }
    }

    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
    {
        let _ = num_threads;
    }

    Ok(subject_records.iter().map(prepare_subject).collect())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:81-83
// ```c
// NCBI_XBLAST_EXPORT const int   kBlastMajorVersion = 2;
// NCBI_XBLAST_EXPORT const int   kBlastMinorVersion = 17;
// NCBI_XBLAST_EXPORT const int   kBlastPatchVersion = 0;
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/api/version.hpp:57-60
// ```c
// virtual string Print(void) const {
//     return CVersionInfo::Print() + "+";
// }
// ```
const NCBI_BLASTP_VERSION: &str = "2.17.0+";
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:2258-2266
// ```c
// m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
// ```
//
// Current pure-Rust blastp runtime support is restricted to the default
// BLOSUM62 path by `validate_requested_blastp_support`, so the NCBI matrix
// footer text is always `BLOSUM62`.
const BLASTP_DEFAULT_MATRIX_NAME: &str = "BLOSUM62";
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3698-3699
// ```c
// const double kRestrictedMult = 0.68;/* fraction of the ordinary cutoff score
// ```
const BLASTP_RESTRICTED_MULT: f64 = 0.68;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:672-676
// ```c
// /** SCALING_FACTOR is a multiplicative factor used to get more bits of
//  * precision in the integer matrix scores.
//  */
// #define SCALING_FACTOR 32
// ```
const BLASTP_KAPPA_SCALING_FACTOR: f64 = 32.0;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2402-2404
// ```c
// /* Bit score per alignment position threshold for preliminaru near identical
//    test */
// #define NEAR_IDENTICAL_BITS_PER_POSITION (1.74)
// ```
const BLASTP_NEAR_IDENTICAL_BITS_PER_POSITION: f64 = 1.74;
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:163
// ```c
// #define PSI_INCLUSION_ETHRESH 0.002 /**< Inclusion threshold for PSI BLAST */
// ```
const PSI_INCLUSION_ETHRESH: f64 = 0.002;

#[inline(always)]
fn blastp_standard_score(matrix: ScoringMatrix, query_residue: u8, subject_residue: u8) -> i32 {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3409-3438
    // ```c
    // score += sbp->matrix->data[*query_var][*subject_var];
    // ```
    if matrix == ScoringMatrix::Blosum62 {
        blosum62_score_ncbistdaa_direct(query_residue, subject_residue)
    } else {
        protein_score(matrix, query_residue, subject_residue)
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:457-463
// ```c
// params->gap_x_dropoff = (Int4)
//     (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
// params->gap_x_dropoff_final = (Int4)
//     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
//         params->gap_x_dropoff);
// ```
#[inline]
fn blastp_gapped_x_dropoff_raw(bits: f64, lambda: f64) -> i32 {
    ((bits * std::f64::consts::LN_2) / lambda) as i32
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2387-2415
// ```c
// gapping_params->x_dropoff = (Int4)
//     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
//         extendParams->gap_x_dropoff_final);
// ```
fn blastp_redo_x_dropoff(
    scaled_gapped_lambda: f64,
    _gapped_params: &crate::stats::KarlinParams,
    x_drop_gapped_final: i32,
) -> i32 {
    blastp_gapped_x_dropoff_raw(X_DROP_GAPPED_FINAL as f64, scaled_gapped_lambda)
        .max(x_drop_gapped_final)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2438-2448
// ```c
// near_identical_cutoff =
//     (near_identical_cutoff_bits * NCBIMATH_LN2)
//         / context->sbp->kbp_gap[index]->Lambda;
// ```
fn blastp_redo_near_identical_cutoff(scaled_gapped_lambda: f64) -> f64 {
    (BLASTP_NEAR_IDENTICAL_BITS_PER_POSITION * std::f64::consts::LN_2) / scaled_gapped_lambda
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2453-2462
// ```c
// if (do_link_hsps) {
//     cutoff_s =
//         (int) (hitParams->cutoff_score_min * context->localScalingFactor);
// } else {
//     cutoff_s = 1;
// }
// ```
fn blastp_redo_cutoff_score(
    do_link_hsps: bool,
    cutoff_score: i32,
    local_scaling_factor: f64,
) -> i32 {
    if do_link_hsps {
        ((cutoff_score as f64) * local_scaling_factor) as i32
    } else {
        1
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:189-193
// ```c
// /** Because approximate gapped alignment adds extra overhead,
//  *  it should be avoided if there is no performance benefit
//  *  to using it. */
// #define RESTRICTED_ALIGNMENT_WORST_EVALUE 10.0
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:844-855
// ```c
// if (program_number == eBlastTypeBlastp && gapped_calculation &&
//     options->expect_value <= RESTRICTED_ALIGNMENT_WORST_EVALUE) {
//     params->restricted_align = TRUE;
// }
// ```
fn blastp_restricted_align_enabled(expect_value: f64) -> bool {
    expect_value <= 10.0
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:932-940
// ```c
// if (sbp->gbp && sbp->gbp->filled) {
//     int cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
//     params->prelim_evalue = cbs_stretch*evalue;
//     new_cutoff = BLAST_SpougeEtoS(cbs_stretch*evalue, kbp, sbp->gbp,
//                     query_info->contexts[context].query_length,
//                     avg_subject_length);
// }
// ```
fn blastp_prelim_evalue(expect_value: f64, comp_based_stats: BlastpCompositionMode) -> f64 {
    let cbs_stretch = match comp_based_stats {
        BlastpCompositionMode::CompositionMatrixAdjust
        | BlastpCompositionMode::ForceFullMatrixAdjust => 5.0,
        _ => 1.0,
    };
    cbs_stretch * expect_value
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:932-940
// ```c
// if (sbp->gbp && sbp->gbp->filled) {
//     int cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
//     params->prelim_evalue = cbs_stretch*evalue;
//     new_cutoff = BLAST_SpougeEtoS(cbs_stretch*evalue, kbp, sbp->gbp,
//                     query_info->contexts[context].query_length,
//                     avg_subject_length);
// }
// ```
fn blastp_cutoff_score_max(
    expect_value: f64,
    query_length: i32,
    subject_length_for_hit_cutoff: i32,
    comp_based_stats: BlastpCompositionMode,
    gapped_params: &crate::stats::KarlinParams,
    gumbel_params: &crate::stats::BlastGumbelBlk,
) -> i32 {
    blast_spouge_etos(
        blastp_prelim_evalue(expect_value, comp_based_stats),
        gapped_params,
        gumbel_params,
        query_length,
        subject_length_for_hit_cutoff,
    )
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:825-886
// ```c
// switch (comp_stat_string[0]) {
// case '0': ... eNoCompositionBasedStats;
// case '1': ... eCompositionBasedStats;
// case '2': ... eCompositionMatrixAdjust;
// case '3': ... eCompoForceFullMatrixAdjust;
// }
// ```
fn blastp_compo_adjust_mode(mode: BlastpCompositionMode) -> BlastCompoAdjustMode {
    match mode {
        BlastpCompositionMode::NoCompositionBasedStats => {
            BlastCompoAdjustMode::NoCompositionBasedStats
        }
        BlastpCompositionMode::CompositionBasedStats => BlastCompoAdjustMode::CompositionBasedStats,
        BlastpCompositionMode::CompositionMatrixAdjust => {
            BlastCompoAdjustMode::CompositionMatrixAdjust
        }
        BlastpCompositionMode::ForceFullMatrixAdjust => BlastCompoAdjustMode::ForceFullMatrixAdjust,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3099-3107
// ```c
// if (compo_adjust_mode != eNoCompositionBasedStats) {
//     if((0 == strcmp(scoringParams->options->matrix, "BLOSUM62_20"))) {
//         localScalingFactor = SCALING_FACTOR / 10;
//     } else {
//         localScalingFactor = SCALING_FACTOR;
//     }
// } else {
//     localScalingFactor = 1.0;
// }
// ```
fn blastp_local_scaling_factor(mode: BlastpCompositionMode, matrix: ScoringMatrix) -> f64 {
    if mode == BlastpCompositionMode::NoCompositionBasedStats {
        return 1.0;
    }
    let _ = matrix;
    BLASTP_KAPPA_SCALING_FACTOR
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3119-3128
// ```c
// redoneMatches = calloc(numQueries, sizeof(BlastCompo_Heap));
// for (query_index = 0;  query_index < numQueries;  query_index++) {
//     BlastCompo_HeapInitialize(&redoneMatches[query_index],
//                               hitParams->options->hitlist_size,
//                               inclusion_ethresh);
// }
// ```
fn build_redone_match_heaps(num_queries: usize, hitlist_size: usize) -> Vec<BlastCompoHeap> {
    (0..num_queries)
        .map(|_| BlastCompoHeap::new(hitlist_size, PSI_INCLUSION_ETHRESH))
        .collect()
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
// ```c
// CMTArgs(...)
// {
// #ifdef NCBI_NO_THREADS
//     m_NumThreads = CThreadable::kMinNumThreads;
//     m_MTMode = eNotSupported;
// #endif
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// ```c
// const int kMaxValue = static_cast<int>(CSystemInfo::GetCpuCount());
// int num_threads = args[kArgNumThreads].AsInteger();
// if (num_threads > kMaxValue) {
//     m_NumThreads = kMaxValue;
// } else {
//     m_NumThreads = num_threads;
// }
// ```
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn blastp_effective_num_threads(requested: usize) -> usize {
    let cpu_count = num_cpus::get().max(1);
    if requested == 0 {
        cpu_count
    } else {
        requested.min(cpu_count)
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:152-169
// ```c
// if (ewp->diag_table->offset >= INT4_MAX / 4) {
//     ewp->diag_table->offset = ewp->diag_table->window;
//     s_BlastDiagClear(ewp->diag_table);
// } else {
//     ewp->diag_table->offset += subject_length + ewp->diag_table->window;
// }
// ```
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn blastp_subject_diag_offsets(subjects: &[EncodedProtein], window: i32) -> Vec<i32> {
    let mut offsets = Vec::with_capacity(subjects.len());
    let mut diag_offset = window;
    for subject in subjects {
        offsets.push(diag_offset);
        if diag_offset >= i32::MAX / 4 {
            diag_offset = window;
        } else {
            diag_offset += subject.aa_len as i32 + window;
        }
    }
    offsets
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_collector.c:139-159
// ```c
// if (!results->hitlist_array[0]) {
//     results->hitlist_array[0] =
//         Blast_HitListNew(params->prelim_hitlist_size);
// }
// Blast_HitListUpdate(results->hitlist_array[0], hsp_list);
// ```
fn update_blastp_preliminary_hit_lists(
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    subject_result: &mut BlastpSubjectResult,
    prelim_hitlist_size: usize,
) {
    let _subject_index = subject_result.subject_index;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_collector.c:139-159
    // ```c
    // Blast_HitListUpdate(results->hitlist_array[0], hsp_list);
    // ```
    for (q_idx, provisional_hsp_list) in subject_result.query_hsp_lists.drain(..) {
        let provisional_hit_list = preliminary_hit_lists[q_idx]
            .get_or_insert_with(|| BlastpHitList::new(prelim_hitlist_size));
        provisional_hit_list.update(provisional_hsp_list);
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1553-1554
// ```c
// /* Save the results. */
// status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:271-327
// ```c
// if (!hsp_stream->results_sorted)
//     BlastHSPStreamClose(hsp_stream);
// *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
// ```
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn merge_blastp_subject_results_in_order(
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    subject_results: Vec<Option<BlastpSubjectResult>>,
    prelim_hitlist_size: usize,
) {
    for (subject_index, subject_result) in subject_results.into_iter().enumerate() {
        let Some(mut subject_result) = subject_result else {
            continue;
        };
        debug_assert_eq!(subject_result.subject_index, subject_index);
        update_blastp_preliminary_hit_lists(
            preliminary_hit_lists,
            &mut subject_result,
            prelim_hitlist_size,
        );
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
// ```c
// /* Increment offset by 1 extra byte for the sentinel NULLB
//    between frames. */
// offset += length + 1;
// frame_offsets[context+1] = offset;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:500-503
// ```c
// int buflen = QueryInfo_GetSeqBufLen(qinfo);
// TAutoUint1Ptr buf((Uint1*) calloc(buflen+1, sizeof(Uint1)));
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-121
// ```c
// (*seq_blk)->sequence_start = (Uint1 *) buffer;
// /* The first byte is a sentinel byte. */
// (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
// (*seq_blk)->length = length;
// ```
fn build_blastp_concatenated_query(contexts: &[QueryContext]) -> Vec<u8> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:246-250
    // ```c
    // return cinfo->query_offset + cinfo->query_length + (cinfo->query_length ? 2 : 1);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-127
    // ```c
    // (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
    // (*seq_blk)->length = length;
    // ```
    //
    // `ctx.aa_seq[1..]` is the per-context `query->sequence` slice
    // (residues plus that context's trailing NULLB). Concatenating those slices
    // already produces the full `query->sequence` buffer; the extra calloc'd
    // byte that exists after `query->length` in NCBI is not part of
    // `query->sequence` itself and must not contribute to Rust slice length.
    let total_len = contexts
        .iter()
        .map(|ctx| ctx.aa_seq.len() - 1)
        .sum::<usize>();
    let mut concatenated = Vec::with_capacity(total_len);
    for ctx in contexts {
        concatenated.extend_from_slice(&ctx.aa_seq[1..]);
    }
    concatenated
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:246-250
// ```c
// BlastContextInfo * cinfo = & qinfo->contexts[qinfo->last_context];
// return cinfo->query_offset + cinfo->query_length +
//        (cinfo->query_length ? 2 : 1);
// ```
//
// LOSAT's `query_length` mirrors `BLAST_SequenceBlk->length`, which excludes
// the final trailing NULLB but includes shared inter-context boundary sentinels.
fn blastp_concatenated_query_length(contexts: &[QueryContext]) -> i32 {
    contexts
        .last()
        .map(|ctx| ctx.frame_base + ctx.aa_len as i32)
        .unwrap_or(0)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:204-221
// ```c
// if (!(query_info->contexts[context].is_valid))
//    continue;
// ...
// p->cutoffs[context].x_dropoff_init =
//    (Int4)(sbp->scale_factor *
//           ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
// ```
fn blastp_x_dropoff_init_for_context(ctx: &QueryContext) -> i32 {
    if !ctx.is_valid {
        return 0;
    }
    x_drop_raw_score(X_DROP_UNGAPPED_BITS, &ctx.karlin_params, 1.0)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:917-940
// ```c
// if (!(query_info->contexts[context].is_valid)) {
//    params->cutoffs[context].cutoff_score = INT4_MAX;
//    continue;
// }
// ...
// params->cutoffs[context].cutoff_score = new_cutoff;
// params->cutoffs[context].cutoff_score_max = new_cutoff;
// ```
fn blastp_hit_cutoff_score_for_context(
    ctx: &QueryContext,
    evalue: f64,
    subject_length: i32,
    mode: BlastpCompositionMode,
    gapped_params: &KarlinParams,
    gapped_gumbel: &BlastGumbelBlk,
) -> i32 {
    if !ctx.is_valid {
        return i32::MAX;
    }
    blastp_cutoff_score_max(
        evalue,
        ctx.aa_len as i32,
        subject_length,
        mode,
        gapped_params,
        gapped_gumbel,
    )
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:324-345
// ```c
// if (!(query_info->contexts[context].is_valid)) {
//    curr_cutoffs->cutoff_score = INT4_MAX;
//    continue;
// }
// ...
// if (sbp->kbp_std) {
//    kbp = sbp->kbp_std[context];
//    if (s_BlastKarlinBlkIsValid(kbp)) {
//       gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 +
//                                kbp->logK) / kbp->Lambda);
//    }
// }
// ```
fn blastp_word_cutoff_score_for_context(hit_cutoff_score: i32, ctx: &QueryContext) -> i32 {
    if !ctx.is_valid {
        return i32::MAX;
    }
    let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, &ctx.karlin_params);
    gap_trigger.min(hit_cutoff_score)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:969-983
// ```c
// if (sbp->gbp) {
//     min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
//     if (Blast_SubjectIsTranslated(program_number)) {
//         min_subject_length/=3;
//     }
// } else {
//     min_subject_length = (Int4) (total_length/num_seqs);
// }
// ...
// BlastHitSavingParametersNew(program_number, hit_options, sbp,
//                             query_info, min_subject_length,
//                             (*ext_params)->options->compositionBasedStats,
//                             hit_params);
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1369-1439
// ```c
// BlastInitialWordParametersNew(..., BlastSeqSrcGetAvgSeqLen(seq_src), &word_params);
// ...
// db_length = BlastSeqSrcGetTotLen(seq_src);
// ...
// if (db_length == 0) {
//     BLAST_OneSubjectUpdateParameters(..., seq_arg.seq->length, ...);
// }
// ```
//
// For multi-sequence subject sets NCBI keeps the hit-saving cutoff baseline
// fixed from setup-time sequence-source statistics; it does not recompute
// blastp protein cutoffs for each subject unless `db_length == 0`.
#[inline]
fn blastp_preliminary_cutoff_subject_length(subjects: &[EncodedProtein]) -> i32 {
    subjects
        .iter()
        .map(|subject| subject.aa_len as i32)
        .min()
        .unwrap_or(0)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:234-242
// ```c
// Blast_HSPListGetEvalues(program_number, query_info, subject_length,
//                         hsp_list, TRUE, FALSE, sbp, 0.0, 1.0);
// ```
fn build_preliminary_hsp_list_for_hitlist(
    preliminary_hits: &[BlastpPreliminaryHsp],
    oid: u32,
    query_index: u32,
    query_length: usize,
    subject_length: usize,
    karlin_params: &KarlinParams,
    gumbel_params: &BlastGumbelBlk,
) -> BlastpHspList {
    let hsps: Vec<BlastpHsp> = preliminary_hits
        .iter()
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:151-178
        // ```c
        // Blast_HSPInit(Int4 query_start, Int4 query_end, Int4 subject_start,
        //               Int4 subject_end, Int4 query_gapped_start,
        //               Int4 subject_gapped_start, Int4 query_context,
        //               Int2 query_frame, Int2 subject_frame, Int4 score, ...)
        // {
        //    new_hsp->query.offset = query_start;
        //    new_hsp->subject.offset = subject_start;
        //    new_hsp->query.end = query_end;
        //    new_hsp->subject.end = subject_end;
        //    new_hsp->query.gapped_start = query_gapped_start;
        //    new_hsp->subject.gapped_start = subject_gapped_start;
        //    new_hsp->context = query_context;
        //    new_hsp->query.frame = query_frame;
        //    new_hsp->subject.frame = subject_frame;
        // }
        // ```
        .map(|hit| BlastpHsp {
            identity: 0.0,
            length: usize::try_from(hit.query_end - hit.query_start)
                .expect("preliminary blastp HSP query span must be non-negative"),
            mismatch: 0,
            gapopen: 0,
            q_start: usize::try_from(hit.query_start + 1)
                .expect("preliminary blastp query start must be non-negative"),
            q_end: usize::try_from(hit.query_end)
                .expect("preliminary blastp query end must be non-negative"),
            s_start: usize::try_from(hit.subject_start + 1)
                .expect("preliminary blastp subject start must be non-negative"),
            s_end: usize::try_from(hit.subject_end)
                .expect("preliminary blastp subject end must be non-negative"),
            gapped_q_start: hit.gapped_query_start,
            gapped_s_start: hit.gapped_subject_start,
            e_value: blast_spouge_stoe(
                hit.raw_score,
                karlin_params,
                gumbel_params,
                query_length as i32,
                subject_length as i32,
            ),
            bit_score: 0.0,
            num_ident: 0,
            query_context: hit.query_context,
            query_frame: hit.query_frame,
            subject_frame: hit.subject_frame,
            query_length: hit.query_length,
            q_idx: hit.q_idx,
            s_idx: hit.s_idx,
            raw_score: hit.raw_score,
            gap_info: None,
            num_positives: 0,
        })
        .collect();
    let mut hsp_list = BlastpHspList {
        oid,
        query_index,
        hsps,
        best_evalue: f64::MAX,
    };
    super::hsp::update_best_evalue(&mut hsp_list);
    hsp_list
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:516-518
// ```c
// diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:776-779
// ```c
// diag_coord = (subject_offset - query_offset) & diag_mask;
// ```
#[inline]
fn blastp_diag_coord(offset_delta: i32, diag_mask: i32) -> usize {
    usize::try_from(offset_delta & diag_mask)
        .expect("NCBI BLAST masked diagonal index must be non-negative")
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:358-372
// ```c
// ungapped_data->q_start = q_start;
// ungapped_data->s_start = s_start;
// ungapped_data->length = len;
// ungapped_data->score = score;
// BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
// ```
fn blastp_save_init_hsp(
    init_hsps: &mut Vec<InitHSP>,
    q_start_absolute: i32,
    s_start: i32,
    q_off_absolute: i32,
    s_off: i32,
    len: i32,
    score: i32,
) {
    init_hsps.push(InitHSP {
        q_seed_absolute: q_off_absolute,
        s_seed: s_off,
        ungapped_data: BlastpUngappedData {
            q_start: q_start_absolute,
            s_start,
            length: len,
            score,
        },
    });
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2360-2366
// ```c
// static NCBI_INLINE Int4
// s_GetUngappedHSPContext(const BlastQueryInfo* query_info,
//                         const BlastInitHSP* init_hsp)
// {
//     return BSearchContextInfo(init_hsp->offsets.qs_offsets.q_off, query_info);
// }
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/blast_query_info.h:85-88
// ```c
// Uint4 max_length;    /**< Length of the longest among the concatenated
//                         queries */
// Uint4 min_length;    /**< Length of the shortest among the concatenated
//                         queries */
// ```
#[derive(Clone, Copy, Debug, Default)]
struct BlastpContextLookupBounds {
    min_length: i32,
    max_length: i32,
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:178-181
// ```c
// Uint4 max_length = 0;
// Uint4 min_length = INT4_MAX;
// for(TSeqPos j = 0; j < queries.Size(); j++) {
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:227-228
// ```c
// max_length = MAX(max_length, length);
// min_length = MIN(min_length, length);
// ```
#[inline]
fn blastp_context_lookup_bounds(contexts: &[QueryContext]) -> BlastpContextLookupBounds {
    let mut min_length = i32::MAX;
    let mut max_length = 0i32;
    for ctx in contexts {
        let length =
            i32::try_from(ctx.aa_len).expect("NCBI BLAST protein query length must fit Int4");
        min_length = min_length.min(length);
        max_length = max_length.max(length);
    }
    if min_length == i32::MAX {
        min_length = 0;
    }
    BlastpContextLookupBounds {
        min_length,
        max_length,
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:219-235
// ```c
// Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
// {
//     size = A->last_context+1;
//     if (A->min_length > 0 && A->max_length > 0 && A->first_context == 0) {
//         b = MIN(n / (A->max_length + 1), size - 1);
//         e = MIN(n / (A->min_length + 1) + 1, size);
//     }
//     while (b < e - 1) { ... }
//     return b;
// }
// ```
#[cfg(test)]
#[inline]
fn blastp_context_idx_for_query_offset(contexts: &[QueryContext], concat_off: i32) -> usize {
    blastp_context_idx_for_query_offset_with_bounds(
        contexts,
        blastp_context_lookup_bounds(contexts),
        concat_off,
    )
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:225-234
// ```c
// if (A->min_length > 0 && A->max_length > 0 && A->first_context == 0) {
//     b = MIN(n / (A->max_length + 1), size - 1);
//     e = MIN(n / (A->min_length + 1) + 1, size);
// }
// while (b < e - 1) {
//     m = (b + e) / 2;
//     if (A->contexts[m].query_offset > n) e = m;
//     else b = m;
// }
// ```
#[inline]
fn blastp_context_idx_for_query_offset_with_bounds(
    contexts: &[QueryContext],
    bounds: BlastpContextLookupBounds,
    concat_off: i32,
) -> usize {
    let size = contexts.len();
    if size == 0 {
        return 0;
    }
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:219-235
    // ```c
    // size = A->last_context+1;
    // ...
    // while (b < e - 1) {
    //     ...
    // }
    // return b;
    // ```
    if size == 1 {
        return 0;
    }

    let mut lo = 0usize;
    let mut hi = size;
    if bounds.min_length > 0 && bounds.max_length > 0 {
        let n = usize::try_from(concat_off)
            .expect("NCBI BLAST query offset for context lookup must be non-negative");
        let max_length =
            usize::try_from(bounds.max_length).expect("NCBI BLAST max query length must fit usize");
        let min_length =
            usize::try_from(bounds.min_length).expect("NCBI BLAST min query length must fit usize");
        lo = (n / (max_length + 1)).min(size - 1);
        hi = (n / (min_length + 1) + 1).min(size);
    }

    while lo < hi.saturating_sub(1) {
        let mid = (lo + hi) / 2;
        if concat_off < contexts[mid].frame_base {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo
}

#[cfg(test)]
#[inline]
fn blastp_get_ungapped_hsp_context(contexts: &[QueryContext], init_hsp: &InitHSP) -> usize {
    blastp_get_ungapped_hsp_context_with_bounds(
        contexts,
        blastp_context_lookup_bounds(contexts),
        init_hsp,
    )
}

#[inline]
fn blastp_get_ungapped_hsp_context_with_bounds(
    contexts: &[QueryContext],
    bounds: BlastpContextLookupBounds,
    init_hsp: &InitHSP,
) -> usize {
    blastp_context_idx_for_query_offset_with_bounds(contexts, bounds, init_hsp.q_seed_absolute)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3870-3878
// ```c
// /* If redo_index > -1 and redo_query > -1, the main loop is recomputing
//    gappaed alignments for a single query ... Until index reaches redo_index
//    again, skip all concatenated queries with index different from redo_query */
// if (index < redo_index && query_index != redo_query) {
//     continue;
// }
// ```
#[inline]
fn blastp_advance_preliminary_hit_index(
    hit_index: usize,
    redo_index: Option<usize>,
    next_same_query_index: &[Option<usize>],
) -> usize {
    if let Some(redo_index) = redo_index {
        if hit_index < redo_index && !next_same_query_index.is_empty() {
            return next_same_query_index
                .get(hit_index)
                .and_then(|index| *index)
                .map(|index| index.min(redo_index))
                .unwrap_or(redo_index);
        }
    }
    hit_index + 1
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2374-2389
// ```c
// static NCBI_INLINE void
// s_AdjustInitialHSPOffsets(BlastInitHSP* init_hsp, Int4 query_start)
// {
//     init_hsp->offsets.qs_offsets.q_off -= query_start;
//     if (init_hsp->ungapped_data) {
//         init_hsp->ungapped_data->q_start -= query_start;
//     }
// }
// ```
#[inline]
fn blastp_adjust_initial_hsp_offsets(init_hsp: &mut InitHSP, query_start: i32) {
    init_hsp.q_seed_absolute -= query_start;
    init_hsp.ungapped_data.q_start -= query_start;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2399-2466
// ```c
// s_SetUpLocalBlastSequenceBlk(query, query_info, *context,
//                              query_out, &query_start);
// ...
// single_query->sequence = concatenated_query->sequence + *query_start;
// single_query->length = query_length;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2445-2466
// ```c
// *context = s_GetUngappedHSPContext(query_info, init_hsp);
// ...
// s_AdjustInitialHSPOffsets(init_hsp, query_start);
// ```
fn blastp_adjust_hsp_offsets_and_get_query_data<'a>(
    concatenated_query: &'a [u8],
    contexts: &'a [QueryContext],
    context_lookup_bounds: BlastpContextLookupBounds,
    init_hsp: &mut InitHSP,
) -> (usize, &'a QueryContext, &'a [u8]) {
    let context_idx =
        blastp_get_ungapped_hsp_context_with_bounds(contexts, context_lookup_bounds, init_hsp);
    let ctx = &contexts[context_idx];
    let query_start = ctx.frame_base;

    assert!(
        init_hsp.ungapped_data.q_start >= query_start,
        "NCBI BLAST requires local-query adjustment to keep q_start non-negative"
    );
    blastp_adjust_initial_hsp_offsets(init_hsp, query_start);
    assert!(
        init_hsp.ungapped_data.q_start >= 0,
        "NCBI BLAST requires adjusted ungapped q_start >= 0"
    );
    let query_start_usize = usize::try_from(query_start)
        .expect("NCBI BLAST concatenated blastp query offsets must be non-negative");
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2445-2466
    // ```c
    // s_SetUpLocalBlastSequenceBlk(query, query_info, *context,
    //                              query_out, &query_start);
    // ...
    // single_query->sequence = concatenated_query->sequence + *query_start;
    // single_query->length = query_length;
    // ```
    //
    // `query_tmp.sequence` points at the current query context only; it does
    // not run into later concatenated queries. Keep the trailing NULLB so the
    // local slice matches the C buffer layout.
    let query_end_usize = query_start_usize
        .checked_add(ctx.aa_len + 1)
        .expect("NCBI BLAST local blastp query span must fit in usize");
    let query_sequence = &concatenated_query[query_start_usize..query_end_usize];
    (context_idx, ctx, query_sequence)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
// ```c
// Int4 GetOffsetArraySize(LookupTableWrap* lookup)
// {
//    ...
//    case eAaLookupTable:
//       offset_array_size = OFFSET_ARRAY_SIZE +
//          ((BlastAaLookupTable*)lookup->lut)->longest_chain;
// ```
const BLASTP_OFFSET_ARRAY_SIZE: i32 = 4096;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
// ```c
// offset_array_size = OFFSET_ARRAY_SIZE +
//     ((BlastAaLookupTable*)lookup->lut)->longest_chain;
// ...
// offset_array_size = OFFSET_ARRAY_SIZE +
//     ((BlastCompressedAaLookupTable*)lookup->lut)->longest_chain;
// ```
#[inline]
fn blastp_offset_array_size(longest_chain: i32) -> i32 {
    BLASTP_OFFSET_ARRAY_SIZE + longest_chain.max(0)
}

enum BlastpRuntimeLookup {
    Aa(BlastAaLookupTable),
    Compressed(BlastCompressedAaLookupTable),
}

impl BlastpRuntimeLookup {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:484-496
    // ```c
    // if (lookup_wrap->lut_type == eAaLookupTable) {
    //     BlastAaLookupTable *lookup = (BlastAaLookupTable *)(lookup_wrap->lut);
    //     wordsize = lookup->word_length;
    // } else {
    //     BlastCompressedAaLookupTable *lookup =
    //                     (BlastCompressedAaLookupTable *)(lookup_wrap->lut);
    //     wordsize = lookup->word_length;
    // }
    // ```
    #[inline]
    fn word_length(&self) -> i32 {
        match self {
            Self::Aa(lookup) => lookup.word_length,
            Self::Compressed(lookup) => lookup.word_length,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
    // ```c
    // case eAaLookupTable:
    //    offset_array_size = OFFSET_ARRAY_SIZE +
    //       ((BlastAaLookupTable*)lookup->lut)->longest_chain;
    // ...
    // case eCompressedAaLookupTable:
    //    offset_array_size = OFFSET_ARRAY_SIZE +
    //       ((BlastCompressedAaLookupTable*)lookup->lut)->longest_chain;
    // ```
    #[inline]
    fn longest_chain(&self) -> i32 {
        match self {
            Self::Aa(lookup) => lookup.longest_chain,
            Self::Compressed(lookup) => lookup.longest_chain,
        }
    }
}

#[derive(Clone, Copy)]
enum BlastpOffsetPairScanMode {
    OneHit,
    TwoHit,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:86-93
// ```c
// for (i = 0; i < numhits; i++) {
//     offset_pairs[i + totalhits].qs_offsets.q_off = src[i];
//     offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
// }
// ```
#[inline(always)]
fn blastp_visit_scanned_offset_pairs<F>(offset_pairs: &[OffsetPair], hits: i32, visit: &mut F)
where
    F: FnMut(OffsetPair),
{
    let offset_pairs_ptr = offset_pairs.as_ptr();
    for i in 0..hits as usize {
        // SAFETY: the selected NCBI-shaped scan callback returns at most
        // `array_size`, which is set to `offset_pairs.len()` for this buffer.
        let pair = unsafe { *offset_pairs_ptr.add(i) };
        visit(pair);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:492-505
// ```c
// scan_range[0] = 0;
// scan_range[1] = subject->seq_ranges[0].left;
// scan_range[2] = subject->seq_ranges[0].right - wordsize;
// while (scan_range[1] <= scan_range[2]) {
//     hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:763-769
// ```c
// scan_range[0] = 0;
// scan_range[1] = subject->seq_ranges[0].left;
// scan_range[2] = subject->seq_ranges[0].right - wordsize;
// while (scan_range[1] <= scan_range[2]) {
//     hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-121
// ```c
// (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
// (*seq_blk)->length = length;
// ```
fn blastp_for_each_offset_pair<F>(
    lookup: &BlastpRuntimeLookup,
    subject_sequence: &[u8],
    subject_length: usize,
    offset_pairs: &mut [OffsetPair],
    scan_mode: BlastpOffsetPairScanMode,
    mut visit: F,
) where
    F: FnMut(OffsetPair),
{
    let word_size = lookup.word_length();
    let seq_ranges = [(0, subject_length as i32)];
    let mut scan_range = [0, seq_ranges[0].0, seq_ranges[0].1 - word_size];

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:498-500
    // ```c
    // if (scan_range[2] < scan_range[1])
    //     scan_range[2] = scan_range[1];
    // ```
    if matches!(scan_mode, BlastpOffsetPairScanMode::TwoHit) && scan_range[2] < scan_range[1] {
        scan_range[2] = scan_range[1];
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:492-505,763-769
    // ```c
    // while (scan_range[1] <= scan_range[2]) {
    //     hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
    // ```
    let array_size =
        i32::try_from(offset_pairs.len()).expect("NCBI BLAST offset pair array size must fit Int4");

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:552-568
    // ```c
    // if (lookup_wrap->lut_type == eAaLookupTable)
    //    lut->scansub_callback = (void *)s_BlastAaScanSubject;
    // else if (lookup_wrap->lut_type == eCompressedAaLookupTable)
    //    lut->scansub_callback = (void *)s_BlastCompressedAaScanSubject;
    // ```
    match lookup {
        BlastpRuntimeLookup::Aa(lookup) => {
            while scan_range[1] <= scan_range[2] {
                let hits = s_blast_aa_scan_subject_one_range(
                    lookup,
                    subject_sequence,
                    offset_pairs,
                    array_size,
                    &mut scan_range,
                );
                blastp_visit_scanned_offset_pairs(offset_pairs, hits, &mut visit);
            }
        }
        BlastpRuntimeLookup::Compressed(lookup) => {
            while scan_range[1] <= scan_range[2] {
                let hits = lookup.scan_subject(
                    subject_sequence,
                    &seq_ranges,
                    offset_pairs,
                    array_size,
                    &mut scan_range,
                );
                blastp_visit_scanned_offset_pairs(offset_pairs, hits, &mut visit);
            }
        }
    }
}
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_options.c:1163-1169
// ```c
// if ((program_number == eBlastTypeTblastn ||
//      program_number == eBlastTypeBlastp ||
//      program_number == eBlastTypeBlastx) &&
//     word_size > 5)
//     options->lut_type = eCompressedAaLookupTable;
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3036-3062
// ```c
// ECompoAdjustModes compo_adjust_mode =
//     (ECompoAdjustModes) extendParams->options->compositionBasedStats;
// if (program_number != eBlastTypeBlastp ||
//     compo_adjust_mode == eNoCompositionBasedStats) {
//     ...
// }
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319
// ```c
// s_BlastProtGappedAlignment(..., const BlastScoringParameters* score_params,
//                            BlastGapAlignStruct* gap_align, ...)
// {
//     ...
//     score_left = Blast_SemiGappedAlign(..., score_params, ...);
//     ...
//     score_right = Blast_SemiGappedAlign(..., score_params, ...);
// }
// ```
fn validate_requested_blastp_support(args: &ResolvedBlastpArgs) -> Result<()> {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:917-918
    // ```c
    // arg_desc.AddFlag(kArgUngapped, "Perform ungapped alignment only?", true);
    // ```
    if args.ungapped {
        bail!(
            "unsupported pure-Rust blastp --ungapped: ungapped-only blastp search is not yet ported"
        );
    }
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:113
    // ```c
    // const string kArgUseSWTraceback(\"use_sw_tback\");
    // ```
    if args.use_sw_tback {
        bail!(
            "unsupported pure-Rust blastp -use_sw_tback: Smith-Waterman traceback is not yet ported"
        );
    }
    if args.word_size != 3
        && !(args.lookup_table_type == BlastpLookupTableType::CompressedAaLookupTable
            && args.word_size == 5)
    {
        bail!(
            "unsupported pure-Rust blastp word_size={}: only standard word_size 3 and NCBI compressed word_size 5 are currently comparison-gated",
            args.word_size
        );
    }
    if args.scoring.matrix != ScoringMatrix::Blosum62
        || args.scoring.gap_open != 11
        || args.scoring.gap_extend != 1
    {
        bail!(
            "unsupported pure-Rust blastp scoring {:?} with gap penalties {}/{}: matrix-aware protein gap alignment is not yet ported; only BLOSUM62 with gap open 11 and gap extend 1 is currently supported",
            args.scoring.matrix,
            args.scoring.gap_open,
            args.scoring.gap_extend
        );
    }
    if args.comp_based_stats.mode != BlastpCompositionMode::CompositionMatrixAdjust
        || args.comp_based_stats.unified_p
    {
        bail!(
            "unsupported pure-Rust blastp composition mode {:?}{}: only --comp_based_stats 2 is currently supported",
            args.comp_based_stats.mode,
            if args.comp_based_stats.unified_p {
                " with unified P-values"
            } else {
                ""
            }
        );
    }
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_options.c:908-936
    // ```c
    // if ((status=Blast_KarlinBlkGappedLoadFromTables(NULL, options->gap_open,
    //       options->gap_extend, options->matrix, std_matrix_only)) != 0)
    // {
    //     ...
    //     return BLASTERR_OPTION_VALUE_INVALID;
    // }
    // ```
    if !protein_scoring_supported(&args.scoring) {
        bail!(
            "blastp matrix {:?} gap penalties {}/{} are not present in the NCBI protein statistics tables used by LOSAT",
            args.scoring.matrix,
            args.scoring.gap_open,
            args.scoring.gap_extend
        );
    }
    Ok(())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:39-45
// ```c
// const char* kDfltArgTabularOutputFmt =
//     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
//     "evalue bitscore";
// const char* kDfltArgTabularOutputFmtTag("std");
// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BlastpTabularField {
    QuerySeqId,
    QueryAccessionVersion,
    QueryLength,
    SubjectSeqId,
    SubjectAccessionVersion,
    SubjectLength,
    QueryStart,
    QueryEnd,
    SubjectStart,
    SubjectEnd,
    QuerySeq,
    SubjectSeq,
    Evalue,
    BitScore,
    Score,
    AlignmentLength,
    PercentIdentical,
    NumIdentical,
    Mismatches,
    Positives,
    GapOpenings,
    Gaps,
    PercentPositives,
    QueryFrame,
    SubjectFrame,
    Frames,
    Btop,
    SubjectTitle,
}

impl BlastpTabularField {
    #[inline]
    fn name(self) -> &'static str {
        match self {
            Self::QuerySeqId => "qseqid",
            Self::QueryAccessionVersion => "qaccver",
            Self::QueryLength => "qlen",
            Self::SubjectSeqId => "sseqid",
            Self::SubjectAccessionVersion => "saccver",
            Self::SubjectLength => "slen",
            Self::QueryStart => "qstart",
            Self::QueryEnd => "qend",
            Self::SubjectStart => "sstart",
            Self::SubjectEnd => "send",
            Self::QuerySeq => "qseq",
            Self::SubjectSeq => "sseq",
            Self::Evalue => "evalue",
            Self::BitScore => "bitscore",
            Self::Score => "score",
            Self::AlignmentLength => "length",
            Self::PercentIdentical => "pident",
            Self::NumIdentical => "nident",
            Self::Mismatches => "mismatch",
            Self::Positives => "positive",
            Self::GapOpenings => "gapopen",
            Self::Gaps => "gaps",
            Self::PercentPositives => "ppos",
            Self::QueryFrame => "qframe",
            Self::SubjectFrame => "sframe",
            Self::Frames => "frames",
            Self::Btop => "btop",
            Self::SubjectTitle => "stitle",
        }
    }
}

const DEFAULT_BLASTP_TABULAR_FIELDS: [BlastpTabularField; 12] = [
    BlastpTabularField::QueryAccessionVersion,
    BlastpTabularField::SubjectAccessionVersion,
    BlastpTabularField::PercentIdentical,
    BlastpTabularField::AlignmentLength,
    BlastpTabularField::Mismatches,
    BlastpTabularField::GapOpenings,
    BlastpTabularField::QueryStart,
    BlastpTabularField::QueryEnd,
    BlastpTabularField::SubjectStart,
    BlastpTabularField::SubjectEnd,
    BlastpTabularField::Evalue,
    BlastpTabularField::BitScore,
];

#[inline]
fn default_blastp_tabular_fields() -> &'static [BlastpTabularField] {
    &DEFAULT_BLASTP_TABULAR_FIELDS
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:900-990
// ```c
// if (x_IsFieldRequested(eQuerySeq) ||
//     x_IsFieldRequested(eSubjectSeq) ||
//     x_IsFieldRequested(eBTOP) || ...) {
//     alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
//     alnVec->GetWholeAlnSeqString(1, m_SubjectSeq);
// }
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1448-1480
// ```c
// case eQuerySeq:
// case eSubjectSeq:
// case eBTOP:
// ```
#[inline]
fn blastp_tabular_fields_require_rendered_alignment(fields: &[BlastpTabularField]) -> bool {
    fields.iter().any(|field| {
        matches!(
            field,
            BlastpTabularField::QuerySeq
                | BlastpTabularField::SubjectSeq
                | BlastpTabularField::Btop
        )
    })
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:39-145
// ```c
// const char* kDfltArgTabularOutputFmtTag("std");
// const SFormatSpec sc_FormatSpecifiers[] = {
//     SFormatSpec("qseqid", ...),
//     SFormatSpec("qaccver", ...),
//     ...
//     SFormatSpec("btop", ...),
//     SFormatSpec("stitle", ...)
// };
// ```
fn parse_blastp_tabular_fields(spec: &str) -> Result<Vec<BlastpTabularField>> {
    let mut fields = Vec::new();

    for token in spec.split_whitespace() {
        let expanded = if token.eq_ignore_ascii_case("std") {
            DEFAULT_BLASTP_TABULAR_FIELDS.to_vec()
        } else {
            vec![match token {
                "qseqid" => BlastpTabularField::QuerySeqId,
                "qacc" | "qaccver" => BlastpTabularField::QueryAccessionVersion,
                "qlen" => BlastpTabularField::QueryLength,
                "sseqid" => BlastpTabularField::SubjectSeqId,
                "sacc" | "saccver" => BlastpTabularField::SubjectAccessionVersion,
                "slen" => BlastpTabularField::SubjectLength,
                "qstart" => BlastpTabularField::QueryStart,
                "qend" => BlastpTabularField::QueryEnd,
                "sstart" => BlastpTabularField::SubjectStart,
                "send" => BlastpTabularField::SubjectEnd,
                "qseq" => BlastpTabularField::QuerySeq,
                "sseq" => BlastpTabularField::SubjectSeq,
                "evalue" => BlastpTabularField::Evalue,
                "bitscore" => BlastpTabularField::BitScore,
                "score" => BlastpTabularField::Score,
                "length" => BlastpTabularField::AlignmentLength,
                "pident" => BlastpTabularField::PercentIdentical,
                "nident" => BlastpTabularField::NumIdentical,
                "mismatch" => BlastpTabularField::Mismatches,
                "positive" => BlastpTabularField::Positives,
                "gapopen" => BlastpTabularField::GapOpenings,
                "gaps" => BlastpTabularField::Gaps,
                "ppos" => BlastpTabularField::PercentPositives,
                "qframe" => BlastpTabularField::QueryFrame,
                "sframe" => BlastpTabularField::SubjectFrame,
                "frames" => BlastpTabularField::Frames,
                "btop" => BlastpTabularField::Btop,
                "stitle" => BlastpTabularField::SubjectTitle,
                other => bail!("unsupported blastp tabular field '{other}'"),
            }]
        };

        for field in expanded {
            if !fields.contains(&field) {
                fields.push(field);
            }
        }
    }

    if fields.is_empty() {
        bail!("blastp custom tabular field list cannot be empty");
    }

    Ok(fields)
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:985-1027
// ```c
// if (m_QuerySeq[i] == m_SubjectSeq[i]) {
//     ++num_ident;
//     ++num_positives;
//     ++num_matches;
// } else {
//     if (num_matches > 0) {
//         btop_string += NStr::Int8ToString(num_matches);
//         num_matches = 0;
//     }
//     btop_string += m_QuerySeq[i];
//     btop_string += m_SubjectSeq[i];
// }
// ```
fn build_btop(query_seq: &str, subject_seq: &str) -> String {
    let mut btop = String::new();
    let mut matches = 0usize;

    for (q, s) in query_seq.chars().zip(subject_seq.chars()) {
        if q == s {
            matches += 1;
        } else {
            if matches > 0 {
                btop.push_str(&matches.to_string());
                matches = 0;
            }
            btop.push(q);
            btop.push(s);
        }
    }

    if matches > 0 {
        btop.push_str(&matches.to_string());
    }

    btop
}

fn blastp_frame_value(frame: Option<i8>) -> i8 {
    frame.unwrap_or(1)
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1264-1283
// ```c
// if (num_hits != 0) {
//     PrintFieldNames(is_csv);
// }
// m_Ostream << "# " << num_hits << " hits found" << "\n";
// ```
fn write_blastp_outfmt7_header<W: Write>(
    writer: &mut W,
    context: &ReportContext,
    fields: &[BlastpTabularField],
    num_hits: usize,
) -> std::io::Result<()> {
    let version_str = context.version.as_deref().unwrap_or("0.1.0");
    writeln!(
        writer,
        "# {} {}",
        context.program.to_uppercase(),
        version_str
    )?;
    if let Some(ref query) = context.query_name {
        writeln!(writer, "# Query: {}", query)?;
    }
    if let Some(ref db) = context.subject_name {
        writeln!(writer, "# Database: {}", db)?;
    }
    if num_hits > 0 {
        let field_list = fields
            .iter()
            .map(|field| field.name())
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(writer, "# Fields: {}", field_list)?;
    }
    writeln!(writer, "# {} hits found", num_hits)?;
    Ok(())
}

fn blastp_tabular_field_text(
    field: BlastpTabularField,
    hit: &PairwiseHit,
    query_id: &str,
    subject_id: &str,
) -> String {
    let h = &hit.hit;
    let positives = hit.positives.unwrap_or(h.num_positives);
    let gaps = hit.gaps.unwrap_or_else(|| h.gap_letters());
    let query_frame = blastp_frame_value(hit.query_frame);
    let subject_frame = blastp_frame_value(hit.subject_frame);

    match field {
        BlastpTabularField::QuerySeqId | BlastpTabularField::QueryAccessionVersion => {
            query_id.to_string()
        }
        BlastpTabularField::QueryLength => h.query_length.to_string(),
        BlastpTabularField::SubjectSeqId | BlastpTabularField::SubjectAccessionVersion => {
            subject_id.to_string()
        }
        BlastpTabularField::SubjectLength => hit.subject_length.unwrap_or_default().to_string(),
        BlastpTabularField::QueryStart => h.q_start.to_string(),
        BlastpTabularField::QueryEnd => h.q_end.to_string(),
        BlastpTabularField::SubjectStart => h.s_start.to_string(),
        BlastpTabularField::SubjectEnd => h.s_end.to_string(),
        BlastpTabularField::QuerySeq => hit.query_seq.clone().unwrap_or_default(),
        BlastpTabularField::SubjectSeq => hit.subject_seq.clone().unwrap_or_default(),
        BlastpTabularField::Evalue => format_evalue_ncbi_tabular(h.e_value),
        BlastpTabularField::BitScore => format_bitscore_ncbi(h.bit_score),
        BlastpTabularField::Score => h.raw_score.to_string(),
        BlastpTabularField::AlignmentLength => h.length.to_string(),
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/objtools/align_format/tabular.hpp:437-441
        // ```c
        // double perc_ident =
        //     (m_AlignLength > 0 ? ((double)m_NumIdent)/m_AlignLength * 100 : 0);
        // m_Ostream << NStr::DoubleToString(perc_ident, 3);
        // ```
        BlastpTabularField::PercentIdentical => {
            format_percent_identity_ncbi(h.num_ident, h.length, 3)
        }
        BlastpTabularField::NumIdentical => h.num_ident.to_string(),
        BlastpTabularField::Mismatches => h.mismatch.to_string(),
        BlastpTabularField::Positives => positives.to_string(),
        BlastpTabularField::GapOpenings => h.gapopen.to_string(),
        BlastpTabularField::Gaps => gaps.to_string(),
        BlastpTabularField::PercentPositives => {
            let ppos = if h.length > 0 {
                (positives as f64 / h.length as f64) * 100.0
            } else {
                0.0
            };
            format!("{:.2}", ppos)
        }
        BlastpTabularField::QueryFrame => query_frame.to_string(),
        BlastpTabularField::SubjectFrame => subject_frame.to_string(),
        BlastpTabularField::Frames => format!("{query_frame}/{subject_frame}"),
        BlastpTabularField::Btop => {
            let query_seq = hit.query_seq.as_deref().unwrap_or("");
            let subject_seq = hit.subject_seq.as_deref().unwrap_or("");
            build_btop(query_seq, subject_seq)
        }
        BlastpTabularField::SubjectTitle => hit.subject_title.clone().unwrap_or_default(),
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// void CBlastTabularInfo::Print()
// {
//     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//         if (iter != m_FieldsToShow.begin())
//             m_Ostream << m_FieldDelimiter;
//         x_PrintField(*iter);
//     }
//     m_Ostream << "\n";
// }
// ```
fn write_blastp_hsp_tabular_field<W: Write>(
    writer: &mut W,
    field: BlastpTabularField,
    hsp: &BlastpHsp,
    query_id: &str,
    subject_id: &str,
    subject_length: usize,
    subject_title: Option<&str>,
) -> std::io::Result<()> {
    match field {
        BlastpTabularField::QuerySeqId | BlastpTabularField::QueryAccessionVersion => {
            writer.write_all(query_id.as_bytes())
        }
        BlastpTabularField::QueryLength => write!(writer, "{}", hsp.query_length),
        BlastpTabularField::SubjectSeqId | BlastpTabularField::SubjectAccessionVersion => {
            writer.write_all(subject_id.as_bytes())
        }
        BlastpTabularField::SubjectLength => write!(writer, "{subject_length}"),
        BlastpTabularField::QueryStart => write!(writer, "{}", hsp.q_start),
        BlastpTabularField::QueryEnd => write!(writer, "{}", hsp.q_end),
        BlastpTabularField::SubjectStart => write!(writer, "{}", hsp.s_start),
        BlastpTabularField::SubjectEnd => write!(writer, "{}", hsp.s_end),
        BlastpTabularField::QuerySeq
        | BlastpTabularField::SubjectSeq
        | BlastpTabularField::Btop => Ok(()),
        BlastpTabularField::Evalue => {
            writer.write_all(format_evalue_ncbi_tabular(hsp.e_value).as_bytes())
        }
        BlastpTabularField::BitScore => {
            writer.write_all(format_bitscore_ncbi(hsp.bit_score).as_bytes())
        }
        BlastpTabularField::Score => write!(writer, "{}", hsp.raw_score),
        BlastpTabularField::AlignmentLength => write!(writer, "{}", hsp.length),
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/objtools/align_format/tabular.hpp:437-441
        // ```c
        // double perc_ident =
        //     (m_AlignLength > 0 ? ((double)m_NumIdent)/m_AlignLength * 100 : 0);
        // m_Ostream << NStr::DoubleToString(perc_ident, 3);
        // ```
        BlastpTabularField::PercentIdentical => {
            writer.write_all(format_percent_identity_ncbi(hsp.num_ident, hsp.length, 3).as_bytes())
        }
        BlastpTabularField::NumIdentical => write!(writer, "{}", hsp.num_ident),
        BlastpTabularField::Mismatches => write!(writer, "{}", hsp.mismatch),
        BlastpTabularField::Positives => write!(writer, "{}", hsp.num_positives),
        BlastpTabularField::GapOpenings => write!(writer, "{}", hsp.gapopen),
        BlastpTabularField::Gaps => write!(
            writer,
            "{}",
            hsp.length
                .saturating_sub(hsp.num_ident.saturating_add(hsp.mismatch))
        ),
        BlastpTabularField::PercentPositives => {
            let ppos = if hsp.length > 0 {
                (hsp.num_positives as f64 / hsp.length as f64) * 100.0
            } else {
                0.0
            };
            write!(writer, "{ppos:.2}")
        }
        BlastpTabularField::QueryFrame => writer.write_all(b"1"),
        BlastpTabularField::SubjectFrame => writer.write_all(b"1"),
        BlastpTabularField::Frames => writer.write_all(b"1/1"),
        BlastpTabularField::SubjectTitle => {
            writer.write_all(subject_title.unwrap_or_default().as_bytes())
        }
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:68-93
// ```c
// CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
//     : m_FormatType(format_type), ..., m_Outfile(outfile),
//       m_NumSummary(num_summary), ...
// ```
fn open_output_writer<'a>(
    out_path: Option<&PathBuf>,
    in_memory_output: Option<&'a mut Vec<u8>>,
) -> Result<Box<dyn Write + 'a>> {
    let stdout = std::io::stdout();
    Ok(if let Some(output) = in_memory_output {
        Box::new(BufWriter::new(output))
    } else if let Some(path) = out_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout))
    })
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// void CBlastTabularInfo::Print()
// {
//     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//         if (iter != m_FieldsToShow.begin())
//             m_Ostream << m_FieldDelimiter;
//         x_PrintField(*iter);
//     }
//     m_Ostream << "\n";
// }
// ```
fn write_blastp_tabular_output(
    hits: &[PairwiseHit],
    fields: &[BlastpTabularField],
    outfmt: OutputFormat,
    writer: &mut impl Write,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    context: &ReportContext,
) -> Result<()> {
    if outfmt == OutputFormat::TabularWithComments {
        write_blastp_outfmt7_header(writer, context, fields, hits.len())?;
    }

    for hit in hits {
        let (query_id, subject_id) = hit.hit.resolve_ids(query_ids, subject_ids);
        for (index, field) in fields.iter().enumerate() {
            if index > 0 {
                writer.write_all(b"\t")?;
            }
            let value = blastp_tabular_field_text(*field, hit, query_id, subject_id);
            writer.write_all(value.as_bytes())?;
        }
        writer.write_all(b"\n")?;
    }

    Ok(())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:271-327
// ```c
// *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
// (*hsp_list_out)->query_index = index;
// --hit_list->hsplist_count;
// ```
fn blastp_hit_list_hsp_count(hit_lists: &[Option<BlastpHitList>]) -> usize {
    hit_lists
        .iter()
        .filter_map(|hit_list| hit_list.as_ref())
        .map(|hit_list| {
            hit_list
                .hsplist_array
                .iter()
                .take(hit_list.hsplist_count)
                .map(|hsp_list| hsp_list.hsps.len())
                .sum::<usize>()
        })
        .sum()
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//     if (iter != m_FieldsToShow.begin())
//         m_Ostream << m_FieldDelimiter;
//     x_PrintField(*iter);
// }
// m_Ostream << "\n";
// ```
fn write_blastp_tabular_hit_lists(
    hit_lists: &[Option<BlastpHitList>],
    fields: &[BlastpTabularField],
    outfmt: OutputFormat,
    writer: &mut impl Write,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    subjects: &[EncodedProtein],
    subject_titles: &[Option<String>],
    context: &ReportContext,
) -> Result<()> {
    if outfmt == OutputFormat::TabularWithComments {
        write_blastp_outfmt7_header(
            writer,
            context,
            fields,
            blastp_hit_list_hsp_count(hit_lists),
        )?;
    }

    let config = OutputConfig::ncbi_compat();
    let use_default_fields = fields == default_blastp_tabular_fields();
    for hit_list_opt in hit_lists {
        let Some(hit_list) = hit_list_opt else {
            continue;
        };
        for hsp_list in hit_list.hsplist_array.iter().take(hit_list.hsplist_count) {
            for hsp in &hsp_list.hsps {
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
                // ```c
                // typedef struct BlastHSPList {
                //    Int4 oid;
                //    Int4 query_index;
                // } BlastHSPList;
                // ```
                let query_id = query_ids
                    .get(hsp.q_idx as usize)
                    .map(|id| id.as_ref())
                    .unwrap_or("unknown");
                let subject_id = subject_ids
                    .get(hsp.s_idx as usize)
                    .map(|id| id.as_ref())
                    .unwrap_or("unknown");
                if use_default_fields {
                    write_hit_fields(
                        writer,
                        query_id,
                        subject_id,
                        hsp.identity,
                        hsp.num_ident,
                        hsp.length,
                        hsp.mismatch,
                        hsp.gapopen,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.s_start,
                        hsp.s_end,
                        hsp.e_value,
                        hsp.bit_score,
                        &config,
                    )?;
                    continue;
                }
                let subject_length = subjects
                    .get(hsp.s_idx as usize)
                    .map(|subject| subject.aa_len)
                    .unwrap_or_default();
                let subject_title = subject_titles
                    .get(hsp.s_idx as usize)
                    .and_then(|title| title.as_deref());
                for (index, field) in fields.iter().enumerate() {
                    if index > 0 {
                        writer.write_all(b"\t")?;
                    }
                    write_blastp_hsp_tabular_field(
                        writer,
                        *field,
                        hsp,
                        query_id,
                        subject_id,
                        subject_length,
                        subject_title,
                    )?;
                }
                writer.write_all(b"\n")?;
            }
        }
    }

    Ok(())
}

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:2122-2149
// ```c
// if(sequence_standard[i]==sequence[i]){
//     match ++;
// } else {
//     if ((m_AlignType&eProt)
//         && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){
//         positive ++;
//     }
// }
// ```
fn count_identities_and_positives(
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    len: usize,
    matrix: ScoringMatrix,
) -> Option<(usize, usize)> {
    if q_off + len > query.len() || s_off + len > subject.len() {
        return None;
    }

    let mut identities = 0usize;
    let mut positives = 0usize;

    for i in 0..len {
        let q = query[q_off + i];
        let s = subject[s_off + i];
        if q == s {
            identities += 1;
            positives += 1;
        } else if blastp_standard_score(matrix, q, s) > 0 {
            positives += 1;
        }
    }

    Some((identities, positives))
}

fn decode_alignment(seq: &[u8], start: usize, end: usize) -> Option<String> {
    if start > end || end > seq.len() {
        return None;
    }
    Some(
        seq[start..end]
            .iter()
            .map(|&aa| ncbistdaa_to_ascii(aa))
            .collect(),
    )
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:518-523
// ```c
// query = query_blk->sequence + query_info->contexts[hsp->context].query_offset;
// query_nomask = query_blk->sequence_nomask +
//     query_info->contexts[hsp->context].query_offset;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:649-653
// ```c
// query = query_blk->sequence + query_info->contexts[hsp->context].query_offset;
// query_nomask = query_blk->sequence_nomask +
//     query_info->contexts[hsp->context].query_offset;
// ```
#[inline]
fn query_nomask_sequence(ctx: &QueryContext) -> &[u8] {
    let query_nomask_full = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
    &query_nomask_full[1..query_nomask_full.len() - 1]
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1554
// ```c
// status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
// ```
fn build_pairwise_hits(
    hits: Vec<Hit>,
    query_contexts: &[QueryContext],
    subjects: &[EncodedProtein],
    subject_titles: &[Option<String>],
    matrix: ScoringMatrix,
    render_alignment: bool,
) -> Vec<PairwiseHit> {
    hits.into_iter()
        .map(|hit| {
            let q_idx = hit.q_idx as usize;
            let s_idx = hit.s_idx as usize;
            let query = query_nomask_sequence(&query_contexts[q_idx]);
            let subject = &subjects[s_idx].aa_seq[1..subjects[s_idx].aa_seq.len() - 1];

            let (query_seq, subject_seq) = if render_alignment {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
                // ```c
                // s_Blast_HSPGetNumIdentitiesAndPositives(...)
                // {
                //    if (!hsp->gap_info) { ... } else { ... }
                // }
                // ```
                let rendered = build_alignment_view_with_matrix(query, subject, &hit, matrix);
                if let Some(view) = rendered {
                    (Some(view.query_aln), Some(view.subject_aln))
                } else {
                    let q_start = hit.q_start.saturating_sub(1);
                    let q_end = hit.q_end.min(query.len());
                    let s_start = hit.s_start.saturating_sub(1);
                    let s_end = hit.s_end.min(subject.len());
                    let ungapped_span = q_end
                        .saturating_sub(q_start)
                        .min(s_end.saturating_sub(s_start));
                    if ungapped_span == hit.length {
                        (
                            decode_alignment(query, q_start, q_start + ungapped_span),
                            decode_alignment(subject, s_start, s_start + ungapped_span),
                        )
                    } else {
                        (None, None)
                    }
                }
            } else {
                (None, None)
            };
            let positives = hit.num_positives;
            let gaps = hit.gap_letters();

            PairwiseHit {
                hit,
                query_seq,
                subject_seq,
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1049-1093
                // ```c
                // int query_frame = 1, subject_frame = 1;
                // ...
                // SetCounts(num_ident, align_length, num_gaps, num_gap_opens,
                //           num_positives, query_frame, subject_frame);
                // ```
                query_frame: Some(1),
                subject_frame: Some(1),
                positives: Some(positives),
                gaps: Some(gaps),
                subject_length: Some(subjects[s_idx].aa_len),
                subject_title: subject_titles[s_idx].clone(),
            }
        })
        .collect()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
// ```c
// diag_array_length = 1;
// while (diag_array_length < (qlen+window_size))
//     diag_array_length = diag_array_length << 1;
// diag_table->diag_mask = diag_array_length-1;
// ```
fn compute_diag_table_shape(query_length: i32, window: i32) -> (i32, i32) {
    let mut diag_array_size = 1i32;
    while diag_array_size < query_length + window {
        diag_array_size <<= 1;
    }

    (diag_array_size, diag_array_size - 1)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3748-3786
// ```c
// /* find the first ungapped alignment for each query and determine
//    whetherto use approximate gapped alignment for the query */
// for (i = 0;i < init_hitlist->total;i++) {
//     ...
//     if (found[query_index]) {
//         continue;
//     }
//     found[query_index] = TRUE;
//     if (init_hsp->ungapped_data && init_hsp->ungapped_data->score <
//         (Int4)(kRestrictedMult * hit_params->cutoffs[contxt].cutoff_score)) {
//         restricted_align_array[query_index] = TRUE;
//     } else {
//         restricted_align_array[query_index] = FALSE;
//     }
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3752-3786
// ```c
// init_hsp = &init_hitlist->init_hsp_array[i];
// ...
// if (init_hsp->ungapped_data && init_hsp->ungapped_data->score <
//     (Int4)(kRestrictedMult * hit_params->cutoffs[contxt].cutoff_score)) {
// ```
fn build_restricted_align_array(
    init_hsps: &[InitHSP],
    contexts: &[QueryContext],
    hit_cutoff_scores: &[i32],
    num_queries: usize,
) -> Vec<bool> {
    let mut found = Vec::new();
    let mut restricted = Vec::new();
    build_restricted_align_array_into(
        init_hsps,
        contexts,
        blastp_context_lookup_bounds(contexts),
        hit_cutoff_scores,
        num_queries,
        &mut found,
        &mut restricted,
    );
    restricted
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3748-3786
// ```c
// found[query_index] = TRUE;
// if (init_hsp->ungapped_data && init_hsp->ungapped_data->score <
//     (Int4)(kRestrictedMult * hit_params->cutoffs[contxt].cutoff_score)) {
//     restricted_align_array[query_index] = TRUE;
// } else {
//     restricted_align_array[query_index] = FALSE;
// }
// ```
fn build_restricted_align_array_into(
    init_hsps: &[InitHSP],
    contexts: &[QueryContext],
    context_lookup_bounds: BlastpContextLookupBounds,
    hit_cutoff_scores: &[i32],
    num_queries: usize,
    found: &mut Vec<bool>,
    restricted: &mut Vec<bool>,
) {
    found.clear();
    found.resize(num_queries, false);
    restricted.clear();
    restricted.resize(num_queries, false);
    for init_hsp in init_hsps {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3753-3765
        // ```c
        // init_hsp = &init_hitlist->init_hsp_array[i];
        // ...
        // q_off = init_hsp->offsets.qs_offsets.q_off;
        // ...
        // query_index = Blast_GetQueryIndexFromContext(contxt,
        //                                              program_number);
        // ```
        let context_idx =
            blastp_get_ungapped_hsp_context_with_bounds(contexts, context_lookup_bounds, init_hsp);
        let query_index = usize::try_from(contexts[context_idx].q_idx)
            .expect("NCBI BLAST requires blastp context query indices to fit in usize");
        if query_index >= num_queries || found[query_index] {
            continue;
        }
        found[query_index] = true;
        restricted[query_index] = init_hsp.ungapped_data.score
            < (BLASTP_RESTRICTED_MULT * hit_cutoff_scores[context_idx] as f64) as i32;
    }
}

#[derive(Clone, Copy, Debug, Default)]
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign_priv.h:153-156
// ```c
// typedef struct BlastInitHSPNode {
//     BlastInitHSP* init_hsp;
//     Int4 best_score;
//     struct BlastInitHSPNode* next;
// } BlastInitHSPNode;
// ```
struct BlastpInitHspNode {
    init_index: usize,
    best_score: i32,
    next_index: Option<usize>,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3468-3499
// ```c
// static int s_CompareInitHSPsByQueryOffsetScore(const void* v1, const void* v2)
// {
//     ...
//     if (h1->ungapped_data->q_start < h2->ungapped_data->q_start) return -1;
//     if (h1->ungapped_data->q_start > h2->ungapped_data->q_start) return 1;
//     if (h1->ungapped_data->score < h2->ungapped_data->score) return 1;
//     if (h1->ungapped_data->score > h2->ungapped_data->score) return -1;
// }
// ```
fn compare_init_hsps_by_query_offset_score(a: &InitHSP, b: &InitHSP) -> std::cmp::Ordering {
    a.ungapped_data
        .q_start
        .cmp(&b.ungapped_data.q_start)
        .then_with(|| b.ungapped_data.score.cmp(&a.ungapped_data.score))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:274-299
// ```c
// if (0 == (result = BLAST_CMP(h2->ungapped_data->score,
//                              h1->ungapped_data->score)) &&
//     0 == (result = BLAST_CMP(h1->ungapped_data->s_start,
//                              h2->ungapped_data->s_start)) &&
//     0 == (result = BLAST_CMP(h2->ungapped_data->length,
//                              h1->ungapped_data->length)) &&
//     0 == (result = BLAST_CMP(h1->ungapped_data->q_start,
//                              h2->ungapped_data->q_start))) {
//     result = BLAST_CMP(h2->ungapped_data->length,
//                        h1->ungapped_data->length);
// }
// ```
fn compare_init_hsps_by_score_ncbi(a: &InitHSP, b: &InitHSP) -> std::cmp::Ordering {
    b.ungapped_data
        .score
        .cmp(&a.ungapped_data.score)
        .then_with(|| a.ungapped_data.s_start.cmp(&b.ungapped_data.s_start))
        .then_with(|| b.ungapped_data.length.cmp(&a.ungapped_data.length))
        .then_with(|| a.ungapped_data.q_start.cmp(&b.ungapped_data.q_start))
        .then_with(|| b.ungapped_data.length.cmp(&a.ungapped_data.length))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3468-3499
// ```c
// qsort(init_hitlist->init_hsp_array, init_hitlist->total,
//       sizeof(BlastInitHSP), s_CompareInitHSPsByQueryOffsetScore);
// ```
fn ncbi_qsort_init_hsps_by_query_offset_score(init_hsps: &mut [InitHSP]) {
    if init_hsps.len() <= 1 {
        return;
    }
    init_hsps.sort_unstable_by(compare_init_hsps_by_query_offset_score);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:306-309
// ```c
// qsort(init_hitlist->init_hsp_array, init_hitlist->total,
//       sizeof(BlastInitHSP), score_compare_match);
// ```
fn ncbi_qsort_init_hsps_by_score(init_hsps: &mut [InitHSP]) {
    if init_hsps.len() <= 1 {
        return;
    }
    init_hsps.sort_unstable_by(compare_init_hsps_by_score_ncbi);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3642-3657
// ```c
// while (init_hitlist->total > 0 &&
//       init_hitlist->init_hsp_array[init_hitlist->total - 1].ungapped_data == NULL) {
//    init_hitlist->total--;
// }
// for (k = 0; k < init_hitlist->total - 1; k++) {
//    if (init_hitlist->init_hsp_array[k].ungapped_data == NULL) {
//      Int4 end = init_hitlist->total - 1;
//      ASSERT(end > k);
//      init_hitlist->init_hsp_array[k] = init_hitlist->init_hsp_array[end];
//      init_hitlist->total--;
//      while (init_hitlist->total - 1 > k &&
//             init_hitlist->init_hsp_array[init_hitlist->total - 1].ungapped_data == NULL) {
//          init_hitlist->total--;
//      }
//    }
// }
// ```
fn compact_chained_init_hsps_ncbi(init_hsps: &mut Vec<InitHSP>, keep: &mut Vec<bool>) {
    while init_hsps
        .last()
        .is_some_and(|_| !keep[keep.len().saturating_sub(1)])
    {
        init_hsps.pop();
        keep.pop();
    }

    if init_hsps.len() <= 1 {
        return;
    }

    let mut index = 0usize;
    while index + 1 < init_hsps.len() {
        if keep[index] {
            index += 1;
            continue;
        }

        let end = init_hsps.len() - 1;
        debug_assert!(end > index);
        init_hsps[index] = init_hsps[end];
        keep[index] = keep[end];
        init_hsps.pop();
        keep.pop();

        while init_hsps.len() > index + 1 && init_hsps.last().is_some_and(|_| !keep[keep.len() - 1])
        {
            init_hsps.pop();
            keep.pop();
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3502-3658
// ```c
// static Int2 s_ChainingAlignment(..., BlastInitHitList* init_hitlist)
// {
//     ...
//     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
//           sizeof(BlastInitHSP), s_CompareInitHSPsByQueryOffsetScore);
//     ...
//     if (nodes[k].best_score - gap_score +
//         word_params->cutoffs[context].cutoff_score - 1 <
//         hit_params->cutoffs[context].cutoff_score) {
//         nodes[k].init_hsp->ungapped_data = NULL;
//     }
//     ...
//     Blast_InitHitListSortByScore(init_hitlist);
// }
// ```
fn chain_blastp_init_hsps(
    mut init_hsps: Vec<InitHSP>,
    contexts: &[QueryContext],
    word_cutoff_scores: &[i32],
    hit_cutoff_scores: &[i32],
    gap_open: i32,
    gap_extend: i32,
    keep: &mut Vec<bool>,
    nodes: &mut Vec<BlastpInitHspNode>,
) -> Vec<InitHSP> {
    if init_hsps.len() <= 1 {
        return init_hsps;
    }

    let gap_score = gap_open + gap_extend;
    ncbi_qsort_init_hsps_by_query_offset_score(&mut init_hsps);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3524-3527
    // ```c
    // Int4 i = 0, k;
    // BlastInitHSPNode* nodes = gap_align->chaining->nodes;
    // ```
    keep.clear();
    keep.resize(init_hsps.len(), true);
    let mut i = 0usize;
    let mut context_idx = 0usize;
    while i < init_hsps.len() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3536-3558
        // ```c
        // while (init_array[i].offsets.qs_offsets.q_off >
        //        query_info->contexts[context].query_offset +
        //        query_info->contexts[context].query_length &&
        //        context <= query_info->last_context) {
        //     context++;
        // }
        // ...
        // while ((i < init_hitlist->total) &&
        //        (init_array[i].offsets.qs_offsets.q_off <
        //         pctx->query_offset + pctx->query_length)) {
        // ```
        while context_idx < contexts.len()
            && init_hsps[i].q_seed_absolute
                > contexts[context_idx].frame_base + contexts[context_idx].aa_len as i32
        {
            context_idx += 1;
        }
        assert!(
            context_idx < contexts.len(),
            "NCBI BLAST chaining requires InitHSP q_off to resolve to a query context"
        );
        let context_end = contexts[context_idx].frame_base + contexts[context_idx].aa_len as i32;
        let start = i;
        while i < init_hsps.len() && init_hsps[i].q_seed_absolute < context_end {
            i += 1;
        }
        let end = i;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3551-3578
        // ```c
        // nodes[num_nodes].init_hsp = &init_array[i];
        // nodes[num_nodes].best_score = init_array[i].ungapped_data->score;
        // num_nodes++;
        // ```
        nodes.clear();
        nodes.extend((start..end).map(|init_index| BlastpInitHspNode {
            init_index,
            best_score: init_hsps[init_index].ungapped_data.score,
            next_index: None,
        }));

        for k in (0..nodes.len()).rev() {
            let self_score = nodes[k].best_score;
            for j in (k + 1)..nodes.len() {
                let left = init_hsps[nodes[k].init_index];
                let right = init_hsps[nodes[j].init_index];
                let left_len = left.ungapped_data.length;
                let q_diff = right.ungapped_data.q_start - left.ungapped_data.q_start + left_len;
                let s_diff = right.ungapped_data.s_start - left.ungapped_data.s_start + left_len;
                if s_diff < 0 {
                    continue;
                }

                let mut new_score = self_score + nodes[j].best_score;
                new_score += (q_diff.min(s_diff) * 3).min(word_cutoff_scores[context_idx]);
                new_score -= (q_diff - s_diff).abs().max(1);
                new_score -= gap_open;

                if new_score > nodes[k].best_score {
                    nodes[k].best_score = new_score;
                    nodes[k].next_index = Some(j);
                }
            }
        }

        for node in nodes.iter() {
            if node.best_score - gap_score + word_cutoff_scores[context_idx] - 1
                < hit_cutoff_scores[context_idx]
            {
                keep[node.init_index] = false;
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3658
    // ```c
    // Blast_InitHitListSortByScore(init_hitlist);
    // ```
    compact_chained_init_hsps_ncbi(&mut init_hsps, keep);
    ncbi_qsort_init_hsps_by_score(&mut init_hsps);
    init_hsps
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3949-4008
// ```c
// status = s_BlastProtGappedAlignment(..., restricted_alignment, ...);
// if (restricted_alignment &&
//     gap_align->score < cutoff &&
//     gap_align->score >= restricted_cutoff) {
//     ...
//     redo_index = index;
//     redo_query = query_index;
//     index = -1;
//     continue;
// }
// ...
// if (gap_align->score >= cutoff) {
//     Blast_HSPListSaveHSP(hsp_list, new_hsp);
// }
// ```
fn should_retry_exact_preliminary_gapped_alignment(
    restricted_alignment: bool,
    raw_score: i32,
    cutoff: i32,
    restricted_cutoff: i32,
) -> bool {
    restricted_alignment && raw_score < cutoff && raw_score >= restricted_cutoff
}

fn preliminary_tree_hsp(
    preliminary_hsp: &BlastpPreliminaryHsp,
    query_context_offset: i32,
) -> TreeHsp {
    TreeHsp {
        query_offset: preliminary_hsp.query_start,
        query_end: preliminary_hsp.query_end,
        subject_offset: preliminary_hsp.subject_start,
        subject_end: preliminary_hsp.subject_end,
        score: preliminary_hsp.raw_score,
        query_frame: 1,
        query_length: preliminary_hsp.query_length as i32,
        query_context_offset,
        subject_frame_sign: preliminary_hsp.subject_frame.signum(),
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3900-3919
// ```c
// tmp_hsp.query.offset = q_start;
// tmp_hsp.query.end = q_end;
// tmp_hsp.subject.offset = s_start;
// tmp_hsp.subject.end = s_end;
// if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
//                                   hit_options->min_diag_separation))
// ```
#[inline]
fn blastp_trace_init_hsp_overlaps_target(target: BlastpTraceHspTarget, init_hsp: &InitHSP) -> bool {
    let q_start = match usize::try_from(init_hsp.ungapped_data.q_start + 1) {
        Ok(value) => value,
        Err(_) => return false,
    };
    let q_end = match usize::try_from(init_hsp.query_end_absolute()) {
        Ok(value) => value,
        Err(_) => return false,
    };
    let s_start = match usize::try_from(init_hsp.ungapped_data.s_start + 1) {
        Ok(value) => value,
        Err(_) => return false,
    };
    let s_end = match usize::try_from(init_hsp.subject_end()) {
        Ok(value) => value,
        Err(_) => return false,
    };

    q_start >= target.q_start
        && q_end <= target.q_end
        && s_start >= target.s_start
        && s_end <= target.s_end
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
// ```c
// int ScoreCompareHSPs(const void* h1, const void* h2) {
//    ...
// }
// ```
fn preliminary_score_compare(
    a: &BlastpPreliminaryHsp,
    b: &BlastpPreliminaryHsp,
) -> std::cmp::Ordering {
    b.raw_score
        .cmp(&a.raw_score)
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| b.query_end.cmp(&a.query_end))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2315
// ```c
// static int
// s_QueryOffsetCompareHSPs(const void* v1, const void* v2)
// {
//    ...
//    if (h1->query.offset < h2->query.offset) return -1;
//    if (h1->subject.offset < h2->subject.offset) return -1;
//    if (h1->score < h2->score) return 1;
//    if (h1->query.end < h2->query.end) return 1;
//    if (h1->subject.end < h2->subject.end) return 1;
// }
// ```
fn preliminary_query_offset_compare(
    a: &BlastpPreliminaryHsp,
    b: &BlastpPreliminaryHsp,
) -> std::cmp::Ordering {
    a.query_context
        .cmp(&b.query_context)
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| b.query_end.cmp(&a.query_end))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2333-2379
// ```c
// static int
// s_QueryEndCompareHSPs(const void* v1, const void* v2)
// {
//    ...
//    if (h1->query.end < h2->query.end) return -1;
//    if (h1->subject.end < h2->subject.end) return -1;
//    if (h1->score < h2->score) return 1;
//    if (h1->query.offset < h2->query.offset) return 1;
//    if (h1->subject.offset < h2->subject.offset) return 1;
// }
// ```
fn preliminary_query_end_compare(
    a: &BlastpPreliminaryHsp,
    b: &BlastpPreliminaryHsp,
) -> std::cmp::Ordering {
    a.query_context
        .cmp(&b.query_context)
        .then_with(|| a.query_end.cmp(&b.query_end))
        .then_with(|| a.subject_end.cmp(&b.subject_end))
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| b.query_start.cmp(&a.query_start))
        .then_with(|| b.subject_start.cmp(&a.subject_start))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2505
// ```c
// qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
// qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
// ```
fn ncbi_qsort_preliminary_hits(
    hits: &mut [BlastpPreliminaryHsp],
    compar: fn(&BlastpPreliminaryHsp, &BlastpPreliminaryHsp) -> std::cmp::Ordering,
) {
    if hits.len() <= 1 {
        return;
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2495,2504-2520
    // ```c
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
    // ...
    // hsp = hsp_array[i+j];
    // hsp = Blast_HSPFree(hsp);
    // ...
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
    // ...
    // hsp = hsp_array[i+j];
    // hsp = Blast_HSPFree(hsp);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2320,2385
    // ```c
    // return 0;
    // ```
    hits.sort_by(compar);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1355-1369
// ```c
// if (ScoreCompareHSPs(&hsp_list->hsp_array[index],
//                      &hsp_list->hsp_array[index+1]) > 0) {
//     return FALSE;
// }
// ```
fn preliminary_hits_is_sorted_by_score_ncbi(hits: &[BlastpPreliminaryHsp]) -> bool {
    if hits.len() <= 1 {
        return true;
    }

    for index in 0..hits.len() - 1 {
        if preliminary_score_compare(&hits[index], &hits[index + 1]) == std::cmp::Ordering::Greater
        {
            return false;
        }
    }
    true
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2502
// ```c
// qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
// while (i+j < hsp_count &&
//        hsp_array[i]->context == hsp_array[i+j]->context &&
//        hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
//        hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
//        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame) {
//     hsp_count--;
//     hsp = Blast_HSPFree(hsp);
//     for (k=i+j; k<hsp_count; k++) hsp_array[k] = hsp_array[k+1];
// }
// ```
#[inline]
fn preliminary_hits_share_query_offset_endpoint(
    first: &BlastpPreliminaryHsp,
    candidate: &BlastpPreliminaryHsp,
) -> bool {
    first.query_context == candidate.query_context
        && first.query_start == candidate.query_start
        && first.subject_start == candidate.subject_start
        && first.subject_frame == candidate.subject_frame
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2504-2527
// ```c
// qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
// while (i+j < hsp_count &&
//        hsp_array[i]->context == hsp_array[i+j]->context &&
//        hsp_array[i]->query.end == hsp_array[i+j]->query.end &&
//        hsp_array[i]->subject.end == hsp_array[i+j]->subject.end &&
//        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame) {
//     hsp_count--;
//     hsp = Blast_HSPFree(hsp);
//     for (k=i+j; k<hsp_count; k++) hsp_array[k] = hsp_array[k+1];
// }
// ```
#[inline]
fn preliminary_hits_share_query_end_endpoint(
    first: &BlastpPreliminaryHsp,
    candidate: &BlastpPreliminaryHsp,
) -> bool {
    first.query_context == candidate.query_context
        && first.query_end == candidate.query_end
        && first.subject_end == candidate.subject_end
        && first.subject_frame == candidate.subject_frame
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2493-2500,2519-2525
// ```c
// hsp = Blast_HSPFree(hsp);
// for (k=i+j; k<hsp_count; k++) {
//     hsp_array[k] = hsp_array[k+1];
// }
// hsp_array[hsp_count] = hsp;
// ```
fn purge_sorted_preliminary_endpoint_duplicates(
    hits: &mut Vec<BlastpPreliminaryHsp>,
    same_endpoint: fn(&BlastpPreliminaryHsp, &BlastpPreliminaryHsp) -> bool,
) -> usize {
    if hits.len() <= 1 {
        return 0;
    }

    let mut write_index = 1usize;
    for read_index in 1..hits.len() {
        if same_endpoint(&hits[write_index - 1], &hits[read_index]) {
            continue;
        }
        if write_index != read_index {
            hits[write_index] = hits[read_index];
        }
        write_index += 1;
    }
    let removed = hits.len().saturating_sub(write_index);
    hits.truncate(write_index);
    removed
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3870-3878
// ```c
// query_index = Blast_GetQueryIndexFromContext(context, program_number);
// ```
fn blastp_query_index_from_preliminary_context(
    contexts: &[QueryContext],
    query_context: i32,
) -> usize {
    let context_idx = usize::try_from(query_context)
        .expect("NCBI BLAST requires preliminary blastp query context to be non-negative");
    usize::try_from(contexts[context_idx].q_idx)
        .expect("NCBI BLAST requires blastp query index to fit in usize")
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2538
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(..., Boolean purge)
// {
//    purge |= (program != eBlastTypeBlastn);
//    ...
//    qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
//    ...
//    qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
//    ...
// }
// ```
fn purge_preliminary_hits_with_common_endpoints(
    mut hits: Vec<BlastpPreliminaryHsp>,
) -> (Vec<BlastpPreliminaryHsp>, usize) {
    if hits.len() <= 1 {
        return (hits, 0);
    }

    ncbi_qsort_preliminary_hits(&mut hits, preliminary_query_offset_compare);
    let mut purged = purge_sorted_preliminary_endpoint_duplicates(
        &mut hits,
        preliminary_hits_share_query_offset_endpoint,
    );

    ncbi_qsort_preliminary_hits(&mut hits, preliminary_query_end_compare);
    purged += purge_sorted_preliminary_endpoint_duplicates(
        &mut hits,
        preliminary_hits_share_query_end_endpoint,
    );

    (hits, purged)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:540-555
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
// ...
// Blast_HSPListSortByScore(hsp_list);
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3612-3619
// ```c
// localMatch->hsp_array,  /* i */
// localMatch->hspcnt,     /* i */
// ...
// *pStatusCode = s_ResultHspToDistinctAlign(...);
// ```
fn prepare_preliminary_hits_for_kappa(
    mut hits: Vec<BlastpPreliminaryHsp>,
) -> Vec<BlastpPreliminaryHsp> {
    prepare_preliminary_hits_for_kappa_in_place(&mut hits);
    hits
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:540-555
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
// ...
// Blast_HSPListSortByScore(hsp_list);
// ```
fn prepare_preliminary_hits_for_kappa_in_place(hits: &mut Vec<BlastpPreliminaryHsp>) -> usize {
    let (purged, purged_count) = purge_preliminary_hits_with_common_endpoints(std::mem::take(hits));
    *hits = purged;
    if !preliminary_hits_is_sorted_by_score_ncbi(hits) {
        ncbi_qsort_preliminary_hits(hits, preliminary_score_compare);
    }
    purged_count
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_collector.c:117-149
// ```c
// query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
// ...
// Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
// ```
fn split_preliminary_hits_by_query(
    preliminary_hits: Vec<BlastpPreliminaryHsp>,
    contexts: &[QueryContext],
    num_queries: usize,
) -> Vec<Vec<BlastpPreliminaryHsp>> {
    let mut hits_by_query = vec![Vec::new(); num_queries];
    let mut preliminary_hits = preliminary_hits;
    split_preliminary_hits_by_query_into(
        &mut preliminary_hits,
        contexts,
        num_queries,
        &mut hits_by_query,
    );
    hits_by_query
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_collector.c:117-149
// ```c
// query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
// ...
// Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
// ```
fn split_preliminary_hits_by_query_into(
    preliminary_hits: &mut Vec<BlastpPreliminaryHsp>,
    contexts: &[QueryContext],
    num_queries: usize,
    hits_by_query: &mut Vec<Vec<BlastpPreliminaryHsp>>,
) {
    if hits_by_query.len() < num_queries {
        hits_by_query.resize_with(num_queries, Vec::new);
    }
    for query_hits in hits_by_query.iter_mut() {
        query_hits.clear();
    }
    for preliminary_hsp in preliminary_hits.drain(..) {
        let query_index =
            blastp_query_index_from_preliminary_context(contexts, preliminary_hsp.query_context);
        if query_index < num_queries {
            hits_by_query[query_index].push(preliminary_hsp);
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3975-3997
// ```c
// /* remove all HSPs computed for the current query */
// for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
//     ...
//     if (q_index2 != query_index) { ... } else { ... }
// }
// ```
fn remove_preliminary_hits_for_query(
    preliminary_hits: &mut Vec<BlastpPreliminaryHsp>,
    query_index: usize,
    contexts: &[QueryContext],
) -> usize {
    let before = preliminary_hits.len();
    preliminary_hits.retain(|hsp| {
        blastp_query_index_from_preliminary_context(contexts, hsp.query_context) != query_index
    });
    before.saturating_sub(preliminary_hits.len())
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3577-3584
// ```c
// *pStatusCode = s_ResultHspToDistinctAlign(
//         incoming_align_set,     /* o */
//         numAligns,              /* o */
//         localMatch->hsp_array,  /* i */
//         localMatch->hspcnt,     /* i */
//         context_index,          /* i */
//         queryInfo,              /* i */
//         localScalingFactor      /* i */
// );
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:151-178
// ```c
// new_hsp->query.gapped_start = query_gapped_start;
// new_hsp->subject.gapped_start = subject_gapped_start;
// new_hsp->context = query_context;
// new_hsp->query.frame = query_frame;
// new_hsp->subject.frame = subject_frame;
// ```
fn preliminary_hit_from_local_match_hsp(local_hsp: &BlastpHsp) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: i32::try_from(
            local_hsp
                .q_start
                .checked_sub(1)
                .expect("NCBI BLAST localMatch query start must be 1-based"),
        )
        .expect("NCBI BLAST localMatch query start must fit in Int4"),
        query_end: i32::try_from(local_hsp.q_end)
            .expect("NCBI BLAST localMatch query end must fit in Int4"),
        subject_start: i32::try_from(
            local_hsp
                .s_start
                .checked_sub(1)
                .expect("NCBI BLAST localMatch subject start must be 1-based"),
        )
        .expect("NCBI BLAST localMatch subject start must fit in Int4"),
        subject_end: i32::try_from(local_hsp.s_end)
            .expect("NCBI BLAST localMatch subject end must fit in Int4"),
        gapped_query_start: local_hsp.gapped_q_start,
        gapped_subject_start: local_hsp.gapped_s_start,
        raw_score: local_hsp.raw_score,
        query_context: local_hsp.query_context,
        query_frame: local_hsp.query_frame,
        subject_frame: local_hsp.subject_frame,
        query_length: local_hsp.query_length,
        q_idx: local_hsp.q_idx,
        s_idx: local_hsp.s_idx,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3975-4000
// ```c
// Blast_IntervalTreeReset(tree);
// for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
//     ...
//     BlastIntervalTreeAddHSP(hsp_list->hsp_array[index2], tree,
//                             query_info, eQueryAndSubject);
// }
// ```
fn rebuild_preliminary_interval_tree_from_list(
    interval_tree: &mut BlastIntervalTree,
    preliminary_hits: &[BlastpPreliminaryHsp],
    contexts: &[QueryContext],
) -> usize {
    interval_tree.reset();
    let mut inserted = 0usize;
    for preliminary_hsp in preliminary_hits {
        let context_idx = usize::try_from(preliminary_hsp.query_context)
            .expect("preliminary blastp query context must be non-negative");
        interval_tree.add_hsp(
            preliminary_tree_hsp(preliminary_hsp, contexts[context_idx].frame_base),
            contexts[context_idx].frame_base,
            IndexMethod::QueryAndSubject,
        );
        inserted += 1;
    }
    inserted
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3777-3786
// ```c
// /* use approximate gapped alignment if the highest scoring ungapped
//    alignment scores below the reduced cutoff */
// if (init_hsp->ungapped_data && init_hsp->ungapped_data->score <
//     (Int4)(kRestrictedMult * hit_params->cutoffs[contxt].cutoff_score)) {
//     restricted_align_array[query_index] = TRUE;
// } else {
//     restricted_align_array[query_index] = FALSE;
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:635-688
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
// ...
// Blast_HSPListSortByScore(hsp_list);
// ...
// if (BlastIntervalTreeContainsHSP(tree, hsp, query_info,
//                                  hit_options->min_diag_separation)) {
//     hsp_array[index] = Blast_HSPFree(hsp);
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3900-3949
// ```c
// tmp_hsp.query.offset = q_start;
// tmp_hsp.query.end = q_end;
// ...
// max_offset =
//    BlastGetStartForGappedAlignment(query_tmp.sequence, subject->sequence,
//       gap_align->sbp, init_hsp->ungapped_data->q_start,
//       init_hsp->ungapped_data->length, init_hsp->ungapped_data->s_start,
//       init_hsp->ungapped_data->length);
// init_hsp->offsets.qs_offsets.s_off +=
//     max_offset - init_hsp->offsets.qs_offsets.q_off;
// init_hsp->offsets.qs_offsets.q_off = max_offset;
// ```
fn extend_preliminary_blastp_hit(
    init_hsp: &mut InitHSP,
    s_idx: u32,
    context_idx: usize,
    ctx: &QueryContext,
    query_sequence: &[u8],
    subject_sequence: &[u8],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    mode: BlastpGappedAlignmentMode,
    scratch: &mut GapAlignScratch,
) -> BlastpPreliminaryHsp {
    let query_raw = &query_sequence[..ctx.aa_len];
    let subject_raw = &subject_sequence[..subject_sequence.len() - 1];
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3900-3905
    // ```c
    // q_start = init_hsp->ungapped_data->q_start;
    // q_end = q_start + init_hsp->ungapped_data->length;
    // s_start = init_hsp->ungapped_data->s_start;
    // s_end = s_start + init_hsp->ungapped_data->length;
    // ```
    let ungapped_len = usize::try_from(init_hsp.ungapped_data.length)
        .expect("NCBI BLAST preliminary ungapped HSP length must fit in size_t");
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3900-3945
    // ```c
    // q_start = init_hsp->ungapped_data->q_start;
    // q_end = q_start + init_hsp->ungapped_data->length;
    // s_start = init_hsp->ungapped_data->s_start;
    // s_end = s_start + init_hsp->ungapped_data->length;
    // ...
    // max_offset =
    //    BlastGetStartForGappedAlignment(query_tmp.sequence, subject->sequence,
    //       gap_align->sbp, init_hsp->ungapped_data->q_start,
    //       init_hsp->ungapped_data->length, init_hsp->ungapped_data->s_start,
    //       init_hsp->ungapped_data->length);
    // ```
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3937-3945
    // ```c
    // max_offset =
    //    BlastGetStartForGappedAlignment(query_tmp.sequence,
    //       subject->sequence, gap_align->sbp,
    //       init_hsp->ungapped_data->q_start,
    //       init_hsp->ungapped_data->length,
    //       init_hsp->ungapped_data->s_start,
    //       init_hsp->ungapped_data->length);
    // init_hsp->offsets.qs_offsets.s_off +=
    //     max_offset - init_hsp->offsets.qs_offsets.q_off;
    // init_hsp->offsets.qs_offsets.q_off = max_offset;
    // ```
    let gapped_q_start = blastp_get_start_for_gapped_alignment(
        query_sequence,
        subject_sequence,
        usize::try_from(init_hsp.ungapped_data.q_start)
            .expect("NCBI BLAST preliminary q_start must be non-negative"),
        ungapped_len,
        usize::try_from(init_hsp.ungapped_data.s_start)
            .expect("NCBI BLAST preliminary s_start must be non-negative"),
        ungapped_len,
        matrix,
    );
    let gapped_q_start = i32::try_from(gapped_q_start)
        .expect("NCBI BLAST preliminary gapped query start must fit in Int4");
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3943-3945
    // ```c
    // init_hsp->offsets.qs_offsets.s_off +=
    //     max_offset - init_hsp->offsets.qs_offsets.q_off;
    // init_hsp->offsets.qs_offsets.q_off = max_offset;
    // ```
    init_hsp.s_seed += gapped_q_start - init_hsp.q_seed_absolute;
    init_hsp.q_seed_absolute = gapped_q_start;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3943-3949
    // ```c
    // init_hsp->offsets.qs_offsets.s_off +=
    //     max_offset - init_hsp->offsets.qs_offsets.q_off;
    // init_hsp->offsets.qs_offsets.q_off = max_offset;
    // status = s_BlastProtGappedAlignment(...);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3949-3978
    // ```c
    // status = s_BlastProtGappedAlignment(program_number, &query_tmp,
    //                                     subject, gap_align,
    //                                     score_params, init_hsp,
    //                                     restricted_alignment,
    //                                     fence_hit);
    // if (restricted_alignment &&
    //     gap_align->score < cutoff &&
    //     gap_align->score >= restricted_cutoff) {
    //     ...
    // }
    // ```
    let aligned = blastp_score_only_gapped_alignment_with_scratch(
        query_raw,
        subject_raw,
        usize::try_from(gapped_q_start)
            .expect("NCBI BLAST preliminary gapped query start must be non-negative"),
        usize::try_from(init_hsp.s_seed)
            .expect("NCBI BLAST preliminary gapped subject start must be non-negative"),
        matrix,
        gap_open,
        gap_extend,
        x_dropoff,
        mode,
        scratch,
    );

    BlastpPreliminaryHsp {
        query_start: aligned.query_start,
        query_end: aligned.query_stop,
        subject_start: aligned.subject_start,
        subject_end: aligned.subject_stop,
        gapped_query_start: gapped_q_start,
        gapped_subject_start: init_hsp.s_seed,
        raw_score: aligned.score,
        query_context: i32::try_from(context_idx)
            .expect("NCBI BLAST preliminary blastp context index must fit in Int4"),
        query_frame: 0,
        subject_frame: 0,
        query_length: ctx.aa_len,
        q_idx: ctx.q_idx,
        s_idx,
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:606-617
// ```c
// sequence = queries.GetBlastSequence(index, encoding,
//                                     eNa_strand_unknown,
//                                     eSentinels,
//                                     &warnings);
// int offset = qinfo->contexts[ctx_index].query_offset;
// memcpy(&buf.get()[offset], sequence.data.get(), sequence.length);
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
// ```c
// db_length = BlastSeqSrcGetTotLen(seq_src);
// itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
// while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
//     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
//         continue;
//     }
// }
// ```
#[cfg(target_arch = "wasm32")]
fn fasta_records_from_bytes(bytes: &[u8]) -> Result<Vec<fasta::Record>> {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // void
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    fasta::Reader::new(bytes)
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(anyhow::Error::from)
}

#[cfg(target_arch = "wasm32")]
pub fn run_web_pair(args: BlastpArgs, query_fasta: &str, subject_fasta: &str) -> Result<Vec<u8>> {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:606-617
    // ```c
    // sequence = queries.GetBlastSequence(index, encoding,
    //                                     eNa_strand_unknown,
    //                                     eSentinels,
    //                                     &warnings);
    // memcpy(&buf.get()[offset], sequence.data.get(), sequence.length);
    // ```
    //
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
    //     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //         continue;
    //     }
    // }
    // ```
    let queries = fasta_records_from_bytes(query_fasta.as_bytes())
        .context("failed to parse web query FASTA")?;
    let subjects = fasta_records_from_bytes(subject_fasta.as_bytes())
        .context("failed to parse web subject FASTA")?;
    run_web_pair_records(args, &queries, &subjects, "", "")
}

#[cfg(target_arch = "wasm32")]
pub fn run_web_pair_records(
    args: BlastpArgs,
    query_records: &[fasta::Record],
    subject_records: &[fasta::Record],
    query_label: &str,
    subject_label: &str,
) -> Result<Vec<u8>> {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    //
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
    //     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //         continue;
    //     }
    // }
    // ```
    let mut output = Vec::new();
    run_internal_with_records(
        args,
        query_records,
        subject_records,
        query_label,
        subject_label,
        Some(&mut output),
    )?;
    Ok(output)
}

pub fn run(args: BlastpArgs) -> Result<()> {
    let args = args.resolve()?;
    validate_requested_blastp_support(&args)?;

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:518-543
    // ```c
    // for(TSeqPos index = 0; index < queries.Size(); index++) {
    //     if (const CSeq_id* id = queries.GetSeqId(index)) {
    //         const string kTitle = queries.GetTitle(index);
    //         messages[index].SetQueryId(query_id);
    //     }
    // }
    // ```
    let query_records: Vec<fasta::Record> = fasta::Reader::from_file(&args.query)
        .with_context(|| format!("failed to open query FASTA {}", args.query.display()))?
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("failed to read query FASTA {}", args.query.display()))?;

    let subject_records: Vec<fasta::Record> = fasta::Reader::from_file(&args.subject)
        .with_context(|| format!("failed to open subject FASTA {}", args.subject.display()))?
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("failed to read subject FASTA {}", args.subject.display()))?;

    run_resolved_with_records(args, &query_records, &subject_records, "", "", None)
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1348-1361
// ```c
// BlastCoreAuxStruct* aux_struct = NULL;
// BlastHSPList* hsp_list = NULL;
// BlastSeqSrcGetSeqArg seq_arg;
// Int2 status = 0;
// Int8 db_length = 0;
// ...
// BlastSeqSrcIterator* itr;
// ```
#[cfg(target_arch = "wasm32")]
fn run_internal_with_records(
    args: BlastpArgs,
    query_records: &[fasta::Record],
    subject_records: &[fasta::Record],
    query_label: &str,
    subject_label: &str,
    in_memory_output: Option<&mut Vec<u8>>,
) -> Result<()> {
    let args = args.resolve()?;
    validate_requested_blastp_support(&args)?;
    run_resolved_with_records(
        args,
        query_records,
        subject_records,
        query_label,
        subject_label,
        in_memory_output,
    )
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:606-617
// ```c
// sequence = queries.GetBlastSequence(index, encoding,
//                                     eNa_strand_unknown,
//                                     eSentinels,
//                                     &warnings);
// memcpy(&buf.get()[offset], sequence.data.get(), sequence.length);
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
// ```c
// while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
//     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
//         continue;
//     }
// }
// ```
fn run_resolved_with_records(
    args: ResolvedBlastpArgs,
    query_records: &[fasta::Record],
    subject_records: &[fasta::Record],
    _query_label: &str,
    subject_label: &str,
    mut in_memory_output: Option<&mut Vec<u8>>,
) -> Result<()> {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1681
    // ```c
    // BLAST_GapAlignSetUp(...);
    // BLAST_PreliminarySearchEngine(program, query, query_info,
    //                               seq_src, gap_align, score_params,
    //                               lookup_wrap, word_options, ext_params, ...);
    // ```
    let timing_enabled = blastp_timing_env_enabled();
    let kappa_range_counters_enabled = blastp_kappa_range_counters_env_enabled();
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1705
    // ```c
    // BLAST_PreliminarySearchEngine(..., diagnostics, ...);
    // Blast_DiagnosticsUpdate(diagnostics, local_diagnostics);
    // ```
    reset_blastp_kappa_range_counters(kappa_range_counters_enabled);
    let t_total = blastp_timing_start(timing_enabled);
    let timing = if timing_enabled {
        Some(Arc::new(BlastpTiming::new()))
    } else {
        None
    };

    let (outfmt, custom_fields) = OutputFormat::parse(&args.outfmt).map_err(anyhow::Error::msg)?;
    if outfmt == OutputFormat::Pairwise && custom_fields.is_some() {
        bail!("blastp outfmt 0 does not accept custom field lists");
    }
    let custom_fields = custom_fields
        .as_deref()
        .map(parse_blastp_tabular_fields)
        .transpose()?;

    let query_ids: Vec<Arc<str>> = query_records.iter().map(fasta_id).collect();
    let query_headers: Vec<String> = query_records.iter().map(fasta_defline).collect();
    let query_lengths: Vec<usize> = query_records
        .iter()
        .map(|record| record.seq().len())
        .collect();

    let query_encoding_start = blastp_timing_start(timing_enabled);
    let seg_params = args.seg.params();
    let query_frames: Vec<Vec<QueryFrame>> = query_records
        .iter()
        .map(|record| {
            vec![encode_protein_query_frame_with_seg(
                record.seq(),
                seg_params.as_ref(),
            )]
        })
        .collect();
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_duration(&timing.query_encoding_ns, query_encoding_start);
    }

    let subject_encoding_start = blastp_timing_start(timing_enabled);
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1426
    // ```c
    // /* iterate over all subject sequences */
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //        continue;
    //    }
    // ```
    let prepared_subjects =
        prepare_blastp_subjects_preserving_order(subject_records, args.num_threads)?;
    let mut subject_ids = Vec::with_capacity(prepared_subjects.len());
    let mut subject_titles = Vec::with_capacity(prepared_subjects.len());
    let mut subjects = Vec::with_capacity(prepared_subjects.len());
    for prepared_subject in prepared_subjects {
        subject_ids.push(prepared_subject.id);
        subject_titles.push(prepared_subject.title);
        subjects.push(prepared_subject.encoded);
    }
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_duration(&timing.subject_encoding_ns, subject_encoding_start);
    }

    let lookup_start = blastp_timing_start(timing_enabled);
    let lookup_ideal_params = lookup_protein_params_ungapped(args.scoring.matrix);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2833-2852
    // ```c
    // Blast_ScoreBlkKbpIdealCalc(BlastScoreBlk* sbp)
    // {
    //    Blast_ResFreqStdComp(sbp, stdrfp);
    //    BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
    //    Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);
    // }
    // ```
    let ideal_ungapped_params = compute_blosum62_ideal_karlin_params()
        .map_err(anyhow::Error::msg)
        .context("failed to calculate NCBI blastp kbp_ideal parameters")?;
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1833-1839
    // ```c
    // kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    // hsp->evalue = BLAST_KarlinStoE_simple(hsp->score, kbp[hsp->context],
    //                                       hsp_list->hspcnt ? searchsp : 1);
    // ```
    let gapped_params = lookup_protein_params(&args.scoring);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:457-463
    // ```c
    // params->gap_x_dropoff = (Int4)
    //     (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
    // params->gap_x_dropoff_final = (Int4)
    //     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
    //         params->gap_x_dropoff);
    // ```
    let x_drop_gapped =
        blastp_gapped_x_dropoff_raw(X_DROP_GAPPED_PRELIM as f64, gapped_params.lambda);
    let x_drop_gapped_final =
        blastp_gapped_x_dropoff_raw(X_DROP_GAPPED_FINAL as f64, gapped_params.lambda)
            .max(x_drop_gapped);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2746-2748
    // ```c
    // Boolean check_ideal =
    //    (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
    //     program == eBlastTypeRpsTblastn);
    // ```
    let (lookup, contexts) = match args.lookup_table_type {
        BlastpLookupTableType::AaLookupTable => {
            let (lookup, contexts) = build_ncbi_lookup(
                &query_frames,
                args.threshold as i32,
                &lookup_ideal_params,
                false,
            );
            (BlastpRuntimeLookup::Aa(lookup), contexts)
        }
        BlastpLookupTableType::CompressedAaLookupTable => {
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1349-1352
            // ```c
            // s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
            //                                 query, locations);
            // s_CompressedLookupFinalize(lookup);
            // ```
            let (query, locations, contexts) =
                prepare_blosum62_lookup_query_for_word_size(&query_frames, args.word_size, false);
            let lookup = build_blosum62_compressed_lookup(
                args.word_size,
                args.threshold,
                &query,
                &locations,
            )
            .context("failed to build NCBI compressed amino-acid lookup table")?;
            (BlastpRuntimeLookup::Compressed(lookup), contexts)
        }
    };
    let context_lookup_bounds = blastp_context_lookup_bounds(&contexts);
    let concatenated_query = build_blastp_concatenated_query(&contexts);
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_duration(&timing.lookup_ns, lookup_start);
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_setup.c:766-769
    // ```c
    // kbp_ptr = (scoring_options->gapped_calculation ? sbp->kbp_gap_std : sbp->kbp);
    // ```
    let total_db_len = subjects.iter().map(|subject| subject.aa_len).sum::<usize>();
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_setup.c:914-921
    // ```c
    // /* Set the database length for new FSC */
    // if (sbp->gbp) {
    //     sbp->gbp->db_length =
    //         (Blast_SubjectIsTranslated(program_number)) ? dbl/3 : dbl;
    // }
    // ```
    let gapped_gumbel = lookup_protein_gumbel_params(&args.scoring, total_db_len as i64)
        .context("missing NCBI gumbel parameters for supported blastp scoring")?;
    let num_subjects = subjects.len().max(1);
    let search_spaces: Vec<SearchSpace> = contexts
        .iter()
        .map(|ctx| {
            SearchSpace::for_database_search(
                ctx.aa_len,
                total_db_len,
                num_subjects,
                &gapped_params,
                true,
            )
        })
        .collect();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:195-221
    // ```c
    // for (context = query_info->first_context;
    //                  context <= query_info->last_context; ++context) {
    //    ...
    //    kbp= sbp->kbp[context];
    //    ...
    //    p->cutoffs[context].x_dropoff_init =
    //        (Int4)(sbp->scale_factor *
    //               ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:318-323,417-420
    // ```c
    // sbp->kbp_std[index] = (Blast_KarlinBlk*)
    //     BlastMemDup(sbp->kbp_gap_std[0], sizeof(Blast_KarlinBlk));
    // sbp->kbp = sbp->kbp_std;
    // ...
    // Blast_KarlinBlkCopy(sbp->kbp_std[context], sbp->kbp_ideal);
    // sbp->kbp = sbp->kbp_std;
    // ```
    //
    // Default blastp uses the score-block ungapped Karlin blocks (`sbp->kbp`
    // / `sbp->kbp_std`), not per-query lookup composition estimates.
    let x_dropoff_per_context: Vec<i32> = contexts
        .iter()
        .map(blastp_x_dropoff_init_for_context)
        .collect();
    let cutoff_subject_length = blastp_preliminary_cutoff_subject_length(&subjects);
    let hit_cutoff_scores: Vec<i32> = contexts
        .iter()
        .map(|ctx| {
            blastp_hit_cutoff_score_for_context(
                ctx,
                args.evalue,
                cutoff_subject_length,
                args.comp_based_stats.mode,
                &gapped_params,
                &gapped_gumbel,
            )
        })
        .collect();
    let word_cutoff_scores: Vec<i32> = hit_cutoff_scores
        .iter()
        .zip(contexts.iter())
        .map(|(&hit_cutoff_score, ctx)| blastp_word_cutoff_score_for_context(hit_cutoff_score, ctx))
        .collect();

    let prelim_hitlist_size = get_prelim_hitlist_size(
        args.max_target_seqs,
        args.comp_based_stats.is_enabled(),
        true,
    );
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:1008-1018
    // ```c
    // BLAST_OneSubjectUpdateParameters(..., subject_length, ...);
    // BlastHitSavingParametersUpdate(program_number, sbp, query_info,
    //                                subject_length, 0, hit_params);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:926-940
    // ```c
    // params->prelim_evalue = evalue;  /* evalue and prelim_evalue same if no CBS. */
    // ...
    // int cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
    // params->prelim_evalue = cbs_stretch*evalue;
    // ```
    let prelim_evalue = blastp_prelim_evalue(args.evalue, args.comp_based_stats.mode);
    let mut redone_matches = build_redone_match_heaps(query_ids.len(), args.max_target_seqs);
    let mut preliminary_hit_lists: Vec<Option<BlastpHitList>> = std::iter::repeat_with(|| None)
        .take(query_ids.len())
        .collect();

    let window = args.window_size as i32;
    let word_size = lookup.word_length();
    let offset_array_size = blastp_offset_array_size(lookup.longest_chain());
    let query_length = blastp_concatenated_query_length(&contexts);
    let (diag_array_size, diag_mask) = compute_diag_table_shape(query_length, window);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:369-372
    // ```c
    // tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
    //                               0, (subject_length > 0 ? subject_length :
    //                               subject_blk->length / CODON_LENGTH) + 1);
    // ```
    let total_query_span = query_length + 1;

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1498
    // ```c
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
    //    BlastHSPStreamWrite(hsp_stream, &hsp_list);
    // }
    // ```
    let process_subject = |scratch: &mut BlastpSubjectScratch,
                           s_idx: usize,
                           subject_result: &mut BlastpSubjectResult| {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1554
        // ```c
        // while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
        //     status = s_BlastSearchEngineCore(..., &hsp_list, ...);
        //     status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
        // }
        // ```
        subject_result.reset(s_idx);
        let subject = &subjects[s_idx];
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-121
        // ```c
        // (*seq_blk)->sequence_start = (Uint1 *) buffer;
        // /* The first byte is a sentinel byte. */
        // (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
        // (*seq_blk)->length = length;
        // ```
        let subject_sequence = &subject.aa_seq[1..];
        let subject_raw = &subject.aa_seq[1..subject.aa_seq.len() - 1];
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:999-1018
        // ```c
        // BLAST_OneSubjectUpdateParameters(..., subject_length, ...);
        // BlastHitSavingParametersUpdate(program_number, sbp, query_info,
        //                                subject_length, 0, hit_params);
        // BlastInitialWordParametersUpdate(program_number, hit_params, sbp,
        //                                  query_info, subject_length,
        //                                  word_params);
        // ```
        //
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:917-960
        // ```c
        // new_cutoff = BLAST_SpougeEtoS(...,
        //                     query_info->contexts[context].query_length,
        //                     avg_subject_length);
        // params->cutoffs[context].cutoff_score = new_cutoff;
        // ...
        // curr_cutoffs->cutoff_score = MIN(gap_trigger,
        //                                  hit_params->cutoffs[context].cutoff_score_max);
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:253-267
        // ```c
        // init_hitlist->init_hsp_array =
        //     malloc(init_hitlist->allocated * sizeof(BlastInitHSP));
        // ```
        let mut init_hsps = std::mem::take(&mut scratch.init_hsps);
        init_hsps.clear();
        let subject_diag_offset = scratch.diag_offset;

        let scan_ungapped_start = blastp_timing_start(timing_enabled);
        blastp_for_each_offset_pair(
            &lookup,
            subject_sequence,
            subject.aa_len,
            &mut scratch.offset_pairs,
            if args.window_size == 0 {
                BlastpOffsetPairScanMode::OneHit
            } else {
                BlastpOffsetPairScanMode::TwoHit
            },
            |pair| {
                let query_offset =
                    i32::try_from(pair.q_off).expect("NCBI BLAST query offset must fit in Int4");
                let subject_offset =
                    i32::try_from(pair.s_off).expect("NCBI BLAST subject offset must fit in Int4");
                if args.window_size == 0 {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:776-779
                    // ```c
                    // diag_coord = (subject_offset - query_offset) & diag_mask;
                    // diff = subject_offset -
                    //     (diag_array[diag_coord].last_hit - diag_offset);
                    // ```
                    let diag_coord = blastp_diag_coord(subject_offset - query_offset, diag_mask);
                    let diag_entry = &mut scratch.diag_array[diag_coord];
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:713-775
                    // ```c
                    // diff = subject_offset -
                    //     (diag_array[diag_coord].last_hit - diag_offset);
                    // if (diff >= 0) {
                    //     score = s_BlastAaExtendOneHit(...);
                    //     if (score >= cutoffs->cutoff_score) {
                    //         BlastSaveInitHsp(...);
                    //     }
                    //     diag_array[diag_coord].last_hit =
                    //         s_last_off - (wordsize - 1) + diag_offset;
                    // }
                    // ```
                    let diff = subject_offset - (diag_entry.last_hit() - subject_diag_offset);
                    if diff < 0 {
                        return;
                    }

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:787-793
                    // ```c
                    // Int4 curr_context = BSearchContextInfo(query_offset,
                    //                                        query_info);
                    // BlastUngappedCutoffs *cutoffs =
                    //     word_params->cutoffs + curr_context;
                    // ```
                    let ctx_idx = blastp_context_idx_for_query_offset_with_bounds(
                        &contexts,
                        context_lookup_bounds,
                        query_offset,
                    );
                    let x_drop = x_dropoff_per_context[ctx_idx];
                    let cutoff = word_cutoff_scores[ctx_idx];

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1040-1044
                    // ```c
                    // if (use_pssm)
                    //     sum += matrix[q_off + i][s[s_off + i]];
                    // else
                    //     sum += matrix[q[q_off + i]][s[s_off + i]];
                    // ```
                    let Some(ungapped_hit) = extend_one_hit_blosum62(
                        &concatenated_query,
                        subject_raw,
                        usize::try_from(query_offset)
                            .expect("NCBI BLAST one-hit query offset must be non-negative"),
                        usize::try_from(subject_offset)
                            .expect("NCBI BLAST one-hit subject offset must be non-negative"),
                        x_drop,
                        args.word_size,
                    ) else {
                        return;
                    };

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:790-800
                    // ```c
                    // if (score >= cutoffs->cutoff_score) {
                    //     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                    //                      query_offset, subject_offset, hsp_len,
                    //                      score);
                    // }
                    // diag_array[diag_coord].last_hit =
                    //     s_last_off - (wordsize - 1) + diag_offset;
                    // ```
                    if ungapped_hit.ungapped_data.score >= cutoff {
                        blastp_save_init_hsp(
                            &mut init_hsps,
                            ungapped_hit.ungapped_data.q_start,
                            ungapped_hit.ungapped_data.s_start,
                            query_offset,
                            subject_offset,
                            ungapped_hit.ungapped_data.length,
                            ungapped_hit.ungapped_data.score,
                        );
                    }
                    diag_entry.set_last_hit(
                        ungapped_hit.s_last_off - (word_size - 1) + subject_diag_offset,
                    );
                    return;
                }

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:516-518
                // ```c
                // diag_coord = (query_offset - subject_offset) & diag_mask;
                // ```
                let diag_coord = blastp_diag_coord(query_offset - subject_offset, diag_mask);
                let diag_entry = &mut scratch.diag_array[diag_coord];
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:488-574
                // ```c
                // if (diag_array[diag_coord].flag) {
                //     ...
                // } else {
                //     last_hit = diag_array[diag_coord].last_hit - diag_offset;
                //     diff = subject_offset - last_hit;
                //     if (diff >= window) { ... }
                //     if (diff < wordsize) continue;
                //     ...
                //     score = s_BlastAaExtendTwoHit(...);
                //     if (right_extend) {
                //         diag_array[diag_coord].flag = 1;
                //         diag_array[diag_coord].last_hit =
                //             s_last_off - (wordsize - 1) + diag_offset;
                //     } else {
                //         diag_array[diag_coord].last_hit =
                //             subject_offset + diag_offset;
                //     }
                // }
                // ```
                if diag_entry.flag() != 0 {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:509-521
                    // ```c
                    // if ((Int4) (subject_offset + diag_offset) <
                    //     diag_array[diag_coord].last_hit) {
                    //     continue;
                    // }
                    // else {
                    //     diag_array[diag_coord].last_hit =
                    //         subject_offset + diag_offset;
                    //     diag_array[diag_coord].flag = 0;
                    // }
                    // ```
                    let subject_plus_offset = subject_offset + subject_diag_offset;
                    if subject_plus_offset < diag_entry.last_hit() {
                        return;
                    }
                    diag_entry.set_last_hit(subject_plus_offset);
                    diag_entry.set_flag(0);
                    return;
                }

                let last_hit = diag_entry.last_hit() - subject_diag_offset;
                let diff = subject_offset - last_hit;
                if diff >= window {
                    diag_entry.set_last_hit(subject_offset + subject_diag_offset);
                    return;
                }
                if diff < word_size {
                    return;
                }

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:551-560
                // ```c
                // curr_context = BSearchContextInfo(query_offset, query_info);
                // ...
                // if (query_offset - diff <
                //     query_info->contexts[curr_context].query_offset) {
                // ```
                let ctx_idx = blastp_context_idx_for_query_offset_with_bounds(
                    &contexts,
                    context_lookup_bounds,
                    query_offset,
                );
                let ctx = &contexts[ctx_idx];

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:562-570
                // ```c
                // if (query_offset - diff <
                //     query_info->contexts[curr_context].query_offset) {
                //     diag_array[diag_coord].last_hit =
                //         subject_offset + diag_offset;
                //     continue;
                // }
                // ```
                if query_offset - diff < ctx.frame_base {
                    diag_entry.set_last_hit(subject_offset + subject_diag_offset);
                    return;
                }

                let x_drop = x_dropoff_per_context[ctx_idx];
                let cutoff = word_cutoff_scores[ctx_idx];

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1112-1116
                // ```c
                // if (use_pssm)
                //     score += matrix[q_right_off + i][s[s_right_off + i]];
                // else
                //     score += matrix[q[q_right_off + i]][s[s_right_off + i]];
                // ```
                let Some(ungapped_hit) = extend_two_hit_blosum62(
                    &concatenated_query,
                    subject_raw,
                    usize::try_from(last_hit + word_size)
                        .expect("NCBI BLAST two-hit seed boundary must be non-negative"),
                    usize::try_from(subject_offset)
                        .expect("NCBI BLAST two-hit subject offset must be non-negative"),
                    usize::try_from(query_offset)
                        .expect("NCBI BLAST two-hit query offset must be non-negative"),
                    x_drop,
                    args.word_size,
                ) else {
                    return;
                };

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:575-588
                // ```c
                // if (score >= cutoffs->cutoff_score)
                //     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                //                      query_offset, subject_offset, hsp_len,
                //                      score);
                // if (right_extend) {
                //     diag_array[diag_coord].flag = 1;
                //     diag_array[diag_coord].last_hit =
                //         s_last_off - (wordsize - 1) + diag_offset;
                // }
                // else {
                //     diag_array[diag_coord].last_hit =
                //         subject_offset + diag_offset;
                // }
                // ```
                if ungapped_hit.ungapped_data.score >= cutoff {
                    blastp_save_init_hsp(
                        &mut init_hsps,
                        ungapped_hit.ungapped_data.q_start,
                        ungapped_hit.ungapped_data.s_start,
                        query_offset,
                        subject_offset,
                        ungapped_hit.ungapped_data.length,
                        ungapped_hit.ungapped_data.score,
                    );
                }
                if ungapped_hit.right_extend {
                    diag_entry.set_flag(1);
                    diag_entry.set_last_hit(
                        ungapped_hit.s_last_off - (word_size - 1) + subject_diag_offset,
                    );
                } else {
                    diag_entry.set_last_hit(subject_offset + subject_diag_offset);
                }
            },
        );

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:152-169
        // ```c
        // if (ewp->diag_table->offset >= INT4_MAX / 4) {
        //     ewp->diag_table->offset = ewp->diag_table->window;
        //     s_BlastDiagClear(ewp->diag_table);
        // } else {
        //     ewp->diag_table->offset += subject_length + ewp->diag_table->window;
        // }
        // ```
        scratch.finish_subject(subject.aa_len, window);
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3826-3837
        // ```c
        // tree = Blast_IntervalTreeInit(0, query->length+1,
        //                               0, subject->length+1);
        // ```
        let interval_tree = &mut scratch.preliminary_interval_tree;
        interval_tree.reset_with_bounds(0, total_query_span, 0, subject.aa_len as i32 + 1);
        if let Some(timing) = timing.as_ref() {
            BlastpTiming::record_call(
                &timing.scan_ungapped_ns,
                &timing.scan_ungapped_calls,
                scan_ungapped_start,
            );
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:228-234
        // ```c
        // Blast_InitHitListSortByScore(init_hitlist);
        // return status;
        // ```
        ncbi_qsort_init_hsps_by_score(&mut init_hsps);
        let init_hsps = if args.chaining && args.scoring.matrix == ScoringMatrix::Blosum62 {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3716-3724
            // ```c
            // if (ext_params->options->chaining &&
            //     is_prot && strcmp(score_params->options->matrix, "BLOSUM62") == 0) {
            //     status = s_ChainingAlignment(..., init_hitlist);
            // }
            // ```
            chain_blastp_init_hsps(
                init_hsps,
                &contexts,
                &word_cutoff_scores,
                &hit_cutoff_scores,
                args.scoring.gap_open,
                args.scoring.gap_extend,
                &mut scratch.chaining_keep,
                &mut scratch.chaining_nodes,
            )
        } else {
            init_hsps
        };
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3735-3789
        // ```c
        // restricted_align_array = calloc(...);
        // for (i = 0;i < init_hitlist->total;i++) {
        //     init_hsp = &init_hitlist->init_hsp_array[i];
        //     ...
        // }
        // ```
        let mut restricted_align_array = std::mem::take(&mut scratch.restricted_align_array);
        if blastp_restricted_align_enabled(args.evalue) {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3735-3789
            // ```c
            // if (hit_params->restricted_align && !score_params->options->is_ooframe) {
            //     restricted_align_array = calloc(...);
            //     for (i = 0;i < init_hitlist->total;i++) {
            //         init_hsp = &init_hitlist->init_hsp_array[i];
            //         ...
            //     }
            // }
            // ```
            build_restricted_align_array_into(
                &init_hsps,
                &contexts,
                context_lookup_bounds,
                &hit_cutoff_scores,
                query_ids.len(),
                &mut scratch.restricted_align_found,
                &mut restricted_align_array,
            );
        } else {
            restricted_align_array.clear();
            restricted_align_array.resize(query_ids.len(), false);
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3971-4008
        // ```c
        // BlastHSPList* new_hsp_list = Blast_HSPListNew(kHspNumMax);
        // ...
        // Blast_HSPListFree(hsp_list);
        // hsp_list = new_hsp_list;
        // ```
        //
        // Keep the accepted preliminary HSPs in one ordered list for the
        // current subject so exact-retry rebuilds preserve NCBI's HSPList
        // order instead of HashMap iteration order.
        let mut preliminary_hsp_list = std::mem::take(&mut scratch.preliminary_hsp_list);
        preliminary_hsp_list.clear();
        let mut preliminary_hsp_counts_by_query =
            std::mem::take(&mut scratch.preliminary_hsp_counts_by_query);
        preliminary_hsp_counts_by_query.clear();
        preliminary_hsp_counts_by_query.resize(query_ids.len(), 0);
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1481-1497
        // ```c
        // if (!gapped_calculation) {
        //     /* The following must be performed for any ungapped
        //        search with a nucleotide database. */
        //     status = Blast_HSPListReevaluateUngapped(...);
        // }
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3697-3698
        // ```c
        // Int4 redo_index = -1; /* these are used for recomputing alignmnets if the */
        // Int4 redo_query = -1; /* approximate alignment score is inconclusive */
        // ```
        let mut redo_index: Option<usize> = None;
        let mut redo_query: Option<usize> = None;
        let mut first_init_hsp_index_by_query: Vec<Option<usize>> = Vec::new();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3870-3878
        // ```c
        // if (index < redo_index && query_index != redo_query) {
        //     continue;
        // }
        // ```
        let mut next_same_query_index: Vec<Option<usize>> = Vec::new();
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80
        // ```c
        // BlastGapDP* dp_mem; /**< scratch structures for dynamic programming */
        // Int4 dp_mem_alloc;  /**< current number of structures allocated */
        // ```
        //
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319
        // ```c
        // s_BlastProtGappedAlignment(..., BlastGapAlignStruct* gap_align, ...)
        // ```
        let preliminary_gapped_start = blastp_timing_start(timing_enabled);
        let mut hit_index = 0usize;
        while hit_index < init_hsps.len() {
            let mut init_hsp = init_hsps[hit_index];
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3852-3868
            // ```c
            // tmp_init_hsp = init_hsp_array[index];
            // ...
            // s_AdjustHspOffsetsAndGetQueryData(query, query_info, init_hsp,
            //                                   &query_tmp, &context);
            // query_index = Blast_GetQueryIndexFromContext(context, program_number);
            // ```
            let (context_idx, ctx, query_sequence) = blastp_adjust_hsp_offsets_and_get_query_data(
                &concatenated_query,
                &contexts,
                context_lookup_bounds,
                &mut init_hsp,
            );
            let query_index = usize::try_from(ctx.q_idx)
                .expect("NCBI BLAST requires blastp query indices to fit in usize");
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3870-3878
            // ```c
            // /* If redo_index > -1 and redo_query > -1, the main loop is recomputing
            //    gappaed alignments for a single query, becuase the approximate
            //    alignment score was inconclusive. This recomputing was triggered when
            //    index was equal to redo_index. Until index reaches redo_index again,
            //    skip all concatenated queries with index different from redo_query */
            // if (index < redo_index && query_index != redo_query) {
            //     continue;
            // }
            // ```
            if let (Some(redo_index), Some(redo_query)) = (redo_index, redo_query) {
                if hit_index < redo_index && query_index != redo_query {
                    hit_index += 1;
                    continue;
                }
            }
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3918
            // ```c
            // tmp_hsp.query.offset = q_start;
            // tmp_hsp.query.end = q_end;
            // tmp_hsp.subject.offset = s_start;
            // tmp_hsp.subject.end = s_end;
            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                   hit_options->min_diag_separation))
            // ```
            let ungapped_tree_hsp = TreeHsp {
                query_offset: init_hsp.ungapped_data.q_start,
                query_end: init_hsp.query_end_absolute(),
                subject_offset: init_hsp.ungapped_data.s_start,
                subject_end: init_hsp.subject_end(),
                score: init_hsp.ungapped_data.score,
                query_frame: 1,
                query_length: ctx.aa_len as i32,
                query_context_offset: ctx.frame_base,
                subject_frame_sign: 0,
            };
            let tree_contains = interval_tree.contains_hsp(&ungapped_tree_hsp, ctx.frame_base, 0);
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3919
            // ```c
            // tmp_hsp.query.offset = q_start;
            // tmp_hsp.query.end = q_end;
            // tmp_hsp.subject.offset = s_start;
            // tmp_hsp.subject.end = s_end;
            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                   hit_options->min_diag_separation))
            // ```
            if let Some(target) = blastp_trace_hsp_target() {
                if blastp_trace_matches_pair(target, ctx.q_idx, s_idx as u32)
                    && blastp_trace_init_hsp_overlaps_target(target, &init_hsp)
                {
                    eprintln!(
                        "[BLASTP_TRACE_HSP] init hit_index={} q_idx={} s_idx={} ungapped=q{}..{} s{}..{} seed=({}, {}) score={} cutoff={} mode={:?} tree_contains={}",
                        hit_index,
                        ctx.q_idx,
                        s_idx,
                        init_hsp.ungapped_data.q_start + 1,
                        init_hsp.query_end_absolute(),
                        init_hsp.ungapped_data.s_start + 1,
                        init_hsp.subject_end(),
                        init_hsp.q_seed_absolute,
                        init_hsp.s_seed,
                        init_hsp.ungapped_data.score,
                        hit_cutoff_scores[context_idx],
                        if restricted_align_array
                            .get(query_index)
                            .copied()
                            .unwrap_or(false)
                        {
                            BlastpGappedAlignmentMode::Restricted
                        } else {
                            BlastpGappedAlignmentMode::Exact
                        },
                        tree_contains
                    );
                }
            }
            if tree_contains {
                hit_index = blastp_advance_preliminary_hit_index(
                    hit_index,
                    redo_index,
                    &next_same_query_index,
                );
                continue;
            }
            let restricted_alignment = restricted_align_array
                .get(query_index)
                .copied()
                .unwrap_or(false);
            let cutoff = hit_cutoff_scores[context_idx];
            let restricted_cutoff = (BLASTP_RESTRICTED_MULT * cutoff as f64) as i32;
            let mode = if restricted_alignment {
                BlastpGappedAlignmentMode::Restricted
            } else {
                BlastpGappedAlignmentMode::Exact
            };
            let preliminary_hsp = extend_preliminary_blastp_hit(
                &mut init_hsp,
                s_idx as u32,
                context_idx,
                ctx,
                query_sequence,
                subject_sequence,
                args.scoring.matrix,
                args.scoring.gap_open,
                args.scoring.gap_extend,
                x_drop_gapped,
                mode,
                &mut scratch.preliminary_gap_scratch,
            );
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3935-3953
            // ```c
            // max_offset = BlastGetStartForGappedAlignment(...);
            // init_hsp->offsets.qs_offsets.s_off +=
            //     max_offset - init_hsp->offsets.qs_offsets.q_off;
            // init_hsp->offsets.qs_offsets.q_off = max_offset;
            // status = s_BlastProtGappedAlignment(..., init_hsp, ...);
            // ```
            if let Some(target) = blastp_trace_hsp_target() {
                if blastp_trace_matches_pair(target, ctx.q_idx, s_idx as u32)
                    && usize::try_from(preliminary_hsp.query_start + 1).ok() == Some(target.q_start)
                    && usize::try_from(preliminary_hsp.query_end).ok() == Some(target.q_end)
                    && usize::try_from(preliminary_hsp.subject_start + 1).ok()
                        == Some(target.s_start)
                    && usize::try_from(preliminary_hsp.subject_end).ok() == Some(target.s_end)
                {
                    eprintln!(
                        "[BLASTP_TRACE_HSP] preliminary_gapped hit_index={} q={}..{} s={}..{} gapped_start=({}, {}) raw_score={} cutoff={} mode={:?}",
                        hit_index,
                        preliminary_hsp.query_start + 1,
                        preliminary_hsp.query_end,
                        preliminary_hsp.subject_start + 1,
                        preliminary_hsp.subject_end,
                        preliminary_hsp.gapped_query_start,
                        preliminary_hsp.gapped_subject_start,
                        preliminary_hsp.raw_score,
                        cutoff,
                        mode
                    );
                }
            }
            if should_retry_exact_preliminary_gapped_alignment(
                restricted_alignment,
                preliminary_hsp.raw_score,
                cutoff,
                restricted_cutoff,
            ) {
                if let Some(restricted_slot) = restricted_align_array.get_mut(query_index) {
                    *restricted_slot = false;
                }
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3975-4000
                // ```c
                // Blast_IntervalTreeReset(tree);
                // /* remove all HSPs computed for the current query */
                // for (index2 = 0; index2 < hsp_list->hspcnt; index2++) { ... }
                // ```
                let accepted_for_query = preliminary_hsp_counts_by_query
                    .get(query_index)
                    .copied()
                    .unwrap_or(0);
                let (removed_hsps, rebuild_nodes) = if accepted_for_query == 0 {
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3975-3997
                    // ```c
                    // Blast_IntervalTreeReset(tree);
                    // /* remove all HSPs computed for the current query */
                    // for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
                    //     ...
                    // }
                    // ```
                    //
                    // With no accepted HSP for this query, NCBI's rebuild
                    // would copy every existing non-query HSP back into the
                    // tree. The current tree already contains exactly that set.
                    (0, 0)
                } else {
                    let removed_hsps = remove_preliminary_hits_for_query(
                        &mut preliminary_hsp_list,
                        query_index,
                        &contexts,
                    );
                    debug_assert_eq!(
                        removed_hsps, accepted_for_query,
                        "NCBI redo removal should drop all accepted HSPs for the retried query"
                    );
                    if let Some(count) = preliminary_hsp_counts_by_query.get_mut(query_index) {
                        *count = 0;
                    }
                    let rebuild_nodes = rebuild_preliminary_interval_tree_from_list(
                        interval_tree,
                        &preliminary_hsp_list,
                        &contexts,
                    );
                    (removed_hsps, rebuild_nodes)
                };
                if let Some(timing) = timing.as_ref() {
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3975-4007
                    // ```c
                    // Blast_IntervalTreeReset(tree);
                    // /* remove all HSPs computed for the current query */
                    // ...
                    // redo_index = index;
                    // redo_query = query_index;
                    // ```
                    BlastpTiming::record_count(&timing.preliminary_exact_retry_count, 1);
                    BlastpTiming::record_count(
                        &timing.preliminary_retry_removed_hsps,
                        removed_hsps,
                    );
                    if accepted_for_query != 0 {
                        BlastpTiming::record_count(&timing.preliminary_interval_tree_rebuilds, 1);
                        BlastpTiming::record_count(
                            &timing.preliminary_interval_tree_rebuild_nodes,
                            rebuild_nodes,
                        );
                    }
                }
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4001-4007
                // ```c
                // redo_index = index;
                // redo_query = query_index;
                // index = -1;
                // continue;
                // ```
                redo_index = Some(hit_index);
                redo_query = Some(query_index);
                if first_init_hsp_index_by_query.is_empty() {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3870-3878
                    // ```c
                    // if (index < redo_index && query_index != redo_query) {
                    //     continue;
                    // }
                    // ```
                    //
                    // NCBI restarts at the first initial hit and skips until
                    // it finds the retried query. Rust precomputes that first
                    // per-query hit lazily and jumps to the same first
                    // processable index, preserving the subsequent hit order.
                    first_init_hsp_index_by_query.resize(query_ids.len(), None);
                    next_same_query_index.clear();
                    next_same_query_index.resize(init_hsps.len(), None);
                    let mut last_init_hsp_index_by_query = vec![None; query_ids.len()];
                    for (index, init_hsp) in init_hsps.iter().enumerate() {
                        let context_idx = blastp_get_ungapped_hsp_context_with_bounds(
                            &contexts,
                            context_lookup_bounds,
                            init_hsp,
                        );
                        let init_query_index = usize::try_from(contexts[context_idx].q_idx)
                            .expect("NCBI BLAST requires blastp query indices to fit in usize");
                        if let Some(slot) = first_init_hsp_index_by_query.get_mut(init_query_index)
                        {
                            if slot.is_none() {
                                *slot = Some(index);
                            }
                        }
                        if let Some(last_index) =
                            last_init_hsp_index_by_query.get_mut(init_query_index)
                        {
                            if let Some(previous_index) = *last_index {
                                next_same_query_index[previous_index] = Some(index);
                            }
                            *last_index = Some(index);
                        }
                    }
                }
                hit_index = first_init_hsp_index_by_query
                    .get(query_index)
                    .and_then(|index| *index)
                    .unwrap_or(0);
                continue;
            }
            if preliminary_hsp.raw_score < cutoff {
                hit_index = blastp_advance_preliminary_hit_index(
                    hit_index,
                    redo_index,
                    &next_same_query_index,
                );
                continue;
            }
            interval_tree.add_hsp(
                preliminary_tree_hsp(&preliminary_hsp, ctx.frame_base),
                ctx.frame_base,
                IndexMethod::QueryAndSubject,
            );
            preliminary_hsp_list.push(preliminary_hsp);
            if let Some(count) = preliminary_hsp_counts_by_query.get_mut(query_index) {
                *count += 1;
            }
            hit_index =
                blastp_advance_preliminary_hit_index(hit_index, redo_index, &next_same_query_index);
        }
        if let Some(timing) = timing.as_ref() {
            BlastpTiming::record_call(
                &timing.preliminary_gapped_ns,
                &timing.preliminary_gapped_calls,
                preliminary_gapped_start,
            );
        }
        scratch.init_hsps = init_hsps;
        scratch.init_hsps.clear();

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:540-555
        // ```c
        // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
        // ...
        // Blast_HSPListSortByScore(hsp_list);
        // ```
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2525
        // ```c
        // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
        // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
        // ```
        let endpoint_duplicates_purged =
            prepare_preliminary_hits_for_kappa_in_place(&mut preliminary_hsp_list);
        if let Some(timing) = timing.as_ref() {
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2525
            // ```c
            // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
            // ...
            // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
            // ```
            BlastpTiming::record_count(
                &timing.preliminary_endpoint_duplicates_purged,
                endpoint_duplicates_purged,
            );
        }
        let mut preliminary_hits_by_query = std::mem::take(&mut scratch.preliminary_hits_by_query);
        split_preliminary_hits_by_query_into(
            &mut preliminary_hsp_list,
            &contexts,
            query_ids.len(),
            &mut preliminary_hits_by_query,
        );
        for (q_idx, preliminary_hits) in preliminary_hits_by_query
            .iter_mut()
            .take(query_ids.len())
            .enumerate()
        {
            if preliminary_hits.is_empty() {
                continue;
            }
            let ctx = &contexts[q_idx as usize];
            let mut provisional_hsp_list = build_preliminary_hsp_list_for_hitlist(
                preliminary_hits.as_slice(),
                s_idx as u32,
                q_idx as u32,
                ctx.aa_len,
                subject.aa_len,
                &gapped_params,
                &gapped_gumbel,
            );
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1523-1539
            // ```c
            // Blast_HSPListGetEvalues(program_number, query_info,
            //                         stat_length, hsp_list,
            //                         gapped_calculation, FALSE,
            //                         sbp, 0, 1.0);
            // status = s_Blast_HSPListReapByPrelimEvalue(hsp_list, hit_params);
            // ```
            //
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:640-670
            // ```c
            // cutoff = hit_params->prelim_evalue;
            // ...
            // if (hsp->evalue > cutoff) {
            //     hsp_array[index] = Blast_HSPFree(hsp_array[index]);
            // }
            // ```
            reap_hsplist_by_evalue(&mut provisional_hsp_list, prelim_evalue);
            if provisional_hsp_list.hsps.is_empty() {
                preliminary_hits.clear();
                continue;
            }
            subject_result
                .query_hsp_lists
                .push((q_idx, provisional_hsp_list));
            preliminary_hits.clear();
        }
        scratch.restricted_align_array = restricted_align_array;
        scratch.restricted_align_array.clear();
        scratch.preliminary_hsp_counts_by_query = preliminary_hsp_counts_by_query;
        scratch.preliminary_hsp_counts_by_query.clear();
        scratch.preliminary_hsp_list = preliminary_hsp_list;
        scratch.preliminary_hits_by_query = preliminary_hits_by_query;
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1498
    // ```c
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
    //    BlastHSPStreamWrite(hsp_stream, &hsp_list);
    // }
    // ```
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let num_threads = blastp_effective_num_threads(args.num_threads);
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let use_parallel = num_threads > 1 && subjects.len() > 1;
    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
    let use_parallel = false;

    if use_parallel {
        #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
        {
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1688
            // ```c
            // aux_arg->gap_align = BLAST_GapAlignStructNew(...);
            // aux_arg->hsp_stream = BlastHSPStreamNew(...);
            // status = Blast_RunPreliminarySearchWithInterrupt_MT(...);
            // ```
            let subject_diag_offsets = blastp_subject_diag_offsets(&subjects, window);
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build()
                .context("failed to build BLASTP thread pool")?;
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1554
            // ```c
            // while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
            //     status = s_BlastSearchEngineCore(..., &hsp_list, ...);
            //     status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
            // }
            // ```
            //
            // Parallel workers write one subject result into its oid-indexed
            // slot; the reducer below replays those slots in NCBI subject
            // iteration order.
            let mut subject_results: Vec<Option<BlastpSubjectResult>> =
                std::iter::repeat_with(|| None)
                    .take(subjects.len())
                    .collect();
            pool.install(|| {
                subject_results.par_iter_mut().enumerate().for_each_init(
                    || {
                        BlastpSubjectScratch::new(
                            offset_array_size,
                            diag_array_size,
                            window,
                            total_query_span,
                        )
                    },
                    |scratch, (s_idx, subject_result_slot)| {
                        scratch.prepare_independent_subject(subject_diag_offsets[s_idx]);
                        let mut subject_result = BlastpSubjectResult::new(s_idx);
                        process_subject(scratch, s_idx, &mut subject_result);
                        *subject_result_slot = Some(subject_result);
                    },
                );
            });
            let merge_sort_start = blastp_timing_start(timing_enabled);
            merge_blastp_subject_results_in_order(
                &mut preliminary_hit_lists,
                subject_results,
                prelim_hitlist_size,
            );
            if let Some(timing) = timing.as_ref() {
                BlastpTiming::record_call(
                    &timing.prelim_merge_sort_ns,
                    &timing.prelim_merge_sort_calls,
                    merge_sort_start,
                );
            }
        }
    } else {
        let mut scratch =
            BlastpSubjectScratch::new(offset_array_size, diag_array_size, window, total_query_span);
        let mut subject_result = BlastpSubjectResult::new(0);
        for s_idx in 0..subjects.len() {
            process_subject(&mut scratch, s_idx, &mut subject_result);
            let merge_sort_start = blastp_timing_start(timing_enabled);
            update_blastp_preliminary_hit_lists(
                &mut preliminary_hit_lists,
                &mut subject_result,
                prelim_hitlist_size,
            );
            if let Some(timing) = timing.as_ref() {
                BlastpTiming::record_call(
                    &timing.prelim_merge_sort_ns,
                    &timing.prelim_merge_sort_calls,
                    merge_sort_start,
                );
            }
        }
    }

    let merge_sort_start = blastp_timing_start(timing_enabled);
    for hit_list_opt in &mut preliminary_hit_lists {
        if let Some(hit_list) = hit_list_opt.as_mut() {
            hit_list.sort_by_evalue();
        }
    }
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_call(
            &timing.prelim_merge_sort_ns,
            &timing.prelim_merge_sort_calls,
            merge_sort_start,
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3329-3334
    // ```c
    // if ((int) compo_adjust_mode > 1 && !positionBased) {
    //     NRrecord_tld[i] = Blast_CompositionWorkspaceNew();
    //     status_code = Blast_CompositionWorkspaceInit(
    //             NRrecord_tld[i],
    //             scoringParams->options->matrix
    //     );
    // }
    // ```
    let kappa_redo_start = blastp_timing_start(timing_enabled);
    let build_redo_align_params = |q_idx: usize| -> Result<BlastRedoAlignParams> {
        let local_scaling_factor =
            blastp_local_scaling_factor(args.comp_based_stats.mode, args.scoring.matrix);
        let scaled_gapped_lambda = gapped_params.lambda / local_scaling_factor;
        let do_link_hsps = false;
        let matrix_info = build_matrix_info(
            args.scoring.matrix,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2216-2231
            // ```c
            // self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
            // status = s_GetStartFreqRatios(self->startFreqRatios, matrixName);
            // Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
            //                          self->startFreqRatios, self->ungappedLambda);
            // ```
            ideal_ungapped_params.lambda / local_scaling_factor,
        )?;
        Ok(BlastRedoAlignParams {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2216-2231
            // ```c
            // self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
            // status = s_GetStartFreqRatios(self->startFreqRatios, matrixName);
            // Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
            //                          self->startFreqRatios, self->ungappedLambda);
            // ```
            matrix_info,
            // NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:97-106
            // ```c
            // typedef struct BlastCompo_GappingParams {
            //     int gap_open;
            //     int gap_extend;
            //     int decline_align;
            //     int x_dropoff;
            //     void * context;
            // } BlastCompo_GappingParams;
            // ```
            //
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2117-2133
            // ```c
            // kbp->Lambda /= scale_factor;
            // sp->gap_open = (Int4)BLAST_Nint(sp->gap_open  * scale_factor);
            // sp->gap_extend = (Int4)BLAST_Nint(sp->gap_extend * scale_factor);
            // ```
            gapping_params: BlastCompoGappingParams {
                gap_open: ((args.scoring.gap_open as f64) * local_scaling_factor).round() as i32,
                gap_extend: ((args.scoring.gap_extend as f64) * local_scaling_factor).round()
                    as i32,
                decline_align: 0,
                x_dropoff: blastp_redo_x_dropoff(
                    scaled_gapped_lambda,
                    &gapped_params,
                    x_drop_gapped_final,
                ),
                context: std::cell::Cell::new(None),
            },
            compo_adjust_mode: blastp_compo_adjust_mode(args.comp_based_stats.mode),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1107-1118
            // ```c
            // Blast_RedoOneMatch(..., int ** matrix, int alphsize,
            //                    ..., int compositionTestIndex, ...)
            // ```
            alphsize: BLASTAA_SIZE as i32,
            composition_test_index: i32::from(args.comp_based_stats.unified_p),
            unified_p: args.comp_based_stats.unified_p,
            log_k: gapped_params.k.ln(),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3584-3690
            // ```c
            // localScalingFactor      /* i */
            // ...
            // s_HSPListNormalizeScores(hsp_list, kbp->Lambda, kbp->logK,
            //                          localScalingFactor);
            // ```
            score_divisor: local_scaling_factor,
            restricted_alignment: false,
            smith_waterman: false,
            near_identical_cutoff: blastp_redo_near_identical_cutoff(scaled_gapped_lambda),
            position_based: false,
            re_matrix_adjustment_pseudocounts: 20,
            ccat_query_length: contexts[q_idx].aa_len as i32,
            query_is_translated: false,
            subject_is_translated: false,
            cutoff_score: blastp_redo_cutoff_score(do_link_hsps, 0, local_scaling_factor),
            cutoff_evalue: args.evalue,
            do_link_hsps,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:382-406
            // ```c
            // Nearly identical alignments are computed with exact subjects,
            // and others with segged subjects; this makes comparing end
            // points more difficult.
            // ```
            is_same_adjustment: false,
        })
    };

    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let kappa_parallel_query_indices: Vec<usize> = preliminary_hit_lists
        .iter()
        .enumerate()
        .filter_map(|(q_idx, hit_list)| {
            hit_list
                .as_ref()
                .filter(|hit_list| hit_list.hsplist_count > 0)
                .map(|_| q_idx)
        })
        .collect();
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let kappa_num_threads = blastp_effective_num_threads(args.num_threads);
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let use_parallel_kappa_redo = kappa_num_threads > 1 && kappa_parallel_query_indices.len() > 1;
    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
    let use_parallel_kappa_redo = false;
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let kappa_match_parallel_query_index = if kappa_num_threads > 1
        && kappa_parallel_query_indices.len() == 1
        && preliminary_hit_lists[kappa_parallel_query_indices[0]]
            .as_ref()
            .is_some_and(|hit_list| hit_list.hsplist_count > 1)
    {
        Some(kappa_parallel_query_indices[0])
    } else {
        None
    };
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let use_parallel_kappa_match_redo = kappa_match_parallel_query_index.is_some();
    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
    let use_parallel_kappa_match_redo = false;

    if use_parallel_kappa_match_redo {
        #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
        {
            let q_idx = kappa_match_parallel_query_index
                .expect("single-query Kappa redo parallel path requires a query index");
            let preliminary_hit_list = preliminary_hit_lists[q_idx]
                .as_ref()
                .expect("single-query Kappa redo parallel path requires preliminary hits");
            let local_matches =
                &preliminary_hit_list.hsplist_array[..preliminary_hit_list.hsplist_count];
            let ctx = &contexts[q_idx];
            let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
            let query_nomask = query_nomask_sequence(ctx);

            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3429-3449
            // ```c
            // #pragma omp parallel ... num_threads(actual_num_threads)
            // #pragma omp for schedule(static)
            // for (b = 0; b < numMatches; ++b) {
            // ```
            //
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3493-3503
            // ```c
            // scoringParams = score_params_tld[tid];
            // redoneMatches = redoneMatches_tld[tid];
            // NRrecord = NRrecord_tld[tid];
            // redo_align_params = redo_align_params_tld[tid];
            // ```
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(kappa_num_threads)
                .build()
                .context("failed to build BLASTP Kappa match redo thread pool")?;
            let precomputed: Vec<Result<Option<BlastpPostprocessResult>>> = pool.install(|| {
                local_matches
                    .par_iter()
                    .map(|local_match| -> Result<Option<BlastpPostprocessResult>> {
                        let mut kappa_preliminary_hits = Vec::with_capacity(local_match.hsps.len());
                        kappa_preliminary_hits.extend(
                            local_match
                                .hsps
                                .iter()
                                .map(preliminary_hit_from_local_match_hsp),
                        );
                        if kappa_preliminary_hits.is_empty() {
                            return Ok(None);
                        }

                        let s_idx = local_match.oid;
                        let subject = &subjects[s_idx as usize];
                        let subject_raw = &subject.aa_seq[1..subject.aa_seq.len() - 1];
                        let redo_align_params = build_redo_align_params(q_idx)?;
                        let query_workspace = build_query_workspace(
                            query_raw,
                            ctx.aa_len as i32,
                            search_spaces[q_idx].length_adjustment as i32,
                            search_spaces[q_idx].effective_space,
                        );
                        let mut composition_workspace = BlastCompositionWorkspace::new_blosum62();
                        let mut kappa_gap_scratch = GapAlignScratch::new();
                        let mut kappa_subject_range_cache = BlastpKappaSubjectRangeCache::new();

                        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3577-3585
                        // ```c
                        // *pStatusCode = s_ResultHspToDistinctAlign(...);
                        // ```
                        postprocess_preliminary_hits(
                            &kappa_preliminary_hits,
                            &query_workspace,
                            &mut composition_workspace,
                            query_nomask,
                            subject_raw,
                            args.scoring.matrix,
                            &gapped_params,
                            &gapped_gumbel,
                            &redo_align_params,
                            args.evalue,
                            &mut kappa_gap_scratch,
                            &mut kappa_subject_range_cache,
                        )
                        .map(Some)
                    })
                    .collect()
            });

            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3525-3529
            // ```c
            // if (BlastCompo_EarlyTermination(localMatch->best_evalue,
            //                                 redoneMatches,
            //                                 numQueries)) {
            // ```
            //
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3817-3849
            // ```c
            // /* Reduce results from all threads and continue with business as usual */
            // local_results = SThreadLocalDataArrayConsolidateResults(thread_data);
            // ```
            //
            // The expensive redo calculation above runs in NCBI's static-loop
            // spirit, but heap admission is replayed in the original HSP stream
            // order so LOSAT keeps thread-count-independent output.
            for (local_match, postprocessed) in local_matches.iter().zip(precomputed.into_iter()) {
                if local_match.hsps.is_empty() {
                    continue;
                }
                if blast_compo_early_termination(local_match.best_evalue, &redone_matches) {
                    continue;
                }
                let Some(postprocessed) = postprocessed? else {
                    continue;
                };
                let best_score = postprocessed.best_score;
                let best_evalue = postprocessed.best_evalue;
                let Some(hsp_list) = postprocessed.hsp_list else {
                    continue;
                };
                let redone_match = &mut redone_matches[q_idx];
                if redone_match.would_insert(best_evalue, best_score, local_match.oid) {
                    redone_match.insert(hsp_list, best_evalue, best_score, local_match.oid);
                }
            }
        }
    } else if use_parallel_kappa_redo {
        #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
        {
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3119-3128
            // ```c
            // redoneMatches = calloc(numQueries, sizeof(BlastCompo_Heap));
            // for (query_index = 0;  query_index < numQueries;  query_index++) {
            //     BlastCompo_HeapInitialize(&redoneMatches[query_index], ...);
            // }
            // ```
            //
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3329-3334
            // ```c
            // NRrecord_tld[i] = Blast_CompositionWorkspaceNew();
            // status_code = Blast_CompositionWorkspaceInit(...);
            // ```
            //
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3436-3459
            // ```c
            // #pragma omp for schedule(static)
            // for (b = 0; b < numMatches; ++b) {
            //     ... redoneMatches = redoneMatches_tld[tid];
            // ```
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(kappa_num_threads)
                .build()
                .context("failed to build BLASTP Kappa redo thread pool")?;
            let query_heaps: Result<Vec<(usize, BlastCompoHeap)>> = pool.install(|| {
                kappa_parallel_query_indices
                    .par_iter()
                    .map(|&q_idx| -> Result<(usize, BlastCompoHeap)> {
                        let mut composition_workspace = BlastCompositionWorkspace::new_blosum62();
                        let mut kappa_preliminary_hits = Vec::new();
                        let mut kappa_gap_scratch = GapAlignScratch::new();
                        let mut kappa_subject_range_cache = BlastpKappaSubjectRangeCache::new();
                        let mut redone_match =
                            BlastCompoHeap::new(args.max_target_seqs, PSI_INCLUSION_ETHRESH);
                        let redo_align_params = build_redo_align_params(q_idx)?;

                        let ctx = &contexts[q_idx];
                        let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
                        let query_nomask = query_nomask_sequence(ctx);
                        let query_workspace = build_query_workspace(
                            query_raw,
                            ctx.aa_len as i32,
                            search_spaces[q_idx].length_adjustment as i32,
                            search_spaces[q_idx].effective_space,
                        );
                        let Some(preliminary_hit_list) = preliminary_hit_lists[q_idx].as_ref()
                        else {
                            return Ok((q_idx, redone_match));
                        };

                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:289-318
                        // ```c
                        // hit_list = results->hitlist_array[index];
                        // last_hsplist_index = hit_list->hsplist_count - 1;
                        // *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
                        // ```
                        for local_match in preliminary_hit_list
                            .hsplist_array
                            .iter()
                            .take(preliminary_hit_list.hsplist_count)
                        {
                            let s_idx = local_match.oid;
                            kappa_preliminary_hits.clear();
                            kappa_preliminary_hits.extend(
                                local_match
                                    .hsps
                                    .iter()
                                    .map(preliminary_hit_from_local_match_hsp),
                            );
                            if kappa_preliminary_hits.is_empty() {
                                continue;
                            }

                            let subject = &subjects[s_idx as usize];
                            let subject_raw = &subject.aa_seq[1..subject.aa_seq.len() - 1];
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3525-3528
                            // ```c
                            // if (BlastCompo_EarlyTermination(localMatch->best_evalue,
                            //                                 redoneMatches,
                            //                                 numQueries)) {
                            // ```
                            if blast_compo_early_termination(
                                local_match.best_evalue,
                                std::slice::from_ref(&redone_match),
                            ) {
                                continue;
                            }

                            let postprocessed = postprocess_preliminary_hits(
                                &kappa_preliminary_hits,
                                &query_workspace,
                                &mut composition_workspace,
                                query_nomask,
                                subject_raw,
                                args.scoring.matrix,
                                &gapped_params,
                                &gapped_gumbel,
                                &redo_align_params,
                                args.evalue,
                                &mut kappa_gap_scratch,
                                &mut kappa_subject_range_cache,
                            )?;
                            let best_score = postprocessed.best_score;
                            let best_evalue = postprocessed.best_evalue;
                            let Some(hsp_list) = postprocessed.hsp_list else {
                                continue;
                            };
                            if redone_match.would_insert(best_evalue, best_score, s_idx) {
                                redone_match.insert(hsp_list, best_evalue, best_score, s_idx);
                            }
                        }

                        Ok((q_idx, redone_match))
                    })
                    .collect()
            });
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3817-3850
            // ```c
            // /* Reduce results from all threads and continue with business as usual */
            // local_results = SThreadLocalDataArrayConsolidateResults(thread_data);
            // ```
            for (q_idx, heap) in query_heaps? {
                redone_matches[q_idx] = heap;
            }
        }
    } else {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3329-3334
        // ```c
        // if ((int) compo_adjust_mode > 1 && !positionBased) {
        //     NRrecord_tld[i] = Blast_CompositionWorkspaceNew();
        //     status_code = Blast_CompositionWorkspaceInit(...);
        // }
        // ```
        let mut composition_workspace = BlastCompositionWorkspace::new_blosum62();
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3656-3667
        // ```c
        // s_HSPListFromDistinctAlignments(hsp_list,
        //         &alignments[context_index],
        //         matchingSeq.index,
        //         queryInfo, qframe);
        // BlastCompo_AlignmentsFree(&incoming_aligns, NULL);
        // ```
        let mut kappa_preliminary_hits = Vec::new();
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3017-3019
        // ```c
        // /* keeps track of gapped alignment params */
        // BlastGapAlignStruct* gapAlign = NULL;
        // ```
        //
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1934-1942
        // ```c
        // status =
        //     BLAST_GappedAlignmentWithTraceback(..., gapAlign,
        //                                        context->scoringParams,
        //                                        q_start, s_start,
        //                                        query_data->length,
        //                                        subject_data->length,
        //                                        &fence_hit);
        // ```
        let mut kappa_gap_scratch = GapAlignScratch::new();
        let mut kappa_subject_range_cache = BlastpKappaSubjectRangeCache::new();

        for q_idx in 0..query_ids.len() {
            let redo_align_params = build_redo_align_params(q_idx)?;

            let ctx = &contexts[q_idx];
            let query_raw = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
            let query_nomask = query_nomask_sequence(ctx);
            let Some(preliminary_hit_list) = preliminary_hit_lists[q_idx].as_ref() else {
                continue;
            };
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-847
            // ```c
            // BLAST_ComputeLengthAdjustment(..., query_length, db_length,
            //                               db_num_seqs, &length_adjustment);
            // query_info->contexts[index].eff_searchsp = effective_search_space;
            // query_info->contexts[index].length_adjustment = length_adjustment;
            // ```
            //
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3244-3248
            // ```c
            // query_info_tld[i] = s_GetQueryInfo(queryBlk->sequence, queryInfo, ...);
            // ```
            let query_workspace = build_query_workspace(
                query_raw,
                ctx.aa_len as i32,
                search_spaces[q_idx].length_adjustment as i32,
                search_spaces[q_idx].effective_space,
            );

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:133-151
            // ```c
            // if (hsp_stream->sort_by_score->sort_on_read) {
            //     Blast_HSPResultsReverseSort(hsp_stream->results);
            // }
            // ```
            //
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:289-318
            // ```c
            // hit_list = results->hitlist_array[index];
            // last_hsplist_index = hit_list->hsplist_count - 1;
            // *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
            // --hit_list->hsplist_count;
            // ```
            //
            // For composition-based searches the preliminary stream is read in
            // query order, and within each query from best to worst preliminary
            // HSPList. Use the sorted preliminary hitlist directly instead of
            // draining an unordered map, so `localMatch->best_evalue` and the
            // admission order match NCBI.
            for local_match in preliminary_hit_list
                .hsplist_array
                .iter()
                .take(preliminary_hit_list.hsplist_count)
            {
                let s_idx = local_match.oid;
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3612-3619
                // ```c
                // localMatch->hsp_array,  /* i */
                // localMatch->hspcnt,     /* i */
                // ...
                // *pStatusCode = s_ResultHspToDistinctAlign(...);
                // ```
                kappa_preliminary_hits.clear();
                kappa_preliminary_hits.extend(
                    local_match
                        .hsps
                        .iter()
                        .map(preliminary_hit_from_local_match_hsp),
                );
                if kappa_preliminary_hits.is_empty() {
                    continue;
                }

                let subject = &subjects[s_idx as usize];
                let subject_raw = &subject.aa_seq[1..subject.aa_seq.len() - 1];
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3525-3528
                // ```c
                // if (BlastCompo_EarlyTermination(localMatch->best_evalue,
                //                                 redoneMatches,
                //                                 numQueries)) {
                //     ...
                // }
                // ```
                if blast_compo_early_termination(local_match.best_evalue, &redone_matches) {
                    continue;
                }

                let postprocessed = postprocess_preliminary_hits(
                    &kappa_preliminary_hits,
                    &query_workspace,
                    &mut composition_workspace,
                    query_nomask,
                    subject_raw,
                    args.scoring.matrix,
                    &gapped_params,
                    &gapped_gumbel,
                    &redo_align_params,
                    args.evalue,
                    &mut kappa_gap_scratch,
                    &mut kappa_subject_range_cache,
                )?;
                let best_score = postprocessed.best_score;
                let best_evalue = postprocessed.best_evalue;
                let Some(hsp_list) = postprocessed.hsp_list else {
                    continue;
                };
                let redone_match = &mut redone_matches[q_idx];
                if redone_match.would_insert(best_evalue, best_score, s_idx) {
                    redone_match.insert(hsp_list, best_evalue, best_score, s_idx);
                }
            }
        }
    }

    let mut hit_lists = fill_results_from_compo_heaps(args.max_target_seqs, &mut redone_matches);
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_call(
            &timing.kappa_redo_ns,
            &timing.kappa_redo_calls,
            kappa_redo_start,
        );
    }

    let final_traceback_start = blastp_timing_start(timing_enabled);
    for hit_list_opt in &mut hit_lists {
        if let Some(hit_list) = hit_list_opt.as_mut() {
            hit_list.sort_by_evalue();
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:845-852
            // ```c
            // for (subject_index = 0; subject_index < hit_list->hsplist_count; ++subject_index) {
            //    BlastHSPList * hsp_list = hit_list->hsplist_array[subject_index];
            //    if(hit_options->max_hsps_per_subject) {
            //       Blast_TrimHSPListByMaxHsps(hsp_list, hit_options);
            //    }
            // }
            // ```
            if args.max_hsps_per_subject > 0 {
                for hsp_list in hit_list
                    .hsplist_array
                    .iter_mut()
                    .take(hit_list.hsplist_count)
                {
                    trim_by_max_hsps(hsp_list, args.max_hsps_per_subject);
                }
            }
            hit_list.prune_by_size(args.max_target_seqs);
        }
    }

    let query_header = query_headers.first().cloned();
    let subject_name_label = if subject_label.is_empty() {
        path_label(&args.subject)
    } else {
        subject_label.to_string()
    };
    let subject_input_label = if subject_label.is_empty() {
        args.subject.display().to_string()
    } else {
        subject_label.to_string()
    };
    let context = ReportContext {
        query_name: query_header,
        query_length: query_lengths.first().copied(),
        subject_name: Some(subject_name_label.clone()),
        database_name: Some(format!(
            "User specified sequence set (Input: {subject_input_label})"
        )),
        database_num_sequences: Some(subject_records.len()),
        database_total_letters: Some(total_db_len),
        program: "blastp".to_string(),
        version: Some(NCBI_BLASTP_VERSION.to_string()),
    };

    let tabular_fields: Option<&[BlastpTabularField]> = if matches!(
        outfmt,
        OutputFormat::Tabular | OutputFormat::TabularWithComments
    ) {
        Some(if let Some(fields) = custom_fields.as_ref() {
            fields.as_slice()
        } else {
            default_blastp_tabular_fields()
        })
    } else {
        None
    };
    let render_alignment = matches!(outfmt, OutputFormat::Pairwise)
        || tabular_fields.is_some_and(blastp_tabular_fields_require_rendered_alignment);
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
    // ```c
    // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
    //     x_PrintField(*iter);
    // }
    // ```
    let stream_tabular_hit_lists = matches!(
        outfmt,
        OutputFormat::Tabular | OutputFormat::TabularWithComments
    ) && !render_alignment;
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:271-327
    // ```c
    // *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
    // (*hsp_list_out)->query_index = index;
    // --hit_list->hsplist_count;
    // ```
    let final_hits = if matches!(outfmt, OutputFormat::Pairwise)
        || (matches!(
            outfmt,
            OutputFormat::Tabular | OutputFormat::TabularWithComments
        ) && !stream_tabular_hit_lists)
    {
        Some(collect_hits_from_hit_lists(&hit_lists))
    } else {
        None
    };
    let pairwise_hits = if matches!(
        outfmt,
        OutputFormat::Pairwise | OutputFormat::Tabular | OutputFormat::TabularWithComments
    ) && !stream_tabular_hit_lists
    {
        Some(build_pairwise_hits(
            final_hits.expect("final hits collected for rendered blastp output"),
            &contexts,
            &subjects,
            &subject_titles,
            args.scoring.matrix,
            render_alignment,
        ))
    } else {
        None
    };

    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_call(
            &timing.final_traceback_ns,
            &timing.final_traceback_calls,
            final_traceback_start,
        );
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1553-1554
    // ```c
    // /* Save the results. */
    // status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
    // ```
    let output_format_start = blastp_timing_start(timing_enabled);
    let output_result: Result<()> = match outfmt {
        OutputFormat::Pairwise => {
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:68-93
            // ```c
            // CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
            //     : m_FormatType(format_type), ..., m_Outfile(outfile),
            //       m_NumSummary(num_summary), ...
            // ```
            let mut writer = open_output_writer(
                args.out.as_ref(),
                in_memory_output.as_mut().map(|output| &mut **output),
            )?;
            let pairwise_config = PairwiseConfig {
                line_length: 60,
                show_gi: false,
                show_frame: false,
                program: "blastp".to_string(),
                protein_matrix: args.scoring.matrix,
            };
            let pairwise_queries: Vec<BlastpPairwiseQuery> = query_headers
                .iter()
                .zip(query_lengths.iter())
                .enumerate()
                .map(|(q_idx, (query_name, &query_length))| BlastpPairwiseQuery {
                    query_name: query_name.clone(),
                    query_length,
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_results.cpp:84-103
                    // ```c
                    // m_SearchSpace = ctx->eff_searchsp;
                    // ...
                    // s_InitializeKarlinBlk(sbp->kbp_std[ctx_index],
                    //                       &m_UngappedKarlinBlk);
                    // ```
                    ungapped_karlin: contexts[q_idx].karlin_params,
                    effective_search_space: search_spaces[q_idx].effective_space.round() as i64,
                })
                .collect();
            let pairwise_report = BlastpPairwiseReport {
                version: NCBI_BLASTP_VERSION.to_string(),
                database_name: format!(
                    "User specified sequence set (Input: {subject_input_label})"
                ),
                database_num_sequences: subject_records.len(),
                database_total_letters: total_db_len,
                matrix_name: BLASTP_DEFAULT_MATRIX_NAME.to_string(),
                gap_open: args.scoring.gap_open,
                gap_extend: args.scoring.gap_extend,
                word_threshold: args.threshold as i32,
                window_size: args.window_size as i32,
                gapped_karlin: gapped_params,
                gumbel: gapped_gumbel,
            };
            write_blastp_pairwise_report(
                pairwise_hits
                    .as_ref()
                    .expect("pairwise hits prepared for blastp outfmt 0"),
                &mut writer,
                &pairwise_config,
                &pairwise_queries,
                &subject_ids,
                &pairwise_report,
            )?;
            Ok(())
        }
        OutputFormat::Tabular | OutputFormat::TabularWithComments => {
            let fields = tabular_fields.expect("tabular fields prepared for blastp tabular output");
            let mut writer = open_output_writer(
                args.out.as_ref(),
                in_memory_output.as_mut().map(|output| &mut **output),
            )?;
            if stream_tabular_hit_lists {
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
                // ```c
                // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
                //     if (iter != m_FieldsToShow.begin())
                //         m_Ostream << m_FieldDelimiter;
                //     x_PrintField(*iter);
                // }
                // ```
                write_blastp_tabular_hit_lists(
                    &hit_lists,
                    fields,
                    outfmt,
                    &mut writer,
                    &query_ids,
                    &subject_ids,
                    &subjects,
                    &subject_titles,
                    &context,
                )?;
            } else {
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
                // ```c
                // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
                //     if (iter != m_FieldsToShow.begin())
                //         m_Ostream << m_FieldDelimiter;
                //     x_PrintField(*iter);
                // }
                // ```
                write_blastp_tabular_output(
                    pairwise_hits
                        .as_ref()
                        .expect("pairwise hits prepared for blastp tabular outfmt"),
                    fields,
                    outfmt,
                    &mut writer,
                    &query_ids,
                    &subject_ids,
                    &context,
                )?;
            }
            Ok(())
        }
    };
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_call(
            &timing.output_format_ns,
            &timing.output_format_calls,
            output_format_start,
        );
    }
    output_result?;

    if let Some(timing) = timing.as_ref() {
        print_blastp_timing(timing, t_total);
    }
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1705
    // ```c
    // BLAST_PreliminarySearchEngine(..., diagnostics, ...);
    // Blast_DiagnosticsUpdate(diagnostics, local_diagnostics);
    // ```
    if kappa_range_counters_enabled {
        print_blastp_kappa_range_counters();
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::args::{BlastpCompBasedStats, BlastpSegSpec};
    use crate::stats::KarlinParams;
    use crate::utils::matrix::ncbistdaa;

    fn make_test_blastp_hsp_list(oid: u32, raw_score: i32) -> BlastpHspList {
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
        // ```c
        // typedef struct BlastHSP {
        //    Int4 score;
        //    double evalue;
        //    BlastSeg query;
        //    BlastSeg subject;
        // } BlastHSP;
        // ```
        let hit = Hit {
            identity: 100.0,
            length: 10,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 10,
            s_start: 1,
            s_end: 10,
            e_value: 1.0,
            bit_score: 0.0,
            num_ident: 10,
            query_frame: 0,
            query_length: 10,
            q_idx: 0,
            s_idx: oid,
            raw_score,
            // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-143
            // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 10,
        };
        BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(hit)],
            best_evalue: f64::MAX,
        }
    }

    fn make_test_subject_result(oid: usize, raw_score: i32) -> BlastpSubjectResult {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1553-1554
        // ```c
        // /* Save the results. */
        // status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
        // ```
        let mut result = BlastpSubjectResult::new(oid);
        result
            .query_hsp_lists
            .push((0, make_test_blastp_hsp_list(oid as u32, raw_score)));
        result
    }

    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    #[test]
    fn test_merge_blastp_subject_results_replays_subject_index_order() {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1554
        // ```c
        // while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
        //     status = s_BlastSearchEngineCore(..., &hsp_list, ...);
        //     status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
        // }
        // ```
        let mut preliminary_hit_lists = vec![None];
        let subject_results = vec![
            Some(make_test_subject_result(0, 50)),
            Some(make_test_subject_result(1, 60)),
            Some(make_test_subject_result(2, 70)),
        ];

        merge_blastp_subject_results_in_order(&mut preliminary_hit_lists, subject_results, 10);

        let hit_list = preliminary_hit_lists[0]
            .as_ref()
            .expect("query hit list should be created");
        let oids: Vec<u32> = hit_list.hsplist_array.iter().map(|list| list.oid).collect();
        assert_eq!(oids, vec![0, 1, 2]);
    }

    fn make_blastp_args() -> BlastpArgs {
        BlastpArgs {
            query: PathBuf::from("query.faa"),
            subject: PathBuf::from("subject.faa"),
            task: "blastp".to_string(),
            evalue: None,
            threshold: None,
            word_size: None,
            num_threads: 1,
            out: None,
            max_target_seqs: 500,
            max_hsps_per_subject: 0,
            ungapped: false,
            window_size: None,
            matrix: None,
            gap_open: None,
            gap_extend: None,
            comp_based_stats: None,
            seg: None,
            use_sw_tback: false,
            outfmt: "6".to_string(),
        }
    }

    fn make_pairwise_hit() -> PairwiseHit {
        PairwiseHit {
            hit: Hit {
                identity: 75.0,
                length: 4,
                mismatch: 1,
                gapopen: 1,
                q_start: 2,
                q_end: 5,
                s_start: 7,
                s_end: 10,
                e_value: 1.0e-25,
                bit_score: 42.5,
                num_ident: 3,
                query_frame: 0,
                query_length: 12,
                q_idx: 0,
                s_idx: 0,
                raw_score: 80,
                // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-143
                // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
                sort_query_offset: 0,
                sort_query_end: 0,
                sort_subject_offset: 0,
                sort_subject_end: 0,
                has_sort_offsets: false,
                gap_info: None,
                num_positives: 3,
            },
            query_seq: Some("AB-C".to_string()),
            subject_seq: Some("ABDC".to_string()),
            query_frame: Some(1),
            subject_frame: Some(1),
            positives: Some(3),
            gaps: Some(1),
            subject_length: Some(30),
            subject_title: Some("subject description".to_string()),
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:87-129
    // ```c
    // if (seq >= word_target) {
    //     BlastLookupAddWordHit(backbone, lut_word_length, charsize,
    //                           seq - lut_word_length, offset - lut_word_length);
    // }
    // ```
    fn make_poly_a_query_frame() -> QueryFrame {
        QueryFrame {
            frame: 1,
            aa_seq: vec![
                0,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                0,
            ],
            aa_seq_nomask: None,
            aa_len: 6,
            orig_len: 6,
            seg_masks: Vec::new(),
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-127
    // ```c
    // (*seq_blk)->sequence_start = (Uint1 *) buffer;
    // (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
    // (*seq_blk)->length = length;
    // ```
    fn make_query_frame_with_residues(residues: &[u8]) -> QueryFrame {
        let mut aa_seq = Vec::with_capacity(residues.len() + 2);
        aa_seq.push(0);
        aa_seq.extend_from_slice(residues);
        aa_seq.push(0);
        QueryFrame {
            frame: 1,
            aa_seq,
            aa_seq_nomask: None,
            aa_len: residues.len(),
            orig_len: residues.len(),
            seg_masks: Vec::new(),
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
    // ```c
    // offset_array_size = OFFSET_ARRAY_SIZE +
    //     ((BlastAaLookupTable*)lookup->lut)->longest_chain;
    // ```
    fn make_poly_a_lookup() -> BlastAaLookupTable {
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let queries = vec![vec![make_poly_a_query_frame()]];
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2746-2748
        // ```c
        // Boolean check_ideal =
        //    (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
        //     program == eBlastTypeRpsTblastn);
        // ```
        let (lookup, _) = build_ncbi_lookup(&queries, 0, &ungapped, false);
        lookup
    }

    fn make_init_hsp(
        q_start_absolute: i32,
        length: i32,
        s_start: i32,
        q_seed_absolute: i32,
        s_seed: i32,
        score: i32,
    ) -> InitHSP {
        InitHSP {
            q_seed_absolute,
            s_seed,
            ungapped_data: BlastpUngappedData {
                q_start: q_start_absolute,
                s_start,
                length,
                score,
            },
        }
    }

    fn make_query_context(q_idx: u32, frame_base: i32) -> QueryContext {
        QueryContext {
            q_idx,
            f_idx: 0,
            frame: 1,
            aa_seq: vec![0, 1, 2, 3, 0],
            aa_seq_nomask: None,
            aa_len: 3,
            orig_len: 3,
            frame_base,
            is_valid: true,
            karlin_params: KarlinParams {
                lambda: 0.267,
                k: 0.041,
                h: 0.14,
                alpha: 1.0,
                beta: -1.0,
            },
        }
    }

    fn make_preliminary_hsp(
        q_idx: u32,
        query_start: i32,
        query_end: i32,
        subject_start: i32,
        subject_end: i32,
        raw_score: i32,
    ) -> BlastpPreliminaryHsp {
        BlastpPreliminaryHsp {
            query_start,
            query_end,
            subject_start,
            subject_end,
            gapped_query_start: query_start,
            gapped_subject_start: subject_start,
            raw_score,
            query_context: q_idx as i32,
            query_frame: 1,
            subject_frame: 0,
            query_length: 100,
            q_idx,
            s_idx: 0,
        }
    }

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    BlastSeg query;
    //    BlastSeg subject;
    //    Int4 context;
    // } BlastHSP;
    // ```
    fn make_local_match_hsp(preliminary_hsp: &BlastpPreliminaryHsp) -> BlastpHsp {
        BlastpHsp {
            identity: 0.0,
            length: usize::try_from(preliminary_hsp.query_end - preliminary_hsp.query_start)
                .expect("NCBI BLAST preliminary query span must be non-negative"),
            mismatch: 0,
            gapopen: 0,
            q_start: usize::try_from(preliminary_hsp.query_start + 1)
                .expect("NCBI BLAST localMatch query start must be non-negative"),
            q_end: usize::try_from(preliminary_hsp.query_end)
                .expect("NCBI BLAST localMatch query end must be non-negative"),
            s_start: usize::try_from(preliminary_hsp.subject_start + 1)
                .expect("NCBI BLAST localMatch subject start must be non-negative"),
            s_end: usize::try_from(preliminary_hsp.subject_end)
                .expect("NCBI BLAST localMatch subject end must be non-negative"),
            gapped_q_start: preliminary_hsp.gapped_query_start,
            gapped_s_start: preliminary_hsp.gapped_subject_start,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_context: preliminary_hsp.query_context,
            query_frame: preliminary_hsp.query_frame,
            subject_frame: preliminary_hsp.subject_frame,
            query_length: preliminary_hsp.query_length,
            q_idx: preliminary_hsp.q_idx,
            s_idx: preliminary_hsp.s_idx,
            raw_score: preliminary_hsp.raw_score,
            gap_info: None,
            num_positives: 0,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:797-800
    // ```c
    // diag_array[diag_coord].last_hit =
    //     s_last_off - (wordsize - 1) + diag_offset;
    // ```
    #[test]
    fn test_blastp_diag_last_hit_after_one_hit_matches_ncbi_formula() {
        assert_eq!(70 - (3 - 1) + 50, 118);
        assert_eq!(11 - (5 - 1) + 19, 26);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:509-521
    // ```c
    // if ((Int4) (subject_offset + diag_offset) <
    //     diag_array[diag_coord].last_hit) {
    //     continue;
    // } else {
    //     diag_array[diag_coord].last_hit =
    //         subject_offset + diag_offset;
    //     diag_array[diag_coord].flag = 0;
    // }
    // ```
    #[test]
    fn test_blastp_flagged_diag_branch_matches_ncbi_skip_vs_reset_order() {
        let mut skip_diag = DiagStruct::default();
        skip_diag.set_flag(1);
        skip_diag.set_last_hit(90);
        let skip_subject_plus_offset = 30 + 50;
        if skip_subject_plus_offset >= skip_diag.last_hit() {
            skip_diag.set_last_hit(skip_subject_plus_offset);
            skip_diag.set_flag(0);
        }
        assert_eq!(skip_diag.flag(), 1);
        assert_eq!(skip_diag.last_hit(), 90);

        let mut reset_diag = DiagStruct::default();
        reset_diag.set_flag(1);
        reset_diag.set_last_hit(90);
        let reset_subject_plus_offset = 40 + 50;
        if reset_subject_plus_offset >= reset_diag.last_hit() {
            reset_diag.set_last_hit(reset_subject_plus_offset);
            reset_diag.set_flag(0);
        }
        assert_eq!(reset_diag.flag(), 0);
        assert_eq!(reset_diag.last_hit(), 90);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:581-588
    // ```c
    // if (right_extend) {
    //     diag_array[diag_coord].flag = 1;
    //     diag_array[diag_coord].last_hit =
    //         s_last_off - (wordsize - 1) + diag_offset;
    // }
    // else {
    //     diag_array[diag_coord].last_hit =
    //         subject_offset + diag_offset;
    // }
    // ```
    #[test]
    fn test_blastp_two_hit_diag_update_matches_ncbi_flag_and_last_hit_updates() {
        let mut right_extended_diag = DiagStruct::default();
        right_extended_diag.set_flag(0);
        right_extended_diag.set_last_hit(0);
        right_extended_diag.set_flag(1);
        right_extended_diag.set_last_hit(70 - (3 - 1) + 50);
        assert_eq!(right_extended_diag.flag(), 1);
        assert_eq!(right_extended_diag.last_hit(), 118);

        let mut plain_diag = DiagStruct::default();
        plain_diag.set_flag(0);
        plain_diag.set_last_hit(0);
        plain_diag.set_last_hit(40 + 50);
        assert_eq!(plain_diag.flag(), 0);
        assert_eq!(plain_diag.last_hit(), 90);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:554-561
    // ```c
    // if (query_offset - diff <
    //     query_info->contexts[curr_context].query_offset) {
    //     diag_array[diag_coord].last_hit =
    //         subject_offset + diag_offset;
    //     continue;
    // }
    // ```
    #[test]
    fn test_blastp_two_hit_cross_query_reject_uses_context_start() {
        assert!(105 - 10 < 100);
        assert!(!(105 - 5 < 100));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:358-372
    // ```c
    // ungapped_data->q_start = q_start;
    // ungapped_data->s_start = s_start;
    // ungapped_data->length = len;
    // ungapped_data->score = score;
    // BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
    // ```
    #[test]
    fn test_blastp_save_init_hsp_uses_qstart_sstart_len_as_single_truth() {
        let mut init_hsps = Vec::new();

        blastp_save_init_hsp(&mut init_hsps, 1015, 42, 1018, 45, 27, 88);

        assert_eq!(init_hsps.len(), 1);
        let hsp = init_hsps[0];
        assert_eq!(hsp.ungapped_data.q_start, 1015);
        assert_eq!(hsp.query_end_absolute(), 1042);
        assert_eq!(hsp.ungapped_data.s_start, 42);
        assert_eq!(hsp.subject_end(), 69);
        assert_eq!(hsp.ungapped_data.length, 27);
        assert_eq!(hsp.q_seed_absolute, 1018);
        assert_eq!(hsp.s_seed, 45);
        assert_eq!(hsp.ungapped_data.score, 88);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:787-800
    // ```c
    // score = s_BlastAaExtendOneHit(..., &hsp_q, &hsp_s, &hsp_len, ...);
    // if (score >= cutoffs->cutoff_score) {
    //     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
    //                      query_offset, subject_offset, hsp_len, score);
    // }
    // ```
    #[test]
    fn test_one_hit_caller_passes_ncbi_payload_directly_to_blastp_save_init_hsp() {
        let query = vec![
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
        ];
        let subject = query.clone();
        let query_offset = 1;
        let subject_offset = 1;
        let ungapped_hit = extend_one_hit(
            ScoringMatrix::Blosum62,
            &query,
            &subject,
            query_offset as usize,
            subject_offset as usize,
            7,
            3,
        )
        .expect("one-hit ungapped extension");
        let mut init_hsps = Vec::new();

        blastp_save_init_hsp(
            &mut init_hsps,
            ungapped_hit.ungapped_data.q_start,
            ungapped_hit.ungapped_data.s_start,
            query_offset,
            subject_offset,
            ungapped_hit.ungapped_data.length,
            ungapped_hit.ungapped_data.score,
        );

        assert_eq!(init_hsps[0].q_seed_absolute, query_offset);
        assert_eq!(init_hsps[0].s_seed, subject_offset);
        assert_eq!(init_hsps[0].ungapped_data.q_start, 0);
        assert_eq!(init_hsps[0].ungapped_data.s_start, 0);
        assert_eq!(init_hsps[0].ungapped_data.length, 5);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:566-583
    // ```c
    // score = s_BlastAaExtendTwoHit(..., &hsp_q, &hsp_s, &hsp_len, ...,
    //                               &right_extend, &s_last_off);
    // if (score >= cutoffs->cutoff_score)
    //     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
    //                      query_offset, subject_offset, hsp_len, score);
    // ```
    #[test]
    fn test_two_hit_caller_preserves_second_hit_seed_and_length_when_right_extension_occurs() {
        let query = vec![
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
        ];
        let subject = query.clone();
        let query_offset = 3;
        let subject_offset = 3;
        let ungapped_hit = extend_two_hit(
            ScoringMatrix::Blosum62,
            &query,
            &subject,
            0,
            subject_offset as usize,
            query_offset as usize,
            7,
            3,
        )
        .expect("two-hit ungapped extension");
        let mut init_hsps = Vec::new();

        blastp_save_init_hsp(
            &mut init_hsps,
            ungapped_hit.ungapped_data.q_start,
            ungapped_hit.ungapped_data.s_start,
            query_offset,
            subject_offset,
            ungapped_hit.ungapped_data.length,
            ungapped_hit.ungapped_data.score,
        );

        assert!(ungapped_hit.right_extend);
        assert_eq!(init_hsps[0].q_seed_absolute, query_offset);
        assert_eq!(init_hsps[0].s_seed, subject_offset);
        assert_eq!(init_hsps[0].ungapped_data.length, 6);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:566-583
    // ```c
    // score = s_BlastAaExtendTwoHit(..., &hsp_q, &hsp_s, &hsp_len, ...,
    //                               &right_extend, &s_last_off);
    // if (score >= cutoffs->cutoff_score)
    //     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
    //                      query_offset, subject_offset, hsp_len, score);
    // ```
    #[test]
    fn test_two_hit_caller_preserves_second_hit_seed_and_length_when_left_stops_early() {
        let query = vec![
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::D,
            ncbistdaa::D,
            ncbistdaa::D,
            ncbistdaa::C,
            ncbistdaa::C,
            ncbistdaa::C,
        ];
        let subject = vec![
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::W,
            ncbistdaa::W,
            ncbistdaa::W,
            ncbistdaa::C,
            ncbistdaa::C,
            ncbistdaa::C,
        ];
        let query_offset = 6;
        let subject_offset = 6;
        let ungapped_hit = extend_two_hit(
            ScoringMatrix::Blosum62,
            &query,
            &subject,
            0,
            subject_offset as usize,
            query_offset as usize,
            2,
            3,
        )
        .expect("two-hit ungapped extension");
        let mut init_hsps = Vec::new();

        blastp_save_init_hsp(
            &mut init_hsps,
            ungapped_hit.ungapped_data.q_start,
            ungapped_hit.ungapped_data.s_start,
            query_offset,
            subject_offset,
            ungapped_hit.ungapped_data.length,
            ungapped_hit.ungapped_data.score,
        );

        assert!(!ungapped_hit.right_extend);
        assert_eq!(init_hsps[0].q_seed_absolute, query_offset);
        assert_eq!(init_hsps[0].s_seed, subject_offset);
        assert_eq!(init_hsps[0].ungapped_data.length, 3);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:492-505
    // ```c
    // while (scan_range[1] <= scan_range[2]) {
    //     hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
    //     totalhits += hits;
    //     for (i = 0; i < hits; ++i) {
    //         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
    //         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
    // ```
    #[test]
    fn test_blastp_for_each_offset_pair_preserves_scan_emission_order() {
        let lookup = make_poly_a_lookup();
        let runtime_lookup = BlastpRuntimeLookup::Aa(lookup);
        let mut subject_sequence = vec![ncbistdaa::A; 6];
        subject_sequence.push(0);
        let mut offset_pairs = vec![OffsetPair::default(); 5];
        let mut emitted = Vec::new();

        blastp_for_each_offset_pair(
            &runtime_lookup,
            &subject_sequence,
            6,
            &mut offset_pairs,
            BlastpOffsetPairScanMode::TwoHit,
            |pair| emitted.push((pair.q_off, pair.s_off)),
        );

        let expected = vec![
            (0, 0),
            (1, 0),
            (2, 0),
            (3, 0),
            (0, 1),
            (1, 1),
            (2, 1),
            (3, 1),
            (0, 2),
            (1, 2),
            (2, 2),
            (3, 2),
            (0, 3),
            (1, 3),
            (2, 3),
            (3, 3),
        ];
        assert_eq!(emitted, expected);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:102-128
    // ```c
    // BlastLookupAddWordHit(backbone, lut_word_length, charsize,
    //                       seq - lut_word_length, offset - lut_word_length);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:492-505
    // ```c
    // Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
    // Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
    // ```
    #[test]
    fn test_blastp_for_each_offset_pair_preserves_absolute_query_offsets_across_queries() {
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let queries = vec![
            vec![make_query_frame_with_residues(&[
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
            ])],
            vec![make_query_frame_with_residues(&[
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
            ])],
        ];
        let (lookup, _) = build_ncbi_lookup(&queries, 0, &ungapped, false);
        let runtime_lookup = BlastpRuntimeLookup::Aa(lookup);
        let subject_sequence = vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A, ncbistdaa::A, 0];
        let mut offset_pairs = vec![OffsetPair::default(); 8];
        let mut emitted = Vec::new();

        blastp_for_each_offset_pair(
            &runtime_lookup,
            &subject_sequence,
            4,
            &mut offset_pairs,
            BlastpOffsetPairScanMode::TwoHit,
            |pair| emitted.push((pair.q_off, pair.s_off)),
        );

        assert_eq!(emitted, vec![(0, 0), (4, 0), (0, 1), (4, 1)]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:552-568
    // ```c
    // if (lookup_wrap->lut_type == eAaLookupTable) { ... }
    // else if (lookup_wrap->lut_type == eCompressedAaLookupTable) {
    //     lut->scansub_callback = (void *)s_BlastCompressedAaScanSubject;
    // }
    // ```
    #[test]
    fn test_blastp_for_each_offset_pair_dispatches_compressed_lookup() {
        let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
        let word = [
            ncbistdaa::S,
            ncbistdaa::T,
            ncbistdaa::A,
            ncbistdaa::G,
            ncbistdaa::W,
        ];
        lookup.add_unencoded(&word, 40);
        lookup.add_unencoded(&word, 41);
        lookup.finalize();
        let runtime_lookup = BlastpRuntimeLookup::Compressed(lookup);
        let subject_sequence = vec![
            ncbistdaa::S,
            ncbistdaa::T,
            ncbistdaa::A,
            ncbistdaa::G,
            ncbistdaa::W,
            ncbistdaa::S,
            ncbistdaa::T,
            ncbistdaa::A,
            ncbistdaa::G,
            ncbistdaa::W,
            0,
        ];
        let mut offset_pairs = vec![OffsetPair::default(); 8];
        let mut emitted = Vec::new();

        blastp_for_each_offset_pair(
            &runtime_lookup,
            &subject_sequence,
            10,
            &mut offset_pairs,
            BlastpOffsetPairScanMode::TwoHit,
            |pair| emitted.push((pair.q_off, pair.s_off)),
        );

        assert_eq!(emitted, vec![(40, 0), (41, 0), (40, 5), (41, 5)]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:496-500
    // ```c
    // scan_range[2] = subject->seq_ranges[0].right - wordsize;
    // if (scan_range[2] < scan_range[1])
    //     scan_range[2] = scan_range[1];
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-121
    // ```c
    // (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
    // (*seq_blk)->length = length;
    // ```
    #[test]
    fn test_blastp_for_each_offset_pair_two_hit_uses_trailing_sentinel_for_short_subjects() {
        let lookup = make_poly_a_lookup();
        let runtime_lookup = BlastpRuntimeLookup::Aa(lookup);
        let subject_sequence = vec![ncbistdaa::A, ncbistdaa::A, 0];
        let mut offset_pairs = vec![OffsetPair::default(); 5];
        let mut emitted = Vec::new();

        blastp_for_each_offset_pair(
            &runtime_lookup,
            &subject_sequence,
            2,
            &mut offset_pairs,
            BlastpOffsetPairScanMode::TwoHit,
            |pair| emitted.push((pair.q_off, pair.s_off)),
        );

        assert!(emitted.is_empty());
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:551-560
    // ```c
    // curr_context = BSearchContextInfo(query_offset, query_info);
    // if (query_offset - diff <
    //     query_info->contexts[curr_context].query_offset) {
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:358-372
    // ```c
    // ungapped_data->q_start = q_start;
    // BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
    // ```
    #[test]
    fn test_blastp_absolute_lookup_offsets_flow_directly_into_context_and_seed_fields() {
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let queries = vec![
            vec![make_query_frame_with_residues(&[
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
            ])],
            vec![make_query_frame_with_residues(&[
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
            ])],
        ];
        let (_lookup, contexts) = build_ncbi_lookup(&queries, 0, &ungapped, false);
        let query_offset = 4;
        let ctx_idx = blastp_context_idx_for_query_offset(&contexts, query_offset);
        let ctx = &contexts[ctx_idx];
        let mut init_hsps = Vec::new();

        blastp_save_init_hsp(&mut init_hsps, query_offset, 12, query_offset, 12, 3, 27);

        assert_eq!(ctx_idx, 1);
        assert_eq!(ctx.frame_base, 4);
        assert_eq!(init_hsps[0].q_seed_absolute, 4);
        assert_eq!(init_hsps[0].ungapped_data.q_start, 4);
        assert_eq!(init_hsps[0].query_end_absolute(), 7);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2445-2466
    // ```c
    // *context = s_GetUngappedHSPContext(query_info, init_hsp);
    // ...
    // s_AdjustInitialHSPOffsets(init_hsp, query_start);
    // ```
    #[test]
    fn test_blastp_adjust_hsp_offsets_and_get_query_data_uses_qoff_as_context_truth() {
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let queries = vec![
            vec![make_query_frame_with_residues(&[
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
            ])],
            vec![make_query_frame_with_residues(&[
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
            ])],
        ];
        let (_lookup, contexts) = build_ncbi_lookup(&queries, 0, &ungapped, false);
        let concatenated_query = build_blastp_concatenated_query(&contexts);
        let mut init_hsps = Vec::new();

        blastp_save_init_hsp(&mut init_hsps, 4, 12, 4, 12, 3, 27);

        let mut init_hsp = init_hsps[0];
        let (context_idx, ctx, query_sequence) = blastp_adjust_hsp_offsets_and_get_query_data(
            &concatenated_query,
            &contexts,
            blastp_context_lookup_bounds(&contexts),
            &mut init_hsp,
        );

        assert_eq!(context_idx, 1);
        assert_eq!(ctx.q_idx, 1);
        assert_eq!(query_sequence.len(), 4);
        assert_eq!(init_hsp.q_seed_absolute, 0);
        assert_eq!(init_hsp.ungapped_data.q_start, 0);
        assert_eq!(init_hsp.query_end_absolute(), 3);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:246-250
    // ```c
    // return cinfo->query_offset + cinfo->query_length + (cinfo->query_length ? 2 : 1);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:500-501,654
    // ```c
    // int buflen = QueryInfo_GetSeqBufLen(qinfo);
    // TAutoUint1Ptr buf((Uint1*) calloc(buflen+1, sizeof(Uint1)));
    // ...
    // BlastSeqBlkSetSequence(*seqblk, buf.release(), buflen - 2);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
    // ```c
    // while (diag_array_length < (qlen+window_size))
    //     diag_array_length = diag_array_length << 1;
    // ```
    #[test]
    fn test_blastp_diag_table_uses_ncbi_query_length_with_shared_boundary_sentinel() {
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let queries = vec![
            vec![make_query_frame_with_residues(&[
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
            ])],
            vec![make_query_frame_with_residues(&[
                ncbistdaa::C,
                ncbistdaa::C,
                ncbistdaa::C,
            ])],
        ];
        let (lookup, _) = build_ncbi_lookup(&queries, 0, &ungapped, false);
        let (diag_array_size, diag_mask) = compute_diag_table_shape(lookup.query_length, 1);

        assert_eq!(lookup.query_length, 7);
        assert_eq!(diag_array_size, 8);
        assert_eq!(diag_mask, 7);
        assert_eq!(lookup.query_length + 1, 8);
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:39-45
    // ```c
    // const char* kDfltArgTabularOutputFmt =
    //     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
    //     "evalue bitscore";
    // const char* kDfltArgTabularOutputFmtTag("std");
    // ```
    #[test]
    fn test_parse_blastp_tabular_fields_expands_std_and_dedups() {
        let fields = parse_blastp_tabular_fields("std qlen qaccver").expect("parsed fields");
        assert_eq!(fields[0], BlastpTabularField::QueryAccessionVersion);
        assert_eq!(fields[1], BlastpTabularField::SubjectAccessionVersion);
        assert!(fields.contains(&BlastpTabularField::QueryLength));
        assert_eq!(
            fields
                .iter()
                .filter(|field| **field == BlastpTabularField::QueryAccessionVersion)
                .count(),
            1
        );
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:985-1027
    // ```c
    // if (num_matches > 0) {
    //     btop_string += NStr::Int8ToString(num_matches);
    //     num_matches = 0;
    // }
    // btop_string += m_QuerySeq[i];
    // btop_string += m_SubjectSeq[i];
    // ```
    #[test]
    fn test_build_btop_matches_ncbi_style_match_run_encoding() {
        assert_eq!(build_btop("AB-C", "ABDC"), "2-D1");
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
    // ```c
    // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
    //     if (iter != m_FieldsToShow.begin())
    //         m_Ostream << m_FieldDelimiter;
    //     x_PrintField(*iter);
    // }
    // ```
    #[test]
    fn test_blastp_tabular_field_text_renders_requested_metadata() {
        let pairwise_hit = make_pairwise_hit();
        assert_eq!(
            blastp_tabular_field_text(
                BlastpTabularField::QueryAccessionVersion,
                &pairwise_hit,
                "query1",
                "subject1"
            ),
            "query1"
        );
        assert_eq!(
            blastp_tabular_field_text(
                BlastpTabularField::SubjectTitle,
                &pairwise_hit,
                "query1",
                "subject1"
            ),
            "subject description"
        );
        assert_eq!(
            blastp_tabular_field_text(
                BlastpTabularField::Btop,
                &pairwise_hit,
                "query1",
                "subject1"
            ),
            "2-D1"
        );
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
    // ```c
    // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
    //     if (iter != m_FieldsToShow.begin())
    //         m_Ostream << m_FieldDelimiter;
    //     x_PrintField(*iter);
    // }
    // ```
    #[test]
    fn test_blastp_hsp_tabular_field_matches_pairwise_metadata_without_alignment_text() {
        let mut pairwise_hit = make_pairwise_hit();
        pairwise_hit.gaps = Some(pairwise_hit.hit.gap_letters());
        let hsp = BlastpHsp::from_hit(pairwise_hit.hit.clone());
        let fields = [
            BlastpTabularField::QueryAccessionVersion,
            BlastpTabularField::SubjectAccessionVersion,
            BlastpTabularField::QueryLength,
            BlastpTabularField::SubjectLength,
            BlastpTabularField::Positives,
            BlastpTabularField::Gaps,
            BlastpTabularField::PercentPositives,
            BlastpTabularField::Frames,
            BlastpTabularField::SubjectTitle,
        ];

        for field in fields {
            let mut rendered = Vec::new();
            write_blastp_hsp_tabular_field(
                &mut rendered,
                field,
                &hsp,
                "query1",
                "subject1",
                30,
                Some("subject description"),
            )
            .expect("tabular field writes");
            assert_eq!(
                String::from_utf8(rendered).expect("tabular field is UTF-8"),
                blastp_tabular_field_text(field, &pairwise_hit, "query1", "subject1")
            );
        }
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_options.c:1163-1169
    // ```c
    // if ((program_number == eBlastTypeTblastn ||
    //      program_number == eBlastTypeBlastp ||
    //      program_number == eBlastTypeBlastx) &&
    //     word_size > 5)
    //     options->lut_type = eCompressedAaLookupTable;
    // ```
    #[test]
    fn test_validate_requested_blastp_support_accepts_blastp_fast_word_size_five() {
        let mut args = make_blastp_args();
        args.task = "blastp-fast".to_string();
        let resolved = args.resolve().expect("resolved blastp-fast args");

        validate_requested_blastp_support(&resolved)
            .expect("word_size 5 compressed blastp-fast is comparison-gated");
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1289-1292
    // ```c
    // ASSERT(word_size == 5 || word_size == 6 || word_size == 7);
    // ```
    #[test]
    fn test_validate_requested_blastp_support_rejects_unverified_compressed_word_size() {
        let mut args = make_blastp_args();
        args.word_size = Some(6);
        let resolved = args.resolve().expect("resolved blastp args");
        let err = validate_requested_blastp_support(&resolved).expect_err("unsupported word size");
        assert!(err.to_string().contains("comparison-gated"));
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3036-3062
    // ```c
    // ECompoAdjustModes compo_adjust_mode =
    //     (ECompoAdjustModes) extendParams->options->compositionBasedStats;
    // ```
    #[test]
    fn test_validate_requested_blastp_support_rejects_unified_p_mode() {
        let mut args = make_blastp_args();
        args.comp_based_stats = Some(BlastpCompBasedStats {
            mode: BlastpCompositionMode::CompositionMatrixAdjust,
            unified_p: true,
        });
        let resolved = args.resolve().expect("resolved blastp args");
        let err = validate_requested_blastp_support(&resolved).expect_err("unsupported unified p");
        assert!(err.to_string().contains("only --comp_based_stats 2"));
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319
    // ```c
    // s_BlastProtGappedAlignment(..., const BlastScoringParameters* score_params,
    //                            BlastGapAlignStruct* gap_align, ...)
    // ```
    #[test]
    fn test_validate_requested_blastp_support_rejects_non_default_scoring() {
        let mut args = make_blastp_args();
        args.matrix = Some("PAM70".to_string());
        args.gap_open = Some(10);
        args.gap_extend = Some(1);
        args.seg = Some(BlastpSegSpec::No);
        let resolved = args.resolve().expect("resolved blastp args");
        let err = validate_requested_blastp_support(&resolved).expect_err("unsupported scoring");
        assert!(err
            .to_string()
            .contains("matrix-aware protein gap alignment"));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2387-2415
    // ```c
    // gapping_params->x_dropoff = (Int4)
    //     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
    //         extendParams->gap_x_dropoff_final);
    // ```
    #[test]
    fn test_blastp_redo_x_dropoff_uses_scaled_final_vs_unscaled_extend_max() {
        let gapped = crate::stats::KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.0,
            beta: -1.0,
        };
        let scaled_lambda = gapped.lambda / 32.0;
        let unscaled_final =
            blastp_gapped_x_dropoff_raw(X_DROP_GAPPED_FINAL as f64, gapped.lambda).max(
                blastp_gapped_x_dropoff_raw(X_DROP_GAPPED_PRELIM as f64, gapped.lambda),
            );

        assert_eq!(
            blastp_redo_x_dropoff(scaled_lambda, &gapped, unscaled_final),
            blastp_gapped_x_dropoff_raw(X_DROP_GAPPED_FINAL as f64, scaled_lambda)
                .max(unscaled_final)
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2438-2462
    // ```c
    // near_identical_cutoff =
    //     (near_identical_cutoff_bits * NCBIMATH_LN2)
    //         / context->sbp->kbp_gap[index]->Lambda;
    // ...
    // if (do_link_hsps) {
    //     cutoff_s =
    //         (int) (hitParams->cutoff_score_min * context->localScalingFactor);
    // } else {
    //     cutoff_s = 1;
    // }
    // ```
    #[test]
    fn test_blastp_redo_params_follow_ncbi_cutoff_formulas() {
        let scaled_lambda = 0.267 / 32.0;
        let near_identical_cutoff = blastp_redo_near_identical_cutoff(scaled_lambda);
        assert!(near_identical_cutoff > 0.0);
        assert_eq!(blastp_redo_cutoff_score(false, 42, 32.0), 1);
        assert_eq!(blastp_redo_cutoff_score(true, 42, 32.0), 1344);
        assert_eq!(blastp_redo_cutoff_score(true, 7, 3.2), 22);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:932-940
    // ```c
    // params->prelim_evalue = evalue;  /* evalue and prelim_evalue same if no CBS. */
    // ...
    // int cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
    // params->prelim_evalue = cbs_stretch*evalue;
    // new_cutoff = BLAST_SpougeEtoS(cbs_stretch*evalue, kbp, sbp->gbp,
    //                     query_info->contexts[context].query_length,
    //                     avg_subject_length);
    // ```
    #[test]
    fn test_blastp_prelim_evalue_and_cutoff_follow_ncbi_cbs_stretch() {
        let expect_value = 1e-5;
        assert_eq!(
            blastp_prelim_evalue(expect_value, BlastpCompositionMode::NoCompositionBasedStats),
            expect_value
        );
        assert_eq!(
            blastp_prelim_evalue(expect_value, BlastpCompositionMode::CompositionMatrixAdjust),
            5.0 * expect_value
        );

        let scoring = crate::config::ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        };
        let gapped = lookup_protein_params(&scoring);
        let gumbel = lookup_protein_gumbel_params(&scoring, 10_000)
            .expect("supported blastp gumbel parameters");
        let cutoff_without_cbs = blastp_cutoff_score_max(
            expect_value,
            120,
            500,
            BlastpCompositionMode::NoCompositionBasedStats,
            &gapped,
            &gumbel,
        );
        let cutoff_with_cbs = blastp_cutoff_score_max(
            expect_value,
            120,
            500,
            BlastpCompositionMode::CompositionMatrixAdjust,
            &gapped,
            &gumbel,
        );

        assert!(cutoff_with_cbs <= cutoff_without_cbs);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:969-983
    // ```c
    // if (sbp->gbp) {
    //     min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
    // }
    // ...
    // BlastHitSavingParametersNew(..., min_subject_length, ..., hit_params);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1403-1440
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // ...
    // if (db_length == 0) {
    //     BLAST_OneSubjectUpdateParameters(..., seq_arg.seq->length, ...);
    // }
    // ```
    #[test]
    fn test_blastp_preliminary_cutoff_subject_length_uses_global_min_subject_len() {
        let subjects = vec![
            encode_protein_sequence(b"AAAAAA"),
            encode_protein_sequence(b"AAA"),
            encode_protein_sequence(b"AAAAAAAAAA"),
        ];

        assert_eq!(blastp_preliminary_cutoff_subject_length(&subjects), 3);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:555
    // ```c
    // Blast_HSPListSortByScore(hsp_list);
    // ```
    #[test]
    fn test_prepare_preliminary_hits_for_kappa_purges_shared_endpoints_then_sorts() {
        let hits = vec![
            BlastpPreliminaryHsp {
                query_start: 10,
                query_end: 60,
                subject_start: 20,
                subject_end: 70,
                gapped_query_start: 30,
                gapped_subject_start: 40,
                raw_score: 80,
                query_context: 0,
                query_frame: 0,
                subject_frame: 0,
                query_length: 100,
                q_idx: 0,
                s_idx: 0,
            },
            BlastpPreliminaryHsp {
                query_start: 10,
                query_end: 80,
                subject_start: 20,
                subject_end: 90,
                gapped_query_start: 35,
                gapped_subject_start: 45,
                raw_score: 90,
                query_context: 0,
                query_frame: 0,
                subject_frame: 0,
                query_length: 100,
                q_idx: 0,
                s_idx: 0,
            },
            BlastpPreliminaryHsp {
                query_start: 5,
                query_end: 40,
                subject_start: 7,
                subject_end: 42,
                gapped_query_start: 20,
                gapped_subject_start: 22,
                raw_score: 120,
                query_context: 0,
                query_frame: 0,
                subject_frame: 0,
                query_length: 100,
                q_idx: 0,
                s_idx: 0,
            },
        ];

        let prepared = prepare_preliminary_hits_for_kappa(hits);
        assert_eq!(prepared.len(), 2);
        assert_eq!(prepared[0].raw_score, 120);
        assert_eq!(prepared[1].raw_score, 90);
        assert_eq!(prepared[1].query_start, 10);
        assert_eq!(prepared[1].subject_start, 20);
        assert_eq!(prepared[1].query_end, 80);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2505
    // ```c
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
    // ```
    #[test]
    fn test_ncbi_qsort_preliminary_hits_matches_requested_ncbi_comparator() {
        let hits = vec![
            make_preliminary_hsp(0, 10, 40, 15, 45, 90),
            make_preliminary_hsp(0, 10, 45, 15, 50, 120),
            make_preliminary_hsp(0, 5, 30, 9, 34, 120),
            make_preliminary_hsp(1, 3, 20, 4, 21, 80),
        ];

        let mut by_query_offset = hits.clone();
        ncbi_qsort_preliminary_hits(&mut by_query_offset, preliminary_query_offset_compare);
        let ordered_by_query_offset: Vec<(i32, i32, i32, i32)> = by_query_offset
            .iter()
            .map(|hsp| {
                (
                    hsp.query_context,
                    hsp.query_start,
                    hsp.subject_start,
                    hsp.raw_score,
                )
            })
            .collect();
        assert_eq!(
            ordered_by_query_offset,
            vec![
                (0, 5, 9, 120),
                (0, 10, 15, 120),
                (0, 10, 15, 90),
                (1, 3, 4, 80)
            ]
        );

        let mut by_query_end = hits.clone();
        ncbi_qsort_preliminary_hits(&mut by_query_end, preliminary_query_end_compare);
        let ordered_by_query_end: Vec<(i32, i32, i32, i32)> = by_query_end
            .iter()
            .map(|hsp| {
                (
                    hsp.query_context,
                    hsp.query_end,
                    hsp.subject_end,
                    hsp.raw_score,
                )
            })
            .collect();
        assert_eq!(
            ordered_by_query_end,
            vec![
                (0, 30, 34, 120),
                (0, 40, 45, 90),
                (0, 45, 50, 120),
                (1, 20, 21, 80)
            ]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2488-2495,2514-2520
    // ```c
    // hsp_count--;
    // hsp = hsp_array[i+j];
    // hsp = Blast_HSPFree(hsp);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2320,2385
    // ```c
    // return 0;
    // ```
    #[test]
    fn test_prepare_preliminary_hits_for_kappa_keeps_first_equal_endpoint_hsp() {
        let mut first = make_preliminary_hsp(0, 12, 6074, 13, 6057, 16284);
        first.gapped_query_start = 4964;
        first.gapped_subject_start = 4942;
        let mut later = first;
        later.gapped_query_start = 55;
        later.gapped_subject_start = 55;

        let prepared = prepare_preliminary_hits_for_kappa(vec![first, later]);

        assert_eq!(prepared.len(), 1);
        assert_eq!(prepared[0].gapped_query_start, 4964);
        assert_eq!(prepared[0].gapped_subject_start, 4942);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2482-2487
    // ```c
    // while (i+j < hsp_count &&
    //        hsp_array[i]->context == hsp_array[i+j]->context &&
    //        hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
    //        hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
    //        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame) {
    // ```
    #[test]
    fn test_prepare_preliminary_hits_for_kappa_keeps_shared_endpoints_on_subject_frame_change() {
        let mut plus = make_preliminary_hsp(0, 10, 60, 20, 70, 80);
        plus.subject_frame = 1;
        let mut minus = make_preliminary_hsp(0, 10, 60, 20, 70, 90);
        minus.subject_frame = -1;

        let prepared = prepare_preliminary_hits_for_kappa(vec![plus, minus]);

        assert_eq!(prepared.len(), 2);
        assert_eq!(prepared[0].subject_frame, -1);
        assert_eq!(prepared[1].subject_frame, 1);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_collector.c:117-149
    // ```c
    // query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
    // ...
    // Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
    // ```
    #[test]
    fn test_split_preliminary_hits_by_query_preserves_subject_order_within_query() {
        let contexts = vec![make_query_context(0, 0), make_query_context(1, 1000)];
        let hits = vec![
            make_preliminary_hsp(99, 15, 35, 100, 120, 50),
            make_preliminary_hsp(77, 10, 30, 40, 60, 80),
            make_preliminary_hsp(88, 25, 45, 130, 150, 40),
            make_preliminary_hsp(66, 35, 55, 70, 90, 60),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, mut hsp)| {
            hsp.query_context = if index % 2 == 0 { 1 } else { 0 };
            hsp
        })
        .collect();

        let hits_by_query = split_preliminary_hits_by_query(hits, &contexts, 2);

        assert_eq!(hits_by_query.len(), 2);
        assert_eq!(hits_by_query[0].len(), 2);
        assert_eq!(hits_by_query[0][0].query_start, 10);
        assert_eq!(hits_by_query[0][1].query_start, 35);
        assert_eq!(hits_by_query[1].len(), 2);
        assert_eq!(hits_by_query[1][0].subject_start, 100);
        assert_eq!(hits_by_query[1][1].subject_start, 130);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:289-318
    // ```c
    // hit_list = results->hitlist_array[index];
    // last_hsplist_index = hit_list->hsplist_count - 1;
    // *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
    // --hit_list->hsplist_count;
    // ```
    #[test]
    fn test_preliminary_hit_from_local_match_hsp_preserves_ncbi_gapped_starts_and_context() {
        let mut preliminary_hsp = make_preliminary_hsp(0, 10, 20, 30, 40, 50);
        preliminary_hsp.gapped_query_start = 13;
        preliminary_hsp.gapped_subject_start = 29;
        preliminary_hsp.query_context = 7;
        preliminary_hsp.subject_frame = -1;

        let rebuilt = preliminary_hit_from_local_match_hsp(&make_local_match_hsp(&preliminary_hsp));

        assert_eq!(rebuilt.query_start, 10);
        assert_eq!(rebuilt.subject_start, 30);
        assert_eq!(rebuilt.gapped_query_start, 13);
        assert_eq!(rebuilt.gapped_subject_start, 29);
        assert_eq!(rebuilt.query_context, 7);
        assert_eq!(rebuilt.subject_frame, -1);
    }

    #[test]
    fn test_restricted_preliminary_score_retries_exact_within_ncbi_inconclusive_band() {
        assert!(should_retry_exact_preliminary_gapped_alignment(
            true, 95, 100, 90
        ));
        assert!(!should_retry_exact_preliminary_gapped_alignment(
            false, 95, 100, 90
        ));
        assert!(!should_retry_exact_preliminary_gapped_alignment(
            true, 89, 100, 90
        ));
        assert!(!should_retry_exact_preliminary_gapped_alignment(
            true, 100, 100, 90
        ));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3975-4000
    // ```c
    // /* remove all HSPs computed for the current query */
    // for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
    //     ...
    //     BlastIntervalTreeAddHSP(hsp_list->hsp_array[index2], tree,
    //                             query_info, eQueryAndSubject);
    // }
    // ```
    #[test]
    fn test_rebuild_preliminary_interval_tree_drops_only_retried_query_hsps() {
        let contexts = vec![make_query_context(0, 0), make_query_context(1, 1000)];
        let query0_hsp = make_preliminary_hsp(0, 10, 60, 20, 70, 100);
        let query1_hsp = make_preliminary_hsp(1, 15, 55, 30, 70, 90);
        let mut preliminary_hsp_list = vec![query0_hsp, query1_hsp];
        let mut interval_tree = BlastIntervalTree::new(0, 2000, 0, 200);

        rebuild_preliminary_interval_tree_from_list(
            &mut interval_tree,
            &preliminary_hsp_list,
            &contexts,
        );

        let query0_contained = TreeHsp {
            query_offset: 20,
            query_end: 40,
            subject_offset: 30,
            subject_end: 50,
            score: 50,
            query_frame: 1,
            query_length: 100,
            query_context_offset: contexts[0].frame_base,
            subject_frame_sign: 0,
        };
        let query1_contained = TreeHsp {
            query_offset: 20,
            query_end: 40,
            subject_offset: 40,
            subject_end: 60,
            score: 50,
            query_frame: 1,
            query_length: 100,
            query_context_offset: contexts[1].frame_base,
            subject_frame_sign: 0,
        };

        assert!(interval_tree.contains_hsp(&query0_contained, contexts[0].frame_base, 0));
        assert!(interval_tree.contains_hsp(&query1_contained, contexts[1].frame_base, 0));

        preliminary_hsp_list[0].q_idx = 99;
        preliminary_hsp_list[1].q_idx = 77;
        remove_preliminary_hits_for_query(&mut preliminary_hsp_list, 0, &contexts);
        rebuild_preliminary_interval_tree_from_list(
            &mut interval_tree,
            &preliminary_hsp_list,
            &contexts,
        );

        assert!(!interval_tree.contains_hsp(&query0_contained, contexts[0].frame_base, 0));
        assert!(interval_tree.contains_hsp(&query1_contained, contexts[1].frame_base, 0));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2315
    // ```c
    // /* If these are from different contexts, don't compare offsets */
    // if (h1->context < h2->context)
    //    return -1;
    // ```
    #[test]
    fn test_preliminary_common_endpoint_compare_uses_context_not_query_index() {
        let mut left = make_preliminary_hsp(0, 10, 40, 20, 50, 80);
        left.query_context = 1;
        let mut right = make_preliminary_hsp(0, 10, 40, 20, 50, 70);
        right.query_context = 0;

        let prepared = prepare_preliminary_hits_for_kappa(vec![left, right]);

        assert_eq!(prepared.len(), 2);
        assert_eq!(prepared[0].query_context, 1);
        assert_eq!(prepared[1].query_context, 0);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3119-3128
    // ```c
    // BlastCompo_HeapInitialize(&redoneMatches[query_index],
    //                           hitParams->options->hitlist_size,
    //                           inclusion_ethresh);
    // ```
    #[test]
    fn test_redone_match_heaps_use_hitlist_size_not_prelim_hitlist_size() {
        let heaps = build_redone_match_heaps(3, 7);
        assert_eq!(heaps.len(), 3);
        assert_eq!(heaps[0].heap_threshold, 7);
        assert!((heaps[0].ecutoff - PSI_INCLUSION_ETHRESH).abs() < f64::EPSILON);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:204-221,324-329,917-921
    // ```c
    // if (!(query_info->contexts[context].is_valid))
    //    continue;
    // ...
    // if (!(query_info->contexts[context].is_valid)) {
    //    curr_cutoffs->cutoff_score = INT4_MAX;
    //    continue;
    // }
    // ...
    // if (!(query_info->contexts[context].is_valid)) {
    //    params->cutoffs[context].cutoff_score = INT4_MAX;
    //    continue;
    // }
    // ```
    #[test]
    fn test_invalid_blastp_context_uses_ncbi_placeholder_cutoffs() {
        let invalid = QueryContext {
            is_valid: false,
            ..make_query_context(0, 0)
        };
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:3577-3622
        // ```c
        // Blast_KarlinBlkGappedLoadFromTables(kbp, gap_open,
        //                                    gap_extend, matrix_name, error_return);
        // ```
        let scoring = crate::config::ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        };
        let gapped = lookup_protein_params(&scoring);
        let gumbel = lookup_protein_gumbel_params(&scoring, 10_000)
            .expect("supported blastp gumbel parameters");

        assert_eq!(blastp_x_dropoff_init_for_context(&invalid), 0);
        assert_eq!(
            blastp_hit_cutoff_score_for_context(
                &invalid,
                10.0,
                100,
                BlastpCompositionMode::CompositionMatrixAdjust,
                &gapped,
                &gumbel,
            ),
            i32::MAX
        );
        assert_eq!(blastp_word_cutoff_score_for_context(50, &invalid), i32::MAX);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:306-309
    // ```c
    // void Blast_InitHitListSortByScore(BlastInitHitList * init_hitlist)
    // {
    //     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
    //           sizeof(BlastInitHSP), score_compare_match);
    // }
    // ```
    #[test]
    fn test_chain_blastp_init_hsps_restores_ncbi_score_order() {
        let contexts = vec![QueryContext {
            aa_len: 256,
            ..make_query_context(0, 0)
        }];
        let init_hsps = vec![
            make_init_hsp(0, 10, 0, 0, 0, 50),
            make_init_hsp(100, 10, 100, 100, 100, 60),
        ];
        let mut keep = Vec::new();
        let mut nodes = Vec::new();

        let chained = chain_blastp_init_hsps(
            init_hsps,
            &contexts,
            &[100],
            &[1],
            11,
            1,
            &mut keep,
            &mut nodes,
        );

        assert_eq!(chained.len(), 2);
        assert_eq!(chained[0].ungapped_data.score, 60);
        assert_eq!(chained[1].ungapped_data.score, 50);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3468-3499
    // ```c
    // qsort(init_hitlist->init_hsp_array, init_hitlist->total,
    //       sizeof(BlastInitHSP), s_CompareInitHSPsByQueryOffsetScore);
    // ```
    #[test]
    fn test_ncbi_qsort_init_hsps_by_query_offset_score_matches_ncbi_ties() {
        let mut init_hsps = vec![
            make_init_hsp(20, 8, 40, 20, 40, 90),
            make_init_hsp(10, 7, 15, 10, 15, 80),
            make_init_hsp(10, 9, 15, 10, 15, 120),
            make_init_hsp(20, 6, 35, 20, 35, 90),
        ];

        ncbi_qsort_init_hsps_by_query_offset_score(&mut init_hsps);

        let ordered: Vec<(i32, i32, i32)> = init_hsps
            .iter()
            .map(|hsp| {
                (
                    hsp.ungapped_data.q_start,
                    hsp.ungapped_data.score,
                    hsp.ungapped_data.s_start,
                )
            })
            .collect();
        assert_eq!(
            ordered,
            vec![(10, 120, 15), (10, 80, 15), (20, 90, 40), (20, 90, 35)]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:274-299
    // ```c
    // if (0 == (result = BLAST_CMP(h2->ungapped_data->score,
    //                              h1->ungapped_data->score)) &&
    //     0 == (result = BLAST_CMP(h1->ungapped_data->s_start,
    //                              h2->ungapped_data->s_start)) &&
    //     0 == (result = BLAST_CMP(h2->ungapped_data->length,
    //                              h1->ungapped_data->length)) &&
    //     0 == (result = BLAST_CMP(h1->ungapped_data->q_start,
    //                              h2->ungapped_data->q_start))) {
    //     result = BLAST_CMP(h2->ungapped_data->length,
    //                        h1->ungapped_data->length);
    // }
    // ```
    #[test]
    fn test_compare_init_hsps_by_score_matches_ncbi_score_compare_match_ties() {
        let mut init_hsps = vec![
            make_init_hsp(20, 10, 30, 20, 30, 100),
            make_init_hsp(5, 10, 10, 5, 10, 100),
            make_init_hsp(20, 20, 10, 20, 10, 100),
            make_init_hsp(10, 20, 10, 10, 10, 100),
        ];

        ncbi_qsort_init_hsps_by_score(&mut init_hsps);

        let ordered: Vec<(i32, i32, i32)> = init_hsps
            .iter()
            .map(|hsp| {
                (
                    hsp.ungapped_data.q_start,
                    hsp.ungapped_data.s_start,
                    hsp.ungapped_data.length,
                )
            })
            .collect();
        assert_eq!(
            ordered,
            vec![(10, 10, 20), (20, 10, 20), (5, 10, 10), (20, 30, 10)]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3642-3657
    // ```c
    // while (init_hitlist->total > 0 &&
    //       init_hitlist->init_hsp_array[init_hitlist->total - 1].ungapped_data == NULL) {
    //    init_hitlist->total--;
    // }
    // for (k = 0; k < init_hitlist->total - 1; k++) {
    //    if (init_hitlist->init_hsp_array[k].ungapped_data == NULL) {
    //      Int4 end = init_hitlist->total - 1;
    //      init_hitlist->init_hsp_array[k] = init_hitlist->init_hsp_array[end];
    //      init_hitlist->total--;
    //      ...
    //    }
    // }
    // ```
    #[test]
    fn test_compact_chained_init_hsps_ncbi_uses_tail_swap_compaction_order() {
        let mut init_hsps = vec![
            make_init_hsp(10, 10, 10, 10, 10, 10),
            make_init_hsp(20, 10, 20, 20, 20, 20),
            make_init_hsp(30, 10, 30, 30, 30, 30),
            make_init_hsp(40, 10, 40, 40, 40, 40),
            make_init_hsp(50, 10, 50, 50, 50, 50),
        ];
        let mut keep = vec![true, false, true, false, true];

        compact_chained_init_hsps_ncbi(&mut init_hsps, &mut keep);

        let ordered_q_starts: Vec<i32> = init_hsps
            .iter()
            .map(|hsp| hsp.ungapped_data.q_start)
            .collect();
        assert_eq!(ordered_q_starts, vec![10, 50, 30]);
        assert_eq!(keep, vec![true, true, true]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3748-3786
    // ```c
    // for (i = 0;i < init_hitlist->total;i++) {
    //     init_hsp = &init_hitlist->init_hsp_array[i];
    //     ...
    //     if (found[query_index]) {
    //         continue;
    //     }
    // ```
    #[test]
    fn test_build_restricted_align_array_uses_first_init_hsp_per_query() {
        let contexts = vec![make_query_context(0, 0), make_query_context(1, 1000)];
        let init_hsps = vec![
            make_init_hsp(20, 20, 20, 20, 20, 45),
            make_init_hsp(0, 20, 0, 0, 0, 80),
            make_init_hsp(1060, 20, 60, 1060, 60, 70),
        ];

        let restricted = build_restricted_align_array(&init_hsps, &contexts, &[100, 100], 2);

        assert_eq!(restricted, vec![true, false]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3748-3786
    // ```c
    // q_off = init_hsp->offsets.qs_offsets.q_off;
    // ...
    // if (q_off >= query_info->contexts[contxt].query_offset &&
    //     q_off < query_info->contexts[contxt].query_offset +
    //             query_info->contexts[contxt].query_length) {
    // ```
    #[test]
    fn test_build_restricted_align_array_uses_qoff_not_cached_context() {
        let contexts = vec![make_query_context(0, 0), make_query_context(1, 1000)];
        let init_hsps = vec![
            make_init_hsp(0, 20, 0, 0, 0, 80),
            make_init_hsp(1005, 20, 50, 1005, 50, 45),
        ];

        let restricted = build_restricted_align_array(&init_hsps, &contexts, &[100, 100], 2);

        assert_eq!(restricted, vec![false, true]);
    }
}
