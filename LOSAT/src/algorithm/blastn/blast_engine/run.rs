//! Main BLASTN run function
//!
//! This module contains the main `run()` function that coordinates the BLASTN
//! search process.

use crate::api::local_blast::{FormatOutput, FormatProbe, OutputSink, ReportOutputs};
use crate::common::GapEditOp;
use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
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
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1674-1783
// ```c
// BLAST_PreliminarySearchEngine(..., hsp_stream, ...);
// ...
// BLAST_ComputeTraceback(..., hsp_stream, ..., results, ...);
// ```
#[cfg(all(
    feature = "parallel",
    not(all(target_arch = "wasm32", feature = "wasm-threads"))
))]
use std::sync::mpsc::channel;
use std::sync::Arc;

use crate::config::NuclScoringSpec;
use crate::core::blast_encoding::{
    encode_iupac_to_blastna, encode_iupac_to_blastna_into, encode_subject_ncbi2na_packed,
    COMPRESSION_RATIO,
};
use crate::stats::length_adjustment::compute_length_adjustment_ncbi;

// Import from parent blastn module (super::super:: because we're in blast_engine/)
use super::super::alignment::{
    blast_get_offsets_for_gapped_alignment,
    blast_get_start_for_gapped_alignment_nucl,
    build_blastna_matrix,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_scratch,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_traceback_with_scratch,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936 (BLAST_GreedyGappedAlignment)
    // (the score-only call goes through `x_greedy_prelim`, on the one-byte copy of the subject)
    greedy_gapped_alignment_with_traceback,
    stats_from_edit_ops,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    GapAlignScratch,
    GreedyAlignScratch,
};
// No NCBI counterpart: dispatch of the preliminary greedy alignment
// (NCBI reference: blast_gapalign.c:4019-4025, BLAST_GreedyGappedAlignment with compressed_subject).
use super::super::alignment::x_fence;
use super::super::alignment::x_greedy_prelim;
use super::super::archive;
use super::super::args::{BlastnArgs, QueryStrand};
use super::super::blast_extend::DiagStruct;
use super::super::coordination::{
    build_lookup_tables, chunk_query_masks, collect_lowercase_masks, configure_task,
    finalize_task_config, prepare_sequence_data, query_masks, reverse_mask_intervals,
    subject_metadata_from_records, subject_trace_id, LookupTables, QueryLetters, SequenceData,
    SubjectMetadata,
};
use super::super::disc_lookup::{choose_disc_scan_subject, disc_word_scan_subject};
use super::super::extension::{
    build_compressed_query, build_nucl_score_table, extend_hit_ungapped_approx_ncbi,
    extend_hit_ungapped_exact_ncbi, type_of_word, SmallNaWord,
};
use super::super::filtering::{
    blast_hsp_test_identity_and_length, hsp_test, purge_hsps_with_common_endpoints,
    purge_hsps_with_common_endpoints_ex, reevaluate_hsp_with_ambiguities_gapped_ex,
    subject_best_hit, subject_best_hit_by, BestHitKey, ReevalParams,
};
use super::super::hsp::{
    evalue_comp, get_prelim_hitlist_size, score_compare_hsps, sort_hsplist_by_evalue,
    sort_hsps_by_score, trim_by_max_hsps, update_best_evalue, BlastnHitList, BlastnHsp,
    BlastnHspList, BlastnOutputFormat, HitList, HitListEntry, NCBI_BLASTN_VERSION,
};
use super::super::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use super::super::lookup::{
    build_unmasked_ranges, encode_reverse_complement_blastna_into, DiscTemplates,
};
use super::super::lookup::{x_dc_csr, x_mb_lut_keep};
use super::super::pairwise::{
    pairwise_hits, pairwise_windows, shown_query_masks, DisplayMasks, Rows as PairwiseRows,
};
use super::super::query_split::{
    calculate_num_chunks, restrict_masks, split_query_batch, QueryChunk, SplitSizes,
};
use super::super::scoring::{
    check_greedy_gap_costs, check_losat_limits, check_scoring_options, context_blocks,
    context_blocks_ungapped, context_ungapped_blocks, gap_x_dropoffs, karlin_error, ContextKarlin,
};
use super::super::seqalign::write_blastn_seqalign;
use super::super::ungapped;
use crate::blastinput::app::{self as app, Program};
use crate::core::blast_hits::reap_by_query_coverage;
use crate::report::record::{QueryNames, SubjectNames};
// EXPERIMENT (LOSAT_X_BNSTORE): the subjects' letters read from the records or from packed subjects.
use super::super::x_subject_store::{
    x_blastna_into, x_bnstore, x_read_subject_stores, SubjectLetters,
};
use crate::blastinput::fasta_reader::{
    read_queries, read_subjects, FastaInputSource, FastaRecord, QueryEnd, QueryRecords, ReadError,
    ReaderConfig,
};
use crate::blastinput::input_files::{
    check_utf8_file_name, is_blank, open_input, write_empty_subject_warnings_by_length,
};
use crate::blastinput::query_batch::{next_query_batch_end, BatchSizeMixer};
use crate::blastinput::seq_range::{
    cut_queries, cut_subjects, next_ranged_batch, parse_optional_range, subjects_read, Placements,
    QueryInput, RangeRole, SequenceRange,
};
use crate::report::pairwise::{
    write_blastn_pairwise_prolog, write_blastn_pairwise_report, BlastnPairwiseQuery,
    BlastnPairwiseReport, PairwiseConfig, PairwiseDatabase, PairwiseHit,
};
use crate::report::query_warnings::query_warning;
use crate::sequence::packed_subject::PackedNucleotide;
use crate::stats::KarlinParams;
// NCBI reference: c++/src/algo/blast/core/blast_parameters.c:342-344,370-374
// gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda);
// new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
use super::super::ncbi_cutoffs::{
    cutoff_score_for_ungapped_extension, cutoff_score_max_from_evalue, gap_trigger_raw_score,
    GAP_TRIGGER_BIT_SCORE_NUCL,
};
use super::super::tracing as blastn_trace;
use crate::utils::dust::MaskedInterval;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1676
// ```c
//     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
//
//         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
// ```
// EXPERIMENT (PAIRPAR, graduated on 2026-10-22: no variable): the one-hit seed stage of one subject
// chunk scanned in strips by the pool threads and folded in scan order by the owner (a child
// module, so that it uses the private types and kernels of this file without changing them).
#[path = "x_pair_seed.rs"]
mod x_pair_seed;

// EXPERIMENT (LOSAT_X_BNFEWSUBJ): a `-subject` search whose longest 5 Mb unit is a large part of
// its letters searches its subjects, and a batch of chunks the same way its chunks, one at a time
// on the whole pool, so that the seed stage scans each chunk in strips (`x_pair_seed`).
#[path = "x_few_subj.rs"]
mod x_few_subj;

// EXPERIMENT (LOSAT_X_BNSTREAM): the subjects of a subject set searched by the pool threads in
// OID order with a bounded window of results that the owner hands over in OID order.
#[path = "x_bnstream.rs"]
mod x_bnstream;

// EXPERIMENT (LOSAT_X_BNLUTDROP): a query batch's lookup tables dropped when its preliminary stage
// ends, before its traceback, as NCBI's traceback stage drops its table.
#[path = "x_lut_drop.rs"]
mod x_lut_drop;

// EXPERIMENT (LOSAT_X_BNQDROP): the letters of a query batch's records dropped once the batch has
// encoded them, where no report reads them.
#[path = "x_query_drop.rs"]
mod x_query_drop;

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:310-318
// ```c++
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
//     ...
//         BLAST_PROF_START( APP.POST );
//         formatter.PrintEpilog(opt);
// ```
// EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of each query batch written when the batch loop of
// `run_in_pool` takes the batch, and the epilogs after the loop (a child module, so that it reads
// the private types of this file).
#[path = "x_rpt_stream.rs"]
mod x_rpt_stream;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
// ```c++
//     if (args.Exist(kArgDb) && args[kArgDb]) {
// 	std::string local_dblist = NStr::TruncateSpaces( args[kArgDb].AsString() );
//
//         m_SearchDb.Reset(new CSearchDatabase( local_dblist,
//                                              mol_type));
// ```
// The search of an NCBI BLAST database (`-db`): its reading and its drivers of the
// per-subject kernel of this file (a child module, so that it uses the kernel's private
// types).
#[path = "database.rs"]
mod database;

// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c (whole file)
// BLASTN `-culling_limit`: the culling writer of the preliminary stage and the culling pipe of
// the traceback (a child module, so that it uses the kernel's private types).
#[path = "culling.rs"]
mod culling;

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-300
// ```c++
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//     ...
//                 results = lcl_blast.Run();
//                 if (!batch_size)
//                     input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()));
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-603 (the query buffer
// of a batch) and c++/src/algo/blast/core/blast_stat.c:2754-2784 (the ungapped Karlin block of
// each context).
// EXPERIMENT (LOSAT_X_QPREPPAR): both steps of `search_query_batch` as ordered parallel maps on
// the batch's search pool (a child module, so that it reads the private serial function and the
// sentinel of this file).
#[path = "x_query_prep_par.rs"]
mod x_query_prep_par;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850 and
// c++/src/algo/blast/core/blast_traceback.c:375 (the two loops over one subject's HSPs).
// EXPERIMENT (LOSAT_X_QGROUPPAR): both loops of one subject run per query on the batch's search
// pool (a child module, beside the other modules of `search_query_batch`).
#[path = "x_qgroup_par.rs"]
mod x_qgroup_par;

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
// ```c++
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
// 	    BLAST_PROF_START( APP.LOOP.PRE );
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
// ```
// EXPERIMENT (LOSAT_X_QSTREAM): the records of each query batch read from the query file when the
// batch loop of `run_in_pool` needs them, after a first reading that keeps only each record's
// length (a child module, so that it reads the private types of this file).
#[path = "x_query_stream.rs"]
mod x_query_stream;

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279 (the batch reading that
// LOSAT_X_QSTREAM's first reading precedes).
// EXPERIMENT (LOSAT_X_QSCANPAR): that first reading cut into parts at record starts and read by
// the search's threads (a child module, beside `x_query_stream`, which calls it).
#[path = "x_query_scan_par.rs"]
mod x_query_scan_par;

/// EXPERIMENT (LOSAT_X_QPREPPAR): the summary line `main` prints at exit. No NCBI
/// counterpart: it prints counters and changes no value NCBI computes.
pub fn x_query_prep_par_print_summary() {
    x_query_prep_par::print_summary();
}

/// EXPERIMENT (LOSAT_X_QGROUPPAR, LOSAT_X_STATS): the summary lines `main` prints at exit.
/// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub fn x_qgroup_par_print_summary() {
    x_qgroup_par::print_summary();
}

/// EXPERIMENT (LOSAT_X_BNTBRANGELOG): the summary line `main` prints at exit. No NCBI
/// counterpart: it prints counters and changes no value NCBI computes.
pub fn x_tb_range_print_summary() {
    database::x_tb_range_print_summary();
}

/// EXPERIMENT (LOSAT_X_BNLUTDROP): the summary line `main` prints at exit (LOSAT_TIMING). No NCBI
/// counterpart: it prints counters and changes no value NCBI computes.
pub fn x_lut_drop_print_summary() {
    x_lut_drop::print_summary();
}

/// EXPERIMENT (LOSAT_X_BNQDROP): the summary line `main` prints at exit (LOSAT_TIMING). No NCBI
/// counterpart: it prints counters and changes no value NCBI computes.
pub fn x_query_drop_print_summary() {
    x_query_drop::print_summary();
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/ncbi_math.h:160-161
// ```c
// #define NCBIMATH_LN2 0.69314718055994530941723212145818
// ```
const NCBIMATH_LN2: f64 = 0.69314718055994530941723212145818;

/// The text of a query's outfmt 0 `Query=` line and outfmt 7 `# Query:` line: the record's
/// title bytes (the local ID is not shown without `-parse_deflines`; the writers apply
/// `TruncateSpaces`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:732-734
/// ```c++
///     string all_id_str = GetSeqIdString(cbs, believe_query);
///     all_id_str += " ";
///     all_id_str = NStr::TruncateSpaces(all_id_str + GetSeqDescrString(cbs));
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:633-635
/// ```c++
///     CRef<CSeq_id> wid = FindBestChoice(ids, CSeq_id::WorstRank);
///
///     if (wid && (wid->Which()!= CSeq_id::e_Local || believe_local_id)){
/// ```
#[inline]
fn record_title(record: &FastaRecord) -> Arc<[u8]> {
    Arc::from(record.title.as_slice())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/gapinfo.h:43-61
// ```c
// typedef enum EGapAlignOpType {
//    eGapAlignDel = 0, /**< Deletion: a gap in query */
//    eGapAlignSub = 3, /**< Substitution */
//    eGapAlignIns = 6, /**< Insertion: a gap in subject */
// } EGapAlignOpType;
// typedef struct GapEditScript {
//    EGapAlignOpType* op_type;
//    Int4* num;
//    Int4 size;
// } GapEditScript;
// ```
fn format_gap_edit_ops_for_trace(ops: &[GapEditOp]) -> String {
    let mut out = String::from("[");
    for (idx, op) in ops.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        match *op {
            GapEditOp::Sub(n) => {
                out.push('S');
                out.push_str(&n.to_string());
            }
            GapEditOp::Del(n) => {
                out.push('D');
                out.push_str(&n.to_string());
            }
            GapEditOp::Ins(n) => {
                out.push('I');
                out.push_str(&n.to_string());
            }
        }
    }
    out.push(']');
    out
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1887-1890
// ```c
// /* Get effective search space from the query information block */
// hsp->evalue =
//     BLAST_KarlinStoE_simple(score, kbp[kbp_context],
//                      query_info->contexts[hsp->context].eff_searchsp);
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1923-1926
// ```c
// hsp->bit_score =
//    (hsp->score*kbp[hsp->context]->Lambda - kbp[hsp->context]->logK) /
//    NCBIMATH_LN2;
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4157-4171
// ```c
// BLAST_KarlinStoE_simple(Int4 S, Blast_KarlinBlk* kbp, Int8 searchsp)
// {
//    ...
//    return (double) searchsp * exp((double)(-Lambda * S) + kbp->logK);
// }
// ```
fn calculate_blastn_context_statistics(
    raw_score: i32,
    params: &crate::stats::KarlinParams,
    eff_searchsp: i64,
    round_down_evalue_score: bool,
) -> (f64, f64) {
    let log_k = params.k.ln();
    let bit_score = (params.lambda * (raw_score as f64) - log_k) / NCBIMATH_LN2;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1864-1870
    // ```c
    // /* Round score down to even number for E-value calculations only. */
    // score = hsp->score;
    // if (hsp_list && hsp_list->hspcnt != 0
    //         && gapped_calculation && sbp->round_down) {
    //     score &= ~1;
    // }
    // ```
    let evalue_score = if round_down_evalue_score {
        raw_score & !1
    } else {
        raw_score
    };
    let e_value = if params.lambda < 0.0 || params.k < 0.0 || params.h < 0.0 {
        -1.0
    } else {
        (eff_searchsp as f64) * (-(params.lambda) * (evalue_score as f64) + log_k).exp()
    };
    (bit_score, e_value)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:662-664
// ```c
// if (hsp->evalue > cutoff) {
//    hsp_array[index] = Blast_HSPFree(hsp_array[index]);
// } else {
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1996-1998
// ```c
// if (hsp->evalue > cutoff) {
//    hsp_array[index] = Blast_HSPFree(hsp_array[index]);
// } else {
// ```
// The HSP survives unless `evalue > cutoff`, the same comparison as C (a NaN
// e-value survives; LOSAT rejects a NaN `-evalue` before the search).
#[allow(clippy::neg_cmp_op_on_partial_ord)]
fn hsp_survives_evalue_reap(evalue: f64, cutoff: f64) -> bool {
    !(evalue > cutoff)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4155-4160
// ```c
// #define MAX_SUBJECT_OFFSET 90000
// #define MAX_TOTAL_GAPS 3000
// ```
const MAX_SUBJECT_OFFSET: i32 = 90000;
const MAX_TOTAL_GAPS: i32 = 3000;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:54
// ```c
// #define MAX_DBSEQ_LEN 5000000
// ```
const MAX_DBSEQ_LEN: usize = 5_000_000;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:192
// ```c
// #define DBSEQ_CHUNK_OVERLAP 100
// ```
const DBSEQ_CHUNK_OVERLAP: usize = 100;
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/lookup_wrap.h:119
// ```c
// #define OFFSET_ARRAY_SIZE 4096
// ```
const OFFSET_ARRAY_SIZE: usize = 4096;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:120-121
// ```c
// const Uint1 kProtSentinel = NULLB;
// const Uint1 kNuclSentinel = 0xF;
// ```
const NUCL_SENTINEL: u8 = 0x0F;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1534-1539
// ```c
// #define OVERLAP_DIAG_CLOSE 10
// ```
const OVERLAP_DIAG_CLOSE: i32 = 10;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:50-54
// ```c
// /** Number of hash buckets in BLAST_DiagHash */
// #define DIAGHASH_NUM_BUCKETS 512
// /** Default hash chain length */
// #define DIAGHASH_CHAIN_LENGTH 256
// ```
const DIAGHASH_NUM_BUCKETS: usize = 512;
const DIAGHASH_CHAIN_LENGTH: usize = 256;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:122-143
// ```c
// typedef struct SubjectSplitStruct {
//    Uint1* sequence;
//    SSeqRange  full_range;
//    SSeqRange* seq_ranges;
//    Int4 num_seq_ranges;
//    Int4 allocated;
//    SSeqRange* hard_ranges;
//    Int4 num_hard_ranges;
//    Int4 hm_index;
//    SSeqRange* soft_ranges;
//    Int4 num_soft_ranges;
//    Int4 sm_index;
//    Int4 offset;
//    Int4 next;
// } SubjectSplitStruct;
// ```
struct SubjectSplitState {
    full_right: i32,
    hard_ranges: [(i32, i32); 1],
    soft_ranges: Vec<(i32, i32)>,
    hm_index: usize,
    sm_index: usize,
    offset: i32,
    next: i32,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:264-276
// ```c
// if (backup->offset == 0 && residual == 0 && backup->next == backup->full_range.right) {
//     subject->seq_ranges = backup->soft_ranges;
//     subject->num_seq_ranges = backup->num_soft_ranges;
//     return SUBJECT_SPLIT_OK;
// }
// /* if soft masking is off */
// if (subject->mask_type != eSoftSubjMasking) {
//     s_AllocateSeqRange(subject, backup, 1);
//     subject->seq_ranges[0].left = residual;
//     subject->seq_ranges[0].right = subject->length;
//     return SUBJECT_SPLIT_OK;
// }
// ```
#[derive(Clone)]
enum SeqRanges {
    Inline { range: (i32, i32) },
    SoftRangesFull,
    SoftRangesSlice { start_idx: usize, count: usize },
}

impl SeqRanges {
    #[inline]
    fn as_slice<'a>(
        &'a self,
        soft_ranges: &'a [(i32, i32)],
        scratch: &'a mut Vec<(i32, i32)>,
        chunk_offset: i32,
        chunk_length: i32,
    ) -> &'a [(i32, i32)] {
        match self {
            SeqRanges::Inline { range } => std::slice::from_ref(range),
            SeqRanges::SoftRangesFull => soft_ranges,
            SeqRanges::SoftRangesSlice { start_idx, count } => {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:184-198
                // ```c
                // if (backup->allocated >= num_seq_ranges) return;
                // if (backup->allocated) {
                //     sfree(subject->seq_ranges);
                // }
                // backup->allocated = num_seq_ranges;
                // subject->seq_ranges = (SSeqRange *) calloc(backup->allocated,
                //                                sizeof(SSeqRange));
                // ```
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:298-307
                // ```c
                // for (i=0; i<len; i++) {
                //     subject->seq_ranges[i].left = backup->soft_ranges[i+start].left - backup->offset;
                //     subject->seq_ranges[i].right = backup->soft_ranges[i+start].right - backup->offset;
                // }
                // if (subject->seq_ranges[0].left < 0)
                //     subject->seq_ranges[0].left = 0;
                // if (subject->seq_ranges[len-1].right > subject->length)
                //     subject->seq_ranges[len-1].right = subject->length;
                // ```
                let start = *start_idx;
                let count = *count;
                scratch.clear();
                if scratch.capacity() < count {
                    scratch.reserve(count.saturating_sub(scratch.capacity()));
                }
                let end_idx = start.saturating_add(count);
                for &(left, right) in soft_ranges[start..end_idx].iter() {
                    scratch.push((left - chunk_offset, right - chunk_offset));
                }
                if let Some(first) = scratch.first_mut() {
                    if first.0 < 0 {
                        first.0 = 0;
                    }
                }
                if let Some(last) = scratch.last_mut() {
                    if last.1 > chunk_length {
                        last.1 = chunk_length;
                    }
                }
                scratch.as_slice()
            }
        }
    }
}

#[derive(Clone)]
struct SubjectChunk {
    offset: usize,
    length: usize,
    residual: usize,
    seq_ranges: SeqRanges,
    masked: bool,
    overlap: usize,
}

enum SubjectChunkStatus {
    Done,
    NoRange,
    Ok(SubjectChunk),
}

impl SubjectSplitState {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:146-184
    // ```c
    // backup->full_range.left = 0;
    // backup->full_range.right = subject->length;
    // backup->hard_ranges = &(backup->full_range);
    // backup->num_hard_ranges = 1;
    // backup->hm_index = 0;
    // backup->soft_ranges = &(backup->full_range);
    // backup->num_soft_ranges = 1;
    // backup->sm_index = 0;
    // backup->offset = backup->hard_ranges[0].left;
    // backup->next = backup->offset;
    // ```
    fn new(subject_len: usize, soft_ranges: Vec<(i32, i32)>) -> Self {
        let full_right = subject_len as i32;
        Self {
            full_right,
            // One range, as NCBI's `&(backup->full_range)`: an array, not a vector per subject.
            hard_ranges: [(0, full_right)],
            soft_ranges,
            hm_index: 0,
            sm_index: 0,
            offset: 0,
            next: 0,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:221-310
    // ```c
    // if (backup->next >= backup->full_range.right) return SUBJECT_SPLIT_DONE;
    // residual = is_nucleotide ?  backup->next % COMPRESSION_RATIO : 0;
    // backup->offset = backup->next - residual;
    // if (backup->offset + MAX_DBSEQ_LEN < backup->hard_ranges[backup->hm_index].right) {
    //     subject->length = MAX_DBSEQ_LEN;
    //     backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
    // } else {
    //     subject->length = backup->hard_ranges[backup->hm_index].right - backup->offset;
    //     backup->hm_index++;
    //     backup->next = (backup->hm_index < backup->num_hard_ranges) ?
    //                    backup->hard_ranges[backup->hm_index].left :
    //                    backup->full_range.right;
    // }
    // if (backup->offset == 0 && residual == 0 && backup->next == backup->full_range.right) {
    //     subject->seq_ranges = backup->soft_ranges;
    //     subject->num_seq_ranges = backup->num_soft_ranges;
    //     return SUBJECT_SPLIT_OK;
    // }
    // if (subject->mask_type != eSoftSubjMasking) {
    //     subject->seq_ranges[0].left = residual;
    //     subject->seq_ranges[0].right = subject->length;
    //     return SUBJECT_SPLIT_OK;
    // }
    // /* soft masking is on, sequence is chunked, must re-allocate and adjust */
    // ASSERT(residual == 0);
    // ...
    // if (len == 0) return SUBJECT_SPLIT_NO_RANGE;
    // ```
    fn next_chunk(&mut self, subject_masked: bool, chunk_overlap: usize) -> SubjectChunkStatus {
        if self.next >= self.full_right {
            return SubjectChunkStatus::Done;
        }

        let residual = (self.next as usize) % COMPRESSION_RATIO;
        self.offset = self.next - residual as i32;
        let offset = self.offset;

        let hard_right = self
            .hard_ranges
            .get(self.hm_index)
            .map(|r| r.1)
            .unwrap_or(self.full_right);

        let length = if offset + (MAX_DBSEQ_LEN as i32) < hard_right {
            self.next = offset + MAX_DBSEQ_LEN as i32 - chunk_overlap as i32;
            MAX_DBSEQ_LEN
        } else {
            let len = (hard_right - offset).max(0) as usize;
            self.hm_index = self.hm_index.saturating_add(1);
            self.next = if self.hm_index < self.hard_ranges.len() {
                self.hard_ranges[self.hm_index].0
            } else {
                self.full_right
            };
            len
        };

        let no_chunking = offset == 0 && residual == 0 && self.next == self.full_right;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:264-276
        // ```c
        // if (backup->offset == 0 && residual == 0 && backup->next == backup->full_range.right) {
        //     subject->seq_ranges = backup->soft_ranges;
        //     subject->num_seq_ranges = backup->num_soft_ranges;
        //     return SUBJECT_SPLIT_OK;
        // }
        // /* if soft masking is off */
        // if (subject->mask_type != eSoftSubjMasking) {
        //     s_AllocateSeqRange(subject, backup, 1);
        //     subject->seq_ranges[0].left = residual;
        //     subject->seq_ranges[0].right = subject->length;
        //     return SUBJECT_SPLIT_OK;
        // }
        // ```
        let seq_ranges = if no_chunking {
            if subject_masked {
                SeqRanges::SoftRangesFull
            } else {
                SeqRanges::Inline {
                    range: (0, length as i32),
                }
            }
        } else if !subject_masked {
            SeqRanges::Inline {
                range: (residual as i32, length as i32),
            }
        } else {
            debug_assert_eq!(residual, 0);
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:279-291
            // ```c
            // start = backup->offset;
            // len = start + subject->length;
            // i = backup->sm_index;
            // while (backup->soft_ranges[i].right < start) ++i;
            // start = i;
            // while (i < backup->num_soft_ranges
            //     && backup->soft_ranges[i].left < len) ++i;
            // len = i - start;
            // backup->sm_index = i - 1;
            // ```
            let start = offset;
            let end = offset + length as i32;
            let mut i = self.sm_index;
            while i < self.soft_ranges.len() && self.soft_ranges[i].1 < start {
                i += 1;
            }
            let start_idx = i;
            while i < self.soft_ranges.len() && self.soft_ranges[i].0 < end {
                i += 1;
            }
            let count = i - start_idx;
            self.sm_index = i.saturating_sub(1);
            if count == 0 {
                return SubjectChunkStatus::NoRange;
            }
            SeqRanges::SoftRangesSlice { start_idx, count }
        };

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:580-584
        // ```c
        // overlap = (backup.offset == backup.hard_ranges[backup.hm_index].left) ?
        //           0 : dbseq_chunk_overlap;
        // ```
        let range_start = if self.hm_index < self.hard_ranges.len() {
            self.hard_ranges[self.hm_index].0
        } else {
            self.full_right
        };
        let overlap = if offset == range_start {
            0
        } else {
            chunk_overlap
        };

        SubjectChunkStatus::Ok(SubjectChunk {
            offset: offset as usize,
            length,
            residual,
            seq_ranges,
            masked: subject_masked,
            overlap,
        })
    }
}
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4163-4191
// ```c
// if (subject_length < MAX_SUBJECT_OFFSET) {
//    *start_shift = 0;
//    return;
// }
// ...
// if (s_offset <= max_extension_left) {
//    *start_shift = 0;
// } else {
//    *start_shift = s_offset - max_extension_left;
//    *subject_offset_ptr = max_extension_left;
// }
// *subject_length_ptr =
//    MIN(subject_length, s_offset + max_extension_right) - *start_shift;
// ```
#[inline]
fn adjust_subject_range(
    subject_offset: &mut i32,
    subject_length: &mut i32,
    query_offset: i32,
    query_length: i32,
) -> i32 {
    let mut start_shift: i32 = 0;
    if *subject_length < MAX_SUBJECT_OFFSET {
        return start_shift;
    }

    let s_offset = *subject_offset;
    let max_extension_left = query_offset + MAX_TOTAL_GAPS;
    let max_extension_right = query_length - query_offset + MAX_TOTAL_GAPS;

    if s_offset <= max_extension_left {
        start_shift = 0;
    } else {
        start_shift = s_offset - max_extension_left;
        *subject_offset = max_extension_left;
    }

    let adjusted_end = std::cmp::min(*subject_length, s_offset + max_extension_right);
    *subject_length = adjusted_end - start_shift;
    start_shift
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:159-182
// ```c
// if (ewp->diag_table->offset >= INT4_MAX / 4) {
//     ewp->diag_table->offset = ewp->diag_table->window;
//     s_BlastDiagClear(ewp->diag_table);
// } else {
//     ewp->diag_table->offset += subject_length + ewp->diag_table->window;
// }
// if (ewp->hash_table->offset >= INT4_MAX / 4) {
//     ewp->hash_table->occupancy = 1;
//     ewp->hash_table->offset = ewp->hash_table->window;
//     memset(ewp->hash_table->backbone, 0,
//            ewp->hash_table->num_buckets * sizeof(Int4));
// } else {
//     ewp->hash_table->offset += subject_length + ewp->hash_table->window;
// }
// ```
#[inline]
fn advance_diag_table_offset(
    diag_table_offset: &mut i32,
    diag_window: i32,
    subject_len: usize,
    use_array_indexing: bool,
    hit_level_array: &mut Vec<DiagStruct>,
    hit_len_array: &mut Vec<u8>,
    diag_hash: &mut DiagHashTable,
) {
    if use_array_indexing {
        if *diag_table_offset >= i32::MAX / 4 {
            *diag_table_offset = diag_window;
            let init_diag = DiagStruct {
                last_hit: -diag_window,
                flag: 0,
            };
            hit_level_array.fill(init_diag);
            if !hit_len_array.is_empty() {
                hit_len_array.fill(0);
            }
        } else {
            *diag_table_offset = diag_table_offset.saturating_add(subject_len as i32 + diag_window);
        }
    } else if diag_hash.offset >= i32::MAX / 4 {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:174-179
        // ```c
        // ewp->hash_table->occupancy = 1;
        // ewp->hash_table->offset = ewp->hash_table->window;
        // memset(ewp->hash_table->backbone, 0,
        //        ewp->hash_table->num_buckets * sizeof(Int4));
        // ```
        diag_hash.occupancy = 1;
        diag_hash.offset = diag_hash.window;
        diag_hash.backbone.fill(0);
        // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): every chain is empty now. No NCBI
        // counterpart (empty while the switch is off).
        diag_hash.max_level.fill(i32::MIN);
        // EXPERIMENT (LOSAT_X_DIAGLV2, `x_diag_lv.rs`): no bucket has an insert now. No NCBI
        // counterpart (empty while the switch is off).
        diag_hash.x_last.fill(x_diag_lv::LvLast::EMPTY);
        #[cfg(debug_assertions)]
        {
            diag_hash.x_last_fresh = true;
        }
    } else {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:180-181
        // ```c
        // ewp->hash_table->offset += subject_length + ewp->hash_table->window;
        // ```
        diag_hash.offset = diag_hash
            .offset
            .saturating_add(subject_len as i32 + diag_hash.window);
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:62-71
// ```c
// typedef struct DiagHashCell {
//    Int4 diag;            /**< This hit's diagonal */
//    signed int level      : 31; /**< This hit's offset in the subject sequence */
//    unsigned int hit_saved : 1;  /**< Whether or not this hit has been saved */
//    Int4  hit_len;        /**< The length of last hit */
//    Uint4 next;           /**< Offset of next element in the chain */
// }  DiagHashCell;
// ```
#[derive(Clone, Copy)]
struct DiagHashCell {
    diag: i32,
    level: i32,
    hit_saved: bool,
    hit_len: i32,
    next: u32,
}

impl Default for DiagHashCell {
    fn default() -> Self {
        Self {
            diag: 0,
            level: 0,
            hit_saved: false,
            hit_len: 0,
            next: 0,
        }
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:97-105
// ```c
// typedef struct BLAST_DiagHash {
//    Uint4 num_buckets;   /**< Number of buckets to be used for storing hit offsets */
//    Uint4 occupancy;     /**< Number of occupied elements */
//    Uint4 capacity;      /**< Total number of elements */
//    Uint4 *backbone;     /**< Array of offsets to heads of chains. */
//    DiagHashCell *chain; /**< Array of data cells. */
//    Int4 offset;         /**< "offset" added to query and subject position so that "last_hit" doesn't have to be zeroed out every time. */
//    Int4 window;         /**< The "window" size, within which two (or more) hits must be found in order to be extended. */
// } BLAST_DiagHash;
// ```
struct DiagHashTable {
    num_buckets: usize,
    occupancy: u32,
    capacity: u32,
    backbone: Vec<u32>,
    chain: Vec<DiagHashCell>,
    offset: i32,
    window: i32,
    // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): per bucket, an upper bound of the levels of
    // its cells (the largest level `insert` wrote since the table was made or cleared;
    // `i32::MIN` for none). Empty while the switch is off. No NCBI counterpart.
    max_level: Vec<i32>,
    // EXPERIMENT (LOSAT_X_DIAGLV2, `x_diag_lv.rs`): per bucket, the latest insert and an upper
    // bound of the levels of the bucket's other cells. Empty unless the switch and
    // LOSAT_X_DIAGLV are on and the table takes two hits. No NCBI counterpart. It is the memo of
    // the cells only from a dc-megablast strips stage's `x_lv2_reset` to the end of that stage
    // (and on a new or cleared table): the inserts of the reference fold (`insert_plain`) do not
    // keep it, so only that stage's folds may read it (`retrieve_below::<true>`) or keep it
    // (`insert_memo`, `insert_inl::<true>`).
    x_last: Vec<x_diag_lv::LvLast>,
    // Debug builds only (LOSAT_X_DIAGLV2): whether `x_last` is the memo of the cells (set by
    // `new`, the clear and `x_lv2_reset`, ended by an insert that does not keep a memo the table
    // has); the memo's reads and upkeep assert it.
    #[cfg(debug_assertions)]
    x_last_fresh: bool,
}

impl DiagHashTable {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:122-134
    // ```c
    // ewp->hash_table->num_buckets = DIAGHASH_NUM_BUCKETS;
    // ewp->hash_table->backbone =
    //     calloc(ewp->hash_table->num_buckets, sizeof(Uint4));
    // ewp->hash_table->capacity = DIAGHASH_CHAIN_LENGTH;
    // ewp->hash_table->chain =
    //     calloc(ewp->hash_table->capacity, sizeof(DiagHashCell));
    // ewp->hash_table->occupancy = 1;
    // ewp->hash_table->window = word_params->options->window_size;
    // ewp->hash_table->offset = word_params->options->window_size;
    // ```
    fn new(window: i32) -> Self {
        let capacity = DIAGHASH_CHAIN_LENGTH as u32;
        // EXPERIMENT (LOSAT_X_DIAGLV): the level bounds when the switch is on.
        let max_level = x_diag_lv::new_levels();
        Self {
            num_buckets: DIAGHASH_NUM_BUCKETS,
            occupancy: 1,
            capacity,
            backbone: vec![0; DIAGHASH_NUM_BUCKETS],
            chain: vec![DiagHashCell::default(); capacity as usize],
            offset: window,
            window,
            max_level,
            // EXPERIMENT (LOSAT_X_DIAGLV2): the memo when the switch is on.
            x_last: x_diag_lv::new_last(window),
            #[cfg(debug_assertions)]
            x_last_fresh: true,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:361-381
    // ```c
    // Uint4 bucket = ((Uint4) diag * 0x9E370001) % DIAGHASH_NUM_BUCKETS;
    // Uint4 index = table->backbone[bucket];
    // while (index) {
    //     if (table->chain[index].diag == diag) {
    //         *level = table->chain[index].level;
    //         *hit_len = table->chain[index].hit_len;
    //         *hit_saved = table->chain[index].hit_saved;
    //         return 1;
    //     }
    //     index = table->chain[index].next;
    // }
    // return 0;
    // ```
    fn retrieve(&self, diag: i32) -> Option<(i32, i32, bool)> {
        let bucket =
            ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        let mut index = self.backbone[bucket];
        while index != 0 {
            let cell = &self.chain[index as usize];
            if cell.diag == diag {
                return Some((cell.level, cell.hit_len, cell.hit_saved));
            }
            index = cell.next;
        }
        None
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:396-449
    // ```c
    // Uint4 bucket = ((Uint4) diag * 0x9E370001) % DIAGHASH_NUM_BUCKETS;
    // Uint4 index = table->backbone[bucket];
    // while (index) {
    //     if (table->chain[index].diag == diag) {
    //         table->chain[index].level = level;
    //         table->chain[index].hit_len = len;
    //         table->chain[index].hit_saved = hit_saved;
    //         return 1;
    //     } else {
    //         if (s_off - table->chain[index].level > window_size) {
    //             table->chain[index].diag = diag;
    //             table->chain[index].level = level;
    //             table->chain[index].hit_len = len;
    //             table->chain[index].hit_saved = hit_saved;
    //             return 1;
    //         }
    //     }
    //     index = table->chain[index].next;
    // }
    // if (table->occupancy == table->capacity) {
    //     table->capacity *= 2;
    //     table->chain =
    //         realloc(table->chain, table->capacity * sizeof(DiagHashCell));
    //     if (table->chain == NULL)
    //         return 0;
    // }
    // cell = table->chain + table->occupancy;
    // cell->diag = diag;
    // cell->level = level;
    // cell->hit_len = len;
    // cell->hit_saved = hit_saved;
    // cell->next = table->backbone[bucket];
    // table->backbone[bucket] = table->occupancy;
    // table->occupancy++;
    // return 1;
    // ```
    // The insert of a dc-megablast strips stage with the memo of LOSAT_X_DIAGLV2 (`x_last`,
    // valid from the stage's `x_lv2_reset`): the statements of `insert_plain` and the memo's
    // upkeep (`insert_inl::<true>`).
    fn insert_memo(
        &mut self,
        diag: i32,
        level: i32,
        hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        window_size: i32,
    ) -> bool {
        self.insert_inl::<true>(diag, level, hit_len, hit_saved, s_off, window_size)
    }

    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:396-449
    /// (`s_BlastDiagHashInsert`, quoted above `insert_memo`).
    /// `insert_memo` without the memo of LOSAT_X_DIAGLV2 (the code before that switch): the reference
    /// fold (`process_offset_pairs`), which never reads the memo (it lives in the strips stages
    /// of dc-megablast only, rebuilt at each stage's start: `x_diag_lv.rs`), the one-hit strips
    /// folds (`x_pair_seed.rs` `fold`, `x_pair_seed_na.rs` `fold_na`), and the strips fold of
    /// dc-megablast in its instance with LOSAT_X_FOLDLEAN off on a table without the memo.
    #[inline(never)]
    pub(super) fn insert_plain(
        &mut self,
        diag: i32,
        level: i32,
        hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        window_size: i32,
    ) -> bool {
        self.insert_inl::<false>(diag, level, hit_len, hit_saved, s_off, window_size)
    }

    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:396-449
    /// (`s_BlastDiagHashInsert`, quoted above `insert_memo`).
    /// The body of `insert_memo` (`MEMO`) and `insert_plain`, inlined into its caller; the strips fold
    /// of dc-megablast inlines it with LOSAT_X_FOLDLEAN (`x_pair_seed_disc.rs` `fold_two_hit`).
    /// The same statements. Without `MEMO` the memo (LOSAT_X_DIAGLV2) is not kept.
    #[inline(always)]
    pub(super) fn insert_inl<const MEMO: bool>(
        &mut self,
        diag: i32,
        level: i32,
        hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        window_size: i32,
    ) -> bool {
        let bucket =
            ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        // (Debug builds: the memo's upkeep needs the memo of the cells; an insert that does not
        // keep a memo the table has ends it until the next stage's `x_lv2_reset`.)
        #[cfg(debug_assertions)]
        {
            if MEMO {
                assert!(
                    self.x_last.is_empty() || self.x_last_fresh,
                    "LOSAT_X_DIAGLV2: an insert keeps a stale memo (no x_lv2_reset since an insert without it)"
                );
            } else if !self.x_last.is_empty() {
                self.x_last_fresh = false;
            }
        }
        // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): this insert writes `level` into one cell
        // of this bucket (the match, a take-over or the append below). No NCBI counterpart;
        // with the switch off `max_level` is empty and this is one length test.
        if let Some(max_level) = self.max_level.get_mut(bucket) {
            if *max_level < level {
                *max_level = level;
            }
            // EXPERIMENT (LOSAT_X_DIAGLV2, `x_diag_lv.rs`): this insert is the bucket's latest.
            // No NCBI counterpart; with the switch off `x_last` is empty and this is one length
            // test (none without `MEMO`: the reference fold's `insert_plain`).
            if MEMO {
                if let Some(last) = self.x_last.get_mut(bucket) {
                    last.note(diag, level, hit_len, hit_saved);
                }
            }
        }
        let mut index = self.backbone[bucket];
        while index != 0 {
            let cell = &mut self.chain[index as usize];
            if cell.diag == diag {
                cell.level = level;
                cell.hit_len = hit_len;
                cell.hit_saved = hit_saved;
                return true;
            }
            if s_off - cell.level > window_size {
                cell.diag = diag;
                cell.level = level;
                cell.hit_len = hit_len;
                cell.hit_saved = hit_saved;
                return true;
            }
            index = cell.next;
        }

        if self.occupancy == self.capacity {
            self.capacity = self.capacity.saturating_mul(2);
            self.chain
                .resize(self.capacity as usize, DiagHashCell::default());
        }

        let cell_index = self.occupancy;
        let cell = &mut self.chain[cell_index as usize];
        cell.diag = diag;
        cell.level = level;
        cell.hit_len = hit_len;
        cell.hit_saved = hit_saved;
        cell.next = self.backbone[bucket];
        self.backbone[bucket] = cell_index;
        self.occupancy = self.occupancy.saturating_add(1);
        true
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:833-839
// ```c
//     rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
//
//     /* if there is no record in hashtable, we set last_hit to be a very negative number */
//     if(!rc)  last_hit = 0;
//
//     /* hit within the explored area should be rejected*/
//     if (s_off_pos < last_hit) return 0;
// ```
// EXPERIMENT (LOSAT_X_DIAGLV): a fold's retrieve skips the walk of
// a bucket whose cells all lie below the level that could change a decision of the seed (a
// child module, so that it uses the private table of this file).
#[path = "x_diag_lv.rs"]
mod x_diag_lv;

#[inline(always)]
fn diag_hash_insert_window(window_size: usize, scan_range: usize, word_length: usize) -> i32 {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:825-858, 939-941
    // ```c
    // Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
    // ...
    // if (Delta < 0) Delta = 0;
    // ...
    // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
    //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
    //                       hit_ready, s_off_pos, window_size + Delta + 1);
    // ```
    // The clamp at line 858 is inside the two-hit off-diagonal block. In BLASTN
    // one-hit mode (window_size == 0) that block is skipped, so the insert sees
    // the original negative Delta.
    let window = window_size as i32;
    let delta = (scan_range as i32).min(window - word_length as i32);
    window + delta + 1
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:135-149
// ```c
// typedef union BlastOffsetPair {
//     struct { Uint4 q_off; Uint4 s_off; } qs_offsets;
// } BlastOffsetPair;
// ```
#[derive(Clone, Copy)]
struct OffsetPair {
    q_off: usize,
    s_off: usize,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:614-621
    // ```c
    // static Int4
    // s_BlastnDiagTableExtendInitialHit(...,
    //                              BlastQueryInfo * query_info,
    //                              Int4 s_range,
    //                              Int4 word_length, Int4 lut_word_length,
    //                              ...)
    // ```
    // Each seed must carry the current subject range boundary used by the
    // original NCBI extension call.
    s_range: usize,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
// ```c
// if (hsp->query.frame != hsp->subject.frame) {
//    *q_end = query_length - hsp->query.offset;
//    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
//    *s_end = hsp->subject.offset + 1;
//    *s_start = hsp->subject.end;
// } else {
//    *q_start = hsp->query.offset + 1;
//    *q_end = hsp->query.end;
//    *s_start = hsp->subject.offset + 1;
//    *s_end = hsp->subject.end;
// }
// ```
fn adjust_blastn_offsets(
    query_offset: usize,
    query_end: usize,
    subject_offset: usize,
    subject_end: usize,
    query_length: usize,
    query_frame: i32,
) -> (usize, usize, usize, usize) {
    if query_frame < 0 {
        let q_end = query_length.saturating_sub(query_offset);
        let q_start = query_length.saturating_sub(query_end).saturating_add(1);
        let s_end = subject_offset + 1;
        let s_start = subject_end;
        (q_start, q_end, s_start, s_end)
    } else {
        let q_start = query_offset + 1;
        let q_end = query_end;
        let s_start = subject_offset + 1;
        let s_end = subject_end;
        (q_start, q_end, s_start, s_end)
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_util.h:51-55
// ```c
// #define NCBI2NA_UNPACK_BASE(x, N) (((x)>>(2*(N))) & NCBI2NA_MASK)
// ```
#[inline(always)]
fn packed_base_at(packed: &[u8], pos: usize) -> u8 {
    let byte = packed[pos / COMPRESSION_RATIO];
    let shift = 2 * (3 - (pos % COMPRESSION_RATIO));
    (byte >> shift) & 0x03
}

#[inline]
fn packed_kmer_at(packed: &[u8], start: usize, k: usize) -> u64 {
    let mut kmer = 0u64;
    for i in 0..k {
        let code = packed_base_at(packed, start + i) as u64;
        kmer = (kmer << 2) | code;
    }
    kmer
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:459-487
// ```c
// Uint1 *s  = subject->sequence + s_off / COMPRESSION_RATIO;
// Int4 shift = 2* (16 - s_off % COMPRESSION_RATIO - lut_word_length);
// Int4 index;
// switch (shift) {
// case  8:
// case 10:
// case 12:
// case 14:
//     index = (s[0] << 24 | s[1] << 16 | s[2] << 8) >> shift;
// break;
// case 16:
// case 18:
// case 20:
// case 22:
//     index = (s[0] << 24 | s[1] << 16 ) >> shift;
// break;
// case 24:
//     index = s[0];
// break;
// default:
//     index = (s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]) >> shift;
// break;
// }
// ```
#[inline(always)]
fn packed_kmer_at_seed_mask(packed: &[u8], start: usize, lut_word_length: usize) -> u64 {
    if lut_word_length == 0 {
        return 0;
    }

    let shift = 2i32 * (16i32 - (start % COMPRESSION_RATIO) as i32 - lut_word_length as i32);
    if shift < 0 || shift > 24 {
        return packed_kmer_at(packed, start, lut_word_length);
    }

    let byte_idx = start / COMPRESSION_RATIO;
    let shift_u = shift as u32;

    match shift {
        8 | 10 | 12 | 14 => {
            if byte_idx + 3 <= packed.len() {
                let s0 = packed[byte_idx] as u32;
                let s1 = packed[byte_idx + 1] as u32;
                let s2 = packed[byte_idx + 2] as u32;
                (((s0 << 24) | (s1 << 16) | (s2 << 8)) >> shift_u) as u64
            } else {
                packed_kmer_at(packed, start, lut_word_length)
            }
        }
        16 | 18 | 20 | 22 => {
            if byte_idx + 2 <= packed.len() {
                let s0 = packed[byte_idx] as u32;
                let s1 = packed[byte_idx + 1] as u32;
                (((s0 << 24) | (s1 << 16)) >> shift_u) as u64
            } else {
                packed_kmer_at(packed, start, lut_word_length)
            }
        }
        24 => {
            if byte_idx < packed.len() {
                packed[byte_idx] as u64
            } else {
                packed_kmer_at(packed, start, lut_word_length)
            }
        }
        _ => {
            if byte_idx + 4 <= packed.len() {
                let s0 = packed[byte_idx] as u32;
                let s1 = packed[byte_idx + 1] as u32;
                let s2 = packed[byte_idx + 2] as u32;
                let s3 = packed[byte_idx + 3] as u32;
                (((s0 << 24) | (s1 << 16) | (s2 << 8) | s3) >> shift_u) as u64
            } else {
                packed_kmer_at(packed, start, lut_word_length)
            }
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1482-1623
// ```c
// if (scan_step % COMPRESSION_RATIO == 0 &&
//     (subject->mask_type == eNoSubjMasking)) {
//     Uint1 *s_end = abs_start + scan_range[1] / COMPRESSION_RATIO;
//     Int4 shift = 2 * (12 - lut_word_length);
//     s = abs_start + scan_range[0] / COMPRESSION_RATIO;
//     scan_step = scan_step / COMPRESSION_RATIO;
//     for ( ; s <= s_end; s += scan_step) {
//         index = s[0] << 16 | s[1] << 8 | s[2];
//         index = index >> shift;
//         ...
//     }
// } else if (lut_word_length > 9) {
//     for (; scan_range[0] <= scan_range[1]; scan_range[0] += scan_step) {
//         Int4 shift = 2*(16 - (scan_range[0] % COMPRESSION_RATIO + lut_word_length));
//         s = abs_start + (scan_range[0] / COMPRESSION_RATIO);
//         index = s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3];
//         index = (index >> shift) & mask;
//         MB_ACCESS_HITS();
//     }
// } else {
//     for (; scan_range[0] <= scan_range[1]; scan_range[0] += scan_step) {
//         Int4 shift = 2*(12 - (scan_range[0] % COMPRESSION_RATIO + lut_word_length));
//         s = abs_start + (scan_range[0] / COMPRESSION_RATIO);
//         index = s[0] << 16 | s[1] << 8 | s[2];
//         index = (index >> shift) & mask;
//         MB_ACCESS_HITS();
//     }
// }
// ```
fn scan_subject_kmers_range_mb_any<F>(
    packed: &[u8],
    subject_len: usize,
    lut_word_length: usize,
    scan_step: usize,
    subject_masked: bool,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64),
{
    if subject_len < lut_word_length || lut_word_length < 9 || lut_word_length > 12 {
        return;
    }

    let packed_len = packed.len();
    if packed_len == 0 {
        return;
    }

    let mask = (1u64 << (2 * lut_word_length)) - 1;

    if scan_step % COMPRESSION_RATIO == 0 && !subject_masked {
        if packed_len >= 3 {
            let shift = 2 * (12 - lut_word_length);
            let step_bytes = scan_step / COMPRESSION_RATIO;
            let mut byte_idx = start / COMPRESSION_RATIO;
            let max_byte_idx = packed_len.saturating_sub(3);
            let end_byte_idx = (end / COMPRESSION_RATIO).min(max_byte_idx);
            while byte_idx <= end_byte_idx {
                let idx = ((packed[byte_idx] as u32) << 16)
                    | ((packed[byte_idx + 1] as u32) << 8)
                    | (packed[byte_idx + 2] as u32);
                let kmer = ((idx >> shift) as u64) & mask;
                on_kmer(byte_idx * COMPRESSION_RATIO, kmer);
                byte_idx = byte_idx.saturating_add(step_bytes);
            }
            let mut pos = byte_idx * COMPRESSION_RATIO;
            while pos <= end {
                let kmer = packed_kmer_at(packed, pos, lut_word_length);
                on_kmer(pos, kmer);
                pos = pos.saturating_add(scan_step);
            }
            return;
        }
    }

    if lut_word_length > 9 {
        if packed_len >= 4 {
            let max_byte_idx = packed_len.saturating_sub(4);
            let max_fast_pos = max_byte_idx * COMPRESSION_RATIO + (COMPRESSION_RATIO - 1);
            let fast_end = end.min(max_fast_pos);
            let mut pos = start;
            while pos <= fast_end {
                let byte_idx = pos / COMPRESSION_RATIO;
                let shift = 2 * (16 - ((pos % COMPRESSION_RATIO) + lut_word_length));
                let idx = ((packed[byte_idx] as u32) << 24)
                    | ((packed[byte_idx + 1] as u32) << 16)
                    | ((packed[byte_idx + 2] as u32) << 8)
                    | (packed[byte_idx + 3] as u32);
                let kmer = ((idx >> shift) as u64) & mask;
                on_kmer(pos, kmer);
                pos = pos.saturating_add(scan_step);
            }
            while pos <= end {
                let kmer = packed_kmer_at(packed, pos, lut_word_length);
                on_kmer(pos, kmer);
                pos = pos.saturating_add(scan_step);
            }
            return;
        }
    } else if packed_len >= 3 {
        let max_byte_idx = packed_len.saturating_sub(3);
        let max_fast_pos = max_byte_idx * COMPRESSION_RATIO + (COMPRESSION_RATIO - 1);
        let fast_end = end.min(max_fast_pos);
        let mut pos = start;
        while pos <= fast_end {
            let byte_idx = pos / COMPRESSION_RATIO;
            let shift = 2 * (12 - ((pos % COMPRESSION_RATIO) + lut_word_length));
            let idx = ((packed[byte_idx] as u32) << 16)
                | ((packed[byte_idx + 1] as u32) << 8)
                | (packed[byte_idx + 2] as u32);
            let kmer = ((idx >> shift) as u64) & mask;
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, lut_word_length);
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mut pos = start;
    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, lut_word_length);
        on_kmer(pos, kmer);
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1639-1701
// ```c
// static Int4 s_MBScanSubject_9_1(...)
// {
//     ...
//     switch (scan_range[0] % COMPRESSION_RATIO) {
//     case 1: init_index = s[0] << 16 | s[1] << 8 | s[2]; goto base_1;
//     case 2: init_index = s[0] << 16 | s[1] << 8 | s[2]; goto base_2;
//     case 3: init_index = s[0] << 16 | s[1] << 8 | s[2]; goto base_3;
//     }
//     while (scan_range[0] <= scan_range[1]) {
//         init_index = s[0] << 16 | s[1] << 8 | s[2];
//         index = init_index >> 6;
//         MB_ACCESS_HITS();
//         scan_range[0]++;
// base_1:
//         ...
// base_2:
//         ...
// base_3:
//         index = init_index & kLutWordMask;
//         s++;
//         MB_ACCESS_HITS();
//         scan_range[0]++;
//     }
// }
// ```
#[inline(always)]
fn scan_subject_kmers_range_mb_9_1<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64),
{
    const LUT_WORD_LENGTH: usize = 9;
    debug_assert_eq!(scan_step, 1);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 3 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut state = pos % COMPRESSION_RATIO;
    let mut init_index: u32 = 0;
    let mut init_valid = false;

    while pos <= end {
        if !init_valid {
            if byte_idx + 2 >= packed_len {
                break;
            }
            init_index = ((packed[byte_idx] as u32) << 16)
                | ((packed[byte_idx + 1] as u32) << 8)
                | (packed[byte_idx + 2] as u32);
            init_valid = true;
        }

        match state {
            0 => {
                let index = ((init_index >> 6) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 1;
            }
            1 => {
                let index = ((init_index >> 4) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 2;
            }
            2 => {
                let index = ((init_index >> 2) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 3;
            }
            _ => {
                let index = (init_index as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 0;
                byte_idx = byte_idx.saturating_add(1);
                init_valid = false;
            }
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        on_kmer(pos, kmer);
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1715-1755
// ```c
// static Int4 s_MBScanSubject_9_2(...)
// {
//     ...
//     if (scan_range[0] % COMPRESSION_RATIO == 2) {
//         init_index = s[0] << 16 | s[1] << 8 | s[2];
//         goto base_2;
//     }
//     while (scan_range[0] <= scan_range[1]) {
//         init_index = s[0] << 16 | s[1] << 8 | s[2];
//         index = init_index >> 6;
//         MB_ACCESS_HITS();
//         scan_range[0] += 2;
// base_2:
//         ...
//         index = (init_index >> 2) & kLutWordMask;
//         s++;
//         MB_ACCESS_HITS();
//         scan_range[0] += 2;
//     }
// }
// ```
#[inline(always)]
fn scan_subject_kmers_range_mb_9_2<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64),
{
    const LUT_WORD_LENGTH: usize = 9;
    debug_assert_eq!(scan_step, 2);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 3 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut use_shift2 = pos % COMPRESSION_RATIO == 2;
    let mut init_index: u32 = 0;
    let mut init_valid = false;

    while pos <= end {
        if !init_valid {
            if byte_idx + 2 >= packed_len {
                break;
            }
            init_index = ((packed[byte_idx] as u32) << 16)
                | ((packed[byte_idx + 1] as u32) << 8)
                | (packed[byte_idx + 2] as u32);
            init_valid = true;
        }

        if !use_shift2 {
            let index = ((init_index >> 6) as u64) & mask;
            on_kmer(pos, index);
            pos = pos.saturating_add(scan_step);
            use_shift2 = true;
        } else {
            let index = ((init_index >> 2) as u64) & mask;
            on_kmer(pos, index);
            pos = pos.saturating_add(scan_step);
            byte_idx = byte_idx.saturating_add(1);
            init_valid = false;
            use_shift2 = false;
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        on_kmer(pos, kmer);
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1768-1832
// ```c
// static Int4 s_MBScanSubject_10_1(...)
// {
//     ...
//     switch (scan_range[0] % COMPRESSION_RATIO) {
//     case 1: init_index = s[0] << 16 | s[1] << 8 | s[2]; goto base_1;
//     case 2: init_index = s[0] << 16 | s[1] << 8 | s[2]; goto base_2;
//     case 3: init_index = s[0] << 16 | s[1] << 8 | s[2]; goto base_3;
//     }
//     while (scan_range[0] <= scan_range[1]) {
//         init_index = s[0] << 16 | s[1] << 8 | s[2];
//         index = init_index >> 4;
//         MB_ACCESS_HITS();
//         scan_range[0]++;
// base_1:
//         ...
// base_2:
//         ...
// base_3:
//         init_index = init_index << 8 | s[3];
//         index = (init_index >> 6) & kLutWordMask;
//         s++;
//         MB_ACCESS_HITS();
//         scan_range[0]++;
//     }
// }
// ```
#[inline(always)]
fn scan_subject_kmers_range_mb_10_1<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64),
{
    const LUT_WORD_LENGTH: usize = 10;
    debug_assert_eq!(scan_step, 1);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 4 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut state = pos % COMPRESSION_RATIO;
    let mut init_index: u32 = 0;
    let mut init_valid = false;

    while pos <= end {
        if !init_valid {
            if byte_idx + 2 >= packed_len {
                break;
            }
            init_index = ((packed[byte_idx] as u32) << 16)
                | ((packed[byte_idx + 1] as u32) << 8)
                | (packed[byte_idx + 2] as u32);
            init_valid = true;
        }

        match state {
            0 => {
                let index = ((init_index >> 4) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 1;
            }
            1 => {
                let index = ((init_index >> 2) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 2;
            }
            2 => {
                let index = (init_index as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 3;
            }
            _ => {
                if byte_idx + 3 >= packed_len {
                    break;
                }
                let s3 = packed[byte_idx + 3] as u32;
                init_index = (init_index << 8) | s3;
                let index = ((init_index >> 6) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 0;
                byte_idx = byte_idx.saturating_add(1);
                init_valid = false;
            }
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        on_kmer(pos, kmer);
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1845-1884
// ```c
// static Int4 s_MBScanSubject_10_2(...)
// {
//     ...
//     if (scan_range[0] % COMPRESSION_RATIO == 2) {
//         init_index = s[0] << 16 | s[1] << 8 | s[2];
//         goto base_2;
//     }
//     while (scan_range[0] <= scan_range[1]) {
//         init_index = s[0] << 16 | s[1] << 8 | s[2];
//         index = init_index >> 4;
//         MB_ACCESS_HITS();
//         scan_range[0] += 2;
// base_2:
//         ...
//         index = init_index & kLutWordMask;
//         s++;
//         MB_ACCESS_HITS();
//         scan_range[0] += 2;
//     }
// }
// ```
#[inline(always)]
fn scan_subject_kmers_range_mb_10_2<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64),
{
    const LUT_WORD_LENGTH: usize = 10;
    debug_assert_eq!(scan_step, 2);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 3 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut use_shift0 = pos % COMPRESSION_RATIO == 2;
    let mut init_index: u32 = 0;
    let mut init_valid = false;

    while pos <= end {
        if !init_valid {
            if byte_idx + 2 >= packed_len {
                break;
            }
            init_index = ((packed[byte_idx] as u32) << 16)
                | ((packed[byte_idx + 1] as u32) << 8)
                | (packed[byte_idx + 2] as u32);
            init_valid = true;
        }

        if !use_shift0 {
            let index = ((init_index >> 4) as u64) & mask;
            on_kmer(pos, index);
            pos = pos.saturating_add(scan_step);
            use_shift0 = true;
        } else {
            let index = (init_index as u64) & mask;
            on_kmer(pos, index);
            pos = pos.saturating_add(scan_step);
            byte_idx = byte_idx.saturating_add(1);
            init_valid = false;
            use_shift0 = false;
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        on_kmer(pos, kmer);
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1898-1964
// ```c
// static Int4 s_MBScanSubject_10_3(...)
// {
//     ...
//     switch (scan_range[0] % COMPRESSION_RATIO) {
//     case 1: init_index = s[0] << 8 | s[1]; s -= 2; goto base_3;
//     case 2: init_index = s[0] << 16 | s[1] << 8 | s[2]; s--; goto base_2;
//     case 3: init_index = s[0] << 16 | s[1] << 8 | s[2]; goto base_1;
//     }
//     while (scan_range[0] <= scan_range[1]) {
//         init_index = s[0] << 16 | s[1] << 8 | s[2];
//         index = init_index >> 4;
//         MB_ACCESS_HITS();
//         scan_range[0] += 3;
// base_1:
//         init_index = init_index << 8 | s[3];
//         index = (init_index >> 6) & kLutWordMask;
//         MB_ACCESS_HITS();
//         scan_range[0] += 3;
// base_2:
//         index = init_index & kLutWordMask;
//         MB_ACCESS_HITS();
//         scan_range[0] += 3;
// base_3:
//         init_index = init_index << 8 | s[4];
//         index = (init_index >> 2) & kLutWordMask;
//         s += 3;
//         MB_ACCESS_HITS();
//         scan_range[0] += 3;
//     }
// }
// ```
#[inline(always)]
fn scan_subject_kmers_range_mb_10_3<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64),
{
    const LUT_WORD_LENGTH: usize = 10;
    debug_assert_eq!(scan_step, 3);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 5 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let packed_len_i32 = packed_len as i32;
    let mut pos = start;
    let mut byte_idx = (pos / COMPRESSION_RATIO) as i32;
    let mut state = 0u8;
    let mut init_index: u32 = 0;
    let mut use_fast = true;

    match pos % COMPRESSION_RATIO {
        1 => {
            if byte_idx + 1 >= packed_len_i32 {
                use_fast = false;
            } else {
                init_index = ((packed[byte_idx as usize] as u32) << 8)
                    | (packed[(byte_idx + 1) as usize] as u32);
                byte_idx -= 2;
                state = 3;
            }
        }
        2 => {
            if byte_idx + 2 >= packed_len_i32 {
                use_fast = false;
            } else {
                init_index = ((packed[byte_idx as usize] as u32) << 16)
                    | ((packed[(byte_idx + 1) as usize] as u32) << 8)
                    | (packed[(byte_idx + 2) as usize] as u32);
                byte_idx -= 1;
                state = 2;
            }
        }
        3 => {
            if byte_idx + 2 >= packed_len_i32 {
                use_fast = false;
            } else {
                init_index = ((packed[byte_idx as usize] as u32) << 16)
                    | ((packed[(byte_idx + 1) as usize] as u32) << 8)
                    | (packed[(byte_idx + 2) as usize] as u32);
                state = 1;
            }
        }
        _ => {
            state = 0;
        }
    }

    if !use_fast {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    while pos <= end {
        match state {
            0 => {
                if byte_idx < 0 || byte_idx + 2 >= packed_len_i32 {
                    break;
                }
                let idx0 = byte_idx as usize;
                let idx1 = (byte_idx + 1) as usize;
                let idx2 = (byte_idx + 2) as usize;
                init_index = ((packed[idx0] as u32) << 16)
                    | ((packed[idx1] as u32) << 8)
                    | (packed[idx2] as u32);
                let index = ((init_index >> 4) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 1;
            }
            1 => {
                let idx3 = byte_idx + 3;
                if idx3 < 0 || idx3 >= packed_len_i32 {
                    break;
                }
                init_index = (init_index << 8) | (packed[idx3 as usize] as u32);
                let index = ((init_index >> 6) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 2;
            }
            2 => {
                let index = (init_index as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                state = 3;
            }
            _ => {
                let idx4 = byte_idx + 4;
                if idx4 < 0 || idx4 >= packed_len_i32 {
                    break;
                }
                init_index = (init_index << 8) | (packed[idx4 as usize] as u32);
                let index = ((init_index >> 2) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                byte_idx += 3;
                state = 0;
            }
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        on_kmer(pos, kmer);
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1978-2041
// ```c
// static Int4 s_MBScanSubject_11_1Mod4(...)
// {
//     ...
//     switch (scan_range[0] % COMPRESSION_RATIO) {
//     case 1: goto base_1;
//     case 2: goto base_2;
//     case 3: goto base_3;
//     }
//     while (scan_range[0] <= scan_range[1]) {
//         index = s[0] << 16 | s[1] << 8 | s[2];
//         index = index >> 2;
//         s += scan_step_byte;
//         MB_ACCESS_HITS();
//         scan_range[0] += scan_step;
// base_1:
//         ...
// base_2:
//         index = s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3];
//         index = (index >> 6) & kLutWordMask;
//         s += scan_step_byte;
//         MB_ACCESS_HITS();
//         scan_range[0] += scan_step;
// base_3:
//         index = s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3];
//         index = (index >> 4) & kLutWordMask;
//         s += scan_step_byte + 1;
//         MB_ACCESS_HITS();
//         scan_range[0] += scan_step;
//     }
// }
// ```
#[inline(always)]
fn scan_subject_kmers_range_mb_11_1mod4<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64),
{
    const LUT_WORD_LENGTH: usize = 11;
    debug_assert_eq!(scan_step % COMPRESSION_RATIO, 1);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 4 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let scan_step_byte = scan_step / COMPRESSION_RATIO;
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut state = pos % COMPRESSION_RATIO;

    while pos <= end {
        match state {
            0 => {
                if byte_idx + 2 >= packed_len {
                    break;
                }
                let idx = ((packed[byte_idx] as u32) << 16)
                    | ((packed[byte_idx + 1] as u32) << 8)
                    | (packed[byte_idx + 2] as u32);
                let index = ((idx >> 2) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                byte_idx = byte_idx.saturating_add(scan_step_byte);
                state = 1;
            }
            1 => {
                if byte_idx + 2 >= packed_len {
                    break;
                }
                let idx = ((packed[byte_idx] as u32) << 16)
                    | ((packed[byte_idx + 1] as u32) << 8)
                    | (packed[byte_idx + 2] as u32);
                let index = (idx as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                byte_idx = byte_idx.saturating_add(scan_step_byte);
                state = 2;
            }
            2 => {
                if byte_idx + 3 >= packed_len {
                    break;
                }
                let idx = ((packed[byte_idx] as u32) << 24)
                    | ((packed[byte_idx + 1] as u32) << 16)
                    | ((packed[byte_idx + 2] as u32) << 8)
                    | (packed[byte_idx + 3] as u32);
                let index = ((idx >> 6) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                byte_idx = byte_idx.saturating_add(scan_step_byte);
                state = 3;
            }
            _ => {
                if byte_idx + 3 >= packed_len {
                    break;
                }
                let idx = ((packed[byte_idx] as u32) << 24)
                    | ((packed[byte_idx + 1] as u32) << 16)
                    | ((packed[byte_idx + 2] as u32) << 8)
                    | (packed[byte_idx + 3] as u32);
                let index = ((idx >> 4) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                byte_idx = byte_idx.saturating_add(scan_step_byte + 1);
                state = 0;
            }
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        on_kmer(pos, kmer);
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2055-2107
// ```c
// static Int4 s_MBScanSubject_11_2Mod4(...)
// {
//     ...
//     if ((scan_range[0] % 2) == 0) { top_shift = 2; bottom_shift = 6; }
//     else { top_shift = 0; bottom_shift = 4; }
//     if ((scan_range[0] % COMPRESSION_RATIO == 2) ||
//         (scan_range[0] % COMPRESSION_RATIO == 3))
//         goto base_23;
//     while (scan_range[0] <= scan_range[1]) {
//         index = s[0] << 16 | s[1] << 8 | s[2];
//         index = (index >> top_shift) & kLutWordMask;
//         s += scan_step_byte;
//         MB_ACCESS_HITS();
//         scan_range[0] += scan_step;
// base_23:
//         index = s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3];
//         index = (index >> bottom_shift) & kLutWordMask;
//         s += scan_step_byte + 1;
//         MB_ACCESS_HITS();
//         scan_range[0] += scan_step;
//     }
// }
// ```
#[inline(always)]
fn scan_subject_kmers_range_mb_11_2mod4<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64),
{
    const LUT_WORD_LENGTH: usize = 11;
    debug_assert_eq!(scan_step % COMPRESSION_RATIO, 2);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 4 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let scan_step_byte = scan_step / COMPRESSION_RATIO;
    let (top_shift, bottom_shift) = if (start % 2) == 0 {
        (2u32, 6u32)
    } else {
        (0u32, 4u32)
    };
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut use_bottom = pos % COMPRESSION_RATIO >= 2;

    while pos <= end {
        if !use_bottom {
            if byte_idx + 2 >= packed_len {
                break;
            }
            let idx = ((packed[byte_idx] as u32) << 16)
                | ((packed[byte_idx + 1] as u32) << 8)
                | (packed[byte_idx + 2] as u32);
            let index = ((idx >> top_shift) as u64) & mask;
            on_kmer(pos, index);
            pos = pos.saturating_add(scan_step);
            byte_idx = byte_idx.saturating_add(scan_step_byte);
            use_bottom = true;
        } else {
            if byte_idx + 3 >= packed_len {
                break;
            }
            let idx = ((packed[byte_idx] as u32) << 24)
                | ((packed[byte_idx + 1] as u32) << 16)
                | ((packed[byte_idx + 2] as u32) << 8)
                | (packed[byte_idx + 3] as u32);
            let index = ((idx >> bottom_shift) as u64) & mask;
            on_kmer(pos, index);
            pos = pos.saturating_add(scan_step);
            byte_idx = byte_idx.saturating_add(scan_step_byte + 1);
            use_bottom = false;
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        on_kmer(pos, kmer);
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2120-2188
// ```c
// static Int4 s_MBScanSubject_11_3Mod4(...)
// {
//     ...
//     switch (scan_range[0] % COMPRESSION_RATIO) {
//     case 1: s -= 2; goto base_3;
//     case 2: s--; goto base_2;
//     case 3: goto base_1;
//     }
//     while (scan_range[0] <= scan_range[1]) {
//         index = s[0] << 16 | s[1] << 8 | s[2];
//         index = index >> 2;
//         s += scan_step_byte;
//         MB_ACCESS_HITS();
//         scan_range[0] += scan_step;
// base_1:
//         index = s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3];
//         index = (index >> 4) & kLutWordMask;
//         s += scan_step_byte;
//         MB_ACCESS_HITS();
//         scan_range[0] += scan_step;
// base_2:
//         index = s[1] << 24 | s[2] << 16 | s[3] << 8 | s[4];
//         index = (index >> 6) & kLutWordMask;
//         s += scan_step_byte;
//         MB_ACCESS_HITS();
//         scan_range[0] += scan_step;
// base_3:
//         index = s[2] << 16 | s[3] << 8 | s[4];
//         index = index & kLutWordMask;
//         s += scan_step_byte + 3;
//         MB_ACCESS_HITS();
//         scan_range[0] += scan_step;
//     }
// }
// ```
#[inline(always)]
fn scan_subject_kmers_range_mb_11_3mod4<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64),
{
    const LUT_WORD_LENGTH: usize = 11;
    debug_assert_eq!(scan_step % COMPRESSION_RATIO, 3);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 5 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            on_kmer(pos, kmer);
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let packed_len_i32 = packed_len as i32;
    let scan_step_byte = scan_step / COMPRESSION_RATIO;
    let mut pos = start;
    let mut byte_idx = (pos / COMPRESSION_RATIO) as i32;
    let mut state = 0u8;

    match pos % COMPRESSION_RATIO {
        1 => {
            byte_idx -= 2;
            state = 3;
        }
        2 => {
            byte_idx -= 1;
            state = 2;
        }
        3 => {
            state = 1;
        }
        _ => {
            state = 0;
        }
    }

    while pos <= end {
        match state {
            0 => {
                if byte_idx < 0 || byte_idx + 2 >= packed_len_i32 {
                    break;
                }
                let idx0 = byte_idx as usize;
                let idx1 = (byte_idx + 1) as usize;
                let idx2 = (byte_idx + 2) as usize;
                let idx = ((packed[idx0] as u32) << 16)
                    | ((packed[idx1] as u32) << 8)
                    | (packed[idx2] as u32);
                let index = ((idx >> 2) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                byte_idx += scan_step_byte as i32;
                state = 1;
            }
            1 => {
                let idx3 = byte_idx + 3;
                if idx3 < 0 || idx3 >= packed_len_i32 {
                    break;
                }
                let idx0 = byte_idx as usize;
                let idx1 = (byte_idx + 1) as usize;
                let idx2 = (byte_idx + 2) as usize;
                let idx3u = idx3 as usize;
                let idx = ((packed[idx0] as u32) << 24)
                    | ((packed[idx1] as u32) << 16)
                    | ((packed[idx2] as u32) << 8)
                    | (packed[idx3u] as u32);
                let index = ((idx >> 4) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                byte_idx += scan_step_byte as i32;
                state = 2;
            }
            2 => {
                let idx4 = byte_idx + 4;
                if idx4 < 0 || idx4 >= packed_len_i32 {
                    break;
                }
                let idx1 = (byte_idx + 1) as usize;
                let idx2 = (byte_idx + 2) as usize;
                let idx3 = (byte_idx + 3) as usize;
                let idx4u = idx4 as usize;
                let idx = ((packed[idx1] as u32) << 24)
                    | ((packed[idx2] as u32) << 16)
                    | ((packed[idx3] as u32) << 8)
                    | (packed[idx4u] as u32);
                let index = ((idx >> 6) as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                byte_idx += scan_step_byte as i32;
                state = 3;
            }
            _ => {
                let idx4 = byte_idx + 4;
                if idx4 < 0 || idx4 >= packed_len_i32 {
                    break;
                }
                let idx2 = (byte_idx + 2) as usize;
                let idx3 = (byte_idx + 3) as usize;
                let idx4u = idx4 as usize;
                let idx = ((packed[idx2] as u32) << 16)
                    | ((packed[idx3] as u32) << 8)
                    | (packed[idx4u] as u32);
                let index = (idx as u64) & mask;
                on_kmer(pos, index);
                pos = pos.saturating_add(scan_step);
                byte_idx += (scan_step_byte + 3) as i32;
                state = 0;
            }
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        on_kmer(pos, kmer);
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1635-1676
// ```c
// scansub = (TNaScanSubjectFunction)lookup->scansub_callback;
// ...
// while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
//     hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
//     ...
// }
// ```
#[derive(Clone, Copy)]
enum MbScanSubjectKind {
    Any,
    Scan9_1,
    Scan9_2,
    Scan10_1,
    Scan10_2,
    Scan10_3,
    Scan11_1Mod4,
    Scan11_2Mod4,
    Scan11_3Mod4,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2602-2677
// ```c
// static void s_MBChooseScanSubject(LookupTableWrap *lookup_wrap)
// {
//     ...
//     switch (mb_lt->lut_word_length) {
//     case 9:  ...
//     case 10: ...
//     case 11: ...
//     case 12:
//     case 16:
//         mb_lt->scansub_callback = (void *)s_MBScanSubject_Any;
//         break;
//     }
// }
// ```
#[inline(always)]
fn choose_mb_scan_subject_kind(lut_word_length: usize, scan_step: usize) -> MbScanSubjectKind {
    match lut_word_length {
        9 => {
            if scan_step == 1 {
                MbScanSubjectKind::Scan9_1
            } else if scan_step == 2 {
                MbScanSubjectKind::Scan9_2
            } else {
                MbScanSubjectKind::Any
            }
        }
        10 => {
            if scan_step == 1 {
                MbScanSubjectKind::Scan10_1
            } else if scan_step == 2 {
                MbScanSubjectKind::Scan10_2
            } else if scan_step == 3 {
                MbScanSubjectKind::Scan10_3
            } else {
                MbScanSubjectKind::Any
            }
        }
        11 => match scan_step % COMPRESSION_RATIO {
            0 => MbScanSubjectKind::Any,
            1 => MbScanSubjectKind::Scan11_1Mod4,
            2 => MbScanSubjectKind::Scan11_2Mod4,
            _ => MbScanSubjectKind::Scan11_3Mod4,
        },
        12 => MbScanSubjectKind::Any,
        _ => MbScanSubjectKind::Any,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2602-2677
// ```c
// static void s_MBChooseScanSubject(LookupTableWrap *lookup_wrap)
// {
//     ...
//     switch (mb_lt->lut_word_length) {
//     case 9:
//         if (scan_step == 1)
//             mb_lt->scansub_callback = (void *)s_MBScanSubject_9_1;
//         if (scan_step == 2)
//             mb_lt->scansub_callback = (void *)s_MBScanSubject_9_2;
//         else
//             mb_lt->scansub_callback = (void *)s_MBScanSubject_Any;
//         break;
//     ...
//     }
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1651-1667
// ```c
// if (subject->mask_type != eNoSubjMasking) {
//     ...
//     scansub = (TNaScanSubjectFunction)
//           BlastChooseNucleotideScanSubjectAny(lookup_wrap);
//     ...
//     scan_range[1] = subject->seq_ranges[0].left + word_length - lut_word_length;
//     scan_range[2] = subject->seq_ranges[0].right - lut_word_length;
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2994-3006
// ```c
// void * BlastChooseNucleotideScanSubjectAny(LookupTableWrap *lookup_wrap)
// {
//     ...
//     return (void *)s_MBScanSubject_Any;
// }
// ```
#[inline(always)]
fn select_mb_scan_kind(
    lut_word_length: usize,
    scan_step: usize,
    subject_masked: bool,
) -> Option<MbScanSubjectKind> {
    let base = if (9..=12).contains(&lut_word_length) {
        Some(choose_mb_scan_subject_kind(lut_word_length, scan_step))
    } else {
        None
    };
    if subject_masked {
        base.map(|_| MbScanSubjectKind::Any)
    } else {
        base
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:56-66
// ```c
// index &= (mb_lt->hashsize-1);
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:88-91
// ```c
// index = lookup->final_backbone[index & lookup->mask];
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:119-123
// ```c
// index &= (lookup->mask);
// ```
#[inline(always)]
fn mask_lookup_index(index: u64, lut_word_length: usize) -> u64 {
    let mask = (1u64 << (2 * lut_word_length)) - 1;
    index & mask
}

fn scan_subject_kmers_range<F>(
    packed: &[u8],
    subject_len: usize,
    lut_word_length: usize,
    scan_step: usize,
    subject_masked: bool,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64),
{
    if subject_len < lut_word_length || lut_word_length == 0 {
        return;
    }
    if start > end || start + lut_word_length > subject_len {
        return;
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:165-176
    // ```c
    // ASSERT(lookup->scan_step > 0);
    // ```
    debug_assert!(scan_step > 0);

    if lut_word_length <= 8 {
        let mask = (1u64 << (2 * lut_word_length)) - 1;
        let packed_len = packed.len();
        if packed_len == 0 {
            return;
        }
        if lut_word_length > 5 {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:172-209
            // ```c
            // if (lut_word_length > 5) {
            //     if (scan_step % COMPRESSION_RATIO == 0 &&
            //         (subject->mask_type == eNoSubjMasking)) {
            //         Uint1 *s_end = abs_start + scan_range[1] / COMPRESSION_RATIO;
            //         Int4 shift = 2 * (8 - lut_word_length);
            //         s = abs_start + scan_range[0] / COMPRESSION_RATIO;
            //         scan_step = scan_step / COMPRESSION_RATIO;
            //         for (; s <= s_end; s += scan_step) {
            //             index = s[0] << 8 | s[1];
            //             index = index >> shift;
            //             ...
            //         }
            //     }
            // }
            // ```
            if scan_step % COMPRESSION_RATIO == 0
                && !subject_masked
                && start % COMPRESSION_RATIO == 0
                && packed_len >= 2
            {
                let shift = 2 * (8 - lut_word_length);
                let step_bytes = scan_step / COMPRESSION_RATIO;
                let mut byte_idx = start / COMPRESSION_RATIO;
                let max_byte_idx = packed_len.saturating_sub(2);
                let end_byte_idx = (end / COMPRESSION_RATIO).min(max_byte_idx);
                while byte_idx <= end_byte_idx {
                    let idx = ((packed[byte_idx] as u16) << 8) | packed[byte_idx + 1] as u16;
                    let kmer = ((idx >> shift) as u64) & mask;
                    on_kmer(byte_idx * COMPRESSION_RATIO, kmer);
                    byte_idx += step_bytes;
                }
                let mut pos = byte_idx * COMPRESSION_RATIO;
                while pos <= end {
                    let kmer = packed_kmer_at(packed, pos, lut_word_length);
                    on_kmer(pos, kmer);
                    pos = pos.saturating_add(scan_step);
                }
                return;
            }

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:210-240
            // ```c
            // for (; scan_range[0] <= scan_range[1]; scan_range[0] += scan_step) {
            //     Int4 shift = 2*(12 - (scan_range[0] % COMPRESSION_RATIO + lut_word_length));
            //     s = abs_start + (scan_range[0] / COMPRESSION_RATIO);
            //     index = s[0] << 16 | s[1] << 8 | s[2];
            //     index = (index >> shift) & mask;
            //     ...
            // }
            // ```
            if packed_len >= 3 {
                let max_byte_idx = packed_len.saturating_sub(3);
                let max_fast_pos = max_byte_idx * COMPRESSION_RATIO + (COMPRESSION_RATIO - 1);
                let fast_end = end.min(max_fast_pos);
                let mut pos = start;
                while pos <= fast_end {
                    let byte_idx = pos / COMPRESSION_RATIO;
                    let shift = 2 * (12 - ((pos % COMPRESSION_RATIO) + lut_word_length));
                    let idx = ((packed[byte_idx] as u32) << 16)
                        | ((packed[byte_idx + 1] as u32) << 8)
                        | (packed[byte_idx + 2] as u32);
                    let kmer = ((idx >> shift) as u64) & mask;
                    on_kmer(pos, kmer);
                    pos = pos.saturating_add(scan_step);
                }
                while pos <= end {
                    let kmer = packed_kmer_at(packed, pos, lut_word_length);
                    on_kmer(pos, kmer);
                    pos = pos.saturating_add(scan_step);
                }
                return;
            }
        } else {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:245-260
            // ```c
            // for (; scan_range[0] <= scan_range[1]; scan_range[0] += scan_step) {
            //     Int4 shift = 2*(8 - (scan_range[0] % COMPRESSION_RATIO + lut_word_length));
            //     s = abs_start + (scan_range[0] / COMPRESSION_RATIO);
            //     index = s[0] << 8 | s[1];
            //     index = (index >> shift) & mask;
            //     ...
            // }
            // ```
            if packed_len >= 2 {
                let max_byte_idx = packed_len.saturating_sub(2);
                let max_fast_pos = max_byte_idx * COMPRESSION_RATIO + (COMPRESSION_RATIO - 1);
                let fast_end = end.min(max_fast_pos);
                let mut pos = start;
                while pos <= fast_end {
                    let byte_idx = pos / COMPRESSION_RATIO;
                    let shift = 2 * (8 - ((pos % COMPRESSION_RATIO) + lut_word_length));
                    let idx = ((packed[byte_idx] as u16) << 8) | packed[byte_idx + 1] as u16;
                    let kmer = ((idx >> shift) as u64) & mask;
                    on_kmer(pos, kmer);
                    pos = pos.saturating_add(scan_step);
                }
                while pos <= end {
                    let kmer = packed_kmer_at(packed, pos, lut_word_length);
                    on_kmer(pos, kmer);
                    pos = pos.saturating_add(scan_step);
                }
                return;
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1598-1605
    // ```c
    // for (; scan_range[0] <= scan_range[1]; scan_range[0] += scan_step) {
    //     Int4 shift = 2*(16 - (scan_range[0] % COMPRESSION_RATIO + lut_word_length));
    //     s = abs_start + (scan_range[0] / COMPRESSION_RATIO);
    //     index = s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3];
    //     index = (index >> shift) & mask;
    //     MB_ACCESS_HITS();
    // }
    // ```
    // Fallback: advance by scan_step and recompute the k-mer at each position.
    let mut pos = start;
    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, lut_word_length);
        on_kmer(pos, kmer);
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/masksubj.inl:43-58
// ```c
// while (range[1] > range[2]) {
//     range[0]++;
//     if (range[0] >= (Int4)subject->num_seq_ranges) return FALSE;
//     range[1] = subject->seq_ranges[range[0]].left + word_length - lut_word_length;
//     range[2] = subject->seq_ranges[range[0]].right - lut_word_length;
// }
// return TRUE;
// ```
#[inline(always)]
fn determine_scanning_offsets(
    seq_ranges: &[(i32, i32)],
    word_length: i32,
    lut_word_length: i32,
    range: &mut [i32; 3],
) -> bool {
    while range[1] > range[2] {
        range[0] += 1;
        if range[0] >= seq_ranges.len() as i32 {
            return false;
        }
        let (left, right) = seq_ranges[range[0] as usize];
        range[1] = left + word_length - lut_word_length;
        range[2] = right - lut_word_length;
    }
    true
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1673-1676
// ```c
// while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
//     hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
// }
// ```
#[inline(always)]
fn scan_subject_kmers_with_offsets<G>(
    seq_ranges: &[(i32, i32)],
    word_length: usize,
    lut_word_length: usize,
    mut scan_range: [i32; 3],
    mut scan_fn: G,
) where
    G: FnMut(usize, usize),
{
    while determine_scanning_offsets(
        seq_ranges,
        word_length as i32,
        lut_word_length as i32,
        &mut scan_range,
    ) {
        let start = scan_range[1];
        let end = scan_range[2];
        if start >= 0 && end >= 0 {
            scan_fn(start as usize, end as usize);
        }
        scan_range[1] = scan_range[2] + 1;
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1647-1674
// ```c
// scan_range[0] = 0;  /* subject seq mask index */
// scan_range[1] = 0;  /* start pos of scan */
// scan_range[2] = subject->length - lut_word_length;
// while (s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
//     hitsfound = scansub(..., &scan_range[1]);
// }
// ```
fn scan_subject_kmers_with_ranges<F>(
    packed: &[u8],
    subject_len: usize,
    word_length: usize,
    lut_word_length: usize,
    scan_step: usize,
    seq_ranges: &[(i32, i32)],
    subject_masked: bool,
    mb_scan_kind: Option<MbScanSubjectKind>,
    mut on_kmer: F,
) where
    F: FnMut(usize, usize, u64),
{
    if seq_ranges.is_empty() {
        return;
    }

    let mut scan_range = [0i32, 0i32, 0i32];
    if subject_masked {
        let (left, right) = seq_ranges[0];
        scan_range[1] = left + word_length as i32 - lut_word_length as i32;
        scan_range[2] = right - lut_word_length as i32;
    } else {
        scan_range[1] = 0;
        scan_range[2] = subject_len as i32 - lut_word_length as i32;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1635-1676
    // ```c
    // scansub = (TNaScanSubjectFunction)lookup->scansub_callback;
    // ...
    // while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
    //     hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
    //     ...
    // }
    // ```
    match mb_scan_kind {
        Some(MbScanSubjectKind::Any) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_any(
                    packed,
                    subject_len,
                    lut_word_length,
                    scan_step,
                    subject_masked,
                    start,
                    end,
                    &mut |kmer_start, kmer| on_kmer(kmer_start, s_range, kmer),
                );
            },
        ),
        Some(MbScanSubjectKind::Scan9_1) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_9_1(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| on_kmer(kmer_start, s_range, kmer),
                );
            },
        ),
        Some(MbScanSubjectKind::Scan9_2) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_9_2(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| on_kmer(kmer_start, s_range, kmer),
                );
            },
        ),
        Some(MbScanSubjectKind::Scan10_1) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_10_1(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| on_kmer(kmer_start, s_range, kmer),
                );
            },
        ),
        Some(MbScanSubjectKind::Scan10_2) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_10_2(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| on_kmer(kmer_start, s_range, kmer),
                );
            },
        ),
        Some(MbScanSubjectKind::Scan10_3) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_10_3(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| on_kmer(kmer_start, s_range, kmer),
                );
            },
        ),
        Some(MbScanSubjectKind::Scan11_1Mod4) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_11_1mod4(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| on_kmer(kmer_start, s_range, kmer),
                );
            },
        ),
        Some(MbScanSubjectKind::Scan11_2Mod4) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_11_2mod4(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| on_kmer(kmer_start, s_range, kmer),
                );
            },
        ),
        Some(MbScanSubjectKind::Scan11_3Mod4) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_11_3mod4(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| on_kmer(kmer_start, s_range, kmer),
                );
            },
        ),
        None => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range(
                    packed,
                    subject_len,
                    lut_word_length,
                    scan_step,
                    subject_masked,
                    start,
                    end,
                    &mut |kmer_start, kmer| on_kmer(kmer_start, s_range, kmer),
                );
            },
        ),
    }
}

/// Scans a subject for the discontiguous words of a discontiguous megablast table, with
/// the template length as both the word length and the lookup word length; calls
/// `on_word(s_off, s_range, index, index2)` for every word in NCBI's order.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1624-1634
/// ```c
///     else if (lookup_wrap->lut_type == eMBLookupTable) {
///         BlastMBLookupTable *lookup =
///                                 (BlastMBLookupTable *) lookup_wrap->lut;
///         if (lookup->discontiguous) {
///             word_length = lookup->template_length;
///             lut_word_length = lookup->template_length;
///         } else {
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1647-1684
/// ```c
///     scan_range[0] = 0;  /* subject seq mask index */
///     scan_range[1] = 0;	/* start pos of scan */
///     scan_range[2] = subject->length - lut_word_length; /*end pos (inclusive) of scan*/
///
///     /* if sequence is masked, fall back to generic scanner and extender */
///     if (subject->mask_type != eNoSubjMasking) {
///         if (lookup_wrap->lut_type == eMBLookupTable &&
///             ((BlastMBLookupTable *) lookup_wrap->lut)->discontiguous) {
///             /* discontiguous scan subs assumes any (non-aligned starting offset */
///         } else {
///         ...
///         }
///         /* generic scanner permits any (non-aligned) starting offset */
///         scan_range[1] = subject->seq_ranges[0].left + word_length - lut_word_length;
///         scan_range[2] = subject->seq_ranges[0].right - lut_word_length;
///     }
///
///     ASSERT(scansub);
///     ASSERT(extend);
///
///     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
///
///         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
///
///         if (hitsfound == 0)
///             continue;
///
///         total_hits += hitsfound;
///         hits_extended += extend(offset_pairs, hitsfound, word_params,
///                                 lookup_wrap, query, subject, matrix,
///                                 query_info, ewp, init_hitlist, scan_range[2] + lut_word_length);
///     }
/// ```
fn scan_subject_disc_words_with_ranges<F>(
    packed: &[u8],
    subject_len: usize,
    disc: &DiscTemplates,
    seq_ranges: &[(i32, i32)],
    subject_masked: bool,
    mut on_word: F,
) where
    F: FnMut(usize, usize, u64, Option<u64>),
{
    if seq_ranges.is_empty() {
        return;
    }
    let template_length = disc.template_length;
    let mut scan_range = [0i32, 0i32, 0i32];
    if subject_masked {
        let (left, right) = seq_ranges[0];
        scan_range[1] = left;
        scan_range[2] = right - template_length as i32;
    } else {
        scan_range[1] = 0;
        scan_range[2] = subject_len as i32 - template_length as i32;
    }
    let kind = choose_disc_scan_subject(disc.template_type, disc.two_templates);
    scan_subject_kmers_with_offsets(
        seq_ranges,
        template_length,
        template_length,
        scan_range,
        |start, end| {
            let s_range = end.saturating_add(template_length);
            disc_word_scan_subject(
                kind,
                packed,
                disc.template_type,
                disc.second_template_type,
                template_length as i32,
                [start as i32, end as i32],
                &mut |s_off, index, index2| {
                    on_word(s_off as usize, s_range, index as u64, index2.map(u64::from))
                },
            );
        },
    );
}

/// Structure to hold ungapped hit data for batch processing
/// NCBI reference: blast_gapalign.c - init_hsp_array is sorted by score descending
/// before gapped extension
#[derive(Clone, Copy)]
struct UngappedHit {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3913
    // ```c
    // tmp_hsp.context = context;
    // tmp_hsp.query.offset = q_start;
    // tmp_hsp.query.end = q_end;
    // tmp_hsp.query.frame = query_info->contexts[context].frame;
    // ```
    context_idx: u32,
    query_idx: u32,
    query_frame: i32,
    query_context_offset: i32,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:135-149
    // ```c
    // typedef union BlastOffsetPair {
    //     struct { Uint4 q_off; Uint4 s_off; } qs_offsets;
    // } BlastOffsetPair;
    // ```
    // Initial word hit offsets (seed) stored in init_hsp->offsets.qs_offsets.
    seed_q_off: usize,
    seed_s_off: usize,
    // Ungapped extension results (0-based coordinates)
    qs: usize,  // query start
    qe: usize,  // query end (exclusive)
    ss: usize,  // subject start in search_seq coordinates
    se: usize,  // subject end in search_seq coordinates
    score: i32, // ungapped score
}

// No NCBI counterpart: only where LOSAT keeps its buffers (no value changes).
/// EXPERIMENT (LOSAT_X_BNMAPPED, a switch of the strict set, `utils/xswitch.rs`: on unless it is
/// `0`): the large per-subject buffers of BLASTN (the
/// blocks and sequences read from a database and their BLASTNA letters, `database.rs`; the
/// ungapped hits of a subject as the seed stages push them, `x_push_ungapped`, and the scratch of
/// their sort by score, outside the per-query groups of `LOSAT_X_QGROUPPAR` and the merge of the
/// folded strips, `x_pair_seed_part.rs`) are kept in buffers that glibc maps
/// from 1 MiB on, so that they are returned when they are dropped instead of staying in the
/// arena of the worker thread that freed them; the database's buffers and the sort's scratch go
/// to shared pools (`utils::x_mapped::X_BYTES`, `X_SORT`) when they are done, so that the next
/// thread that needs one takes it with its pages resident. The hits are not pooled: a thread
/// keeps its scratch, and with it its hit buffer, from one subject to the next (as NCBI keeps its
/// `init_hitlist`), so idle hit buffers in a pool only added to what was in use (Sakai x EDL933
/// `-subject` blastn at 4 threads, max RSS 298-327 MiB with a pool of 64 MiB, 245-264 without,
/// 238-261 with the switch off).
fn x_bn_mapped() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNMAPPED"))
}

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_MBHITPF, in the strict set: on unless it is `0`) the
/// first extension of a hit reads the query at a place the hit's cell chose, a cache miss for a
/// query of megabases; the loop over the offset pairs prefetches the line of the hit `X_HIT_AHEAD`
/// places on. It changes no value.
fn x_mb_hit_prefetch() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| {
        let on = crate::utils::xswitch::on("LOSAT_X_MBHITPF");
        crate::utils::xreach::x_reach!("LOSAT_X_MBHITPF", on);
        on
    })
}
const X_HIT_AHEAD: usize = 16;

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNPVONCE, in the strict set: on unless it is `0`) the
/// scan of a table laid out as packed CSR tests each word once, in the head's presence bits, when
/// those fit the second-level cache (4^11 cells or fewer); a larger head keeps the small presence
/// vector in front. The words it keeps are the words with hits either way (the presence vector
/// only skips words without hits). It changes no value.
fn x_bn_pv_once(head: &x_dc_csr::XPacked<'_>) -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    let on = *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNPVONCE"));
    let once = on && head.x_bits_bytes() <= X_PV_ONCE_MAX_BITS;
    crate::utils::xreach::x_reach!("LOSAT_X_BNPVONCE", once);
    once
}
const X_PV_ONCE_MAX_BITS: usize = 1 << 20;

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNTBRANGE) the traceback of a database search decodes
/// a sequence of at least `x_tb_range_min` bases in the ranges of its HSPs only, as NCBI's partial
/// fetch does (`database::traceback_listed`).
pub(super) fn x_bn_tb_range() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNTBRANGE"))
}

/// No NCBI counterpart: the shortest sequence that LOSAT_X_BNTBRANGE decodes in ranges, 1 Mi
/// bases (below it the whole BLASTNA is small, and an allocation below glibc's mmap threshold
/// may come from an arena that writes it); LOSAT_X_BNTBRANGEMIN=<bases> sets it (0 = every
/// sequence, for the tests: the fixtures' sequences are short).
pub(super) fn x_tb_range_min() -> usize {
    static MIN: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *MIN.get_or_init(|| {
        std::env::var("LOSAT_X_BNTBRANGEMIN")
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .unwrap_or(crate::utils::x_mapped::X_MAPPED_FROM)
    })
}

/// No NCBI counterpart: records for `LOSAT_X_REACH` whether a buffer of `bytes` takes
/// LOSAT_X_BNMAPPED's own path (mapped, and pooled where it has a pool) or stays as `Vec` keeps
/// it (below 1 MiB, or not glibc).
pub(super) fn x_bn_mapped_reach(bytes: usize) {
    crate::utils::xreach::x_reach!(
        "LOSAT_X_BNMAPPED",
        crate::utils::x_mapped::x_is_mapped(bytes)
    );
}

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNMAPPED) the pool of the scratch of the hits' sort
/// by score (`x_mapped::x_sort_by_mapped`), whose idle buffers may hold `x_sort_idle_bytes`.
/// Without it each sort maps its scratch again and touches every page of it (Sakai x EDL933
/// `-subject` blastn, where the hits of one chunk are sorted at a time: x1.009-1.027 of the time
/// with the switch off).
static X_SORT: crate::utils::x_mapped::XPool<UngappedHit> =
    crate::utils::x_mapped::XPool::new(x_sort_idle_bytes);

/// No NCBI counterpart: the built-in bound of `x_sort_idle_bytes`, in MiB: one scratch of the
/// largest hit lists measured (Sakai x EDL933 `-subject` blastn sorts 10 lists of 1.0-11.4 MiB;
/// with 8 each of the 5 above 9 MiB is mapped again, x1.009 of the time).
const X_SORT_IDLE_MIB: usize = 16;

/// No NCBI counterpart: the idle bound of `X_SORT` (`LOSAT_X_SORTIDLE` MiB, a tuning variable,
/// 0 to 65,536; built-in `X_SORT_IDLE_MIB`).
fn x_sort_idle_bytes() -> usize {
    static IDLE: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *IDLE.get_or_init(|| crate::utils::x_mapped::x_idle_bytes("LOSAT_X_SORTIDLE", X_SORT_IDLE_MIB))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:254-268
// ```c
//     if (init_hitlist->total >= init_hitlist->allocated) {
//         Int4 new_size = MAX(100, 2 * init_hitlist->allocated);
//         ...
//         init_hsp_array = (BlastInitHSP *)
//             realloc(init_hitlist->init_hsp_array, sizeof(BlastInitHSP) * new_size);
// ```
/// `hits.push(hit)`, growing as `Vec` does; with LOSAT_X_BNMAPPED with the capacity of
/// `x_mapped::x_reserve_mapped`.
#[inline(always)]
fn x_push_ungapped(hits: &mut Vec<UngappedHit>, hit: UngappedHit) {
    if hits.len() == hits.capacity() {
        x_grow_ungapped(hits);
    }
    hits.push(hit);
}

#[cold]
#[inline(never)]
fn x_grow_ungapped(hits: &mut Vec<UngappedHit>) {
    if x_bn_mapped() {
        crate::utils::x_mapped::x_reserve_mapped(hits, 1);
    } else {
        hits.reserve(1);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:271-300
// ```c
// static int score_compare_match(const void *v1, const void *v2)
// {
//     ...
//     if (0 == (result = BLAST_CMP(h2->ungapped_data->score,
//                                  h1->ungapped_data->score)) &&
//         0 == (result = BLAST_CMP(h1->ungapped_data->s_start,
//                                  h2->ungapped_data->s_start)) &&
//         0 == (result = BLAST_CMP(h2->ungapped_data->length,
//                                  h1->ungapped_data->length)) &&
//         0 == (result = BLAST_CMP(h1->ungapped_data->q_start,
//                                  h2->ungapped_data->q_start))) {
//         result = BLAST_CMP(h2->ungapped_data->length,
//                            h1->ungapped_data->length);
//     }
// }
// ```
// NCBI's `ungapped_data->q_start` is an offset in the concatenated query (the
// scan and the ungapped extension run on the whole query block), so hits of
// different contexts never tie on it; `qs` here is context-local, and the
// context offset is added back for the comparison.
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:206-207
// ```c
// ungapped_data->q_start = (Int4)(q_beg - query->sequence);
// ungapped_data->s_start = s_off - (q_off - ungapped_data->q_start);
// ```
fn score_compare_ungapped_hits(a: &UngappedHit, b: &UngappedHit) -> std::cmp::Ordering {
    let a_len = a.qe.saturating_sub(a.qs);
    let b_len = b.qe.saturating_sub(b.qs);
    let a_q_start = i64::from(a.query_context_offset) + a.qs as i64;
    let b_q_start = i64::from(b.query_context_offset) + b.qs as i64;

    b.score
        .cmp(&a.score)
        .then_with(|| a.ss.cmp(&b.ss))
        .then_with(|| b_len.cmp(&a_len))
        .then_with(|| a_q_start.cmp(&b_q_start))
        .then_with(|| b_len.cmp(&a_len))
}

/// Preliminary gapped HSP data for greedy traceback.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4031
/// ```c
/// } else if (is_greedy) {
///    if (init_hsp->ungapped_data) {
///        init_hsp->offsets.qs_offsets.q_off =
///            init_hsp->ungapped_data->q_start + init_hsp->ungapped_data->length/2;
///        init_hsp->offsets.qs_offsets.s_off =
///            init_hsp->ungapped_data->s_start + init_hsp->ungapped_data->length/2;
///    }
///    status = BLAST_GreedyGappedAlignment(..., (Boolean) TRUE, FALSE, fence_hit);
///    init_hsp->offsets.qs_offsets.q_off = gap_align->greedy_query_seed_start;
///    init_hsp->offsets.qs_offsets.s_off = gap_align->greedy_subject_seed_start;
/// }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4058-4076
/// ```c
/// if (gap_align->score >= cutoff) {
///    status = Blast_HSPInit(gap_align->query_start,
///              gap_align->query_stop, gap_align->subject_start,
///              gap_align->subject_stop,
///              init_hsp->offsets.qs_offsets.q_off,
///              init_hsp->offsets.qs_offsets.s_off, context,
///              query_frame, subject->frame, gap_align->score,
///              &(gap_align->edit_script), &new_hsp);
/// }
/// ```
#[derive(Clone, Debug)]
struct PrelimHit {
    context_idx: u32,
    query_idx: u32,
    query_frame: i32,
    query_context_offset: i32,
    prelim_qs: usize,
    prelim_qe: usize,
    prelim_ss: usize,
    prelim_se: usize,
    prelim_score: i32,
    seed_qs: usize,
    seed_ss: usize,
    /// The e-value of the preliminary stage (set before its reap); a merge of two HSPs keeps
    /// the first one's (`s_BlastMergeTwoHSPs` does not change `evalue`).
    prelim_evalue: f64,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
// ```c
// if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
//     0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
//     0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
//     0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
//     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
// }
// ```
fn score_compare_prelim_hits(a: &PrelimHit, b: &PrelimHit) -> std::cmp::Ordering {
    b.prelim_score
        .cmp(&a.prelim_score)
        .then_with(|| a.prelim_ss.cmp(&b.prelim_ss))
        .then_with(|| b.prelim_se.cmp(&a.prelim_se))
        .then_with(|| a.prelim_qs.cmp(&b.prelim_qs))
        .then_with(|| b.prelim_qe.cmp(&a.prelim_qe))
}

type PrelimHitCompare = fn(&PrelimHit, &PrelimHit) -> std::cmp::Ordering;

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2321
// ```c
// if (h1->context < h2->context) return -1;
// if (h1->query.offset < h2->query.offset) return -1;
// if (h1->subject.offset < h2->subject.offset) return -1;
// if (h1->score < h2->score) return 1;
// if (h1->query.end < h2->query.end) return 1;
// if (h1->subject.end < h2->subject.end) return 1;
// ```
fn query_offset_compare_prelim_hits(a: &PrelimHit, b: &PrelimHit) -> std::cmp::Ordering {
    a.context_idx
        .cmp(&b.context_idx)
        .then_with(|| a.prelim_qs.cmp(&b.prelim_qs))
        .then_with(|| a.prelim_ss.cmp(&b.prelim_ss))
        .then_with(|| b.prelim_score.cmp(&a.prelim_score))
        .then_with(|| b.prelim_qe.cmp(&a.prelim_qe))
        .then_with(|| b.prelim_se.cmp(&a.prelim_se))
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2327-2383
// ```c
// if (h1->context < h2->context) return -1;
// if (h1->query.end < h2->query.end) return -1;
// if (h1->subject.end < h2->subject.end) return -1;
// if (h1->score < h2->score) return 1;
// if (h1->query.offset < h2->query.offset) return 1;
// if (h1->subject.offset < h2->subject.offset) return 1;
// ```
fn query_end_compare_prelim_hits(a: &PrelimHit, b: &PrelimHit) -> std::cmp::Ordering {
    a.context_idx
        .cmp(&b.context_idx)
        .then_with(|| a.prelim_qe.cmp(&b.prelim_qe))
        .then_with(|| a.prelim_se.cmp(&b.prelim_se))
        .then_with(|| b.prelim_score.cmp(&a.prelim_score))
        .then_with(|| b.prelim_qs.cmp(&a.prelim_qs))
        .then_with(|| b.prelim_ss.cmp(&a.prelim_ss))
}

fn qsort_prelim_hits_by(hsps: &mut [PrelimHit], compare: PrelimHitCompare) {
    if hsps.len() <= 1 {
        return;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
    // ```c
    // qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //       ScoreCompareHSPs);
    // ```
    // Keep the NCBI comparator fields and call timing while using the same
    // stable Rust sort on every target.
    hsps.sort_by(compare);
}

fn sort_prelim_hits_by_score(hsps: &mut [PrelimHit]) {
    if hsps.len() <= 1 {
        return;
    }

    let sorted = hsps
        .windows(2)
        .all(|pair| score_compare_prelim_hits(&pair[0], &pair[1]) != std::cmp::Ordering::Greater);
    if sorted {
        return;
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:555
    // ```c
    // Blast_HSPListSortByScore(hsp_list);
    // ```
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1374-1382
    // ```c
    // void Blast_HSPListSortByScore(BlastHSPList* hsp_list)
    // {
    //     if (!hsp_list || hsp_list->hspcnt <= 1)
    //         return;
    //     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
    //         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //               ScoreCompareHSPs);
    //     }
    // }
    // ```
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:358-365
    // ```c
    // #ifdef _DEBUG
    // { Blast_HSPListSortByScore(hsp_list); }
    // #endif
    // ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
    // ```
    // Preliminary HSPs drive traceback interval-tree insertion order. The
    // target-neutral stable sort preserves incoming order for comparator ties.
    qsort_prelim_hits_by(hsps, score_compare_prelim_hits);
}

fn purge_prelim_hits_with_common_endpoints(mut hsps: Vec<PrelimHit>) -> Vec<PrelimHit> {
    if hsps.len() <= 1 {
        return hsps;
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:540-557
    // ```c
    // if (aux_struct->GetGappedScore) {
    //     /* Removes redundant HSPs. */
    //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    // }
    // Blast_HSPListSortByScore(hsp_list);
    // ```
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
    // ```c
    // purge |= (program != eBlastTypeBlastn);
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
    // ...
    // if (!purge && (hsp->query.end > hsp_array[i]->query.end)) {
    //     s_CutOffGapEditScript(...);
    // } else {
    //     hsp = Blast_HSPFree(hsp);
    // }
    // ...
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
    // ```
    // Preliminary BLASTN HSPs use purge=TRUE, so common-start/common-end
    // duplicates are deleted, not trimmed. The full traceback stage still uses
    // the later purge=FALSE/TRUE passes over final HSPs with edit scripts.
    qsort_prelim_hits_by(&mut hsps, query_offset_compare_prelim_hits);
    // In place: `dedup_by` drops an HSP that matches the last one kept, as the copy into a
    // second vector did.
    hsps.dedup_by(|hsp, prev| {
        prev.context_idx == hsp.context_idx
            && prev.prelim_qs == hsp.prelim_qs
            && prev.prelim_ss == hsp.prelim_ss
    });
    let mut kept = hsps;

    if kept.len() <= 1 {
        return kept;
    }

    qsort_prelim_hits_by(&mut kept, query_end_compare_prelim_hits);
    kept.dedup_by(|hsp, prev| {
        prev.context_idx == hsp.context_idx
            && prev.prelim_qe == hsp.prelim_qe
            && prev.prelim_se == hsp.prelim_se
    });

    kept
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits_priv.h:66-72
// ```c
// #define CONTAINED_IN_HSP(a,b,c,d,e,f) \
//     (((a <= c && b >= c) && (d <= f && e >= f)) ? TRUE : FALSE)
// ```
#[inline]
fn contained_in_hsp(a: usize, b: usize, c: usize, d: usize, e: usize, f: usize) -> bool {
    a <= c && b >= c && d <= f && e >= f
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1488-1533
// ```c
// static Boolean s_BlastMergeTwoHSPs(BlastHSP* hsp1, BlastHSP* hsp2, Boolean allow_gap)
// ```
fn merge_two_prelim_hits(hsp1: &mut PrelimHit, hsp2: &PrelimHit, allow_gap: bool) -> bool {
    if !allow_gap
        && (hsp1.prelim_ss as isize - hsp2.prelim_ss as isize - hsp1.prelim_qs as isize
            + hsp2.prelim_qs as isize)
            != 0
    {
        return false;
    }

    // BLASTN subject frame is always +1, so no frame mismatch handling needed here.

    if contained_in_hsp(
        hsp1.prelim_qs,
        hsp1.prelim_qe,
        hsp2.prelim_qs,
        hsp1.prelim_ss,
        hsp1.prelim_se,
        hsp2.prelim_ss,
    ) || contained_in_hsp(
        hsp1.prelim_qs,
        hsp1.prelim_qe,
        hsp2.prelim_qe,
        hsp1.prelim_ss,
        hsp1.prelim_se,
        hsp2.prelim_se,
    ) {
        let len1 = hsp1.prelim_qe.saturating_sub(hsp1.prelim_qs) as f64;
        let len2 = hsp2.prelim_qe.saturating_sub(hsp2.prelim_qs) as f64;
        let score_density = (hsp1.prelim_score as f64 + hsp2.prelim_score as f64) / (len1 + len2);

        hsp1.prelim_qs = hsp1.prelim_qs.min(hsp2.prelim_qs);
        hsp1.prelim_ss = hsp1.prelim_ss.min(hsp2.prelim_ss);
        hsp1.prelim_qe = hsp1.prelim_qe.max(hsp2.prelim_qe);
        hsp1.prelim_se = hsp1.prelim_se.max(hsp2.prelim_se);

        if hsp2.prelim_score > hsp1.prelim_score {
            hsp1.seed_qs = hsp2.seed_qs;
            hsp1.seed_ss = hsp2.seed_ss;
            hsp1.prelim_score = hsp2.prelim_score;
        }

        let new_len = hsp1.prelim_qe.saturating_sub(hsp1.prelim_qs) as f64;
        let score_from_density = (score_density * new_len) as i32;
        if score_from_density > hsp1.prelim_score {
            hsp1.prelim_score = score_from_density;
        }
        return true;
    }

    false
}

/// NCBI's preliminary subject-best-hit filter (`-subject_besthit`) over the combined HSP
/// list of a subject, applied after each chunk is merged, before traceback.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:586-591
/// ```c
///         if((hit_params->options->hsp_filt_opt != NULL) &&
///            (hit_params->options->hsp_filt_opt->subject_besthit_opts != NULL)) {
///            	Blast_HSPListSubjectBestHit(program_number,
///            								hit_params->options->hsp_filt_opt->subject_besthit_opts,
///            								query_info, combined_hsp_list);
///         }
/// ```
fn prelim_subject_best_hit(combined: &mut Vec<PrelimHit>, query_lengths: &[usize]) {
    subject_best_hit_by(combined, |hit| BestHitKey {
        context: hit.context_idx,
        query_frame: hit.query_frame,
        query_offset: hit.prelim_qs,
        query_end: hit.prelim_qe,
        query_length: query_lengths[hit.query_idx as usize],
    });
}

/// Where the second of two preliminary HSP lists continues the first: the start of a subject
/// chunk, or the offsets of a query chunk's contexts (plus and minus strand; `Blast_HSPListsMerge`
/// takes `contexts_per_query < 0` for the subject).
#[derive(Clone, Copy)]
enum HspListSplit {
    Subject { offset: usize },
    Query { offsets: [i32; 2] },
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2857-2862
// ```c
// Int2 Blast_HSPListsMerge(BlastHSPList** hsp_list_ptr,
//                    BlastHSPList** combined_hsp_list_ptr,
//                    Int4 hsp_num_max, Int4 *split_offsets,
//                    Int4 contexts_per_query, Int4 chunk_overlap_size,
//                    Boolean allow_gap, Boolean short_reads)
// ```
// `hsp_num_max` is `INT4_MAX` for BLASTN (`BlastHspNumMax`), so every HSP is kept.
fn merge_prelim_hit_lists(
    combined: &mut Vec<PrelimHit>,
    mut incoming: Vec<PrelimHit>,
    split: HspListSplit,
    chunk_overlap_size: usize,
    allow_gap: bool,
) -> Vec<PrelimHit> {
    if incoming.is_empty() {
        return incoming;
    }
    if combined.is_empty() {
        combined.extend(incoming.drain(..));
        return incoming;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2890-2945
    // ```c
    //    if (contexts_per_query < 0) {      /* subject seq is split */
    //       for (index1 = 0; index1 < combined_hsp_list->hspcnt; index1++) {
    //          hsp1 = combined_hsp_list->hsp_array[index1];
    //          if (hsp1->subject.end > split_offsets[0]) {
    //             /* At least part of this HSP lies in the overlap strip. */
    //             hsp_var = combined_hsp_list->hsp_array[hspcnt1];
    //             combined_hsp_list->hsp_array[hspcnt1] = hsp1;
    //             combined_hsp_list->hsp_array[index1] = hsp_var;
    //             ++hspcnt1;
    //          }
    //       }
    //       ...
    //          if (hsp2->subject.offset < split_offsets[0] + chunk_overlap_size) {
    //       ...
    //    else {            /* query seq is split */
    //       ...
    //          offset_idx = hsp1->context % contexts_per_query;
    //          if (split_offsets[offset_idx] < 0) continue;
    //          if ((hsp1->query.frame >= 0 && hsp1->query.end >
    //                          split_offsets[offset_idx]) ||
    //              (hsp1->query.frame < 0 && hsp1->query.offset <
    //                          split_offsets[offset_idx] + chunk_overlap_size)) {
    //       ...
    //          if ((hsp2->query.frame < 0 && hsp2->query.end >
    //                          split_offsets[offset_idx]) ||
    //              (hsp2->query.frame >= 0 && hsp2->query.offset <
    //                          split_offsets[offset_idx] + chunk_overlap_size)) {
    // ```
    // The swaps put the HSPs of the overlap strip first, in order, and reorder the others;
    // the stable sort by score below keeps that order for equal HSPs.
    let overlap = chunk_overlap_size as i64;
    let in_strip = |hsp: &PrelimHit, combined_side: bool| match split {
        HspListSplit::Subject { offset } => {
            if combined_side {
                hsp.prelim_se > offset
            } else {
                hsp.prelim_ss < offset.saturating_add(chunk_overlap_size)
            }
        }
        HspListSplit::Query { offsets } => {
            let split_offset = offsets[(hsp.context_idx % 2) as usize] as i64;
            if split_offset < 0 {
                return false;
            }
            let ends_after = hsp.prelim_qe as i64 > split_offset;
            let starts_before = (hsp.prelim_qs as i64) < split_offset + overlap;
            if (hsp.query_frame >= 0) == combined_side {
                ends_after
            } else {
                starts_before
            }
        }
    };
    let mut hspcnt1 = 0usize;
    for index1 in 0..combined.len() {
        if in_strip(&combined[index1], true) {
            combined.swap(hspcnt1, index1);
            hspcnt1 += 1;
        }
    }
    let mut hspcnt2 = 0usize;
    for index2 in 0..incoming.len() {
        if in_strip(&incoming[index2], false) {
            incoming.swap(hspcnt2, index2);
            hspcnt2 += 1;
        }
    }

    if hspcnt1 > 0 && hspcnt2 > 0 {
        let mut deleted = vec![false; incoming.len()];
        for index1 in 0..hspcnt1 {
            let hsp1_context = combined[index1].context_idx;
            for index2 in 0..hspcnt2 {
                if deleted[index2] || incoming[index2].context_idx != hsp1_context {
                    continue;
                }
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2982-2990
                // ```c
                // if (contexts_per_query < 0 || hsp1->query.frame >= 0) {
                //    end_diag = s_HSPEndDiag(hsp1);
                //    start_diag = s_HSPStartDiag(hsp2);
                // }
                // else {
                //    end_diag = s_HSPEndDiag(hsp2);
                //    start_diag = s_HSPStartDiag(hsp1);
                // }
                // ```
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1465-1477
                // ```c
                //     return hsp->query.offset - hsp->subject.offset;
                //     return hsp->query.end - hsp->subject.end;
                // ```
                // `hsp1` may already have been expanded by an earlier merge in this
                // inner loop, so recompute its diagonals every iteration.
                let end_diag = |hsp: &PrelimHit| hsp.prelim_qe as i64 - hsp.prelim_se as i64;
                let start_diag = |hsp: &PrelimHit| hsp.prelim_qs as i64 - hsp.prelim_ss as i64;
                let hsp1_first = matches!(split, HspListSplit::Subject { .. })
                    || combined[index1].query_frame >= 0;
                let (end, start) = if hsp1_first {
                    (end_diag(&combined[index1]), start_diag(&incoming[index2]))
                } else {
                    (end_diag(&incoming[index2]), start_diag(&combined[index1]))
                };
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2991-2995
                // ```c
                //             if (ABS(end_diag - start_diag) < OVERLAP_DIAG_CLOSE) {
                //                if (s_BlastMergeTwoHSPs(hsp1, hsp2, allow_gap)) {
                //                   /* Free the second HSP. */
                //                   hspp2[index2] = Blast_HSPFree(hsp2);
                //                }
                // ```
                if (end - start).abs() < OVERLAP_DIAG_CLOSE as i64
                    && merge_two_prelim_hits(&mut combined[index1], &incoming[index2], allow_gap)
                {
                    deleted[index2] = true;
                }
            }
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3000-3001
        // ```c
        // /* Purge the nulled out HSPs from the new HSP list */
        // Blast_HSPListPurgeNullHSPs(hsp_list);
        // ```
        let mut delete_idx = 0usize;
        incoming.retain(|_| {
            let keep = !deleted[delete_idx];
            delete_idx += 1;
            keep
        });
    }

    append_prelim_hit_list(combined, incoming)
}

/// NCBI's `s_BlastHSPListsCombineByScore` when every HSP is kept (`hsp_num_max` is
/// `INT4_MAX`): the HSPs of the second list follow those of the first, then the list is sorted
/// by score. `Blast_HSPListAppend` is this alone.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2757-2766
/// ```c
///    if (new_hspcnt >= hsp_list->hspcnt + combined_hsp_list->hspcnt) {
///       /* All HSPs from both arrays are saved */
///       for (index=combined_hsp_list->hspcnt, index1=0;
///            index1<hsp_list->hspcnt; index1++) {
///          if (hsp_list->hsp_array[index1] != NULL)
///             combined_hsp_list->hsp_array[index++] = hsp_list->hsp_array[index1];
///       }
///       combined_hsp_list->hspcnt = new_hspcnt;
///       Blast_HSPListSortByScore(combined_hsp_list);
/// ```
fn append_prelim_hit_list(
    combined: &mut Vec<PrelimHit>,
    mut incoming: Vec<PrelimHit>,
) -> Vec<PrelimHit> {
    if incoming.is_empty() {
        return incoming;
    }
    if combined.is_empty() {
        combined.extend(incoming.drain(..));
        return incoming;
    }
    combined.extend(incoming.drain(..));
    sort_prelim_hits_by_score(combined);
    incoming
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1414-1435
// ```c
// s_EvalueCompareHSPs(const void* v1, const void* v2)
// {
//    ...
//    if ((retval = s_EvalueComp(h1->evalue, h2->evalue)) != 0)
//       return retval;
//
//    return ScoreCompareHSPs(v1, v2);
// }
// ```
fn evalue_compare_prelim_hits(a: &PrelimHit, b: &PrelimHit) -> std::cmp::Ordering {
    evalue_comp(a.prelim_evalue, b.prelim_evalue).then_with(|| score_compare_prelim_hits(a, b))
}

/// A subject's preliminary HSPs of one query: NCBI's `BlastHSPList` in the HSP stream of the
/// preliminary stage, kept in the query's hit list of `prelim_hitlist_size` subjects.
struct PrelimHspList {
    oid: u32,
    hsps: Vec<PrelimHit>,
    best_evalue: f64,
}

impl HitListEntry for PrelimHspList {
    fn oid(&self) -> u32 {
        self.oid
    }

    fn hsp_count(&self) -> usize {
        self.hsps.len()
    }

    fn best_evalue(&self) -> f64 {
        self.best_evalue
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1740-1748
    // ```c
    // s_BlastGetBestEvalue(const BlastHSPList* hsp_list)
    // {
    //     int index = 0;
    //     double best_evalue = (double) INT4_MAX;
    //
    //     for (index=0; index<hsp_list->hspcnt; index++)
    //        best_evalue = MIN(hsp_list->hsp_array[index]->evalue, best_evalue);
    //
    //     return best_evalue;
    // ```
    fn update_best_evalue(&mut self) {
        let mut best = i32::MAX as f64;
        for hsp in &self.hsps {
            if hsp.prelim_evalue < best {
                best = hsp.prelim_evalue;
            }
        }
        self.best_evalue = best;
    }

    fn first_score(&self) -> Option<i32> {
        self.hsps.first().map(|hsp| hsp.prelim_score)
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1442-1455
    // ```c
    //     if (hsp_list->hspcnt > 1) {
    //         Int4 index;
    //         BlastHSP** hsp_array = hsp_list->hsp_array;
    //         /* First check if HSP array is already sorted. */
    //         for (index = 0; index < hsp_list->hspcnt - 1; ++index) {
    //             if (s_EvalueCompareHSPs(&hsp_array[index], &hsp_array[index+1]) > 0) {
    //                 break;
    //             }
    //         }
    //         /* Sort the HSP array if it is not sorted yet. */
    //         if (index < hsp_list->hspcnt - 1) {
    //             qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //                   s_EvalueCompareHSPs);
    //         }
    // ```
    fn sort_by_evalue(&mut self) {
        let sorted = self.hsps.windows(2).all(|pair| {
            evalue_compare_prelim_hits(&pair[0], &pair[1]) != std::cmp::Ordering::Greater
        });
        if !sorted {
            qsort_prelim_hits_by(&mut self.hsps, evalue_compare_prelim_hits);
        }
    }
}

/// The hit lists of the preliminary stage, one per query.
type PrelimHitLists = Vec<Option<HitList<PrelimHspList>>>;

/// NCBI's HSP collector of the preliminary stage: each subject's HSPs (in subject order) are
/// split into one list per query, in their order, and saved in the query's hit list, which
/// holds at most `prelim_hitlist_size` subjects (`Blast_HitListUpdate` drops the list with
/// the worst preliminary e-value).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_collector.c:116-161
/// ```c
///       for (index = 0; index < hsp_list->hspcnt; index++) {
///          Int4 query_index;
///          hsp = hsp_list->hsp_array[index];
///          query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
///          ...
///          Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
///       ...
///       for (index = 0; index < results->num_queries; index++) {
///          if (hsp_list_array[index]) {
///             if (!results->hitlist_array[index]) {
///                results->hitlist_array[index] =
///                   Blast_HitListNew(params->prelim_hitlist_size);
///             }
///             Blast_HitListUpdate(results->hitlist_array[index],
///                                 hsp_list_array[index]);
///          }
///       }
///       ...
///    } else if (hsp_list->hspcnt > 0) {
///       /* Single query; save the HSP list directly into the results
///          structure */
///       if (!results->hitlist_array[0]) {
///          results->hitlist_array[0] =
///             Blast_HitListNew(params->prelim_hitlist_size);
///       }
///       Blast_HitListUpdate(results->hitlist_array[0], hsp_list);
/// ```
fn collect_prelim_hit_lists(
    subject_hits: Vec<Vec<PrelimHit>>,
    num_queries: usize,
    prelim_hitlist_size: usize,
) -> PrelimHitLists {
    let mut collector = PrelimHitCollector::new(num_queries, prelim_hitlist_size);
    for (oid, hits) in subject_hits.into_iter().enumerate() {
        collector.add(oid, hits);
    }
    collector.hit_lists
}

/// EXPERIMENT (LOSAT_X_BNSTREAM, a strict-set switch: on unless it is `0`)
/// : with two or more threads, the subjects of a
/// query batch are claimed by the threads in OID order and their results go to the hit lists in
/// OID order while the threads search the next subjects (`x_bnstream::search_in_order`); a
/// thread claims a subject only while it is fewer than `X_BNSTREAM_SUBJECTS_PER_THREAD`
/// subjects per thread past the last one handed over, so that the HSPs of that window are held
/// at a time instead of those of every subject. The hit lists are the same.
// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes.
fn x_bnstream() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNSTREAM"))
}

/// No NCBI counterpart: the subjects per thread of the window of `LOSAT_X_BNSTREAM`.
const X_BNSTREAM_SUBJECTS_PER_THREAD: usize = 64;

/// EXPERIMENT (LOSAT_X_BNHITGROUP, a strict-set switch: on unless it is `0`): a subject's HSPs
/// go to lists kept by the collector and only the queries of those HSPs are visited, instead of a
/// new list per query of the batch, each visited: that cost a 24-byte cell and a visit per query
/// for every subject with an HSP, on the thread that takes the subjects in OID order (reads 1e6 x
/// 16S, 33,000 queries per batch x 28,036 subjects: 53 % of the cycles at 8 threads). The hit
/// lists are the same.
// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes. NCBI
// visits every query per subject (hspfilter_collector.c, above); AGENTS.md rule 9 allows less work.
fn x_bnhitgroup() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNHITGROUP"))
}

/// The preliminary hit lists of a batch filled one subject at a time, in OID order
/// (`collect_prelim_hit_lists`), so that only the HSP lists the hit lists keep stay in memory
/// instead of those of every subject.
struct PrelimHitCollector {
    hit_lists: PrelimHitLists,
    prelim_hitlist_size: usize,
    /// LOSAT_X_BNHITGROUP: one empty list per query from the first subject with an HSP, and the
    /// queries of the subject being added.
    x_slots: Vec<Vec<PrelimHit>>,
    x_touched: Vec<u32>,
}

impl PrelimHitCollector {
    fn new(num_queries: usize, prelim_hitlist_size: usize) -> Self {
        let mut hit_lists: PrelimHitLists = Vec::with_capacity(num_queries);
        hit_lists.resize_with(num_queries, || None);
        Self {
            hit_lists,
            prelim_hitlist_size,
            x_slots: Vec::new(),
            x_touched: Vec::new(),
        }
    }

    /// The HSPs of subject `oid`; the subjects must come in OID order.
    fn add(&mut self, oid: usize, hits: Vec<PrelimHit>) {
        if hits.is_empty() {
            return;
        }
        if x_bnhitgroup() {
            self.x_add_grouped(oid, hits);
        } else {
            self.add_per_query(oid, hits);
        }
    }

    /// `LOSAT_X_BNHITGROUP=0`: a list per query of the batch, as NCBI's `hsp_list_array`.
    fn add_per_query(&mut self, oid: usize, hits: Vec<PrelimHit>) {
        let mut per_query: Vec<Vec<PrelimHit>> = vec![Vec::new(); self.hit_lists.len()];
        for hit in hits {
            per_query[hit.query_idx as usize].push(hit);
        }
        for (query, hsps) in per_query.into_iter().enumerate() {
            if !hsps.is_empty() {
                self.update(query, oid, hsps);
            }
        }
    }

    /// EXPERIMENT (LOSAT_X_BNHITGROUP): the same lists, made in the collector's slots: only the
    /// queries of the subject's HSPs are visited. A query's hit list does not read another's, so
    /// that the queries may come in the order of their first HSPs; the HSPs of a query keep their
    /// order. A slot gives its list away, so that the lists grow as those of `add_per_query`.
    fn x_add_grouped(&mut self, oid: usize, hits: Vec<PrelimHit>) {
        if self.x_slots.is_empty() {
            self.x_slots.resize_with(self.hit_lists.len(), Vec::new);
        }
        let mut touched = std::mem::take(&mut self.x_touched);
        for hit in hits {
            let slot = &mut self.x_slots[hit.query_idx as usize];
            if slot.is_empty() {
                touched.push(hit.query_idx);
            }
            slot.push(hit);
        }
        for &query in &touched {
            let hsps = std::mem::take(&mut self.x_slots[query as usize]);
            self.update(query as usize, oid, hsps);
        }
        touched.clear();
        self.x_touched = touched;
    }

    fn update(&mut self, query: usize, oid: usize, hsps: Vec<PrelimHit>) {
        let prelim_hitlist_size = self.prelim_hitlist_size;
        self.hit_lists[query]
            .get_or_insert_with(|| HitList::new(prelim_hitlist_size))
            .update(PrelimHspList {
                oid: oid as u32,
                hsps,
                best_evalue: 0.0,
            });
    }
}

/// The preliminary HSPs that the traceback reads, for each subject: the kept lists of the
/// queries, a query after another. NCBI's `BlastHSPStreamClose` sorts the kept lists by
/// subject, and the traceback reads each subject's lists (`BlastHSPStreamBatchRead`); the
/// traceback of one query's list does not depend on the lists of the other queries.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:196-203
/// ```c
///    /* sort in order of decreasing subject OID. HSPLists will be
///       read out from the end of hsplist_array later */
///
///    hsp_stream->num_hsplists = num_hsplists;
///    if (num_hsplists > 1) {
///       qsort(hsp_stream->sorted_hsplists, num_hsplists,
///                     sizeof(BlastHSPList *), s_SortHSPListByOid);
///    }
/// ```
fn prelim_hits_by_subject(hit_lists: PrelimHitLists, num_subjects: usize) -> Vec<Vec<PrelimHit>> {
    let mut by_subject: Vec<Vec<PrelimHit>> = Vec::with_capacity(num_subjects);
    by_subject.resize_with(num_subjects, Vec::new);
    for hit_list in hit_lists.into_iter().flatten() {
        for list in hit_list.hsplist_array {
            by_subject[list.oid as usize].extend(list.hsps);
        }
    }
    by_subject
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:196-203
// ```c
//    /* sort in order of decreasing subject OID. HSPLists will be
//       read out from the end of hsplist_array later */
//
//    hsp_stream->num_hsplists = num_hsplists;
//    if (num_hsplists > 1) {
//       qsort(hsp_stream->sorted_hsplists, num_hsplists,
//                     sizeof(BlastHSPList *), s_SortHSPListByOid);
//    }
// ```
/// `prelim_hits_by_subject` for a database: the ordinal ids of the kept lists, ascending,
/// each with its kept HSPs in the same order (the lists of the queries, a query after
/// another), so that the traceback reads only these sequences
/// (`database::traceback_listed`).
fn listed_prelim_hits(hit_lists: PrelimHitLists) -> Vec<(u32, Vec<PrelimHit>)> {
    let mut by_subject: std::collections::BTreeMap<u32, Vec<PrelimHit>> =
        std::collections::BTreeMap::new();
    for hit_list in hit_lists.into_iter().flatten() {
        for list in hit_list.hsplist_array {
            by_subject.entry(list.oid).or_default().extend(list.hsps);
        }
    }
    by_subject.into_iter().collect()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
// ```c
// hsp_list = Blast_HSPListFree(hsp_list);
// BlastInitHitListReset(init_hitlist);
// ```
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
struct SubjectScratch {
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    double evalue;
    //    BlastSeg query;
    //    BlastSeg subject;
    //    Int4 context;
    // } BlastHSP;
    // ```
    hits: Vec<BlastnHsp>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4076
    // ```c
    // status = BLAST_GreedyGappedAlignment(..., (Boolean) TRUE, FALSE, fence_hit);
    // init_hsp->offsets.qs_offsets.q_off = gap_align->greedy_query_seed_start;
    // init_hsp->offsets.qs_offsets.s_off = gap_align->greedy_subject_seed_start;
    // status = Blast_HSPInit(gap_align->query_start,
    //           gap_align->query_stop, gap_align->subject_start,
    //           gap_align->subject_stop,
    //           init_hsp->offsets.qs_offsets.q_off,
    //           init_hsp->offsets.qs_offsets.s_off, context,
    //           query_frame, subject->frame, gap_align->score,
    //           &(gap_align->edit_script), &new_hsp);
    // ```
    prelim_hits: Vec<PrelimHit>,
    ungapped_hits: Vec<UngappedHit>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:356-357
    // ```c
    // gap_align->fwd_prelim_tback = GapPrelimEditBlockNew();
    // gap_align->rev_prelim_tback = GapPrelimEditBlockNew();
    // ```
    greedy_align_scratch: GreedyAlignScratch,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:356-357
    // ```c
    //    gap_align->fwd_prelim_tback = GapPrelimEditBlockNew();
    //    gap_align->rev_prelim_tback = GapPrelimEditBlockNew();
    // ```
    // Every look-ahead helper has its own scratch, as each NCBI call has the `gap_align` structure
    // that holds its working memory. A helper never uses the caller's scratch. Reuse only.
    // No NCBI counterpart: scratch memory of the look-ahead helpers.
    x_ahead_scratch: Vec<(GapAlignScratch, GreedyAlignScratch)>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:509-513
    // ```c
    //             BLAST_GappedAlignmentWithTraceback(program_number,
    //                 query, adjusted_subject, gap_align, score_params,
    // ```
    // No NCBI counterpart: the scratch of the speculative traceback batches, one pair per thread
    // of the pool, kept with the rest of this scratch so that each subject's traceback reuses it
    // instead of making it again. Reuse only, as `x_ahead_scratch`.
    x_spec_scratch: Vec<(GapAlignScratch, GreedyAlignScratch)>,
    cutoff_scores: Vec<i32>,
    x_dropoff_scores: Vec<i32>,
    reduced_cutoff_scores: Vec<i32>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:184-198
    // ```c
    // if (backup->allocated >= num_seq_ranges) return;
    // if (backup->allocated) {
    //     sfree(subject->seq_ranges);
    // }
    // backup->allocated = num_seq_ranges;
    // subject->seq_ranges = (SSeqRange *) calloc(backup->allocated,
    //                                sizeof(SSeqRange));
    // ```
    seq_ranges_scratch: Vec<(i32, i32)>,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:77-81
    // ```c
    // typedef struct BLAST_DiagTable {
    //    DiagStruct* hit_level_array;
    //    Uint1* hit_len_array;
    //    Int4 diag_array_length;
    // } BLAST_DiagTable;
    // ```
    hit_level_array: Vec<DiagStruct>,
    hit_len_array: Vec<u8>,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:97-105
    // ```c
    // typedef struct BLAST_DiagHash {
    //    Uint4 num_buckets;
    //    Uint4 occupancy;
    //    Uint4 capacity;
    //    Uint4 *backbone;
    //    DiagHashCell *chain;
    //    Int4 offset;
    //    Int4 window;
    // } BLAST_DiagHash;
    // ```
    diag_hash: DiagHashTable,
    // EXPERIMENT (LOSAT_X_PAIRPARFOLD, `x_pair_seed_part.rs`): the allocations of the partition
    // tables and buffers of the partitioned fold of the strips paths. Every partitioned stage
    // copies `diag_hash` into the partition tables, so no content carries over from one stage
    // to the next. Empty while the switch is off. No NCBI counterpart.
    x_pair_fold: x_pair_seed::PartScratch,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:159-176
    // ```c
    // if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //     ewp->diag_table->offset = ewp->diag_table->window;
    //     s_BlastDiagClear(ewp->diag_table);
    // } else {
    //     ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    // }
    // ```
    diag_table_offset: i32,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
    // ```c
    // Int4 GetOffsetArraySize(LookupTableWrap* lookup)
    // {
    //    Int4 offset_array_size;
    //    switch (lookup->lut_type) {
    //    case eMBLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastMBLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    ...
    //    default:
    //       offset_array_size = OFFSET_ARRAY_SIZE;
    //       break;
    //    }
    //    return offset_array_size;
    // }
    // ```
    offset_pairs: Vec<OffsetPair>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:146-170
    // ```c
    // Blast_IntervalTreeInit(Int4 q_start, Int4 q_end,
    //                        Int4 s_start, Int4 s_end)
    // ...
    //     tree->s_min = s_start;
    //     tree->s_max = s_end;
    // ...
    //     s_IntervalRootNodeInit(tree, q_start, q_end, &retval);
    // ```
    // The interval tree of the subject's preliminary gapped stage and of its traceback, kept with
    // this scratch and given each subject's ranges with `reset_with_bounds` (the same state as a
    // new tree) instead of being made again for every subject. Reuse only.
    x_tree: BlastIntervalTree,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
// ```c
// hsp_list = Blast_HSPListFree(hsp_list);
// BlastInitHitListReset(init_hitlist);
// ```
impl SubjectScratch {
    fn new(query_count: usize, offset_array_size: usize, window_size: usize) -> Self {
        Self {
            hits: Vec::new(),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4076
            // ```c
            // status = BLAST_GreedyGappedAlignment(..., (Boolean) TRUE, FALSE, fence_hit);
            // init_hsp->offsets.qs_offsets.q_off = gap_align->greedy_query_seed_start;
            // init_hsp->offsets.qs_offsets.s_off = gap_align->greedy_subject_seed_start;
            // status = Blast_HSPInit(gap_align->query_start,
            //           gap_align->query_stop, gap_align->subject_start,
            //           gap_align->subject_stop,
            //           init_hsp->offsets.qs_offsets.q_off,
            //           init_hsp->offsets.qs_offsets.s_off, context,
            //           query_frame, subject->frame, gap_align->score,
            //           &(gap_align->edit_script), &new_hsp);
            // ```
            prelim_hits: Vec::new(),
            ungapped_hits: Vec::new(),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:356-357
            // ```c
            // gap_align->fwd_prelim_tback = GapPrelimEditBlockNew();
            // gap_align->rev_prelim_tback = GapPrelimEditBlockNew();
            // ```
            greedy_align_scratch: GreedyAlignScratch::new(),
            // No NCBI counterpart: the helper scratch list starts empty and is filled when helpers are used;
            // it does not change any value NCBI computes.
            x_ahead_scratch: Vec::new(),
            x_spec_scratch: Vec::new(),
            cutoff_scores: Vec::with_capacity(query_count),
            x_dropoff_scores: Vec::with_capacity(query_count),
            reduced_cutoff_scores: Vec::with_capacity(query_count),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:184-198
            // ```c
            // if (backup->allocated >= num_seq_ranges) return;
            // if (backup->allocated) {
            //     sfree(subject->seq_ranges);
            // }
            // backup->allocated = num_seq_ranges;
            // subject->seq_ranges = (SSeqRange *) calloc(backup->allocated,
            //                                sizeof(SSeqRange));
            // ```
            seq_ranges_scratch: Vec::new(),
            hit_level_array: Vec::new(),
            hit_len_array: Vec::new(),
            diag_hash: DiagHashTable::new(window_size as i32),
            // No NCBI counterpart: empty until a partitioned stage runs (LOSAT_X_PAIRPARFOLD).
            x_pair_fold: x_pair_seed::PartScratch::default(),
            diag_table_offset: window_size as i32,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:991-1041
            // ```c
            // Int4 offset_array_size = GetOffsetArraySize(lookup_wrap);
            // ...
            // aux_struct->offset_pairs =
            //   (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
            // ```c
            // Int4 GetOffsetArraySize(LookupTableWrap* lookup)
            // {
            //    ...
            //    offset_array_size = OFFSET_ARRAY_SIZE +
            //       ((BlastMBLookupTable*)lookup->lut)->longest_chain;
            //    ...
            // }
            // ```
            offset_pairs: vec![
                OffsetPair {
                    q_off: 0,
                    s_off: 0,
                    s_range: 0,
                };
                offset_array_size
            ],
            // No NCBI counterpart: every user sets the subject's ranges first (`reset_with_bounds`).
            x_tree: BlastIntervalTree::new(0, 1, 0, 1),
        }
    }
}

#[derive(Clone)]
struct QueryContext {
    query_idx: u32,
    frame: i32,
    query_offset: i32,
    /// The context's letter count. The letters are held once, in BLASTNA, in the query buffer.
    len: usize,
    masks: Vec<MaskedInterval>,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-512
/// ```c
///         } else if (kGreedyTraceback) {
///             BLAST_GreedyGappedAlignment(query, adjusted_subject,
/// ...
///         } else {
///           BLAST_GappedAlignmentWithTraceback(program_number, query,
/// ```
/// The result of the traceback call of NCBI, of the dynamic-programming kind or the greedy
/// kind, computed before the ordered loop reaches the HSP.
///
/// No NCBI counterpart: a traceback computed ahead of the ordered
/// loop, by the dynamic-programming kernel or by the greedy one.
#[derive(PartialEq)]
enum XSpecTraceback<D, G> {
    Dp(D),
    Greedy(G),
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:356-357
/// ```c
///    gap_align->fwd_prelim_tback = GapPrelimEditBlockNew();
///    gap_align->rev_prelim_tback = GapPrelimEditBlockNew();
/// ```
/// The working memory of one gapped alignment call (the `gap_align` of NCBI, and the greedy aligner's
/// auxiliary memory), borrowed from the caller or from a helper.
///
/// No NCBI counterpart: the scratch memory one gapped alignment may use.
struct XAheadScratch<'a> {
    gap: &'a mut GapAlignScratch,
    greedy: &'a mut GreedyAlignScratch,
}

#[path = "x_ctx_o1.rs"]
mod x_ctx_o1;

struct QueryContextIndex {
    offsets: Vec<usize>,
    min_length: usize,
    max_length: usize,
    // Direct mapping from query offset to context index for fast lookup.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:219-238 (BSearchContextInfo)
    // ```c
    // if (A->contexts[m].query_offset > n) e = m;
    // else b = m;
    // ...
    // return b;
    // ```
    direct_map: Vec<u32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243
    // ```c
    // Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
    // ...
    //     while (b < e - 1) {
    // 	m = (b + e) / 2;
    // 	if (A->contexts[m].query_offset > n)
    // 	    e = m;
    // 	else
    // 	    b = m;
    //     }
    //     return b;
    // ```
    // NCBI finds the context of a query offset by this binary search over the `query_offset` of the
    // contexts. `direct_map` is a table of its answer for every offset. The runs below hold the same
    // answer as intervals: for an offset inside the sequence of context `c` it is `c`, else 0.
    // No NCBI counterpart: what `direct_map` holds, as the list of
    // its non-zero runs `(start, end, context)` sorted by `start`, plus the
    // length it would have. Four bytes per query base become a few words.
    x_ranges: Vec<(usize, usize, u32)>,
    x_total_len: usize,
    // EXPERIMENT (LOSAT_X_CTXO1): what `direct_map` holds, in constant time from the layout of the
    // contexts (`x_ctx_o1.rs`). When it is built it answers the offsets the table covers, the runs
    // only those of its crowded blocks, and the table is not built.
    x_o1: Option<x_ctx_o1::CtxO1>,
}

/// No NCBI counterpart: reads the LOSAT_X_CTXO1 switch once (`x_ctx_o1.rs`).
fn x_ctx_o1_on() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_CTXO1"))
}

impl QueryContextIndex {
    fn new(contexts: &[QueryContext], one_strand: bool) -> Self {
        let offsets: Vec<usize> = contexts
            .iter()
            .map(|ctx| ctx.query_offset.max(0) as usize)
            .collect();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:229-243
        // ```c
        //                 case eNa_strand_plus:
        //                     s_QueryInfo_SetContext(query_info, ctx_index, length);
        //                     s_QueryInfo_SetContext(query_info, ctx_index+1, 0);
        //                     // the missing strand is present in query_info as
        //                     // zero-lenghth context
        //                     min_length = 0;
        // ```
        // With one strand searched `min_length` is 0, so `bsearch_context_info` searches all
        // contexts (the bounds of the other case assume every context has a length).
        let min_length = if one_strand {
            0
        } else {
            contexts
                .iter()
                .map(|ctx| ctx.len)
                .filter(|&len| len > 0)
                .min()
                .unwrap_or(0)
        };
        let max_length = contexts.iter().map(|ctx| ctx.len).max().unwrap_or(0);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:219-238 (BSearchContextInfo)
        // ```c
        // if (A->contexts[m].query_offset > n) e = m;
        // else b = m;
        // ...
        // return b;
        // ```
        let total_len = contexts
            .iter()
            .map(|ctx| ctx.query_offset.max(0) as usize + ctx.len)
            .max()
            .unwrap_or(0);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:225-227,235-241
        // ```c
        //     if (A->min_length > 0 && A->max_length > 0 && A->first_context == 0) {
        //         b = MIN(n / (A->max_length + 1), size - 1);
        //         e = MIN(n / (A->min_length + 1) + 1, size);
        // ...
        //     while (b < e - 1) {
        // ```
        // Builds the runs (start, end, context) from the same context offsets and lengths the table is
        // filled from. If two sequences would overlap, the table's "later context wins" rule would not be
        // a list of disjoint runs, so the runs are dropped and the table is used.
        // No NCBI counterpart: the runs the loop below would fill.
        // Later contexts overwrite earlier ones in `direct_map`, so the list
        // only stands in for it when no two runs overlap; context 0 needs no
        // run because the map starts out as zeros.
        // `x_mode` 1: the runs answer and the table is not built; 0: the table answers (overlapping
        // runs).
        let mut x_mode = 1u8;
        // EXPERIMENT (LOSAT_X_CTXO1): the constant-time form, where the layout allows it. It takes
        // the place of the runs (kept only for the offsets of its crowded blocks) and of the table.
        let x_o1 = if x_ctx_o1_on() {
            x_ctx_o1::CtxO1::build(contexts, total_len)
        } else {
            None
        };
        // the runs answer what the form leaves (offsets in a crowded block of `x_ctx_o1.rs`)
        let x_o1_alone = x_o1.as_ref().is_some_and(|o1| !o1.needs_runs());
        let mut x_ranges: Vec<(usize, usize, u32)> = Vec::new();
        if x_mode != 0 && !x_o1_alone {
            for (idx, ctx) in contexts.iter().enumerate() {
                let start = ctx.query_offset.max(0) as usize;
                let end = start.saturating_add(ctx.len);
                if start < end && end <= total_len && idx != 0 {
                    x_ranges.push((start, end, idx as u32));
                }
            }
            let mut all: Vec<(usize, usize)> = contexts
                .iter()
                .map(|ctx| {
                    let start = ctx.query_offset.max(0) as usize;
                    (start, start.saturating_add(ctx.len))
                })
                .filter(|&(start, end)| start < end && end <= total_len)
                .collect();
            all.sort_unstable();
            if all.windows(2).any(|w| w[0].1 > w[1].0) {
                x_mode = 0;
                x_ranges.clear();
            }
            x_ranges.sort_unstable();
        }
        let mut direct_map = Vec::new();
        if x_mode != 1 {
            direct_map = vec![0u32; total_len];
            for (idx, ctx) in contexts.iter().enumerate() {
                let start = ctx.query_offset.max(0) as usize;
                let end = start.saturating_add(ctx.len);
                if start < end && end <= direct_map.len() {
                    direct_map[start..end].fill(idx as u32);
                }
            }
        }
        Self {
            offsets,
            min_length,
            max_length,
            direct_map,
            x_ranges,
            x_total_len: if x_mode == 0 || x_o1_alone {
                0
            } else {
                total_len
            },
            x_o1,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:235-242
    /// ```c
    ///     while (b < e - 1) {
    /// 	m = (b + e) / 2;
    /// 	if (A->contexts[m].query_offset > n)
    /// 	    e = m;
    /// 	else
    /// 	    b = m;
    ///     }
    ///     return b;
    /// ```
    /// The context of offset `n` from the runs: the last run that starts at or before `n`, if `n` is
    /// inside it. This is the value the table holds at `n`.
    ///
    /// `direct_map[n]` for `n < x_total_len`, from the runs.
    #[inline]
    fn x_direct(&self, n: usize) -> usize {
        // the last run that starts at or before `n`
        let after = self.x_ranges.partition_point(|&(start, _, _)| start <= n);
        if after == 0 {
            return 0;
        }
        let (_, end, context) = self.x_ranges[after - 1];
        if n < end {
            context as usize
        } else {
            0
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219
    /// ```c
    /// Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
    /// ```
    /// The context of query offset `n`. EXPERIMENT (LOSAT_X_CTXO1): the constant-time form answers
    /// here, inlined into the seed loops; any other form is the out-of-line `x_context_for_offset_ref`.
    #[inline(always)]
    fn context_for_offset(&self, n: usize) -> usize {
        if let Some(o1) = &self.x_o1 {
            if let Some(context) = o1.get(n) {
                return context;
            }
        }
        self.x_context_for_offset_ref(n)
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:219-238 (BSearchContextInfo)
    // ```c
    // if (A->contexts[m].query_offset > n) e = m;
    // else b = m;
    // ...
    // return b;
    // ```
    // The forms besides LOSAT_X_CTXO1, out of line as `context_for_offset` was.
    #[inline(never)]
    fn x_context_for_offset_ref(&self, n: usize) -> usize {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243
        // ```c
        // Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
        // ...
        //     while (b < e - 1) {
        // 	m = (b + e) / 2;
        // 	if (A->contexts[m].query_offset > n)
        // 	    e = m;
        // 	else
        // 	    b = m;
        //     }
        //     return b;
        // ```
        // The code after this block (the table, then `bsearch_context_info`) ports `BSearchContextInfo`.
        // With the runs the same context is found without the table.
        if n < self.x_total_len {
            return self.x_direct(n);
        }
        if n < self.direct_map.len() {
            return self.direct_map[n] as usize;
        }
        bsearch_context_info(n, self)
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:219-238
// ```c
// Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
// {
//     Int4 m=0, b=0, e=0, size=0;
//     size = A->last_context+1;
//     if (A->min_length > 0 && A->max_length > 0 && A->first_context == 0) {
//         b = MIN(n / (A->max_length + 1), size - 1);
//         e = MIN(n / (A->min_length + 1) + 1, size);
//     } else { b = 0; e = size; }
//     while (b < e - 1) {
//         m = (b + e) / 2;
//         if (A->contexts[m].query_offset > n) e = m;
//         else b = m;
//     }
//     return b;
// }
// ```
fn bsearch_context_info(n: usize, index: &QueryContextIndex) -> usize {
    let size = index.offsets.len();
    if size == 0 {
        return 0;
    }
    let mut b = 0usize;
    let mut e = size;
    if index.min_length > 0 && index.max_length > 0 {
        b = (n / (index.max_length + 1)).min(size - 1);
        e = (n / (index.min_length + 1) + 1).min(size);
    }
    while b + 1 < e {
        let m = (b + e) / 2;
        if index.offsets[m] > n {
            e = m;
        } else {
            b = m;
        }
    }
    b
}

/// The query length that picks the diagonal container (hash above 8000 letters, array
/// otherwise): NCBI's `query_offset + query_length` of the LAST context. With `-strand plus`
/// the last context is the empty minus context of the last query, one past the end of the plus
/// context, so this is the searched block's length + 1; when the last context has letters it
/// is the searched block's length (`query_concat_length`).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:225-229
/// ```c
/// (query_info->contexts[query_info->last_context].query_offset +
///  query_info->contexts[query_info->last_context].query_length) > kQueryLenForHashTable)
///     p->container_type = eDiagHash;
/// ```
fn diag_container_query_length(query_contexts: &[QueryContext]) -> usize {
    query_contexts
        .last()
        .map_or(0, |last| last.query_offset.max(0) as usize + last.len)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:500-603
// ```c
// int buflen = QueryInfo_GetSeqBufLen(qinfo);
// TAutoUint1Ptr buf((Uint1*) calloc(buflen+1, sizeof(Uint1)));
// ...
// sequence = queries.GetBlastSequence(index, encoding, strand, eSentinels);
// ...
// memcpy(&buf.get()[offset], sequence.data.get(), sequence.length);
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:120-121
// ```c
// const Uint1 kProtSentinel = NULLB;
// const Uint1 kNuclSentinel = 0xF;
// ```
/// The query buffer with its sentinels (`query_concat_length + 2` bytes), every context
/// written straight into it in BLASTNA: for query `q` with context offsets
/// `context_offsets[2 * q]` (plus) and `context_offsets[2 * q + 1]` (minus), the plus
/// context is the BLASTNA code of the letters and the minus context that of their reverse
/// complement. The searched block (`query_concat_length` bytes) is `[1..n + 1]` of it
/// (LOSAT's memory: no second copy of the block, no IUPAC copy of either strand).
fn build_query_blastna_buffer<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
) -> Vec<u8> {
    debug_assert_eq!(2 * query_letters.count(), context_offsets.len());

    let mut with_sentinels = vec![NUCL_SENTINEL; query_concat_length + 2];

    for (q, offsets) in context_offsets.chunks_exact(2).enumerate() {
        query_letters.read(q, |letters| {
            let length = letters.len();
            // A context that is empty or does not end inside the block stays sentinels.
            // (Twin: `place` and the two encoder calls below are copied in x_query_prep_par.rs
            // (LOSAT_X_QPREPPAR) and, with the strand switches, in
            // `build_strand_query_blastna_buffer`; change them there too.)
            let place = |offset: i32| -> Option<std::ops::Range<usize>> {
                let start = offset.max(0) as usize;
                let end = start.saturating_add(length);
                (start < end && end <= query_concat_length).then(|| start + 1..end + 1)
            };
            if let Some(range) = place(offsets[0]) {
                encode_iupac_to_blastna_into(letters, &mut with_sentinels[range]);
            }
            if let Some(range) = place(offsets[1]) {
                encode_reverse_complement_blastna_into(letters, &mut with_sentinels[range]);
            }
        });
    }

    with_sentinels
}

/// The query buffer of a search of one strand (`-strand plus` or `minus`): as
/// `build_query_blastna_buffer`, but a context is written only when its length is above 0,
/// with its own length (the strand that is not searched has a context of length 0 and takes no
/// room, so it must not be written over the next query).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:71-81
/// ```c
///         Uint4 shift = prev_len ? prev_len + 1 : 0;
///
///         qinfo->contexts[index].query_offset = prev_loc + shift;
/// ```
fn build_strand_query_blastna_buffer<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    strand: QueryStrand,
) -> Vec<u8> {
    debug_assert_eq!(2 * query_letters.count(), context_offsets.len());
    let (plus_on, minus_on) = strand.searched();

    let mut with_sentinels = vec![NUCL_SENTINEL; query_concat_length + 2];

    for (q, offsets) in context_offsets.chunks_exact(2).enumerate() {
        query_letters.read(q, |letters| {
            let place = |offset: i32, length: usize| -> Option<std::ops::Range<usize>> {
                let start = offset.max(0) as usize;
                let end = start.saturating_add(length);
                (start < end && end <= query_concat_length).then(|| start + 1..end + 1)
            };
            if plus_on {
                if let Some(range) = place(offsets[0], letters.len()) {
                    encode_iupac_to_blastna_into(letters, &mut with_sentinels[range]);
                }
            }
            if minus_on {
                if let Some(range) = place(offsets[1], letters.len()) {
                    encode_reverse_complement_blastna_into(letters, &mut with_sentinels[range]);
                }
            }
        });
    }

    with_sentinels
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_util.c:190-203
// ```c
// Int4 num_entries = 0;
// Int4 curr_max = 0;
// ...
// num_entries += loc->ssr->right - loc->ssr->left;
// curr_max = MAX(curr_max, loc->ssr->right);
// ...
// *max_off = curr_max;
// return num_entries;
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:402-406
// ```c
// BlastLookupIndexQueryExactMatches(thin_backbone,
//                                   lookup->word_length,
//                                   BITS_PER_NUC,
//                                   lookup->lut_word_length,
//                                   query, locations);
// ```
fn compute_lookup_query_stats(
    contexts: &[QueryContext],
    context_masks: &[Vec<MaskedInterval>],
) -> (usize, usize) {
    debug_assert_eq!(contexts.len(), context_masks.len());
    let mut approx_entries = 0i64;
    let mut max_q_off = 0usize;

    // NCBI computes the lookup segments (BLAST_ComplementMaskLocations, which skips the
    // contexts marked invalid) before its ungapped blocks mark contexts invalid
    // (blast_setup.c:633-653), so every context of a query with letters counts.
    for (ctx, masks) in contexts.iter().zip(context_masks.iter()) {
        if ctx.len == 0 {
            continue;
        }
        let ctx_offset = ctx.query_offset.max(0) as usize;
        let ranges = build_unmasked_ranges(ctx.len, masks);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1081-1115
        // ```c
        //          if (first) {
        //             last_interval_open = TRUE;
        //             first = FALSE;
        //
        //             if (filter_start > start_offset) {
        //                /* beginning of sequence not filtered */
        //                left = start_offset;
        //             } else {
        //                /* beginning of sequence filtered */
        //                left = filter_end + 1;
        //                continue;
        //             }
        //          }
        // ...
        //       if (last_interval_open) {
        //          /* Need to finish SSeqRange* for last interval. */
        //          right = end_offset;
        // ```
        // A context masked from its first letter to its last gets the range
        // (end_offset + 1, end_offset): one entry less, and its end in `max_off`.
        if ranges.is_empty() {
            approx_entries -= 1;
            max_q_off = max_q_off.max(ctx_offset + ctx.len - 1);
            continue;
        }
        for (start, end) in ranges {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_util.c:190-203
            // ```c
            // num_entries += loc->ssr->right - loc->ssr->left;
            // curr_max = MAX(curr_max, loc->ssr->right);
            // ```
            // NCBI's `right` is inclusive (the last unmasked letter), LOSAT's `end` is not.
            approx_entries += (end - start) as i64 - 1;
            max_q_off = max_q_off.max(ctx_offset + end - 1);
        }
    }

    // A negative estimate (every context masked) selects the smallest table, as 0 does.
    (approx_entries.max(0) as usize, max_q_off)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:739-776
// ```c
// SSeqRange ranges2scan[] = { {0, 501}, {700, 1001} , {subject_bases, subject_bases}};
// ...
// if ( s_off >= (Uint4)ranges2scan[j].left &&
//      s_off <  (Uint4)ranges2scan[j].right ) {
//     hit_found = TRUE;
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1647-1674
// ```c
// scan_range[1] = subject->seq_ranges[0].left + word_length - lut_word_length;
// scan_range[2] = subject->seq_ranges[0].right - lut_word_length;
// while (s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
//     hitsfound = scansub(..., &scan_range[1]);
// }
// ```
// The ranges of a soft-masked subject (bl2seq `-lcase_masking`) come from the inclusive
// (from, to) pairs of its lowercase masks.
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:700-707
// ```c
//         ITERATE(CPacked_seqint::Tdata, itr, slp->GetPacked_int().Get()) {
//     	    CSeqDB::TOffsetPair p;
//             p.first = ((*itr)->GetFrom() > offset)? (*itr)->GetFrom() - offset : 0;
//             p.second = MIN((*itr)->GetTo() - offset, length-1);
//
//             if ((*itr)->GetTo() >= offset && p.first < length) {
//                 output.push_back(p);
//             }
// ```
// The pairs are stored from `_data[1]`, but `get_data` returns `_data`, which
// `SetupSubjects_OMF` passes as `size + 1` ranges.
// NCBI reference: ncbi-blast/c++/include/objtools/blast/seqdb_reader/seqdb.hpp:280-282
// ```c
//         value_type& operator[](size_type i) { return (value_type &)_data[1+ 2*i]; }
//
//         value_type * get_data() const { return (value_type *) _data; }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:819-826
// ```c
//             s_SeqLoc2MaskedSubjRanges(masks, &*range, length,  masked_ranges);
//             if ( !masked_ranges.empty() ) {
//                 ...
//                 BlastSeqBlkSetSeqRanges(subj, (SSeqRange*) masked_ranges.get_data(),
//                                     static_cast<Uint4>(masked_ranges.size()) + 1, true, eSoftSubjMasking);
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:211-213
// ```c
//     // Fill out the boundary of the sequence to compliment the masks
//     tmp[0].left = 0;
//     tmp[num_seq_ranges - 1].right = seq_blk->length;
// ```
// So range k is {`to` of mask k-1, `from` of mask k}: it starts at the last lowercase
// letter of the previous mask (a word may start there) and ends before the next mask.
fn build_subject_seq_ranges_from_masks(
    masks: &[MaskedInterval],
    subject_len: usize,
) -> Vec<(i32, i32)> {
    let mut ranges = Vec::with_capacity(masks.len() + 1);
    let mut left = 0i32;
    for mask in masks.iter().filter(|mask| mask.start < subject_len) {
        ranges.push((left, mask.start as i32));
        left = (mask.end.min(subject_len) - 1) as i32;
    }
    ranges.push((left, subject_len as i32));
    ranges
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2049-2067
// ```c
// Int2 Blast_TrimHSPListByMaxHsps(BlastHSPList* hsp_list,
//                                const BlastHitSavingOptions* hit_options)
// {
//    if ((hsp_list == NULL) ||
//        (hit_options->max_hsps_per_subject == 0) ||
//        (hsp_list->hspcnt <= hit_options->max_hsps_per_subject))
//       return 0;
//    ...
//    hsp_list->hspcnt = hsp_max;
//    return 0;
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:877-892
// ```c
// for (query_index = 0; query_index < results->num_queries; ++query_index) {
//    if (!(hit_list = results->hitlist_array[query_index])) continue;
//    for (subject_index = hitlist_size;
//         subject_index < hit_list->hsplist_count; ++subject_index) {
//       hit_list->hsplist_array[subject_index] =
//       Blast_HSPListFree(hit_list->hsplist_array[subject_index]);
//    }
//    hit_list->hsplist_count = MIN(hit_list->hsplist_count, hitlist_size);
// }
// ```
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
// ```c
// typedef struct BlastHSPList {
//    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
//    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
//                       Set to 0 if not applicable */
//    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
//    Int4 hspcnt; /**< Number of HSPs saved */
//    ...
// } BlastHSPList;
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:493-555
// ```c
// if (aux_struct->WordFinder) {
//     aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
//                            word_params, aux_struct->ewp,
//                            aux_struct->offset_pairs,
//                            kScanSubjectOffsetArraySize,
//                            init_hitlist, ungapped_stats);
//     if (init_hitlist->total == 0) continue;
// }
// ...
// if (aux_struct->GetGappedScore) {
//     status = aux_struct->GetGappedScore(program_number, query,
//             query_info, subject, gap_align, score_params, ext_params,
//             hit_params, word_params, init_hitlist, &hsp_list,
//             gapped_stats, NULL);
// }
// ...
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
// Blast_HSPListSortByScore(hsp_list);
// ```
struct BlastnTiming {
    scan_ns: std::sync::atomic::AtomicU64,
    scan_calls: std::sync::atomic::AtomicU64,
    ungapped_ns: std::sync::atomic::AtomicU64,
    ungapped_calls: std::sync::atomic::AtomicU64,
    gapped_ns: std::sync::atomic::AtomicU64,
    gapped_calls: std::sync::atomic::AtomicU64,
    traceback_ns: std::sync::atomic::AtomicU64,
    traceback_calls: std::sync::atomic::AtomicU64,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:358-613
    // ```c
    // ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
    // tree = Blast_IntervalTreeInit(...);
    // for (index=0; index < num_initial_hsps; index++) {
    //    if (... !BlastIntervalTreeContainsHSP(tree, hsp, query_info, ...)) {
    //       BlastGetOffsetsForGappedAlignment(...);
    //       BLAST_GreedyGappedAlignment(...);
    //       Blast_HSPUpdateWithTraceback(gap_align, hsp);
    //       BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
    //    }
    // }
    // ```
    traceback_sort_prelim_ns: std::sync::atomic::AtomicU64,
    traceback_tree_precheck_ns: std::sync::atomic::AtomicU64,
    traceback_start_offsets_ns: std::sync::atomic::AtomicU64,
    traceback_alignment_ns: std::sync::atomic::AtomicU64,
    traceback_alignment_dp_ns: std::sync::atomic::AtomicU64,
    traceback_alignment_greedy_ns: std::sync::atomic::AtomicU64,
    traceback_hsp_build_ns: std::sync::atomic::AtomicU64,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692
    // ```c
    // extra_start = Blast_HSPListPurgeHSPsWithCommonEndpoints(..., FALSE);
    // delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(...);
    // Blast_HSPListPurgeHSPsWithCommonEndpoints(..., TRUE);
    // Blast_HSPListSortByScore(hsp_list);
    // Blast_IntervalTreeReset(tree);
    // for (index = 0; index < hsp_list->hspcnt; index++) {
    //    if (BlastIntervalTreeContainsHSP(...)) ...
    //    else BlastIntervalTreeAddHSP(...);
    // }
    // ```
    purge_common_endpoint_pass1_ns: std::sync::atomic::AtomicU64,
    reevaluate_ambiguities_ns: std::sync::atomic::AtomicU64,
    identity_length_test_ns: std::sync::atomic::AtomicU64,
    purge_common_endpoint_pass2_ns: std::sync::atomic::AtomicU64,
    traceback_score_resort_ns: std::sync::atomic::AtomicU64,
    traceback_final_tree_contains_ns: std::sync::atomic::AtomicU64,
    traceback_final_tree_add_ns: std::sync::atomic::AtomicU64,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:836-892
    // ```c
    // Blast_TrimHSPListByMaxHsps(hsp_list, hit_options);
    // s_BlastPruneExtraHits(results, hitlist_size);
    // ```
    traceback_postprocess_ns: std::sync::atomic::AtomicU64,
    traceback_prelim_hsps: std::sync::atomic::AtomicU64,
    traceback_tree_precheck_skipped_hsps: std::sync::atomic::AtomicU64,
    traceback_full_traceback_hsps: std::sync::atomic::AtomicU64,
    traceback_deleted_evalue_cutoff_hsps: std::sync::atomic::AtomicU64,
    traceback_deleted_reevaluation_hsps: std::sync::atomic::AtomicU64,
    traceback_deleted_identity_length_hsps: std::sync::atomic::AtomicU64,
    traceback_removed_endpoint_pass1_hsps: std::sync::atomic::AtomicU64,
    traceback_removed_endpoint_pass2_hsps: std::sync::atomic::AtomicU64,
    traceback_removed_final_tree_hsps: std::sync::atomic::AtomicU64,
    traceback_edit_script_lengths: std::sync::Mutex<Vec<u32>>,
    traceback_alignment_lengths: std::sync::Mutex<Vec<u32>>,
    format_ns: std::sync::atomic::AtomicU64,
    format_calls: std::sync::atomic::AtomicU64,
}

impl BlastnTiming {
    fn new() -> Self {
        Self {
            scan_ns: std::sync::atomic::AtomicU64::new(0),
            scan_calls: std::sync::atomic::AtomicU64::new(0),
            ungapped_ns: std::sync::atomic::AtomicU64::new(0),
            ungapped_calls: std::sync::atomic::AtomicU64::new(0),
            gapped_ns: std::sync::atomic::AtomicU64::new(0),
            gapped_calls: std::sync::atomic::AtomicU64::new(0),
            traceback_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_calls: std::sync::atomic::AtomicU64::new(0),
            traceback_sort_prelim_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_tree_precheck_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_start_offsets_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_alignment_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_alignment_dp_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_alignment_greedy_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_hsp_build_ns: std::sync::atomic::AtomicU64::new(0),
            purge_common_endpoint_pass1_ns: std::sync::atomic::AtomicU64::new(0),
            reevaluate_ambiguities_ns: std::sync::atomic::AtomicU64::new(0),
            identity_length_test_ns: std::sync::atomic::AtomicU64::new(0),
            purge_common_endpoint_pass2_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_score_resort_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_final_tree_contains_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_final_tree_add_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_postprocess_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_prelim_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_tree_precheck_skipped_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_full_traceback_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_deleted_evalue_cutoff_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_deleted_reevaluation_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_deleted_identity_length_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_removed_endpoint_pass1_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_removed_endpoint_pass2_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_removed_final_tree_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_edit_script_lengths: std::sync::Mutex::new(Vec::new()),
            traceback_alignment_lengths: std::sync::Mutex::new(Vec::new()),
            format_ns: std::sync::atomic::AtomicU64::new(0),
            format_calls: std::sync::atomic::AtomicU64::new(0),
        }
    }

    fn record_duration(counter: &std::sync::atomic::AtomicU64, start: std::time::Instant) {
        counter.fetch_add(
            start.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
    }

    fn record_ns(counter: &std::sync::atomic::AtomicU64, elapsed_ns: u64) {
        counter.fetch_add(elapsed_ns, std::sync::atomic::Ordering::Relaxed);
    }

    fn record_count(counter: &std::sync::atomic::AtomicU64, count: u64) {
        counter.fetch_add(count, std::sync::atomic::Ordering::Relaxed);
    }

    fn record_traceback_lengths(&self, edit_script_lengths: Vec<u32>, alignment_lengths: Vec<u32>) {
        if !edit_script_lengths.is_empty() {
            if let Ok(mut samples) = self.traceback_edit_script_lengths.lock() {
                samples.extend(edit_script_lengths);
            }
        }
        if !alignment_lengths.is_empty() {
            if let Ok(mut samples) = self.traceback_alignment_lengths.lock() {
                samples.extend(alignment_lengths);
            }
        }
    }
}

fn timing_seconds(counter: &std::sync::atomic::AtomicU64) -> f64 {
    counter.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9
}

fn timing_count(counter: &std::sync::atomic::AtomicU64) -> u64 {
    counter.load(std::sync::atomic::Ordering::Relaxed)
}

fn timing_length_stats(samples: &std::sync::Mutex<Vec<u32>>) -> (usize, u32, f64) {
    let Ok(samples) = samples.lock() else {
        return (0, 0, 0.0);
    };
    if samples.is_empty() {
        return (0, 0, 0.0);
    }
    let mut sorted = samples.clone();
    sorted.sort_unstable();
    let count = sorted.len();
    let max = sorted[count - 1];
    let median = if count % 2 == 0 {
        let upper = sorted[count / 2] as f64;
        let lower = sorted[(count / 2) - 1] as f64;
        (lower + upper) / 2.0
    } else {
        sorted[count / 2] as f64
    };
    (count, max, median)
}

// No NCBI counterpart: LOSAT_TIMING lines of one query batch (diagnostics only). The batch
// returns them as text, which `run_in_pool` writes to standard error when it takes the batch.
fn blastn_timing_text(
    timing: &BlastnTiming,
    t_search_start: Option<std::time::Instant>,
    t_total: Option<std::time::Instant>,
) -> String {
    use std::fmt::Write as _;
    let mut text = String::new();
    // One line of the block, as `eprintln!` wrote it before.
    macro_rules! timing_line {
        ($($arg:tt)*) => {{
            let _ = writeln!(text, $($arg)*);
        }};
    }
    let scan_s = timing_seconds(&timing.scan_ns);
    let scan_n = timing_count(&timing.scan_calls);
    let ungapped_s = timing_seconds(&timing.ungapped_ns);
    let ungapped_n = timing_count(&timing.ungapped_calls);
    let gapped_s = timing_seconds(&timing.gapped_ns);
    let gapped_n = timing_count(&timing.gapped_calls);
    let traceback_s = timing_seconds(&timing.traceback_ns);
    let traceback_n = timing_count(&timing.traceback_calls);
    let format_s = timing_seconds(&timing.format_ns);
    let format_n = timing_count(&timing.format_calls);

    timing_line!("[TIMING] scan_lookup: {:.3}s (calls={})", scan_s, scan_n);
    timing_line!(
        "[TIMING] ungapped_extend: {:.3}s (calls={})",
        ungapped_s,
        ungapped_n
    );
    timing_line!(
        "[TIMING] gapped_extend: {:.3}s (calls={})",
        gapped_s,
        gapped_n
    );
    timing_line!(
        "[TIMING] traceback_prune: {:.3}s (calls={})",
        traceback_s,
        traceback_n
    );
    timing_line!(
        "[TIMING] traceback_sort_prelim: {:.3}s",
        timing_seconds(&timing.traceback_sort_prelim_ns)
    );
    timing_line!(
        "[TIMING] traceback_tree_precheck: {:.3}s",
        timing_seconds(&timing.traceback_tree_precheck_ns)
    );
    timing_line!(
        "[TIMING] traceback_start_offsets: {:.3}s",
        timing_seconds(&timing.traceback_start_offsets_ns)
    );
    timing_line!(
        "[TIMING] traceback_alignment: {:.3}s (dp={:.3}s greedy={:.3}s)",
        timing_seconds(&timing.traceback_alignment_ns),
        timing_seconds(&timing.traceback_alignment_dp_ns),
        timing_seconds(&timing.traceback_alignment_greedy_ns)
    );
    timing_line!(
        "[TIMING] traceback_hsp_build: {:.3}s",
        timing_seconds(&timing.traceback_hsp_build_ns)
    );
    timing_line!(
        "[TIMING] purge_common_endpoint_pass1: {:.3}s",
        timing_seconds(&timing.purge_common_endpoint_pass1_ns)
    );
    timing_line!(
        "[TIMING] reevaluate_ambiguities: {:.3}s",
        timing_seconds(&timing.reevaluate_ambiguities_ns)
    );
    timing_line!(
        "[TIMING] identity_length_test: {:.3}s",
        timing_seconds(&timing.identity_length_test_ns)
    );
    timing_line!(
        "[TIMING] purge_common_endpoint_pass2: {:.3}s",
        timing_seconds(&timing.purge_common_endpoint_pass2_ns)
    );
    timing_line!(
        "[TIMING] traceback_score_resort: {:.3}s",
        timing_seconds(&timing.traceback_score_resort_ns)
    );
    timing_line!(
        "[TIMING] traceback_final_tree_contains: {:.3}s",
        timing_seconds(&timing.traceback_final_tree_contains_ns)
    );
    timing_line!(
        "[TIMING] traceback_final_tree_add: {:.3}s",
        timing_seconds(&timing.traceback_final_tree_add_ns)
    );
    timing_line!(
        "[TIMING] traceback_postprocess: {:.3}s",
        timing_seconds(&timing.traceback_postprocess_ns)
    );
    timing_line!(
        "[TIMING] traceback_counts: prelim={} tree_skipped={} full_traceback={} deleted_evalue_cutoff={} deleted_reevaluation={} deleted_identity_length={} endpoint_pass1_removed={} endpoint_pass2_removed={} final_tree_removed={}",
        timing_count(&timing.traceback_prelim_hsps),
        timing_count(&timing.traceback_tree_precheck_skipped_hsps),
        timing_count(&timing.traceback_full_traceback_hsps),
        timing_count(&timing.traceback_deleted_evalue_cutoff_hsps),
        timing_count(&timing.traceback_deleted_reevaluation_hsps),
        timing_count(&timing.traceback_deleted_identity_length_hsps),
        timing_count(&timing.traceback_removed_endpoint_pass1_hsps),
        timing_count(&timing.traceback_removed_endpoint_pass2_hsps),
        timing_count(&timing.traceback_removed_final_tree_hsps)
    );
    let (edit_count, edit_max, edit_median) =
        timing_length_stats(&timing.traceback_edit_script_lengths);
    let (aln_count, aln_max, aln_median) = timing_length_stats(&timing.traceback_alignment_lengths);
    timing_line!(
        "[TIMING] traceback_lengths: edit_script_samples={} edit_script_max={} edit_script_median={:.1} alignment_samples={} alignment_max={} alignment_median={:.1}",
        edit_count, edit_max, edit_median, aln_count, aln_max, aln_median
    );
    timing_line!(
        "[TIMING] format_output: {:.3}s (calls={})",
        format_s,
        format_n
    );
    if let Some(t_search_start) = t_search_start {
        let search_total = t_search_start.elapsed();
        // Diagnostic (LOSAT_TIMING): the same interval feeds the common line (utils/xtiming.rs):
        // one query batch's search, from after its lookup table to its end (the batches' sum
        // over all threads). No NCBI counterpart.
        crate::utils::xtiming::add_duration(crate::utils::xtiming::Phase::Search, search_total);
        timing_line!("[TIMING] search_total: {:.3}s", search_total.as_secs_f64());
    }
    if let Some(t_total) = t_total {
        timing_line!("[TIMING] total: {:.3}s", t_total.elapsed().as_secs_f64());
    }
    text
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_collector.c:121-161
// ```c
// if (!(tmp_hsp_list = hsp_list_array[query_index])) {
//    hsp_list_array[query_index] = tmp_hsp_list =
//       Blast_HSPListNew(params->hsp_num_max);
//    tmp_hsp_list->oid = hsp_list->oid;
// }
// Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
// ...
// if (!results->hitlist_array[index]) {
//    results->hitlist_array[index] =
//       Blast_HitListNew(params->prelim_hitlist_size);
// }
// Blast_HitListUpdate(results->hitlist_array[index],
//                     hsp_list_array[index]);
// ```
fn update_hitlists_with_subject_hits(
    hit_lists: &mut [Option<BlastnHitList>],
    hits: Vec<BlastnHsp>,
    prelim_hitlist_size: usize,
) {
    if hits.is_empty() {
        return;
    }

    // Each query's HSPs get a vector of their exact count (LOSAT's memory: a vector grown by
    // `push` holds four HSPs from the first one, and most subjects give a query one HSP).
    let mut counts = vec![0usize; hit_lists.len()];
    for hit in &hits {
        if let Some(count) = counts.get_mut(hit.q_idx as usize) {
            *count += 1;
        }
    }
    let mut per_query: Vec<Vec<BlastnHsp>> = counts.into_iter().map(Vec::with_capacity).collect();
    for hit in hits {
        let q_idx = hit.q_idx as usize;
        if q_idx < per_query.len() {
            per_query[q_idx].push(hit);
        }
    }

    for (q_idx, mut hsps) in per_query.into_iter().enumerate() {
        if hsps.is_empty() {
            continue;
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:555
        // ```c
        // Blast_HSPListSortByScore(hsp_list);
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1374-1382
        // ```c
        // void Blast_HSPListSortByScore(BlastHSPList* hsp_list)
        // {
        //     if (hsp_list->hspcnt > 1) {
        //         qsort(hsp_list->hsp_array, hsp_list->hspcnt,
        //               sizeof(BlastHSP*), ScoreCompareHSPs);
        //     }
        // }
        // ```
        sort_hsps_by_score(&mut hsps);

        let oid = hsps[0].s_idx;
        let hsp_list = BlastnHspList {
            oid,
            query_index: q_idx as u32,
            hsps,
            best_evalue: i32::MAX as f64,
        };

        let hit_list =
            hit_lists[q_idx].get_or_insert_with(|| BlastnHitList::new(prelim_hitlist_size));
        hit_list.update(hsp_list);
    }
}

/// The SAM document that a report call adds its queries to (outfmt 17): the stream's, which
/// the last call writes, or the document of a call of its own (the off path), which the call
/// writes at its end. NCBI writes the document when its formatter is destroyed.
struct SamTarget<'a> {
    document: &'a mut crate::report::sam::SamDocument,
    /// Whether the call writes the document after adding its queries.
    finish: bool,
}

/// What the pairwise report, the hit records and the query warnings need besides the
/// final hit list.
struct BlastnReportInputs<'a> {
    queries: &'a [FastaRecord],
    /// The subject records and their letters (theirs, or with LOSAT_X_BNSTORE their packed
    /// subjects').
    subjects: SubjectLetters<'a>,
    /// Where the searched letters lie in the records (`-query_loc`, `-subject_loc`).
    ranges: &'a ReportRanges<'a>,
    /// Per query, plus strand: DUST and, with `-lcase_masking`, the input lowercase.
    query_masks: &'a [Vec<MaskedInterval>],
    /// The DUST-only masks of every query (outfmt 11's `LCaseMask`), or with
    /// `-window_masker_db` the WindowMasker masks; empty without outfmt 11.
    archive_dust: &'a [Vec<MaskedInterval>],
    lcase_masking: bool,
    /// Per query context (the plus strand of query `q` is context `2 * q`).
    query_eff_searchsp: &'a [i64],
    /// The first valid context of a query, as its offset in the query's two contexts: 1 for
    /// `-strand minus`, else 0 (`context 2 * q + first_context`).
    first_context: usize,
    megablast: bool,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    /// The task is megablast or blastn (`BlastnPairwiseReport::zero_gap_extension_formula`).
    zero_gap_extension_formula: bool,
    /// The two-hit window of the task, printed in the epilog when not 0.
    window_size: usize,
    /// An ungapped (`-ungapped`) search: one Karlin block, no gap costs in the reports, and the
    /// hit list of the end of the traceback is the preliminary one (`finish_query_hit_list`).
    ungapped: bool,
    /// Per query context (the plus strand of query `q` is context `2 * q`); `None` for
    /// the contexts of an invalid query (not searched by NCBI) and of the strand that
    /// `-strand` does not search.
    context_karlin: &'a [Option<ContextKarlin>],
    /// The numbers of descriptions and alignments of the pairwise report
    /// (`BlastnArgs::report_counts`).
    num_descriptions: usize,
    num_alignments: usize,
    /// `-culling_limit` (0: no culling), read when the hit lists are finished.
    culling_limit: i32,
    db_num_seqs: usize,
    db_len_total: usize,
    /// The `Database:` lines of a database search's outfmt 0 report.
    database: Option<&'a PairwiseDatabase>,
    /// The tabular subject ids of a database search's reported sequences, by rank
    /// (`database::ReportedSubjects::tabular_subject_ids`); `None` for a subject set.
    database_subject_ids: Option<&'a [Arc<[u8]>]>,
    /// The ordinal ID of each reported sequence of a database search, by rank
    /// (`database::ReportedSubjects::ordinals`); `None` for a subject set.
    database_ordinals: Option<&'a [u32]>,
    /// The opened database of a database search, which the XML2 and JSON formats read the taxids
    /// of the reported sequences from; `None` for a subject set.
    database_reader: Option<&'a database::BlastnDatabase>,
    /// The reported sequences with parsed Seq-ids of a database search, by rank
    /// (`database::ReportedSubjects::parsed`; empty without them): their ids and titles come from
    /// their definition lines.
    database_parsed: &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    /// Whether the hit lists were finished after each query batch's traceback (a database
    /// search, `run_in_pool`), not by `post_process_hits_and_write`.
    finished_hit_lists: bool,
    /// Whether NCBI's epilog (`PrintEpilog`: the outfmt 0 footer, outfmt 7's last line)
    /// follows the reports; an error in a later batch skips it.
    epilog: bool,
    /// Whether the outfmt 0 report starts with the prolog, which is otherwise written
    /// before the search (`write_pairwise_prologs`).
    prolog: bool,
    /// The queries whose reports earlier calls wrote, which outfmt 7's epilog counts too
    /// (NCBI's `m_QueriesFormatted`): with LOSAT_X_RPTSTREAM the last call writes only the
    /// epilogs, after the reports of every batch (`x_rpt_stream`); 0 otherwise.
    queries_before: usize,
    /// EXPERIMENT (LOSAT_X_RPTSTREAM): the report values of a subject set's subjects, made once
    /// for the reports of all batches (`x_rpt_stream::SubjectCache`); `None`: made here.
    subject_cache: Option<&'a x_rpt_stream::SubjectCache>,
    /// EXPERIMENT (LOSAT_X_PAIRWISEPAR): the pool that formats the outfmt 0 rows and text
    /// (`report::x_pairwise_par`); `None`: the serial writer.
    x_pairwise_pool: Option<&'a crate::utils::threading::SearchPool<'a>>,
    /// EXPERIMENT (LOSAT_X_BNTABPAR): the pool that formats the outfmt 6 / 7 rows
    /// (`algorithm::blastn::x_tab_par`); `None`: the serial writer.
    x_tab_pool: Option<&'a crate::utils::threading::SearchPool<'a>>,
    /// `-out`: the base of the side files of outfmt 13 and 14 and of the names their master
    /// file lists.
    out: Option<&'a std::path::Path>,
    /// The e-value threshold of the search (BLAST XML's `Parameters_expect`).
    expect: f64,
    /// The DUST options as NCBI's filter string prints them (BLAST XML's `Parameters_filter`);
    /// `None`: no DUST.
    dust: Option<(i32, i32, i32)>,
    /// `-window_masker_db` as given (the filter string's `W -d <name>;`).
    window_masker_db: Option<&'a str>,
    /// The database length behind the effective lengths (`-dbsize`, or the total of the
    /// subjects), for the length adjustment that BLAST XML's `Statistics_hsp-len` shows.
    eff_db_length: i64,
    /// The run's BLAST archive (outfmt 11), and the queries of each batch of this call
    /// (`BatchResults::batches`).
    archive: Option<&'a crate::report::archive::ArchiveHandle>,
    batch_spans: &'a [std::ops::Range<usize>],
    /// Whether both strands of each query are searched (`-strand both`).
    both_strands: bool,
    minus_strand: bool,
    /// Whether `taxdb.bti` is on NCBI's search path, for the custom tabular reports
    /// (`ShownSubjects::tabular_taxdb_found`).
    tabular_taxdb_found: bool,
}

impl BlastnReportInputs<'_> {
    /// The number of sequences the reports print: the database's own (with a taxid list, the
    /// whole database's: `SeqDb::printed_num_seqs`), or the subject set's.
    fn shown_num_seqs(&self) -> usize {
        self.database
            .map_or(self.db_num_seqs, |database| database.num_sequences as usize)
    }

    /// The total length the reports print (as `shown_num_seqs`).
    fn shown_total_letters(&self) -> usize {
        self.database.map_or(self.db_len_total, |database| {
            usize::try_from(database.total_letters).unwrap_or(usize::MAX)
        })
    }
}

/// The per-query and run-level data of the pairwise report.
// Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
// shared post-processing makes every BLASTN run compile it in Wasm hosts.
#[inline(never)]
fn blastn_pairwise_report(
    report: &BlastnReportInputs<'_>,
    unsearched: Vec<bool>,
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
) -> Result<(Vec<BlastnPairwiseQuery>, BlastnPairwiseReport)> {
    let context_karlin = report.context_karlin;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_results.cpp:82-104
    // ```c
    //     // find the first valid context corresponding to this query
    // ```
    // The first valid context of a query is its minus context with `-strand minus`.
    let first_context = report.first_context;
    let queries = report
        .queries
        .iter()
        .enumerate()
        .map(|(q_idx, query)| BlastnPairwiseQuery {
            query_name: query_titles[q_idx].to_vec(),
            // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:742-744
            // ```c++
            //         if(cbs.IsSetInst() && cbs.GetInst().CanGetLength()){
            //             out << "\nLength=";
            //             out << cbs.GetInst().GetLength() <<"\n";
            // ```
            // The length of the query's record, also where `-query_loc` searches part of it.
            query_length: report
                .ranges
                .queries
                .placements
                .length(q_idx, query.seq().len()),
            karlin: context_karlin[2 * q_idx + first_context]
                .map(|blocks| (blocks.ungapped, blocks.gapped)),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_results.cpp:82-104
            // ```c
            //     // find the first valid context corresponding to this query
            //     ...
            //     m_SearchSpace = ctx->eff_searchsp;
            // ```
            // An invalid query has no valid context, and its search space stays 0.
            effective_search_space: if context_karlin[2 * q_idx + first_context].is_some() {
                report.query_eff_searchsp[2 * q_idx + first_context]
            } else {
                0
            },
        })
        .collect();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2910-2928
    // ```c
    //     m_NumDescriptions = m_DfltNumDescriptions;
    //     m_NumAlignments = m_DfltNumAlignments;
    //     ...
    //     if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //         m_NumDescriptions = args[kArgMaxTargetSequences].AsInteger();
    //         m_NumAlignments = args[kArgMaxTargetSequences].AsInteger();
    // ```
    // The defaults are 500 descriptions and 250 alignments (format_flags.cpp:219,221).
    let (num_descriptions, num_alignments) = (report.num_descriptions, report.num_alignments);
    let pairwise_report = BlastnPairwiseReport {
        version: NCBI_BLASTN_VERSION.to_string(),
        megablast: report.megablast,
        database_name: subject_title.to_string(),
        database_num_sequences: report.shown_num_seqs(),
        database_total_letters: report.shown_total_letters(),
        database: report.database.cloned(),
        reward: report.reward,
        penalty: report.penalty,
        gap_open: report.gap_open,
        gap_extend: report.gap_extend,
        zero_gap_extension_formula: report.zero_gap_extension_formula,
        window_size: report.window_size,
        gapped: !report.ungapped,
        num_descriptions,
        num_alignments,
        unsearched,
        epilog: report.epilog,
        prolog: report.prolog,
        parsed_subjects: report.database_parsed.to_vec(),
    };
    Ok((queries, pairwise_report))
}

/// The values of the first query of a batch (`first`, an index of the call) that the archive
/// takes (`archive::FirstContext`): the blocks and search statistics of its first valid
/// context (`CBlastAncillaryData`, blast_results.cpp:72-115; the plus strand's, or the minus
/// strand's with `-strand minus`), and the search space of its context 0 (the option
/// `SplitQuery_SetEffectiveSearchSpace` sets).
fn archive_first_context(
    report: &BlastnReportInputs<'_>,
    query_lengths: &[usize],
    first: usize,
) -> Result<archive::FirstContext> {
    let context0_space = report.query_eff_searchsp[2 * first];
    let Some(blocks) = report.context_karlin[2 * first + report.first_context] else {
        return Ok(archive::FirstContext {
            karlin: None,
            search_space: 0,
            length_adjustment: 0,
            context0_space,
        });
    };
    let search_space = report.query_eff_searchsp[2 * first + report.first_context];
    let length_adjustment = compute_length_adjustment_ncbi(
        query_lengths[first] as i64,
        report.eff_db_length,
        report.db_num_seqs as i64,
        &blocks.gapped,
    )
    .length_adjustment;
    Ok(archive::FirstContext {
        karlin: Some(blocks),
        search_space,
        length_adjustment,
        context0_space,
    })
}

/// The `Statistics` of each query of a BLAST XML report (outfmt 5): those of its first valid
/// context, which NCBI's `CBlastAncillaryData` reads, from the gapped Karlin-Altschul block
/// (blastxml_format.cpp:628-650, data4xmlformat.cpp:236-298): the length adjustment, which the
/// search keeps no field for, from `context_lengths` with the search's own inputs, and the
/// effective search space the search kept. A query without a valid context has `XmlStats::NONE`.
fn blastn_xml_stats(
    report: &BlastnReportInputs<'_>,
    query_lengths: &[usize],
) -> Result<Vec<crate::report::blast_xml::XmlStats>> {
    (0..report.queries.len())
        .map(|q_idx| {
            let Some(blocks) = report.context_karlin[2 * q_idx + report.first_context] else {
                return Ok(crate::report::blast_xml::XmlStats::NONE);
            };
            let eff_space = report.query_eff_searchsp[2 * q_idx + report.first_context];
            let length_adjustment = compute_length_adjustment_ncbi(
                query_lengths[q_idx] as i64,
                report.eff_db_length,
                report.db_num_seqs as i64,
                &blocks.gapped,
            )
            .length_adjustment;
            Ok(crate::report::blast_xml::XmlStats {
                hsp_len: length_adjustment as i32,
                eff_space,
                kappa: blocks.gapped.k,
                lambda: blocks.gapped.lambda,
                entropy: blocks.gapped.h,
            })
        })
        .collect()
}

/// Moves the reported coordinates of every HSP from the searched letters to the records:
/// by the start of the query's interval (`RemapToQueryLoc`) and of the subject's interval
/// (`s_RemapToSubjectLoc`), on both strands. Without ranges nothing moves.
///
/// NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1520-1534
/// ```c++
/// void RemapToQueryLoc(CRef<CSeq_align> sar, const CSeq_loc & query)
/// {
///     _ASSERT(sar);
///     const int query_row = 0;
///
///     TSeqPos q_shift = 0;
///
///     if (query.IsInt()) {
///         q_shift = query.GetInt().GetFrom();
///     }
///
///     if (q_shift > 0) {
///         sar->OffsetRow(query_row, q_shift);
///     }
/// }
/// ```
/// NCBI reference: c++/src/objmgr/util/seq_align_util.cpp:71-78
/// ```c++
///     // Create source seq-loc
///     CSeq_loc src_loc(*id, 0, GetLength(loc, scope) - 1);
///     ENa_strand strand = loc.GetStrand();
///     if (strand != eNa_strand_unknown) {
///         src_loc.SetStrand(strand);
///     }
///     CSeq_loc_Mapper mapper(src_loc, loc, scope);
///     return mapper.Map(align, row);
/// ```
/// A subject interval has the strand both, so its mapping is the shift by its start.
fn shift_to_records(hit_lists: &mut [Option<BlastnHitList>], ranges: &ReportRanges<'_>) {
    let (queries, subjects) = (&ranges.queries.placements, ranges.subjects.placements);
    if !queries.ranged() && !subjects.ranged() {
        return;
    }
    for hit_list in hit_lists.iter_mut().flatten() {
        for hsp_list in &mut hit_list.hsplist_array {
            for hsp in &mut hsp_list.hsps {
                let q_shift = queries.offset(hsp.q_idx as usize);
                let s_shift = subjects.offset(hsp.s_idx as usize);
                hsp.q_start += q_shift;
                hsp.q_end += q_shift;
                hsp.s_start += s_shift;
                hsp.s_end += s_shift;
            }
        }
    }
}

/// The subject ID of the tabular formats (`sseqid`, `sacc`, `saccver`) of every subject, by
/// subject index (`report::defline::tabular_subject_id`: the title's first word, or the
/// local ID, and for `lcl|Subject_...` the first word of `GenerateDefline`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:828-841
/// ```c++
///     bool setSubjectId = (x_IsFieldRequested(eSubjectSeqId) ||
///     		 	 	 	 x_IsFieldRequested(eSubjectGi) ||
///     		 	 	 	 x_IsFieldRequested(eSubjectAccession) ||
///     		 	 	 	 x_IsFieldRequested(eSubjAccessionVersion));
///
///     if(setSubjectIds || setSubjectTaxInfo || setSubjectTaxInfoAll || setSubjectTitle ||
///        x_IsFieldRequested(eSubjectStrand) || setSubjectId)
///     {
///         try {
///        		const CBioseq_Handle& subject_bh =
///                 scope.GetBioseqHandle(align.GetSeq_id(1));
///             if(setSubjectId) {
///             	SetSubjectId(subject_bh);
///             }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:862-869
/// ```c++
///         } catch (const CException&) {
///             list<CRef<CSeq_id> > subject_ids;
///             CRef<CSeq_id> id(new CSeq_id());
///             id->Assign(align.GetSeq_id(1));
///             subject_ids.push_back(id);
///             SetSubjectId(subject_ids);
///             bioseqs_found = false;
///         }
/// ```
/// A subject with a row whose `GenerateDefline` throws (a `Subject_` first word and a title
/// in no encoding that `CUtf8::GuessEncoding` recognises) makes NCBI 2.17.0 write part of the
/// row and stop with a `CCoreException` that names its build's source files (exit 255;
/// `AUTHORITY.md` §J-6, §K-1, RP-20). No result can be defined for it, so LOSAT rejects such
/// a subject explicitly before it writes any format, pending the PD-LOSAT-NCBI-DEFECTS
/// question; a subject without rows is never named.
fn tabular_subject_ids(
    hit_lists: &[Option<BlastnHitList>],
    subjects: &[FastaRecord],
) -> Result<Vec<Arc<[u8]>>> {
    let ids: Vec<Option<Arc<[u8]>>> = subjects
        .iter()
        .map(|subject| {
            crate::report::defline::tabular_subject_id(
                subject.local_id.as_bytes(),
                &subject.title,
                false,
            )
            .ok()
            .map(Arc::from)
        })
        .collect();
    for hit_list in hit_lists.iter().flatten() {
        for hsp_list in &hit_list.hsplist_array {
            let Some(hsp) = hsp_list.hsps.first() else {
                continue;
            };
            let s_idx = hsp.s_idx as usize;
            if ids.get(s_idx).is_some_and(Option::is_none) {
                anyhow::bail!(
                    "subject record {} ({}) has a title whose first word starts with 'Subject_' and whose non-UTF-8 bytes are in no encoding that NCBI BLAST+ recognises; NCBI BLAST+ writes part of a tabular row for such a subject and stops with an exception that names its build's source files (exit 255), which is not supported by LOSAT's BLASTN",
                    s_idx + 1,
                    String::from_utf8_lossy(subjects[s_idx].shown_id())
                );
            }
        }
    }
    Ok(ids
        .into_iter()
        .map(|id| id.unwrap_or_else(|| Arc::from(&b""[..])))
        .collect())
}

/// The end of the traceback for the final hit list of one query (`hit_list`, of a query of
/// `query_len` letters, the length of each of its contexts): the HSP lists cut to
/// `max_hsps_per_subject`, their HSPs below `query_cov_hsp_perc` percent of the query
/// removed (`-qcov_hsp_perc`; a list left empty is removed) and their subject best hits
/// (`-subject_besthit`), the lists sorted by e-value and cut to `hitlist_size`, and the HSPs
/// of each list sorted for the Seq-aligns (`post_process_hits_and_write`; a database search
/// does it after each query batch's traceback, `run_in_pool`).
fn finish_query_hit_list(
    hit_list: &mut BlastnHitList,
    query_len: usize,
    hitlist_size: usize,
    max_hsps_per_subject: usize,
    query_cov_hsp_perc: f64,
    subject_besthit: bool,
    culling_limit: i32,
    ungapped: bool,
) {
    // NCBI's culling pipe is the first thing the traceback does with the consolidated results
    // (`-culling_limit`; the mask level of 101, the default, is no step).
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:1721-1723,1763-1766
    // ```c
    //         /* post-traceback pipes */
    //         BlastHSPStreamTBackClose(hsp_stream, results);
    // ...
    //     	s_FilterBlastResults(results, hit_params->options, query_info, program_number);
    // ```
    if culling_limit > 0 {
        culling::cull_traceback(
            hit_list,
            query_len,
            culling_limit,
            usize::try_from(get_prelim_hitlist_size(hitlist_size, false, !ungapped)).unwrap_or(0),
        );
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:836-867
    // ```c
    // for (subject_index = 0; subject_index < hit_list->hsplist_count; ++subject_index) {
    //     BlastHSPList * hsp_list = hit_list->hsplist_array[subject_index];
    //     if (hit_options->max_hsps_per_subject) {
    //         Blast_TrimHSPListByMaxHsps(hsp_list, hit_options);
    //     }
    //     if(hit_options->query_cov_hsp_perc) {
    //         Blast_HSPListReapByQueryCoverage(hsp_list, hit_options, query_info, program_number);
    //         if(hsp_list->hspcnt == 0){
    //             hit_list->hsplist_array[subject_index] =  Blast_HSPListFree(hsp_list);
    //         }
    //     }
    //     if ((hit_options->hsp_filt_opt != NULL) &&
    //         (hit_options->hsp_filt_opt->subject_besthit_opts != NULL)) {
    //         Blast_HSPListSubjectBestHit(program_number,
    //             hit_options->hsp_filt_opt->subject_besthit_opts,
    //             query_info, hsp_list);
    //     }
    // }
    // if(hit_options->query_cov_hsp_perc) {
    //     Blast_HitListPurgeNullHSPLists(hit_list);
    // }
    // ```
    for hsp_list in &mut hit_list.hsplist_array {
        if max_hsps_per_subject > 0 {
            trim_by_max_hsps(hsp_list, max_hsps_per_subject);
        }
        // Both contexts of a query have its length, and the offsets are those of the HSP's
        // own context.
        if query_cov_hsp_perc > 0.0
            && reap_by_query_coverage(&mut hsp_list.hsps, query_cov_hsp_perc, |hsp| {
                (
                    hsp.internal_q_offset_0 as i64,
                    hsp.internal_q_end_0 as i64,
                    query_len as i64,
                )
            })
        {
            update_best_evalue(hsp_list);
        }
        // NCBI calls Blast_HSPListSubjectBestHit on a list the coverage left empty after
        // freeing it (blast_traceback.c:853-861); LOSAT reads it as the empty list it was.
        if subject_besthit {
            subject_best_hit(&mut hsp_list.hsps, query_len);
        }
    }
    if query_cov_hsp_perc > 0.0 {
        hit_list
            .hsplist_array
            .retain(|hsp_list| !hsp_list.hsps.is_empty());
        hit_list.hsplist_count = hit_list.hsplist_array.len();
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3383-3400
    // ```c
    // Int2 Blast_HSPResultsSortByEvalue(BlastHSPResults* results)
    // {
    //    for (index = 0; index < results->num_queries; ++index) {
    //       hit_list = results->hitlist_array[index];
    //       if (hit_list != NULL && hit_list->hsplist_count > 1) {
    //          qsort(hit_list->hsplist_array, hit_list->hsplist_count,
    //                sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
    //       }
    //       s_BlastHitListPurge(hit_list);
    //    }
    // }
    // ```
    hit_list.sort_by_evalue();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:877-892
    // ```c
    // static void s_BlastPruneExtraHits(BlastHSPResults* results, Int4 hitlist_size)
    // {
    //    for (subject_index = hitlist_size;
    //         subject_index < hit_list->hsplist_count; ++subject_index) {
    //       hit_list->hsplist_array[subject_index] =
    //       Blast_HSPListFree(hit_list->hsplist_array[subject_index]);
    //    }
    //    hit_list->hsplist_count = MIN(hit_list->hsplist_count, hitlist_size);
    // }
    // ```
    hit_list.prune_by_size(hitlist_size);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:1571-1577
    // ```c
    //         // Sort HSPs with e-values as first priority and scores as
    //         // tie-breakers, since that is the order we want to see them in
    //         // in Seq-aligns.
    //         Blast_HSPListSortByEvalue(hsp_list);
    // ```
    // Every format shows the HSPs of a subject in this order. It is the score order,
    // except where the two strands of a query have different gapped blocks.
    for hsp_list in &mut hit_list.hsplist_array {
        sort_hsplist_by_evalue(hsp_list);
    }
}

/// The HSPs of score 0 dropped from the finished hit list `hit_list`, and the HSP lists
/// left without HSPs (`post_process_hits_and_write`).
fn drop_hsps_of_score_zero(hit_list: &mut BlastnHitList) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:672-674
    // ```c
    //     if (hsp->score == 0) {
    //         return CRef<CSeq_align>();
    //     }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:1485-1490
    // ```c
    //             seqalign =
    //                 s_BlastHSP2SeqAlign(program, hsp, query_id, subject_id,
    //                                     query_length, subject_length);
    //         }
    //
    //         if (seqalign.Empty()) continue;
    // ```
    // An HSP of score 0 (the traceback reads the ambiguity codes that the preliminary
    // search read as random bases) has no Seq-align, so no format shows it, and a
    // subject without other HSPs has no alignments.
    for hsp_list in &mut hit_list.hsplist_array {
        hsp_list.hsps.retain(|hsp| hsp.raw_score != 0);
    }
    hit_list
        .hsplist_array
        .retain(|hsp_list| !hsp_list.hsps.is_empty());
    hit_list.hsplist_count = hit_list.hsplist_array.len();
}

fn post_process_hits_and_write(
    mut hit_lists: Vec<Option<BlastnHitList>>,
    hitlist_size: usize,
    max_hsps_per_subject: usize,
    query_cov_hsp_perc: f64,
    subject_besthit: bool,
    query_lengths: &[usize],
    outputs: &mut ReportOutputs<'_>,
    verbose: bool,
    query_ids: &[Arc<[u8]>],
    output_formats: &[BlastnOutputFormat],
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
    timing: Option<&BlastnTiming>,
    report: &BlastnReportInputs<'_>,
    unsearched: Vec<bool>,
    query_warnings: &[Vec<u8>],
    mut sam: Vec<Option<SamTarget<'_>>>,
) -> Result<()> {
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:183-187
    // ```c
    // typedef struct BlastHSPResults {
    //    Int4 num_queries;
    //    BlastHitList** hitlist_array;
    // } BlastHSPResults;
    // ```
    if verbose {
        let total_hits: usize = hit_lists
            .iter()
            .filter_map(|list| list.as_ref())
            .map(|list| {
                list.hsplist_array
                    .iter()
                    .map(|l| l.hsps.len())
                    .sum::<usize>()
            })
            .sum();
        eprintln!(
            "[INFO] Received {} raw hits total, starting post-processing...",
            total_hits
        );
    }
    let chain_start = std::time::Instant::now();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:836-867
    // ```c
    // for (subject_index = 0; subject_index < hit_list->hsplist_count; ++subject_index) {
    //     BlastHSPList * hsp_list = hit_list->hsplist_array[subject_index];
    //     if (hit_options->max_hsps_per_subject) {
    //         Blast_TrimHSPListByMaxHsps(hsp_list, hit_options);
    //     }
    //     if ((hit_options->hsp_filt_opt != NULL) &&
    //         (hit_options->hsp_filt_opt->subject_besthit_opts != NULL)) {
    //         Blast_HSPListSubjectBestHit(program_number,
    //             hit_options->hsp_filt_opt->subject_besthit_opts,
    //             query_info, hsp_list);
    //     }
    // }
    // ```
    // A database search has finished the hit lists of each query batch after its traceback
    // (`run_in_pool`: NCBI reads the definition lines of the batch's final hit lists between
    // the two steps), so they are not finished again.
    if !report.finished_hit_lists {
        for (q_idx, hit_list_opt) in hit_lists.iter_mut().enumerate() {
            let hit_list = match hit_list_opt {
                Some(value) => value,
                None => continue,
            };
            finish_query_hit_list(
                hit_list,
                query_lengths.get(q_idx).copied().unwrap_or_default(),
                hitlist_size,
                max_hsps_per_subject,
                query_cov_hsp_perc,
                subject_besthit,
                report.culling_limit,
                report.ungapped,
            );
            drop_hsps_of_score_zero(hit_list);
        }
    }

    let total_hits: usize = hit_lists
        .iter()
        .filter_map(|list| list.as_ref())
        .map(|list| {
            list.hsplist_array
                .iter()
                .map(|l| l.hsps.len())
                .sum::<usize>()
        })
        .sum();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:877-892
    // ```c
    // static void s_BlastPruneExtraHits(BlastHSPResults* results, Int4 hitlist_size)
    // ```
    if let Some(timing) = timing {
        let post_elapsed = chain_start.elapsed();
        let post_elapsed_ns = post_elapsed.as_nanos() as u64;
        BlastnTiming::record_ns(&timing.traceback_ns, post_elapsed_ns);
        BlastnTiming::record_ns(&timing.traceback_postprocess_ns, post_elapsed_ns);
        timing
            .traceback_calls
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if verbose {
            eprintln!(
                "[INFO] Post-processing done in {:.2}s, {} hits after filtering, writing output...",
                post_elapsed.as_secs_f64(),
                total_hits
            );
        }
    } else if verbose {
        eprintln!(
            "[INFO] Post-processing done in {:.2}s, {} hits after filtering, writing output...",
            chain_start.elapsed().as_secs_f64(),
            total_hits
        );
    }
    let write_start = std::time::Instant::now();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3077-3106
    // ```c
    // static int s_EvalueCompareHSPLists(const void* v1, const void* v2) { ... }
    // ```
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-467
    // ```c
    // CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
    // formatter.PrintProlog();
    // ...
    // ITERATE(CSearchResultSet, result, *results) {
    //     ...
    //         formatter.PrintOneResultSet(**result, queries);
    //     ...
    // }
    // ```
    // Each requested format prints the same final result without searching again.
    // `unsearched` marks the queries of the batches without a valid context, which NCBI
    // does not search (outfmt 0 and 7 show them without results).
    // Rendered hits are needed by outfmt 0 and by the caller's hit records; both get the
    // same final hit list, built once.
    // NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1630-1637
    // ```c++
    // 	if (seqinfo_src->CanReturnPartialSequence() == true)
    // 	{
    //         	CConstRef<CSeq_loc> subj_loc = seqinfo_src->GetSeqLoc(kOid);
    //         	NON_CONST_ITERATE(vector<CRef<CSeq_align > >, iter, hit_align) {
    //              	   RemapToQueryLoc(*iter, query_loc);
    //              	   if ( !is_ooframe )
    //                    	s_RemapToSubjectLoc(*iter, *subj_loc);
    // ```
    // The reports print record coordinates: the HSPs move by the start of the query's and
    // the subject's interval (`shift_to_records`); their internal offsets stay those of
    // the searched letters, from which the shown rows are read.
    shift_to_records(&mut hit_lists, report.ranges);
    let wants_rows =
        outputs.hits.is_some() || output_formats.contains(&BlastnOutputFormat::Pairwise);
    // EXPERIMENT (LOSAT_X_RPTSTREAM): the same masks, made once for every batch.
    let subject_masks = (wants_rows && report.lcase_masking).then(|| {
        // The records as read with `-subject_loc` (never packed), else the searched ones.
        let subject_records = match report.ranges.subjects.records {
            Some(records) => SubjectLetters::new(records, None),
            None => report.subjects,
        };
        match report
            .subject_cache
            .and_then(|cache| cache.masks.as_deref())
        {
            Some(masks) => std::borrow::Cow::Borrowed(masks),
            None => std::borrow::Cow::Owned(
                (0..subject_records.records().len())
                    .map(|index| subject_records.lowercase_masks(index))
                    .collect::<Vec<_>>(),
            ),
        }
    });
    // The shown masks of every query and subject, made once for every window of queries.
    let rows = wants_rows.then(|| {
        PairwiseRows::new(
            report.queries,
            report.subjects,
            &DisplayMasks {
                query: report.query_masks,
                query_offsets: &report.ranges.queries.placements,
                subject: subject_masks.as_deref(),
            },
        )
    });
    // The rendered hits of the queries `queries`.
    let rows_of = |queries: std::ops::Range<usize>| -> Result<Vec<PairwiseHit>> {
        let rows = rows
            .as_ref()
            .expect("the rows' inputs are made for outfmt 0");
        let mut hits = pairwise_hits(&hit_lists[queries], rows, report.x_pairwise_pool)?;
        // NCBI reference: c++/src/objtools/align_format/showalign.cpp:2479
        // ```c++
        //     out<<kLengthString<<bsp_handle.GetBioseqLength()<<"\n";
        // ```
        // The length of the subject's record, also where `-subject_loc` searches part
        // of it.
        for hit in &mut hits {
            let s_idx = hit.hit.s_idx as usize;
            hit.subject_length = hit
                .subject_length
                .map(|searched| report.ranges.subjects.placements.length(s_idx, searched));
        }
        Ok(hits)
    };
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1411-1585
    // ```c++
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ...
    //     display.DisplaySeqalign(m_Outfile);
    // ```
    // NCBI makes and prints the alignments of one query's result set at a time; outfmt 0 makes
    // the rows of a window of queries, writes them and drops them before the next window's
    // (`pairwise::pairwise_windows`, line DBID, ID4f). The rows' errors are LOSAT's internal
    // invariants (`pairwise_hits`), so making them while the report is written changes no
    // run's output. The caller's hit records and a format observer (whose HSP numbers count
    // over the whole report) read every query's rows: those are made once, before any format
    // is written.
    let pairwise_hits = if outputs.hits.is_some()
        || (outputs.observer.is_some() && output_formats.contains(&BlastnOutputFormat::Pairwise))
    {
        Some(rows_of(0..hit_lists.len())?)
    } else {
        None
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1411
    // ```c
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ```
    // The caller receives the same final result that every formatter prints.
    if let (Some(hits_sink), Some(hits)) = (outputs.hits.as_mut(), pairwise_hits.as_ref()) {
        hits_sink(hits);
    }
    // Built before any format is written, so that a report that cannot be made fails the
    // run without output.
    let pairwise_data = if output_formats.contains(&BlastnOutputFormat::Pairwise) {
        Some(blastn_pairwise_report(
            report,
            unsearched.clone(),
            query_titles,
            subject_title,
        )?)
    } else {
        None
    };
    // EXPERIMENT (LOSAT_X_RPTSTREAM): a subject set's ids made once for every batch; the
    // switch acts only when no subject's id is refused (`x_rpt_stream::SubjectCache`), so
    // they are the ids made here.
    let subject_ids: std::borrow::Cow<'_, [Arc<[u8]>]> =
        if output_formats.iter().any(|format| format.is_tabular()) {
            match (report.database_subject_ids, report.subject_cache) {
                (Some(ids), _) => std::borrow::Cow::Borrowed(ids),
                (None, Some(cache)) => std::borrow::Cow::Borrowed(cache.ids.as_slice()),
                (None, None) => std::borrow::Cow::Owned(tabular_subject_ids(
                    &hit_lists,
                    report.subjects.records(),
                )?),
            }
        } else {
            std::borrow::Cow::Owned(Vec::new())
        };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-498
    // ```c++
    //     //get defline
    //     sdl->defline = CDeflineGenerator().GenerateDefline(m_ScopeRef->GetBioseqHandle(*(sdl->id)), sequence::CDeflineGenerator::fLeavePrefixSuffix);
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2273
    // ```c++
    // 	alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
    // ```
    // The outfmt 0 descriptions and headings are made from the subjects' title bytes
    // (`report::defline::generate_defline`), by subject index.
    // EXPERIMENT (LOSAT_X_RPTSTREAM): the same titles, made once for every batch.
    let subject_titles: std::borrow::Cow<'_, [Arc<[u8]>]> =
        if output_formats.contains(&BlastnOutputFormat::Pairwise) {
            match report
                .subject_cache
                .and_then(|cache| cache.titles.as_deref())
            {
                Some(titles) => std::borrow::Cow::Borrowed(titles),
                None => report.subjects.records().iter().map(record_title).collect(),
            }
        } else {
            std::borrow::Cow::Owned(Vec::new())
        };
    let observer = &mut outputs.observer;
    // The warnings are written once, with the first format (`QueryWarnings`).
    let mut warning_lines = Some(query_warnings);
    for (format_index, (format, &output_format)) in
        outputs.formats.iter_mut().zip(output_formats).enumerate()
    {
        // The warnings that a format other than the archive writes are not the archive's
        // messages (a search with several formats, `run_local`).
        if let Some(handle) = report.archive {
            handle.set_muted(output_format != BlastnOutputFormat::Archive);
        }
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        // The taxonomy database's messages of XML2 and JSON go to the diagnostics at every
        // format.
        // So do those of the formatters of a tabular field list that shows a taxonomy name
        // (`x_CheckTaxDB`, `algorithm::blastn::tabular`).
        let mut format_warnings = (warning_lines.is_some()
            || matches!(
                output_format,
                BlastnOutputFormat::Xml2(_) | BlastnOutputFormat::Archive
            )
            || output_format
                .tabular_spec()
                .is_some_and(|spec| spec.checks_taxdb()))
        .then(|| crate::report::query_warnings::QueryWarnings {
            before: warning_lines.take().unwrap_or(&[]),
            sink: &mut *outputs.diagnostics,
        });
        // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:68-93
        // ```c
        // CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
        //     : m_FormatType(format_type), ..., m_Outfile(outfile),
        // ```
        let mut writer = format.sink.open()?;
        if output_format == BlastnOutputFormat::Pairwise {
            let (queries, pairwise_report) = pairwise_data
                .as_ref()
                .expect("the pairwise report data is built for outfmt 0");
            let config = PairwiseConfig {
                program: "blastn".to_string(),
                show_frame: false,
                ..PairwiseConfig::default()
            };
            let mut write = |hits: &[PairwiseHit], range: std::ops::Range<usize>| {
                write_blastn_pairwise_report(
                    hits,
                    &mut writer,
                    &config,
                    queries,
                    range,
                    &subject_titles,
                    pairwise_report,
                    probe.as_mut(),
                    format_warnings.as_mut(),
                    report.x_pairwise_pool,
                )
            };
            match pairwise_hits.as_deref() {
                Some(hits) => write(hits, 0..queries.len())?,
                None => {
                    // The windows of the parallel writer (`x_pairwise_par`).
                    let windows =
                        pairwise_windows(&hit_lists, queries.len(), report.x_pairwise_pool);
                    if windows.is_empty() {
                        // No query: the prolog and the epilog alone.
                        write(&[], 0..0)?;
                    }
                    for range in windows {
                        let lists =
                            range.start.min(hit_lists.len())..range.end.min(hit_lists.len());
                        let hits = rows_of(lists)?;
                        write(&hits, range)?;
                    }
                }
            }
        } else if let BlastnOutputFormat::Seqalign(serial) = output_format {
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
            // ```c++
            //     if (m_FormatType == CFormattingArgs::eAsnText || m_FormatType == CFormattingArgs::eJsonSeqalign) {
            //         if (results.HasAlignments()) {
            //     	    CRef<CSeq_align_set> aln_set (new CSeq_align_set);
            //             CBlastFormatUtil::PruneSeqalign(*(results.GetSeqAlign()), *aln_set, m_HitlistSize);
            //             if(m_FormatType == CFormattingArgs::eAsnText)
            //             	m_Outfile << MSerial_AsnText << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
            //             else
            //             	m_Outfile << MSerial_Json << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
            //         }
            // ```
            // The Seq-annot of each query with alignments (`seqalign`).
            let subjects = match (report.database_ordinals, report.database) {
                (Some(ordinals), Some(database)) => SubjectNames::Database {
                    ordinals,
                    title: &database.definition,
                    parsed: report.database_parsed,
                },
                _ => SubjectNames::Records(report.subjects.records()),
            };
            write_blastn_seqalign(
                serial,
                &hit_lists,
                &QueryNames {
                    records: report.queries,
                    searched_lengths: query_lengths,
                    placements: &report.ranges.queries.placements,
                },
                subjects,
                report.ranges.subjects.placements,
                &mut writer,
                format_warnings.as_mut(),
            )?;
        } else if output_format == BlastnOutputFormat::Archive {
            // The BLAST archive (outfmt 11): one document per query batch of the call.
            let handle = report
                .archive
                .expect("a BLAST archive request has its handle");
            let subjects = match (report.database_ordinals, report.database) {
                (Some(ordinals), Some(database)) => SubjectNames::Database {
                    ordinals,
                    title: &database.definition,
                    parsed: report.database_parsed,
                },
                _ => SubjectNames::Records(report.subjects.records()),
            };
            let call = archive::ArchiveCall {
                handle,
                queries: report.queries,
                query_lengths,
                query_placements: &report.ranges.queries.placements,
                subject_records: report.database.is_none().then(|| {
                    report
                        .ranges
                        .subjects
                        .records
                        .unwrap_or_else(|| report.subjects.records())
                }),
                subjects,
                subject_placements: report.ranges.subjects.placements,
                hit_lists: &hit_lists,
                query_masks: report.query_masks,
                dust_masks: report.archive_dust,
                unsearched: &unsearched,
                lcase_masking: report.lcase_masking,
                window_masker: report.window_masker_db.is_some(),
                both_strands: report.both_strands,
                minus_strand: report.minus_strand,
                spans: report.batch_spans,
            };
            archive::write_blastn_archives(
                &call,
                format_warnings.as_mut(),
                &|first| archive_first_context(report, query_lengths, first),
                &mut writer,
            )?;
        } else if let BlastnOutputFormat::Sam(_) = output_format {
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:750-755
            // ```c++
            //     if (m_FormatType == CFormattingArgs::eSAM) {
            //         if (results.HasAlignments()) {
            //             m_SamFormatter->Print(*results.GetSeqAlign());
            //         }
            //         return;
            //     }
            // ```
            // The queries' records are added to the run's document, which NCBI writes when its
            // formatter is destroyed (`SamDocument::finish`): after the last call, or in a call
            // that is the run's only one.
            let target = sam[format_index]
                .as_mut()
                .expect("a SAM document is made for outfmt 17");
            if let Some(warnings) = format_warnings.as_mut() {
                for q_idx in 0..hit_lists.len() {
                    warnings.before_query(q_idx, &mut writer)?;
                }
            }
            let names = QueryNames {
                records: report.queries,
                searched_lengths: query_lengths,
                placements: &report.ranges.queries.placements,
            };
            match (
                report.database_ordinals,
                report.database,
                report.database_reader,
            ) {
                (Some(ordinals), Some(database), Some(reader)) => {
                    // A database's sequences are searched whole; the lengths come from its
                    // index (no letters).
                    let lengths = reader.sam_lengths(ordinals)?;
                    let cache = database::SamLetterCache::default();
                    let letters_of = |rank: usize, from: usize, to: usize| {
                        reader.sam_letters(&cache, ordinals[rank], from, to)
                    };
                    super::super::sam::add_blastn_sam(
                        target.document,
                        &hit_lists,
                        &names,
                        report.queries_before,
                        SubjectNames::Database {
                            ordinals,
                            title: &database.definition,
                            parsed: report.database_parsed,
                        },
                        report.ranges.subjects.placements,
                        &super::super::sam::SubjectLengths {
                            searched: &|rank| lengths[rank],
                            whole: &|rank| lengths[rank],
                        },
                        &mut super::super::sam::BlastnSamLetters {
                            queries: report.queries,
                            query_placements: &report.ranges.queries.placements,
                            subjects: super::super::sam::SubjectSource::Database(&letters_of),
                        },
                    )?;
                }
                _ => {
                    // The records as read with `-subject_loc` (never packed), else the
                    // searched ones: the letters are asked for in record coordinates.
                    let whole_records = match report.ranges.subjects.records {
                        Some(records) => SubjectLetters::new(records, None),
                        None => report.subjects,
                    };
                    super::super::sam::add_blastn_sam(
                        target.document,
                        &hit_lists,
                        &names,
                        report.queries_before,
                        SubjectNames::Records(report.subjects.records()),
                        report.ranges.subjects.placements,
                        &super::super::sam::SubjectLengths {
                            searched: &|s_idx| report.subjects.len(s_idx),
                            whole: &|s_idx| {
                                report
                                    .ranges
                                    .subjects
                                    .placements
                                    .length(s_idx, report.subjects.len(s_idx))
                            },
                        },
                        &mut super::super::sam::BlastnSamLetters {
                            queries: report.queries,
                            query_placements: &report.ranges.queries.placements,
                            subjects: super::super::sam::SubjectSource::Records(whole_records),
                        },
                    )?;
                }
            }
            if target.finish {
                target.document.finish(&mut writer)?;
            }
        } else if output_format.is_structured_xml() {
            // The XML2 / JSON documents of outfmt 15 and 16 take the path of outfmt 5: the same
            // per-query records, written by another writer.
            let xml2_kind = output_format.xml2_kind();
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1441
            // ```c++
            //     if (m_FormatType == CFormattingArgs::eAsnText
            //       || m_FormatType == CFormattingArgs::eAsnBinary
            //       || m_FormatType == CFormattingArgs::eXml
            //       ...
            //     {
            //         x_PrintStructuredReport(results, queries);
            //         return;
            //     }
            // ```
            // The CBlastOutput document of the queries (`blast_xml`).
            let is_database = report.database.is_some();
            let xml_stats = blastn_xml_stats(report, query_lengths)?;
            let xml_messages: Vec<Vec<String>> = report
                .queries
                .iter()
                .enumerate()
                .map(|(q_idx, query)| {
                    let mut messages = Vec::new();
                    if report.context_karlin[2 * q_idx + report.first_context].is_none() {
                        messages.push(if query.seq().is_empty() {
                            NO_DATA_MESSAGE.to_string()
                        } else {
                            crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string()
                        });
                    }
                    if !unsearched[q_idx]
                        && report
                            .database_reader
                            .is_some_and(|database| database.filtered_empty())
                    {
                        messages.push(
                            crate::report::query_warnings::EMPTY_DATABASE_MESSAGE.to_string(),
                        );
                    }
                    messages
                })
                .collect();
            let version = format!("BLASTN {NCBI_BLASTN_VERSION}");
            let filter = crate::report::blast_xml::filter::filter_string(
                report.dust,
                None,
                report
                    .window_masker_db
                    .map(crate::report::blast_xml::filter::WindowMasker::Database),
                true,
            );
            let subject_records = report.subjects.records();
            let database_length = |rank: usize| report.subjects.len(rank);
            let record_length = |s_idx: usize| {
                report
                    .ranges
                    .subjects
                    .placements
                    .length(s_idx, report.subjects.len(s_idx))
            };
            // The taxids of the reported sequences and the taxonomy database's names, read for
            // the XML2 and JSON formats only.
            let (xml_taxids, xml_taxonomy, xml_parsed_taxonomy, xml_taxdb_notes): (
                Vec<i32>,
                Vec<crate::report::blast_xml::HitTaxonomy>,
                Vec<Vec<crate::report::blast_xml::HitTaxonomy>>,
                crate::algorithm::common::database::TaxdbNotes,
            ) = match (xml2_kind, report.database_ordinals, report.database_reader) {
                (Some(_), Some(ordinals), Some(database)) if is_database => {
                    database.read_taxids(ordinals, report.database_parsed, &hit_lists)?
                }
                _ => (Vec::new(), Vec::new(), Vec::new(), Default::default()),
            };
            // The masked regions of the queries in record coordinates and the queries of the
            // batches that were not searched, which only XML2 shows (outfmt 5 shows neither).
            let (xml_masks, xml_searched): (
                Vec<Vec<crate::report::blast_xml::QueryMask>>,
                Vec<bool>,
            ) = if xml2_kind.is_some() {
                (
                    (0..report.queries.len())
                        .map(|q_idx| {
                            shown_query_masks(
                                &report.query_masks[q_idx],
                                query_lengths[q_idx],
                                report.ranges.queries.placements.offset(q_idx),
                            )
                            .iter()
                            .map(|mask| crate::report::blast_xml::QueryMask {
                                from: mask.start as i32,
                                to: mask.end as i32 - 1,
                                frame: 0,
                            })
                            .collect()
                        })
                        .collect(),
                    unsearched.iter().map(|&skipped| !skipped).collect(),
                )
            } else {
                (Vec::new(), Vec::new())
            };
            let xml_subjects = match report.database_ordinals {
                Some(ordinals) if is_database => crate::report::blast_xml::XmlSubjects::Database {
                    ordinals,
                    titles: crate::report::blast_xml::XmlTitles::Records(subject_records),
                    // Only XML2 shows the taxid.
                    taxids: &xml_taxids,
                    taxonomy: &xml_taxonomy,
                    length: &database_length,
                    parsed: report.database_parsed,
                    parsed_taxonomy: &xml_parsed_taxonomy,
                },
                _ => crate::report::blast_xml::XmlSubjects::Records {
                    records: subject_records,
                    length: &record_length,
                    protein: false,
                },
            };
            let xml = super::super::blast_xml::BlastnXml {
                run: crate::report::blast_xml::XmlRun {
                    program: "blastn",
                    version: &version,
                    reference: if report.megablast {
                        crate::report::blast_xml::REFERENCE_MEGABLAST
                    } else {
                        crate::report::blast_xml::REFERENCE_GAPPED_BLAST
                    },
                    db_name: if is_database {
                        subject_title.as_bytes()
                    } else {
                        b""
                    },
                    matrix: None,
                    expect: report.expect,
                    gap_open: report.gap_open,
                    gap_extend: report.gap_extend,
                    reward_penalty: Some((report.reward, report.penalty)),
                    filter: &filter,
                    db_num: if is_database {
                        report.shown_num_seqs() as i32
                    } else {
                        0
                    },
                    db_len: if is_database {
                        report.shown_total_letters() as i64
                    } else {
                        0
                    },
                    positives: None,
                    mask_char: b'N',
                },
                queries: report.queries,
                searched_lengths: query_lengths,
                stats: xml_stats,
                messages: xml_messages,
                masks: xml_masks,
                searched: xml_searched,
                subjects: xml_subjects,
                taxdb_notes: xml_taxdb_notes,
                letters: report.subjects,
                epilog: report.epilog,
            };
            match xml2_kind {
                None => super::super::blast_xml::write_blastn_xml(
                    &xml,
                    &mut crate::report::blast_xml::BlastXmlWriter::new(report.queries_before),
                    &hit_lists,
                    &mut writer,
                    format_warnings.as_mut(),
                )?,
                Some(kind) => {
                    // The gap costs are written for a gapped search (NCBI reads them from
                    // the options when `GetGappedMode()`); `cbs` and the genetic codes
                    // are not those of a nucleotide search (blastxml2_format.cpp:531-546).
                    let params = crate::report::blast_xml2::Xml2Run {
                        target: if is_database {
                            crate::report::blast_xml2::Xml2Target::Database
                        } else {
                            crate::report::blast_xml2::Xml2Target::Subjects(subject_records.len())
                        },
                        gapped: !report.ungapped,
                        cbs: 0,
                        query_gencode: None,
                        db_gencode: None,
                    };
                    let mut document =
                        crate::report::blast_xml2::Xml2Writer::for_run(kind, params, report.out);
                    super::super::blast_xml::write_blastn_xml(
                        &xml,
                        &mut document,
                        &hit_lists,
                        &mut writer,
                        format_warnings.as_mut(),
                    )?;
                }
            }
        } else if let Some(spec) = output_format.tabular_spec() {
            // The tabular reports with a field list or a delimiter, and outfmt 10 and 20
            // (`algorithm::blastn::tabular`); never the default 6 and 7 below.
            write_custom_tabular(
                spec,
                report,
                &hit_lists,
                query_ids,
                query_titles,
                query_lengths,
                subject_title,
                &subject_ids,
                &unsearched,
                &mut writer,
                format_warnings.as_mut(),
                probe.as_mut(),
            )?;
        } else {
            crate::algorithm::blastn::x_tab_par::write_output(
                report.x_tab_pool,
                &hit_lists,
                &mut writer,
                query_ids,
                &subject_ids,
                output_format,
                query_titles,
                subject_title,
                &unsearched,
                report
                    .epilog
                    .then_some(report.queries_before + hit_lists.len()),
                probe.as_mut(),
                format_warnings.as_mut(),
            )?;
        }
        writer.flush()?;
    }
    // Without an output format the warnings are written all at once.
    if let Some(handle) = report.archive {
        handle.set_muted(false);
    }
    if let Some(before) = warning_lines {
        for query_warnings in before {
            outputs.diagnostics.write_all(query_warnings)?;
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-782
    // ```c
    // // tabular formatting just prints each alignment in turn
    // if (m_FormatType == CFormattingArgs::eTabular ||
    //     m_FormatType == CFormattingArgs::eTabularWithComments ||
    //     m_FormatType == CFormattingArgs::eCommaSeparatedValues ||
    //     m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader) {
    //   CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
    // ```
    if let Some(timing) = timing {
        let write_elapsed = write_start.elapsed();
        timing.format_ns.fetch_add(
            write_elapsed.as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
        timing
            .format_calls
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if verbose {
            eprintln!(
                "[INFO] Output written in {:.2}s",
                write_elapsed.as_secs_f64()
            );
        }
    } else if verbose {
        eprintln!(
            "[INFO] Output written in {:.2}s",
            write_start.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

/// The custom tabular report of the queries of `hit_lists` (`algorithm::blastn::tabular`): the
/// subjects' values as the report shows them (a database's stored lengths and, where a field
/// shows them, its taxids), outfmt 20's header line with the run's first report and outfmt 7's
/// last line with the epilog. Kept out of line: only a field list reaches it.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-846
/// ```c++
///        		const CBioseq_Handle& subject_bh =
///                 scope.GetBioseqHandle(align.GetSeq_id(1));
///             if(setSubjectId) {
///             	SetSubjectId(subject_bh);
///             }
///             m_SubjectLength = subject_bh.GetBioseqLength();
///
///             if(setSubjectIds || setSubjectTaxInfo || setSubjectTitle || setSubjectTaxInfoAll) {
///             	CRef<CBlast_def_line_set> bdlRef =
///             			CSeqDB::ExtractBlastDefline(subject_bh);
/// ```
/// `slen` is the length of the subject's Bioseq: the whole record with `-subject_loc`, as the
/// outfmt 0 `Length=` line.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn write_custom_tabular(
    spec: &crate::report::tabular::TabularSpec,
    report: &BlastnReportInputs<'_>,
    hit_lists: &[Option<BlastnHitList>],
    query_ids: &[Arc<[u8]>],
    query_titles: &[Arc<[u8]>],
    query_lengths: &[usize],
    subject_title: &str,
    subject_ids: &[Arc<[u8]>],
    unsearched: &[bool],
    writer: &mut dyn std::io::Write,
    warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
    probe: Option<&mut FormatProbe<'_>>,
) -> Result<()> {
    use super::super::tabular::{write_blastn_tabular, BlastnTabular};
    use crate::algorithm::common::tabular::TabularSubjects;
    let needs = spec.needs();
    let record_length = |s_idx: usize| {
        report
            .ranges
            .subjects
            .placements
            .length(s_idx, report.subjects.len(s_idx))
    };
    let (lengths, taxonomy) = match (report.database_ordinals, report.database_reader) {
        (Some(ordinals), Some(database)) => (
            if spec.has(crate::report::tabular::TabularField::SubjectLength) {
                database.sam_lengths(ordinals)?
            } else {
                Vec::new()
            },
            if needs.taxonomy {
                database.tabular_taxonomy(ordinals, report.database_parsed, hit_lists)?
            } else {
                Vec::new()
            },
        ),
        _ => (Vec::new(), Vec::new()),
    };
    let names = match report.database_reader {
        Some(database) => database.tabular_names(spec, &taxonomy)?,
        None => None,
    };
    let subjects = match report.database_ordinals {
        Some(ordinals) if report.database_reader.is_some() => TabularSubjects::Database {
            ids: subject_ids,
            ordinals,
            records: report.subjects.records(),
            parsed: report.database_parsed,
            lengths: &lengths,
            taxonomy: &taxonomy,
        },
        _ => TabularSubjects::Records {
            ids: subject_ids,
            length: &record_length,
        },
    };
    let tab = BlastnTabular {
        spec,
        queries: report.queries,
        query_ids,
        query_titles,
        searched_lengths: query_lengths,
        range_length: report.ranges.queries.range_length,
        database: subject_title,
        subjects,
        letters: report.subjects,
        unsearched,
        prolog: report.queries_before == 0,
        epilog: report
            .epilog
            .then_some(report.queries_before + hit_lists.len()),
        taxdb_found: report.tabular_taxdb_found,
        names: names.as_ref(),
    };
    write_blastn_tabular(&tab, hit_lists, writer, warnings, probe)
}

/// `settings` are the NCBI application settings that LOSAT reproduces
/// (`ncbi_environment::check_ncbi_application_settings`): the data loaders decide whether
/// the readers of the subjects and the queries try a first line as a Seq-id
/// (`ReaderConfig`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:67-72
/// ```c++
///     CNcbiApplication* app = CNcbiApplication::Instance();
///     if (app) {
///         const CNcbiRegistry& registry = app->GetConfig();
///         x_LoadDataLoadersConfig(registry);
///         x_LoadBlastDbDataLoaderConfig(registry);
///     }
/// ```
pub fn run(
    mut args: BlastnArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:128-133
    // ```c
    // if(RecoverSearchStrategy(args, m_CmdLineArgs)) {
    // 	m_OptsHndl.Reset(&*m_CmdLineArgs->SetOptionsForSavedStrategy(args));
    // }
    // else {
    // 	m_OptsHndl.Reset(&*m_CmdLineArgs->SetOptions(args));
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
    // ```c
    //     if (GetExportSearchStrategyStream(args) ||
    //            m_FormattingArgs->ArchiveFormatRequested(args)) {
    //         locality = CBlastOptions::eBoth;
    //     }
    // ```
    // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
    // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
    // there too (`cli_output_formats`).
    let output_formats = cli_output_formats(&args.outfmt, args.out.as_deref())?;
    // The BLAST archive (outfmt 11): its request, and the refusals of what LOSAT cannot
    // write, before anything is read or written.
    let archive_requested = output_formats.contains(&BlastnOutputFormat::Archive);
    // The error archives follow the `-outfmt` text (`writes_error_archive`).
    args.archive = archive::handle_for(&args, args.dust_filtering.as_deref(), archive_requested)?
        .map(|handle| handle.with_outfmt_text(&args.outfmt));
    // LOSAT's thread capability is checked before any input or output, as before.
    crate::utils::threading::validate_threads(args.num_threads)?;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3631-3636
    // ```c
    //     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
    //         (*arg)->ExtractAlgorithmOptions(args, opts);
    //     }
    //
    //     m_IsUngapped = !opts.GetGappedMode();
    //     try { retval->Validate(); }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blastn_args.cpp:63-70
    // ```c
    //     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
    //     m_BlastDbArgs->SetDatabaseMaskingSupport(true);
    //     arg.Reset(m_BlastDbArgs);
    //     m_Args.push_back(arg);
    //
    //     m_StdCmdLineArgs.Reset(new CStdCmdLineArgs);
    //     arg.Reset(m_StdCmdLineArgs);
    //     m_Args.push_back(arg);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2553-2557
    // ```c
    //         CRef<blast::CBlastQueryVector> subjects;
    //         m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(),
    //                                        subj_range, parse_deflines,
    //                                        use_lcase_masks, subjects, m_IsMapper);
    //         m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects));
    // ```
    // The first handler opens and reads the subjects (an empty subject set fails there),
    // before the query and the output are opened and the options are checked. Files are
    // opened when a handler asks for them, and each once: a named pipe gives its bytes to
    // one reader only.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2558-2562
    // ```c
    //     } else if (!m_IsIgBlast){
    //         // IgBlast permits use of germline database
    //         NCBI_THROW(CInputException, eInvalidInput,
    //            "Either a BLAST database or subject sequence(s) must be specified");
    //     }
    // ```
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:172-176
    // ```c
    //     catch (const blast::CInputException& e) {                               \
    //         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
    //         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
    //     }                                                                       \
    // ```
    // The handler of the database arguments raises it, before the query and the output
    // are opened (`CStdCmdLineArgs` comes after it).
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
    // ```c++
    //     if (args.Exist(kArgDb) && args[kArgDb]) {
    // 	std::string local_dblist = NStr::TruncateSpaces( args[kArgDb].AsString() );
    //
    //         m_SearchDb.Reset(new CSearchDatabase( local_dblist,
    //                                              mol_type));
    // ```
    // With `-db` the handler of the database arguments only keeps the name: the database
    // is opened after every option handler has run (`search_cli`). With `-subject` it
    // reads the subjects here, before the query and the output are opened.
    let archive = args.archive.clone();
    let mut stderr = std::io::stderr();
    let mut archive_tee;
    let diagnostics: &mut (dyn std::io::Write + Send) = match &archive {
        Some(handle) => {
            archive_tee = crate::report::archive::DiagTee::new(&mut stderr, handle.clone());
            &mut archive_tee
        }
        None => &mut stderr,
    };
    let cli_subjects = if args.db.is_some() {
        None
    } else {
        match read_cli_subjects(&args, &settings, &mut *diagnostics, archive.is_some()) {
            Ok(read) => Some(read),
            Err(error) => {
                // NCBI reads the `-subject` file inside `Run`, after the diagnostics handler
                // that saves the messages is set: an error archive follows the error.
                write_run_error_archive(&args, archive.as_ref(), &error);
                return Err(error);
            }
        }
    };
    if archive.is_some() {
        if let Some(read) = &cli_subjects {
            archive::check_titles(&read.subjects, "subject")?;
            if let Some(full) = &read.records {
                archive::check_titles(full, "subject")?;
            }
        }
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
    // ```c
    //     if (args.Exist(kArgQuery) && args[kArgQuery].HasValue() &&
    //         m_InputStream == NULL) {
    // ...
    //         else {
    //             m_InputStream = &args[kArgQuery].AsInputFile();
    //         }
    //     }
    // ...
    //     else {
    //         m_OutputStream = &args[kArgOutput].AsOutputFile();
    //     }
    // ```
    // NCBI reference: ncbi-blast/c++/src/corelib/ncbiargs.cpp:777-781
    // ```c
    //     if (AsString() == "-") {
    // #if defined(NCBI_OS_MSWIN)
    //         NcbiSys_setmode(NcbiSys_fileno(stdout), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
    // #endif
    //         m_Ios = &cout;
    // ```
    // The output file is created here, before the options are checked, so a run that
    // stops at a check leaves an empty file; `-` is standard output.
    check_utf8_file_name(&args.query, "query", "BLASTN")?;
    // NCBI opens the query file inside `Run` too (`CStdCmdLineArgs` in its handler of the
    // arguments does not): its error is saved and an error archive follows.
    let query_file = match open_input(&args.query, "query", "BLASTN") {
        Ok(file) => file,
        Err(error) => {
            write_run_error_archive(&args, archive.as_ref(), &error);
            return Err(error);
        }
    };
    // The output file is opened once, as NCBI's stream (a named pipe gives one reader one
    // end of file).
    let out_file = open_out_file(&args)?;
    let outfmt = args.outfmt.clone();
    let pairwise = output_formats.contains(&BlastnOutputFormat::Pairwise);
    let mut report = cli_report_stream(out_file);
    let result = {
        let mut outputs =
            ReportOutputs::single(&outfmt, OutputSink::Writer(&mut report), &mut *diagnostics);
        search_cli(
            args,
            query_file,
            cli_subjects.as_ref().map(|read| {
                (
                    read.subjects(),
                    SubjectRanges {
                        placements: &read.placements,
                        records: read.records.as_deref(),
                    },
                )
            }),
            &mut outputs,
            output_formats,
            settings,
        )
    };
    // NCBI writes an error archive after the archives of the batches before when the
    // application stops with one of its own errors, and when the run ends with messages that
    // no archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status:
    // blastn_app.cpp:331-334).
    if let Some(handle) = &archive {
        match &result {
            Err(error) => {
                if let Some(lines) = archive::error_lines(error, false) {
                    let _ =
                        crate::report::archive::write_error_archive(handle, &lines, &mut report);
                }
            }
            Ok(()) => {
                if handle.has_messages() {
                    let _ = crate::report::archive::write_error_archive(handle, &[], &mut report);
                }
            }
        }
    }
    // What was written before an error (such as the outfmt 0 prolog) stays in the file.
    let flushed = std::io::Write::flush(&mut report);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:118-119
    // ```c
    // {
    //     m_Outfile.exceptions(NcbiBadbit);
    // ```
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:252-255
    // ```c
    //     catch (const std::ios::failure&) {                                      \
    //         LOG_POST(Error << "BLAST failed to write output");                  \
    //         exit_code = BLAST_OUTPUT_ERROR;                                     \
    //     }                                                                       \
    // ```
    // The outfmt 0 formatter's stream throws when a write fails (oracle: `-out /dev/full`,
    // exit 6). For outfmt 6/7 NCBI aborts instead, and LOSAT reports the error
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES).
    if report.failed && pairwise {
        return Err(crate::cli::NativeError {
            exit: 6,
            message: "BLAST failed to write output\n".to_string(),
        }
        .into());
    }
    result?;
    flushed.context("failed to write the output")
}

/// The subjects of `-subject` as the handler of the database arguments reads them (`run`):
/// the searched records (their intervals with `-subject_loc`), where their letters lie in
/// the records, and the records as read (with `-subject_loc`).
/// The error archive of a run that stopped in `Run` with `error` before its search wrote
/// anything, when NCBI writes one (`PrintErrorArchive`): into `-out`, or nothing when that
/// cannot be opened (NCBI's call swallows that).
fn write_run_error_archive(
    args: &BlastnArgs,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    error: &anyhow::Error,
) {
    let Some(handle) = archive else {
        return;
    };
    let Some(lines) = archive::error_lines(error, true) else {
        return;
    };
    if let Ok(out_file) = open_out_file(args) {
        let mut report = cli_report_stream(out_file);
        let _ = crate::report::archive::write_error_archive(handle, &lines, &mut report);
        let _ = std::io::Write::flush(&mut report);
    }
}

/// The `-out` file, created (`-`: standard output).
fn open_out_file(args: &BlastnArgs) -> Result<Option<std::io::BufWriter<std::fs::File>>> {
    if let Some(path) = args.out.as_deref() {
        check_utf8_file_name(path, "out", "BLASTN")?;
    }
    Ok(
        match args.out.as_deref().filter(|path| path.as_os_str() != "-") {
            Some(path) => Some(std::io::BufWriter::new(
                std::fs::File::create(path).map_err(|_| crate::cli::inaccessible("out", path))?,
            )),
            None => None,
        },
    )
}

/// The report stream of the command line: the `-out` file, or standard output.
fn cli_report_stream(
    out_file: Option<std::io::BufWriter<std::fs::File>>,
) -> crate::cli::ReportStream {
    crate::cli::ReportStream {
        inner: match out_file {
            Some(file) => Box::new(file) as Box<dyn std::io::Write + Send>,
            None => crate::cli::report_standard_output(),
        },
        failed: false,
    }
}

struct CliSubjects {
    subjects: Vec<FastaRecord>,
    placements: Placements,
    records: Option<Vec<FastaRecord>>,
    /// EXPERIMENT (LOSAT_X_BNSTORE): the letters of `subjects`, packed (`x_read_subject_stores`;
    /// the records keep none).
    stores: Option<Vec<PackedNucleotide>>,
}

impl CliSubjects {
    /// The searched subjects.
    fn subjects(&self) -> Subjects<'_> {
        Subjects::Records(SubjectLetters::new(&self.subjects, self.stores.as_deref()))
    }
}

/// Reads the subjects of `-subject` where NCBI's handler of the database arguments reads
/// them (`run`; NCBI references there), with NCBI's error when `-subject` is missing (a
/// search without `-db`).
fn read_cli_subjects(
    args: &BlastnArgs,
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    diagnostics: &mut (dyn std::io::Write + Send),
    keep_deflines: bool,
) -> Result<CliSubjects> {
    let Some(subject_path) = args.subject.as_deref() else {
        return Err(crate::cli::NativeError {
            exit: 1,
            message: "BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified\nPlease refer to the BLAST+ user manual.\n".to_string(),
        }
        .into());
    };
    check_utf8_file_name(subject_path, "subject", "BLASTN")?;
    let subject_file = open_input(subject_path, "subject", "BLASTN")?;
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
    // ```c++
    //             subj_input_stream = &args[kArgSubject].AsInputFile();
    //         }
    //
    //         TSeqRange subj_range;
    //         if (args.Exist(kArgSubjectLocation) && args[kArgSubjectLocation]) {
    //             subj_range =
    //                 ParseSequenceRange(args[kArgSubjectLocation].AsString(),
    //                             "Invalid specification of subject location");
    //         }
    // ```
    // The subject range is read after the subject file is opened, before it is read.
    let subject_range =
        parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "BLASTN")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
    // ```c++
    //     CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
    //     CRef<CBlastInput> input(new CBlastInput(fasta));
    //     CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
    //     sequences = input->GetAllSeqs(*scope);
    //     return scope;
    // ```
    // NCBI reads the records one at a time, writes each one's messages as it reads it, and
    // checks each one's range after reading it: a range that starts past the end of a
    // record stops the reading there (`read_subjects`).
    let mut subject_source = FastaInputSource::from_argument(
        subject_path,
        subject_file,
        ReaderConfig::subject("BLASTN", false, settings.data_loaders).with_defline(keep_deflines),
    );
    // EXPERIMENT (LOSAT_X_BNSTORE): without a subject range, each record's letters are packed as
    // soon as it has been read, and the records keep none (`x_read_subject_stores`).
    // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs); with
    // LOSAT_X_BNSTORE the packing of each record is in this interval too. No NCBI counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let (read_subjects, stores) = {
        let mut warn = |message: &[u8]| std::io::Write::write_all(&mut *diagnostics, message);
        // The archive writes the subjects' letters from the records, not from a packed store.
        if x_bnstore() && subject_range.is_none() && !keep_deflines {
            x_read_subject_stores(&mut subject_source, &mut warn)?
        } else {
            (
                read_subjects(&mut subject_source, subject_range.as_ref(), &mut warn)?,
                None,
            )
        }
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    if x_bnstore() {
        crate::utils::xreach::x_reach!("LOSAT_X_BNSTORE", stores.is_some());
    }
    drop(subject_source);
    // NCBI reads records without residues (and intervals without letters) without a message
    // and reports them when it sets up the search (`search`).
    let ranged_subjects = cut_subjects(&read_subjects, subject_range.as_ref())?;
    let (subjects, subject_placements, full_subjects) = match ranged_subjects {
        Some((cut, placements)) => (cut, placements, Some(read_subjects)),
        None => (read_subjects, Placements::default(), None),
    };
    check_subjects_not_empty(&subjects)?;
    Ok(CliSubjects {
        subjects,
        placements: subject_placements,
        records: full_subjects,
        stores,
    })
}

/// The part of `run` after the output is opened: NCBI's filtering handler (`-dust`,
/// `-window_masker_db`), the
/// processing of the options, the opening of the database of `-db`, `Query is Empty!`, the
/// reading of the queries and the search. `subjects` are those of `-subject` (`None` with
/// `-db`).
fn search_cli(
    mut args: BlastnArgs,
    query_file: std::fs::File,
    subjects: Option<(Subjects<'_>, SubjectRanges<'_>)>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    args.resolve_dust()?;
    args.resolve_window_masker()?;
    let query_range = parse_query_range(&args)?;
    process_options(&mut args, outputs)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:191-197
    // ```c++
    //         /*** Initialize the database/subject ***/
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope;
    //         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    //         _ASSERT(db_adapter && scope);
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:213-219
    // ```c++
    //         try {
    //             // Try to open the BLAST database even for remote searches, as if
    //             // it is available locally, it will be better to fetch the
    //             // sequence data for formatting from this (local) source
    //             CRef<CSeqDB> seqdb = search_db->GetSeqDb();
    //             db_adapter.Reset(new CLocalDbAdapter(*search_db));
    //             dl = RegisterOMDataLoader(seqdb);
    // ```
    // The database of `-db` is opened here: after the option handlers (the query file is
    // open and the output file exists) and before the test for an empty query.
    // The taxid or seqid list (line DI-8) is read first: NCBI's handler of the database
    // arguments reads it (`s_GetTaxIDList`, `CSeqDBFileGiList`) before the database opens.
    let database = match args.database_name() {
        Some(name) => {
            let list = args.db_filter.resolve(&settings, outputs.diagnostics)?;
            let database = database::BlastnDatabase::open_restricted(
                name,
                &settings,
                list.as_ref(),
                outputs.diagnostics,
            )?
            .showing_gis(args.show_gis);
            crate::algorithm::common::database::note_archive_database_filter(
                args.archive.as_ref(),
                &args.db_filter,
                list.as_ref(),
                database.filtered_empty(),
                args.searchsp.is_some_and(|space| space > 0),
                args.dbsize.is_some_and(|size| size > 0),
                "BLASTN",
            )?;
            Some(database)
        }
        None => None,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104
    // ```c++
    // void CBlastTabularInfo::x_CheckTaxDB()
    // {
    //    	if( x_IsFieldRequested(eSubjectSciNames) ||
    //     ...
    //    		string resolved = SeqDB_ResolveDbPath("taxdb.bti");
    // ```
    // A tabular field list that shows a name of the taxonomy database asks whether `taxdb.bti`
    // is on the search path: its formatters warn without it (the plural names), and LOSAT
    // refuses a name of a subject with a taxid with it (`algorithm::blastn::tabular`).
    if output_formats
        .iter()
        .filter_map(|format| format.tabular_spec())
        .any(|spec| spec.needs().taxonomy_names)
    {
        args.tabular_taxdb = Some(match database.as_ref() {
            Some(database) => database.taxdb_found(),
            None => super::super::tabular::subject_taxdb_found(&settings)?,
        });
    }
    // LOSAT's refusals of a database for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let Some(database) = database.as_ref() {
        for serial in output_formats.iter().filter_map(|format| format.serial()) {
            crate::report::serial::seqalign::check_database(
                serial,
                database.name(),
                &database.pairwise_database().definition,
                Program::Blastn,
            )?;
        }
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:201-214
    // ```c++
    //         SDataLoaderConfig dlconfig =
    //             InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
    //                                                    db_adapter);
    //         CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
    //                                      query_opts->UseLowercaseMasks(),
    //                                      query_opts->GetParseDeflines(),
    //                                      query_opts->GetRange());
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())) {
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    //         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    //         CBlastInput input(&fasta);
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
    // ```c++
    // 	char c;
    // 	CNcbiStreampos orig_p = in.tellg();
    // 	// Piped input
    // 	if(orig_p < 0)
    // 		return false;
    // ```
    // The query's data loaders are those of the subjects (`SDataLoaderConfig`; the BLAST
    // database of a database search is added only to the loaders' databases). The position
    // is taken on the opened file, before it is read (`stream_is_empty`). When the subjects
    // were read from standard input too, `cin` has reached its end (a failed stream), so
    // NCBI gets no position for the query either, and its reader is at its end: no query
    // batch is read.
    let shared_standard_input =
        args.query.as_os_str() == "-" && args.subject_path().as_os_str() == "-";
    // EXPERIMENT (LOSAT_X_QSTREAM): what its rule reads (`x_query_stream::Rule`): the task and the
    // largest query batch. EXPERIMENT (LOSAT_X_BNQPACK) reads the query file's size and that batch
    // too (`x_query_drop::packs`).
    let x_file_bytes = x_query_file_bytes(&args.query);
    let x_batch_letters = x_query_batch_letters(&args.task);
    let x_rule = x_query_stream::Rule {
        megablast: args.task == "megablast",
        batch_letters: x_batch_letters,
    };
    // EXPERIMENT (LOSAT_X_QSTREAM): whether the query file is read a second time, batch by batch
    // (a second handle on it, at the position its reader starts from), where the reports may be
    // written batch by batch (LOSAT_X_RPTSTREAM).
    let x_opening = x_query_stream::Opening::new(
        &args.query,
        &query_file,
        shared_standard_input,
        x_rpt_stream::switch_on() && !args.parse_deflines,
        x_rule,
    );
    // EXPERIMENT (LOSAT_X_BNQDROP): whether the records read before the search move into a store
    // that the batch loop takes them out of (`x_query_drop::Rule`), decided before the reading so
    // that EXPERIMENT (LOSAT_X_BNQPACK) can pack their letters as they are read.
    let x_drop_rule = !x_opening.acts()
        && x_query_drop::on()
        && x_query_drop::Rule {
            formats: &output_formats,
            query_range: query_range.is_some(),
            hit_records: outputs.hits.is_some(),
        }
        .acts();
    let x_pack = x_query_drop::packs(x_drop_rule, x_file_bytes, x_batch_letters);
    if x_drop_rule && x_query_drop::pack_on() {
        crate::utils::xreach::x_reach!("LOSAT_X_BNQPACK", x_pack);
    }
    crate::blastinput::fasta_reader::believed::refuse_with_subject(
        args.parse_deflines,
        args.subject.is_some(),
        "BLASTN",
    )?;
    let query_config = ReaderConfig::query("BLASTN", false, settings.data_loaders)
        .with_defline(args.archive.is_some())
        .with_parse_deflines(args.parse_deflines)
        .with_strand(&args.strand);
    let mut query_source = FastaInputSource::from_argument(&args.query, query_file, query_config);
    if !shared_standard_input && query_source.stream_is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [blastn] Query is Empty!\n")?;
        // NCBI returns from `Run` here, before its `PrintErrorArchive`: the message is not an
        // archive's.
        if let Some(handle) = &args.archive {
            handle.take_messages();
        }
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
    // ```c++
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
    // 	    BLAST_PROF_START( APP.LOOP.PRE );
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
    // The records are read here, and their messages are written batch by batch
    // (`run_in_pool`). What ends the reading (the end of the input, an `eEOF`, a reader
    // error or LOSAT's rejection of a Seq-id line) comes in the batch that NCBI reads it in,
    // after the reports of the batches before (`run_in_pool`).
    // Diagnostic (LOSAT_TIMING): the query FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    // EXPERIMENT (LOSAT_X_QSTREAM): with the switch acting this reading keeps only each record's
    // length and how the reading ends; the batches read their records again (`x_query_stream`).
    let mut x_packed = Vec::new();
    let (records, end, x_store) = if x_opening.acts() {
        let (store, end) = x_query_stream::Store::read(
            x_opening,
            &mut query_source,
            query_config,
            query_range.as_ref(),
            args.num_threads,
        )?;
        (Vec::new(), end, store)
    } else if x_pack {
        let (QueryRecords { records, end }, packed) =
            crate::blastinput::fasta_reader::x_read_queries_packed(&mut query_source);
        x_packed = packed;
        (records, end, x_query_stream::Store::off(x_opening))
    } else {
        let QueryRecords { records, end } = read_queries(&mut query_source);
        (records, end, x_query_stream::Store::off(x_opening))
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(query_source);
    if args.archive.is_some() {
        archive::check_titles(&records, "query")?;
    }
    let reading = match end {
        QueryEnd::Input => QueryReading::End,
        QueryEnd::EofError => QueryReading::BlankLines,
        QueryEnd::Error { error, warnings } => QueryReading::Error { error, warnings },
    };
    // A database has no subject ranges (`-subject_loc` excludes `-db`).
    let no_placements = Placements::default();
    let (subjects, subject_ranges) = match (subjects, database.as_ref()) {
        (Some((subjects, ranges)), _) => (subjects, ranges),
        (None, Some(database)) => (
            Subjects::Database(database),
            SubjectRanges {
                placements: &no_placements,
                records: None,
            },
        ),
        (None, None) => unreachable!("`run` reads the subjects of a search without -db"),
    };
    // EXPERIMENT (LOSAT_X_BNQDROP): where no report reads the queries' letters, the records read
    // before the search move into a store that the batch loop takes them out of, so a batch can drop
    // its letters once it has encoded them (`x_query_drop.rs`).
    let x_taken = !x_store.acts() && x_query_drop::on() && {
        crate::utils::xreach::x_reach!("LOSAT_X_BNQDROP", x_drop_rule);
        x_drop_rule
    };
    let x_taken_store;
    let queries = if x_store.acts() {
        SearchQueries::Streamed(&x_store)
    } else if x_taken {
        x_taken_store = x_query_drop::Store::new(records, x_packed);
        SearchQueries::Taken(&x_taken_store)
    } else {
        SearchQueries::Records(&records)
    };
    search(
        args,
        queries,
        reading,
        Ranges {
            query: query_range,
            subjects: subject_ranges,
        },
        subjects,
        outputs,
        output_formats,
    )
}

/// How the reading of the queries of a search ended, after its records.
#[derive(Debug)]
pub(crate) enum QueryReading {
    /// The records of `run_local`: an input without records is NCBI's `Query is Empty!`
    /// (the inputs of the ABI are read as regular files).
    Records,
    /// The CLI's reader reached the end of its input (`CBlastInput::End()`); without
    /// records, no query batch is read.
    End,
    /// The CLI's reader found only blank and comment lines where it expected a record
    /// (`eEOF`), which ends the batch being read.
    BlankLines,
    /// The CLI's reader failed on the record after the last one read (a reader error, or
    /// LOSAT's explicit rejection of a line that NCBI would fetch as a Seq-id), after
    /// writing `warnings`.
    Error { error: ReadError, warnings: Vec<u8> },
}

impl QueryReading {
    /// Whether NCBI's reader reads on after the last record (`End()` is false): an `eEOF`
    /// or an error then comes in the batch being read.
    pub(crate) fn reads_past_records(&self) -> bool {
        matches!(self, QueryReading::BlankLines | QueryReading::Error { .. })
    }
}

/// The query range, read by NCBI's query options handler: after the filtering handler
/// (`-dust`), before the formatting handler and the check of the options.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blastn_args.cpp:100-102
/// ```c++
///     m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
///     arg.Reset(m_QueryOptsArgs);
///     m_Args.push_back(arg);
/// ```
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
/// ```c++
///     // set the sequence range
///     if (args.Exist(kArgQueryLocation) && args[kArgQueryLocation]) {
///         m_Range = ParseSequenceRange(args[kArgQueryLocation].AsString(),
///                                      "Invalid specification of query location");
///     }
/// ```
fn parse_query_range(args: &BlastnArgs) -> Result<Option<SequenceRange>> {
    parse_optional_range(args.query_loc.as_deref(), RangeRole::Query, "BLASTN")
}

/// The ranges of a search: the query range as given, and the subjects' ranges.
struct Ranges<'a> {
    query: Option<SequenceRange>,
    subjects: SubjectRanges<'a>,
}

/// Where each subject's searched letters lie in its record, and the records as read
/// (with `-subject_loc`; the report shows their lowercase letters).
#[derive(Clone, Copy)]
struct SubjectRanges<'a> {
    placements: &'a Placements,
    records: Option<&'a [FastaRecord]>,
}

/// The ranges as the report uses them: the query input and the subjects' ranges.
struct ReportRanges<'a> {
    queries: &'a QueryInput,
    subjects: SubjectRanges<'a>,
}

/// NCBI's error for a subject file without records, raised when it reads the subjects
/// (the first option handler, `run`).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
/// ```c
/// CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
///     : m_QueryVector(& queries)
/// {
///     if (queries.Empty()) {
///         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:225-227
/// ```c
///             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
///             exit_code = BLAST_ENGINE_ERROR;                                 \
/// ```
pub(crate) fn check_subjects_not_empty<R>(subjects: &[R]) -> Result<()> {
    if subjects.is_empty() {
        return Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
        }
        .into());
    }
    Ok(())
}

/// The report of the CLI's `-outfmt`, read as NCBI reads it before its option handlers run
/// (`run`; `blastinput/app.rs` `parse_formatting_string`, `report_format`). A format that
/// NCBI itself rejects later has no report: the run stops at NCBI's error, the formatting
/// handler's in `process_options` (AIRR rearrangement, FASTA) or the XInclude formats' in
/// `search` (13 and 14 to standard output).
fn cli_output_formats(
    outfmt: &str,
    out: Option<&std::path::Path>,
) -> Result<Vec<BlastnOutputFormat>> {
    let choice = app::parse_formatting_string(outfmt)?;
    app::report_format(&choice, Program::Blastn, out)?
        .map(|format| BlastnOutputFormat::of(format, &choice))
        .into_iter()
        .collect()
}

/// The requested output formats of `run_local`, each read as the CLI reads its single
/// `-outfmt` (`cli_output_formats`); a format that NCBI rejects later fails here with NCBI's
/// error, as in the other programs' `run_local` (they have no output file).
fn parse_output_formats<'a>(
    formats: impl IntoIterator<Item = &'a str>,
) -> Result<Vec<BlastnOutputFormat>> {
    formats
        .into_iter()
        .map(|outfmt| {
            let choice = app::parse_formatting_string(outfmt)?;
            match app::report_format(&choice, Program::Blastn, None)? {
                Some(format) => BlastnOutputFormat::of(format, &choice),
                None => {
                    app::formatting_handler_check(&choice, Program::Blastn)?;
                    app::xinclude_check(&choice)?;
                    anyhow::bail!(
                        "output format {} is not supported by LOSAT's BLASTN",
                        choice.number
                    )
                }
            }
        })
        .collect()
}

/// NCBI's processing of the options after the query and the output are opened
/// (`SetOptions`): the formatting handler (its errors for the formats of other programs,
/// then its few-matches warning), and the check of the options.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886
/// ```c
///     ParseFormattingString(args, m_OutputFormat, m_CustomOutputFormatSpec,m_CustomDelim);
///     if((m_OutputFormat == eSAM) && !(m_FormatFlags & eIsSAM) ){
/// ...
///     if((m_OutputFormat == eAirrRearrangement) && !(m_FormatFlags & eIsAirrRearrangement) ){
///         NCBI_THROW(CInputException, eInvalidInput,
///                    "AIRR rearrangement format is only applicable to igblastn" );
///     }
///     if (m_OutputFormat == eFasta) {
///         NCBI_THROW(CInputException, eInvalidInput,
///                    "FASTA output format is only applicable to magicblast");
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
/// ```c
///     if(hitlist_size < 5){
///    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3636
/// ```c
///     try { retval->Validate(); }
/// ```
/// The hit list size is the -max_target_seqs value, or 500 when it is omitted.
fn process_options(args: &mut BlastnArgs, outputs: &mut ReportOutputs<'_>) -> Result<()> {
    // The formatting handler comes after the query options handler (blastn_args.cpp:100-105).
    let mut outfmts = Vec::with_capacity(outputs.formats.len());
    for format in &outputs.formats {
        let choice = app::parse_formatting_string(format.outfmt)?;
        app::formatting_handler_check(&choice, Program::Blastn)?;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:1397-1440
        // (`BLASTUngappedHspListToSeqAlign`: a subject's Seq-align of Dense-diags, and
        // objtools/align_format/sam.cpp:66-80 for SAM's one line per subject)
        // O4a-1 writes the formats that need no Seq-align of Dense-diags: 8, 9, 11, 12 (the
        // Seq-align formats) and 17 (SAM) of an ungapped search are refused (O4a-2).
        if args.ungapped && matches!(choice.number, 8 | 9 | 11 | 12 | 17) {
            anyhow::bail!(
                "output format {} with -ungapped (an ungapped search) is not supported by LOSAT's BLASTN",
                choice.number
            );
        }
        outfmts.push(choice.number);
    }
    // The hit list size and the counts of the pairwise report, with the formatting
    // handler's warnings (`ExtractAlgorithmOptions`), once for the run.
    let (counts, warnings) = crate::algorithm::blastn::scoring::report_counts_of(args, &outfmts)?;
    args.report_counts = counts;
    outputs.diagnostics.write_all(&warnings)?;
    check_scoring_options(args)
}

// NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-467
// ```c
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// ```
// NCBI formats one result set without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// Runs one BLASTN search over already parsed records and writes every requested
/// output format from the same result (the shared entry of the CLI, web ABI v1 and
/// v2).
///
/// Each requested `-outfmt` is parsed as on the command line, before the search starts.
/// No BLASTN search option depends on the output formats it supports (6 and 7). The
/// `-query` and `-subject` values of `args` are used only as display names. The hit
/// records (`ReportOutputs::hits`) are not produced yet.
///
/// The BLAST archive (outfmt 11) writes each record's DefLine as the reader kept it
/// (`ReaderConfig::with_defline(true)`). Records without a kept defline get `>` and their
/// title, which differs from NCBI for a defline with extra spaces or control bytes.
pub fn run_local(
    mut args: BlastnArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The BLAST archive (outfmt 11) among the formats: its request, the refusals of what
    // LOSAT cannot write, and a diagnostics stream that saves the messages.
    let archive_requested = outputs.formats.iter().any(|format| {
        app::parse_formatting_string(format.outfmt)
            .is_ok_and(|choice| choice.number == app::output_format::ARCHIVE_FORMAT)
    });
    if archive_requested && args.db.is_none() {
        QueryStrand::parse(&args.strand)?;
        args.archive = archive::handle_for(&args, args.dust_filtering.as_deref(), true)?;
        archive::check_titles(query_records, "query")?;
        archive::check_titles(subject_records, "subject")?;
        if let Some(handle) = args.archive.clone() {
            return outputs.with_archive_diagnostics(handle, |view| {
                run_local_formats(args, query_records, subject_records, view)
            });
        }
    }
    run_local_formats(args, query_records, subject_records, outputs)
}

fn run_local_formats(
    mut args: BlastnArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The callers of `run_local` (LOSAT Web's adapter, ABI v1) give subject records; a
    // database is read from files by the command line only.
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported by LOSAT Web yet");
    }
    // The counts file of `-window_masker_db` comes from the caller's counts source (LOSAT
    // Web's registered inputs); a caller that gives none has no file to read.
    if args.window_masker_db.is_some() && args.window_masker_counts.is_none() {
        anyhow::bail!(
            "-window_masker_db (a WindowMasker counts file) is not supported by LOSAT Web yet"
        );
    }
    // The command line has checked `-strand` at parsing (clap); a caller of `run_local` has not.
    QueryStrand::parse(&args.strand)?;
    let output_formats = parse_output_formats(outputs.formats.iter().map(|format| format.outfmt))?;
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range =
        parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "BLASTN")?;
    for record in &subject_records[..subjects_read(subject_records, subject_range.as_ref())] {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = cut_subjects(subject_records, subject_range.as_ref())?;
    let (searched_subjects, subject_placements) = match &ranged_subjects {
        Some((records, placements)) => (records.as_slice(), placements.clone()),
        None => (subject_records, Placements::default()),
    };
    check_subjects_not_empty(searched_subjects)?;
    args.resolve_dust()?;
    args.resolve_window_masker()?;
    let query_range = parse_query_range(&args)?;
    process_options(&mut args, outputs)?;
    search(
        args,
        SearchQueries::Records(query_records),
        QueryReading::Records,
        Ranges {
            query: query_range,
            subjects: SubjectRanges {
                placements: &subject_placements,
                records: ranged_subjects.is_some().then_some(subject_records),
            },
        },
        Subjects::Records(SubjectLetters::new(searched_subjects, None)),
        outputs,
        output_formats,
    )
}

/// The queries of a search: records read before it, or EXPERIMENT (LOSAT_X_QSTREAM) a query file
/// whose batches read their records when the batch loop needs them (`x_query_stream`).
#[derive(Clone, Copy)]
enum SearchQueries<'a> {
    Records(&'a [FastaRecord]),
    Streamed(&'a x_query_stream::Store),
    /// EXPERIMENT (LOSAT_X_BNQDROP): records read before the search that the batch loop takes out
    /// of their store (`x_query_drop.rs`).
    Taken(&'a x_query_drop::Store),
}

fn search(
    args: BlastnArgs,
    queries: SearchQueries<'_>,
    reading: QueryReading,
    ranges: Ranges<'_>,
    subjects: Subjects<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
) -> Result<()> {
    // EXPERIMENT (LOSAT_X_QSTREAM): a streamed search has no records here (its reading is never
    // `Records`), and its query input comes from the lengths of the first reading; EXPERIMENT
    // (LOSAT_X_BNQDROP) a taken one has its records in its store, with its query input.
    let query_records: &[FastaRecord] = match queries {
        SearchQueries::Records(records) => records,
        SearchQueries::Streamed(_) | SearchQueries::Taken(_) => &[],
    };
    // NCBI reads the subjects before the queries (`run`); the records are those of NCBI's
    // reader (`fasta_reader`; `run_local`'s callers read them so). The CLI has already
    // reported an empty query file and read its queries; the records of the other callers
    // get NCBI's warning for an input without records.
    if query_records.is_empty() && matches!(reading, QueryReading::Records) {
        outputs
            .diagnostics
            .write_all(b"Warning: [blastn] Query is Empty!\n")?;
        return Ok(());
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:129-136
    // ```c
    //     if(m_IsDbScan) {
    // 	int num_seqs=0;
    //         int total_length=0;
    // 	if (!is_remote_search)
    //         {
    //                 BlastSeqSrc* seqsrc = db_adapter.MakeSeqSrc();
    //                 num_seqs=BlastSeqSrcGetNumSeqs(seqsrc);
    //                 total_length=static_cast<int>(BlastSeqSrcGetTotLen(seqsrc));
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
    // ```c
    //         catch(CBlastException & e ) {
    //         	// Skip bad subject sequence
    //         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
    //         		seqblk_vec->push_back(subj);
    //         ...
    //         		warning += "Subject sequence contains no data";
    //         		ERR_POST(Warning << warning);
    //         		continue;
    // ```
    // The formatter, made after `Query is Empty!` and before the outfmt 0 prolog, sets up
    // the subjects: a subject without letters (a record without residues, or an interval
    // that starts just past its record's end) gets its warning, stays in the database
    // statistics with no letters, and never hits. (A database is not set up so.)
    if let Subjects::Records(subject_set) = subjects {
        write_empty_subject_warnings_by_length(
            subject_set.records(),
            |index| subject_set.len(index),
            "blastn",
            outputs.diagnostics,
        )?;
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:255-258
    // ```c
    //         if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    //         }
    //         formatter.PrintProlog();
    // ```
    // The XInclude formats written to standard output fail after the formatter is made (a
    // format 13 or 14 to a file takes the file name as the base of its side files).
    for format in &outputs.formats {
        app::xinclude_check_for(
            &app::parse_formatting_string(format.outfmt)?,
            args.out.as_deref(),
        )?;
    }
    // NCBI reads the queries after `Query is Empty!`, one batch at a time, with its
    // reader's warnings (`run_in_pool`).
    // NCBI's x_CleanAndCompress reads past the end of some titles of punctuation (NCBI
    // crashes when it writes the title of such a subject with hits); LOSAT stops at the end
    // of the title (`report/defline.rs`, approved exception 2 of PD-LOSAT-NCBI-DEFECTS).
    // LOSAT's limits come where NCBI starts the search, after its checks and its
    // `Query is Empty!` success.
    check_losat_limits(&args)?;
    check_unsupported_environment(matches!(subjects, Subjects::Database(_)))?;
    let batching = query_batching_from_environment(
        super::super::coordination::task_uses_megablast_chunks(&args.task),
    )?;
    // The query range applies to every query record as NCBI's batch reader reads it
    // (`cut_queries`): a record is searched cut, and one whose interval starts just past
    // its end is a query without data (`run_in_pool`).
    let ranged_queries = ranges
        .query
        .as_ref()
        .filter(|_| matches!(queries, SearchQueries::Records(_)))
        .map(|range| cut_queries(query_records, range));
    let input_records = query_records;
    let (query_records, query_input) = match (ranged_queries, queries) {
        (Some(ranged), _) => (
            std::borrow::Cow::Owned(ranged.records),
            std::borrow::Cow::Owned(ranged.input),
        ),
        (None, SearchQueries::Records(_)) => (
            std::borrow::Cow::Borrowed(query_records),
            std::borrow::Cow::Owned(QueryInput::whole(query_records)),
        ),
        (None, SearchQueries::Streamed(store)) => (
            std::borrow::Cow::Borrowed(query_records),
            std::borrow::Cow::Borrowed(store.query_input()),
        ),
        (None, SearchQueries::Taken(store)) => (
            std::borrow::Cow::Borrowed(query_records),
            std::borrow::Cow::Borrowed(store.query_input()),
        ),
    };
    let queries = match queries {
        SearchQueries::Records(_) => x_query_stream::Queries::Records {
            input: input_records,
            searched: &query_records,
        },
        SearchQueries::Streamed(store) => x_query_stream::Queries::Streamed(store),
        SearchQueries::Taken(store) => x_query_stream::Queries::Taken(store),
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:129-138
    // ```c
    // 	int num_seqs=0;
    //         int total_length=0;
    // 	if (!is_remote_search)
    //         {
    //                 BlastSeqSrc* seqsrc = db_adapter.MakeSeqSrc();
    //                 num_seqs=BlastSeqSrcGetNumSeqs(seqsrc);
    //                 total_length=static_cast<int>(BlastSeqSrcGetTotLen(seqsrc));
    //         }
    // ```
    // The report and the query batches (`GetDbTotalLength`) take the total through an
    // `int`, which LOSAT does not reproduce beyond its range. (A database's total is the
    // `Int8` of `GetBlastDbInfo`, blast_format.cpp:140-154.)
    if let Subjects::Records(subject_set) = subjects {
        let subject_letters: usize = (0..subject_set.records().len())
            .map(|index| subject_set.len(index))
            .sum();
        if subject_letters > i32::MAX as usize {
            anyhow::bail!(
                "the subjects have {subject_letters} letters; NCBI BLAST+ reports a total of 2^31 letters or more through a 32-bit int, which is not supported by LOSAT's BLASTN"
            );
        }
    }
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    let x_result = crate::utils::threading::with_search_pool(args.num_threads, "blastn", |pool| {
        run_in_pool(
            args,
            queries,
            reading,
            &query_input,
            subjects,
            ranges.subjects,
            outputs,
            output_formats,
            pool,
            batching,
        )
    });
    // No NCBI counterpart: the PAIRPAR counters with LOSAT_TIMING; nothing without it.
    x_pair_seed::summary();
    x_result
}

/// Rejects the environment variables that put NCBI BLAST+ into a search mode or report that
/// LOSAT does not reproduce.
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:204-210
/// ```c
///     if ( (subjects = db_args->GetSubjects(scope)) ) {
///         _ASSERT(search_db.Empty());
/// 	char* bl2seq_legacy = getenv("BL2SEQ_LEGACY");
/// 	if (bl2seq_legacy)
///         	db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, false));
/// 	else
///         	db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, true));
/// ```
/// `BL2SEQ_LEGACY` (any value, also an empty one) turns off the scan of the subjects as a
/// database: `CLocalBlast` then searches each subject on its own (local_blast.cpp:189,289)
/// and `CBlastFormat` writes the legacy bl2seq report (`m_IsBl2Seq && !m_IsDbScan`). NCBI
/// reads it in that branch only: a database search (`database`) does not read it.
fn check_unsupported_environment(database: bool) -> Result<()> {
    if !database && std::env::var_os("BL2SEQ_LEGACY").is_some() {
        anyhow::bail!(
            "the environment variable BL2SEQ_LEGACY, which makes NCBI BLAST+ search each subject on its own and write its legacy bl2seq report, is not supported by LOSAT's BLASTN"
        );
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:732-737
    // ```c
    // 		char * pre_fetch_limit_str = getenv("PRE_FETCH_SEQS_LIMIT");
    // 		if (pre_fetch_limit_str) {
    // 			int pre_fetch_limit = NStr::StringToInt(pre_fetch_limit_str);
    // 			if(pre_fetch_limit == 0) {
    // 				return false;
    // 			}
    // ```
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:761-766
    // ```c
    //     if (results.size() == 0) {
    //         return;
    //     }
    //     if(!s_PreFetchSeqs(results, format_type)){
    //     	return;
    //     }
    // ```
    // An integer only decides whether the report's sequences are fetched ahead, which
    // changes no output. Before the report of each query batch NCBI converts the value, and
    // one that `NStr::StringToInt` cannot convert (an empty one too) stops it with a
    // `CStringException` whose text names its build's source files (exit 255).
    if let Some(value) = std::env::var_os("PRE_FETCH_SEQS_LIMIT") {
        if crate::blastinput::query_batch::ncbi_string_to_int(&value).is_none() {
            anyhow::bail!(
                "the environment variable PRE_FETCH_SEQS_LIMIT has the value {:?}, which NCBI BLAST+ cannot convert to an int (it stops with a CStringException that names its build's source files); this is not supported by LOSAT's BLASTN",
                value.to_string_lossy()
            );
        }
    }
    Ok(())
}

/// An error of a query batch after the first that does not depend on the batch's queries
/// (the Karlin-Altschul error of the scoring options, or the average subject length;
/// `search_query_batch`): the first batch with a valid query meets it, and NCBI raises it
/// after the reports of the batches before, whose queries were all invalid.
#[derive(Debug)]
struct ErrorAfterInvalidBatches(anyhow::Error);

impl std::fmt::Display for ErrorAfterInvalidBatches {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for ErrorAfterInvalidBatches {}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:295-312
// ```c++
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
// 		        lcl_blast.SetBatchNumber( batch_num );
//                 results = lcl_blast.Run();
//     ...
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
// ```
/// A database error of a query batch's search (`database::search_blocks`,
/// `database::traceback_listed`) as NCBI reports it: a sequence file that the preliminary
/// search cannot read is NCBI's `BLAST engine error` (`database::PreliminaryReadError`,
/// local_blast.cpp:240-251), one that the traceback cannot read its `BLAST Database
/// error`. NCBI raises it in the batch's search, after the outfmt 0 prolog and the reports
/// of the batches before (`ErrorAfterInvalidBatches`; for the first batch, the prolog of a
/// file sink is written here, the others were written before the search). LOSAT's
/// refusals are errors of their own.
fn database_search_error(
    error: anyhow::Error,
    first_batch: bool,
    outputs: BatchOutputs<'_, '_>,
    megablast: bool,
    shown: &ShownSubjects,
) -> anyhow::Error {
    let error = match error.downcast::<database::PreliminaryReadError>() {
        Ok(read_error) => crate::blastinput::app::engine_error(read_error.0),
        Err(error) => error,
    };
    if !error.is::<crate::cli::NativeError>() {
        return error;
    }
    if !first_batch {
        return ErrorAfterInvalidBatches(error).into();
    }
    let (outputs, output_formats) = outputs.expect(FIRST_BATCH_OUTPUTS);
    match write_pairwise_prologs(outputs, output_formats, megablast, shown, false) {
        Ok(()) => error,
        Err(write_error) => write_error,
    }
}

/// NCBI's query batch size (`GetQueryBatchSize`, 0 for the adaptive `CBatchSizeMixer`) and
/// query chunk sizes, as a query batch of the search uses them.
#[derive(Clone, Copy, Debug)]
struct QueryBatching {
    /// The `BATCH_SIZE` value; 0 when it is not set.
    fixed_batch_size: i32,
    /// The size of the current batch (`TSeqPos`).
    batch_size: u32,
    split: SplitSizes,
}

/// Reads the environment variables that set NCBI's query batches and chunks, in NCBI's
/// order: `CHUNK_SIZE` (the `CBatchSizeMixer` maximum), `BATCH_SIZE`, then
/// `OVERLAP_CHUNK_SIZE` (the query splitter of the first batch).
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:261-262
/// ```c
///         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
///         int batch_size = m_CmdLineArgs->GetQueryBatchSize();
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:55-59
/// ```c
///     m_ChunkSize = SplitQuery_GetChunkSize(m_Options->GetProgram());
///     m_LocalQueryData = m_QueryFactory->MakeLocalQueryData(m_Options);
///     m_TotalQueryLength = m_LocalQueryData->GetSumOfSequenceLengths();
///     m_NumChunks = SplitQuery_CalculateNumChunks(m_Options->GetProgramType(),
///         &m_ChunkSize, m_TotalQueryLength, m_LocalQueryData->GetNumQueries());
/// ```
/// A value that `NStr::StringToInt` cannot convert makes NCBI stop with its
/// `CStringException`, whose text names the source path of the oracle's build (exit 255);
/// LOSAT rejects it. A negative `CHUNK_SIZE` is a `size_t` near 2^64; LOSAT rejects it where
/// it splits a query batch (`search_query_chunks`), and otherwise searches as NCBI.
fn query_batching_from_environment(megablast: bool) -> Result<QueryBatching> {
    let read = |variable: &str| -> Result<Option<i32>> {
        match std::env::var_os(variable) {
            Some(value) if !is_blank(value.as_encoded_bytes()) => {
                match crate::blastinput::query_batch::ncbi_string_to_int(&value) {
                    Some(number) => Ok(Some(number)),
                    None => anyhow::bail!(
                        "the environment variable {variable} has the value {:?}, which NCBI BLAST+ cannot convert to an int (it stops with a CStringException that names its build's source files); this is not supported by LOSAT's BLASTN",
                        value.to_string_lossy()
                    ),
                }
            }
            _ => Ok(None),
        }
    };
    let chunk_size = read("CHUNK_SIZE")?;
    let fixed_batch_size = crate::blastinput::query_batch::get_query_batch_size(
        std::env::var_os("BATCH_SIZE").as_deref(),
    )
    .map_err(|value| {
        anyhow::anyhow!(
            "the environment variable BATCH_SIZE has the value {value:?}, which NCBI BLAST+ cannot convert to an int (it stops with a CStringException that names its build's source files); this is not supported by LOSAT's BLASTN"
        )
    })?;
    let overlap = read("OVERLAP_CHUNK_SIZE")?;
    let split = SplitSizes::new(megablast, chunk_size, overlap);
    Ok(QueryBatching {
        fixed_batch_size,
        batch_size: fixed_batch_size as u32,
        split,
    })
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads());
// (*thread)->Run(); (*thread)->Join(&result);
/// The subjects of a search, read and encoded once for all query batches, as NCBI's
/// `CLocalDbAdapter` of the subject set serves every batch. Only the packed (ncbi2na)
/// letters that the preliminary stage reads are kept for every subject; the traceback
/// makes the blastna letters of a subject with preliminary HSPs from its record.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
/// ```c
/// if (subj_is_na) {
///     BlastSeqBlkSetSequence(subj, sequence.data.release(),
///        ((sentinels == eSentinels) ? sequence.length - 2 :
///         sequence.length));
///     ...
///     SBlastSequence compressed_seq =
///         subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
///                                   eNa_strand_plus, eNoSentinels);
///     BlastSeqBlkSetCompressedSequence(subj,
///                               compressed_seq.data.release());
/// ```
struct PreparedSubjects<'a> {
    records: &'a [FastaRecord],
    /// EXPERIMENT (LOSAT_X_BNSTORE): the letters of the records, packed (the records keep
    /// none); the preliminary stage reads their packed bytes.
    stores: Option<&'a [PackedNucleotide]>,
    metadata: SubjectMetadata,
    /// The packed letters of each record; empty when `stores` holds them.
    packed: Vec<Vec<u8>>,
}

impl<'a> PreparedSubjects<'a> {
    /// The records and their letters.
    fn letters(&self) -> SubjectLetters<'a> {
        SubjectLetters::new(self.records, self.stores)
    }

    /// Subject `index` as the per-subject kernel reads it.
    fn subject(&self, index: usize) -> BlastnSubjectRef<'_> {
        match self.stores {
            Some(stores) => BlastnSubjectRef::packed_record(&stores[index]),
            None => BlastnSubjectRef::record(&self.records[index], &self.packed[index]),
        }
    }
}

/// The subjects of a search as its caller gives them.
#[derive(Clone, Copy)]
enum Subjects<'a> {
    /// A FASTA subject set (`-subject`, or the records of `run_local`): the searched
    /// records (their intervals with `-subject_loc`), and their letters (with
    /// LOSAT_X_BNSTORE, packed).
    Records(SubjectLetters<'a>),
    /// An NCBI BLAST database (`-db`).
    Database(&'a database::BlastnDatabase),
}

/// The steps of NCBI's preliminary stage that finish the HSP list of one subject of an ungapped
/// (`-ungapped`) search (`Blast_HSPListReevaluateUngapped`, the e-values, the reap by the
/// preliminary e-value and by query coverage): the HSPs are re-evaluated with the ambiguity
/// codes of the subject (`s_blastna`, whole, without sentinel bytes), the ones below the
/// word cutoff of their context or the `-perc_identity` are deleted, the list is sorted by
/// score again, and the ones above the e-value (`evalue_of`) are deleted. The reap by the
/// query coverage (`-qcov_hsp_perc`) follows it, where the subject's HSPs are given to the
/// collector (`ungapped_reap_query_coverage`).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2609-2737
/// ```c
///       if (kNucleotideSubject) {
///          delete_hsp = Blast_HSPReevaluateWithAmbiguitiesUngapped(hsp, query_start,
///                 subject_start, word_params, sbp, kTranslateSubject);
///       }
///       if (!delete_hsp) {
///           const Uint1* query_nomask = query_blk->sequence_nomask +
///               query_info->contexts[context].query_offset;
///           Int4 align_length = 0;
///           Blast_HSPGetNumIdentitiesAndPositives(query_nomask, subject_start, hsp,
///                                score_params->options, &align_length, sbp);
///           delete_hsp = Blast_HSPTest(hsp, hit_params->options, align_length);
///       }
///    ...
///    /* Sort the HSP array by score (scores may have changed!) */
///    Blast_HSPListSortByScore(hsp_list);
///    Blast_HSPListAdjustOddBlastnScores(hsp_list, FALSE, sbp);
/// ```
#[allow(clippy::too_many_arguments)]
fn ungapped_subject_stage(
    hits: &mut Vec<PrelimHit>,
    s_blastna: &[u8],
    queries_blastna: &[&[u8]],
    score_matrix: &[i32; 256],
    word_cutoffs: &[i32],
    percent_identity: f64,
    evalue_of: impl Fn(&PrelimHit) -> f64,
    evalue_threshold: f64,
) {
    hits.retain_mut(|hit| {
        let context_idx = hit.context_idx as usize;
        let query = queries_blastna[context_idx];
        let length = hit.prelim_qe - hit.prelim_qs;
        let Some(reevaluated) = ungapped::reevaluate_with_ambiguities_ungapped(
            query,
            s_blastna,
            hit.prelim_qs,
            hit.prelim_ss,
            length,
            word_cutoffs[context_idx],
            score_matrix,
        ) else {
            return false;
        };
        hit.prelim_score = reevaluated.score;
        hit.prelim_qs = reevaluated.q_offset;
        hit.prelim_qe = reevaluated.q_end;
        hit.prelim_ss = reevaluated.s_offset;
        hit.prelim_se = reevaluated.s_end;
        let align_length = hit.prelim_qe - hit.prelim_qs;
        let identities = ungapped::count_identities_ungapped(
            query,
            s_blastna,
            hit.prelim_qs,
            hit.prelim_ss,
            align_length,
        );
        // `min_hit_length` is 0 for BLASTN.
        !ungapped::blast_hsp_test(identities, align_length, percent_identity, 0)
    });
    sort_prelim_hits_by_score(hits);
    hits.retain_mut(|hit| {
        hit.prelim_evalue = evalue_of(hit);
        hsp_survives_evalue_reap(hit.prelim_evalue, evalue_threshold)
    });
}

/// The last step of the preliminary stage's finish of a subject's HSP list in an ungapped
/// (`-ungapped`) search: the HSPs below `-qcov_hsp_perc` percent of their query are deleted
/// (`Blast_HSPListReapByQueryCoverage`). It does not depend on the subject, so the collector's
/// caller makes it, where the subject's HSPs are given over.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1482-1541
/// ```c
///             Blast_HSPListReapByQueryCoverage(hsp_list,hit_params->options, query_info, program_number);
/// ```
fn ungapped_reap_query_coverage(
    hits: &mut Vec<PrelimHit>,
    query_cov_hsp_perc: f64,
    context_length: impl Fn(usize) -> usize,
) {
    reap_by_query_coverage(hits, query_cov_hsp_perc, |hit| {
        (
            hit.prelim_qs as i64,
            hit.prelim_qe as i64,
            context_length(hit.context_idx as usize) as i64,
        )
    });
}

/// LOSAT's refusal of a word cutoff whose arithmetic NCBI does in a way that overflows (memory
/// "no NCBI UB / overflow"): a context of 2^30 letters or more, or a search space beyond `Int8`.
fn ungapped_cutoff_limit_error(limit: ungapped::CutoffLimit) -> anyhow::Error {
    let what = match limit {
        ungapped::CutoffLimit::QueryLength => {
            "a query of 2^30 letters or more, whose doubled length overflows NCBI BLAST+'s 32-bit int in the ungapped word cutoff"
        }
        ungapped::CutoffLimit::SearchSpace => {
            "an average subject length whose ungapped word cutoff search space overflows NCBI BLAST+'s 64-bit int"
        }
    };
    anyhow::anyhow!("{what} is not supported by LOSAT's BLASTN with -ungapped")
}

/// What the query batches search (`search_query_batch`): the prepared subject set, or the
/// database, which is read in blocks for each batch (`database::search_blocks`).
#[derive(Clone, Copy)]
enum SubjectSource<'a> {
    Prepared(&'a PreparedSubjects<'a>),
    Database(&'a database::BlastnDatabase),
}

impl<'a> SubjectSource<'a> {
    /// The number of sequences and the total length (the real ones, also with `-dbsize`).
    fn metadata(self) -> &'a SubjectMetadata {
        match self {
            SubjectSource::Prepared(prepared) => &prepared.metadata,
            SubjectSource::Database(database) => database.metadata(),
        }
    }

    /// The prepared subject set; none for a database.
    fn prepared(self) -> Option<&'a PreparedSubjects<'a>> {
        match self {
            SubjectSource::Prepared(prepared) => Some(prepared),
            SubjectSource::Database(_) => None,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:129-153
// ```c++
//     if(m_IsDbScan) {
// 	int num_seqs=0;
//         int total_length=0;
// 	if (!is_remote_search)
//         {
//                 BlastSeqSrc* seqsrc = db_adapter.MakeSeqSrc();
//                 num_seqs=BlastSeqSrcGetNumSeqs(seqsrc);
//                 total_length=static_cast<int>(BlastSeqSrcGetTotLen(seqsrc));
//         }
// 	CBlastFormatUtil::FillScanModeBlastDbInfo(m_DbInfo, m_DbIsAA,
//                 num_seqs, total_length, m_SubjectTag);
//     } else {
//     ...
//         CBlastFormatUtil::GetBlastDbInfo(m_DbInfo, m_DbName, m_DbIsAA,
//                                    dbfilt_algorithm, is_remote_search);
//     }
// ```
/// What the reports show of the subjects (NCBI's `m_DbInfo` and `m_DbName`): the name of
/// outfmt 7's `# Database:` line and of the outfmt 0 `Database:` lines of a subject set
/// (`User specified sequence set (Input: PATH)`; the `-db` value for a database), the
/// counts, and for a database the lines of `GetBlastDbInfo` (`PairwiseDatabase`).
struct ShownSubjects {
    name: String,
    num_seqs: usize,
    total_letters: usize,
    database: Option<PairwiseDatabase>,
    /// Whether `taxdb.bti` is on NCBI's search path, for the custom tabular reports
    /// (`BlastnArgs::tabular_taxdb`; not found where it was not asked).
    tabular_taxdb_found: bool,
}

impl ShownSubjects {
    fn new(args: &BlastnArgs, source: SubjectSource<'_>) -> Self {
        let metadata = source.metadata();
        match source {
            SubjectSource::Prepared(_) => ShownSubjects {
                name: format!(
                    "User specified sequence set (Input: {})",
                    args.subject_path().display()
                ),
                num_seqs: metadata.db_num_seqs,
                total_letters: metadata.db_len_total,
                database: None,
                tabular_taxdb_found: args.tabular_taxdb.unwrap_or(false),
            },
            SubjectSource::Database(database) => ShownSubjects {
                name: database.name().to_string(),
                num_seqs: metadata.db_num_seqs,
                total_letters: metadata.db_len_total,
                database: Some(database.pairwise_database()),
                tabular_taxdb_found: args.tabular_taxdb.unwrap_or(false),
            },
        }
    }
}

/// One subject as the per-subject kernel (`process_subject` of `search_query_batch`)
/// reads it: a record of a subject set, or a sequence of a database.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
/// ```c
/// if (subj_is_na) {
///     BlastSeqBlkSetSequence(subj, sequence.data.release(),
///        ((sentinels == eSentinels) ? sequence.length - 2 :
///         sequence.length));
///     ...
///     SBlastSequence compressed_seq =
///         subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
///                                   eNa_strand_plus, eNoSentinels);
///     BlastSeqBlkSetCompressedSequence(subj,
///                               compressed_seq.data.release());
/// }
/// ```
#[derive(Clone, Copy)]
struct BlastnSubjectRef<'a> {
    /// The number of bases.
    length: usize,
    /// The packed ncbi2na letters that the preliminary stage scans (four bases per byte,
    /// the count of the last byte's bases in its two low bits): the subject set's
    /// (`encode_subject_ncbi2na_packed`) or the database's stored data.
    packed: &'a [u8],
    /// The letters of a subject set's record: its lowercase masks (`-lcase_masking`) and
    /// the traceback's BLASTNA letters come from them. `None` for a database sequence.
    letters: Option<&'a [u8]>,
    /// The BLASTNA letters of a database sequence in the traceback
    /// (`database::traceback_listed`), without sentinel bytes; `None` otherwise.
    blastna: Option<&'a [u8]>,
    /// EXPERIMENT (LOSAT_X_BNSTORE): the packed subject of a subject set's record that keeps
    /// no letters: its lowercase masks and the traceback's BLASTNA letters come from it.
    store: Option<&'a PackedNucleotide>,
    /// EXPERIMENT (LOSAT_X_BNTBRANGE): with `blastna` decoded in ranges only, set when a
    /// traceback alignment met a fence (`alignment::x_fence`): the traceback stops, and
    /// `database::traceback_listed` traces the sequence again whole. `None` otherwise.
    x_fenced: Option<&'a std::sync::atomic::AtomicBool>,
}

impl<'a> BlastnSubjectRef<'a> {
    /// The record `record` of a subject set, with its packed letters `packed`.
    fn record(record: &'a FastaRecord, packed: &'a [u8]) -> Self {
        BlastnSubjectRef {
            length: record.seq().len(),
            packed,
            letters: Some(record.seq()),
            blastna: None,
            store: None,
            x_fenced: None,
        }
    }

    /// EXPERIMENT (LOSAT_X_BNSTORE): a record of a subject set held as the packed subject
    /// `store`; the preliminary stage scans its own packed bytes.
    fn packed_record(store: &'a PackedNucleotide) -> Self {
        BlastnSubjectRef {
            length: store.len(),
            packed: store.packed(),
            letters: None,
            blastna: None,
            store: Some(store),
            x_fenced: None,
        }
    }
}

/// The results of one query batch, with the values of its query contexts that the reports
/// use.
struct QueryBatch {
    hit_lists: Vec<Option<BlastnHitList>>,
    good_init_extends: u64,
    searched: bool,
    query_masks: Vec<Vec<MaskedInterval>>,
    context_karlin: Vec<Option<ContextKarlin>>,
    query_eff_searchsp: Vec<i64>,
    /// The preliminary hit lists of a query chunk, one per query part
    /// (`BatchStage::ChunkPrelim`).
    chunk_prelim_lists: PrelimHitLists,
    /// The batch's `LOSAT_TIMING` lines (`blastn_timing_text`; empty without it), which
    /// `run_in_pool` writes to standard error when it takes the batch, where the search wrote
    /// them before.
    timing: String,
}

/// Which of NCBI's stages `search_query_batch` runs for its queries.
#[derive(Clone, Copy)]
enum BatchStage<'a> {
    /// The preliminary search and the traceback of a query batch (`CLocalBlast::Run`). A
    /// batch that NCBI splits has its preliminary search run in query chunks
    /// (`search_query_chunks`), which need where each query's searched letters start in
    /// its record (`-query_loc`).
    ///
    /// `query_masks` are the masks of the batch's queries when the loop has made them
    /// (`-window_masker_db`, whose errors stop the batch before its search); `None`: made here.
    Search {
        query_offsets: &'a [usize],
        query_masks: Option<&'a [Vec<MaskedInterval>]>,
    },
    /// The preliminary search of one query chunk of a split batch
    /// (`SplitQuery_CreateChunkData`), with the masks of its query parts and the search
    /// spaces of the batch's contexts.
    ChunkPrelim {
        query_masks: &'a [Vec<MaskedInterval>],
        batch_eff_searchsp: &'a [i64],
    },
}

/// Searches the queries in NCBI's query batches and writes the results of all of them.
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:261-300
/// ```c++
///         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
///         int batch_size = m_CmdLineArgs->GetQueryBatchSize();
///         if (batch_size) {
///             input.SetBatchSize(batch_size);
///         } else {
///             Int8 total_len = formatter.GetDbTotalLength();
///             if (total_len > 0) {
///                 /* the optimal hits per batch scales with total db size */
///                 mixer.SetTargetHits(total_len / 3000);
///             }
///             input.SetBatchSize(mixer.GetBatchSize());
///         }
///         ...
///         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
///             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
///             ...
///                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
///                 ...
///                 results = lcl_blast.Run();
///                 if (!batch_size)
///                     input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()));
/// ```
/// `GetQueryBatchSize` is the `BATCH_SIZE` value (0 when it is not set) and the chunk size
/// that of `CHUNK_SIZE` or the task (`query_batching_from_environment`).
#[allow(clippy::too_many_arguments)]
fn run_in_pool(
    args: BlastnArgs,
    queries: x_query_stream::Queries<'_>,
    reading: QueryReading,
    query_input: &QueryInput,
    subjects: Subjects<'_>,
    subject_ranges: SubjectRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    mut batching: QueryBatching,
) -> Result<()> {
    let prepared: PreparedSubjects<'_>;
    let source = match subjects {
        Subjects::Records(subject_set) => {
            let subject_records = subject_set.records();
            let metadata =
                subject_metadata_from_records(subject_records, |index| subject_set.len(index));
            if metadata.db_num_seqs == 0 {
                // EXPERIMENT (LOSAT_X_QSTREAM): the records are read again for their warnings.
                queries.write_input_warnings(outputs.diagnostics)?;
                return Ok(());
            }
            // Diagnostic (LOSAT_TIMING): the subjects' packing for the common line
            // (utils/xtiming.rs). No NCBI counterpart.
            let t_encode = crate::utils::xtiming::clock();
            // The packed subjects hold the packed letters (LOSAT_X_BNSTORE).
            let packed = if subject_set.stores().is_some() {
                Vec::new()
            } else {
                subject_records
                    .iter()
                    .map(|record| encode_subject_ncbi2na_packed(record.seq()))
                    .collect()
            };
            crate::utils::xtiming::add(crate::utils::xtiming::Phase::Encode, t_encode);
            prepared = PreparedSubjects {
                records: subject_records,
                stores: subject_set.stores(),
                packed,
                metadata,
            };
            SubjectSource::Prepared(&prepared)
        }
        // A database without sequences fails in the first searched batch
        // (`search_query_batch`).
        Subjects::Database(database) => SubjectSource::Database(database),
    };
    // The subjects as the reports show them, and the stored titles of the database
    // sequences that the final hit lists of the batches name (`database::read_titles`).
    let shown = ShownSubjects::new(&args, source);
    let mut titles = database::DatabaseTitles::new();
    let megablast = args.task == "megablast";
    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:258-261
    // ```c
    //         formatter.PrintProlog();
    //
    //         /*** Process the input ***/
    //         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
    // ```
    // NCBI writes the outfmt 0 prolog before it reads the first query batch, and its stream
    // throws at the prolog's first flush when the write fails (blast_format.cpp:118-119),
    // before any query is read or searched.
    write_pairwise_prologs(outputs, &output_formats, megablast, &shown, true)?;
    let mut mixer = BatchSizeMixer::new(batching.split.mixer_max_batch_size());
    if batching.fixed_batch_size == 0 {
        // A database's total is the `Int8` of its index files (`GetDbTotalLength` sums
        // `m_DbInfo`, `GetBlastDbInfo`: the whole database's, also with a list).
        let total_length = shown
            .database
            .as_ref()
            .map_or(source.metadata().db_len_total as i64, |database| {
                database.total_letters as i64
            });
        if total_length > 0 {
            mixer.set_target_hits((total_length / 3000) as i32);
        }
        // `input.SetBatchSize(mixer.GetBatchSize())` converts the `Int4` to `TSeqPos`.
        batching.batch_size = mixer.batch_size(None) as u32;
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:277-280
    // ```c
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
    // 	    BLAST_PROF_START( APP.LOOP.PRE );
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    //             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // ```c
    // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
    //     : m_QueryVector(& queries)
    // {
    //     if (queries.Empty()) {
    //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
    //     }
    // ```
    // A reader at the end of its input reads no batch (`!input.End()`): the report of no
    // query, the prolog and the epilog (an empty pipe, or `-query -` after `-subject -`).
    if queries.input_count() == 0 && matches!(reading, QueryReading::End) {
        return write_batch_reports(
            &args,
            &[],
            ReportSubjects {
                source,
                shown: &shown,
                titles: &titles,
            },
            &ReportRanges {
                queries: query_input,
                subjects: subject_ranges,
            },
            outputs,
            &output_formats,
            BatchResults {
                hit_lists: Vec::new(),
                warnings: &[],
                lengths: &[],
                unsearched: Vec::new(),
                query_masks: &[],
                context_karlin: &[],
                query_eff_searchsp: &[],
                batches: &[],
            },
            true,
            None,
        );
    }
    // A batch size of 0 (`CHUNK_SIZE` 1000 without `BATCH_SIZE`: the mixer's maximum is
    // 0) reads no query, so the first batch fails after the prolog (outfmt 0).
    if batching.batch_size == 0 {
        write_pairwise_prologs(outputs, &output_formats, megablast, &shown, false)?;
        return Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
        }
        .into());
    }

    // The lengths of the searched queries (their intervals with `-query_loc`; with
    // LOSAT_X_QSTREAM from the first reading's lengths).
    let lengths: Vec<usize> = queries.searched_lengths();
    let report_ranges = ReportRanges {
        queries: query_input,
        subjects: subject_ranges,
    };
    let query_offsets: Vec<usize> = (0..lengths.len())
        .map(|index| query_input.placements.offset(index))
        .collect();
    // EXPERIMENT (LOSAT_X_RPTSTREAM): whether the loop writes each batch's reports when it takes
    // the batch (`x_rpt_stream.rs`); the summary lines come when it is dropped, at the end.
    let mut stream = x_rpt_stream::Stream::start(
        &args,
        outputs,
        &output_formats,
        source,
        subject_ranges,
        query_input,
        &lengths,
        batching,
        parallel_pool,
    );
    // EXPERIMENT (LOSAT_X_MBLUTKEEP): the megablast table this loop keeps between its batches
    // (`lookup::x_mb_lut_keep`), freed when the loop ends.
    let x_lut_keep = x_mb_lut_keep::Keep::for_loop();
    let mut hit_lists: Vec<Option<BlastnHitList>> = Vec::with_capacity(lengths.len());
    let mut batch_spans: Vec<std::ops::Range<usize>> = Vec::new();
    let mut unsearched = Vec::with_capacity(lengths.len());
    let mut query_masks = Vec::with_capacity(lengths.len());
    let mut context_karlin = Vec::with_capacity(2 * lengths.len());
    // The warnings written before each query's report (`QueryWarnings`).
    let mut warnings: Vec<Vec<u8>> = vec![Vec::new(); lengths.len()];
    let mut query_eff_searchsp = Vec::with_capacity(2 * lengths.len());
    // EXPERIMENT (LOSAT_X_QSTREAM, LOSAT_X_BNQDROP): the searched records of the batches taken,
    // kept where a writer reads them after their batch: the reports after the loop and the
    // error paths that write the reports of the batches before (LOSAT_X_RPTSTREAM not acting).
    // Without the switches the records read before the search are those.
    let keep = queries.per_batch() && !stream.acting();
    let mut kept: Vec<FastaRecord> = Vec::new();
    // A batch is read from the input records (`input_start..input_end`); its searched
    // queries are `start..end` (all of them without `-query_loc`). A batch size counts the
    // whole records, and a record whose range starts past its end is skipped
    // (`next_ranged_batch`).
    let mut reading = reading;
    // Whether what ended the reading after the records (`reading`) has come.
    let mut reading_ended = false;
    let mut input_start = 0;
    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:277-279
    // ```c
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
    // 	    BLAST_PROF_START( APP.LOOP.PRE );
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
    // The reader is at the end of its input after the last record, unless an `eEOF` or an
    // error ends the reading there; that comes in the batch being read when the records
    // run out before the batch reaches its size, and otherwise in a batch of its own.
    while input_start < query_input.input_lengths.len()
        || (reading.reads_past_records() && !reading_ended)
    {
        let (input_end, reached_size) = next_ranged_batch(
            &query_input.input_lengths,
            &query_input.skipped,
            input_start,
            batching.batch_size,
        );
        let std::ops::Range { start, end } = query_input.searched_in(input_start..input_end);
        let reading_ends_here = input_end == query_input.input_lengths.len()
            && !reached_size
            && reading.reads_past_records();
        // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:277-279
        // ```c
        //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
        // 	    BLAST_PROF_START( APP.LOOP.PRE );
        //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
        // ```
        // `CFastaReader` writes its messages about a batch's queries (their lines and
        // titles) when it reads them, after the report of the batch before: before the
        // report of the batch's first query (`QueryWarnings`). Skipped records are read too.
        // EXPERIMENT (LOSAT_X_QSTREAM): the batch's records, read now from the query file when
        // the switch acts (an error here is the file changed between the two readings).
        let mut batch_records =
            queries.batch_in(input_start..input_end, start..end, parallel_pool)?;
        let batch_title_warnings: Vec<u8> = batch_records
            .input()
            .iter()
            .flat_map(|record| record.warnings.iter().copied())
            .collect();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input.cpp:145-152
        // ```c
        //         try { q.Reset(m_Source->GetNextSequence(scope)); }
        //         catch (const CObjReaderParseException& e) {
        //             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
        //                 break;
        //             }
        //             throw;
        //         }
        // ```
        // An `eEOF` ends the batch with the records read so far. A reader error (or LOSAT's
        // rejection of a Seq-id line) stops the run while the batch is read: after the
        // reports of the batches before (without the epilog) and the messages of the
        // batch's records, which are neither searched nor reported.
        if reading_ends_here {
            reading_ended = true;
            if let QueryReading::Error {
                error,
                warnings: error_warnings,
            } = std::mem::replace(&mut reading, QueryReading::End)
            {
                // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were
                // written when the loop took them.
                let before = |outputs: &mut ReportOutputs<'_>| {
                    write_before_failed_batch(
                        &args,
                        queries.before(&kept, start),
                        ReportSubjects {
                            source,
                            shown: &shown,
                            titles: &titles,
                        },
                        &report_ranges,
                        outputs,
                        &output_formats,
                        BatchResults {
                            hit_lists,
                            warnings: &warnings[..start],
                            lengths: &lengths[..start],
                            unsearched,
                            query_masks: &query_masks,
                            context_karlin: &context_karlin,
                            query_eff_searchsp: &query_eff_searchsp,
                            batches: &batch_spans,
                        },
                    )
                };
                if stream.acting() {
                    stream.finish_sam(outputs)?;
                } else {
                    before(outputs)?;
                }
                outputs.diagnostics.write_all(&batch_title_warnings)?;
                outputs.diagnostics.write_all(&error_warnings)?;
                return Err(error.into_app_error());
            }
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
        // ```c
        // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
        //     : m_QueryVector(& queries)
        // {
        //     if (queries.Empty()) {
        //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
        //     }
        // ```
        // A batch of skipped records only (with `-query_loc`), or a batch that an `eEOF`
        // ends before its first record, fails after the reports of the batches before,
        // without the epilog.
        if start == end {
            // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were written
            // when the loop took them.
            let before = |outputs: &mut ReportOutputs<'_>| {
                write_before_failed_batch(
                    &args,
                    queries.before(&kept, start),
                    ReportSubjects {
                        source,
                        shown: &shown,
                        titles: &titles,
                    },
                    &report_ranges,
                    outputs,
                    &output_formats,
                    BatchResults {
                        hit_lists,
                        warnings: &warnings[..start],
                        lengths: &lengths[..start],
                        unsearched,
                        query_masks: &query_masks,
                        context_karlin: &context_karlin,
                        query_eff_searchsp: &query_eff_searchsp,
                        batches: &batch_spans,
                    },
                )
            };
            if stream.acting() {
                stream.finish_sam(outputs)?;
            } else {
                before(outputs)?;
            }
            outputs.diagnostics.write_all(&batch_title_warnings)?;
            return Err(crate::cli::NativeError {
                exit: 3,
                message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
            }
            .into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:518-530
        // ```c
        //     for(TSeqPos index = 0; index < queries.Size(); index++) {
        //         ENa_strand strand = eNa_strand_unknown;
        //
        //         try {
        //
        //             strand = s_BlastSetup_GetStrand(queries.GetStrand(index), prog,
        //                                           strand_opt);
        //             ...
        //             CRef<CBlastQueryFilteredFrames> frame_to_bsl =
        //                 s_GetRestrictedBlastSeqLocs(queries, index, qinfo, prog);
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:149-152
        // ```c
        //     // 1. Initialize the query data (borrow it from the factory)
        //     CRef<ILocalQueryData> query_data(qf->MakeLocalQueryData(&*options));
        //     retval->m_InternalData->m_Queries = query_data->GetSequenceBlk();
        //     retval->m_InternalData->m_QueryInfo = query_data->GetQueryInfo();
        // ```
        // With `-window_masker_db` the masks come first in each query's set-up: a counts
        // file NCBI cannot read, or a bad index of its hash table, fails every query of the
        // batch (`window_masker::batch_error`), also those without letters. The preliminary
        // search sets up the whole batch (`BlastSetupPreliminarySearchEx`, called by
        // `CBlastPrelimSearch::x_Init`) before it splits the batch into chunks, so a split
        // batch fails the same way; its chunks then make their own masks in the search, and
        // an error there (the whole queries' masks were made) is LOSAT's rejection
        // (`window_masker::rejection`). The masks made here are those the search uses.
        let mut window_masker_masks = None;
        if let Some(window_masker) = &args.window_masker {
            let searched = batch_records.searched();
            let masks = match window_masker.build_error() {
                Some(error) => Err(error.clone()),
                None => window_masker.batch_masks(
                    searched.iter().map(|record| record.seq()),
                    QueryStrand::from_option(&args.strand) == QueryStrand::Minus,
                ),
            };
            if let Err(error) = &masks {
                let Some((masker_warnings, failure)) =
                    super::super::window_masker::batch_error(error, searched.len())
                else {
                    return Err(window_masker_rejection(&args, error));
                };
                // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were
                // written when the loop took them.
                let before = |outputs: &mut ReportOutputs<'_>| {
                    write_before_failed_batch(
                        &args,
                        queries.before(&kept, start),
                        ReportSubjects {
                            source,
                            shown: &shown,
                            titles: &titles,
                        },
                        &report_ranges,
                        outputs,
                        &output_formats,
                        BatchResults {
                            hit_lists,
                            warnings: &warnings[..start],
                            lengths: &lengths[..start],
                            unsearched,
                            query_masks: &query_masks,
                            context_karlin: &context_karlin,
                            query_eff_searchsp: &query_eff_searchsp,
                            batches: &batch_spans,
                        },
                    )
                };
                if stream.acting() {
                    stream.finish_sam(outputs)?;
                } else {
                    before(outputs)?;
                }
                outputs.diagnostics.write_all(&batch_title_warnings)?;
                outputs.diagnostics.write_all(&masker_warnings)?;
                return Err(failure.into());
            }
            window_masker_masks = masks.ok();
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
        // ```c
        //         } catch (const CException& e) {
        //             ...
        //             CRef<CSearchMessage> m
        //                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
        //             messages[index].push_back(m);
        //             s_InvalidateQueryContexts(qinfo, index);
        //         }
        //     ...
        //     // Validate that at least one query context is valid
        //     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
        //         NCBI_THROW(CBlastException, eSetup, messages.ToString());
        //     }
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_aux.cpp:1013-1025
        // ```c
        // TSearchMessages::ToString() const
        // {
        //     string retval;
        //     ITERATE(vector<TQueryMessages>, qm, *this) {
        //         if (qm->empty()) {
        //             continue;
        //         }
        //         ITERATE(TQueryMessages, msg, *qm) {
        //             retval += (*msg)->GetMessage() + " ";
        //         }
        //     }
        //     return retval;
        // }
        // ```
        // A query without letters (`Sequence contains no data`, blast_setup.hpp:190-198) has
        // no valid context, and the strands of every other BLASTN query are valid when the
        // queries are set up: a batch of such queries only stops NCBI with one message per
        // query (`GetMessage` puts the severity first), after the reports of the batches
        // before (`CATCH_ALL`: `BLAST engine error: `, exit 3).
        if batch_records
            .searched()
            .iter()
            .all(|record| record.seq().is_empty())
        {
            // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were written
            // when the loop took them.
            let before = |outputs: &mut ReportOutputs<'_>| {
                write_before_failed_batch(
                    &args,
                    queries.before(&kept, start),
                    ReportSubjects {
                        source,
                        shown: &shown,
                        titles: &titles,
                    },
                    &report_ranges,
                    outputs,
                    &output_formats,
                    BatchResults {
                        hit_lists,
                        warnings: &warnings[..start],
                        lengths: &lengths[..start],
                        unsearched,
                        query_masks: &query_masks,
                        context_karlin: &context_karlin,
                        query_eff_searchsp: &query_eff_searchsp,
                        batches: &batch_spans,
                    },
                )
            };
            if stream.acting() {
                stream.finish_sam(outputs)?;
            } else {
                before(outputs)?;
            }
            outputs.diagnostics.write_all(&batch_title_warnings)?;
            let mut message = String::from("BLAST engine error: ");
            for _ in start..end {
                message.push_str("Warning: ");
                message.push_str(NO_DATA_MESSAGE);
                message.push(' ');
            }
            message.push('\n');
            return Err(crate::cli::NativeError { exit: 3, message }.into());
        }
        warnings[start].extend(batch_title_warnings);
        check_batch_window(&args, batch_records.batch_queries().records())?;
        // The batch's search, on the whole pool. EXPERIMENT (LOSAT_X_MBLUTKEEP): the loop keeps
        // its megablast table from batch to batch.
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:929
        // ```c
        //     opt.SetGappedMode(!args[kArgUngapped]);
        // ```
        // The instance of the batch search is chosen here, once for the batch.
        let searched = if args.ungapped {
            search_query_batch::<true>(
                &args,
                batch_records.batch_queries(),
                source,
                Some((&mut *outputs, &output_formats)),
                parallel_pool,
                batching,
                input_start == 0,
                BatchStage::Search {
                    query_offsets: &query_offsets[start..end],
                    query_masks: window_masker_masks.as_deref(),
                },
                x_lut_keep.as_ref(),
            )
        } else {
            search_query_batch::<false>(
                &args,
                batch_records.batch_queries(),
                source,
                Some((&mut *outputs, &output_formats)),
                parallel_pool,
                batching,
                input_start == 0,
                BatchStage::Search {
                    query_offsets: &query_offsets[start..end],
                    query_masks: window_masker_masks.as_deref(),
                },
                x_lut_keep.as_ref(),
            )
        };
        let batch = match searched {
            Ok(batch) => batch,
            Err(error) => match error.downcast::<ErrorAfterInvalidBatches>() {
                // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:295-318
                // ```c
                //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
                //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
                // 		        lcl_blast.SetBatchNumber( batch_num );
                //                 results = lcl_blast.Run();
                //     ...
                //                 ITERATE(CSearchResultSet, result, *results) {
                //                     formatter.PrintOneResultSet(**result, query_batch);
                //                 }
                //     ...
                //         }
                //         BLAST_PROF_START( APP.POST );
                //         formatter.PrintEpilog(opt);
                // ```
                // The reports of the batches before are written, then the error, which
                // skips `PrintEpilog` (the outfmt 0 footer and outfmt 7's last line).
                Ok(ErrorAfterInvalidBatches(error)) => {
                    // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were
                    // written when the loop took them.
                    let before = |outputs: &mut ReportOutputs<'_>| {
                        write_batch_reports(
                            &args,
                            queries.before(&kept, start),
                            ReportSubjects {
                                source,
                                shown: &shown,
                                titles: &titles,
                            },
                            &report_ranges,
                            outputs,
                            &output_formats,
                            BatchResults {
                                hit_lists,
                                warnings: &warnings[..start],
                                lengths: &lengths[..start],
                                unsearched,
                                query_masks: &query_masks,
                                context_karlin: &context_karlin,
                                query_eff_searchsp: &query_eff_searchsp,
                                batches: &batch_spans,
                            },
                            false,
                            None,
                        )
                    };
                    if stream.acting() {
                        stream.finish_sam(outputs)?;
                    } else {
                        before(outputs)?;
                    }
                    // The failing batch was read (its title warnings) before its search.
                    outputs.diagnostics.write_all(&warnings[start])?;
                    return Err(error);
                }
                // The warnings of the batches read so far, in NCBI's order (the outfmt 0
                // prolog of a Karlin-Altschul error of the first batch is written).
                // EXPERIMENT (LOSAT_X_RPTSTREAM): those of the batches before were written
                // with their reports; the off path writes no report here (the switch acts only
                // where LOSAT's refusals cannot come after the first batch, `x_rpt_stream`).
                Err(error) => {
                    let written = if stream.acting() { start } else { 0 };
                    for query_warnings in &warnings[written..=start] {
                        outputs.diagnostics.write_all(query_warnings)?;
                    }
                    return Err(error);
                }
            },
        };
        // The batch's LOSAT_TIMING block, where its search wrote it before.
        if !batch.timing.is_empty() {
            eprint!("{}", batch.timing);
        }
        let mut batch_hit_lists = batch.hit_lists;
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:264-275
        // ```c++
        //     CRef<ILocalQueryData> qdata = m_QueryFactory->MakeLocalQueryData(m_Options);
        //
        //     vector<TSeqLocInfoVector> subj_masks;
        //     TSeqAlignVector aligns =
        //         LocalBlastResults2SeqAlign(hsp_results,
        //                                    *qdata,
        //                                    *m_SeqInfoSrc,
        //                                    m_OptsMemento->m_ProgramType,
        //                                    m_Options->GetGappedMode(),
        //                                    m_Options->GetOutOfFrameMode(),
        //                                    subj_masks,
        //                                    m_ResultType);
        // ```
        // A database: the batch's hit lists are finished now (`finish_query_hit_list`, the
        // end of NCBI's traceback), and the definition lines of the sequences they name are
        // read in NCBI's order where NCBI makes the batch's Seq-aligns
        // (`database::read_titles`), before the HSPs of score 0 are dropped. A header file
        // that NCBI cannot read fails this batch, after the reports of the batches before
        // and the messages of this batch's reading.
        if let SubjectSource::Database(database) = source {
            let hitlist_size = args.report_counts.hitlist_size;
            for (q_idx, hit_list) in batch_hit_lists.iter_mut().enumerate() {
                if let Some(hit_list) = hit_list {
                    finish_query_hit_list(
                        hit_list,
                        lengths[start + q_idx],
                        hitlist_size,
                        args.max_hsps_per_subject.unwrap_or(0),
                        args.qcov_hsp_perc.unwrap_or(0.0),
                        args.subject_besthit,
                        args.culling_limit_value(),
                        args.ungapped,
                    );
                }
            }
            if let Err(error) = database::read_titles(database, &batch_hit_lists, &mut titles) {
                let error = match error.downcast::<database::HeaderReadError>() {
                    Ok(failure) => failure.error,
                    // LOSAT's refusal: the warnings of the batches read so far, as for the
                    // other errors of a batch (below).
                    // EXPERIMENT (LOSAT_X_RPTSTREAM): those of the batches before were
                    // written with their reports; the off path writes no report here (the
                    // switch acts only where every title is read before the search,
                    // `x_rpt_stream`).
                    Err(error) => {
                        let written = if stream.acting() { start } else { 0 };
                        for query_warnings in &warnings[written..=start] {
                            outputs.diagnostics.write_all(query_warnings)?;
                        }
                        return Err(error);
                    }
                };
                // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were
                // written when the loop took them.
                let before = |outputs: &mut ReportOutputs<'_>| {
                    write_before_failed_batch(
                        &args,
                        queries.before(&kept, start),
                        ReportSubjects {
                            source,
                            shown: &shown,
                            titles: &titles,
                        },
                        &report_ranges,
                        outputs,
                        &output_formats,
                        BatchResults {
                            hit_lists,
                            warnings: &warnings[..start],
                            lengths: &lengths[..start],
                            unsearched,
                            query_masks: &query_masks,
                            context_karlin: &context_karlin,
                            query_eff_searchsp: &query_eff_searchsp,
                            batches: &batch_spans,
                        },
                    )
                };
                if stream.acting() {
                    stream.finish_sam(outputs)?;
                } else {
                    before(outputs)?;
                }
                outputs.diagnostics.write_all(&warnings[start])?;
                return Err(error);
            }
            for hit_list in batch_hit_lists.iter_mut().flatten() {
                drop_hsps_of_score_zero(hit_list);
            }
            // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:304-312
            // ```c++
            //             if (isArchiveFormat) {
            //                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
            //                 m_Bah.ResetMessages();
            //             } else {
            //                 BlastFormatter_PreFetchSequenceData(*results, scope,
            //                 			                        fmt_args->GetFormattedOutputChoice());
            //                 ITERATE(CSearchResultSet, result, *results) {
            //                     formatter.PrintOneResultSet(**result, query_batch);
            //                 }
            // ```
            // Before the reports of a batch NCBI fetches every subject that the batch reports
            // (those with Seq-aligns: the HSPs of score 0 are dropped), and with `taxdb.bti`
            // on the search path the taxonomy database writes its messages there
            // (`BlastnDatabase::prefetch_messages`: the texts of opening it at the run's first
            // lookup, a volume's lookup of a taxid once per run): after the reader's messages
            // of the batch and before the warnings of its first query, which are added below.
            // They go with the warnings, where the first format is outfmt 0, 5, 6, 7 or 17
            // (`BlastnOutputFormat::prefetches_subjects`), so every writer of the batch's
            // reports writes them there (the reports after the loop, LOSAT_X_RPTSTREAM's
            // batch writes, the error paths that write the reports of the batches before).
            // Without `taxdb.bti` nothing is read. A refusal (a taxonomy database or a taxid
            // that LOSAT cannot read as NCBI does) is written as the titles' refusals above;
            // LOSAT_X_RPTSTREAM acts only where none can come (`x_rpt_stream::check`).
            if output_formats
                .first()
                .is_some_and(|format| format.prefetches_subjects())
                && database.taxdb_found()
            {
                // A custom tabular report with a plural name looks every taxid of its rows up,
                // taxid 0 too, which the prefetch does not: the rows of a batch with hits may open
                // the taxonomy database before a later batch's prefetch
                // (`TaxdbRun::open_for_rows`).
                let rows_open = matches!(
                    output_formats.first(),
                    Some(BlastnOutputFormat::CustomTabular(spec)) if spec.checks_taxdb()
                ) && batch_hit_lists
                    .iter()
                    .flatten()
                    .any(|hit_list| !hit_list.hsplist_array.is_empty());
                let prefetch = database
                    .prefetch_messages(&titles, &batch_hit_lists)
                    .and_then(|messages| {
                        if rows_open {
                            database.open_for_tabular_rows()?;
                        }
                        Ok(messages)
                    });
                match prefetch {
                    Ok(messages) => warnings[start].extend(messages),
                    Err(error) => {
                        let written = if stream.acting() { start } else { 0 };
                        for query_warnings in &warnings[written..=start] {
                            outputs.diagnostics.write_all(query_warnings)?;
                        }
                        return Err(error);
                    }
                }
            }
        }
        // The batch numbers its queries from 0.
        let first = start as u32;
        hit_lists.extend(batch_hit_lists.into_iter().map(|list| {
            list.map(|mut list| {
                for hsp_list in &mut list.hsplist_array {
                    hsp_list.query_index += first;
                    for hsp in &mut hsp_list.hsps {
                        hsp.q_idx += first;
                    }
                }
                list
            })
        }));
        batch_spans.push(start..end);
        unsearched.extend(std::iter::repeat_n(!batch.searched, end - start));
        query_masks.extend(batch.query_masks);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1450-1452
        // ```c
        //     if (results.HasWarnings()) {
        //         ERR_POST(Warning << results.GetWarningStrings());
        //     }
        // ```
        // The warning of an invalid query comes with its report (`PrintOneResultSet`),
        // before its preamble (`QueryWarnings`): a query without letters has the message of
        // its set-up (`SetupQueries_OMF`, above), any other invalid query that of its
        // Karlin-Altschul parameters.
        // A structured format (outfmt 8, 12) writes no result warnings
        // (`ReportFormat::writes_query_warnings`); the warnings go with the first format.
        let result_warnings = output_formats
            .first()
            .is_none_or(|format| format.writes_query_warnings());
        let filtered_empty = match source {
            SubjectSource::Database(database) => database.filtered_empty(),
            _ => false,
        };
        // The query's first valid context (`-strand minus`: its minus context).
        let first_context = QueryStrand::from_option(&args.strand).first_context();
        for (index, karlin) in batch
            .context_karlin
            .iter()
            .skip(first_context)
            .step_by(2)
            .enumerate()
        {
            if !result_warnings {
                continue;
            }
            let mut messages = Vec::new();
            if karlin.is_none() {
                messages.push(if lengths[start + index] == 0 {
                    NO_DATA_MESSAGE.to_string()
                } else {
                    crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string()
                });
            }
            if batch.searched && filtered_empty {
                messages.push(crate::report::query_warnings::EMPTY_DATABASE_MESSAGE.to_string());
            }
            warnings[start + index].extend(query_warning(
                "blastn",
                query_input.ordinal(start + index),
                &batch_records.searched()[index],
                &messages,
            ));
        }
        context_karlin.extend(batch.context_karlin);
        query_eff_searchsp.extend(batch.query_eff_searchsp);
        // NCBI's `Int4` count of the batch's successful initial extensions.
        if batching.fixed_batch_size == 0 {
            batching.batch_size =
                mixer.batch_size(Some(batch.good_init_extends as u32 as i32)) as u32;
        }
        // EXPERIMENT (LOSAT_X_RPTSTREAM): the batch's reports, written now, after the mixer
        // has its gie, as NCBI's loop prints a batch's result sets (`x_rpt_stream.rs`). The
        // batch's lists and values are taken back from the end of what the loop keeps, numbered
        // from 0 again, and written with the batch's part of the query ranges.
        if stream.acting() {
            let count = end - start;
            let at = hit_lists.len() - count;
            fn take<T>(kept: &mut Vec<T>, at: usize) -> Vec<T> {
                kept.drain(at..).collect()
            }
            let mut batch_lists = take(&mut hit_lists, at);
            for hsp_list in batch_lists
                .iter_mut()
                .flatten()
                .flat_map(|list| list.hsplist_array.iter_mut())
            {
                hsp_list.query_index -= first;
                for hsp in &mut hsp_list.hsps {
                    hsp.q_idx -= first;
                }
            }
            let batch_unsearched = take(&mut unsearched, at);
            let batch_masks = take(&mut query_masks, at);
            let batch_karlin = take(&mut context_karlin, 2 * at);
            let batch_searchsp = take(&mut query_eff_searchsp, 2 * at);
            let batch_input = x_rpt_stream::batch_input(query_input, start..end);
            let batch_ranges = ReportRanges {
                queries: &batch_input,
                subjects: subject_ranges,
            };
            // EXPERIMENT (LOSAT_X_PAIRWISEPAR): the pool formats the batch's outfmt 0 report while
            // it is idle (between the loop's searches).
            let x_pool = Some(parallel_pool);
            stream.write(outputs, count, false, |outputs, written| {
                write_reports(
                    &args,
                    batch_records.searched(),
                    ReportSubjects {
                        source,
                        shown: &shown,
                        titles: &titles,
                    },
                    &batch_ranges,
                    outputs,
                    &output_formats,
                    BatchResults {
                        hit_lists: batch_lists,
                        warnings: &warnings[start..end],
                        lengths: &lengths[start..end],
                        unsearched: batch_unsearched,
                        query_masks: &batch_masks,
                        context_karlin: &batch_karlin,
                        query_eff_searchsp: &batch_searchsp,
                        batches: &[],
                    },
                    false,
                    written,
                    x_pool,
                )
            })?;
            // The written queries' warnings are not read again (the error paths write those
            // of the current batch).
            for query_warnings in &mut warnings[start..end] {
                *query_warnings = Vec::new();
            }
        }
        // EXPERIMENT (LOSAT_X_QSTREAM, LOSAT_X_BNQDROP): the batch's records are kept for a later
        // writer (moved when the loop took them out of their store), or dropped.
        if keep {
            let searched = batch_records.into_searched();
            queries.note_kept(searched.len());
            kept.extend(searched);
        } else {
            drop(batch_records);
        }
        queries.done(input_end);
        input_start = input_end;
    }
    // EXPERIMENT (LOSAT_X_MBLUTKEEP): the loop's kept table is freed with the loop.
    drop(x_lut_keep);

    let reports = |outputs: &mut ReportOutputs<'_>| {
        write_batch_reports(
            &args,
            queries.all(&kept),
            ReportSubjects {
                source,
                shown: &shown,
                titles: &titles,
            },
            &report_ranges,
            outputs,
            &output_formats,
            BatchResults {
                hit_lists,
                warnings: &warnings,
                lengths: &lengths,
                unsearched,
                query_masks: &query_masks,
                context_karlin: &context_karlin,
                query_eff_searchsp: &query_eff_searchsp,
                batches: &batch_spans,
            },
            true,
            Some(parallel_pool),
        )
    };
    if !stream.acting() {
        return reports(outputs);
    }
    // EXPERIMENT (LOSAT_X_RPTSTREAM): every batch's reports are written; the epilogs remain
    // (NCBI's `PrintEpilog` after its loop), written by the report code for no query, with
    // the queries of every batch counted (outfmt 7's last line).
    stream.write(outputs, 0, true, |outputs, written| {
        write_reports(
            &args,
            &[],
            ReportSubjects {
                source,
                shown: &shown,
                titles: &titles,
            },
            &report_ranges,
            outputs,
            &output_formats,
            BatchResults {
                hit_lists: Vec::new(),
                warnings: &[],
                lengths: &[],
                unsearched: Vec::new(),
                query_masks: &[],
                context_karlin: &[],
                query_eff_searchsp: &[],
                batches: &[],
            },
            true,
            written,
            None,
        )
    })?;
    Ok(())
}

// No NCBI counterpart (when LOSAT reads the records only): it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BNQPACK): the length of the query file in bytes, as NCBI's `CFile` gives it
/// (`stat`, through symbolic links); `None` for standard input (`-`) or a path that cannot be
/// read.
fn x_query_file_bytes(query: &std::path::Path) -> Option<u64> {
    if query == std::path::Path::new("-") {
        return None;
    }
    std::fs::metadata(query).ok().map(|metadata| metadata.len())
}

// No NCBI counterpart (when LOSAT reads the records only): it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_QSTREAM): the letters of the largest query batch of a search of `task`, as
/// `run_in_pool` reads its batches (`query_batching_from_environment`): `BATCH_SIZE` when it is
/// set (a `TSeqPos`), else the `CBatchSizeMixer` maximum (`CHUNK_SIZE`, or the task's chunk size,
/// less 1000: 999,000 letters for blastn, 4,999,000 for megablast and dc-megablast; as the loop's
/// `TSeqPos`, so a maximum below 0 is a batch of about 4.29e9 letters); 0 when the variables cannot
/// be read (the search stops at them before its first batch).
fn x_query_batch_letters(task: &str) -> u64 {
    match query_batching_from_environment(super::super::coordination::task_uses_megablast_chunks(
        task,
    )) {
        Ok(batching) if batching.fixed_batch_size != 0 => u64::from(batching.batch_size),
        Ok(batching) => u64::from(batching.split.mixer_max_batch_size() as u32),
        Err(_) => 0,
    }
}

// No NCBI counterpart: LOSAT's verbose, debug and trace output (it changes no value NCBI computes).
/// EXPERIMENT (LOSAT_X_RPTSTREAM): whether a search writes verbose, debug or trace output while its
/// batches are searched (`-verbose`, `BLEMIR_DEBUG`, `LOSAT_DEBUG_BLASTN`,
/// `LOSAT_DEBUG_COORDS[_START]`, the blastn trace, the thread diagnostics). The switch then keeps
/// the reports written after the loop, where that output comes before them.
fn x_debug_output(args: &BlastnArgs) -> bool {
    args.verbose
        || [
            "BLEMIR_DEBUG",
            "LOSAT_DEBUG_BLASTN",
            "LOSAT_DEBUG_COORDS",
            "LOSAT_DEBUG_COORDS_START",
        ]
        .iter()
        .any(|name| std::env::var_os(name).is_some())
        || blastn_trace::enabled()
        || crate::utils::threading::diagnostics_enabled()
}

/// Whether the outfmt 0 prolog of `sink` is written before the search, as NCBI writes it
/// (`write_pairwise_prologs`). A file sink is created when its report is written
/// (`OutputSink::open`), so its prolog stays at the start of the report.
fn prolog_before_search(sink: &OutputSink<'_>) -> bool {
    !matches!(sink, OutputSink::File(_))
}

/// Whether the prolog of `format` comes with its reports (or, for a run that stops before any
/// report, from `write_pairwise_prologs`) and not before the search: the outfmt 0 prolog of a
/// file sink (`prolog_before_search`), and the start of the XML2 and JSON documents of outfmt
/// 16 and 15 for every sink. NCBI writes that start before the search too (`PrintProlog`), but
/// LOSAT's refusals of what NCBI could write (a title outside US-ASCII, DECISIONS 20) come
/// after the search and before the first byte, so the start is written with the first report, or
/// by the path of an NCBI error that stops the run before one.
///
/// Outfmt 20's line of the field names comes with the first report too (or from
/// `write_pairwise_prologs` for a run that stops before one), after LOSAT's refusals of the
/// custom tabular report; the taxonomy database's warning of its formatter, which a run that
/// can refuse never writes, comes before the search (`tabular::write_tabular_prolog_warning`).
fn prolog_with_report(format: &FormatOutput<'_>, output_format: BlastnOutputFormat) -> bool {
    match output_format {
        BlastnOutputFormat::Pairwise => !prolog_before_search(&format.sink),
        BlastnOutputFormat::Xml2(_) => true,
        BlastnOutputFormat::CustomTabular(spec) => {
            spec.kind == crate::report::tabular::TabularKind::CsvWithHeader
        }
        _ => false,
    }
}

/// Writes and flushes the outfmt 0 prolog of every outfmt 0 sink whose prolog comes before the
/// search (`before_search`) or, for a run that stops before any report, the prolog of the
/// others (`prolog_with_report`: the file sinks of outfmt 0 and every sink of outfmt 15 and 16).
fn write_pairwise_prologs(
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    megablast: bool,
    shown: &ShownSubjects,
    before_search: bool,
) -> Result<()> {
    for (format, &output_format) in outputs.formats.iter_mut().zip(output_formats) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-357
        // ```c++
        //     if (m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader) {
        //         const CBlastTabularInfo::EFieldDelimiter kDelim = CBlastTabularInfo::eComma;
        //         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
        //     ...
        //         tabinfo.PrintFieldSpecs();
        // ```
        // Outfmt 20: the warning of its formatter before the search, the line of the field
        // names with the first report or here before an NCBI error (`prolog_with_report`).
        if let (Some(spec), true) = (output_format.tabular_spec(), before_search) {
            super::super::tabular::write_tabular_prolog_warning(
                spec,
                shown.tabular_taxdb_found,
                &mut *outputs.diagnostics,
            )?;
        }
        let writes_prolog = match output_format {
            BlastnOutputFormat::Pairwise => prolog_before_search(&format.sink) == before_search,
            BlastnOutputFormat::Xml2(_) => !before_search,
            BlastnOutputFormat::CustomTabular(_) => {
                !before_search && prolog_with_report(format, output_format)
            }
            _ => false,
        };
        if writes_prolog {
            write_format_prolog(format, output_format, megablast, shown)?;
        }
    }
    Ok(())
}

/// Writes and flushes the prolog of one sink, for the formats that have one (outfmt 0, 15, 16,
/// 20).
fn write_format_prolog(
    format: &mut FormatOutput<'_>,
    output_format: BlastnOutputFormat,
    megablast: bool,
    shown: &ShownSubjects,
) -> Result<()> {
    if output_format == BlastnOutputFormat::Pairwise {
        let mut writer = format.sink.open()?;
        write_blastn_pairwise_prolog(
            &mut writer,
            NCBI_BLASTN_VERSION,
            megablast,
            &shown.name,
            shown.num_seqs,
            shown.total_letters,
            shown.database.as_ref(),
        )?;
        writer.flush()?;
    } else if let Some(kind) = output_format.xml2_kind() {
        let mut writer = format.sink.open()?;
        writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
        writer.flush()?;
    } else if let Some(spec) = output_format.tabular_spec() {
        // Outfmt 20's line of the field names (nothing for the other tabular kinds).
        if spec.kind == crate::report::tabular::TabularKind::CsvWithHeader {
            let mut writer = format.sink.open()?;
            super::super::tabular::write_tabular_prolog(spec, &mut writer)?;
            writer.flush()?;
        }
    }
    Ok(())
}

/// The results of the query batches searched so far, indexed by query (two contexts per
/// query).
struct BatchResults<'a> {
    hit_lists: Vec<Option<BlastnHitList>>,
    /// The warnings written before each query's report (`QueryWarnings`).
    warnings: &'a [Vec<u8>],
    lengths: &'a [usize],
    unsearched: Vec<bool>,
    query_masks: &'a [Vec<MaskedInterval>],
    context_karlin: &'a [Option<ContextKarlin>],
    query_eff_searchsp: &'a [i64],
    /// The searched queries of each query batch of the call, in order, for the reports of
    /// every batch written at once (the BLAST archive writes one document per batch); empty
    /// when the call is one batch.
    batches: &'a [std::ops::Range<usize>],
}

/// NCBI's message for a query or subject without letters.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup.hpp:190-198
/// ```c
///     TSeqPos size() const {
///         TSeqPos retval = x_Size();
///         if (retval == 0) {
///             NCBI_THROW(CBlastException, eInvalidArgument,
///                        "Sequence contains no data");
///         }
///         return retval;
///     }
/// ```
pub(crate) const NO_DATA_MESSAGE: &str = "Sequence contains no data";

/// Writes what NCBI has written when a query batch fails, before its error: the reports of
/// the batches before (`before`, the queries searched so far; without the epilog), or, for
/// the first batch, the outfmt 0 prolog of the sinks that write it with the reports (the
/// others wrote it before the search, `write_pairwise_prologs`).
#[allow(clippy::too_many_arguments)]
fn write_before_failed_batch(
    args: &BlastnArgs,
    before: &[FastaRecord],
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
) -> Result<()> {
    if before.is_empty() {
        write_pairwise_prologs(
            outputs,
            output_formats,
            args.task == "megablast",
            subjects.shown,
            false,
        )
    } else {
        write_batch_reports(
            args,
            before,
            subjects,
            ranges,
            outputs,
            output_formats,
            results,
            false,
            None,
        )
    }
}

/// The subjects as the reports read them: what the batches searched, how the reports show
/// it, and for a database the stored titles of the sequences that the final hit lists
/// name (`database::read_titles`; empty for a subject set).
#[derive(Clone, Copy)]
struct ReportSubjects<'a> {
    source: SubjectSource<'a>,
    shown: &'a ShownSubjects,
    titles: &'a database::DatabaseTitles,
}

/// Writes the reports of `query_records` (the queries of the batches searched), with
/// NCBI's epilog (`PrintEpilog`) when `epilog`; `x_pool` as `write_reports`'.
#[allow(clippy::too_many_arguments)]
fn write_batch_reports(
    args: &BlastnArgs,
    query_records: &[FastaRecord],
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
    epilog: bool,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
    write_reports(
        args,
        query_records,
        subjects,
        ranges,
        outputs,
        output_formats,
        results,
        epilog,
        x_rpt_stream::Written::default(),
        x_pool,
    )
}

/// `write_batch_reports` with what earlier calls wrote (`written`: nothing, except for
/// LOSAT_X_RPTSTREAM's calls, `x_rpt_stream`).
///
/// EXPERIMENT (LOSAT_X_PAIRWISEPAR): `x_pool` is the search's pool, whose threads are idle between
/// the loop's searches; `None` keeps the serial writer.
#[allow(clippy::too_many_arguments)]
fn write_reports(
    args: &BlastnArgs,
    query_records: &[FastaRecord],
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
    epilog: bool,
    written: x_rpt_stream::Written<'_>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
    // Diagnostic (LOSAT_TIMING): the report for the common line (utils/xtiming.rs); the time of
    // every call adds up. With streaming (`x_rpt_stream`) the calls come between the query
    // batches, while other batches are searched: that part is inside `search` too. No NCBI
    // counterpart.
    let t_output = crate::utils::xtiming::clock();
    let result = write_reports_untimed(
        args,
        query_records,
        subjects,
        ranges,
        outputs,
        output_formats,
        results,
        epilog,
        written,
        x_pool,
    );
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Output, t_output);
    result
}

#[allow(clippy::too_many_arguments)]
fn write_reports_untimed(
    args: &BlastnArgs,
    query_records: &[FastaRecord],
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
    epilog: bool,
    mut written: x_rpt_stream::Written<'_>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
    let megablast = args.task == "megablast";
    let config = configure_task(args);
    let mut hit_lists = results.hit_lists;
    // The SAM documents (outfmt 17), one per SAM format: the stream's, which only its last call
    // writes (the stream acts for one format), else documents of this call's own, which it
    // writes at its end (the off path's single call, and the calls that write the document of
    // the batches before an NCBI error).
    let mut own_sams: Vec<Option<crate::report::sam::SamDocument>> = Vec::new();
    let mut sam_targets: Vec<Option<SamTarget<'_>>> = Vec::new();
    // No SAM format: no target is made (the vectors stay empty and are never indexed).
    let has_sam = output_formats
        .iter()
        .any(|format| format.sam_spec().is_some());
    match written.sam.take() {
        Some(document) => sam_targets.push(Some(SamTarget {
            document,
            finish: epilog,
        })),
        None if !has_sam => {}
        None => {
            own_sams.extend(output_formats.iter().map(|format| {
                format
                    .sam_spec()
                    .map(|spec| x_rpt_stream::sam_document(args, spec))
            }));
            sam_targets.extend(own_sams.iter_mut().map(|document| {
                document.as_mut().map(|document| SamTarget {
                    document,
                    finish: true,
                })
            }));
        }
    }
    // A database: the reported sequences, named by rank from here on, with their stored
    // titles, and their stored letters where a report shows them (outfmt 0 and the hit
    // records; `database::reported_subjects`).
    let reported = match subjects.source {
        SubjectSource::Prepared(_) => None,
        SubjectSource::Database(database) => {
            let letters = outputs.hits.is_some()
                || output_formats.contains(&BlastnOutputFormat::Pairwise)
                || output_formats
                    .iter()
                    .any(|format| format.is_structured_xml())
                || output_formats
                    .iter()
                    .filter_map(|format| format.tabular_spec())
                    .any(|spec| spec.needs().rows);
            let pairwise = output_formats.contains(&BlastnOutputFormat::Pairwise);
            let shown = crate::algorithm::common::database::ParsedShown {
                tabular: output_formats.iter().any(|format| format.is_tabular()),
                pairwise,
                joined: pairwise
                    || output_formats
                        .iter()
                        .any(|format| format.is_structured_xml()),
                sam: output_formats
                    .iter()
                    .any(|format| matches!(format, BlastnOutputFormat::Sam(_))),
            };
            Some(database::reported_subjects(
                database,
                subjects.titles,
                &mut hit_lists,
                letters,
                shown,
                written.queries,
            )?)
        }
    };
    let database_subject_ids = match &reported {
        Some(reported) if output_formats.iter().any(|format| format.is_tabular()) => {
            Some(reported.tabular_subject_ids()?)
        }
        _ => None,
    };
    let subject_letters = match (&reported, subjects.source) {
        (Some(reported), _) => SubjectLetters::new(&reported.records, None),
        (None, SubjectSource::Prepared(prepared)) => prepared.letters(),
        (None, SubjectSource::Database(_)) => SubjectLetters::new(&[], None),
    };
    // The archive's `LCaseMask` holds the DUST intervals of every query, valid or not, whether
    // or not its batch was searched; with `-window_masker_db` those of WindowMasker (empty for a
    // query of a batch whose masks failed: its batch writes no archive).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:456-458
    // ```c
    //     } else if (m_QueryVector.NotEmpty()) {
    //         for (CBlastQueryVector::size_type i = 0; i < m_QueryVector->Size(); i++)
    //             retval.push_back(m_QueryVector->GetMaskedRegions(i));
    // ```
    // (`ExtractUserSpecifiedMasks`, which the archive's request reads, remote_blast.cpp:1358-1359:
    // the query vector's masks after `x_CalculateMasks`, both frames of every interval.)
    let archive_dust = if output_formats.contains(&BlastnOutputFormat::Archive) {
        match &args.window_masker {
            Some(window_masker) => {
                let minus = QueryStrand::from_option(&args.strand) == QueryStrand::Minus;
                query_records
                    .iter()
                    .map(|record| window_masker.masks(record.seq(), minus).unwrap_or_default())
                    .collect()
            }
            None => {
                let serial = crate::utils::threading::SearchPool::x_serial();
                super::super::coordination::apply_dust_masking(
                    args,
                    query_records,
                    x_pool.unwrap_or(&serial),
                )
            }
        }
    } else {
        Vec::new()
    };
    let report_inputs = BlastnReportInputs {
        queries: query_records,
        subjects: subject_letters,
        ranges,
        query_masks: results.query_masks,
        archive_dust: &archive_dust,
        lcase_masking: args.lcase_masking,
        query_eff_searchsp: results.query_eff_searchsp,
        first_context: QueryStrand::from_option(&args.strand).first_context(),
        megablast,
        reward: config.reward,
        penalty: config.penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2272
        // ```c
        //         if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
        // ```
        zero_gap_extension_formula: matches!(args.task.as_str(), "megablast" | "blastn"),
        window_size: config.window_size,
        ungapped: args.ungapped,
        context_karlin: results.context_karlin,
        num_descriptions: args.report_counts.num_descriptions,
        num_alignments: args.report_counts.num_alignments,
        culling_limit: args.culling_limit_value(),
        db_num_seqs: subjects.shown.num_seqs,
        db_len_total: subjects.shown.total_letters,
        database: subjects.shown.database.as_ref(),
        database_subject_ids: database_subject_ids.as_deref(),
        database_ordinals: reported
            .as_ref()
            .map(|reported| reported.ordinals.as_slice()),
        database_reader: match subjects.source {
            SubjectSource::Database(database) => Some(database),
            SubjectSource::Prepared(_) => None,
        },
        database_parsed: reported
            .as_ref()
            .map_or(&[][..], |reported| reported.parsed.as_slice()),
        finished_hit_lists: reported.is_some(),
        epilog,
        prolog: outputs
            .formats
            .iter()
            .zip(output_formats)
            .any(|(format, &output_format)| {
                output_format == BlastnOutputFormat::Pairwise && !prolog_before_search(&format.sink)
            }),
        queries_before: written.queries,
        subject_cache: written.subjects,
        x_pairwise_pool: x_pool.and_then(crate::report::x_pairwise_par::pool),
        x_tab_pool: x_pool.and_then(crate::algorithm::blastn::x_tab_par::pool),
        out: args.out.as_deref(),
        expect: super::super::coordination::determine_evalue(args),
        dust: args.dust.params().map(|(level, window, linker)| {
            (level as i32, window as u32 as i32, linker as u32 as i32)
        }),
        window_masker_db: args.window_masker_db.as_deref(),
        eff_db_length: if args.db_length() > 0 {
            args.db_length()
        } else {
            subjects.shown.total_letters as i64
        },
        archive: args.archive.as_ref(),
        batch_spans: results.batches,
        both_strands: args.strand == "both",
        minus_strand: QueryStrand::from_option(&args.strand) == QueryStrand::Minus,
        tabular_taxdb_found: subjects.shown.tabular_taxdb_found,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:482-496
    // ```c++
    //         if (sid_in->IsLocal()) {
    //             string id_token;
    //             vector<string> title_tokens;
    //             title_tokens =
    //                 NStr::Split(CAlignFormatUtil::GetTitle(bh), " ", title_tokens);
    //             if(title_tokens.empty()){
    //                 id_token = NcbiEmptyString;
    //             } else {
    //                 id_token = title_tokens[0];
    //             }
    //
    //             if (id_token == NcbiEmptyString || parse_local) {
    //                 const CObject_id& obj_id = sid_in->GetLocal();
    //                 if (obj_id.IsStr())
    //                     id_token = obj_id.GetStr();
    // ```
    // The query ID of the tabular formats is the title's first word, or the local ID
    // (`FastaRecord::shown_id`).
    let query_ids: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.shown_id()))
        .collect();
    let query_titles: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.query_line().as_ref()))
        .collect();
    // `User specified sequence set (Input: PATH)` for a subject set, the `-db` value for a
    // database (`ShownSubjects`).
    let subject_title = subjects.shown.name.as_str();
    let hitlist_size = args.report_counts.hitlist_size;
    post_process_hits_and_write(
        hit_lists,
        hitlist_size,
        args.max_hsps_per_subject.unwrap_or(0),
        args.qcov_hsp_perc.unwrap_or(0.0),
        args.subject_besthit,
        results.lengths,
        outputs,
        args.verbose,
        &query_ids,
        output_formats,
        &query_titles,
        subject_title,
        None,
        &report_inputs,
        results.unsearched,
        results.warnings,
        sam_targets,
    )
}

/// The outputs that a query batch's search writes to: the outfmt 0 prolog of a file sink
/// before an error of the first batch (`write_pairwise_prologs`), and nothing otherwise.
/// `run_in_pool` gives them to the batches it searches (a split batch to its query chunks).
type BatchOutputs<'o, 'r> = Option<(&'o mut ReportOutputs<'r>, &'o [BlastnOutputFormat])>;

/// Why a search of the first query batch has its outputs (`BatchOutputs`).
const FIRST_BATCH_OUTPUTS: &str =
    "run_in_pool searches the first query batch itself, with the outputs";

/// `outputs` again, for one call.
fn reborrow_outputs<'o, 'r>(outputs: &'o mut BatchOutputs<'_, 'r>) -> BatchOutputs<'o, 'r> {
    outputs
        .as_mut()
        .map(|(outputs, formats)| (&mut **outputs, &**formats))
}

/// The length adjustment and the effective search space of a valid query context of `query_len`
/// letters against a database of `db_len` letters in `db_num_seqs` sequences
/// (`BLAST_CalcEffLengths`), with the context's gapped Karlin-Altschul block. The search space
/// is the context's `eff_searchsp` when no `-searchsp` replaces it. The reports compute the
/// length adjustment again, which the search does not keep.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-846
/// ```c
/// BLAST_ComputeLengthAdjustment(..., query_length, db_length, db_num_seqs, &length_adjustment);
/// effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
/// if (effective_db_length <= 0) effective_db_length = 1;
/// effective_search_space = effective_db_length * (query_length - length_adjustment);
/// ```
fn context_lengths(
    query_len: i64,
    db_len: i64,
    db_num_seqs: i64,
    gapped: &KarlinParams,
) -> Result<(i64, i64)> {
    let result = compute_length_adjustment_ncbi(query_len, db_len, db_num_seqs, gapped);
    let length_adjustment = result.length_adjustment;
    // NCBI's `Int8` arithmetic, which overflows where the database length (a
    // `-dbsize`, or the total of a large database) times the query length passes
    // 2^63 (never for a subject set without `-dbsize`: its total and the query
    // lengths are below 2^31); LOSAT rejects such a search instead of reproducing
    // the overflow.
    let overflow = || {
        anyhow::anyhow!(
            "a database length of {db_len} letters with a query of {query_len} letters overflows NCBI BLAST+'s 64-bit effective search space, which is not supported by LOSAT's BLASTN"
        )
    };
    let effective_db_length = db_num_seqs
        .checked_mul(length_adjustment)
        .and_then(|adjusted| db_len.checked_sub(adjusted))
        .ok_or_else(overflow)?
        .max(1);
    let effective_query_length = (query_len - length_adjustment).max(1);
    let effective_search_space = effective_db_length
        .checked_mul(effective_query_length)
        .ok_or_else(overflow)?;
    Ok((length_adjustment, effective_search_space))
}

/// LOSAT's rejection of an error of the masks of `-window_masker_db` that the batch loop does
/// not report as NCBI does (`window_masker::rejection`).
fn window_masker_rejection(
    args: &BlastnArgs,
    error: &crate::utils::winmask::WinMaskError,
) -> anyhow::Error {
    super::super::window_masker::rejection(
        args.window_masker_db.as_deref().unwrap_or_default(),
        error,
    )
}

/// Searches one query batch against the subjects, with the query block, lookup table and
/// diagonal table of the batch, as NCBI's `CLocalBlast` does for each batch that
/// `run_in_pool` reads. `x_keep`: the batch loop's keep of LOSAT_X_MBLUTKEEP, which the batch's
/// megablast table is built on and goes back to (`lookup::x_mb_lut_keep`).
/// `blastinput::app::check_diag_table_window` for a query batch (DI-6, `-window_size`), before
/// `search_query_batch`, whose code a check inside would change (its inlining; the default
/// search's code must stay the same). The concatenated length is that of the batch's contexts
/// (`search_query_batch`): each searched strand of a query with letters takes its letters and
/// a sentinel, and the searched block leaves out the last sentinel.
#[inline(never)]
fn check_batch_window(args: &BlastnArgs, records: &[FastaRecord]) -> Result<()> {
    let window = super::super::coordination::determine_window_size(args);
    let (plus_on, minus_on) = QueryStrand::from_option(&args.strand).searched();
    let strands = usize::from(plus_on) + usize::from(minus_on);
    let rooms: usize = records
        .iter()
        .map(|record| record.sequence.len())
        .filter(|&length| length > 0)
        .map(|length| strands * (length + 1))
        .sum();
    crate::blastinput::app::check_diag_table_window(
        i32::try_from(rooms.saturating_sub(1)).unwrap_or(i32::MAX),
        i32::try_from(window).unwrap_or(i32::MAX),
        "BLASTN",
    )
}

fn search_query_batch<const UNGAPPED: bool>(
    args: &BlastnArgs,
    mut query_records: x_query_drop::BatchQueries<'_>,
    subjects: SubjectSource<'_>,
    mut outputs: BatchOutputs<'_, '_>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batching: QueryBatching,
    first_batch: bool,
    stage: BatchStage<'_>,
    x_keep: Option<&Arc<x_mb_lut_keep::Keep>>,
) -> Result<QueryBatch> {
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
    // ```c
    // const int kMaxValue = static_cast<int>(CSystemInfo::GetCpuCount());
    // ...
    // int num_threads = args[kArgNumThreads].AsInteger();
    // if (num_threads > kMaxValue) {
    //     m_NumThreads = kMaxValue;
    // } else {
    //     m_NumThreads = num_threads;
    // }
    // ```
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    let num_threads = parallel_pool.threads();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:493-555
    // ```c
    // if (aux_struct->WordFinder) {
    //     aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
    //                            word_params, aux_struct->ewp,
    //                            aux_struct->offset_pairs,
    //                            kScanSubjectOffsetArraySize,
    //                            init_hitlist, ungapped_stats);
    //     if (init_hitlist->total == 0) continue;
    // }
    // ...
    // if (aux_struct->GetGappedScore) {
    //     status = aux_struct->GetGappedScore(program_number, query,
    //             query_info, subject, gap_align, score_params, ext_params,
    //             hit_params, word_params, init_hitlist, &hsp_list,
    //             gapped_stats, NULL);
    // }
    // ...
    // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    // Blast_HSPListSortByScore(hsp_list);
    // ```
    let timing_enabled = std::env::var_os("LOSAT_TIMING").is_some();
    let t_total = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    let timing = if timing_enabled {
        Some(Arc::new(BlastnTiming::new()))
    } else {
        None
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:478-536
    // ```c
    // while (TRUE) {
    //     status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
    //                                    dbseq_chunk_overlap);
    //     if (status == SUBJECT_SPLIT_DONE) break;
    //     if (status == SUBJECT_SPLIT_NO_RANGE) continue;
    //     ...
    //     if (aux_struct->WordFinder) {
    //         aux_struct->WordFinder(...);
    //         if (init_hitlist->total == 0) continue;
    //     }
    //     ...
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let requested_parallel = num_threads > 1;
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let requested_parallel = false;

    // Configure task-specific parameters (initial configuration)
    let mut config = configure_task(&args);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:509-513
    // ```c
    // BLAST_GappedAlignmentWithTraceback(program_number, query,
    //     adjusted_subject, gap_align, score_params, q_start, s_start,
    //     query_length, adjusted_s_length, fence_hit);
    // ```
    // Preserve ordered debug-coordinate logging; speculative DP is otherwise pure.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,503-507
    // ```c
    //    for (index=0; index < num_initial_hsps; index++) {
    // ...
    //         } else if (kGreedyTraceback) {
    //             BLAST_GreedyGappedAlignment(query, adjusted_subject,
    //                 query_length, adjusted_s_length, gap_align,
    //                 score_params, q_start, s_start, FALSE, TRUE,
    //                 fence_hit);
    // ```
    // The dynamic-programming traceback and, in a parallel search, the greedy traceback of megablast
    // (the call above) are scheduled ahead of the ordered loop (the rule of `LOSAT_X_GREEDYSPEC`,
    // kept as its `x_reach!` key). The ordered loop still tests the tree and adds to it in the NCBI
    // order. This is scheduling only.
    // No NCBI counterpart: the same scheduling for the greedy
    // traceback of megablast. A greedy alignment is a function of its
    // arguments only: every cell it reads from its scratch was written
    // earlier in the same call, and a larger `max_dist` left by an earlier
    // call only saves retries.
    let x_greedy_spec = requested_parallel && !config.use_dp;
    let speculative_traceback = requested_parallel
        && (config.use_dp || x_greedy_spec)
        && std::env::var_os("LOSAT_DEBUG_COORDS").is_none()
        && std::env::var_os("LOSAT_DEBUG_COORDS_START").is_none();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3918-3919
    // ```c
    //    for (index=0; index<init_hitlist->total; index++)
    // ...
    //       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
    //                                         hit_options->min_diag_separation))
    // ```
    // Look-ahead for the preliminary loop (above) and the traceback loop
    // (`Blast_TracebackFromHSPList`, next block). The helper threads only evaluate the alignment call
    // of an index ahead of the loop. The loops, their tree tests and their tree insertions stay on the
    // calling thread in the NCBI order.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,404-405
    // ```c
    //    for (index=0; index < num_initial_hsps; index++) {
    // ...
    //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
    //                                hit_options->min_diag_separation)) {
    // ```
    // The same for the loop of `Blast_TracebackFromHSPList`.
    // No NCBI counterpart: the other threads of the pool
    // evaluate the gapped alignments of the two ordered loops (preliminary
    // and traceback) ahead of the loop; see utils/xahead.rs. Same conditions
    // as the speculative traceback above.
    let x_ahead_window = if speculative_traceback {
        Some(crate::utils::xahead::WINDOW)
    } else {
        None
    };
    if x_ahead_window.is_none() {
        crate::utils::xreach::x_reach!("LOSAT_X_AHEAD", false);
    }
    // Read sequences
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    // The search reads the query records of its caller (`prepare_sequence_data`).
    let queries = query_records.records();
    // The records' letters for the readers below (`x_query_drop::Letters`).
    let letters = query_records.letters();
    let query_ids = queries
        .iter()
        .map(|record| String::from_utf8_lossy(record.shown_id()).into_owned())
        .collect();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1312-1325
    // ```c
    // if (lookup_options->db_filter) {
    //    s_FillPV(query, location, mb_lt, lookup_options);
    //    s_ScanSubjectForWordCounts(seqsrc, mb_lt, counts,
    //                               lookup_options->max_db_word_count);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
    // ```c
    // BlastSeqBlkSetSequence(subj, sequence.data.release(),
    //    ((sentinels == eSentinels) ? sequence.length - 2 :
    //     sequence.length));
    // ...
    // SBlastSequence compressed_seq =
    //     subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
    //                               eNa_strand_plus, eNoSentinels);
    // BlastSeqBlkSetCompressedSequence(subj, compressed_seq.data.release());
    // ```
    // Keep one fetched subject record per oid so metadata, lookup filtering,
    // and search reuse the same subject material instead of reparsing FASTA.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1401-1427
    // ```c
    // memset((void*) &seq_arg, 0, sizeof(seq_arg));
    // seq_arg.encoding = eBlastEncodingProtein;
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //         continue;
    //     }
    // }
    // ```
    let subject_records: Option<&[FastaRecord]> =
        subjects.prepared().map(|prepared| prepared.records);
    let subject_metadata = subjects.metadata().clone();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1475
    // ```c
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    status = s_BlastSearchEngineCore(...);
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:478-500
    // ```c
    // while (TRUE) {
    //     status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
    //                                    dbseq_chunk_overlap);
    //     if (status == SUBJECT_SPLIT_DONE) break;
    //     if (status == SUBJECT_SPLIT_NO_RANGE) continue;
    //     if (aux_struct->WordFinder) {
    //         aux_struct->WordFinder(...);
    // ```
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    // The subjects as tasks of the pool; `use_parallel` below, after LOSAT_X_BNFEWSUBJ's rule.
    let subjects_parallel = requested_parallel && subject_metadata.db_num_seqs > 1;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
    // ```c
    // if (subj_is_na) {
    //     BlastSeqBlkSetSequence(subj, sequence.data.release(),
    //        ((sentinels == eSentinels) ? sequence.length - 2 :
    //         sequence.length));
    //     ...
    //     SBlastSequence compressed_seq =
    //         subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
    //                                   eNa_strand_plus, eNoSentinels);
    //     BlastSeqBlkSetCompressedSequence(subj,
    //                               compressed_seq.data.release());
    // }
    // ```
    // The packed letters of each record for the lookup build, which reads them only for
    // `limit_lookup` (not a BLASTN option, always off); `None` with LOSAT_X_BNSTORE.
    let subject_packed_cache_ref = subjects
        .prepared()
        .filter(|prepared| prepared.stores.is_none())
        .map(|prepared| prepared.packed.as_slice());

    // Prepare sequence data (including DUST masking)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1401-1405
    // ```c
    // /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
    //   sequences are retieved in ncbistdaa/ncbi2na encodings respectively. */
    // seq_arg.encoding = eBlastEncodingProtein;
    // ```
    // Diagnostic (LOSAT_TIMING): the batch's DUST for the common line (utils/xtiming.rs); a
    // chunk takes the batch's masks. No NCBI counterpart.
    let t_seg = crate::utils::xtiming::clock();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:445-446
    // ```c
    //     TMaskedQueryRegions mqr =
    //         queries.GetMaskedRegions(query_index);
    // ```
    // The batch loop made the masks of `-window_masker_db` when it set the batch up (the
    // errors stop the batch there); the search takes them.
    let masks = match stage {
        BatchStage::Search {
            query_masks: Some(masks),
            ..
        } => masks.to_vec(),
        BatchStage::Search { .. } => query_masks(args, &letters, parallel_pool)
            .map_err(|error| window_masker_rejection(args, &error))?,
        BatchStage::ChunkPrelim { query_masks, .. } => query_masks.to_vec(),
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Seg, t_seg);
    let seq_data = prepare_sequence_data(queries, query_ids, masks, subject_metadata);

    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1264-1284
    // ```c
    // x_PrintQueryAndDbNames(program_version, bioseq, dbname, rid, iteration,
    //                        subj_bioseq);
    // ```
    let query_titles_arc = Arc::new(
        seq_data
            .queries
            .iter()
            .map(record_title)
            .collect::<Vec<Arc<[u8]>>>(),
    );
    // The subjects as the reports show them (the outfmt 0 prolog of an error below).
    let shown = ShownSubjects::new(args, subjects);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2589-2593
    // ```c
    // curr_context = hsp_array[i]->context;
    // qlen = query_info->contexts[curr_context].query_length;
    // target_context = (hsp_array[i]->query.frame > 0) ?
    //                  curr_context + 1 : curr_context - 1;
    // ```
    let query_lengths = Arc::new(
        (0..letters.count())
            .map(|index| letters.length(index))
            .collect::<Vec<usize>>(),
    );

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to. */
    //    BlastHSP** hsp_array;
    //    Int4 hspcnt;
    //    ...
    // } BlastHSPList;
    // ```
    let query_ids_arc = Arc::new(
        seq_data
            .query_ids
            .iter()
            .map(|id| Arc::<str>::from(id.as_str()))
            .collect::<Vec<Arc<str>>>(),
    );

    // Build per-context query sequences (plus/minus) for blastn.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:839-847
    // ```c
    // if (context_number % NUM_STRANDS == 0) frame = 1;
    // else frame = -1;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:166-174
    // ```c
    // query_info->contexts[0].query_offset = 0;
    // query_info->contexts[1].query_offset = kStrandLength + 1;
    // ```
    let mut query_contexts: Vec<QueryContext> = Vec::new();
    let mut query_context_masks: Vec<Vec<MaskedInterval>> = Vec::new();
    let mut query_base_offsets: Vec<usize> = Vec::new();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:166-174
    // ```c
    // query_info->contexts[0].query_offset = 0;
    // query_info->contexts[1].query_offset = kStrandLength + 1;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:71-81
    // ```c
    //     if (index) {
    //         Uint4 prev_loc = qinfo->contexts[index-1].query_offset;
    //         Uint4 prev_len = qinfo->contexts[index-1].query_length;
    //
    //         Uint4 shift = prev_len ? prev_len + 1 : 0;
    //
    //         qinfo->contexts[index].query_offset = prev_loc + shift;
    //         qinfo->contexts[index].query_length = length;
    //         if (length == 0)
    //            qinfo->contexts[index].is_valid = false;
    //     } else {
    // ```
    // Each context starts after the context before and its sentinel; the contexts of a
    // query without letters take no room (and are invalid, `context_ungapped` below).
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:229-243
    // ```c
    //                 case eNa_strand_plus:
    //                     s_QueryInfo_SetContext(query_info, ctx_index, length);
    //                     s_QueryInfo_SetContext(query_info, ctx_index+1, 0);
    //                     // the missing strand is present in query_info as
    //                     // zero-lenghth context
    //                     min_length = 0;
    // ```
    // `-strand` sets the strand of the query's Seq-loc: the strand that is not searched has a
    // context of length 0 (invalid, without masks and without room in the buffer). Chosen once
    // here; `both` keeps the plain two-strand code below.
    let strand = QueryStrand::from_option(&args.strand);
    let (plus_on, minus_on) = strand.searched();
    let mut query_concat_offset: usize = 0;
    let context_room = |length: usize| if length == 0 { 0 } else { length + 1 };
    for q_idx in 0..seq_data.queries.len() {
        let q_len = letters.length(q_idx);
        let plus_len = if plus_on { q_len } else { 0 };
        let minus_len = if minus_on { q_len } else { 0 };
        let plus_masks = if plus_on {
            seq_data.query_masks[q_idx].clone()
        } else {
            Vec::new()
        };
        let minus_masks = if minus_on {
            reverse_mask_intervals(&seq_data.query_masks[q_idx], q_len)
        } else {
            Vec::new()
        };
        let plus_offset = query_concat_offset;
        let minus_offset = plus_offset + context_room(plus_len);

        query_base_offsets.push(plus_offset);
        query_contexts.push(QueryContext {
            query_idx: q_idx as u32,
            frame: 1,
            query_offset: plus_offset as i32,
            len: plus_len,
            masks: plus_masks.clone(),
        });
        query_contexts.push(QueryContext {
            query_idx: q_idx as u32,
            frame: -1,
            query_offset: minus_offset as i32,
            len: minus_len,
            masks: minus_masks.clone(),
        });

        query_context_masks.push(plus_masks);
        query_context_masks.push(minus_masks);
        query_concat_offset = minus_offset + context_room(minus_len);
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:245-250
    // ```c
    // Uint4
    // QueryInfo_GetSeqBufLen(const BlastQueryInfo* qinfo)
    // {
    //     BlastContextInfo * cinfo = & qinfo->contexts[qinfo->last_context];
    //     return cinfo->query_offset + cinfo->query_length + (cinfo->query_length ? 2 : 1);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:654
    // ```c
    //     BlastSeqBlkSetSequence(*seqblk, buf.release(), buflen - 2);
    // ```
    // The searched block leaves out the first and the last sentinel of the buffer. A batch
    // has a query with letters (`run_in_pool` stops on a batch without one).
    let query_concat_length = query_contexts.last().map_or(0, |last| {
        let length = last.len;
        (last.query_offset.max(0) as usize + length + if length == 0 { 1 } else { 2 })
            .saturating_sub(2)
    });
    let diag_container_query_length = diag_container_query_length(&query_contexts);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:166-174
    // ```c
    // query_info->contexts[0].query_offset = 0;
    // query_info->contexts[1].query_offset = kStrandLength + 1;
    // ```
    let query_context_offsets: Vec<i32> =
        query_contexts.iter().map(|ctx| ctx.query_offset).collect();
    let query_context_index = QueryContextIndex::new(&query_contexts, strand != QueryStrand::Both);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:498-604
    // ```c
    // EBlastEncoding encoding = GetQueryEncoding(prog);
    // ...
    // sequence = queries.GetBlastSequence(index, encoding, strand, eSentinels);
    // ...
    // memcpy(&buf.get()[offset], sequence.data.get(),
    //        sequence.length);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:500-603
    // ```c
    // int buflen = QueryInfo_GetSeqBufLen(qinfo);
    // TAutoUint1Ptr buf((Uint1*) calloc(buflen+1, sizeof(Uint1)));
    // ...
    // memcpy(&buf.get()[offset], sequence.data.get(), sequence.length);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:120-121
    // ```c
    // const Uint1 kProtSentinel = NULLB;
    // const Uint1 kNuclSentinel = 0xF;
    // ```
    // The buffer is the one BLASTNA copy of the query contexts; the searched block is the
    // buffer without its first and last sentinel, and a context is a slice of it.
    // EXPERIMENT (LOSAT_X_QPREPPAR): with the switch on and a search pool of two or more threads
    // the queries are encoded by the pool, each into its own part of the buffer; the buffer is
    // the same (`x_query_prep_par.rs`).
    // Diagnostic (LOSAT_TIMING): the batch's query encoding for the common line
    // (utils/xtiming.rs). No NCBI counterpart.
    let t_encode = crate::utils::xtiming::clock();
    // `-strand plus` or `minus` builds the buffer serially from the contexts that are searched
    // (`build_strand_query_blastna_buffer`); `both` is the code of the search without the option.
    let encoded_query_buffer = if strand == QueryStrand::Both {
        x_query_prep_par::query_buffer(
            &letters,
            &query_context_offsets,
            query_concat_length,
            parallel_pool,
        )
    } else {
        build_strand_query_blastna_buffer(
            &letters,
            &query_context_offsets,
            query_concat_length,
            strand,
        )
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Encode, t_encode);
    let encoded_query_concat_blastna_with_sentinels: &[u8] = &encoded_query_buffer;
    let encoded_query_concat_blastna: &[u8] = &encoded_query_buffer[1..query_concat_length + 1];
    let encoded_queries_blastna: Vec<&[u8]> = query_contexts
        .iter()
        .map(|ctx| {
            let start = ctx.query_offset.max(0) as usize;
            // The offset of an empty context can be past the block.
            if ctx.len == 0 {
                &[][..]
            } else {
                &encoded_query_concat_blastna[start..start + ctx.len]
            }
        })
        .collect();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:609-619
    // ```c
    // mask_at_hash = SBlastFilterOptionsMaskAtHash(filter_options);
    // ...
    // if (!mask_at_hash) {
    //     BlastSetUp_MaskQuery(query_blk, query_info, filter_maskloc,
    //                          program_number);
    // }
    // ```
    // LOSATN currently exposes the NCBI soft-query-masking path: DUST/lcase
    // masks restrict lookup-table entries through query_context_masks below,
    // while the encoded traceback/scoring sequence remains unmodified.

    // Finalize configuration with query-dependent parameters (adaptive lookup table selection)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:46-47
    // ```c
    // BlastChooseNaLookupTable(const LookupTableOptions* lookup_options,
    //                          Int4 approx_table_entries, Int4 max_q_off,
    //                          Int4 *lut_width)
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:62-67
    // ```c
    //    /* Discontiguous megablast must always use a megablast table */
    //
    //    if (lookup_options->mb_template_length > 0) {
    //         *lut_width = lookup_options->word_size;
    //         return eMBLookupTable;
    //    }
    // ```
    let discontig_template = config.mb_template_length > 0;
    let scoring_spec = NuclScoringSpec {
        reward: config.reward,
        penalty: config.penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
    };
    // EXPERIMENT (LOSAT_X_QPREPPAR): with the switch on and a search pool of two or more threads
    // the contexts' blocks are computed by the pool, in context order; the blocks are the same
    // (`x_query_prep_par.rs`).
    let context_ungapped =
        x_query_prep_par::ungapped_blocks(&encoded_queries_blastna, &scoring_spec, parallel_pool);
    let (approx_table_entries, max_q_off) =
        compute_lookup_query_stats(&query_contexts, &query_context_masks);
    let max_query_length = max_q_off.saturating_add(1).max(1);
    finalize_task_config(
        &mut config,
        approx_table_entries,
        max_query_length,
        discontig_template,
    );

    if args.verbose {
        eprintln!(
            "[INFO] Adaptive lookup: approx_entries={}, max_q_off={}, lut_word_length={}, two_stage={}, scan_step={}",
            approx_table_entries,
            max_q_off,
            config.lut_word_length,
            config.use_two_stage,
            config.scan_step
        );
    }
    // The word extension of NCBI's small-query lookup table reads a compressed query in
    // which ambiguity codes are bases (`build_compressed_query`); the other tables compare
    // the query letters themselves.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:237-240
    // ```c
    //     /* compute a compressed representation of the query, used
    //        for computing ungapped extensions */
    //
    //     BlastCompressBlastnaSequence(query);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1806-1819
    // ```c
    //     else if (lookup_wrap->lut_type == eSmallNaLookupTable) {
    //         ...
    //         if (lut->lut_word_length == lut->word_length)
    //             lut->extend_callback = (void *)s_BlastNaExtendDirect;
    //         else if (lut->lut_word_length % COMPRESSION_RATIO == 0 &&
    //                  lut->scan_step % COMPRESSION_RATIO == 0 &&
    //                  lut->word_length - lut->lut_word_length <= 4)
    //             lut->extend_callback = (void *)s_BlastSmallNaExtendAlignedOneByte;
    //         else
    //             lut->extend_callback = (void *)s_BlastSmallNaExtend;
    //     }
    // ```
    let small_na_compressed_query =
        if config.small_na_lookup && config.effective_word_size > config.lut_word_length {
            build_compressed_query(&encoded_query_concat_blastna)
        } else {
            Vec::new()
        };
    let small_na_aligned_one_byte = config.lut_word_length % COMPRESSION_RATIO == 0
        && config.scan_step % COMPRESSION_RATIO == 0
        && config.effective_word_size - config.lut_word_length <= 4;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1002-1003
    // ```c
    //     if ((status = BlastExtendWordNew(query->length, word_params,
    //                                     &aux_struct->ewp)) != 0)
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
    // ```c
    //                 diag_array_length = 1;
    //                 /* What power of 2 is just longer than the query? */
    //                 while (diag_array_length < (qlen+window_size))
    //                 {
    //                         diag_array_length = diag_array_length << 1;
    //                 }
    //                 /* These are used in the word finders to shift and mask
    //                 rather than dividing and taking the remainder. */
    //                 diag_table->diag_array_length = diag_array_length;
    //                 diag_table->diag_mask = diag_array_length-1;
    // ```
    // `query->length` is the whole query block (both strands and the sentinel, the
    // offsets of the word finder), not one strand: a table sized from one strand made the
    // diagonals of the two strands share entries and lost HSPs (S07+ independent audit).
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:133-142
    // ```c
    //         ewp->hash_table->window = word_params->options->window_size;
    //         ewp->hash_table->offset = word_params->options->window_size;
    // ...
    //         Boolean multiple_hits = (word_params->options->window_size > 0);
    //         ...
    //             s_BlastDiagTableNew(query_length, multiple_hits,
    //                               word_params->options->window_size);
    // ```
    // The two-hit window of the task or of `-window_size` (`configure_task`; a window whose sum
    // with the query would pass 2^30 is refused before this search, `check_batch_window`).
    let window_size = config.window_size;
    let diag_array_length = {
        let mut diag_array_length = 1usize;
        while diag_array_length < query_concat_length + window_size {
            diag_array_length <<= 1;
        }
        diag_array_length
    };
    let diag_mask = diag_array_length - 1;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2758-2792
    // ```c
    //       Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
    //       sbp->sfp[context] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
    //       BlastScoreFreqCalc(sbp, sbp->sfp[context], rfp, stdrfp);
    //       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
    //       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    //       if (loop_status) {
    //           contexts[context].is_valid = FALSE;
    // ```
    // Every query context (strand) has the ungapped block of its composition (kbp_std:
    // the gap trigger and the ungapped X-drop) and a gapped block (kbp_gap: e-values,
    // bit scores, cutoffs and the length adjustment); `scoring.rs` has the NCBI
    // references. NCBI raises an unsupported scoring system when it sets up the first
    // batch with a valid query, after the outfmt 0 prolog.
    let megablast = args.task == "megablast";
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:513-529
    // ```c
    //     if (scoring_options->gapped_calculation) { ... s_BlastScoreBlkGappedFill ... }
    // ```
    // An ungapped (`-ungapped`) search has the ungapped blocks only (`context_blocks_ungapped`).
    let context_blocks_result = if UNGAPPED {
        Ok(context_blocks_ungapped(&context_ungapped, &scoring_spec))
    } else {
        context_blocks(&context_ungapped, &scoring_spec)
    };
    let (context_karlin, round_down_evalue_score) = match context_blocks_result {
        Ok(blocks) => blocks,
        Err(message) => {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/setup_factory.cpp:170-172
            // ```c
            //     Blast_Message2TSearchMessages(blast_msg.Get(), query_info, search_messages);
            //     if (status != 0 &&
            //         (!(blast_msg.Get()) || (blast_msg.Get() && blast_msg.Get()->severity == eBlastSevError)))
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2784-2791
            // ```c
            //           contexts[context].is_valid = FALSE;
            //     ...
            //              Blast_MessageWrite(blast_message, eBlastSevWarning, context,
            //              kBlastErrMsg_CantCalculateUngappedKAParams);
            // ```
            // NCBI raises the error only when the first message of the batch is it: an
            // invalid context puts its warning first, and NCBI then goes on without the
            // gapped blocks (it crashes; with a first batch of invalid queries only, a later
            // batch raises it). LOSAT reproduces the error only for a first batch of valid
            // queries (the error does not depend on the queries, so the first batch with a valid
            // query meets it; NCBI writes the results and warnings of the batches before it).
            // The contexts of a query without letters are invalid from their set-up and put
            // no message here (`context_ungapped_blocks`).
            if query_contexts
                .iter()
                .zip(&context_ungapped)
                .any(|(context, block)| context.len != 0 && block.is_none())
            {
                anyhow::bail!(
                    "these scoring options have no Karlin-Altschul values, and the query batch (up to {} residues) has an invalid query; NCBI BLAST+ goes on without the gapped values and crashes, which is not supported by LOSAT's BLASTN",
                    batching.batch_size
                );
            }
            // A batch after the first: every batch before it had only invalid queries (one
            // with a valid query would have raised the error), and NCBI has written their
            // reports (`run_in_pool`).
            if !first_batch {
                return Err(ErrorAfterInvalidBatches(karlin_error(
                    &message,
                    seq_data.queries.len(),
                ))
                .into());
            }
            // The prolog of a file sink (the others were written before the search).
            let (outputs, output_formats) = outputs.expect(FIRST_BATCH_OUTPUTS);
            for (format, &output_format) in outputs.formats.iter_mut().zip(output_formats) {
                if prolog_with_report(format, output_format) {
                    write_format_prolog(format, output_format, megablast, &shown)?;
                }
            }
            return Err(karlin_error(&message, seq_data.queries.len()));
        }
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:324-332
    // ```c
    //        if (!(query_info->contexts[context].is_valid)) {
    //           /* either this context was never valid, or it was
    //              valid at the beginning of the search but is not
    //              valid now. The latter means that ungapped
    //              alignments can still occur to this context,
    //              so we set the cutoff score to be infinite */
    //           curr_cutoffs->cutoff_score = INT4_MAX;
    //           continue;
    //        }
    // ```
    // An invalid context has an infinite cutoff (below), so it gets no hit. Its slot holds
    // the blocks of a valid context (or placeholder blocks when none is valid), which no
    // hit reads.
    let context_valid: Vec<bool> = context_karlin.iter().map(Option::is_some).collect();
    let any_valid_context = context_valid.iter().any(|&valid| valid);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:908-939
    // ```c
    //    if (seq_src) {
    //       total_length = BlastSeqSrcGetTotLenStats(seq_src);
    //       if (total_length <= 0)
    //           total_length = BlastSeqSrcGetTotLen(seq_src);
    //       ...
    //       if (total_length > 0) {
    //           num_seqs = BlastSeqSrcGetNumSeqsStats(seq_src);
    //           if (num_seqs <= 0)
    //               num_seqs = BlastSeqSrcGetNumSeqs(seq_src);
    //       } else {
    //           /* Not a database search; each subject sequence is considered
    //              individually */
    //           Int4 oid = 0;  /* Get length of first sequence. */
    //           if ( (total_length = BlastSeqSrcGetSeqLen(seq_src, (void*) &oid)) < 0) {
    //               total_length = -1;
    //               num_seqs = -1;
    //           }
    //           num_seqs = 1;
    //       }
    //    }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:969-980
    // ```c
    //    if (sbp->gbp) {
    //        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
    //        if (Blast_SubjectIsTranslated(program_number)) {
    //            min_subject_length/=3;
    //        }
    //    } else {
    //        min_subject_length = (Int4) (total_length/num_seqs);
    //    }
    //
    //    if(min_subject_length <=0) {
    // 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
    //    }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_message.c:215-219
    // ```c
    //     case BLASTERR_SUBJECT_LENGTH_INVALID:
    //         new_msg->message = strdup("The average subject length is too short");
    //         new_msg->severity = eBlastSevFatal;
    //         new_msg->context = context;
    //         break;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:310-314
    // ```c
    //             retval = CPrelimSearchRunner(*m_InternalData, opts_memento.get())();
    //             if (retval) {
    //                 NCBI_THROW(CBlastException, eCoreBlastError,
    //                            BlastErrorCode2String(retval));
    //             }
    // ```
    // The subject set has no BLAST database statistics, and BLASTN has no Gumbel block
    // (blast_setup.c:475-480): a batch with a valid context stops NCBI when the subjects
    // have fewer letters than subjects (subjects without letters count with none). A batch
    // without a valid context is not searched (`CheckInternalData`), so it does not.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:448-458
    // ```c++
    // static Int4
    // s_SeqDbGetSeqLen(void* seqdb_handle, void* args)
    // {
    //     Int4* oid = (Int4*) args;
    //
    //     if (!seqdb_handle || !oid)
    //        return BLAST_SEQSRC_ERROR;
    //
    //     CSeqDB & seqdb = **(TSeqDBData *) seqdb_handle;
    //     return seqdb.GetSeqLength(*oid);
    // }
    // ```
    // A database without sequences (makeblastdb builds none; a volume whose index holds
    // none) has a total length of 0, so the setup above asks its reader for the length of
    // ordinal id 0, which it does not have: NCBI's search stage reports the reader's error
    // as a database file that cannot be read (local_blast.cpp:240-251). Shown with the
    // oracle (campaign dbsearch-20261012, `$BUILD_ROOT/dbsearch/db4/smoke/empty`: the
    // nucleotide fixture's index with no sequence): the outfmt 0 prolog, then `BLAST engine
    // error: Database memory map file error`, exit status 3, for outfmt 0, 6 and 7.
    if any_valid_context {
        if let SubjectSource::Database(database) = subjects {
            if database.metadata().db_num_seqs == 0 {
                return Err(database_search_error(
                    database::PreliminaryReadError("Database memory map file error").into(),
                    first_batch,
                    outputs,
                    megablast,
                    &shown,
                ));
            }
        }
    }
    if any_valid_context {
        let total_length = seq_data.db_len_total as i64;
        let (total_length, num_seqs) = if total_length > 0 {
            (total_length, seq_data.db_num_seqs as i64)
        } else {
            // A database with a total length of 0 holds sequences without letters only:
            // the length of its first sequence (`BlastSeqSrcGetSeqLen` of ordinal id 0) is 0.
            (
                subjects.prepared().map_or(0, |prepared| {
                    if prepared.records.is_empty() {
                        0
                    } else {
                        prepared.letters().len(0)
                    }
                }) as i64,
                1,
            )
        };
        if total_length / num_seqs <= 0 {
            let error = anyhow::Error::from(crate::cli::NativeError {
                exit: 3,
                message: "BLAST engine error: The average subject length is too short\n"
                    .to_string(),
            });
            // A batch after the first: the batches before had no valid query (one with a
            // valid query would have raised the error), and NCBI has written their reports
            // (`run_in_pool`).
            if !first_batch {
                return Err(ErrorAfterInvalidBatches(error).into());
            }
            let (outputs, output_formats) = outputs.expect(FIRST_BATCH_OUTPUTS);
            write_pairwise_prologs(outputs, output_formats, megablast, &shown, false)?;
            return Err(error);
        }
    }
    // Gap costs that the tables accept but LOSAT's greedy extension does not (`scoring.rs`);
    // with no valid context, NCBI extends nothing.
    if any_valid_context && !UNGAPPED {
        check_greedy_gap_costs(&args)?;
    }
    let search_karlin: Vec<ContextKarlin> = {
        let fallback = context_karlin
            .iter()
            .flatten()
            .next()
            .copied()
            .unwrap_or(ContextKarlin {
                ungapped: KarlinParams::default(),
                gapped: KarlinParams::default(),
            });
        context_karlin
            .iter()
            .map(|blocks| blocks.unwrap_or(fallback))
            .collect()
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-846
    // ```c
    // BLAST_ComputeLengthAdjustment(..., query_length, db_length, db_num_seqs, &length_adjustment);
    // effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
    // if (effective_db_length <= 0) effective_db_length = 1;
    // effective_search_space = effective_db_length * (query_length - length_adjustment);
    // query_info->contexts[index].eff_searchsp = effective_search_space;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-721
    // ```c
    //    /* use overriding value from effective lengths options or the real value
    //       from effective lengths parameters. */
    //    if (eff_len_options->db_length > 0)
    //       db_length = eff_len_options->db_length;
    //    else
    //       db_length = eff_len_params->real_db_length;
    // ```
    // A `-dbsize` above 0 replaces the total length in the effective lengths (with `-db`
    // and with `-subject`); the number of sequences, the average subject length above, the
    // batch mixer and the reports keep the real numbers.
    let db_len_total_i64 = if args.db_length() > 0 {
        args.db_length()
    } else {
        seq_data.db_len_total as i64
    };
    let db_num_seqs_i64 = seq_data.db_num_seqs as i64;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:778-786
    // ```c
    //       Int8 effective_search_space =
    //           s_GetEffectiveSearchSpaceForContext(eff_len_options, index,
    //                                               blast_message);
    //     ...
    //       if (query_info->contexts[index].is_valid &&
    //           ((query_length = query_info->contexts[index].query_length) > 0) ) {
    // ```
    // An invalid context keeps the search space of the options, which is 0, except in a query
    // chunk: NCBI sets the options to the search spaces of the batch's contexts
    // (`SplitQuery_SetEffectiveSearchSpace`), and a chunk's context takes the one at its own
    // index in the chunk, not its context in the batch (0, for an invalid context of the
    // batch, has the chunk compute its own).
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:177-182
    // ```c
    //     vector<Int8> eff_searchsp;
    //     for (size_t index = 0; index <= (size_t)qinfo->last_context; index++) {
    //         eff_searchsp.push_back(calc.GetEffSearchSpaceForContext(index));
    //     }
    //     options->SetEffectiveSearchSpace(eff_searchsp);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:690-693
    // ```c
    //     } else if (eff_len_options->num_searchspaces > 1) {
    //         ASSERT(context_index < eff_len_options->num_searchspaces);
    //         retval = eff_len_options->searchsp_eff[context_index];
    //     } else {
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:826-846
    // ```c
    //          if (effective_search_space == 0) {
    //          ...
    //              effective_search_space = effective_db_length *
    //                              (query_length - length_adjustment);
    //          }
    //       }
    //       query_info->contexts[index].eff_searchsp = effective_search_space;
    // ```
    let query_eff_searchsp: Vec<i64> = query_contexts
        .iter()
        .zip(&search_karlin)
        .zip(&context_valid)
        .enumerate()
        .map(|(index, ((ctx, blocks), &valid))| {
            let options_searchsp = match stage {
                BatchStage::Search { .. } => args.options_searchsp(),
                BatchStage::ChunkPrelim {
                    batch_eff_searchsp, ..
                } => batch_eff_searchsp[index],
            };
            if !valid || options_searchsp != 0 {
                return Ok(options_searchsp);
            }
            context_lengths(
                ctx.len as i64,
                db_len_total_i64,
                db_num_seqs_i64,
                &blocks.gapped,
            )
            .map(|(_, effective_search_space)| effective_search_space)
        })
        .collect::<Result<_>>()?;

    // NCBI searches a batch of at least two query chunks chunk by chunk in its preliminary
    // stage, after it has set up the whole batch (the errors and warnings above are the
    // batch's), and runs the traceback on the batch. A batch without a valid context is not
    // searched.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:232-236
    // ```c
    //     CRef<CQuerySplitter> query_splitter = setup_data->m_QuerySplitter;
    //     if (query_splitter->IsQuerySplit()) {
    //
    //         CRef<CSplitQueryBlk> split_query_blk = query_splitter->Split();
    // ```
    let split_prelim_lists = match stage {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:60-61
        // ```c
        //     if (!options->GetGappedMode()) m_NumChunks = 1;
        // ```
        BatchStage::Search { query_offsets, .. } if any_valid_context && !UNGAPPED => {
            search_query_chunks(
                args,
                &letters,
                &seq_data.query_masks,
                query_offsets,
                &query_eff_searchsp,
                &query_context_offsets,
                subjects,
                reborrow_outputs(&mut outputs),
                parallel_pool,
                batching,
                first_batch,
                x_keep,
            )?
        }
        _ => None,
    };
    // EXPERIMENT (LOSAT_X_BNQDROP): the encoding and the query chunks were the last readers of the
    // batch's letters; records held by the caller alone drop them here, before the lookup tables
    // are built (`x_query_drop.rs`). The search reads only the records' count from here on.
    let SequenceData {
        queries: _,
        query_ids,
        query_masks,
        db_len_total,
        db_num_seqs,
    } = seq_data;
    query_records.after_encoding();
    let seq_data = SequenceData {
        queries: query_records.records(),
        query_ids,
        query_masks,
        db_len_total,
        db_num_seqs,
    };

    // Build lookup tables
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1270-1306
    // ```c
    // pv_size = (Int4)(mb_lt->hashsize >> PV_ARRAY_BTS);
    // mb_lt->pv_array_bts = ilog2(mb_lt->hashsize / pv_size);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1312-1325
    // ```c
    // if (lookup_options->db_filter) {
    //    s_FillPV(query, location, mb_lt, lookup_options);
    //    s_ScanSubjectForWordCounts(seqsrc, mb_lt, counts,
    //                               lookup_options->max_db_word_count);
    // }
    // ```
    let subjects_for_lookup: &[FastaRecord] = subject_records.as_deref().unwrap_or(&[]);
    // A split batch has no lookup table of its own; its chunks have theirs.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_aux_priv.cpp:206-207
    // ```c
    //     // 5. Create the lookup table
    //     if ( !retval->m_QuerySplitter->IsQuerySplit() ) {
    // ```
    let (lookup_tables, scan_step) = if split_prelim_lists.is_some() {
        (
            LookupTables {
                two_stage_lookup: None,
                pv_direct_lookup: None,
                na_lookup: None,
            },
            config.scan_step,
        )
    } else {
        // Diagnostic (LOSAT_TIMING): the batch's lookup table for the common line
        // (utils/xtiming.rs). No NCBI counterpart.
        let t_build_lookup = crate::utils::xtiming::clock();
        let tables = build_lookup_tables(
            &config,
            &args,
            &encoded_queries_blastna,
            &query_context_masks,
            &query_context_offsets,
            subjects_for_lookup,
            subject_packed_cache_ref,
            approx_table_entries,
            x_keep,
            Some(subjects.metadata().db_len_total as u64),
        );
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::BuildLookup, t_build_lookup);
        tables
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:991-1041
    // ```c
    // Int4 offset_array_size = GetOffsetArraySize(lookup_wrap);
    // ...
    // aux_struct->offset_pairs =
    //   (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
    // ```c
    // switch (lookup->lut_type) {
    // case eMBLookupTable:
    //    offset_array_size = OFFSET_ARRAY_SIZE +
    //       ((BlastMBLookupTable*)lookup->lut)->longest_chain;
    //    break;
    // case eNaLookupTable:
    //    offset_array_size = OFFSET_ARRAY_SIZE +
    //       ((BlastNaLookupTable*)lookup->lut)->longest_chain;
    //    break;
    // ...
    // default:
    //    offset_array_size = OFFSET_ARRAY_SIZE;
    //    break;
    // }
    // ```
    let offset_array_size = if let Some(two_stage) = lookup_tables.two_stage_lookup.as_ref() {
        OFFSET_ARRAY_SIZE + two_stage.longest_chain()
    } else if let Some(pv_direct) = lookup_tables.pv_direct_lookup.as_ref() {
        OFFSET_ARRAY_SIZE + pv_direct.longest_chain()
    } else if let Some(na_lookup) = lookup_tables.na_lookup.as_ref() {
        OFFSET_ARRAY_SIZE + na_lookup.longest_chain()
    } else {
        OFFSET_ARRAY_SIZE
    };

    // NCBI reference: blast_traceback.c:654 (cutoff_score_min computed per subject)
    // Cutoff scores are computed per subject in the main loop; no global map needed.

    if args.verbose {
        eprintln!("Searching...");
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:478-501
    // ```c
    // while (TRUE) {
    //     status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
    //                                    dbseq_chunk_overlap);
    //     if (status == SUBJECT_SPLIT_DONE) break;
    //     if (status == SUBJECT_SPLIT_NO_RANGE) continue;
    //     ...
    // }
    // ```
    let t_search_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1409
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
    let subject_count = seq_data.db_num_seqs as u64;
    let progress_bar = if subjects_parallel {
        let bar = ProgressBar::new(subject_count);
        bar.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
                .unwrap(),
        );
        Some(bar)
    } else {
        None
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2960-2968
    // ```c
    // if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //    hitlist_size = args[kArgMaxTargetSequences].AsInteger();
    // }
    // m_NumDescriptions = hitlist_size;
    // m_NumAlignments = hitlist_size;
    // ```
    let hitlist_size = args.report_counts.hitlist_size;
    // NCBI blast_args.cpp:317-318: if (... args[kArgMaxHSPsPerSubject]) {
    // opt.SetMaxHspsPerSubject(args[kArgMaxHSPsPerSubject].AsInteger()); }
    // CLI omission maps to the engine's unlimited sentinel.
    let max_hsps_per_subject = args.max_hsps_per_subject.unwrap_or(0);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:43-70 (GetPrelimHitlistSize)
    // ```c
    // Int4
    // GetPrelimHitlistSize(Int4 hitlist_size, Int4 compositionBasedStats,
    //                      Boolean gapped_calculation)
    // {
    //     ...
    //     else if (gapped_calculation) {
    //          prelim_hitlist_size = MIN(MAX(2 * prelim_hitlist_size, 10),
    //                                   prelim_hitlist_size + 50);
    //     }
    //     return prelim_hitlist_size;
    // }
    // ```
    let prelim_hitlist_size =
        usize::try_from(get_prelim_hitlist_size(hitlist_size, false, !UNGAPPED))
            .map_err(|_| anyhow::anyhow!("the preliminary hit list size overflows"))?;

    // Debug mode: set BLEMIR_DEBUG=1 to enable, BLEMIR_DEBUG_WINDOW="q_start-q_end,s_start-s_end" to focus on a region
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();
    // BLASTN-specific debug mode: set LOSAT_DEBUG_BLASTN=1 to enable detailed hit loss diagnostics
    let blastn_debug = std::env::var("LOSAT_DEBUG_BLASTN").is_ok();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:350-692
    // ```c
    // ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
    // ...
    // Blast_HSPListPurgeHSPsWithCommonEndpoints(..., FALSE);
    // ...
    // Blast_IntervalTreeReset(tree);
    // ```
    let blastn_trace_enabled = blastn_trace::enabled();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3267-3282
    // ```c
    // arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
    //                  true);
    // ...
    // m_DebugOutput = static_cast<bool>(args["verbose"]);
    // ```
    let verbose = args.verbose;
    let debug_window: Option<(usize, usize, usize, usize)> =
        std::env::var("BLEMIR_DEBUG_WINDOW").ok().and_then(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                let q_parts: Vec<usize> =
                    parts[0].split('-').filter_map(|x| x.parse().ok()).collect();
                let s_parts: Vec<usize> =
                    parts[1].split('-').filter_map(|x| x.parse().ok()).collect();
                if q_parts.len() == 2 && s_parts.len() == 2 {
                    return Some((q_parts[0], q_parts[1], s_parts[0], s_parts[1]));
                }
            }
            None
        });

    // NCBI BLAST does NOT use mask_array/mask_hash for diagonal suppression
    // REMOVED: disable_mask variable (no longer needed)

    if debug_mode {
        // Build marker to verify correct code is running
        eprintln!(
            "[DEBUG] BLEMIR build: 2024-12-24-v7 (adaptive banding with MAX_WINDOW_SIZE=50000)"
        );
        eprintln!(
            "[DEBUG] Task: {}, Scoring: reward={}, penalty={}, gap_open={}, gap_extend={}",
            args.task, config.reward, config.penalty, config.gap_open, config.gap_extend
        );
        if let Some((q_start, q_end, s_start, s_end)) = debug_window {
            eprintln!(
                "[DEBUG] Focusing on window: query {}-{}, subject {}-{}",
                q_start, q_end, s_start, s_end
            );
        }
    }

    // Pass queries for use in the closure (the lookup tables are passed by the preliminary
    // stage's calls, `process_subject`'s `lookups`)
    let two_stage_lookup_ref = lookup_tables.two_stage_lookup.as_ref();
    let queries_ref = &seq_data.queries;
    let query_contexts_ref = &query_contexts;
    let query_base_offsets_ref = &query_base_offsets;
    let query_eff_searchsp_ref = &query_eff_searchsp;

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
    //    Int4 hspcnt; /**< Number of HSPs saved */
    //    ...
    // } BlastHSPList;
    // ```
    // The records whose titles name the subjects in LOSAT's debug and trace lines (none for a
    // database: its sequences are "unknown" there).
    let trace_subject_records: &[FastaRecord] =
        subjects.prepared().map_or(&[], |prepared| prepared.records);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2547-2556
    // ```c
    // const bool use_lcase_masks = args.Exist(kArgUseLCaseMasking)
    //     ? bool(args[kArgUseLCaseMasking])
    //     : kDfltArgUseLCaseMasking;
    // m_Scope = ReadSequencesToBlast(... use_lcase_masks, subjects, ...);
    // ```
    let lcase_masking = args.lcase_masking;

    // Capture config values for use in closure
    let effective_word_size = config.effective_word_size;
    let min_ungapped_score = config.min_ungapped_score;
    let use_dp = config.use_dp;
    let use_direct_lookup = config.use_direct_lookup;
    let reward = config.reward;
    let penalty = config.penalty;
    let gap_open = config.gap_open;
    let gap_extend = config.gap_extend;
    let search_karlin_ref = &search_karlin;
    let context_valid_ref = &context_valid;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:455-463
    // ```c
    //       double min_lambda = s_BlastFindSmallestLambda(sbp->kbp_gap, query_info, NULL);
    //       params->gap_x_dropoff = (Int4)
    //           (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
    //     ...
    //       params->gap_x_dropoff_final = (Int4)
    //           MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda, params->gap_x_dropoff);
    // ```
    let (x_drop_gapped, x_drop_final) =
        gap_x_dropoffs(&context_karlin, config.x_drop_gapped, config.x_drop_final);
    let scan_range = config.scan_range; // For off-diagonal hit detection
    let min_diag_separation = config.min_diag_separation; // For MB_HSP_CLOSE containment check

    // EXPERIMENT (LOSAT_X_BNFEWSUBJ, x_few_subj.rs): whether the seed stage's strips would take an
    // unmasked chunk of this `-subject` search that one subject's owner searches on this thread's
    // pool, and with it whether the search takes the shape of a one-subject search (`use_parallel`
    // false: the subjects one at a time). A batch split into query chunks has no lookup table
    // here (its traceback only), so it keeps its shape; each query chunk's call decides for itself.
    // `-subject` only: a database (`-db`) keeps its shapes (plan v4 line D).
    let x_strips = x_few_subj::on()
        && subjects.prepared().is_some()
        && x_pair_seed::takes_owner_chunks(
            two_stage_lookup_ref,
            diag_container_query_length <= 8000,
            window_size,
            scan_range,
            !small_na_compressed_query.is_empty(),
            debug_mode || blastn_debug || blastn_trace_enabled || debug_window.is_some(),
        );
    let x_kind = x_few_subj::Kind::of(
        two_stage_lookup_ref.is_some_and(|two_stage| two_stage.disc().is_some()),
        &args.task,
    );
    let x_serial_subjects = subjects_parallel
        && x_strips
        && !lcase_masking
        && subjects.prepared().is_some_and(|prepared| {
            let letters = prepared.letters();
            x_few_subj::serial_subjects(
                (0..prepared.records.len()).map(|index| letters.len(index)),
                num_threads,
                x_kind,
            )
        });
    if x_few_subj::on() && subjects_parallel && subjects.prepared().is_some() {
        crate::utils::xreach::x_reach!("LOSAT_X_BNFEWSUBJ", x_serial_subjects);
    }
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    let use_parallel = subjects_parallel && !x_serial_subjects;
    crate::utils::threading::report_stage("blastn", "subjects", seq_data.db_num_seqs, use_parallel);
    let db_len_total = seq_data.db_len_total;
    let db_num_seqs = seq_data.db_num_seqs;
    let evalue_threshold = super::super::coordination::determine_evalue(&args);
    let subject_besthit = args.subject_besthit;
    // EXPERIMENT (LOSAT_X_QGROUPPAR): whether the gapped stages of a subject searched on this
    // thread may run per query on the pool (x_qgroup_par.rs). Not on the diagnostic paths, which
    // print in the loop's order, and not while the pool searches subjects in parallel
    // (`use_parallel`; several subjects one at a time under LOSAT_X_BNFEWSUBJ take it).
    let x_qgroup_batch = !UNGAPPED
        && x_qgroup_par::on()
        && !use_parallel
        && parallel_pool.enabled()
        && num_threads > 1
        && !verbose
        && !debug_mode
        && !blastn_debug
        && !blastn_trace_enabled
        && debug_window.is_none()
        && std::env::var_os("LOSAT_DEBUG_COORDS").is_none()
        && std::env::var_os("LOSAT_DEBUG_COORDS_START").is_none();
    if x_qgroup_par::on() && !x_qgroup_batch {
        crate::utils::xreach::x_reach!("LOSAT_X_QGROUPPAR", false);
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
    // ```c
    // return ((hsp->num_ident * 100.0 <
    //         align_length * hit_options->percent_identity) ||
    //         align_length < hit_options->min_hit_length) ;
    // ```
    let percent_identity = args.percent_identity_value();
    let min_hit_length = args.min_hit_length;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1052-1127 (BlastScoreBlkNuclMatrixCreate)
    let score_matrix = build_blastna_matrix(reward, penalty);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:236-259
    // ```c
    // for (i = 0; i < 256; i++) {
    //    Int4 score = 0;
    //    if (i & 3) score += penalty; else score += reward;
    //    if ((i >> 2) & 3) score += penalty; else score += reward;
    //    if ((i >> 4) & 3) score += penalty; else score += reward;
    //    if (i >> 6) score += penalty; else score += reward;
    //    table[i] = score;
    // }
    // ```
    let nucl_score_table = build_nucl_score_table(reward, penalty);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:218-221
    // ```c
    // p->cutoffs[context].x_dropoff_init =
    //     (Int4)(sbp->scale_factor *
    //            ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
    // ```
    // NCBI reference: c++/src/algo/blast/api/blast_options_local_priv.cpp:49-54
    // m_InitWordOpts.Reset((BlastInitialWordOptions*)calloc(1, sizeof(BlastInitialWordOptions)));
    // NCBI reference: c++/src/algo/blast/api/blast_nucl_options.cpp:163-174
    // SetInitialWordOptionsDefaults() { SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_NUCL); ... }
    // SetMBInitialWordOptionsDefaults() { SetWindowSize(BLAST_WINDOW_SIZE_NUCL); }
    // NCBI reference: c++/src/algo/blast/api/disc_nucl_options.cpp:66-73
    // SetMBInitialWordOptionsDefaults() { SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_NUCL); ... }
    // Traditional megablast retains the zero-initialized X-drop option (0 bits, so 0 here):
    // its raw X-drop is the subject's word cutoff in ParametersUpdate below. `-xdrop_ungap`
    // replaces the task's bits (`configure_task`).
    let x_drop_ungapped_bits = config.x_drop_ungapped;
    let x_dropoff_init: Vec<i32> = search_karlin
        .iter()
        .map(|blocks| {
            if x_drop_ungapped_bits == 0.0 {
                0
            } else {
                ((x_drop_ungapped_bits * NCBIMATH_LN2) / blocks.ungapped.lambda).ceil() as i32
            }
        })
        .collect();
    let x_dropoff_init_ref = &x_dropoff_init;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:478-536
    // ```c
    // while (TRUE) {
    //     status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
    //                                    dbseq_chunk_overlap);
    //     if (status == SUBJECT_SPLIT_DONE) break;
    //     if (status == SUBJECT_SPLIT_NO_RANGE) continue;
    //     ...
    //     if (aux_struct->WordFinder) {
    //         aux_struct->WordFinder(...);
    //         if (init_hitlist->total == 0) continue;
    //     }
    //     ...
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:301-310
    // ```c
    // Int4 CLocalBlast::GetNumExtensions()
    // {
    //     Int4 retv = 0;
    //     if (m_InternalData) {
    //         BlastDiagnostics * diag = m_InternalData->m_Diagnostics->GetPointer();
    //         if (diag && diag->ungapped_stat) {
    //              retv = diag->ungapped_stat->good_init_extends;
    //         }
    //     }
    //     return retv;
    // }
    // ```
    // The initial hits that the word finder saved in every subject chunk of the batch; they
    // size the next batch.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:336-377
    // ```c
    // if (!gapped_calculation || sbp->matrix_only_scoring) {
    //     ...
    // } else {
    //     new_cutoff = gap_trigger;
    // }
    // new_cutoff = MIN(new_cutoff,
    //                  hit_params->cutoffs[context].cutoff_score_max);
    // curr_cutoffs->cutoff_score = new_cutoff;
    // ```
    // NCBI updates these per subject (`BlastInitialWordParametersUpdate`); LOSAT's values do not
    // depend on the subject, so they are computed once for the batch instead of in every subject
    // (`process_subject`): the cutoff scores, the hit-saving cutoff scores, the X-dropoffs and the
    // reduced cutoff scores of each context.
    let x_batch_cutoffs = {
        // The subject length that NCBI's word cutoff of an ungapped search takes: the average
        // length of the subjects (a `-subject` set's `(Uint4)(total / count)`, a database's
        // `(Int4)(total / MAX(1, count))`); the dbscan mode of the command line holds the total
        // of a `-subject` set, so `BlastInitialWordParametersUpdate` runs once for the batch.
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1434-1445,
        // c++/src/algo/blast/api/seqsrc_multiseq.cpp:245-266, seqsrc_seqdb.cpp:203-209
        let x_avg_subject_length: u64 = if UNGAPPED {
            match subjects {
                SubjectSource::Database(_) => ungapped::average_subject_length_database(
                    seq_data.db_len_total as u64,
                    seq_data.db_num_seqs as u64,
                ),
                SubjectSource::Prepared(_) => ungapped::average_subject_length_subject_set(
                    seq_data.db_len_total as u64,
                    seq_data.db_num_seqs as u64,
                ),
            }
        } else {
            0
        };
        let mut cutoff_scores: Vec<i32> = Vec::with_capacity(query_contexts.len());
        let mut hit_saving_cutoff_scores: Vec<i32> = Vec::with_capacity(query_contexts.len());
        // NCBI reference: c++/src/app/blast/blast_app_util.cpp:206-211;
        // c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181;
        // c++/src/algo/blast/core/blast_engine.c:1434-1445
        // db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, true));
        // if (dbscan_mode) { ... m_iTotalLength += (Int8) (*iter)->length; }
        // if (db_length == 0) { BLAST_OneSubjectUpdateParameters(...); }
        // The CLI subject set has a nonzero total length. Its context search
        // space is retained across subjects, just as for output statistics.
        for context_idx in 0..query_contexts.len() {
            // An invalid context has an infinite cutoff and, unset, a zero X-drop and
            // reduced cutoff (blast_parameters.c:324-332 above).
            if !context_valid_ref[context_idx] {
                cutoff_scores.push(i32::MAX);
                hit_saving_cutoff_scores.push(i32::MAX);
                continue;
            }
            // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:925-946
            // searchsp = query_info->contexts[context].eff_searchsp;
            // BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
            // params->cutoffs[context].cutoff_score_max = new_cutoff;
            // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:876-882
            // ```c
            //    if (options->cutoff_score > 0) {
            //       Int4 new_cutoff = options->cutoff_score * (Int4) sbp->scale_factor;
            //       for (context = query_info->first_context;
            //                           context <= query_info->last_context; ++context) {
            //          params->cutoffs[context].cutoff_score = new_cutoff;
            //          params->cutoffs[context].cutoff_score_max = new_cutoff;
            // ```
            // `-min_raw_gapped_score` replaces the e-value's cutoff (the scale factor of
            // BLASTN is 1); the word cutoff, the reduced cutoff and megablast's ungapped
            // X-drop below follow it.
            let hit_saving_cutoff = match config.min_raw_gapped_score {
                Some(score) => score,
                None => cutoff_score_max_from_evalue(
                    evalue_threshold,
                    query_eff_searchsp[context_idx],
                    &search_karlin_ref[context_idx].gapped,
                ),
            };
            // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:342-374
            // gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda);
            // new_cutoff = gap_trigger;
            // new_cutoff *= (Int4)sbp->scale_factor;
            // new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
            let gap_trigger = gap_trigger_raw_score(
                GAP_TRIGGER_BIT_SCORE_NUCL,
                &search_karlin_ref[context_idx].ungapped,
            );
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:346-374
            // ```c
            //       if (!gapped_calculation || sbp->matrix_only_scoring) {
            //          ...
            //          BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp,
            //                        MIN((Uint8)subj_length, (Uint8)query_length)*((Uint8)subj_length),
            //                        TRUE, gap_decay_rate);
            //       } else {
            //          new_cutoff = gap_trigger;
            //       }
            //       new_cutoff *= (Int4)sbp->scale_factor;
            //       new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
            // ```
            let cutoff = if UNGAPPED {
                ungapped::word_cutoff_score(
                    query_contexts[context_idx].len,
                    x_avg_subject_length,
                    &search_karlin_ref[context_idx].ungapped,
                    hit_saving_cutoff,
                )
                .map_err(|limit| ungapped_cutoff_limit_error(limit))?
            } else {
                cutoff_score_for_ungapped_extension(gap_trigger, hit_saving_cutoff, 1.0)
            };
            cutoff_scores.push(cutoff);
            hit_saving_cutoff_scores.push(hit_saving_cutoff);
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:380-383
        // ```c
        // if (curr_cutoffs->x_dropoff_init == 0)
        //    curr_cutoffs->x_dropoff = new_cutoff;
        // else
        //    curr_cutoffs->x_dropoff = curr_cutoffs->x_dropoff_init;
        // ```
        let mut x_dropoff_scores: Vec<i32> = Vec::with_capacity(cutoff_scores.len());
        for ((cutoff, &x_dropoff_init), &valid) in cutoff_scores
            .iter()
            .zip(x_dropoff_init_ref)
            .zip(context_valid_ref)
        {
            let x_dropoff = if !valid {
                0
            } else if x_dropoff_init == 0 {
                *cutoff
            } else {
                x_dropoff_init
            };
            x_dropoff_scores.push(x_dropoff);
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:408-412
        // ```c
        // curr_cutoffs->reduced_nucl_cutoff_score = (Int4)(0.8 * new_cutoff);
        // ```
        let mut reduced_cutoff_scores: Vec<i32> = Vec::with_capacity(cutoff_scores.len());
        for (cutoff, &valid) in cutoff_scores.iter().zip(context_valid_ref) {
            reduced_cutoff_scores.push(if valid {
                (0.8 * (*cutoff as f64)) as i32
            } else {
                0
            });
        }
        (
            cutoff_scores,
            hit_saving_cutoff_scores,
            x_dropoff_scores,
            reduced_cutoff_scores,
        )
    };
    let good_init_extends = std::sync::atomic::AtomicU64::new(0);
    let process_subject = |s_idx: usize,
                           subject: BlastnSubjectRef<'_>,
                           gap_scratch: &mut GapAlignScratch,
                           subject_scratch: &mut SubjectScratch,
                           lookups: Option<&LookupTables>,
                           prelim_source: Option<&Vec<PrelimHit>>,
                           subject_hits: &mut Option<Vec<BlastnHsp>>,
                           prelim_out: &mut Vec<PrelimHit>| {
        // The batch's lookup tables: the preliminary stage's calls pass them, the traceback's
        // (`prelim_source` set) pass `None` and do not scan (LOSAT_X_BNLUTDROP, `x_lut_drop.rs`).
        let two_stage_lookup_ref = lookups.and_then(|tables| tables.two_stage_lookup.as_ref());
        let pv_direct_lookup_ref = lookups.and_then(|tables| tables.pv_direct_lookup.as_ref());
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:131-156
        // ```c
        // typedef struct BlastNaLookupTable {
        //     ...
        // } BlastNaLookupTable;
        // ```
        let na_lookup_ref = lookups.and_then(|tables| tables.na_lookup.as_ref());
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:177-180
        // ```c
        //     int status = m_PrelimSearch->CheckInternalData();
        //     if (status != 0)
        //     {
        //          // Search was not run, but we send back an empty CSearchResultSet.
        // ```
        // With no valid context, NCBI does not search.
        if !any_valid_context {
            return;
        }
        let queries = queries_ref;
        let query_contexts = query_contexts_ref;
        let query_base_offsets = query_base_offsets_ref;
        let query_eff_searchsp = query_eff_searchsp_ref;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:493-526
        // ```c
        // if (aux_struct->WordFinder) {
        //     aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
        //                            word_params, aux_struct->ewp,
        //                            aux_struct->offset_pairs,
        //                            kScanSubjectOffsetArraySize,
        //                            init_hitlist, ungapped_stats);
        //     if (init_hitlist->total == 0) continue;
        // }
        // if (aux_struct->GetGappedScore) {
        //     status = aux_struct->GetGappedScore(program_number, query,
        //             query_info, subject, gap_align, score_params, ext_params,
        //             hit_params, word_params, init_hitlist, &hsp_list,
        //             gapped_stats, NULL);
        // }
        // ```
        let timing_ref = timing.as_deref();
        let timing_enabled = timing_ref.is_some();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
        // ```c
        // hsp_list = Blast_HSPListFree(hsp_list);
        // BlastInitHitListReset(init_hitlist);
        // ```
        let s_len_full = subject.length;
        // EXPERIMENT (LOSAT_X_BNTBRANGE): set when a traceback alignment met a fence.
        let x_fenced = subject.x_fenced;
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to. */
        //    BlastHSP** hsp_array;
        //    Int4 hspcnt;
        //    ...
        // } BlastHSPList;
        // ```
        // Made only for the debug and trace lines, which alone read it.
        let s_id_trace = if blastn_debug || blastn_trace_enabled {
            trace_subject_records.get(s_idx).map_or_else(
                || "unknown".to_string(),
                |record| subject_trace_id(record.title()),
            )
        } else {
            String::new()
        };
        let s_id = s_id_trace.as_str();

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1647-1662
        // ```c
        // scan_range[0] = 0;  /* subject seq mask index */
        // scan_range[1] = 0;  /* start pos of scan */
        // scan_range[2] = subject->length - lut_word_length;
        // if (subject->mask_type != eNoSubjMasking) {
        //     scan_range[1] = subject->seq_ranges[0].left + word_length - lut_word_length;
        //     scan_range[2] = subject->seq_ranges[0].right - lut_word_length;
        // }
        // ```
        // NCBI reference: ncbi-blast/c++/src/objtools/readers/fasta.cpp:856-874
        // ```c
        // case 'a': case 'b': case 'c': case 'd':
        // case 'g': case 'h':
        // ...
        //     char_type = eCharType_MaskedNonGap;
        //     break;
        // ```
        // NCBI reference: ncbi-blast/c++/src/objtools/readers/fasta.cpp:1079-1089
        // ```c
        // void CFastaReader::x_OpenMask(void)
        // {
        //     m_MaskRangeStart = GetCurrentPos(ePosWithGapsAndSegs);
        // }
        // void CFastaReader::x_CloseMask(void)
        // {
        //     m_CurrentMask->SetPacked_int().AddInterval(...);
        // }
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:812-826
        // ```c
        // if (subjects.GetMask(i).NotEmpty()) {
        //     s_SeqLoc2MaskedSubjRanges(...);
        //     BlastSeqBlkSetSeqRanges(..., eSoftSubjMasking);
        // }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2550-2556
        // ```c++
        //         const bool use_lcase_masks = args.Exist(kArgUseLCaseMasking)
        // 	    ? bool(args[kArgUseLCaseMasking])
        //             : kDfltArgUseLCaseMasking;
        //         CRef<blast::CBlastQueryVector> subjects;
        //         m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(),
        //                                        subj_range, parse_deflines,
        //                                        use_lcase_masks, subjects, m_IsMapper);
        // ```
        // The lowercase masks of the subjects are those of the FASTA records of
        // `-subject`; a database sequence has none (its masks would come from
        // `-db_soft_mask`, which LOSAT does not implement).
        let subject_masks = match (subject.letters, subject.store) {
            (Some(letters), _) if lcase_masking => collect_lowercase_masks(letters),
            (None, Some(store)) if lcase_masking => store.lowercase().to_vec(),
            _ => Vec::new(),
        };
        let subject_masked = !subject_masks.is_empty();
        // Unmasked, every chunk carries its range inline (`SeqRanges::Inline`) and no one reads
        // the soft ranges, so none are made (NCBI points them at `full_range`).
        let soft_ranges = if subject_masked {
            build_subject_seq_ranges_from_masks(subject_masks.as_slice(), s_len_full)
        } else {
            Vec::new()
        };

        if s_len_full < effective_word_size {
            return;
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
        // ```c
        // if (subj_is_na) {
        //     BlastSeqBlkSetSequence(subj, sequence.data.release(),
        //        ((sentinels == eSentinels) ? sequence.length - 2 :
        //         sequence.length));
        //     ...
        //     SBlastSequence compressed_seq =
        //         subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
        //                                   eNa_strand_plus, eNoSentinels);
        //     BlastSeqBlkSetCompressedSequence(subj,
        //                               compressed_seq.data.release());
        // }
        // ```
        // The packed letters of the subject set's record (`PreparedSubjects::packed`) or the
        // stored data of the database sequence (`BlastnSubjectRef::packed`).
        let s_seq_packed_full = subject.packed;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4019-4025
        // ```c
        //             status = BLAST_GreedyGappedAlignment(
        //                          query_tmp.sequence, subject->sequence,
        // ...
        //                          (Boolean) TRUE, FALSE, fence_hit);
        // ```
        // No NCBI counterpart: the one-byte copy of this subject's packed sequence
        // that the preliminary greedy alignment reads instead (made on first use, shared by every
        // chunk and worker of this subject). Off unless the search uses the
        // greedy aligner and this is the preliminary stage.
        let x_prelim_subject = x_greedy_prelim::XPrelimSubject::new(
            !use_dp && prelim_source.is_none(),
            s_seq_packed_full,
            s_len_full,
        );

        let dbseq_chunk_overlap = DBSEQ_CHUNK_OVERLAP;
        let mut split_state = SubjectSplitState::new(s_len_full, soft_ranges);

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:336-377
        // ```c
        // if (!gapped_calculation || sbp->matrix_only_scoring) {
        //     double cutoff_e = s_GetCutoffEvalue(program_number);
        //     Int4 query_length = query_info->contexts[context].query_length;
        //     if (program_number == eBlastTypeBlastn ||
        //         program_number == eBlastTypeMapping)
        //         query_length *= 2;
        //     BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp,
        //                   MIN((Uint8)subj_length,
        //                       (Uint8)query_length) * ((Uint8)subj_length),
        //                   TRUE, gap_decay_rate);
        // } else {
        //     new_cutoff = gap_trigger;
        // }
        // new_cutoff = MIN(new_cutoff,
        //                  hit_params->cutoffs[context].cutoff_score_max);
        // curr_cutoffs->cutoff_score = new_cutoff;
        // ```
        // NCBI uses UNGAPPED params (kbp_std) for gap_trigger calculation
        // and GAPPED params (kbp_gap) for cutoff_score_max calculation
        // The values depend on the context only (the gapped calculation's branch above: gap_trigger
        // capped by cutoff_score_max), so the batch computes them once (`x_batch_cutoffs`).
        let (cutoff_scores, hit_saving_cutoff_scores, x_dropoff_scores, reduced_cutoff_scores) = (
            &x_batch_cutoffs.0,
            &x_batch_cutoffs.1,
            &x_batch_cutoffs.2,
            &x_batch_cutoffs.3,
        );

        // BLASTN debug: Log cutoff scores for this query-subject pair
        if blastn_debug {
            eprintln!(
                "[BLASTN_DEBUG] Subject {} (len={}): cutoff_scores={:?}",
                s_id, s_len_full, cutoff_scores
            );
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
        // ```c
        // /* Delete if not done in last loop iteration to prevent memory leak. */
        // hsp_list = Blast_HSPListFree(hsp_list);
        //
        // BlastInitHitListReset(init_hitlist);
        // ```
        let collect_prelim_hits_for_chunk = |chunk: &SubjectChunk,
                                             soft_ranges: &[(i32, i32)],
                                             gap_scratch: &mut GapAlignScratch,
                                             subject_scratch: &mut SubjectScratch,
                                             reuse_prelim_hits: bool|
         -> Vec<PrelimHit> {
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
            let mut dbg_two_hit_failed = 0usize;
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
                // diagonal hash in `process_offset_pairs` below skip the walk of a bucket whose
                // cells all lie below the level bound, unless diagnostics print what they read.
                let x_lv_ref = diag_hash.x_lv_on()
                    && !(debug_enabled || blastn_trace_enabled || debug_window.is_some());
                if diag_hash.x_lv_on() {
                    crate::utils::xreach::x_reach!("LOSAT_X_DIAGLV", use_diag_hash && x_lv_ref);
                    // EXPERIMENT (LOSAT_X_DIAGLV2): the reference fold below never reads the memo
                    // (the strips stage of dc-megablast records its own side).
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
                    // Dispatch point of PAIRPAR: `x_pair_seed::run` runs this block (scan,
                    // word extension, diagonal hash, ungapped extension; one-hit diagonal-hash
                    // scope only) with the subject scan in strips on the pool threads and the
                    // pairs folded by this thread in the scan order; `continue` ends the
                    // one-iteration strand loop. Otherwise the code below runs.
                    if x_pair_seed::run(
                        &x_pair_seed::Inputs {
                            two_stage,
                            search_seq_packed,
                            s_len,
                            subject_seq_ranges: &subject_seq_ranges,
                            subject_masked,
                            scan_step,
                            query_context_index: &query_context_index,
                            query_contexts,
                            encoded_query_concat_blastna: &encoded_query_concat_blastna,
                            encoded_query_concat_blastna_with_sentinels:
                                &encoded_query_concat_blastna_with_sentinels,
                            cutoff_scores: &cutoff_scores,
                            x_dropoff_scores: &x_dropoff_scores,
                            reduced_cutoff_scores: &reduced_cutoff_scores,
                            score_matrix: &score_matrix,
                            nucl_score_table: &nucl_score_table,
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
                            diagnostics: debug_enabled
                                || blastn_trace_enabled
                                || debug_window.is_some(),
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
                    let small_na_word =
                        (!small_na_compressed_query.is_empty()).then(|| SmallNaWord {
                            compressed_query: &small_na_compressed_query,
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
                    let offset_pairs = subject_scratch.offset_pairs.as_mut_slice();
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
                    let mut process_offset_pairs =
                        |offset_pairs: &[OffsetPair],
                         offset_pairs_len: &mut usize,
                         scan_pause_ns: &mut u64,
                         track_scan_pause: bool| {
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
                                let query_idx = ctx.query_idx as usize;
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
                                if q_off0 + lut_word_length > q_context_end
                                    || kmer_start + lut_word_length > s_len
                                {
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
                                    (kmer_start as isize + diag_array_length - q_off0 as isize)
                                        & diag_mask
                                } else {
                                    kmer_start as isize - q_off0 as isize
                                };

                                // Check if this seed is in the debug window
                                let s_pos = kmer_start + chunk.offset;
                                let in_window =
                                    if let Some((q_start, q_end, s_start, s_end)) = debug_window {
                                        q_pos_usize >= q_start
                                            && q_pos_usize <= q_end
                                            && s_pos >= s_start
                                            && s_pos <= s_end
                                    } else {
                                        false
                                    };

                                if debug_enabled && in_window {
                                    dbg_window_seeds += 1;
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
                                let (q_ext_start, s_ext_start, word_ext_left, word_ext_right) =
                                    if let Some(small_na_word) = &small_na_word {
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
                                            small_na_word.extend(
                                                q_off0,
                                                kmer_start,
                                                q_context_start,
                                                q_context_end,
                                                s_range,
                                            )
                                        };
                                        let Some((q_word, s_word)) = word else {
                                            continue;
                                        };
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
                                        let max_ext_left =
                                            ext_to.min(kmer_start).min(q_off0.saturating_add(1));
                                        let mut ext_left = 0usize;
                                        let mut q_left = q_off0 + 1;
                                        let mut s_off = kmer_start;
                                        let mut s_idx = s_off / COMPRESSION_RATIO;
                                        while ext_left < max_ext_left {
                                            if debug_enabled {
                                                dbg_left_ext_iters += 1;
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
                                            let s_byte =
                                                unsafe { *search_seq_packed.get_unchecked(s_idx) };
                                            let s_base = ((s_byte
                                                << (2 * (s_off % COMPRESSION_RATIO)))
                                                >> 6)
                                                as u8;
                                            // SAFETY: q_left is bounded by the outer-sentinel query buffer.
                                            let q_base = unsafe {
                                                *encoded_query_concat_blastna_with_sentinels
                                                    .get_unchecked(q_left)
                                            };
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
                                                continue;
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
                                            if q_right + (ext_to - ext_left)
                                                > encoded_query_concat_blastna_with_sentinels.len()
                                            {
                                                continue;
                                            }
                                            let mut s_idx = s_off / COMPRESSION_RATIO;

                                            // NCBI BLAST: for (; ext_right < ext_to - ext_left; ++ext_right)
                                            // Reference: na_ungapped.c:1128-1136
                                            while ext_right < (ext_to - ext_left) {
                                                if debug_enabled {
                                                    dbg_right_ext_iters += 1;
                                                }
                                                // NCBI BLAST: if (base mismatch) break;
                                                // Reference: na_ungapped.c:1129-1131
                                                // SAFETY: s_idx tracks s_off / COMPRESSION_RATIO within bounds by s_len check above.
                                                let s_byte = unsafe {
                                                    *search_seq_packed.get_unchecked(s_idx)
                                                };
                                                let s_base = ((s_byte
                                                    << (2 * (s_off % COMPRESSION_RATIO)))
                                                    >> 6)
                                                    as u8;
                                                // SAFETY: q_right is bounded by the outer-sentinel query buffer.
                                                let q_base = unsafe {
                                                    *encoded_query_concat_blastna_with_sentinels
                                                        .get_unchecked(q_right)
                                                };
                                                if s_base != q_base {
                                                    break;
                                                }
                                                ext_right += 1;
                                                q_right += 1;
                                                s_off += 1;
                                                if s_off % COMPRESSION_RATIO == 0 {
                                                    s_idx += 1;
                                                }
                                            }

                                            // NCBI BLAST: if (ext_left + ext_right < ext_to) continue;
                                            // Reference: na_ungapped.c:1125
                                            if ext_left + ext_right < ext_to {
                                                // Word_length match failed, skip this seed
                                                continue;
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
                                    (s_ext_start as isize + diag_array_length
                                        - q_ext_start as isize)
                                        & diag_mask
                                } else {
                                    s_ext_start as isize - q_ext_start as isize
                                };
                                let diag_idx = if use_array_indexing { diag as usize } else { 0 };
                                let (last_hit, hit_saved) =
                                    if use_array_indexing && diag_idx < diag_array_size {
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
                                        packed_kmer_at_seed_mask(
                                            search_seq_packed,
                                            s_pos,
                                            lut_word_length,
                                        ),
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
                                        // NCBI reference: na_ungapped.c:658, 692
                                        // Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
                                        // if (Delta < 0) Delta = 0;
                                        let delta_calc =
                                            window_size as isize - word_length as isize;
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
                                        let s_a = s_off_pos as isize + word_length as isize
                                            - window_size as isize;
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
                                                        let off_entry =
                                                            &hit_level_array[off_diag_idx];
                                                        (
                                                            off_entry.last_hit,
                                                            hit_len_array[off_diag_idx],
                                                        )
                                                    } else {
                                                        (0, 0)
                                                    }
                                                };
                                                // NCBI: off_s_end - delta >= s_a (signed comparison)
                                                // Convert to signed for comparison to match NCBI behavior
                                                if off_s_l > 0
                                                    && (off_s_end as isize - delta) >= s_a
                                                    && (off_s_end as isize - off_s_l as isize)
                                                        <= s_b
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
                                                        let off_entry =
                                                            &hit_level_array[off_diag_idx];
                                                        (
                                                            off_entry.last_hit,
                                                            hit_len_array[off_diag_idx],
                                                        )
                                                    } else {
                                                        (0, 0)
                                                    }
                                                };
                                                // NCBI: off_s_end >= s_a (signed comparison)
                                                // Convert to signed for comparison to match NCBI behavior
                                                if off_s_l > 0
                                                    && (off_s_end as isize) >= s_a
                                                    && (off_s_end as isize - off_s_l as isize
                                                        + delta)
                                                        <= s_b
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
                                                if let Some((off_s_end, off_s_l, _)) =
                                                    diag_hash.retrieve(diag_i32 + delta as i32)
                                                {
                                                    if off_s_l > 0
                                                        && (off_s_end as isize - delta) >= s_a
                                                        && (off_s_end as isize - off_s_l as isize)
                                                            <= s_b
                                                    {
                                                        off_found = true;
                                                        break;
                                                    }
                                                }
                                                if let Some((off_s_end, off_s_l, _)) =
                                                    diag_hash.retrieve(diag_i32 - delta as i32)
                                                {
                                                    if off_s_l > 0
                                                        && (off_s_end as isize) >= s_a
                                                        && (off_s_end as isize - off_s_l as isize
                                                            + delta)
                                                            <= s_b
                                                    {
                                                        off_found = true;
                                                        break;
                                                    }
                                                }
                                            }
                                        }

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
                                    dbg_ungapped_ext_calls += 1;
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
                                        &score_matrix,
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
                                        &nucl_score_table,
                                        reduced_cutoff,
                                        &score_matrix,
                                    )
                                };
                                if let Some(ungapped_start) = ungapped_start {
                                    let elapsed_ns = ungapped_start.elapsed().as_nanos() as u64;
                                    if let Some(timing) = timing_ref {
                                        timing.ungapped_ns.fetch_add(
                                            elapsed_ns,
                                            std::sync::atomic::Ordering::Relaxed,
                                        );
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
                                    let trace_qs =
                                        ungapped.q_start as isize - q_context_start as isize;
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
                                        dbg_ungapped_low += 1;
                                    }
                                    if in_window && debug_mode {
                                        eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}", q_pos_usize, s_pos, ungapped_score, min_ungapped_score);
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
                                        hit_level_array[diag_idx].flag =
                                            if hit_ready { 1 } else { 0 };
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
                                    dbg_gapped_attempted += 1;
                                }

                                if in_window && debug_mode {
                                    eprintln!("[DEBUG WINDOW] Seed at q={}, s={} -> COLLECT UNGAPPED (score={}, len={})", q_pos_usize, s_pos, ungapped_score, qe - qs);
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
                                        seed_q_off: q_off - q_context_start,
                                        seed_s_off: s_off,
                                        qs,
                                        qe,
                                        ss,
                                        se: ungapped_se,
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
                                hit_ready = true;
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
                        };

                    let mut on_word =
                        |kmer_start: usize,
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
                                dbg_total_s_positions += 1;
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
                                process_offset_pairs(
                                    offset_pairs,
                                    &mut offset_pairs_len,
                                    &mut scan_pause_ns,
                                    true,
                                );
                            }

                            // For each match, add to the offset_pairs buffer.
                            two_stage.for_each_hit(current_lut_kmer, |q_off_1| {
                                if debug_enabled {
                                    dbg_seeds_found += 1;
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
                                        dbg_seeds_found += 1;
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
                                process_offset_pairs(
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
                            let mut x_on_word_csr = |(kmer_start, s_range, index): (
                                usize,
                                usize,
                                u64,
                            )| {
                                if debug_enabled {
                                    dbg_total_s_positions += 1;
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
                                    process_offset_pairs(
                                        offset_pairs,
                                        &mut offset_pairs_len,
                                        &mut scan_pause_ns,
                                        true,
                                    );
                                }
                                for &q_off_1 in run {
                                    if debug_enabled {
                                        dbg_seeds_found += 1;
                                    }
                                    offset_pairs[offset_pairs_len] = OffsetPair {
                                        q_off: q_off_1 as usize - 1,
                                        s_off: kmer_start,
                                        s_range,
                                    };
                                    offset_pairs_len += 1;
                                }
                                if offset_pairs_len >= offset_array_size {
                                    process_offset_pairs(
                                        offset_pairs,
                                        &mut offset_pairs_len,
                                        &mut scan_pause_ns,
                                        true,
                                    );
                                }
                            };
                            let mut ring: x_dc_csr::Ring<'_, (usize, usize, u64)> =
                                x_dc_csr::Ring::new(csr);
                            scan_subject_disc_words_with_ranges(
                                search_seq_packed,
                                s_len,
                                disc,
                                &subject_seq_ranges,
                                subject_masked,
                                |s_off, s_range, index, index2| {
                                    debug_assert!(index2.is_none());
                                    ring.push(
                                        (s_off, s_range, index),
                                        |w| w.2 as usize,
                                        &mut x_on_word_csr,
                                    )
                                },
                            );
                            ring.flush(&mut x_on_word_csr);
                        } else {
                            scan_subject_disc_words_with_ranges(
                                search_seq_packed,
                                s_len,
                                disc,
                                &subject_seq_ranges,
                                subject_masked,
                                |s_off, s_range, index, index2| {
                                    on_word(s_off, s_range, index, index2)
                                },
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
                        macro_rules! x_scan_csr {
                            ($head:expr, $pv:expr) => {{
                                use x_dc_csr::{XHead, XRing};
                                let head = $head;
                            let mut x_on_word_csr =
                                |(kmer_start, s_range, index): (usize, usize, u64)| {
                                    if offset_pairs_len >= OFFSET_ARRAY_SIZE {
                                        process_offset_pairs(
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
                                            dbg_seeds_found += 1;
                                        }
                                        offset_pairs[offset_pairs_len] = OffsetPair {
                                            q_off: q_off_1 as usize - 1,
                                            s_off: kmer_start,
                                            s_range,
                                        };
                                        offset_pairs_len += 1;
                                    }
                                    if offset_pairs_len >= offset_array_size {
                                        process_offset_pairs(
                                            offset_pairs,
                                            &mut offset_pairs_len,
                                            &mut scan_pause_ns,
                                            true,
                                        );
                                    }
                                };
                            let mut ring = head.ring::<(usize, usize, u64)>();
                            let mb_scan_kind =
                                select_mb_scan_kind(lut_word_length, scan_step, subject_masked);
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
                                        dbg_total_s_positions += 1;
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
                        let mb_scan_kind =
                            select_mb_scan_kind(lut_word_length, scan_step, subject_masked);
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
                        process_offset_pairs(
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
                    // Dispatch point of LOSAT_X_PAIRPARNA: `x_pair_seed::run_na` runs this block
                    // (scan, diagonal hash, ungapped extension; one-hit diagonal-hash scope only)
                    // with the subject scan in strips on the pool threads and the pairs folded in
                    // the scan order; `continue` ends the one-iteration strand loop. Otherwise
                    // the code below runs.
                    if x_pair_seed::run_na(
                        &x_pair_seed::NaInputs {
                            pv_direct_lookup: pv_direct_lookup_ref,
                            na_lookup: na_lookup_ref,
                            use_direct_lookup,
                            word_length: safe_k,
                            search_seq_packed,
                            s_len,
                            subject_seq_ranges: &subject_seq_ranges,
                            subject_masked,
                            scan_step,
                            query_context_index: &query_context_index,
                            query_contexts,
                            encoded_queries_blastna: &encoded_queries_blastna,
                            cutoff_scores: &cutoff_scores,
                            x_dropoff_scores: &x_dropoff_scores,
                            reduced_cutoff_scores: &reduced_cutoff_scores,
                            score_matrix: &score_matrix,
                            nucl_score_table: &nucl_score_table,
                            diag_offset,
                            diag_hash_window: diag_hash_insert_window(
                                window_size,
                                scan_range,
                                safe_k,
                            ),
                            use_array_indexing,
                            window_size,
                            scan_range,
                            diagnostics: debug_enabled
                                || blastn_trace_enabled
                                || debug_window.is_some(),
                            one_subject_owner: !use_parallel && reuse_prelim_hits,
                            timing: timing_ref,
                        },
                        diag_hash,
                        ungapped_hits,
                        &mut subject_scratch.x_pair_fold,
                    ) {
                        continue;
                    }
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
                        &subject_seq_ranges,
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
                                dbg_total_s_positions += 1;
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
                                    dbg_seeds_found += 1;
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
                                let query_idx = ctx.query_idx as usize;

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
                                    (kmer_start as isize + diag_array_length - q_off0 as isize)
                                        & diag_mask
                                } else {
                                    kmer_start as isize - q_off0 as isize
                                };

                                // Check if this seed is in the debug window
                                let s_pos = kmer_start + chunk.offset;
                                let in_window =
                                    if let Some((q_start, q_end, s_start, s_end)) = debug_window {
                                        q_pos_usize >= q_start
                                            && q_pos_usize <= q_end
                                            && s_pos >= s_start
                                            && s_pos <= s_end
                                    } else {
                                        false
                                    };

                                if debug_enabled && in_window {
                                    dbg_window_seeds += 1;
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
                                let (last_hit, hit_saved) =
                                    if use_array_indexing && diag_idx < diag_array_size {
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
                                let mut word_type = 1u8; // Default: single word (when word_length == lut_word_length)
                                let mut extended = 0usize;
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

                                    // NCBI reference: na_ungapped.c:852-881 - off-diagonal search (DiagHash version)
                                    if word_type == 1 {
                                        // NCBI reference: na_ungapped.c:858
                                        // Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
                                        // if (Delta < 0) Delta = 0;
                                        let delta_calc = window_size as isize - safe_k as isize;
                                        let delta_max = if delta_calc < 0 {
                                            0
                                        } else {
                                            scan_range.min(delta_calc as usize) as isize
                                        };

                                        // NCBI reference: na_ungapped.c:855-856
                                        // Int4 s_a = s_off_pos + word_length - window_size;
                                        // Int4 s_b = s_end_pos - 2 * word_length;
                                        // CRITICAL: NCBI uses signed arithmetic (Int4), so s_a and s_b can be negative
                                        // LOSAT must use signed arithmetic to match NCBI behavior
                                        let s_a = s_off_pos as isize + safe_k as isize
                                            - window_size as isize;
                                        let s_b = s_end_pos as isize - 2 * safe_k as isize;

                                        if use_array_indexing {
                                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:688-694
                                            // ```c
                                            // orig_diag = real_diag + diag_table->diag_array_length;
                                            // off_diag  = (orig_diag + delta) & diag_table->diag_mask;
                                            // ```
                                            let orig_diag = diag + diag_array_length;
                                            // NCBI reference: na_ungapped.c:859
                                            // for (delta = 1; delta <= Delta; ++delta) {
                                            for delta in 1..=delta_max {
                                                // NCBI reference: na_ungapped.c:860-871
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
                                                        let off_entry =
                                                            &hit_level_array[off_diag_idx];
                                                        (
                                                            off_entry.last_hit,
                                                            hit_len_array[off_diag_idx],
                                                        )
                                                    } else {
                                                        (0, 0)
                                                    }
                                                };
                                                if off_s_l > 0
                                                    && (off_s_end as isize - delta) >= s_a
                                                    && (off_s_end as isize - off_s_l as isize)
                                                        <= s_b
                                                {
                                                    off_found = true;
                                                    break;
                                                }

                                                // NCBI reference: na_ungapped.c:872-880
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
                                                        let off_entry =
                                                            &hit_level_array[off_diag_idx];
                                                        (
                                                            off_entry.last_hit,
                                                            hit_len_array[off_diag_idx],
                                                        )
                                                    } else {
                                                        (0, 0)
                                                    }
                                                };
                                                if off_s_l > 0
                                                    && (off_s_end as isize) >= s_a
                                                    && (off_s_end as isize - off_s_l as isize
                                                        + delta)
                                                        <= s_b
                                                {
                                                    off_found = true;
                                                    break;
                                                }
                                            }
                                        } else {
                                            let diag_i32 = diag as i32;
                                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:860-874
                                            // ```c
                                            // off_rc = s_BlastDiagHashRetrieve(hash_table, diag + delta,
                                            //           &off_s_end, &off_s_l, &off_hit_saved);
                                            // ...
                                            // off_rc = s_BlastDiagHashRetrieve(hash_table, diag - delta,
                                            //           &off_s_end, &off_s_l, &off_hit_saved);
                                            // ```
                                            for delta in 1..=delta_max {
                                                if let Some((off_s_end, off_s_l, _)) =
                                                    diag_hash.retrieve(diag_i32 + delta as i32)
                                                {
                                                    if off_s_l > 0
                                                        && (off_s_end as isize - delta) >= s_a
                                                        && (off_s_end as isize - off_s_l as isize)
                                                            <= s_b
                                                    {
                                                        off_found = true;
                                                        break;
                                                    }
                                                }
                                                if let Some((off_s_end, off_s_l, _)) =
                                                    diag_hash.retrieve(diag_i32 - delta as i32)
                                                {
                                                    if off_s_l > 0
                                                        && (off_s_end as isize) >= s_a
                                                        && (off_s_end as isize - off_s_l as isize
                                                            + delta)
                                                            <= s_b
                                                    {
                                                        off_found = true;
                                                        break;
                                                    }
                                                }
                                            }
                                        }

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
                                        &score_matrix,
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
                                        &nucl_score_table,
                                        reduced_cutoff,
                                        &score_matrix,
                                    )
                                };
                                if let Some(ungapped_start) = ungapped_start {
                                    let elapsed_ns = ungapped_start.elapsed().as_nanos() as u64;
                                    if let Some(timing) = timing_ref {
                                        timing.ungapped_ns.fetch_add(
                                            elapsed_ns,
                                            std::sync::atomic::Ordering::Relaxed,
                                        );
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
                                        hit_level_array[diag_idx].flag =
                                            if hit_ready { 1 } else { 0 };
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
                                        dbg_ungapped_low += 1;
                                    }
                                    if in_window && debug_mode {
                                        eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}", q_pos_usize, s_pos, ungapped_score, min_ungapped_score);
                                    }
                                    continue;
                                }

                                if debug_enabled {
                                    dbg_gapped_attempted += 1;
                                }

                                if in_window && debug_mode {
                                    eprintln!("[DEBUG WINDOW] Seed at q={}, s={} -> COLLECT UNGAPPED (score={}, len={})", q_pos_usize, s_pos, ungapped_score, qe - qs);
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
                                        seed_q_off: q_off,
                                        seed_s_off: s_off,
                                        qs,
                                        qe,
                                        ss,
                                        se: ungapped_se,
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
                                hit_ready = true;
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
                } // end of if/else for two-stage vs original lookup
            } // end of strand loop

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:556-563
            // ```c
            //       if (aux_struct->GetGappedScore) {
            //          status = aux_struct->GetGappedScore(program_number, query, query_info,
            //                        subject, gap_align, score_params, ext_params, hit_params,
            //                        init_hitlist, &hsp_list, gapped_stats, fence_hit);
            //       } else {
            //          status = BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
            //                        hit_params->options, &hsp_list);
            //       }
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4719-4775
            // ```c
            // for (index = 0; index < init_hitlist->total; ++index) {
            //     init_hsp = &init_hsp_array[index];
            //     ...
            //     Blast_HSPInit(ungapped_data->q_start, ungapped_data->q_start + length,
            //                   ungapped_data->s_start, ungapped_data->s_start + length,
            //                   init_hsp->offsets.qs_offsets.q_off, init_hsp->offsets.qs_offsets.s_off,
            //                   context, query_frame, subject->frame, ungapped_data->score, NULL, &new_hsp);
            //     Blast_HSPListSaveHSP(hsp_list, new_hsp);
            // }
            // /* Sort by score: */
            // Blast_HSPListSortByScore(hsp_list);
            // ```
            // An ungapped (`-ungapped`) search sorts the initial hits by score first, as a gapped
            // search does (`BlastNaWordFinder` ends with `if (word_params->ungapped_extension)
            // Blast_InitHitListSortByScore(init_hitlist)`, na_ungapped.c:1691-1692, and
            // `ungapped_extension` is TRUE for every program but PHI: blast_parameters.c:189-192),
            // with the same comparator (`score_compare_ungapped_hits` = `score_compare_match`:
            // score descending, s_start, length descending, the concatenated query's q_start;
            // stable, as the oracle's glibc qsort). The order of ties (a plus and a minus hit
            // with the same score and context-relative ranges) then comes out as NCBI's, and
            // `Blast_HSPListSortByScore` below keeps it. Then the HSP list of the chunk is made
            // from those hits (no gapped extension), the subject chunk's offset is taken
            // (blast_engine.c:564-580), and there is no purge of common end points (a gapped
            // search's step).
            if UNGAPPED {
                // The same sort as the gapped branch below (review 2 S5): with LOSAT_X_BNMAPPED
                // (default ON) the scratch of the stable sort is a mapped buffer, not a plain
                // `sort_by` scratch of up to 8 MB left in the thread's arena. The comparator and
                // the stability are the same, so the order is `sort_by`'s. LOSAT_X_QGROUPPAR is
                // not used for an ungapped search.
                if x_bn_mapped() {
                    x_bn_mapped_reach(
                        ungapped_hits
                            .len()
                            .saturating_mul(std::mem::size_of::<UngappedHit>()),
                    );
                    crate::utils::x_mapped::x_sort_by_mapped(
                        ungapped_hits,
                        score_compare_ungapped_hits,
                        &X_SORT,
                    )
                } else {
                    ungapped_hits.sort_by(score_compare_ungapped_hits)
                }
                good_init_extends.fetch_add(
                    ungapped_hits.len() as u64,
                    std::sync::atomic::Ordering::Relaxed,
                );
                let mut chunk_prelim_hits: Vec<PrelimHit> = ungapped_hits
                    .iter()
                    .map(|uh| PrelimHit {
                        context_idx: uh.context_idx,
                        query_idx: uh.query_idx,
                        query_frame: uh.query_frame,
                        query_context_offset: uh.query_context_offset,
                        prelim_qs: uh.qs,
                        prelim_qe: uh.qe,
                        prelim_ss: uh.ss.saturating_add(chunk.offset),
                        prelim_se: uh.se.saturating_add(chunk.offset),
                        prelim_score: uh.score,
                        seed_qs: uh.seed_q_off,
                        seed_ss: uh.seed_s_off.saturating_add(chunk.offset),
                        prelim_evalue: 0.0,
                    })
                    .collect();
                sort_prelim_hits_by_score(&mut chunk_prelim_hits);
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

            // ============================================================
            // NCBI ARCHITECTURE: Batch process ungapped hits in score order
            // ============================================================
            // NCBI reference: blast_gapalign.c:3824 - ASSERT(Blast_InitHitListIsSortedByScore(init_hitlist))
            // NCBI processes init_hsp_array sorted by score DESCENDING
            // High-score HSPs are processed first, their GAPPED HSPs added to interval tree
            // Lower-score HSPs are then checked for containment against these GAPPED HSPs

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:306-310
            // ```c
            // void Blast_InitHitListSortByScore(BlastInitHitList * init_hitlist)
            // {
            //     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
            //           sizeof(BlastInitHSP), score_compare_match);
            // }
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:307-309
            // ```c
            // qsort(init_hsp_array, init_hitlist->total,
            //       sizeof(BlastInitHSP), score_compare_match);
            // ```
            // The oracle's qsort (glibc 2.39) is a stable merge sort: hits that
            // compare equal keep the order in which they were saved. `sort_by`
            // is stable too.
            // EXPERIMENT (LOSAT_X_QGROUPPAR): with the stage run per query, each query's hits are
            // taken in the order saved above and sorted by the task (the comparator and the sort
            // are those below, so each query's hits come in the relative order the batch's sort
            // gives them).
            let x_counts = (x_qgroup_batch && reuse_prelim_hits).then(|| {
                x_qgroup_par::counts(ungapped_hits.iter().map(|hit| hit.query_idx), queries.len())
            });
            let x_order = x_counts
                .as_deref()
                .and_then(|counts| x_qgroup_par::plan(counts, x_qgroup_par::Stage::Prelim));
            match (&x_order, &x_counts) {
                (Some(_), Some(counts)) => {
                    x_qgroup_par::group_in_place(ungapped_hits, |hit| hit.query_idx, counts)
                }
                // EXPERIMENT (LOSAT_X_BNMAPPED): `x_mapped::x_sort_by_mapped` is a stable sort
                // with the same comparator, so the order is that of `sort_by`.
                _ if x_bn_mapped() => {
                    x_bn_mapped_reach(
                        ungapped_hits
                            .len()
                            .saturating_mul(std::mem::size_of::<UngappedHit>()),
                    );
                    crate::utils::x_mapped::x_sort_by_mapped(
                        ungapped_hits,
                        score_compare_ungapped_hits,
                        &X_SORT,
                    )
                }
                _ => ungapped_hits.sort_by(score_compare_ungapped_hits),
            }

            // Debug counters for containment analysis
            let mut dbg_containment_skipped = 0usize;
            let total_ungapped = ungapped_hits.len();
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1688-1689
            // ```c
            //     Blast_UngappedStatsUpdate(ungapped_stats, total_hits, hits_extended,
            //                               init_hitlist->total);
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_diagnostics.c:106-112
            // ```c
            //    if (!ungapped_stats || total_hits == 0)
            //       return;
            //    ...
            //    ungapped_stats->good_init_extends += saved_hits;
            // ```
            good_init_extends
                .fetch_add(total_ungapped as u64, std::sync::atomic::Ordering::Relaxed);
            let gapped_start_time = std::time::Instant::now();
            let mut dbg_gapped_calls = 0usize;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2786-2794
            // ```c
            //    if (!compressed_subject) {
            //       s = subject + s_off;
            //       rem = 4; /* Special value to indicate that sequence is already
            //                   uncompressed */
            //    } else {
            //       s = subject + s_off/4;
            //       rem = s_off % 4;
            // ```
            // No NCBI counterpart: the bases of `s_seq_packed` one per byte, for the
            // greedy calls below; `None` where no score-only greedy call happens (no ungapped
            // hits, or the dynamic-programming extension).
            let x_prelim_s = if total_ungapped > 0 {
                x_prelim_subject.chunk(chunk.offset)
            } else {
                None
            };

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2959-2963,4012-4025,4041-4046
            // ```c
            //    Int4 x_dropoff =  gap_align->gap_x_dropoff;
            //
            //    if(init_hsp->ungapped_data->score < x_dropoff) {
            //           x_dropoff = init_hsp->ungapped_data->score;
            //    }
            // ...
            //          } else if (is_greedy) {
            //             if (init_hsp->ungapped_data) {
            //                 init_hsp->offsets.qs_offsets.q_off =
            //                     init_hsp->ungapped_data->q_start + init_hsp->ungapped_data->length/2;
            // ...
            //             status = BLAST_GreedyGappedAlignment(
            // ...
            //             if (s_end >= (Int4)init_hsp->offsets.qs_offsets.s_off + 8) {
            //                init_hsp->offsets.qs_offsets.s_off += 3;
            // ...
            //             status = s_BlastDynProgNtGappedAlignment(&query_tmp, subject,
            // ```
            // This closure is what each branch of the ordered loop below computes for hit `x_index`, as a
            // function of the index alone: the seed shift by 3, the X-drop limited by the ungapped score, and the
            // aligner call, for the dynamic-programming branch and the greedy branch. The branches below make
            // the same calls on the same arguments, so the closure gives the same value. The tree test and the
            // tree insertion are not part of it.
            // No NCBI counterpart: the gapped extension of the loop below
            // as a function of the hit's index (the two branches of the loop,
            // without their trace logging).
            let x_prelim_extend =
                |uh: &UngappedHit,
                 x_scratch: &mut XAheadScratch<'_>|
                 -> Option<(usize, usize, usize, usize, i32, usize, usize)> {
                    let q_seq_blastna = encoded_queries_blastna[uh.context_idx as usize];
                    if use_dp {
                        let mut seed_qs = uh.seed_q_off;
                        let mut seed_ss = uh.seed_s_off;
                        if uh.se >= uh.seed_s_off.saturating_add(8) {
                            seed_qs = seed_qs.saturating_add(3);
                            seed_ss = seed_ss.saturating_add(3);
                        }
                        let x_drop_score_only = x_drop_gapped.min(uh.score);
                        let (p_qs, p_qe, p_ss, p_se, p_score, _, _, _, _, _) =
                            extend_gapped_heuristic_with_scratch(
                                q_seq_blastna,
                                s_seq_packed,
                                s_len,
                                seed_qs,
                                seed_ss,
                                1,
                                reward,
                                penalty,
                                &score_matrix,
                                gap_open,
                                gap_extend,
                                x_drop_score_only,
                                &mut *x_scratch.gap,
                                true,
                            );
                        Some((p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss))
                    } else {
                        let ungapped_len = uh.qe.saturating_sub(uh.qs);
                        let seed_qs = uh.qs + ungapped_len / 2;
                        let seed_ss = uh.ss + ungapped_len / 2;
                        x_greedy_prelim::x_greedy_gapped_alignment_score_only(
                            x_prelim_s.expect(
                                "the one-byte copy exists wherever a score-only greedy call runs",
                            ),
                            q_seq_blastna,
                            s_len,
                            seed_qs,
                            seed_ss,
                            reward,
                            penalty,
                            gap_open,
                            gap_extend,
                            x_drop_gapped,
                            &mut *x_scratch.greedy,
                        )
                    }
                };
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850
            // ```c
            //    for (index=0; index<init_hitlist->total; index++)
            // ```
            // Look-ahead for the preliminary loop. When the conditions below do not hold,
            // `x_ahead` is `None` and the loop runs as ported.
            // Only the caller's own (single-chunk) loop gets helpers; chunk
            // batches already occupy the pool.
            // The preliminary gapped loop over one list of ungapped hits in score order: the
            // chunk's on the owner (with the look-ahead helpers, whose scratch
            // `x_helper_scratch` holds), or one query's (LOSAT_X_QGROUPPAR, no helpers). Returns
            // the alignments run and the hits the tree held.
            let x_prelim_list = |ungapped_hits: &[UngappedHit],
                                 interval_tree: &mut BlastIntervalTree,
                                 gap_scratch: &mut GapAlignScratch,
                                 greedy_align_scratch: &mut GreedyAlignScratch,
                                 x_helper_scratch: Option<
                &mut Vec<(GapAlignScratch, GreedyAlignScratch)>,
            >,
                                 prelim_hits: &mut Vec<PrelimHit>|
             -> (usize, usize) {
                let mut dbg_containment_skipped = 0usize;
                let mut dbg_gapped_calls = 0usize;
                let x_ungapped: &[UngappedHit] = ungapped_hits;
                let x_prelim_compute = |x_index: usize, x_scratch: &mut XAheadScratch<'_>| {
                    x_prelim_extend(&x_ungapped[x_index], x_scratch)
                };
                let x_ahead = match x_ahead_window {
                    Some(x_window)
                        if x_helper_scratch.is_some()
                            && reuse_prelim_hits
                            && !blastn_trace_enabled
                            && x_ungapped.len() >= 2 * x_window =>
                    {
                        Some(crate::utils::xahead::Ahead::new(x_window))
                    }
                    _ => None,
                };
                if x_ahead_window.is_some() {
                    crate::utils::xreach::x_reach!("LOSAT_X_AHEAD", x_ahead.is_some());
                }
                if x_greedy_spec && x_ahead_window.is_some() {
                    crate::utils::xreach::x_reach!("LOSAT_X_GREEDYSPEC", x_ahead.is_some());
                }
                let mut x_no_helpers = Vec::new();
                let x_helper_scratch = x_helper_scratch.unwrap_or(&mut x_no_helpers);
                if x_ahead.is_some() {
                    while x_helper_scratch.len() + 1 < num_threads {
                        x_helper_scratch.push((GapAlignScratch::new(), GreedyAlignScratch::new()));
                    }
                }
                let x_helpers = x_helper_scratch
                    .iter_mut()
                    .map(|(gap, greedy)| XAheadScratch { gap, greedy });

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,4087-4088
                // ```c
                //    for (index=0; index<init_hitlist->total; index++)
                // ...
                //             status = BlastIntervalTreeAddHSP(new_hsp, tree, query_info,
                //                                     eQueryAndSubject);
                // ```
                // `with_helpers` runs the closure on this thread. The closure is the port of this NCBI loop, in
                // NCBI order; the other threads of the pool only evaluate values ahead (`utils/xahead.rs`).
                // Process each ungapped hit in score order
                crate::utils::xahead::with_helpers(
                    x_ahead.as_ref(),
                    x_helpers,
                    &x_prelim_compute,
                    || {
                        for (idx, uh) in ungapped_hits.iter().enumerate() {
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3918-3919
                            // ```c
                            //    for (index=0; index<init_hitlist->total; index++)
                            // ...
                            //       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
                            //                                         hit_options->min_diag_separation))
                            // ```
                            // The owner offers later indexes to the helpers, leaving out the ones the tree already contains.
                            // The tree is only read here. A wrong guess costs time only, because the loop decides each index
                            // itself.
                            // No NCBI counterpart: offer the next hits to the helpers,
                            // leaving out those the tree already contains.
                            if let Some(x_ahead) = x_ahead.as_ref() {
                                x_ahead.announce(idx, x_ungapped.len(), |x_index| {
                                    let u = &x_ungapped[x_index];
                                    let h = TreeHsp {
                                        query_offset: u.qs as i32,
                                        query_end: u.qe as i32,
                                        subject_offset: u.ss as i32,
                                        subject_end: u.se as i32,
                                        score: u.score,
                                        query_frame: u.query_frame,
                                        query_length: query_contexts[u.context_idx as usize].len
                                            as i32,
                                        query_context_offset: u.query_context_offset,
                                        subject_frame_sign: 1,
                                    };
                                    interval_tree
                                        .containing_hsp(
                                            &h,
                                            u.query_context_offset,
                                            min_diag_separation,
                                        )
                                        .is_none()
                                });
                            }
                            // NCBI reference: blast_gapalign.c:3886-4090 (no per-HSP logging)
                            if verbose && idx % 100 == 0 {
                                eprintln!("[INFO] Gapped extension: {}/{}", idx, total_ungapped);
                            }
                            // Get query sequence for this hit (context-specific)
                            let ctx = &query_contexts[uh.context_idx as usize];
                            let q_seq_blastna = encoded_queries_blastna[uh.context_idx as usize];
                            let subject_len = s_len;

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2793
                            // ```c
                            // if (!compressed_subject) {
                            //    s = subject + s_off;
                            //    rem = 4;
                            // } else {
                            //    s = subject + s_off/4;
                            //    rem = s_off % 4;
                            // }
                            // ```
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:503-507 (traceback uses uncompressed subject)
                            // The dynamic-programming score-only call below reads the packed subject
                            // (`s_seq_score`); the greedy one runs on the one-byte copy (`x_prelim_s`),
                            // NCBI's `compressed_subject = FALSE` form.
                            let s_seq_score = s_seq_packed;

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3913
                            // ```c
                            // tmp_hsp.query.offset = q_start;
                            // tmp_hsp.query.end = q_end;
                            // tmp_hsp.query.frame = query_info->contexts[context].frame;
                            // tmp_hsp.subject.offset = s_start;
                            // tmp_hsp.subject.end = s_end;
                            // ```
                            // NCBI uses 0-based coordinates internally for the tree
                            let subject_frame_sign = 1i32;
                            let ungapped_tree_hsp = TreeHsp {
                                query_offset: uh.qs as i32,
                                query_end: uh.qe as i32,
                                subject_offset: uh.ss as i32,
                                subject_end: uh.se as i32,
                                score: uh.score,
                                query_frame: uh.query_frame,
                                query_length: ctx.len as i32,
                                query_context_offset: uh.query_context_offset,
                                subject_frame_sign,
                            };

                            // NCBI reference: blast_gapalign.c:3918 BlastIntervalTreeContainsHSP
                            // Check if UNGAPPED HSP is contained in existing GAPPED HSPs
                            let containing_hsp = interval_tree.containing_hsp(
                                &ungapped_tree_hsp,
                                uh.query_context_offset,
                                min_diag_separation,
                            );
                            let is_contained = containing_hsp.is_some();

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3919
                            // ```c
                            // tmp_hsp.query.offset = q_start;
                            // tmp_hsp.query.end = q_end;
                            // tmp_hsp.subject.offset = s_start;
                            // tmp_hsp.subject.end = s_end;
                            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
                            //                                   hit_options->min_diag_separation))
                            // ```
                            if blastn_trace_enabled
                                && blastn_trace::should_trace_range(
                                    "prelim",
                                    uh.context_idx,
                                    s_idx,
                                    s_id,
                                    uh.qs,
                                    uh.qe,
                                    uh.ss.saturating_add(chunk.offset),
                                    uh.se.saturating_add(chunk.offset),
                                    ctx.len,
                                    ctx.frame,
                                )
                            {
                                blastn_trace::log(
                            "prelim",
                            format!(
                                "subject={}({}) context={} ungapped=q{}..{} s{}..{} raw_score={} tree_contains={} containing_hsp={:?} min_diag_separation={}",
                                s_id,
                                s_idx,
                                uh.context_idx,
                                uh.qs,
                                uh.qe,
                                uh.ss.saturating_add(chunk.offset),
                                uh.se.saturating_add(chunk.offset),
                                uh.score,
                                is_contained,
                                containing_hsp,
                                min_diag_separation
                            ),
                        );
                            }

                            if is_contained {
                                // NCBI: Skip gapped extension if ungapped HSP is contained
                                dbg_containment_skipped += 1;
                                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3918-3921
                                // ```c
                                //       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
                                //                                         hit_options->min_diag_separation))
                                //       {
                                // ```
                                // When the tree contains the hit, NCBI skips the aligner call for this index. The helpers are
                                // told that nobody needs its value.
                                if let Some(x_ahead) = x_ahead.as_ref() {
                                    x_ahead.skip(idx);
                                }
                                continue;
                            }

                            // Select gapped-start seed within the ungapped HSP.
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4046
                            //
                            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3924-3927,4058
                            // ```c
                            // if (is_rpsblast)
                            //    cutoff = hit_params->cutoffs[rps_cutoff_index].cutoff_score;
                            // else
                            //    cutoff = hit_params->cutoffs[context].cutoff_score;
                            // ...
                            // if (gap_align->score >= cutoff) {
                            // ```
                            // The score-only gapped HSP enters the preliminary HSP list and
                            // interval tree only if it reaches the hit-saving cutoff, not
                            // the lower BlastInitialWordParameters cutoff used for saving
                            // ungapped HSPs.
                            let cutoff_score = hit_saving_cutoff_scores[uh.context_idx as usize];
                            let (
                                prelim_qs,
                                prelim_qe,
                                prelim_ss,
                                prelim_se,
                                prelim_score,
                                seed_qs,
                                seed_ss,
                            ) = if let Some(x_ahead) = x_ahead.as_ref() {
                                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4025,4041-4046
                                // ```c
                                //          } else if (is_greedy) {
                                //             if (init_hsp->ungapped_data) {
                                //                 init_hsp->offsets.qs_offsets.q_off =
                                //                     init_hsp->ungapped_data->q_start + init_hsp->ungapped_data->length/2;
                                // ...
                                //             status = BLAST_GreedyGappedAlignment(
                                // ...
                                //             if (s_end >= (Int4)init_hsp->offsets.qs_offsets.s_off + 8) {
                                //                init_hsp->offsets.qs_offsets.s_off += 3;
                                // ...
                                //             status = s_BlastDynProgNtGappedAlignment(&query_tmp, subject,
                                // ```
                                // Look-ahead: the value of the aligner call comes from a helper or is
                                // computed here with `x_prelim_compute`. Without helpers the `else if use_dp` branch and the greedy
                                // branch below port these calls.
                                // No NCBI counterpart: the value of the branches
                                // below, taken from a helper or evaluated here.
                                let mut x_scratch = XAheadScratch {
                                    gap: &mut *gap_scratch,
                                    greedy: &mut *greedy_align_scratch,
                                };
                                let x_value =
                                    match x_ahead.take(idx, &mut x_scratch, &x_prelim_compute) {
                                        Some(x_value) => x_value,
                                        None => x_prelim_compute(idx, &mut x_scratch),
                                    };
                                match x_value {
                                    Some(x_value) => x_value,
                                    None => continue,
                                }
                            } else if use_dp {
                                // DP seed selection (blastn)
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4033-4045
                                // ```c
                                // if (s_end >= (Int4)init_hsp->offsets.qs_offsets.s_off + 8) {
                                //    init_hsp->offsets.qs_offsets.s_off += 3;
                                //    init_hsp->offsets.qs_offsets.q_off += 3;
                                // }
                                // status = s_BlastDynProgNtGappedAlignment(...);
                                // ```
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4071-4076
                                // ```c
                                // status = Blast_HSPInit(...,
                                //        init_hsp->offsets.qs_offsets.q_off,
                                //        init_hsp->offsets.qs_offsets.s_off, ...);
                                // ```
                                let mut seed_qs = uh.seed_q_off;
                                let mut seed_ss = uh.seed_s_off;
                                if uh.se >= uh.seed_s_off.saturating_add(8) {
                                    seed_qs = seed_qs.saturating_add(3);
                                    seed_ss = seed_ss.saturating_add(3);
                                }

                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2959-2963 (x_dropoff limited by ungapped score)
                                let x_drop_score_only = x_drop_gapped.min(uh.score);
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2973-3003
                                // ```c
                                // offset_adjustment = COMPRESSION_RATIO -
                                //      (init_hsp->offsets.qs_offsets.s_off % COMPRESSION_RATIO);
                                // q_length = init_hsp->offsets.qs_offsets.q_off + offset_adjustment;
                                // s_length = init_hsp->offsets.qs_offsets.s_off + offset_adjustment;
                                // score_left = s_BlastAlignPackedNucl(query, subject, q_length, s_length, ...);
                                // score_right = s_BlastAlignPackedNucl(query+q_length-1,
                                //    subject+(s_length+3)/COMPRESSION_RATIO - 1, ...);
                                // ```
                                if blastn_trace_enabled
                                    && blastn_trace::should_trace_seed(
                                        "prelim",
                                        uh.context_idx,
                                        s_idx,
                                        s_id,
                                        seed_qs,
                                        seed_ss.saturating_add(chunk.offset),
                                    )
                                {
                                    let offset_adjustment =
                                        COMPRESSION_RATIO - (seed_ss % COMPRESSION_RATIO);
                                    let mut score_q_anchor =
                                        seed_qs.saturating_add(offset_adjustment);
                                    let mut score_s_anchor =
                                        seed_ss.saturating_add(offset_adjustment);
                                    if score_q_anchor > q_seq_blastna.len()
                                        || score_s_anchor > subject_len
                                    {
                                        score_q_anchor =
                                            score_q_anchor.saturating_sub(COMPRESSION_RATIO);
                                        score_s_anchor =
                                            score_s_anchor.saturating_sub(COMPRESSION_RATIO);
                                    }
                                    blastn_trace::log(
                                    "prelim",
                                    format!(
                                        "subject={}({}) context={} score_only_seed=({}, {}) local_seed_s={} offset_adjustment={} score_anchor=({}, {}) global_score_anchor_s={} x_drop_score_only={} ungapped_score={}",
                                        s_id,
                                        s_idx,
                                        uh.context_idx,
                                        seed_qs,
                                        seed_ss.saturating_add(chunk.offset),
                                        seed_ss,
                                        offset_adjustment,
                                        score_q_anchor,
                                        score_s_anchor,
                                        score_s_anchor.saturating_add(chunk.offset),
                                        x_drop_score_only,
                                        uh.score
                                    ),
                                );
                                }

                                // Preliminary DP gapped extension (score-only)
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:176-183
                                let (p_qs, p_qe, p_ss, p_se, p_score, _, _, _, _, _) =
                                    extend_gapped_heuristic_with_scratch(
                                        q_seq_blastna,
                                        s_seq_score,
                                        subject_len,
                                        seed_qs,
                                        seed_ss,
                                        1,
                                        reward,
                                        penalty,
                                        &score_matrix,
                                        gap_open,
                                        gap_extend,
                                        x_drop_score_only,
                                        gap_scratch,
                                        true,
                                    );
                                (p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss)
                            } else {
                                // Greedy seed selection (megablast)
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4017
                                let ungapped_len = uh.qe.saturating_sub(uh.qs);
                                let seed_qs = uh.qs + ungapped_len / 2;
                                let seed_ss = uh.ss + ungapped_len / 2;

                                let prelim =
                                    match x_greedy_prelim::x_greedy_gapped_alignment_score_only(
                                        x_prelim_s.expect("the one-byte copy exists wherever a score-only greedy call runs"),
                                        q_seq_blastna,
                                        subject_len,
                                        seed_qs,
                                        seed_ss,
                                        reward,
                                        penalty,
                                        gap_open,
                                        gap_extend,
                                        x_drop_gapped,
                                        greedy_align_scratch,
                                    ) {
                                        Some(value) => value,
                                        None => {
                                            continue;
                                        }
                                    };

                                let (p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss) = prelim;
                                (p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss)
                            };

                            dbg_gapped_calls += 1;
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4058-4091
                            // ```c
                            // if (gap_align->score >= cutoff) {
                            //     status = Blast_HSPInit(gap_align->query_start,
                            //                           gap_align->query_stop,
                            //                           gap_align->subject_start,
                            //                           gap_align->subject_stop, ...);
                            //     status = BlastIntervalTreeAddHSP(new_hsp, tree, query_info,
                            //                                      eQueryAndSubject);
                            // }
                            // ```
                            let trace_prelim_seed = blastn_trace_enabled
                                && blastn_trace::should_trace_seed(
                                    "prelim",
                                    uh.context_idx,
                                    s_idx,
                                    s_id,
                                    seed_qs,
                                    seed_ss.saturating_add(chunk.offset),
                                );
                            if blastn_trace_enabled
                                && (trace_prelim_seed
                                    || blastn_trace::should_trace_range(
                                        "prelim",
                                        uh.context_idx,
                                        s_idx,
                                        s_id,
                                        prelim_qs,
                                        prelim_qe,
                                        prelim_ss.saturating_add(chunk.offset),
                                        prelim_se.saturating_add(chunk.offset),
                                        ctx.len,
                                        ctx.frame,
                                    ))
                            {
                                blastn_trace::log(
                            "prelim",
                            format!(
                                "subject={}({}) context={} seed=({}, {}) prelim=q{}..{} s{}..{} raw_score={} cutoff={} x_drop_score_only={} accepted={}",
                                s_id,
                                s_idx,
                                uh.context_idx,
                                seed_qs,
                                seed_ss.saturating_add(chunk.offset),
                                prelim_qs,
                                prelim_qe,
                                prelim_ss.saturating_add(chunk.offset),
                                prelim_se.saturating_add(chunk.offset),
                                prelim_score,
                                cutoff_score,
                                if use_dp { x_drop_gapped.min(uh.score) } else { x_drop_gapped },
                                prelim_score >= cutoff_score
                            ),
                        );
                            }
                            if prelim_score < cutoff_score {
                                continue;
                            }

                            // Add preliminary GAPPED HSP to interval tree for containment checks
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3913
                            // ```c
                            // tmp_hsp.query.offset = q_start;
                            // tmp_hsp.query.end = q_end;
                            // tmp_hsp.query.frame = query_info->contexts[context].frame;
                            // tmp_hsp.subject.offset = s_start;
                            // tmp_hsp.subject.end = s_end;
                            // ```
                            let gapped_tree_hsp = TreeHsp {
                                query_offset: prelim_qs as i32,
                                query_end: prelim_qe as i32,
                                subject_offset: prelim_ss as i32,
                                subject_end: prelim_se as i32,
                                score: prelim_score,
                                query_frame: uh.query_frame,
                                query_length: ctx.len as i32,
                                query_context_offset: uh.query_context_offset,
                                subject_frame_sign: 1,
                            };
                            interval_tree.add_hsp(
                                gapped_tree_hsp,
                                uh.query_context_offset,
                                IndexMethod::QueryAndSubject,
                            );

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4031
                            // ```c
                            // if (init_hsp->ungapped_data) {
                            //    init_hsp->offsets.qs_offsets.q_off =
                            //        init_hsp->ungapped_data->q_start + init_hsp->ungapped_data->length/2;
                            //    init_hsp->offsets.qs_offsets.s_off =
                            //        init_hsp->ungapped_data->s_start + init_hsp->ungapped_data->length/2;
                            // }
                            // ```
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4058-4076
                            // ```c
                            // if (gap_align->score >= cutoff) {
                            //    status = Blast_HSPInit(gap_align->query_start,
                            //              gap_align->query_stop, gap_align->subject_start,
                            //              gap_align->subject_stop,
                            //              init_hsp->offsets.qs_offsets.q_off,
                            //              init_hsp->offsets.qs_offsets.s_off, context,
                            //              query_frame, subject->frame, gap_align->score,
                            //              &(gap_align->edit_script), &new_hsp);
                            // }
                            // ```
                            prelim_hits.push(PrelimHit {
                                context_idx: uh.context_idx,
                                query_idx: uh.query_idx,
                                query_frame: uh.query_frame,
                                query_context_offset: uh.query_context_offset,
                                prelim_qs,
                                prelim_qe,
                                prelim_ss,
                                prelim_se,
                                prelim_score,
                                seed_qs,
                                seed_ss,
                                prelim_evalue: 0.0,
                            });
                            continue;
                        }
                    },
                );
                // No NCBI counterpart: LOSAT_X_STATS work counters (xstats feature); they do not change any
                // value NCBI computes.
                if let Some(x_ahead) = x_ahead.as_ref() {
                    let (x_evaluated, x_taken) = x_ahead.counts();
                    crate::utils::xstats::AHEAD_PRELIM_EVALUATED
                        .fetch_add(x_evaluated, std::sync::atomic::Ordering::Relaxed);
                    crate::utils::xstats::AHEAD_PRELIM_TAKEN
                        .fetch_add(x_taken, std::sync::atomic::Ordering::Relaxed);
                }

                (dbg_gapped_calls, dbg_containment_skipped)
            };
            // Dispatch point of LOSAT_X_QGROUPPAR for the preliminary loop: per query on the
            // pool, each query's purged list concatenated in query order (the serial purge's
            // output, see x_qgroup_par.rs), or the loop over the chunk's list.
            let mut x_purged = false;
            if let Some(order) = x_order.as_ref() {
                let x_starts = x_qgroup_par::starts(x_counts.as_deref().unwrap_or(&[]));
                let mut x_rest: &mut [UngappedHit] = ungapped_hits;
                let mut x_slices = Vec::with_capacity(queries.len());
                for query in 0..queries.len() {
                    let (head, tail) = x_rest.split_at_mut(x_starts[query + 1] - x_starts[query]);
                    x_slices.push(std::sync::Mutex::new(head));
                    x_rest = tail;
                }
                let x_tree_q_max = (query_concat_length + 1) as i32;
                let x_tree_s_max = (s_len + 1) as i32;
                let x_results = x_qgroup_par::run(
                    parallel_pool,
                    x_qgroup_par::Stage::Prelim,
                    order,
                    queries.len(),
                    || {
                        (
                            BlastIntervalTree::new(0, x_tree_q_max, 0, x_tree_s_max),
                            GapAlignScratch::new(),
                            GreedyAlignScratch::new(),
                        )
                    },
                    |(tree, gap, greedy), query| {
                        let mut hits = x_slices[query as usize]
                            .lock()
                            .expect("each query's hits are taken by one task");
                        hits.sort_by(score_compare_ungapped_hits);
                        tree.reset();
                        let mut out = Vec::new();
                        let dbg = x_prelim_list(&hits, tree, gap, greedy, None, &mut out);
                        (purge_prelim_hits_with_common_endpoints(out), dbg)
                    },
                );
                let mut x_hits: Vec<PrelimHit> = Vec::new();
                for (hits, (calls, skipped)) in x_results.into_iter().flatten() {
                    dbg_gapped_calls += calls;
                    dbg_containment_skipped += skipped;
                    x_hits.extend(hits);
                }
                *prelim_hits = x_hits;
                x_purged = true;
            } else {
                let (calls, skipped) = x_prelim_list(
                    ungapped_hits,
                    interval_tree,
                    gap_scratch,
                    greedy_align_scratch,
                    Some(&mut subject_scratch.x_ahead_scratch),
                    prelim_hits,
                );
                dbg_gapped_calls = calls;
                dbg_containment_skipped = skipped;
            }

            // DEBUG: Log gapped extension stats
            let gapped_elapsed = gapped_start_time.elapsed();
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3917-3922
            // ```c
            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                   hit_options->min_diag_separation))
            // {
            //    BlastHSP* new_hsp;
            // ```
            if let Some(timing) = timing_ref {
                timing.gapped_ns.fetch_add(
                    gapped_elapsed.as_nanos() as u64,
                    std::sync::atomic::Ordering::Relaxed,
                );
                timing.gapped_calls.fetch_add(
                    dbg_gapped_calls as u64,
                    std::sync::atomic::Ordering::Relaxed,
                );
            }
            if verbose {
                eprintln!(
                    "[INFO] Gapped extension took {:?}: calls={}, skipped={}, total_ungapped={}",
                    gapped_elapsed, dbg_gapped_calls, dbg_containment_skipped, total_ungapped
                );
            }

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
            // ```c
            // /* Delete if not done in last loop iteration to prevent memory leak. */
            // hsp_list = Blast_HSPListFree(hsp_list);
            //
            // BlastInitHitListReset(init_hitlist);
            // ```
            let mut chunk_prelim_hits: Vec<PrelimHit> = if reuse_prelim_hits {
                std::mem::take(prelim_hits)
            } else {
                prelim_hits.drain(..).collect()
            };
            if !chunk_prelim_hits.is_empty() {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4058-4085
                // ```c
                // if (gap_align->score >= cutoff) {
                //     status = Blast_HSPInit(..., &new_hsp);
                //     status = Blast_HSPListSaveHSP(hsp_list, new_hsp);
                //     status = BlastIntervalTreeAddHSP(new_hsp, tree, query_info,
                //                                      eQueryAndSubject);
                // }
                // ```
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:540-557
                // ```c
                // if (aux_struct->GetGappedScore) {
                //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
                // }
                // Blast_HSPListSortByScore(hsp_list);
                // ```
                // NCBI purges redundant preliminary score-only gapped HSPs
                // before subject chunk offset adjustment and final traceback.
                // With LOSAT_X_QGROUPPAR the queries' lists were purged by their tasks, and the
                // same stable sort runs on the pool.
                if x_purged {
                    x_qgroup_par::sort_stable(
                        parallel_pool,
                        &mut chunk_prelim_hits,
                        score_compare_prelim_hits,
                    );
                } else {
                    chunk_prelim_hits = purge_prelim_hits_with_common_endpoints(chunk_prelim_hits);
                    sort_prelim_hits_by_score(&mut chunk_prelim_hits);
                }
                for hit in chunk_prelim_hits.iter_mut() {
                    hit.prelim_ss = hit.prelim_ss.saturating_add(chunk.offset);
                    hit.prelim_se = hit.prelim_se.saturating_add(chunk.offset);
                    hit.seed_ss = hit.seed_ss.saturating_add(chunk.offset);
                }
            }
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
                dbg_containment_skipped,
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
        };

        let mut combined_prelim_hits: Vec<PrelimHit> = Vec::new();
        if let Some(kept) = prelim_source {
            // The traceback of the HSPs that the preliminary stage kept (`search_subjects`).
            combined_prelim_hits.clone_from(kept);
        } else {
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
                                let hits = collect_prelim_hits_for_chunk(
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
                                    prelim_subject_best_hit(
                                        &mut combined_prelim_hits,
                                        &query_lengths,
                                    );
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
                                            SubjectScratch::new(
                                                queries_ref.len(),
                                                offset_array_size,
                                                window_size,
                                            ),
                                        )
                                    },
                                    |state, chunk| {
                                        let (gap_scratch, subject_scratch) = state;
                                        let hits = collect_prelim_hits_for_chunk(
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
                                    prelim_subject_best_hit(
                                        &mut combined_prelim_hits,
                                        &query_lengths,
                                    );
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
                            crate::utils::threading::report_stage(
                                "blastn",
                                "subject_chunks",
                                1,
                                false,
                            );
                            let soft_ranges = split_state.soft_ranges.as_slice();
                            let hits = collect_prelim_hits_for_chunk(
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
                                prelim_subject_best_hit(&mut combined_prelim_hits, &query_lengths);
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

            if !UNGAPPED && !combined_prelim_hits.is_empty() {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:885-899
                // ```c
                // Blast_HSPListGetEvalues(program_number, query_info,
                //                                  stat_length, hsp_list_out,
                //                                  score_options->gapped_calculation,
                //                                  isRPS, gap_align->sbp, 0, scale_factor);
                // ...
                // status = s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
                // ```
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:642-671
                // ```c
                // cutoff = hit_params->prelim_evalue;
                // ...
                // if (hsp->evalue > cutoff) {
                //     hsp_array[index] = Blast_HSPFree(hsp_array[index]);
                // }
                // ```
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:868,933
                // ```c
                // params->prelim_evalue = options->expect_value;
                // ```
                // This reap happens after preliminary gapped scoring and before the
                // traceback stream is consumed, so low preliminary scores cannot survive
                // solely because traceback later finds a high-scoring alignment.
                combined_prelim_hits.retain_mut(|prelim| {
                    let eff_searchsp = query_eff_searchsp
                        .get(prelim.context_idx as usize)
                        .copied()
                        .unwrap_or(0);
                    let (_, prelim_evalue) = calculate_blastn_context_statistics(
                        prelim.prelim_score,
                        &search_karlin_ref[prelim.context_idx as usize].gapped,
                        eff_searchsp,
                        round_down_evalue_score,
                    );
                    prelim.prelim_evalue = prelim_evalue;
                    hsp_survives_evalue_reap(prelim_evalue, evalue_threshold)
                });
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1482-1541
            // ```c
            //       if (hsp_list && hsp_list->hspcnt > 0) {
            //          if (!gapped_calculation) {
            //             status = Blast_HSPListReevaluateUngapped(program_number, hsp_list, query,
            //                         seq_arg.seq, word_params, hit_params, query_info, sbp,
            //                         score_params, seq_src, seq_arg.seq->gen_code_string);
            //             ...
            //                Blast_HSPListGetEvalues(program_number, query_info, stat_length,
            //                                        hsp_list, gapped_calculation, FALSE, sbp, 0, 1.0);
            //             ...
            //         	   status = s_Blast_HSPListReapByPrelimEvalue(hsp_list, hit_params);
            //             Blast_HSPListReapByQueryCoverage(hsp_list,hit_params->options, query_info, program_number);
            //             Blast_HSPListGetBitScores(hsp_list, gapped_calculation, sbp);
            // ```
            // An ungapped (`-ungapped`) search finishes each subject's HSP list here, in the
            // preliminary stage: the traceback only keeps what this leaves.
            if UNGAPPED && !combined_prelim_hits.is_empty() {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2645-2668
                // ```c
                //       	seq_arg.encoding = (kTranslateSubject ? eBlastEncodingNcbi4na : eBlastEncodingNucleotide);
                //       	BlastSeqSrcReleaseSequence(seq_src, &seq_arg);
                //       	if (( BLAST_SEQSRC_EXCLUDED == BlastSeqSrcGetSequence(seq_src, &seq_arg))) {
                // ```
                // The subject again, unpacked (BLASTNA, without sentinel bytes). A database
                // sequence has none here: the preliminary search's caller reads it again, for
                // the sequences that keep an HSP, and finishes the list itself
                // (`x_ungapped_database_stage`, lazily as NCBI does).
                let mut x_subject_vec: Vec<u8>;
                let s_blastna: Option<&[u8]> =
                    match (subject.blastna, subject.store, subject.letters) {
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
                    // (`percent_identity` and the cutoffs are the values the traceback below
                    // reads already; the query coverage is reaped by the collector's caller.)
                    ungapped_subject_stage(
                        &mut combined_prelim_hits,
                        s_blastna,
                        encoded_queries_blastna.as_slice(),
                        &score_matrix,
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
                        evalue_threshold,
                    );
                }
            }
        }

        // The preliminary stage gives the subject's HSPs to the collector (`search_subjects`).
        if prelim_source.is_none() {
            *prelim_out = combined_prelim_hits;
            return;
        }
        if combined_prelim_hits.is_empty() {
            return;
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:358-373
        // ```c
        // /* Make sure the HSPs in the HSP list are sorted by score, as they should be. */
        // ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
        // tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
        //                               0, subject_length + 1);
        // ```
        let traceback_start = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };

        let hits = &mut subject_scratch.hits;
        hits.clear();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1528, 1567-1573
        // ```c
        //             const EBlastEncoding encoding = Blast_TracebackGetEncoding(program_number);
        // ...
        //                 seq_arg.oid = batch->hsplist_array[0]->oid;
        //                 seq_arg.encoding = encoding;
        // ...
        //                 if (BlastSeqSrcGetSequence(seqsrc, &seq_arg) < 0) {
        // ```
        // The traceback fetches the uncompressed (blastna) subject of each subject with
        // preliminary HSPs; the preliminary stage reads only the packed one. The blastna letters
        // are made here, for this subject only, instead of being kept for every subject of the
        // search. A database sequence comes with its blastna letters, read with its
        // ambiguity codes (`database::traceback_listed`).
        let mut s_seq_blastna_vec: Vec<u8>;
        let s_seq_blastna: &[u8] = match (subject.blastna, subject.store) {
            (Some(blastna), _) => blastna,
            // EXPERIMENT (LOSAT_X_BNSTORE): made from the packed subject.
            (None, Some(store)) => {
                s_seq_blastna_vec = Vec::new();
                x_blastna_into(store, &mut s_seq_blastna_vec);
                s_seq_blastna_vec.as_slice()
            }
            (None, None) => {
                s_seq_blastna_vec = encode_iupac_to_blastna(
                    subject
                        .letters
                        .expect("a subject reaches the traceback with its letters"),
                );
                s_seq_blastna_vec.as_slice()
            }
        };

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:1526,1690-1706
        // ```c
        //     const Boolean perform_traceback = score_params->options->gapped_calculation;
        //     ...
        //             if (perform_traceback) { ... }
        //             /* Update bit scores */
        //             Blast_HSPListGetBitScores(hsp_list, kGapped, sbp);
        //             ...
        //             status = Blast_HSPResultsInsertHSPList(results, hsp_list, hit_options->hitlist_size);
        // ```
        // The traceback of an ungapped (`-ungapped`) search aligns nothing and purges nothing:
        // the HSPs the preliminary stage kept (re-evaluated, with their e-values) become the
        // subject's hits in the kept order, with the identities counted again on the subject's
        // BLASTNA letters (the ranges the traceback decoded include the HSPs').
        if UNGAPPED {
            let mut final_hits: Vec<BlastnHsp> = Vec::with_capacity(combined_prelim_hits.len());
            for prelim in combined_prelim_hits.iter() {
                let context_idx = prelim.context_idx as usize;
                let length = prelim.prelim_qe - prelim.prelim_qs;
                let matches = ungapped::count_identities_ungapped(
                    encoded_queries_blastna[context_idx],
                    s_seq_blastna,
                    prelim.prelim_qs,
                    prelim.prelim_ss,
                    length,
                );
                let eff_searchsp = query_eff_searchsp[context_idx];
                let (bit_score, eval) = calculate_blastn_context_statistics(
                    prelim.prelim_score,
                    &search_karlin_ref[context_idx].gapped,
                    eff_searchsp,
                    round_down_evalue_score,
                );
                let identity = if length > 0 {
                    ((matches as f64 / length as f64) * 100.0).min(100.0)
                } else {
                    0.0
                };
                let query_length = query_lengths[prelim.query_idx as usize];
                let (hit_q_start, hit_q_end, hit_s_start, hit_s_end) = adjust_blastn_offsets(
                    prelim.prelim_qs,
                    prelim.prelim_qe,
                    prelim.prelim_ss,
                    prelim.prelim_se,
                    query_length,
                    prelim.query_frame,
                );
                final_hits.push(BlastnHsp {
                    identity,
                    length,
                    mismatch: length - matches,
                    gapopen: 0,
                    q_start: hit_q_start,
                    q_end: hit_q_end,
                    s_start: hit_s_start,
                    s_end: hit_s_end,
                    e_value: eval,
                    bit_score,
                    num_ident: matches,
                    query_frame: prelim.query_frame,
                    query_length,
                    q_idx: prelim.query_idx,
                    s_idx: s_idx as u32,
                    raw_score: prelim.prelim_score,
                    internal_q_offset_0: prelim.prelim_qs,
                    internal_q_end_0: prelim.prelim_qe,
                    internal_s_offset_0: prelim.prelim_ss,
                    internal_s_end_0: prelim.prelim_se,
                    internal_query_context_offset: prelim.query_context_offset,
                    gap_info: None,
                    num_positives: matches,
                });
            }
            if !final_hits.is_empty() {
                *subject_hits = Some(final_hits);
            }
            if let Some(bar) = progress_bar.as_ref() {
                bar.inc(1);
            }
            return;
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:371-373
        // ```c
        // tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
        //                               0, subject_length + 1);
        // ```
        // The tree is the scratch's (`x_tree`), given this subject's ranges.
        let interval_tree = &mut subject_scratch.x_tree;
        interval_tree.reset_with_bounds(
            0,
            (query_concat_length + 1) as i32,
            0,
            (s_len_full + 1) as i32,
        );
        let mut traceback_edit_script_lengths: Vec<u32> = Vec::new();
        let mut traceback_alignment_lengths: Vec<u32> = Vec::new();

        // The traceback loop over one list of kept HSPs: the subject's on the owner (with the
        // look-ahead helpers, whose scratch `x_helper_scratch` holds, or the speculative
        // batches), or one query's (LOSAT_X_QGROUPPAR, neither). The HSPs go to `hits`.
        let x_trace_list = |prelim_hits: &[PrelimHit],
                            interval_tree: &mut BlastIntervalTree,
                            gap_scratch: &mut GapAlignScratch,
                            greedy_align_scratch: &mut GreedyAlignScratch,
                            x_helper_scratch: Option<
            &mut Vec<(GapAlignScratch, GreedyAlignScratch)>,
        >,
                            x_spec_scratch: Option<
            &mut Vec<(GapAlignScratch, GreedyAlignScratch)>,
        >,
                            hits: &mut Vec<BlastnHsp>,
                            traceback_edit_script_lengths: &mut Vec<u32>,
                            traceback_alignment_lengths: &mut Vec<u32>| {
            let speculative_traceback = speculative_traceback && x_helper_scratch.is_some();
            #[cfg(not(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            )))]
            let _ = &x_spec_scratch;
            if x_greedy_spec {
                crate::utils::xreach::x_reach!("LOSAT_X_GREEDYSPEC", speculative_traceback);
            }
            if !prelim_hits.is_empty() {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:358-365
                // ```c
                // /* Make sure the HSPs in the HSP list are sorted by score, as they
                //    should be. */
                // #ifdef _DEBUG
                // {
                //     Blast_HSPListSortByScore(hsp_list);
                // }
                // #endif
                // ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
                // ```
                // The oracle is a release build: the HSPs are traced in the order of the
                // stored list (`kept`, each query's list for this subject in turn). That
                // is score order for a list kept as it came (sorted per subject chunk,
                // blast_engine.c:555, and merged by score), and e-value order for a list
                // that `Blast_HitListUpdate` sorted when the hit list became a heap
                // (blast_hits.c:3272-3284).
                if let Some(timing) = timing_ref {
                    BlastnTiming::record_count(
                        &timing.traceback_prelim_hsps,
                        prelim_hits.len() as u64,
                    );
                }

                interval_tree.reset();
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:66-70
                // ```c
                // if (tree->num_used == tree->num_alloc) {
                //     tree->num_alloc = 2 * tree->num_alloc;
                //     tree->nodes = (SIntervalNode *)realloc(tree->nodes, tree->num_alloc *
                //                                                  sizeof(SIntervalNode));
                // }
                // ```
                // Reserve from the preliminary HSP count before the hot add/contains loop.
                // This preserves NCBI node order while avoiding repeated Vec growth.
                interval_tree.reserve_nodes_for_hsps(prelim_hits.len());

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405,436-472,583-612
                // ```c
                // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, hit_options->min_diag_separation)) {
                //     BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                //     AdjustSubjectRange(&s_start, &adjusted_s_length, q_start, query_length, &start_shift);
                //     /* traceback, identity test, then BlastIntervalTreeAddHSP */
                // }
                // ```
                // N02 scheduling only: speculative values have no externally visible effects;
                // the original ordered contains/materialize/add sequence remains authoritative.
                let prepare_traceback = |prelim: &PrelimHit| {
                    let q_seq_blastna = encoded_queries_blastna[prelim.context_idx as usize];
                    let mut trace_q_start = prelim.seed_qs;
                    let mut trace_s_start = prelim.seed_ss;

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:436-460
                    // ```c
                    // if (!kIsOutOfFrame && hsp->query.gapped_start == 0 &&
                    //                       hsp->subject.gapped_start == 0) {
                    //    Boolean retval =
                    //       BlastGetOffsetsForGappedAlignment(query, subject, sbp,
                    //           hsp, &q_start, &s_start);
                    //    if (!retval) { ... }
                    //    hsp->query.gapped_start = q_start;
                    //    hsp->subject.gapped_start = s_start;
                    // } else {
                    //    ...
                    //    BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                    //    q_start = hsp->query.gapped_start;
                    //    s_start = hsp->subject.gapped_start;
                    // }
                    // ```
                    if trace_q_start == 0 && trace_s_start == 0 {
                        let (q_start, s_start) = match blast_get_offsets_for_gapped_alignment(
                            q_seq_blastna,
                            s_seq_blastna,
                            prelim.prelim_qs,
                            prelim.prelim_qe,
                            prelim.prelim_ss,
                            prelim.prelim_se,
                            &score_matrix,
                        ) {
                            Some(value) => value,
                            None => return None,
                        };
                        trace_q_start = q_start;
                        trace_s_start = s_start;
                    } else {
                        let (q_start, s_start) = blast_get_start_for_gapped_alignment_nucl(
                            q_seq_blastna,
                            s_seq_blastna,
                            prelim.prelim_qs,
                            prelim.prelim_qe,
                            prelim.prelim_ss,
                            prelim.prelim_se,
                            trace_q_start,
                            trace_s_start,
                        );
                        trace_q_start = q_start;
                        trace_s_start = s_start;
                    }

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:466-472
                    // ```c
                    // AdjustSubjectRange(&s_start, &adjusted_s_length, q_start,
                    //                    query_length, &start_shift);
                    // adjusted_subject = subject + start_shift;
                    // hsp->subject.gapped_start = s_start;
                    // ```
                    let mut adjusted_s_len = s_seq_blastna.len();
                    let mut s_start_i32 = trace_s_start as i32;
                    let mut s_len_i32 = adjusted_s_len as i32;
                    let start_shift_i32 = adjust_subject_range(
                        &mut s_start_i32,
                        &mut s_len_i32,
                        trace_q_start as i32,
                        q_seq_blastna.len() as i32,
                    );
                    let start_shift = start_shift_i32 as usize;
                    adjusted_s_len = s_len_i32 as usize;
                    let trace_s_start_adj = s_start_i32 as usize;

                    Some((
                        trace_q_start,
                        trace_s_start,
                        start_shift,
                        adjusted_s_len,
                        trace_s_start_adj,
                    ))
                };
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405,436-472,583-612
                // ```c
                // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, hit_options->min_diag_separation)) {
                //     BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                //     AdjustSubjectRange(&s_start, &adjusted_s_length, q_start, query_length, &start_shift);
                //     /* traceback, identity test, then BlastIntervalTreeAddHSP */
                // }
                // ```
                // N02 scheduling only: speculative values have no externally visible effects;
                // the original ordered contains/materialize/add sequence remains authoritative.
                #[cfg(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                ))]
                let mut x_no_spec = Vec::new();
                #[cfg(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                ))]
                let x_spec_scratch = x_spec_scratch.unwrap_or(&mut x_no_spec);
                #[cfg(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                ))]
                let mut speculative_scratch: Vec<_> = if speculative_traceback {
                    // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:509-513
                    // BLAST_GappedAlignmentWithTraceback(...);
                    // Each independent scratch slot retains its intermediate results.
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-507
                    // ```c
                    //         } else if (kGreedyTraceback) {
                    //             BLAST_GreedyGappedAlignment(query, adjusted_subject,
                    //                 query_length, adjusted_s_length, gap_align,
                    // ```
                    // Each slot of the batch also has the scratch of the greedy aligner, so that
                    // the greedy traceback call can run on a slot instead of on the ordered loop.
                    // Scratch only, made once per subject scratch (`x_spec_scratch`).
                    while x_spec_scratch.len() < num_threads {
                        x_spec_scratch.push((GapAlignScratch::new(), GreedyAlignScratch::new()));
                    }
                    x_spec_scratch
                        .iter_mut()
                        .take(num_threads)
                        .map(|(gap, greedy)| (gap, greedy, Vec::new()))
                        .collect()
                } else {
                    Vec::new()
                };
                #[cfg(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                ))]
                let mut speculative_results = Vec::new();
                // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:403-405,598-601
                // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, ...)) {
                //     BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
                // }
                // Reuse only batch storage; live ordered containment remains below.
                #[cfg(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                ))]
                let mut speculative_jobs = Vec::new();
                // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:403-405,509-513,598-601
                // ```c
                // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, ...)) {
                //     BLAST_GappedAlignmentWithTraceback(...);
                //     ...
                //     BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
                // }
                // ```
                // Scheduling only; ordered containment and insertion stay unchanged.
                // Serial Wasm compiles out the speculative path and retains batch 8.
                const X_DEFAULT_SPECULATIVE_TRACEBACK_BATCH_SIZE: usize =
                    if cfg!(all(target_arch = "wasm32", not(feature = "wasm-threads"))) {
                        8
                    } else {
                        16
                    };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375
                // ```c
                //    for (index=0; index < num_initial_hsps; index++) {
                // ```
                // The number of HSPs whose tracebacks are computed ahead of this loop in one batch. Scheduling
                // only: the loop still visits the HSPs in this order and decides each one itself.
                // EXPERIMENT (LOSAT_X_SPECBATCH=n): how many HSPs are traced ahead of
                // the ordered loop at a time. Scheduling only, like the batch itself.
                #[allow(non_snake_case)]
                let SPECULATIVE_TRACEBACK_BATCH_SIZE: usize = {
                    use std::sync::OnceLock;
                    static SIZE: OnceLock<usize> = OnceLock::new();
                    *SIZE.get_or_init(|| {
                        std::env::var("LOSAT_X_SPECBATCH")
                            .ok()
                            .and_then(|v| v.parse::<usize>().ok())
                            .filter(|&n| (1..=4096).contains(&n))
                            .unwrap_or(X_DEFAULT_SPECULATIVE_TRACEBACK_BATCH_SIZE)
                    })
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-440,457,466-469,503-509
                // ```c
                //          if (!kIsOutOfFrame && hsp->query.gapped_start == 0 &&
                //                                hsp->subject.gapped_start == 0) {
                //             Boolean retval =
                //                BlastGetOffsetsForGappedAlignment(query, subject, sbp,
                // ...
                //                BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                // ...
                //             AdjustSubjectRange(&s_start, &adjusted_s_length, q_start,
                //                                query_length, &start_shift);
                // ...
                //         } else if (kGreedyTraceback) {
                //             BLAST_GreedyGappedAlignment(query, adjusted_subject,
                // ...
                //         } else {
                //           BLAST_GappedAlignmentWithTraceback(program_number, query,
                // ```
                // This closure is the part of the loop body that depends on the HSP alone: the start offsets and
                // the subject range (`prepare_traceback`) and the traceback call, dynamic-programming or greedy.
                // Index in, value out. The tree test and the tree insertion stay in the ordered loop.
                // No NCBI counterpart: the traceback of the loop below as a
                // function of the HSP's index (the job of the speculative batch).
                let x_trace_compute = |x_index: usize, x_scratch: &mut XAheadScratch<'_>| {
                    let p = &prelim_hits[x_index];
                    let prepared = prepare_traceback(p)?;
                    let (qs, _, shift, slen, ss) = prepared;
                    if x_fence::x_fenced(x_fenced) {
                        return None;
                    }
                    let x_mark = x_fence::x_fence_mark(x_fenced);
                    if !use_dp {
                        let result = greedy_gapped_alignment_with_traceback(
                            &encoded_queries_blastna[p.context_idx as usize],
                            &s_seq_blastna[shift..shift + slen],
                            slen,
                            qs,
                            ss,
                            reward,
                            penalty,
                            gap_open,
                            gap_extend,
                            x_drop_final,
                            &mut *x_scratch.greedy,
                        );
                        x_fence::x_fence_note(x_fenced, x_mark);
                        return Some((prepared, XSpecTraceback::Greedy(result)));
                    }
                    let result = extend_gapped_heuristic_with_traceback_with_scratch(
                        &encoded_queries_blastna[p.context_idx as usize],
                        &s_seq_blastna[shift..shift + slen],
                        qs,
                        ss,
                        1,
                        reward,
                        penalty,
                        &score_matrix,
                        gap_open,
                        gap_extend,
                        x_drop_final,
                        &mut *x_scratch.gap,
                    );
                    x_fence::x_fence_note(x_fenced, x_mark);
                    Some((prepared, XSpecTraceback::Dp(result)))
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375
                // ```c
                //    for (index=0; index < num_initial_hsps; index++) {
                // ```
                // Look-ahead for the traceback loop. Without helpers `x_ahead` is `None` and
                // the loop runs as ported.
                let x_ahead = match x_ahead_window {
                    Some(x_window)
                        if x_helper_scratch.is_some()
                            && !blastn_trace_enabled
                            && prelim_hits.len() >= 2 * x_window =>
                    {
                        Some(crate::utils::xahead::Ahead::new(x_window))
                    }
                    _ => None,
                };
                if x_ahead_window.is_some() {
                    crate::utils::xreach::x_reach!("LOSAT_X_AHEAD", x_ahead.is_some());
                }
                let mut x_no_helpers = Vec::new();
                let x_helper_scratch = x_helper_scratch.unwrap_or(&mut x_no_helpers);
                if x_ahead.is_some() {
                    while x_helper_scratch.len() + 1 < num_threads {
                        x_helper_scratch.push((GapAlignScratch::new(), GreedyAlignScratch::new()));
                    }
                }
                let x_helpers = x_helper_scratch
                    .iter_mut()
                    .map(|(gap, greedy)| XAheadScratch { gap, greedy });
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,600-601
                // ```c
                //    for (index=0; index < num_initial_hsps; index++) {
                // ...
                //               status = BlastIntervalTreeAddHSP(hsp, tree, query_info,
                //                                          eQueryAndSubject);
                // ```
                // `with_helpers` runs the closure on this thread. The closure is the port of this NCBI loop, in NCBI
                // order; the other threads of the pool only evaluate tracebacks ahead.
                crate::utils::xahead::with_helpers(
                    x_ahead.as_ref(),
                    x_helpers,
                    &x_trace_compute,
                    || {
                        for prelim_index in 0..prelim_hits.len() {
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,404-405
                            // ```c
                            //    for (index=0; index < num_initial_hsps; index++) {
                            // ...
                            //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                            //                                hit_options->min_diag_separation)) {
                            // ```
                            // The owner offers later HSPs to the helpers, leaving out the ones the tree already contains. The
                            // tree is only read here.
                            // No NCBI counterpart: offer the next HSPs to the helpers,
                            // leaving out those the tree already contains.
                            if let Some(x_ahead) = x_ahead.as_ref() {
                                x_ahead.announce(prelim_index, prelim_hits.len(), |x_index| {
                                    let p = &prelim_hits[x_index];
                                    let h = TreeHsp {
                                        query_offset: p.prelim_qs as i32,
                                        query_end: p.prelim_qe as i32,
                                        subject_offset: p.prelim_ss as i32,
                                        subject_end: p.prelim_se as i32,
                                        score: p.prelim_score,
                                        query_frame: p.query_frame,
                                        query_length: query_contexts[p.context_idx as usize].len
                                            as i32,
                                        query_context_offset: p.query_context_offset,
                                        subject_frame_sign: 1,
                                    };
                                    interval_tree
                                        .containing_hsp(
                                            &h,
                                            p.query_context_offset,
                                            min_diag_separation,
                                        )
                                        .is_none()
                                });
                            }
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405,436-472,583-612
                            // ```c
                            // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, hit_options->min_diag_separation)) {
                            //     BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                            //     AdjustSubjectRange(&s_start, &adjusted_s_length, q_start, query_length, &start_shift);
                            //     /* traceback, identity test, then BlastIntervalTreeAddHSP */
                            // }
                            // ```
                            // N02 scheduling only: speculative values have no externally visible effects;
                            // the original ordered contains/materialize/add sequence remains authoritative.
                            #[cfg(all(
                                feature = "parallel",
                                any(not(target_arch = "wasm32"), feature = "wasm-threads")
                            ))]
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,509-512
                            // ```c
                            //    for (index=0; index < num_initial_hsps; index++) {
                            // ...
                            //           BLAST_GappedAlignmentWithTraceback(program_number, query,
                            //                 adjusted_subject, gap_align, score_params, q_start, s_start,
                            //                 query_length, adjusted_s_length,
                            //                 fence_hit);
                            // ```
                            // Look-ahead: with helpers the batch below is not used; the tracebacks come
                            // from `x_ahead.take`. Both ways only compute the traceback call of this loop ahead of its turn.
                            if speculative_traceback
                                && x_ahead.is_none()
                                && prelim_index % SPECULATIVE_TRACEBACK_BATCH_SIZE == 0
                            {
                                let batch_end = (prelim_index + SPECULATIVE_TRACEBACK_BATCH_SIZE)
                                    .min(prelim_hits.len());
                                speculative_jobs.clear();
                                speculative_jobs.extend((prelim_index..batch_end).filter(
                                    |&index| {
                                        let p = &prelim_hits[index];
                                        let h = TreeHsp {
                                            query_offset: p.prelim_qs as i32,
                                            query_end: p.prelim_qe as i32,
                                            subject_offset: p.prelim_ss as i32,
                                            subject_end: p.prelim_se as i32,
                                            score: p.prelim_score,
                                            query_frame: p.query_frame,
                                            query_length: query_contexts[p.context_idx as usize].len
                                                as i32,
                                            query_context_offset: p.query_context_offset,
                                            subject_frame_sign: 1,
                                        };
                                        interval_tree
                                            .containing_hsp(
                                                &h,
                                                p.query_context_offset,
                                                min_diag_separation,
                                            )
                                            .is_none()
                                    },
                                ));
                                let jobs = &speculative_jobs;
                                // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:509-513
                                // BLAST_GappedAlignmentWithTraceback(...);
                                crate::utils::threading::report_stage(
                                    "blastn",
                                    "dp_traceback",
                                    jobs.len(),
                                    jobs.len() > 1,
                                );
                                speculative_results.clear();
                                speculative_results.resize_with(batch_end - prelim_index, || None);
                                let pool = parallel_pool;
                                pool.install(|| {
                                    speculative_scratch.par_iter_mut().enumerate().for_each(
                                        |(slot, (scratch, x_greedy_scratch, results))| {
                                            results.clear();
                                            results.extend(
                                                jobs.iter().skip(slot).step_by(num_threads).filter_map(
                                                    |&index| {
                                                        let p = &prelim_hits[index];
                                                        // NCBI reference: core/blast_traceback.c:436-472,509-513
                                                        // BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                                                        // AdjustSubjectRange(&s_start, &adjusted_s_length,
                                                        //                    q_start, query_length, &start_shift);
                                                        // BLAST_GappedAlignmentWithTraceback(...);
                                                        // Keep the exact preparation with its speculative DP.
                                                        let prepared = prepare_traceback(p)?;
                                                        let (qs, _, shift, slen, ss) = prepared;
                                                        if x_fence::x_fenced(x_fenced) {
                                                            return None;
                                                        }
                                                        let x_mark = x_fence::x_fence_mark(x_fenced);
                                                        if !use_dp {
                                                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-507
                                                            // ```c
                                                            //         } else if (kGreedyTraceback) {
                                                            //             BLAST_GreedyGappedAlignment(query, adjusted_subject,
                                                            //                 query_length, adjusted_s_length, gap_align,
                                                            //                 score_params, q_start, s_start, FALSE, TRUE,
                                                            //                 fence_hit);
                                                            // ```
                                                            // The greedy traceback call of the ordered loop, made ahead on this slot's
                                                            // scratch. The call depends on its arguments only, so the value is the one the ordered loop would
                                                            // compute.
                                                            // No NCBI counterpart: the call of
                                                            // the ordered loop below, on this slot's scratch.
                                                            let result =
                                                                greedy_gapped_alignment_with_traceback(
                                                                    &encoded_queries_blastna
                                                                        [p.context_idx as usize],
                                                                    &s_seq_blastna[shift..shift + slen],
                                                                    slen,
                                                                    qs,
                                                                    ss,
                                                                    reward,
                                                                    penalty,
                                                                    gap_open,
                                                                    gap_extend,
                                                                    x_drop_final,
                                                                    x_greedy_scratch,
                                                                );
                                                            x_fence::x_fence_note(x_fenced, x_mark);
                                                            return Some((
                                                                index,
                                                                (
                                                                    prepared,
                                                                    XSpecTraceback::Greedy(result),
                                                                ),
                                                            ));
                                                        }
                                                        let result =
                                                    extend_gapped_heuristic_with_traceback_with_scratch(
                                                        &encoded_queries_blastna
                                                            [p.context_idx as usize],
                                                        &s_seq_blastna[shift..shift + slen],
                                                        qs,
                                                        ss,
                                                        1,
                                                        reward,
                                                        penalty,
                                                        &score_matrix,
                                                        gap_open,
                                                        gap_extend,
                                                        x_drop_final,
                                                        scratch,
                                                    );
                                                        x_fence::x_fence_note(x_fenced, x_mark);
                                                        Some((
                                                            index,
                                                            (prepared, XSpecTraceback::Dp(result)),
                                                        ))
                                                    },
                                                ),
                                            );
                                        },
                                    )
                                });
                                // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:583-612
                                // Blast_HSPUpdateWithTraceback(gap_align, hsp);
                                // BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
                                // One owner restores each result to its original batch index.
                                for (_, _, results) in &mut speculative_scratch {
                                    for (index, result) in results.drain(..) {
                                        speculative_results[index - prelim_index] = Some(result);
                                    }
                                }
                            }

                            let prelim = &prelim_hits[prelim_index];
                            let ctx = &query_contexts[prelim.context_idx as usize];
                            let q_seq_blastna =
                                encoded_queries_blastna[prelim.context_idx as usize];
                            let q_seq_nomask_blastna =
                                encoded_queries_blastna[prelim.context_idx as usize];

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405
                            // ```c
                            // if (program_number == eBlastTypeRpsBlast ||
                            //     !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                            //                          hit_options->min_diag_separation)) {
                            // ```
                            let tree_precheck_start = if timing_enabled {
                                Some(std::time::Instant::now())
                            } else {
                                None
                            };
                            let subject_frame_sign = 1i32;
                            let prelim_tree_hsp = TreeHsp {
                                query_offset: prelim.prelim_qs as i32,
                                query_end: prelim.prelim_qe as i32,
                                subject_offset: prelim.prelim_ss as i32,
                                subject_end: prelim.prelim_se as i32,
                                score: prelim.prelim_score,
                                query_frame: prelim.query_frame,
                                query_length: ctx.len as i32,
                                query_context_offset: prelim.query_context_offset,
                                subject_frame_sign,
                            };
                            let traceback_containing_hsp = interval_tree.containing_hsp(
                                &prelim_tree_hsp,
                                prelim.query_context_offset,
                                min_diag_separation,
                            );
                            let prelim_traceback_contained = traceback_containing_hsp.is_some();
                            let trace_traceback_seed = blastn_trace_enabled
                                && blastn_trace::should_trace_seed(
                                    "traceback",
                                    prelim.context_idx,
                                    s_idx,
                                    s_id,
                                    prelim.seed_qs,
                                    prelim.seed_ss,
                                );
                            if let (Some(timing), Some(tree_precheck_start)) =
                                (timing_ref, tree_precheck_start)
                            {
                                BlastnTiming::record_duration(
                                    &timing.traceback_tree_precheck_ns,
                                    tree_precheck_start,
                                );
                            }
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405
                            // ```c
                            // if (program_number == eBlastTypeRpsBlast ||
                            //     !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                            //                              hit_options->min_diag_separation)) {
                            // ```
                            if blastn_trace_enabled
                                && (trace_traceback_seed
                                    || blastn_trace::should_trace_range(
                                        "traceback",
                                        prelim.context_idx,
                                        s_idx,
                                        s_id,
                                        prelim.prelim_qs,
                                        prelim.prelim_qe,
                                        prelim.prelim_ss,
                                        prelim.prelim_se,
                                        ctx.len,
                                        ctx.frame,
                                    ))
                            {
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4071-4076
                                // ```c
                                // status = Blast_HSPInit(gap_align->query_start,
                                //               gap_align->query_stop, gap_align->subject_start,
                                //               gap_align->subject_stop,
                                //               init_hsp->offsets.qs_offsets.q_off,
                                //               init_hsp->offsets.qs_offsets.s_off, context,
                                //               query_frame, subject->frame, gap_align->score,
                                //               &(gap_align->edit_script), &new_hsp);
                                // ```
                                // The init_hsp q_off/s_off values become the traceback
                                // gapped_start seed that BlastGetStartForGappedAlignmentNucl
                                // may keep or move.
                                blastn_trace::log(
                            "traceback",
                            format!(
                                "subject={}({}) context={} prelim=q{}..{} s{}..{} seed=({}, {}) raw_score={} tree_contains={} containing_hsp={:?}",
                                s_id,
                                s_idx,
                                prelim.context_idx,
                                prelim.prelim_qs,
                                prelim.prelim_qe,
                                prelim.prelim_ss,
                                prelim.prelim_se,
                                prelim.seed_qs,
                                prelim.seed_ss,
                                prelim.prelim_score,
                                prelim_traceback_contained,
                                traceback_containing_hsp
                            ),
                        );
                            }
                            if prelim_traceback_contained {
                                if let Some(timing) = timing_ref {
                                    BlastnTiming::record_count(
                                        &timing.traceback_tree_precheck_skipped_hsps,
                                        1,
                                    );
                                }
                                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:404-405,607-610
                                // ```c
                                //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                                //                                hit_options->min_diag_separation)) {
                                // ...
                                //       } else {
                                //          /* Contained within another HSP, delete. */
                                // ```
                                // When the tree contains the HSP, NCBI deletes it and makes no traceback call. The helpers are
                                // told that nobody needs its value.
                                if let Some(x_ahead) = x_ahead.as_ref() {
                                    x_ahead.skip(prelim_index);
                                }
                                continue;
                            }

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4017
                            // ```c
                            // if (init_hsp->ungapped_data) {
                            //     init_hsp->offsets.qs_offsets.q_off =
                            //         init_hsp->ungapped_data->q_start + init_hsp->ungapped_data->length/2;
                            //     init_hsp->offsets.qs_offsets.s_off =
                            //         init_hsp->ungapped_data->s_start + init_hsp->ungapped_data->length/2;
                            // }
                            // ```
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:451-460
                            // ```c
                            // BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                            // q_start = hsp->query.gapped_start;
                            // s_start = hsp->subject.gapped_start;
                            // ```
                            // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:403-405,436-472
                            // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, ...)) {
                            //     BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                            //     AdjustSubjectRange(&s_start, &adjusted_s_length,
                            //                        q_start, query_length, &start_shift);
                            // }
                            // Consume only after the live ordered containment check. An HSP
                            // skipped at batch start can become eligible after endpoint
                            // replacement, so missing results still prepare and run here.
                            #[cfg(all(
                                feature = "parallel",
                                any(not(target_arch = "wasm32"), feature = "wasm-threads")
                            ))]
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:403-405,509-512
                            // ```c
                            //       if (program_number == eBlastTypeRpsBlast ||
                            //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                            //                                hit_options->min_diag_separation)) {
                            // ...
                            //           BLAST_GappedAlignmentWithTraceback(program_number, query,
                            //                 adjusted_subject, gap_align, score_params, q_start, s_start,
                            //                 query_length, adjusted_s_length,
                            //                 fence_hit);
                            // ```
                            // Look-ahead. This is reached after the tree test said "not contained". The
                            // value of the traceback call comes from a helper, or the call is made below as before.
                            let precomputed = if let Some(x_ahead) = x_ahead.as_ref() {
                                // No NCBI counterpart: a helper's value, if one has
                                // started on this HSP; otherwise it is evaluated below.
                                let mut x_scratch = XAheadScratch {
                                    gap: &mut *gap_scratch,
                                    greedy: &mut *greedy_align_scratch,
                                };
                                let x_taken =
                                    x_ahead.take(prelim_index, &mut x_scratch, &x_trace_compute);
                                x_taken.flatten()
                            } else if speculative_traceback {
                                speculative_results[prelim_index % SPECULATIVE_TRACEBACK_BATCH_SIZE]
                                    .take()
                            } else {
                                None
                            };
                            #[cfg(not(all(
                                feature = "parallel",
                                any(not(target_arch = "wasm32"), feature = "wasm-threads")
                            )))]
                            let precomputed: Option<(
                                (usize, usize, usize, usize, usize),
                                _,
                            )> = None;
                            let start_offsets_start = timing_enabled.then(std::time::Instant::now);
                            let prepared = precomputed
                                .as_ref()
                                .map(|(prepared, _)| *prepared)
                                .or_else(|| prepare_traceback(prelim));
                            if let (Some(timing), Some(start)) = (timing_ref, start_offsets_start) {
                                BlastnTiming::record_duration(
                                    &timing.traceback_start_offsets_ns,
                                    start,
                                );
                            }
                            let Some((
                                trace_q_start,
                                trace_s_start,
                                start_shift,
                                adjusted_s_len,
                                trace_s_start_adj,
                            )) = prepared
                            else {
                                continue;
                            };
                            let adjusted_subject =
                                &s_seq_blastna[start_shift..start_shift + adjusted_s_len];
                            let x_drop_trace = x_drop_final;
                            if let Some(timing) = timing_ref {
                                BlastnTiming::record_count(
                                    &timing.traceback_full_traceback_hsps,
                                    1,
                                );
                            }
                            let alignment_start = if timing_enabled {
                                Some(std::time::Instant::now())
                            } else {
                                None
                            };
                            let (
                                final_qs,
                                final_qe,
                                mut final_ss,
                                mut final_se,
                                score,
                                matches,
                                mismatches,
                                gaps,
                                gap_letters,
                                aln_len,
                                edit_ops,
                            ) = if use_dp {
                                let (
                                    final_qs,
                                    final_qe,
                                    final_ss,
                                    final_se,
                                    score,
                                    matches,
                                    mismatches,
                                    gaps,
                                    gap_letters,
                                    edit_ops,
                                ) = {
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405,436-472,583-612
                                    // ```c
                                    // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, hit_options->min_diag_separation)) {
                                    //     BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                                    //     AdjustSubjectRange(&s_start, &adjusted_s_length, q_start, query_length, &start_shift);
                                    //     /* traceback, identity test, then BlastIntervalTreeAddHSP */
                                    // }
                                    // ```
                                    // N02 scheduling only: speculative values have no externally visible effects;
                                    // the original ordered contains/materialize/add sequence remains authoritative.
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:273-286,558-585
                                    // ```c
                                    // if (in_hsp->score > tree_hsp->score) return in_hsp;
                                    // /* equal scores: pick the shorter HSP */
                                    // if (index_method == eQueryAndSubject) { /* check common endpoints */ }
                                    // ```
                                    // Endpoint replacement can invalidate batch-start containment. A newly eligible
                                    // HSP therefore runs the unchanged DP here, at its original sequential position.
                                    if x_fence::x_fenced(x_fenced) {
                                        continue;
                                    }
                                    let x_mark = x_fence::x_fence_mark(x_fenced);
                                    let x_result = match precomputed {
                                        Some((_, XSpecTraceback::Dp(result))) => result,
                                        _ => extend_gapped_heuristic_with_traceback_with_scratch(
                                            q_seq_blastna,
                                            adjusted_subject,
                                            trace_q_start,
                                            trace_s_start_adj,
                                            1,
                                            reward,
                                            penalty,
                                            &score_matrix,
                                            gap_open,
                                            gap_extend,
                                            x_drop_trace,
                                            gap_scratch,
                                        ),
                                    };
                                    // EXPERIMENT (LOSAT_X_BNTBRANGE): an alignment that met a
                                    // fence (here or ahead) makes no HSP; the sequence is traced
                                    // again whole (NCBI breaks its loop, blast_traceback.c:516-520).
                                    x_fence::x_fence_note(x_fenced, x_mark);
                                    if x_fence::x_fenced(x_fenced) {
                                        continue;
                                    }
                                    x_result
                                };
                                let aln_len = matches + mismatches + gap_letters;
                                (
                                    final_qs,
                                    final_qe,
                                    final_ss,
                                    final_se,
                                    score,
                                    matches,
                                    mismatches,
                                    gaps,
                                    gap_letters,
                                    aln_len,
                                    edit_ops,
                                )
                            } else {
                                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-507
                                // ```c
                                //         } else if (kGreedyTraceback) {
                                //             BLAST_GreedyGappedAlignment(query, adjusted_subject,
                                //                 query_length, adjusted_s_length, gap_align,
                                //                 score_params, q_start, s_start, FALSE, TRUE,
                                //                 fence_hit);
                                // ```
                                // The `_` arm is the ported call. The other arm uses the
                                // value of the same call made ahead on another scratch.
                                // No NCBI counterpart: the same call, already made
                                // for this HSP on another scratch, or made here.
                                if x_fence::x_fenced(x_fenced) {
                                    continue;
                                }
                                let x_mark = x_fence::x_fence_mark(x_fenced);
                                let x_greedy_result = match precomputed {
                                    Some((_, XSpecTraceback::Greedy(result))) => result,
                                    _ => greedy_gapped_alignment_with_traceback(
                                        q_seq_blastna,
                                        adjusted_subject,
                                        adjusted_subject.len(),
                                        trace_q_start,
                                        trace_s_start_adj,
                                        reward,
                                        penalty,
                                        gap_open,
                                        gap_extend,
                                        x_drop_trace,
                                        &mut *greedy_align_scratch,
                                    ),
                                };
                                // EXPERIMENT (LOSAT_X_BNTBRANGE): as for the DP traceback above.
                                x_fence::x_fence_note(x_fenced, x_mark);
                                if x_fence::x_fenced(x_fenced) {
                                    continue;
                                }
                                match x_greedy_result {
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:583-597
                                    // ```c
                                    // Blast_HSPUpdateWithTraceback(gap_align, hsp);
                                    //
                                    // if (!delete_hsp && !kGreedyTraceback) {
                                    //     Int4 align_length = 0;
                                    //     Blast_HSPGetNumIdentitiesAndPositives(..., &align_length, ...);
                                    //     delete_hsp = Blast_HSPTest(hsp, hit_options, align_length);
                                    // }
                                    // ```
                                    // Greedy traceback keeps score/coordinates/edit-script only here.
                                    // Identity statistics are recomputed later during the NCBI-style
                                    // reevaluation and Blast_HSPTestIdentityAndLength pass.
                                    Some((
                                        final_qs,
                                        final_qe,
                                        final_ss,
                                        final_se,
                                        score,
                                        aln_len,
                                        edit_ops,
                                    )) => (
                                        final_qs, final_qe, final_ss, final_se, score, 0usize,
                                        0usize, 0usize, 0usize, aln_len, edit_ops,
                                    ),
                                    None => {
                                        if let (Some(timing), Some(alignment_start)) =
                                            (timing_ref, alignment_start)
                                        {
                                            let alignment_elapsed_ns =
                                                alignment_start.elapsed().as_nanos() as u64;
                                            BlastnTiming::record_ns(
                                                &timing.traceback_alignment_ns,
                                                alignment_elapsed_ns,
                                            );
                                            BlastnTiming::record_ns(
                                                &timing.traceback_alignment_greedy_ns,
                                                alignment_elapsed_ns,
                                            );
                                        }
                                        continue;
                                    }
                                }
                            };
                            if let (Some(timing), Some(alignment_start)) =
                                (timing_ref, alignment_start)
                            {
                                let alignment_elapsed_ns =
                                    alignment_start.elapsed().as_nanos() as u64;
                                BlastnTiming::record_ns(
                                    &timing.traceback_alignment_ns,
                                    alignment_elapsed_ns,
                                );
                                if use_dp {
                                    BlastnTiming::record_ns(
                                        &timing.traceback_alignment_dp_ns,
                                        alignment_elapsed_ns,
                                    );
                                } else {
                                    BlastnTiming::record_ns(
                                        &timing.traceback_alignment_greedy_ns,
                                        alignment_elapsed_ns,
                                    );
                                }
                            }

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:598-600
                            // ```c
                            // Blast_HSPAdjustSubjectOffset(hsp, start_shift);
                            // ```
                            let hsp_build_start = if timing_enabled {
                                Some(std::time::Instant::now())
                            } else {
                                None
                            };
                            if start_shift != 0 {
                                final_ss = final_ss.saturating_add(start_shift);
                                final_se = final_se.saturating_add(start_shift);
                            }

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:583-597
                            // ```c
                            // Blast_HSPUpdateWithTraceback(gap_align, hsp);
                            // ...
                            // Blast_HSPGetNumIdentitiesAndPositives(query_nomask,
                            //        adjusted_subject, hsp, score_options, &align_length, sbp);
                            // ```
                            // Traceback scoring uses query_blk->sequence; identity/length
                            // reporting uses query_blk->sequence_nomask. In LOSATN's
                            // current soft-query-masking path these slices alias.
                            let (matches, mismatches, gaps, gap_letters, aln_len) = if use_dp {
                                let (matches, mismatches, gaps, gap_letters) = stats_from_edit_ops(
                                    q_seq_nomask_blastna,
                                    s_seq_blastna,
                                    final_qs,
                                    final_ss,
                                    &edit_ops,
                                );
                                (
                                    matches,
                                    mismatches,
                                    gaps,
                                    gap_letters,
                                    matches + mismatches + gap_letters,
                                )
                            } else {
                                (matches, mismatches, gaps, gap_letters, aln_len)
                            };

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:436-472
                            // ```c
                            // BlastGetOffsetsForGappedAlignment(..., &q_start, &s_start);
                            // ...
                            // BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
                            // ...
                            // AdjustSubjectRange(&s_start, &adjusted_s_length, q_start,
                            //                    query_length, &start_shift);
                            // ```
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:583-600
                            // ```c
                            // Blast_HSPUpdateWithTraceback(gap_align, hsp);
                            // ...
                            // Blast_HSPAdjustSubjectOffset(hsp, start_shift);
                            // ```
                            if blastn_trace_enabled
                                && (trace_traceback_seed
                                    || blastn_trace::should_trace_range(
                                        "traceback",
                                        prelim.context_idx,
                                        s_idx,
                                        s_id,
                                        final_qs,
                                        final_qe,
                                        final_ss,
                                        final_se,
                                        ctx.len,
                                        ctx.frame,
                                    ))
                            {
                                blastn_trace::log(
                            "traceback",
                            format!(
                                "subject={}({}) context={} start=({}, {}) adjusted_start=({}, {}) start_shift={} adjusted_s_len={} final=q{}..{} s{}..{} raw_score={} x_drop={} aln_len={} identities={} mismatches={} gaps={} gap_letters={} edit_ops_len={} edit_ops={}",
                                s_id,
                                s_idx,
                                prelim.context_idx,
                                trace_q_start,
                                trace_s_start,
                                trace_q_start,
                                trace_s_start_adj,
                                start_shift,
                                adjusted_s_len,
                                final_qs,
                                final_qe,
                                final_ss,
                                final_se,
                                score,
                                x_drop_trace,
                                aln_len,
                                matches,
                                mismatches,
                                gaps,
                                gap_letters,
                                edit_ops.len(),
                                format_gap_edit_ops_for_trace(&edit_ops)
                            ),
                        );
                            }

                            if timing_enabled {
                                traceback_edit_script_lengths
                                    .push(edit_ops.len().min(u32::MAX as usize) as u32);
                                traceback_alignment_lengths
                                    .push(aln_len.min(u32::MAX as usize) as u32);
                            }

                            if use_dp
                                && hsp_test(matches, aln_len, percent_identity, min_hit_length)
                            {
                                if let Some(timing) = timing_ref {
                                    BlastnTiming::record_count(
                                        &timing.traceback_deleted_identity_length_hsps,
                                        1,
                                    );
                                }
                                if let (Some(timing), Some(hsp_build_start)) =
                                    (timing_ref, hsp_build_start)
                                {
                                    BlastnTiming::record_duration(
                                        &timing.traceback_hsp_build_ns,
                                        hsp_build_start,
                                    );
                                }
                                continue;
                            }

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:583-612
                            // ```c
                            // Blast_HSPUpdateWithTraceback(gap_align, hsp);
                            // if (!delete_hsp && !kGreedyTraceback) {
                            //     Blast_HSPGetNumIdentitiesAndPositives(..., &align_length, ...);
                            //     delete_hsp = Blast_HSPTest(hsp, hit_options, align_length);
                            // }
                            // if (!delete_hsp) {
                            //     Blast_HSPAdjustSubjectOffset(hsp, start_shift);
                            //     status = BlastIntervalTreeAddHSP(hsp, tree, query_info,
                            //                                eQueryAndSubject);
                            // }
                            // ```
                            let final_tree_hsp = TreeHsp {
                                query_offset: final_qs as i32,
                                query_end: final_qe as i32,
                                subject_offset: final_ss as i32,
                                subject_end: final_se as i32,
                                score,
                                query_frame: prelim.query_frame,
                                query_length: ctx.len as i32,
                                query_context_offset: prelim.query_context_offset,
                                subject_frame_sign,
                            };
                            interval_tree.add_hsp(
                                final_tree_hsp,
                                prelim.query_context_offset,
                                IndexMethod::QueryAndSubject,
                            );

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:234-250
                            // ```c
                            // Blast_HSPListGetEvalues(program_number, query_info, subject_length,
                            //                         hsp_list, kGapped, FALSE, sbp, 0,
                            //                         scale_factor);
                            // ...
                            // Blast_HSPListGetBitScores(hsp_list, kGapped, sbp);
                            // ```
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1887-1890
                            // ```c
                            // hsp->evalue =
                            //     BLAST_KarlinStoE_simple(score, kbp[kbp_context],
                            //                      query_info->contexts[hsp->context].eff_searchsp);
                            // ```
                            let eff_searchsp = query_eff_searchsp[prelim.context_idx as usize];
                            let (bit_score, eval) = calculate_blastn_context_statistics(
                                score,
                                &search_karlin_ref[prelim.context_idx as usize].gapped,
                                eff_searchsp,
                                round_down_evalue_score,
                            );
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:234-246
                            // ```c
                            // Blast_HSPListGetEvalues(program_number, query_info, subject_length,
                            //                         hsp_list, kGapped, FALSE, sbp, 0,
                            //                         scale_factor);
                            // Blast_HSPListReapByEvalue(hsp_list, hit_params->options);
                            // Blast_HSPListGetBitScores(hsp_list, kGapped, sbp);
                            // ```
                            // NCBI reaps by E-value only after common-endpoint purging,
                            // ambiguity re-evaluation, score resort, and final containment.
                            // Keep this traceback HSP in the subject list so it can still
                            // participate in those survivor-set decisions.

                            let identity = if aln_len > 0 {
                                ((matches as f64 / aln_len as f64) * 100.0).min(100.0)
                            } else {
                                0.0
                            };

                            let query_length = query_lengths[prelim.query_idx as usize];
                            let (hit_q_start, hit_q_end, hit_s_start, hit_s_end) =
                                adjust_blastn_offsets(
                                    final_qs,
                                    final_qe,
                                    final_ss,
                                    final_se,
                                    query_length,
                                    prelim.query_frame,
                                );

                            let gap_info = if edit_ops.is_empty() {
                                None
                            } else {
                                Some(edit_ops)
                            };
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4058-4077
                            // ```c
                            // if (gap_align->score >= cutoff) {
                            //     ...
                            //     status = Blast_HSPInit(gap_align->query_start,
                            //                   gap_align->query_stop, gap_align->subject_start,
                            //                   gap_align->subject_stop,
                            //                   init_hsp->offsets.qs_offsets.q_off,
                            //                   init_hsp->offsets.qs_offsets.s_off, context,
                            //                   query_frame, subject->frame, gap_align->score,
                            //                   &(gap_align->edit_script), &new_hsp);
                            // }
                            // ```
                            hits.push(BlastnHsp {
                                identity,
                                length: aln_len,
                                mismatch: mismatches,
                                gapopen: gaps,
                                q_start: hit_q_start,
                                q_end: hit_q_end,
                                s_start: hit_s_start,
                                s_end: hit_s_end,
                                e_value: eval,
                                bit_score,
                                num_ident: matches,
                                query_frame: prelim.query_frame,
                                query_length,
                                q_idx: prelim.query_idx,
                                s_idx: s_idx as u32,
                                raw_score: score,
                                internal_q_offset_0: final_qs,
                                internal_q_end_0: final_qe,
                                internal_s_offset_0: final_ss,
                                internal_s_end_0: final_se,
                                internal_query_context_offset: prelim.query_context_offset,
                                gap_info,
                                num_positives: matches,
                            });
                            if let (Some(timing), Some(hsp_build_start)) =
                                (timing_ref, hsp_build_start)
                            {
                                BlastnTiming::record_duration(
                                    &timing.traceback_hsp_build_ns,
                                    hsp_build_start,
                                );
                            }
                        }
                    },
                );
                // No NCBI counterpart: LOSAT_X_STATS work counters (xstats feature); they do not change any
                // value NCBI computes.
                if let Some(x_ahead) = x_ahead.as_ref() {
                    let (x_evaluated, x_taken) = x_ahead.counts();
                    crate::utils::xstats::AHEAD_TRACE_EVALUATED
                        .fetch_add(x_evaluated, std::sync::atomic::Ordering::Relaxed);
                    crate::utils::xstats::AHEAD_TRACE_TAKEN
                        .fetch_add(x_taken, std::sync::atomic::Ordering::Relaxed);
                }
            }
        };

        // =================================================================
        // NCBI POST-GAPPED PROCESSING
        // Reference: blast_traceback.c:633-692
        // =================================================================

        // Step 1: Extract hits
        // The passes after the traceback loop up to the second purge, on one list: the
        // subject's, or one query's (LOSAT_X_QGROUPPAR). Returns the list and its lengths after
        // each step (for the debug line).
        let x_post_purge = |local_hits: Vec<BlastnHsp>| {
            // Step 2: Endpoint purging pass 1 (trim, purge=false)
            // NCBI reference: blast_traceback.c:637-638
            // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
            let hits_step1 = local_hits.len();
            let purge_pass1_start = if timing_enabled {
                Some(std::time::Instant::now())
            } else {
                None
            };
            let (mut local_hits, mut extra_start) =
                purge_hsps_with_common_endpoints_ex(local_hits, false);
            if let (Some(timing), Some(purge_pass1_start)) = (timing_ref, purge_pass1_start) {
                BlastnTiming::record_duration(
                    &timing.purge_common_endpoint_pass1_ns,
                    purge_pass1_start,
                );
                BlastnTiming::record_count(
                    &timing.traceback_removed_endpoint_pass1_hsps,
                    hits_step1.saturating_sub(local_hits.len()) as u64,
                );
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-647
            // ```c
            // Int4 extra_start =
            //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
            // ...
            // for (index=extra_start; index < hsp_list->hspcnt; index++) {
            // ```
            if blastn_trace_enabled
                && blastn_trace::should_trace_subject("purge", None, s_idx, s_id)
            {
                blastn_trace::log(
                    "purge",
                    format!(
                        "subject={}({}) pass=common_endpoint_trim before={} after={} extra_start={}",
                        s_id,
                        s_idx,
                        hits_step1,
                        local_hits.len(),
                        extra_start
                    ),
                );
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:640-644
            // ```c
            // /* Low level greedy algorithm ignores ambiguities, so the score
            //  * needs to be reevaluated. */
            // if (kGreedyTraceback) {
            //    extra_start = 0;
            // }
            // ```
            if !use_dp {
                extra_start = 0;
            }
            let hits_step2 = local_hits.len();

            // Step 3: Re-evaluate trimmed HSPs
            // NCBI reference: blast_traceback.c:647-665
            // The remaining part of the hsp may be extended further

            for hit in local_hits.iter_mut().skip(extra_start) {
                // Get sequences for re-evaluation
                let context_idx =
                    (hit.q_idx as usize) * 2 + if hit.query_frame < 0 { 1 } else { 0 };
                let cutoff = hit_saving_cutoff_scores
                    .get(context_idx)
                    .copied()
                    .unwrap_or_else(|| cutoff_scores.get(context_idx).copied().unwrap_or(0));

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1887-1891
                // ```c
                // hsp->evalue =
                //   BLAST_KarlinStoE_simple(score, kbp[kbp_context],
                //        query_info->contexts[hsp->context].eff_searchsp);
                // ```
                let eff_searchsp = query_eff_searchsp.get(context_idx).copied().unwrap_or(0);
                let reeval_params = ReevalParams {
                    lambda: search_karlin_ref[context_idx].gapped.lambda,
                    k: search_karlin_ref[context_idx].gapped.k,
                    eff_searchsp,
                    db_len: db_len_total,
                    db_num_seqs,
                    round_down_evalue_score,
                };

                // NCBI reference: blast_traceback.c:653-665 (reevaluate with blastna sequences)
                let q_seq_blastna = encoded_queries_blastna[context_idx];
                let q_seq_nomask_blastna = encoded_queries_blastna[context_idx];
                let s_seq_eval = s_seq_blastna;
                let reeval_start = if timing_enabled {
                    Some(std::time::Instant::now())
                } else {
                    None
                };
                let delete = reevaluate_hsp_with_ambiguities_gapped_ex(
                    hit,
                    q_seq_blastna,
                    s_seq_eval,
                    reward,
                    penalty,
                    gap_open,
                    gap_extend,
                    cutoff,
                    &score_matrix,
                    Some(&reeval_params),
                );
                if let (Some(timing), Some(reeval_start)) = (timing_ref, reeval_start) {
                    BlastnTiming::record_duration(&timing.reevaluate_ambiguities_ns, reeval_start);
                }
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:647-665
                // ```c
                // delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(...);
                // if (!delete_hsp)
                //     delete_hsp = Blast_HSPTestIdentityAndLength(...);
                // if (delete_hsp)
                //     hsp_array[index] = Blast_HSPFree(hsp);
                // ```
                let trace_reeval = blastn_trace_enabled
                    && blastn_trace::should_trace_range(
                        "purge",
                        context_idx as u32,
                        s_idx,
                        s_id,
                        hit.internal_q_offset_0,
                        hit.internal_q_end_0,
                        hit.internal_s_offset_0,
                        hit.internal_s_end_0,
                        hit.query_length,
                        hit.query_frame,
                    );
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:656-663
                // ```c
                // delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(...);
                // if (!delete_hsp)
                //     delete_hsp = Blast_HSPTestIdentityAndLength(program_number, hsp,
                //                                                 query_nomask, subject,
                //                                                 score_options, hit_options);
                // ```
                if delete {
                    if let Some(timing) = timing_ref {
                        BlastnTiming::record_count(&timing.traceback_deleted_reevaluation_hsps, 1);
                    }
                    if trace_reeval {
                        blastn_trace::log(
                            "purge",
                            format!(
                                "subject={}({}) context={} pass=reeval q{}..{} s{}..{} deleted=true reason=score_or_gap_info cutoff={}",
                                s_id,
                                s_idx,
                                context_idx,
                                hit.internal_q_offset_0,
                                hit.internal_q_end_0,
                                hit.internal_s_offset_0,
                                hit.internal_s_end_0,
                                cutoff
                            ),
                        );
                    }
                    hit.raw_score = i32::MIN; // Mark for removal
                    continue;
                }
                let identity_length_start = if timing_enabled {
                    Some(std::time::Instant::now())
                } else {
                    None
                };
                let delete_identity_length = blast_hsp_test_identity_and_length(
                    hit,
                    q_seq_nomask_blastna,
                    s_seq_eval,
                    percent_identity,
                    min_hit_length,
                );
                if let (Some(timing), Some(identity_length_start)) =
                    (timing_ref, identity_length_start)
                {
                    BlastnTiming::record_duration(
                        &timing.identity_length_test_ns,
                        identity_length_start,
                    );
                }
                if delete_identity_length {
                    if let Some(timing) = timing_ref {
                        BlastnTiming::record_count(
                            &timing.traceback_deleted_identity_length_hsps,
                            1,
                        );
                    }
                    if trace_reeval {
                        blastn_trace::log(
                            "purge",
                            format!(
                                "subject={}({}) context={} pass=identity_length q{}..{} s{}..{} deleted=true raw_score={} identities={} aln_len={} identity={:.6}",
                                s_id,
                                s_idx,
                                context_idx,
                                hit.internal_q_offset_0,
                                hit.internal_q_end_0,
                                hit.internal_s_offset_0,
                                hit.internal_s_end_0,
                                hit.raw_score,
                                hit.num_ident,
                                hit.length,
                                hit.identity
                            ),
                        );
                    }
                    hit.raw_score = i32::MIN; // Mark for removal
                } else if trace_reeval {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:647-665
                    // ```c
                    // delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(...);
                    // ```
                    blastn_trace::log(
                        "purge",
                        format!(
                            "subject={}({}) context={} pass=reeval q{}..{} s{}..{} deleted=false raw_score={} evalue={:.12e} bit_score={:.12} identities={} aln_len={} identity={:.6} gap_info={}",
                            s_id,
                            s_idx,
                            context_idx,
                            hit.internal_q_offset_0,
                            hit.internal_q_end_0,
                            hit.internal_s_offset_0,
                            hit.internal_s_end_0,
                            hit.raw_score,
                            hit.e_value,
                            hit.bit_score,
                            hit.num_ident,
                            hit.length,
                            hit.identity,
                            hit.gap_info
                                .as_deref()
                                .map(format_gap_edit_ops_for_trace)
                                .unwrap_or_else(|| "[]".to_string())
                        ),
                    );
                }
            }
            local_hits.retain(|h| h.raw_score != i32::MIN);
            let hits_step3 = local_hits.len();

            // Step 4: Endpoint purging pass 2 (delete, purge=true) - BLASTN only
            // NCBI reference: blast_traceback.c:667-668
            // if(program_number == eBlastTypeBlastn) {
            //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
            // }
            let purge_pass2_start = if timing_enabled {
                Some(std::time::Instant::now())
            } else {
                None
            };
            let mut local_hits = purge_hsps_with_common_endpoints(local_hits);
            let hits_step4 = local_hits.len();
            if let (Some(timing), Some(purge_pass2_start)) = (timing_ref, purge_pass2_start) {
                BlastnTiming::record_duration(
                    &timing.purge_common_endpoint_pass2_ns,
                    purge_pass2_start,
                );
                BlastnTiming::record_count(
                    &timing.traceback_removed_endpoint_pass2_hsps,
                    hits_step3.saturating_sub(hits_step4) as u64,
                );
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:666-668
            // ```c
            // Blast_HSPListPurgeNullHSPs(hsp_list);
            // if(program_number == eBlastTypeBlastn) {
            //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
            // }
            // ```
            if blastn_trace_enabled
                && blastn_trace::should_trace_subject("purge", None, s_idx, s_id)
            {
                blastn_trace::log(
                    "purge",
                    format!(
                        "subject={}({}) pass=post_reeval before={} after_reeval={} after_delete_purge={}",
                        s_id, s_idx, hits_step2, hits_step3, hits_step4
                    ),
                );
            }

            (local_hits, hits_step1, hits_step2, hits_step3, hits_step4)
        };
        // The passes after the second purge on one list: the score sort, the final tree and the
        // e-values. Returns the HSPs kept and their number after the tree (for the debug line).
        let x_final_list = |mut local_hits: Vec<BlastnHsp>,
                            interval_tree: &mut BlastIntervalTree| {
            // Step 5: Re-sort by gapped score
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:671-672
            // ```c
            // /* Sort HSPs by score again, as the scores might have changed. */
            // Blast_HSPListSortByScore(hsp_list);
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353 (ScoreCompareHSPs)
            // ```c
            // if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
            //     0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
            //     0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
            //     0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
            //     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
            // }
            // ```
            let score_resort_start = if timing_enabled {
                Some(std::time::Instant::now())
            } else {
                None
            };
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:671-672
            // ```c
            // /* Sort HSPs by score again, as the scores might have changed. */
            // Blast_HSPListSortByScore(hsp_list);
            // ```
            sort_hsps_by_score(&mut local_hits);
            if let (Some(timing), Some(score_resort_start)) = (timing_ref, score_resort_start) {
                BlastnTiming::record_duration(
                    &timing.traceback_score_resort_ns,
                    score_resort_start,
                );
            }

            // Step 6: Phase 2 - Tree reset and second containment pass
            // NCBI reference: blast_traceback.c:678
            // Blast_IntervalTreeReset(tree);
            interval_tree.reset();
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:66-70
            // ```c
            // if (tree->num_used == tree->num_alloc) {
            //     tree->num_alloc = 2 * tree->num_alloc;
            //     tree->nodes = (SIntervalNode *)realloc(tree->nodes, tree->num_alloc *
            //                                                  sizeof(SIntervalNode));
            // }
            // ```
            // The final containment pass adds at most one tree HSP for each surviving
            // HSP. Reserving here only removes allocator churn; containment and order
            // are still controlled by BlastIntervalTreeContainsHSP/AddHSP.
            interval_tree.reserve_nodes_for_hsps(local_hits.len());

            // NCBI reference: blast_traceback.c:679-692
            // Remove any HSPs that are contained within other HSPs.
            // Since the list is sorted by score already, any HSP
            // contained by a previous HSP is guaranteed to have a
            // lower score, and may be purged.
            let mut final_hits = Vec::with_capacity(local_hits.len());
            for hit in local_hits {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:537-550
                // ```c
                // query_start = query_context_offset;
                // region_start = query_start + hsp->query.offset;
                // region_end = query_start + hsp->query.end;
                // ```
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3913
                // ```c
                // tmp_hsp.query.offset = q_start;
                // tmp_hsp.query.end = q_end;
                // tmp_hsp.query.frame = query_info->contexts[context].frame;
                // tmp_hsp.subject.offset = s_start;
                // tmp_hsp.subject.end = s_end;
                // ```
                // NCBI uses canonical coordinates: subject.offset < subject.end always
                let query_context_offset = hit.internal_query_context_offset;
                let tree_hsp = TreeHsp {
                    query_offset: hit.internal_q_offset_0 as i32,
                    query_end: hit.internal_q_end_0 as i32,
                    subject_offset: hit.internal_s_offset_0 as i32,
                    subject_end: hit.internal_s_end_0 as i32,
                    score: hit.raw_score,
                    query_frame: hit.query_frame,
                    query_length: hit.query_length as i32,
                    query_context_offset,
                    subject_frame_sign: 1,
                };

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:678-688
                // ```c
                // Blast_IntervalTreeReset(tree);
                // for (index = 0; index < hsp_list->hspcnt; index++) {
                //     if (BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                //                                  hit_options->min_diag_separation)) {
                //         hsp_array[index] = Blast_HSPFree(hsp);
                //     } else {
                //         BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
                //     }
                // }
                // ```
                let final_tree_contains_start = if timing_enabled {
                    Some(std::time::Instant::now())
                } else {
                    None
                };
                let final_tree_containing_hsp = interval_tree.containing_hsp(
                    &tree_hsp,
                    query_context_offset,
                    min_diag_separation,
                );
                let final_tree_contains = final_tree_containing_hsp.is_some();
                if let (Some(timing), Some(final_tree_contains_start)) =
                    (timing_ref, final_tree_contains_start)
                {
                    BlastnTiming::record_duration(
                        &timing.traceback_final_tree_contains_ns,
                        final_tree_contains_start,
                    );
                }
                if blastn_trace_enabled
                    && blastn_trace::should_trace_range(
                        "purge",
                        (hit.q_idx * 2 + if hit.query_frame < 0 { 1 } else { 0 }) as u32,
                        s_idx,
                        s_id,
                        hit.internal_q_offset_0,
                        hit.internal_q_end_0,
                        hit.internal_s_offset_0,
                        hit.internal_s_end_0,
                        hit.query_length,
                        hit.query_frame,
                    )
                {
                    blastn_trace::log(
                        "purge",
                        format!(
                            "subject={}({}) context={} pass=final_interval_tree q{}..{} s{}..{} raw_score={} contained={} containing_hsp={:?} min_diag_separation={}",
                            s_id,
                            s_idx,
                            hit.q_idx * 2 + if hit.query_frame < 0 { 1 } else { 0 },
                            hit.internal_q_offset_0,
                            hit.internal_q_end_0,
                            hit.internal_s_offset_0,
                            hit.internal_s_end_0,
                            hit.raw_score,
                            final_tree_contains,
                            final_tree_containing_hsp,
                            min_diag_separation
                        ),
                    );
                }
                if !final_tree_contains {
                    let final_tree_add_start = if timing_enabled {
                        Some(std::time::Instant::now())
                    } else {
                        None
                    };
                    interval_tree.add_hsp(
                        tree_hsp,
                        query_context_offset,
                        IndexMethod::QueryAndSubject,
                    );
                    if let (Some(timing), Some(final_tree_add_start)) =
                        (timing_ref, final_tree_add_start)
                    {
                        BlastnTiming::record_duration(
                            &timing.traceback_final_tree_add_ns,
                            final_tree_add_start,
                        );
                    }
                    final_hits.push(hit);
                } else if let Some(timing) = timing_ref {
                    BlastnTiming::record_count(&timing.traceback_removed_final_tree_hsps, 1);
                }
                // else: HSP is contained within another, skip (implicit delete)
            }

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:678-692
            // ```c
            // Blast_IntervalTreeReset(tree);
            // for (index = 0; index < hsp_list->hspcnt; index++) {
            //    if (BlastIntervalTreeContainsHSP(...)) {
            //       hsp_array[index] = Blast_HSPFree(hsp);
            //    } else {
            //       BlastIntervalTreeAddHSP(...);
            //    }
            // }
            // ```
            // Keep the final-containment count separate from the later E-value reap.
            let final_tree_hit_count = final_hits.len();

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:234-246
            // ```c
            // Blast_HSPListGetEvalues(program_number, query_info, subject_length,
            //                         hsp_list, kGapped, FALSE, sbp, 0,
            //                         scale_factor);
            // Blast_HSPListReapByEvalue(hsp_list, hit_params->options);
            // Blast_HSPListGetBitScores(hsp_list, kGapped, sbp);
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1983-2003
            // ```c
            // cutoff = hit_options->expect_value;
            // for (index = 0; index < hsp_list->hspcnt; index++) {
            //    hsp = hsp_array[index];
            //    if (hsp->evalue > cutoff) {
            //       hsp_array[index] = Blast_HSPFree(hsp_array[index]);
            //    } else {
            //       if (index > hsp_cnt)
            //          hsp_array[hsp_cnt] = hsp_array[index];
            //       hsp_cnt++;
            //    }
            // }
            // ```
            // Apply final gapped E-value assignment and reaping after all
            // post-traceback survivor-set operations, matching
            // s_HSPListPostTracebackUpdate.
            let before_evalue_reap = final_hits.len();
            for hit in final_hits.iter_mut() {
                let context_idx =
                    (hit.q_idx as usize) * 2 + if hit.query_frame < 0 { 1 } else { 0 };
                let eff_searchsp = query_eff_searchsp.get(context_idx).copied().unwrap_or(0);
                let (bit_score, eval) = calculate_blastn_context_statistics(
                    hit.raw_score,
                    &search_karlin_ref[context_idx].gapped,
                    eff_searchsp,
                    round_down_evalue_score,
                );
                hit.bit_score = bit_score;
                hit.e_value = eval;
            }
            final_hits.retain(|hit| hsp_survives_evalue_reap(hit.e_value, evalue_threshold));
            if let Some(timing) = timing_ref {
                BlastnTiming::record_count(
                    &timing.traceback_deleted_evalue_cutoff_hsps,
                    before_evalue_reap.saturating_sub(final_hits.len()) as u64,
                );
            }

            (final_hits, final_tree_hit_count)
        };
        // Dispatch point of LOSAT_X_QGROUPPAR for the traceback: the kept HSPs are the queries'
        // hit lists one after another, so each query's are one slice; per query on the pool, the
        // final lists concatenated in query order and sorted by score (stable, so ties stay in
        // query order, as the serial sort leaves them: its input, the second purge's output, is
        // ordered by context first; see x_qgroup_par.rs), or the loop over the subject's list.
        let prelim_hits = &mut combined_prelim_hits;
        let x_trace_counts = (x_qgroup_batch
            && x_qgroup_par::in_query_order(prelim_hits.iter().map(|hit| hit.query_idx)))
        .then(|| x_qgroup_par::counts(prelim_hits.iter().map(|hit| hit.query_idx), queries.len()));
        let x_trace_order = x_trace_counts
            .as_deref()
            .and_then(|counts| x_qgroup_par::plan(counts, x_qgroup_par::Stage::Traceback));
        let (final_hits, [hits_step1, hits_step2, hits_step3, hits_step4, final_tree_hit_count]) =
            if let Some(order) = x_trace_order.as_ref() {
                let x_starts = x_qgroup_par::starts(x_trace_counts.as_deref().unwrap_or(&[]));
                let x_kept: &[PrelimHit] = prelim_hits;
                let x_tree_q_max = (query_concat_length + 1) as i32;
                let x_tree_s_max = (s_len_full + 1) as i32;
                let x_results = x_qgroup_par::run(
                    parallel_pool,
                    x_qgroup_par::Stage::Traceback,
                    order,
                    queries.len(),
                    || {
                        (
                            BlastIntervalTree::new(0, x_tree_q_max, 0, x_tree_s_max),
                            GapAlignScratch::new(),
                            GreedyAlignScratch::new(),
                        )
                    },
                    |(tree, gap, greedy), query| {
                        let query = query as usize;
                        let mut edit_lengths = Vec::new();
                        let mut alignment_lengths = Vec::new();
                        let mut hits = Vec::new();
                        x_trace_list(
                            &x_kept[x_starts[query]..x_starts[query + 1]],
                            tree,
                            gap,
                            greedy,
                            None,
                            None,
                            &mut hits,
                            &mut edit_lengths,
                            &mut alignment_lengths,
                        );
                        let (list, s1, s2, s3, s4) = x_post_purge(hits);
                        let (list, s5) = x_final_list(list, tree);
                        (list, [s1, s2, s3, s4, s5], edit_lengths, alignment_lengths)
                    },
                );
                let mut fast: Vec<BlastnHsp> = Vec::new();
                let mut steps = [0usize; 5];
                for (list, step_counts, edit_lengths, alignment_lengths) in
                    x_results.into_iter().flatten()
                {
                    fast.extend(list);
                    for (sum, step) in steps.iter_mut().zip(step_counts) {
                        *sum += step;
                    }
                    traceback_edit_script_lengths.extend(edit_lengths);
                    traceback_alignment_lengths.extend(alignment_lengths);
                }
                x_qgroup_par::sort_stable(parallel_pool, &mut fast, score_compare_hsps);
                (fast, steps)
            } else {
                x_trace_list(
                    prelim_hits,
                    interval_tree,
                    gap_scratch,
                    &mut subject_scratch.greedy_align_scratch,
                    Some(&mut subject_scratch.x_ahead_scratch),
                    Some(&mut subject_scratch.x_spec_scratch),
                    hits,
                    &mut traceback_edit_script_lengths,
                    &mut traceback_alignment_lengths,
                );
                let (list, s1, s2, s3, s4) = x_post_purge(hits.drain(..).collect());
                let (list, s5) = x_final_list(list, interval_tree);
                (list, [s1, s2, s3, s4, s5])
            };
        prelim_hits.clear();

        // Phase 2 debug output
        if debug_mode || blastn_debug {
            let phase2_tree_removed = hits_step4 - final_tree_hit_count;
            let evalue_removed = final_tree_hit_count.saturating_sub(final_hits.len());
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:234-246
            // ```c
            // Blast_HSPListReapByEvalue(hsp_list, hit_params->options);
            // ```
            // Report tree containment and the later E-value reap as separate
            // phases, preserving the same ordering as NCBI traceback cleanup.
            eprintln!(
                    "[DEBUG] Phase 2 breakdown: step1(extract)={}, step2(purge1)={} (-{}), step3(reeval)={} (-{}), step4(purge2)={} (-{}), step6(tree)={} (-{}), evalue_reap={} (-{})",
                    hits_step1,
                    hits_step2, hits_step1 - hits_step2,
                    hits_step3, hits_step2 - hits_step3,
                    hits_step4, hits_step3 - hits_step4,
                    final_tree_hit_count, phase2_tree_removed,
                    final_hits.len(), evalue_removed
                );
        }

        if !final_hits.is_empty() {
            // NCBI reference: blast_traceback.c:633-692 (post-gapped processing is per-subject)
            *subject_hits = Some(final_hits);
        }
        if let Some(timing) = timing_ref {
            timing.record_traceback_lengths(
                std::mem::take(&mut traceback_edit_script_lengths),
                std::mem::take(&mut traceback_alignment_lengths),
            );
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:358-373
        // ```c
        // /* Make sure the HSPs in the HSP list are sorted by score, as they should be. */
        // ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
        // tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
        //                               0, subject_length + 1);
        // ```
        if let Some(timing) = timing_ref {
            if let Some(traceback_start) = traceback_start {
                let elapsed_ns = traceback_start.elapsed().as_nanos() as u64;
                timing
                    .traceback_ns
                    .fetch_add(elapsed_ns, std::sync::atomic::Ordering::Relaxed);
                timing
                    .traceback_calls
                    .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
        }
        if let Some(bar) = progress_bar.as_ref() {
            bar.inc(1);
        }
    };

    // The batch's results, read after its subjects are searched.
    let batch_result = |hit_lists: Vec<Option<BlastnHitList>>,
                        chunk_prelim_lists: PrelimHitLists| QueryBatch {
        hit_lists,
        chunk_prelim_lists,
        good_init_extends: good_init_extends.load(std::sync::atomic::Ordering::Relaxed),
        searched: any_valid_context,
        query_masks: seq_data.query_masks.clone(),
        context_karlin: context_karlin.clone(),
        query_eff_searchsp: query_eff_searchsp.clone(),
        timing: String::new(),
    };
    // The records of a subject set and their letters (none for a database, which has
    // drivers of its own below).
    let prepared_ref = subjects.prepared();
    let subject_records_ref: &[FastaRecord] = prepared_ref.map_or(&[], |prepared| prepared.records);
    let prepared_subject = |s_idx: usize| {
        prepared_ref
            .expect("the subjects of a subject set")
            .subject(s_idx)
    };
    // One pass over the subjects in subject order: the preliminary stage gives each subject's
    // preliminary HSPs (`kept` is `None`), the traceback saves the final HSPs of the kept
    // preliminary HSPs in the final hit lists of `hitlist_size` subjects
    // (`Blast_HSPResultsInsertHSPList`).
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1427
    // ```c
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //        continue;
    //    }
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:1700-1707
    // ```c
    //                 batch->hsplist_array[hsplist_itr] = NULL;
    //                 if (hsp_list->hspcnt == 0) {
    //                     hsp_list = Blast_HSPListFree(hsp_list);
    //                 }
    //                 else {
    //                     Blast_HSPResultsInsertHSPList(thread_data->tld[tid]->results, hsp_list,
    //                                   hit_params->options->hitlist_size);
    //                 }
    // ```
    // Each subject's results go to `sink` in OID order, which fills the hit lists as NCBI's
    // `Blast_HitListUpdate` does, one HSP list at a time. One thread hands each subject over as
    // soon as it is searched, so only the HSPs the hit lists keep stay in memory. Two or more
    // threads search all subjects of the batch and then hand them over, or, with
    // `LOSAT_X_BNSTREAM`, hand each subject over in OID order while they search the next ones.
    let search_subjects = |kept: Option<&[Vec<PrelimHit>]>,
                           lookups: Option<&LookupTables>,
                           sink: &mut (dyn FnMut(usize, Option<Vec<BlastnHsp>>, Vec<PrelimHit>)
                                     + Send)| {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:991-1041
        // ```c
        // Int4 offset_array_size = GetOffsetArraySize(lookup_wrap);
        // ...
        // aux_struct->offset_pairs =
        //   (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
        // ```
        let new_state = || {
            (
                GapAlignScratch::new(),
                SubjectScratch::new(queries_ref.len(), offset_array_size, window_size),
            )
        };
        let search_one = |state: &mut (GapAlignScratch, SubjectScratch), s_idx: usize| {
            let (gap_scratch, subject_scratch) = state;
            let mut subject_hits: Option<Vec<BlastnHsp>> = None;
            let mut prelim_hits: Vec<PrelimHit> = Vec::new();
            process_subject(
                s_idx,
                prepared_subject(s_idx),
                gap_scratch,
                subject_scratch,
                lookups,
                kept.map(|kept| &kept[s_idx]),
                &mut subject_hits,
                &mut prelim_hits,
            );
            (subject_hits, prelim_hits)
        };
        let num_subjects = subject_records_ref.len();
        if use_parallel {
            #[cfg(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))]
            {
                if x_bnstream() {
                    crate::utils::xreach::x_reach!("LOSAT_X_BNSTREAM", true);
                    let threads = parallel_pool.threads().max(1);
                    let window = threads * X_BNSTREAM_SUBJECTS_PER_THREAD;
                    parallel_pool.install(|| {
                        x_bnstream::search_in_order(
                            num_subjects,
                            window,
                            false,
                            &new_state,
                            &search_one,
                            &mut |s_idx, (hits, prelim)| sink(s_idx, hits, prelim),
                        )
                    });
                } else {
                    let results: Vec<(Option<Vec<BlastnHsp>>, Vec<PrelimHit>)> = parallel_pool
                        .install(|| {
                            (0..num_subjects)
                                .into_par_iter()
                                .map_init(new_state, search_one)
                                .collect()
                        });
                    for (s_idx, (hits, prelim)) in results.into_iter().enumerate() {
                        sink(s_idx, hits, prelim);
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
            if x_bnstream() {
                crate::utils::xreach::x_reach!("LOSAT_X_BNSTREAM", false);
            }
            let mut state = new_state();
            for s_idx in 0..num_subjects {
                let (hits, prelim) = search_one(&mut state, s_idx);
                sink(s_idx, hits, prelim);
            }
        }
    };

    // NCBI's preliminary stage keeps, for each query, the preliminary HSPs of at most
    // `prelim_hitlist_size` subjects (a split batch merges those of its query chunks), and
    // its traceback reads only those.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1674-1783
    // ```c
    // BLAST_PreliminarySearchEngine(program_number, query, query_info,
    //    seq_src, gap_align, score_params, lookup_wrap, word_options,
    //    ext_params, hit_params, eff_len_params, psi_options,
    //    db_options, hsp_stream, diagnostics, interrupt_search,
    //    progress_info);
    // ...
    // BLAST_ComputeTraceback(program_number, hsp_stream, query, query_info,
    //    seq_src, gap_align, score_params, ext_params, hit_params,
    //    eff_len_params, db_options, psi_options, rps_info, pattern_blk,
    //    results, interrupt_search, progress_info);
    // ```
    // A database has its own drivers of the two passes (`database::search_blocks`,
    // `database::traceback_listed`), with a scratch per worker made as the subject set's
    // parallel pass makes them.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:991-1041
    // ```c
    // Int4 offset_array_size = GetOffsetArraySize(lookup_wrap);
    // ...
    // aux_struct->offset_pairs =
    //   (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
    // ```
    let new_scratch = || {
        (
            GapAlignScratch::new(),
            SubjectScratch::new(queries_ref.len(), offset_array_size, window_size),
        )
    };
    // An ungapped (`-ungapped`) search only (the other instance never calls them): the steps of
    // NCBI's preliminary stage that finish a subject's HSP list and that `process_subject` leaves
    // to its caller. The reap by the query coverage is made as the HSPs are given to the
    // collector (it takes nothing of the subject). A database sequence is read again, unpacked
    // (BLASTNA), when its HSP list is not empty only, and re-evaluated there.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1482
    // ```c
    //           if (hsp_list && hsp_list->hspcnt > 0) {
    //              ...  Blast_HSPListReevaluateUngapped(hsp_list, query, subject, ...
    // ```
    let x_ungapped_query_coverage = |hits: &mut Vec<PrelimHit>| {
        ungapped_reap_query_coverage(hits, args.qcov_hsp_perc.unwrap_or(0.0), |context_idx| {
            query_contexts_ref[context_idx].len
        });
    };
    let x_ungapped_database_stage =
        |database: &database::BlastnDatabase, ordinal: u32, hits: &mut Vec<PrelimHit>| {
            if hits.is_empty() {
                return Ok(());
            }
            database.with_blastna(ordinal, |s_blastna| {
                ungapped_subject_stage(
                    hits,
                    s_blastna,
                    encoded_queries_blastna.as_slice(),
                    &score_matrix,
                    &x_batch_cutoffs.0,
                    percent_identity,
                    |prelim| {
                        let eff_searchsp = query_eff_searchsp_ref
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
                    evalue_threshold,
                );
                x_ungapped_query_coverage(hits);
            })
        };
    // The callback of the database's preliminary pass in an ungapped search: the search of the
    // sequence as for a default search, then the finish of its HSP list when it has HSPs left
    // (`x_ungapped_database_stage`). A read error goes to the cell `$error` (the one of the
    // smallest ordinal stays) and the sequence has no HSPs. The default search has its own
    // callback (a different closure, so that the default instance of this function holds no
    // capture of this one).
    macro_rules! ungapped_database_callback {
        ($database:expr, $error:expr) => {
            |state: &mut (GapAlignScratch, SubjectScratch),
             ordinal: u32,
             subject: BlastnSubjectRef<'_>| {
                let (gap_scratch, subject_scratch) = state;
                let mut subject_hits: Option<Vec<BlastnHsp>> = None;
                let mut prelim_hits: Vec<PrelimHit> = Vec::new();
                process_subject(
                    ordinal as usize,
                    subject,
                    gap_scratch,
                    subject_scratch,
                    Some(&lookup_tables),
                    None,
                    &mut subject_hits,
                    &mut prelim_hits,
                );
                if let Err(error) = x_ungapped_database_stage($database, ordinal, &mut prelim_hits)
                {
                    let mut first = $error.lock().expect("the error cell is not poisoned");
                    if first
                        .as_ref()
                        .is_none_or(|(smallest, _)| ordinal < *smallest)
                    {
                        *first = Some((ordinal, error));
                    }
                    return None;
                }
                (!prelim_hits.is_empty()).then_some(prelim_hits)
            }
        };
    }
    let prelim_lists = match split_prelim_lists {
        Some(lists) => lists,
        None => match subjects {
            // `-culling_limit`: NCBI's culling writer takes the place of the HSP collector
            // (`BlastHSPWriter` of `CreateHspWriter`), for a batch and for a query chunk alike
            // (a chunk's queries are its parts, so the lengths are the chunk's). Chosen once;
            // the default arms below are unchanged.
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/setup_factory.cpp:318-342
            // ```c
            //     else if (filt_opts) {
            //         if (filt_opts->best_hit && (filt_opts->best_hit_stage & ePrelimSearch))
            //         ...
            //         else if (filt_opts->culling_opts &&
            //                  (filt_opts->culling_stage & ePrelimSearch))
            //         {
            // ```
            SubjectSource::Prepared(_) if args.culling_limit_value() > 0 => {
                let mut writer = culling::PrelimCulling::new(
                    &query_lengths,
                    prelim_hitlist_size,
                    args.culling_limit_value(),
                );
                if UNGAPPED {
                    search_subjects(None, Some(&lookup_tables), &mut |oid, _, mut prelim| {
                        x_ungapped_query_coverage(&mut prelim);
                        writer.add(oid, prelim)
                    });
                } else {
                    search_subjects(None, Some(&lookup_tables), &mut |oid, _, prelim| {
                        writer.add(oid, prelim)
                    });
                }
                writer.finish()
            }
            SubjectSource::Database(database) if args.culling_limit_value() > 0 => {
                // The same drivers as the default arm below, with the culling writer as the
                // sink (a second instance of `database::search_blocks`).
                let mut writer = culling::PrelimCulling::new(
                    &query_lengths,
                    prelim_hitlist_size,
                    args.culling_limit_value(),
                );
                let mut add = |ordinal: u32, prelim_hits: Vec<PrelimHit>| {
                    writer.add(ordinal as usize, prelim_hits)
                };
                if any_valid_context {
                    let searched = if UNGAPPED {
                        // The first error of reading a database sequence again for the
                        // ungapped re-evaluation (`BlastnDatabase::with_blastna`).
                        let x_ungapped_error: std::sync::Mutex<Option<(u32, anyhow::Error)>> =
                            std::sync::Mutex::new(None);
                        database::search_blocks(
                            database,
                            parallel_pool,
                            use_parallel,
                            &new_scratch,
                            &ungapped_database_callback!(database, x_ungapped_error),
                            &mut add,
                        )
                        .and_then(|()| {
                            match x_ungapped_error
                                .lock()
                                .expect("the error cell is not poisoned")
                                .take()
                            {
                                Some((_, error)) => Err(error),
                                None => Ok(()),
                            }
                        })
                    } else {
                        database::search_blocks(
                            database,
                            parallel_pool,
                            use_parallel,
                            &new_scratch,
                            &|state: &mut (GapAlignScratch, SubjectScratch),
                              ordinal: u32,
                              subject: BlastnSubjectRef<'_>| {
                                let (gap_scratch, subject_scratch) = state;
                                let mut subject_hits: Option<Vec<BlastnHsp>> = None;
                                let mut prelim_hits: Vec<PrelimHit> = Vec::new();
                                process_subject(
                                    ordinal as usize,
                                    subject,
                                    gap_scratch,
                                    subject_scratch,
                                    Some(&lookup_tables),
                                    None,
                                    &mut subject_hits,
                                    &mut prelim_hits,
                                );
                                (!prelim_hits.is_empty()).then_some(prelim_hits)
                            },
                            &mut add,
                        )
                    };
                    searched.map_err(|error| {
                        database_search_error(
                            error,
                            first_batch,
                            reborrow_outputs(&mut outputs),
                            megablast,
                            &shown,
                        )
                    })?;
                }
                writer.finish()
            }
            SubjectSource::Prepared(_) => {
                let mut collector = PrelimHitCollector::new(queries_ref.len(), prelim_hitlist_size);
                if UNGAPPED {
                    search_subjects(None, Some(&lookup_tables), &mut |oid, _, mut prelim| {
                        x_ungapped_query_coverage(&mut prelim);
                        collector.add(oid, prelim)
                    });
                } else {
                    search_subjects(None, Some(&lookup_tables), &mut |oid, _, prelim| {
                        collector.add(oid, prelim)
                    });
                }
                collector.hit_lists
            }
            SubjectSource::Database(database) => {
                // The subjects come one at a time in the order of their ordinal ids
                // (`database::search_blocks`), as the subject set's come in OID order.
                let mut collector = PrelimHitCollector::new(queries_ref.len(), prelim_hitlist_size);
                // NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:177-180
                // ```c
                //     int status = m_PrelimSearch->CheckInternalData();
                //     if (status != 0)
                //     {
                //          // Search was not run, but we send back an empty CSearchResultSet.
                // ```
                // With no valid context NCBI does not search, and reads no sequence.
                let mut add = |ordinal: u32, prelim_hits: Vec<PrelimHit>| {
                    collector.add(ordinal as usize, prelim_hits)
                };
                if any_valid_context {
                    let searched = if UNGAPPED {
                        // The first error of reading a database sequence again for the
                        // ungapped re-evaluation (`BlastnDatabase::with_blastna`).
                        let x_ungapped_error: std::sync::Mutex<Option<(u32, anyhow::Error)>> =
                            std::sync::Mutex::new(None);
                        database::search_blocks(
                            database,
                            parallel_pool,
                            use_parallel,
                            &new_scratch,
                            &ungapped_database_callback!(database, x_ungapped_error),
                            &mut add,
                        )
                        .and_then(|()| {
                            match x_ungapped_error
                                .lock()
                                .expect("the error cell is not poisoned")
                                .take()
                            {
                                Some((_, error)) => Err(error),
                                None => Ok(()),
                            }
                        })
                    } else {
                        database::search_blocks(
                            database,
                            parallel_pool,
                            use_parallel,
                            &new_scratch,
                            &|state: &mut (GapAlignScratch, SubjectScratch),
                              ordinal: u32,
                              subject: BlastnSubjectRef<'_>| {
                                let (gap_scratch, subject_scratch) = state;
                                let mut subject_hits: Option<Vec<BlastnHsp>> = None;
                                let mut prelim_hits: Vec<PrelimHit> = Vec::new();
                                process_subject(
                                    ordinal as usize,
                                    subject,
                                    gap_scratch,
                                    subject_scratch,
                                    Some(&lookup_tables),
                                    None,
                                    &mut subject_hits,
                                    &mut prelim_hits,
                                );
                                (!prelim_hits.is_empty()).then_some(prelim_hits)
                            },
                            &mut add,
                        )
                    };
                    searched.map_err(|error| {
                        database_search_error(
                            error,
                            first_batch,
                            reborrow_outputs(&mut outputs),
                            megablast,
                            &shown,
                        )
                    })?;
                }
                collector.hit_lists
            }
        },
    };
    if let BatchStage::ChunkPrelim { .. } = stage {
        return Ok(batch_result(Vec::new(), prelim_lists));
    }
    // Dispatch point of LOSAT_X_BNLUTDROP: the traceback reads the kept HSPs only, so the tables
    // are dropped here (NCBI's traceback stage, traceback_stage.cpp:208-211); with `=0` they are
    // held until the batch returns (`x_lut_drop.rs`).
    let _x_kept_lookup_tables = x_lut_drop::after_prelim(lookup_tables);
    let hit_lists = match subjects {
        SubjectSource::Prepared(_) => {
            let kept = prelim_hits_by_subject(prelim_lists, subject_records_ref.len());
            let mut hit_lists: Vec<Option<BlastnHitList>> = Vec::with_capacity(queries_ref.len());
            hit_lists.resize_with(queries_ref.len(), || None);
            search_subjects(Some(&kept), None, &mut |_, hits, _| {
                if let Some(hits) = hits {
                    update_hitlists_with_subject_hits(&mut hit_lists, hits, hitlist_size);
                }
            });
            hit_lists
        }
        SubjectSource::Database(database) => {
            let listed = listed_prelim_hits(prelim_lists);
            let found = database::traceback_listed(
                database,
                parallel_pool,
                use_parallel,
                &listed,
                &new_scratch,
                &|state: &mut (GapAlignScratch, SubjectScratch),
                  ordinal: u32,
                  subject: BlastnSubjectRef<'_>,
                  kept: &Vec<PrelimHit>| {
                    let (gap_scratch, subject_scratch) = state;
                    let mut subject_hits: Option<Vec<BlastnHsp>> = None;
                    let mut prelim_hits: Vec<PrelimHit> = Vec::new();
                    process_subject(
                        ordinal as usize,
                        subject,
                        gap_scratch,
                        subject_scratch,
                        None,
                        Some(kept),
                        &mut subject_hits,
                        &mut prelim_hits,
                    );
                    subject_hits
                },
                &|kept: &Vec<PrelimHit>| {
                    kept.iter()
                        .map(|hit| (hit.prelim_ss, hit.prelim_se))
                        .collect()
                },
            )
            .map_err(|error| {
                database_search_error(
                    error,
                    first_batch,
                    reborrow_outputs(&mut outputs),
                    megablast,
                    &shown,
                )
            })?;
            // The final HSPs of the listed sequences in the order of their ordinal ids, as
            // the subject set's traceback saves those of its subjects
            // (`Blast_HSPResultsInsertHSPList`, `search_subjects`).
            let mut hit_lists: Vec<Option<BlastnHitList>> = Vec::with_capacity(queries_ref.len());
            hit_lists.resize_with(queries_ref.len(), || None);
            for hits in found {
                update_hitlists_with_subject_hits(&mut hit_lists, hits, hitlist_size);
            }
            hit_lists
        }
    };

    if let Some(bar) = progress_bar.as_ref() {
        bar.finish();
    }
    let mut batch = batch_result(hit_lists, Vec::new());
    if let Some(timing) = timing.as_ref() {
        batch.timing = blastn_timing_text(timing.as_ref(), t_search_start, t_total);
    }
    Ok(batch)
}

/// NCBI's preliminary search of a split query batch: each query chunk is searched as a batch
/// of its query parts (`BatchStage::ChunkPrelim`), and its preliminary hit lists are mapped
/// onto the batch's contexts and merged with those of the chunks before it
/// (`merge_query_chunk`). The merged hit lists of the batch's queries; `None` when NCBI does
/// not split the batch.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:237-289
/// ```c
///         for (Uint4 i = 0; i < query_splitter->GetNumberOfChunks(); i++) {
///             try {
///                 CRef<IQueryFactory> chunk_qf =
///                     query_splitter->GetQueryFactoryForChunk(i);
///                 ...
///                 CRef<SInternalData> chunk_data =
///                     SplitQuery_CreateChunkData(chunk_qf, m_Options,
///                                                m_InternalData,
///                                                GetNumberOfThreads());
///                 ...
///                     retval =
///                         CPrelimSearchRunner(*chunk_data, opts_memento.get())();
///                 ...
///                 BlastHSPStreamMerge(split_query_blk->GetCStruct(), i,
///                                 chunk_data->m_HspStream->GetPointer(),
///                                 m_InternalData->m_HspStream->GetPointer());
///                 ...
///             } catch (const CBlastException& e) {
///                 // This error message is safe to ignore for a given chunk,
///                 // because the chunks might end up producing a region of
///                 // the query for which ungapped Karlin-Altschul blocks
///                 // cannot be calculated
///                 ...
///                 if (e.GetMsg().find(err_msg1) == NPOS && e.GetMsg().find(err_msg2) == NPOS) {
///                     throw;
///                 }
/// ```
/// A chunk has its own query block, masks, Karlin-Altschul blocks, lookup table and
/// diagnostics (a split batch gives the next batch size 0 initial hits), and NCBI drops its
/// warnings. A chunk without a valid context finds nothing, as NCBI ignores its error.
#[allow(clippy::too_many_arguments)]
fn search_query_chunks(
    args: &BlastnArgs,
    queries: &x_query_drop::Letters<'_>,
    query_masks: &[Vec<MaskedInterval>],
    query_offsets: &[usize],
    batch_eff_searchsp: &[i64],
    context_offsets: &[i32],
    subjects: SubjectSource<'_>,
    mut outputs: BatchOutputs<'_, '_>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batching: QueryBatching,
    first_batch: bool,
    x_keep: Option<&Arc<x_mb_lut_keep::Keep>>,
) -> Result<Option<PrelimHitLists>> {
    let lengths: Vec<usize> = (0..queries.count())
        .map(|index| queries.length(index))
        .collect();
    // A negative CHUNK_SIZE splits a batch only with a negative OVERLAP_CHUNK_SIZE below it
    // (`calculate_num_chunks`); NCBI's `size_t` chunk ranges then wrap (split_query_cxx.cpp:
    // 145-171): NCBI stops with a CCoreException that names its build's files where a chunk
    // would be split again, and otherwise searches chunk ranges with gaps between them.
    if batching.split.negative_chunk_size()
        && calculate_num_chunks(batching.split, lengths.iter().sum()).0 > 1
    {
        anyhow::bail!(
            "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this query batch with NCBI BLAST+'s size_t chunk ranges wrapped (NCBI stops with a CCoreException that names its build's source files, or searches chunks with gaps between them), which is not supported by LOSAT's BLASTN"
        );
    }
    let Some(chunks) = split_query_batch(
        &lengths,
        batching.split,
        QueryStrand::from_option(&args.strand),
    ) else {
        return Ok(None);
    };
    // With an overlap close to the chunk size NCBI's setup of a chunk would split it again
    // and stop with a CCoreException that names its build's files
    // (`query_split::chunk_would_be_split`). LOSAT searches each chunk once, as NCBI does
    // for every overlap that does not split a chunk again (approved exception 1 of
    // PD-LOSAT-NCBI-DEFECTS).
    let mut merged: PrelimHitLists = Vec::with_capacity(queries.count());
    merged.resize_with(queries.count(), || None);
    for chunk in &chunks {
        // NCBI fails to make the query factory of a chunk without a query ("Empty
        // CBlastQueryVector"); `split_query_batch` makes none.
        anyhow::ensure!(
            !chunk.queries.is_empty(),
            "a query chunk without a query is not supported by LOSAT's BLASTN"
        );
        let parts: Vec<FastaRecord> = chunk
            .queries
            .iter()
            .map(|part| queries.cut(part.query, part.from, part.to))
            .collect();
        let restricted = chunk
            .queries
            .iter()
            .map(|part| restrict_masks(&query_masks[part.query], part, query_offsets[part.query]))
            .collect();
        let part_masks = chunk_query_masks(args, &parts, restricted, parallel_pool)
            .map_err(|error| window_masker_rejection(args, &error))?;
        let batch = search_query_batch::<false>(
            args,
            x_query_drop::BatchQueries::Shared(&parts),
            subjects,
            reborrow_outputs(&mut outputs),
            parallel_pool,
            batching,
            first_batch,
            BatchStage::ChunkPrelim {
                query_masks: &part_masks,
                batch_eff_searchsp,
            },
            x_keep,
        )?;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:174-178
        // ```c
        //     const size_t kOverlap =
        //         Blast_QueryIsTranslated(m_Options->GetProgramType())
        //         ? kOverlapSize / CODON_LENGTH : kOverlapSize;
        //     m_SplitBlk->SetChunkOverlapSize(kOverlap);
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:519-522
        // ```c
        //        Blast_HitListMerge(results1->hitlist_array + i,
        //                           results2->hitlist_array + global_query,
        //                           contexts_per_query, split_points,
        //                           (Int4)SplitQueryBlk_GetChunkOverlapSize(squery_blk),
        // ```
        merge_query_chunk(
            &mut merged,
            batch.chunk_prelim_lists,
            chunk,
            &lengths,
            context_offsets,
            batching.split.overlap as u32 as i32,
        )?;
    }
    Ok(Some(merged))
}

/// NCBI's `BlastHSPStreamMerge` of a query chunk: the HSPs of each query part's hit list
/// move to the part's contexts in the batch, and the list merges with the query's hit list of
/// the chunks before (`merge_prelim_hit_list`).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:476-522
/// ```c
///        for (j = 0; j < contexts_per_query; j++) {
///            split_points[j] = -1;
///        }
///
///        for (j = 0; j < contexts_per_query; j++) {
///            Int4 local_context = i * contexts_per_query + j;
///            if (context_list[local_context] >= 0) {
///                split_points[context_list[local_context] % contexts_per_query] =
///                                 offset_list[local_context];
///            }
///        }
///        ...
///                hsp->context = context_list[local_context];
///                hsp->query.offset += offset_list[local_context];
///                hsp->query.end += offset_list[local_context];
///                hsp->query.gapped_start += offset_list[local_context];
///                hsp->query.frame = BLAST_ContextToFrame(stream2->program,
///                                                        hsp->context);
///            }
///
///            hsplist->query_index = global_query;
///        }
///
///        Blast_HitListMerge(results1->hitlist_array + i,
///                           results2->hitlist_array + global_query,
///                           contexts_per_query, split_points,
///                           (Int4)SplitQueryBlk_GetChunkOverlapSize(squery_blk),
/// ```
/// The overlap is 100 and gaps are allowed (`CSplitQueryBlk` of a gapped search).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:527-534
/// ```c
///    for (i = 0; i < results2->num_queries; i++) {
///        BlastHitList *hitlist = results2->hitlist_array[i];
///        if (hitlist == NULL)
///            continue;
///
///        for (j = 0; j < hitlist->hsplist_count; j++)
///            Blast_HSPListSortByScore(hitlist->hsplist_array[j]);
///    }
/// ```
fn merge_query_chunk(
    merged: &mut PrelimHitLists,
    chunk_lists: PrelimHitLists,
    chunk: &QueryChunk,
    query_lengths: &[usize],
    context_offsets: &[i32],
    chunk_overlap_size: i32,
) -> Result<()> {
    for (local_query, hit_list) in chunk_lists.into_iter().enumerate() {
        let Some(mut hit_list) = hit_list else {
            continue;
        };
        let part = &chunk.queries[local_query];
        let query_length = query_lengths[part.query] as i64;
        for list in &mut hit_list.hsplist_array {
            for hsp in &mut list.hsps {
                let local_context = hsp.context_idx as usize;
                // The chunk's HSPs are in the contexts that it searches.
                let context = chunk.contexts[local_context].ok_or_else(|| {
                    anyhow::anyhow!(
                        "internal error: an HSP of a query chunk is in a context that is not searched"
                    )
                })?;
                let offset = chunk.context_offsets[local_context] as i64;
                let query_offset = hsp.prelim_qs as i64 + offset;
                let query_end = hsp.prelim_qe as i64 + offset;
                let gapped_start = hsp.seed_qs as i64 + offset;
                // The offsets are the starts of the parts (`query_split::context_offsets`).
                anyhow::ensure!(
                    query_offset >= 0 && query_end <= query_length && gapped_start < query_length,
                    "internal error: an HSP of a query chunk maps outside its query"
                );
                hsp.context_idx = context as u32;
                hsp.query_idx = part.query as u32;
                hsp.query_frame = if context % 2 == 0 { 1 } else { -1 };
                hsp.query_context_offset = context_offsets[context];
                hsp.prelim_qs = query_offset as usize;
                hsp.prelim_qe = query_end as usize;
                hsp.seed_qs = gapped_start as usize;
            }
        }
        // A context that the chunk does not have (a strand that is not searched) keeps -1.
        let split_point = |strand: usize| {
            let position = 2 * local_query + strand;
            if chunk.contexts[position].is_some() {
                chunk.context_offsets[position]
            } else {
                -1
            }
        };
        let split_points = [split_point(0), split_point(1)];
        merge_prelim_hit_list(
            hit_list,
            &mut merged[part.query],
            split_points,
            chunk_overlap_size,
        );
    }
    for hit_list in merged.iter_mut().flatten() {
        for list in &mut hit_list.hsplist_array {
            sort_prelim_hits_by_score(&mut list.hsps);
        }
    }
    Ok(())
}

/// NCBI's `Blast_HitListMerge`: a query chunk's hit list of a query merges with the query's
/// hit list of the chunks before into a new hit list of the same size, subject by subject in
/// subject order (`Blast_HitListUpdate` keeps the best `prelim_hitlist_size` subjects by
/// their preliminary e-values). The HSP lists of a subject in both merge with
/// `Blast_HSPListsMerge` when a split point is positive, else `Blast_HSPListAppend`.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2132-2217
/// ```c
///     if (hitlist1 == NULL)
///         return 0;
///     if (hitlist2 == NULL) {
///         *combined_hit_list_ptr = hitlist1;
///         *old_hit_list_ptr = NULL;
///         return 0;
///     }
///     num_hsplists1 = hitlist1->hsplist_count;
///     num_hsplists2 = hitlist2->hsplist_count;
///     new_hitlist = Blast_HitListNew(hitlist1->hsplist_max);
///
///     /* sort the lists of HSPs by oid */
///
///     if (num_hsplists1 > 1) {
///         qsort(hitlist1->hsplist_array, num_hsplists1,
///               sizeof(BlastHSPList*), s_SortHSPListByOid);
///     }
///     ...
///     query_is_split = FALSE;
///     for (i = 0; i < contexts_per_query; i++) {
///         if (split_offsets[i] > 0) {
///             query_is_split = TRUE;
///     ...
///         if (hsplist1->oid < hsplist2->oid) {
///             Blast_HitListUpdate(new_hitlist, hsplist1);
///             i++;
///         }
///         else if (hsplist1->oid > hsplist2->oid) {
///             Blast_HitListUpdate(new_hitlist, hsplist2);
///             j++;
///         }
///         else {
///             ...
///             if (query_is_split) {
///                 Blast_HSPListsMerge(hitlist1->hsplist_array + i,
///                                     hitlist2->hsplist_array + j,
///                                     hsplist2->hsp_max, split_offsets,
///                                     contexts_per_query,
///                                     chunk_overlap_size,
///                                     allow_gap, FALSE);
///             }
///             else {
///                 Blast_HSPListAppend(hitlist1->hsplist_array + i,
///                                     hitlist2->hsplist_array + j,
///                                     hsplist2->hsp_max);
///             }
///             Blast_HitListUpdate(new_hitlist, hitlist2->hsplist_array[j]);
///     ...
///     *old_hit_list_ptr = NULL;
///     *combined_hit_list_ptr = new_hitlist;
/// ```
fn merge_prelim_hit_list(
    mut hitlist1: HitList<PrelimHspList>,
    combined: &mut Option<HitList<PrelimHspList>>,
    split_offsets: [i32; 2],
    chunk_overlap_size: i32,
) {
    let Some(mut hitlist2) = combined.take() else {
        *combined = Some(hitlist1);
        return;
    };
    let mut new_hitlist = HitList::new(hitlist1.hsplist_max);
    // The subjects of a hit list are distinct, so the sort order is total.
    hitlist1.hsplist_array.sort_by_key(|list| list.oid);
    hitlist2.hsplist_array.sort_by_key(|list| list.oid);
    let query_is_split = split_offsets.iter().any(|&offset| offset > 0);
    let mut lists1 = std::mem::take(&mut hitlist1.hsplist_array)
        .into_iter()
        .peekable();
    let mut lists2 = std::mem::take(&mut hitlist2.hsplist_array)
        .into_iter()
        .peekable();
    while let (Some(list1), Some(list2)) = (lists1.peek(), lists2.peek()) {
        if list1.oid < list2.oid {
            new_hitlist.update(lists1.next().unwrap());
        } else if list1.oid > list2.oid {
            new_hitlist.update(lists2.next().unwrap());
        } else {
            let list1 = lists1.next().unwrap();
            let mut list2 = lists2.next().unwrap();
            if query_is_split {
                merge_prelim_hit_lists(
                    &mut list2.hsps,
                    list1.hsps,
                    HspListSplit::Query {
                        offsets: split_offsets,
                    },
                    chunk_overlap_size as usize,
                    true,
                );
            } else {
                append_prelim_hit_list(&mut list2.hsps, list1.hsps);
            }
            new_hitlist.update(list2);
        }
    }
    for list in lists1.chain(lists2) {
        new_hitlist.update(list);
    }
    *combined = Some(new_hitlist);
}

#[cfg(test)]
mod tests {
    use super::*;

    // QueryContextIndex answers the value of `direct_map` at every offset, through whichever form
    // the switches choose (by default LOSAT_X_CTXO1: division, or the block table, here with the
    // runs for its crowded blocks in the second case: a query of 3 letters among 500 of 300).
    #[test]
    fn the_context_index_answers_the_table_at_every_offset() {
        let mut crowded = vec![300usize; 1000];
        crowded[200] = 3;
        crowded[201] = 3;
        crowded.extend([0, 0]);
        for (case, lengths) in [
            vec![24usize; 40],
            crowded,
            vec![0, 0, 18, 18, 30, 30, 25, 25, 19, 19],
            vec![5000, 5000],
        ]
        .into_iter()
        .enumerate()
        {
            let mut contexts = Vec::new();
            let mut offset = 0usize;
            for (i, &len) in lengths.iter().enumerate() {
                contexts.push(QueryContext {
                    query_idx: (i / 2) as u32,
                    frame: if i % 2 == 0 { 1 } else { -1 },
                    query_offset: offset as i32,
                    len,
                    masks: Vec::new(),
                });
                offset += if len == 0 { 0 } else { len + 1 };
            }
            let total = contexts
                .iter()
                .map(|c| c.query_offset as usize + c.len)
                .max()
                .unwrap();
            let mut table = vec![0usize; total];
            for (idx, c) in contexts.iter().enumerate() {
                let start = c.query_offset as usize;
                table[start..start + c.len].fill(idx);
            }
            let index = QueryContextIndex::new(&contexts, false);
            if case == 1 {
                if let Some(o1) = &index.x_o1 {
                    assert!(o1.needs_runs() && !index.x_ranges.is_empty());
                }
            }
            for (n, &want) in table.iter().enumerate() {
                assert_eq!(
                    index.context_for_offset(n),
                    want,
                    "offset {n} of {lengths:?}"
                );
            }
            // past the table: NCBI's search, as before the forms
            for n in total..total + 40 {
                assert_eq!(index.context_for_offset(n), bsearch_context_info(n, &index));
            }
        }
    }

    /// The records of a FASTA text as NCBI's reader reads them (`config`).
    fn reader_records(fasta: &str, config: ReaderConfig) -> Vec<FastaRecord> {
        crate::blastinput::fasta_reader::read_all(
            &mut FastaInputSource::from_bytes(fasta.as_bytes(), config),
            &mut |_| Ok(()),
        )
        .unwrap()
    }

    /// The query buffer as it was built from per-context IUPAC copies: every context encoded
    /// on its own, then copied into the block and into the block with sentinels.
    fn query_buffers_from_context_copies(
        queries: &[&[u8]],
        offsets: &[i32],
        concat_len: usize,
    ) -> (Vec<u8>, Vec<u8>) {
        let mut logical = vec![NUCL_SENTINEL; concat_len];
        let mut with_sentinels = vec![NUCL_SENTINEL; concat_len + 2];
        for (q, seq) in queries.iter().enumerate() {
            let contexts = [
                encode_iupac_to_blastna(seq),
                encode_iupac_to_blastna(&super::super::super::lookup::reverse_complement(seq)),
            ];
            for (encoded, &offset) in contexts.iter().zip(&offsets[2 * q..2 * q + 2]) {
                let start = offset.max(0) as usize;
                let end = start.saturating_add(encoded.len());
                if start >= end || end > logical.len() {
                    continue;
                }
                logical[start..end].copy_from_slice(encoded);
                with_sentinels[start + 1..end + 1].copy_from_slice(encoded);
            }
        }
        (logical, with_sentinels)
    }

    // One buffer written context by context equals the two buffers of the IUPAC copies,
    // with an empty query (no room) between queries and letters of both cases.
    #[test]
    fn query_buffer_equals_the_buffers_of_the_context_copies() {
        let queries: [&[u8]; 5] = [
            b"ACGTTGCAAC",
            b"",
            b"acgtnRYKMSWBDHVNu-acgt",
            b"A",
            b"GGGGCCCCTTTTAAAAccgg*0",
        ];
        let mut offsets = Vec::new();
        let mut at = 0usize;
        let room = |length: usize| if length == 0 { 0 } else { length + 1 };
        for seq in &queries {
            offsets.push(at as i32);
            at += room(seq.len());
            offsets.push(at as i32);
            at += room(seq.len());
        }
        let concat_len = at - 1;
        let (logical, with_sentinels) =
            query_buffers_from_context_copies(&queries, &offsets, concat_len);
        let buffer = build_query_blastna_buffer(&queries[..], &offsets, concat_len);
        assert_eq!(buffer, with_sentinels);
        assert_eq!(&buffer[1..concat_len + 1], logical.as_slice());
        assert_eq!(
            (buffer[0], buffer[concat_len + 1]),
            (NUCL_SENTINEL, NUCL_SENTINEL)
        );
    }

    /// The contexts of `lengths` (query letters) searched on `strand`: offset and length of
    /// each of the two contexts of a query, laid out as NCBI's `s_QueryInfo_SetContext` does.
    fn strand_contexts(lengths: &[usize], strand: QueryStrand) -> Vec<QueryContext> {
        let (plus_on, minus_on) = strand.searched();
        let room = |length: usize| if length == 0 { 0 } else { length + 1 };
        let mut contexts = Vec::new();
        let mut at = 0usize;
        for (q, &length) in lengths.iter().enumerate() {
            for (frame, on) in [(1, plus_on), (-1, minus_on)] {
                let len = if on { length } else { 0 };
                contexts.push(QueryContext {
                    query_idx: q as u32,
                    frame,
                    query_offset: at as i32,
                    len,
                    masks: Vec::new(),
                });
                at += room(len);
            }
        }
        contexts
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:71-81,229-243
    // ```c
    //         Uint4 shift = prev_len ? prev_len + 1 : 0;
    //     ...
    //                     s_QueryInfo_SetContext(query_info, ctx_index+1, 0);
    // ```
    // The buffer of a search of one strand holds only the contexts of that strand, each
    // followed by a sentinel, and equals the one laid out naively; queries of 10, 3, 10
    // letters and an empty query.
    #[test]
    fn a_single_strand_query_buffer_holds_only_the_strand_it_searches() {
        let queries: [&[u8]; 5] = [b"ACGTTGCAAC", b"ACN", b"", b"GGRYCCTTAA", b"T"];
        for strand in [QueryStrand::Plus, QueryStrand::Minus] {
            let lengths: Vec<usize> = queries.iter().map(|q| q.len()).collect();
            let contexts = strand_contexts(&lengths, strand);
            let offsets: Vec<i32> = contexts.iter().map(|c| c.query_offset).collect();
            let last = contexts.last().unwrap();
            let concat_len =
                (last.query_offset as usize + last.len + if last.len == 0 { 1 } else { 2 })
                    .saturating_sub(2);
            // the naive buffer: a sentinel, then each searched context and a sentinel
            let mut naive = vec![NUCL_SENTINEL];
            for letters in &queries {
                if letters.is_empty() {
                    continue;
                }
                let encoded = if strand == QueryStrand::Plus {
                    encode_iupac_to_blastna(letters)
                } else {
                    encode_iupac_to_blastna(&super::super::super::lookup::reverse_complement(
                        letters,
                    ))
                };
                naive.extend(encoded);
                naive.push(NUCL_SENTINEL);
            }
            assert_eq!(naive.len(), concat_len + 2, "{strand:?}");
            let buffer =
                build_strand_query_blastna_buffer(&queries[..], &offsets, concat_len, strand);
            assert_eq!(buffer, naive, "{strand:?}");
        }
        // an empty query alone, and a single query, build a block of sentinels / its letters
        let buffer =
            build_strand_query_blastna_buffer(&[&b"AC"[..]][..], &[0, 3], 2, QueryStrand::Plus);
        assert_eq!(buffer, [NUCL_SENTINEL, 0, 1, NUCL_SENTINEL]);
        let buffer =
            build_strand_query_blastna_buffer(&[&b"AC"[..]][..], &[0, 0], 2, QueryStrand::Minus);
        assert_eq!(buffer, [NUCL_SENTINEL, 2, 3, NUCL_SENTINEL]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:219-243
    // ```c
    //     while (b < e - 1) {
    // 	m = (b + e) / 2;
    // 	if (A->contexts[m].query_offset > n)
    // 	    e = m;
    // 	else
    // 	    b = m;
    //     }
    //     return b;
    // ```
    // With one strand searched the context of every offset inside a valid context is that
    // context, through the index (whichever form the environment chooses), the table, the
    // runs, the constant-time form and `bsearch_context_info` (with `min_length` 0).
    #[test]
    fn a_single_strand_layout_maps_every_offset_to_its_context() {
        for lengths in [
            vec![10usize, 3, 10],
            vec![10, 3, 0, 10, 1],
            vec![0, 7, 7, 7],
            vec![5000, 17, 5000, 2],
            vec![1],
        ] {
            for strand in [QueryStrand::Plus, QueryStrand::Minus, QueryStrand::Both] {
                let contexts = strand_contexts(&lengths, strand);
                let index = QueryContextIndex::new(&contexts, strand != QueryStrand::Both);
                let total = contexts
                    .iter()
                    .map(|c| c.query_offset as usize + c.len)
                    .max()
                    .unwrap();
                let o1 = x_ctx_o1::CtxO1::build(&contexts, total);
                for (idx, ctx) in contexts.iter().enumerate() {
                    for n in ctx.query_offset as usize..ctx.query_offset as usize + ctx.len {
                        let want = (idx, format!("{lengths:?} {strand:?} offset {n}"));
                        assert_eq!(index.context_for_offset(n), want.0, "{}", want.1);
                        // (With both strands NCBI's bounds of `min_length` need every query to
                        // have letters; the one-strand layouts have `min_length` 0.)
                        if strand != QueryStrand::Both || lengths.iter().all(|&len| len > 0) {
                            assert_eq!(bsearch_context_info(n, &index), want.0, "{}", want.1);
                        }
                        if !index.direct_map.is_empty() {
                            assert_eq!(index.direct_map[n] as usize, want.0, "{}", want.1);
                        }
                        if index.x_total_len > n {
                            assert_eq!(index.x_direct(n), want.0, "{}", want.1);
                        }
                        if let Some(o1) = &o1 {
                            if let Some(context) = o1.get(n) {
                                assert_eq!(context, want.0, "{}", want.1);
                            }
                        }
                    }
                }
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:225-229
    // ```c
    // (query_info->contexts[query_info->last_context].query_offset +
    //  query_info->contexts[query_info->last_context].query_length) > kQueryLenForHashTable)
    //     p->container_type = eDiagHash;
    // ```
    // The diagonal hash starts above 8000 on the LAST context's end: with `-strand plus` that
    // is the empty minus context, one past the plus context, so a query of exactly 8000 letters
    // takes the hash; with `minus` the last context is the searched one; with both strands the
    // value is odd (two contexts and the gap between them).
    #[test]
    fn the_diagonal_hash_starts_at_the_end_of_the_last_context() {
        let hash = |lengths: &[usize], strand| {
            diag_container_query_length(&strand_contexts(lengths, strand)) > 8000
        };
        // plus
        assert!(!hash(&[7999], QueryStrand::Plus));
        assert!(hash(&[8000], QueryStrand::Plus));
        assert!(!hash(&[4000, 3998], QueryStrand::Plus));
        assert!(hash(&[4000, 3999], QueryStrand::Plus));
        // minus
        assert!(!hash(&[8000], QueryStrand::Minus));
        assert!(hash(&[8001], QueryStrand::Minus));
        assert!(!hash(&[4000, 3999], QueryStrand::Minus));
        assert!(hash(&[4000, 4000], QueryStrand::Minus));
        // both
        assert!(!hash(&[3999], QueryStrand::Both));
        assert!(hash(&[4000], QueryStrand::Both));
        // The value is the searched block's length (the query buffer without its sentinels)
        // plus 1 for `plus`, and that length otherwise.
        for (lengths, strand, extra) in [
            (vec![10usize, 3, 10], QueryStrand::Plus, 1),
            (vec![10, 3, 10], QueryStrand::Minus, 0),
            (vec![10, 3, 10], QueryStrand::Both, 0),
        ] {
            let contexts = strand_contexts(&lengths, strand);
            let last = contexts.last().unwrap();
            let concat_len =
                (last.query_offset as usize + last.len + if last.len == 0 { 1 } else { 2 })
                    .saturating_sub(2);
            assert_eq!(
                diag_container_query_length(&contexts),
                concat_len + extra,
                "{strand:?}"
            );
        }
        assert_eq!(diag_container_query_length(&[]), 0);
    }

    #[test]
    fn evalue_reap_keeps_unless_greater_than_cutoff() {
        assert!(hsp_survives_evalue_reap(10.0, 10.0));
        assert!(!hsp_survives_evalue_reap(1e-180, 0.0));
        assert!(hsp_survives_evalue_reap(0.0, 0.0));
        assert!(!hsp_survives_evalue_reap(10.000001, 10.0));
        assert!(hsp_survives_evalue_reap(f64::NAN, 10.0));
        assert!(hsp_survives_evalue_reap(10.0, f64::NAN));
        assert!(hsp_survives_evalue_reap(f64::INFINITY, f64::INFINITY));
        assert!(!hsp_survives_evalue_reap(f64::INFINITY, f64::MAX));
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
    fn make_equal_prelim_hit(query_idx: u32, seed_qs: usize, seed_ss: usize) -> PrelimHit {
        PrelimHit {
            context_idx: 0,
            query_idx,
            query_frame: 1,
            query_context_offset: 0,
            prelim_qs: 10,
            prelim_qe: 30,
            prelim_ss: 20,
            prelim_se: 40,
            prelim_score: 50,
            seed_qs,
            seed_ss,
            prelim_evalue: 0.0,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn prelim_hit(
        context_idx: u32,
        query_frame: i32,
        prelim_qs: usize,
        prelim_qe: usize,
        prelim_ss: usize,
        prelim_se: usize,
        prelim_score: i32,
    ) -> PrelimHit {
        PrelimHit {
            context_idx,
            query_idx: 0,
            query_frame,
            query_context_offset: 0,
            prelim_qs,
            prelim_qe,
            prelim_ss,
            prelim_se,
            prelim_score,
            seed_qs: prelim_qs,
            seed_ss: prelim_ss,
            prelim_evalue: 0.0,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2912-2945
    // ```c
    //    else {            /* query seq is split */
    //          if ((hsp1->query.frame >= 0 && hsp1->query.end >
    //                          split_offsets[offset_idx]) ||
    //              (hsp1->query.frame < 0 && hsp1->query.offset <
    //                          split_offsets[offset_idx] + chunk_overlap_size)) {
    //          if ((hsp2->query.frame < 0 && hsp2->query.end >
    //                          split_offsets[offset_idx]) ||
    //              (hsp2->query.frame >= 0 && hsp2->query.offset <
    //                          split_offsets[offset_idx] + chunk_overlap_size)) {
    // ```
    #[test]
    fn a_query_split_merges_hsps_across_the_split_on_either_strand() {
        let split = HspListSplit::Query {
            offsets: [1000, 5000],
        };
        // Plus strand: the earlier chunk's HSP ends after the split point, the later
        // chunk's starts before its overlap ends.
        let mut combined = vec![prelim_hit(0, 1, 900, 1050, 100, 250, 150)];
        let incoming = vec![prelim_hit(0, 1, 1000, 1200, 200, 400, 200)];
        merge_prelim_hit_lists(&mut combined, incoming, split, 100, true);
        assert_eq!(combined.len(), 1);
        assert_eq!((combined[0].prelim_qs, combined[0].prelim_qe), (900, 1200));
        // The score density of the two (350 / 350) over the merged 300 residues.
        assert_eq!(combined[0].prelim_score, 300);
        // Minus strand: the later chunk comes first on the strand, so the tests swap.
        let mut combined = vec![prelim_hit(1, -1, 4950, 5200, 1000, 1250, 150)];
        let incoming = vec![prelim_hit(1, -1, 4800, 5050, 850, 1100, 150)];
        merge_prelim_hit_lists(&mut combined, incoming, split, 100, true);
        assert_eq!(combined.len(), 1);
        assert_eq!((combined[0].prelim_qs, combined[0].prelim_qe), (4800, 5200));
        // The same pair on the plus strand is not in the strips (the earlier HSP does not
        // end after the split point): both stay.
        let mut combined = vec![prelim_hit(0, 1, 800, 950, 1000, 1150, 150)];
        let incoming = vec![prelim_hit(0, 1, 700, 850, 900, 1050, 150)];
        merge_prelim_hit_lists(&mut combined, incoming, split, 100, true);
        assert_eq!(combined.len(), 2);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
    // ```c
    // if (0 == (result = BLAST_CMP(hsp2->score, hsp1->score)) && ...) {
    //     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
    // }
    // return result;
    // ```
    #[test]
    fn stable_prelim_sort_preserves_equal_complete_record_order() {
        let mut hsps = vec![
            make_equal_prelim_hit(7, 101, 201),
            make_equal_prelim_hit(3, 102, 202),
            make_equal_prelim_hit(9, 103, 203),
        ];

        qsort_prelim_hits_by(&mut hsps, score_compare_prelim_hits);

        assert_eq!(
            hsps.iter().map(|hsp| hsp.query_idx).collect::<Vec<_>>(),
            vec![7, 3, 9]
        );
        assert_eq!(
            hsps.iter()
                .map(|hsp| (hsp.seed_qs, hsp.seed_ss))
                .collect::<Vec<_>>(),
            vec![(101, 201), (102, 202), (103, 203)]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2532
    // ```c
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
    // ...
    // hsp = Blast_HSPFree(hsp);
    // ```
    #[test]
    fn prelim_purge_retains_the_first_complete_equal_record() {
        let first = make_equal_prelim_hit(7, 101, 201);
        let second = make_equal_prelim_hit(3, 102, 202);

        let kept = purge_prelim_hits_with_common_endpoints(vec![first, second]);

        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].query_idx, 7);
        assert_eq!((kept[0].seed_qs, kept[0].seed_ss), (101, 201));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:689-706
    // ```c
    // p.first = (slp->GetInt().GetFrom() > offset)? slp->GetInt().GetFrom() - offset : 0;
    // p.second = MIN(slp->GetInt().GetTo() - offset, length-1);
    // if (slp->GetInt().GetTo() >= offset && p.first < length) {
    //     output.push_back(p);
    // }
    // ```
    // Each range starts at the last letter of the previous mask, as NCBI's ranges read
    // from `TSequenceRanges::get_data` (seqdb.hpp:280-282).
    #[test]
    fn test_build_subject_seq_ranges_from_masks() {
        let masks = vec![MaskedInterval::new(2, 5), MaskedInterval::new(8, 12)];
        let ranges = build_subject_seq_ranges_from_masks(&masks, 10);
        assert_eq!(ranges, vec![(0, 2), (4, 8), (9, 10)]);
        let masks = vec![MaskedInterval::new(0, 5)];
        assert_eq!(
            build_subject_seq_ranges_from_masks(&masks, 32),
            vec![(0, 0), (4, 32)]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
    // ```c
    // if (hsp->query.frame != hsp->subject.frame) {
    //    *q_end = query_length - hsp->query.offset;
    //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
    //    *s_end = hsp->subject.offset + 1;
    //    *s_start = hsp->subject.end;
    // }
    // ```
    #[test]
    fn test_adjust_blastn_offsets_minus_query_keeps_internal_subject_order() {
        let adjusted = adjust_blastn_offsets(10, 20, 30, 40, 100, -1);
        assert_eq!(adjusted, (81, 90, 40, 31));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1887-1890
    // ```c
    // hsp->evalue =
    //     BLAST_KarlinStoE_simple(score, kbp[kbp_context],
    //                      query_info->contexts[hsp->context].eff_searchsp);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1923-1926
    // ```c
    // hsp->bit_score =
    //    (hsp->score*kbp[hsp->context]->Lambda - kbp[hsp->context]->logK) /
    //    NCBIMATH_LN2;
    // ```
    #[test]
    fn test_blastn_context_statistics_uses_supplied_eff_searchsp() {
        let params = crate::stats::KarlinParams {
            lambda: 0.625,
            k: 0.41,
            h: 0.78,
            alpha: 1.0,
            beta: -1.0,
        };
        let eff_searchsp = 44_573_947_800_i64;
        let raw_score = 46;

        let (bit_score, evalue) =
            calculate_blastn_context_statistics(raw_score, &params, eff_searchsp, false);

        let expected_bit_score = (params.lambda * raw_score as f64 - params.k.ln()) / NCBIMATH_LN2;
        let expected_evalue =
            (eff_searchsp as f64) * (-(params.lambda) * raw_score as f64 + params.k.ln()).exp();

        assert_eq!(bit_score, expected_bit_score);
        assert_eq!(evalue, expected_evalue);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1864-1870
    // ```c
    // score = hsp->score;
    // if (hsp_list && hsp_list->hspcnt != 0
    //         && gapped_calculation && sbp->round_down) {
    //     score &= ~1;
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1923-1926
    // ```c
    // hsp->bit_score =
    //    (hsp->score*kbp[hsp->context]->Lambda - kbp[hsp->context]->logK) /
    //    NCBIMATH_LN2;
    // ```
    #[test]
    fn test_blastn_context_statistics_rounds_evalue_score_only() {
        let params = crate::stats::KarlinParams {
            lambda: 0.625,
            k: 0.41,
            h: 0.78,
            alpha: 0.8,
            beta: -2.0,
        };
        let eff_searchsp = 86_850_803_746_i64;
        let raw_score = 45;

        let (bit_score, evalue) =
            calculate_blastn_context_statistics(raw_score, &params, eff_searchsp, true);

        let expected_bit_score = (params.lambda * raw_score as f64 - params.k.ln()) / NCBIMATH_LN2;
        let expected_evalue_score = raw_score & !1;
        let expected_evalue = (eff_searchsp as f64)
            * (-(params.lambda) * expected_evalue_score as f64 + params.k.ln()).exp();

        assert_eq!(bit_score, expected_bit_score);
        assert_eq!(evalue, expected_evalue);
    }
    // NCBI reference: c++/src/app/blast/blast_app_util.cpp:204-210;
    // c++/src/algo/blast/core/blast_parameters.c:925-946
    // db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, true));
    // searchsp = query_info->contexts[context].eff_searchsp;
    // BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
    // NCBI BLAST+ 2.17.0+ raw output is frozen for three unequal subjects.
    // Their combined search space excludes weak word7 seeds which could grow
    // into extra gapped HSPs under an incorrectly lowered per-subject cutoff.
    #[test]
    fn word7_uses_subject_set_cutoffs() {
        use clap::Parser;
        #[derive(Parser)]
        struct Options {
            #[command(flatten)]
            blastn: BlastnArgs,
        }
        let sequence: String = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fasta/small_test.fasta"
        ))
        .lines()
        .filter(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect();
        let query: String = (0..3)
            .map(|index| format!(">seq{index}\n{}\n", &sequence[..900]))
            .collect();
        let subject: String = [899, 900, 903]
            .into_iter()
            .enumerate()
            .map(|(index, length)| format!(">seq{index}\n{}\n", &sequence[..length]))
            .collect();
        let expected = include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/unit/helpers/ncbi_reference_data/blastn_word7_subject_set.out"
        ));
        let threads: &[usize] = if cfg!(feature = "parallel") {
            &[1, 4]
        } else {
            &[1]
        };
        for n in threads {
            let options = crate::cli::try_parse_from::<Options, _, _>([
                "test",
                "-query",
                "memory-query",
                "-subject",
                "memory-subject",
                "-task",
                "blastn",
                "-word_size",
                "7",
                "-outfmt",
                "6",
                "-num_threads",
                &n.to_string(),
            ])
            .unwrap();
            // NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
            // BlastSeqBlkSetSequence(subj, sequence.data.release(), ...);
            let mut output = Vec::new();
            let mut stderr = std::io::stderr();
            run_local(
                options.blastn,
                &reader_records(&query, ReaderConfig::query("BLASTN", false, false)),
                &reader_records(&subject, ReaderConfig::subject("BLASTN", false, false)),
                &mut ReportOutputs::single("6", OutputSink::Writer(&mut output), &mut stderr),
            )
            .unwrap();
            assert_eq!(output, expected.as_slice(), "n{n}");
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-406,509-513,583-612
    // ```c
    // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info,
    //                                 hit_options->min_diag_separation)) {
    //   BLAST_GappedAlignmentWithTraceback(...);
    //   Blast_HSPUpdateWithTraceback(gap_align, hsp);
    //   status = BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
    // }
    // ```
    // Repeated, unequal-length gapped matches generate dependent prelim HSPs
    // across batch boundaries. Raw tabular bytes guard accepted HSP order and
    // endpoints; the DP check below compares scores and complete edit scripts.
    #[cfg(feature = "parallel")]
    #[test]
    fn speculative_traceback_preserves_ordered_gapped_hsps() {
        use clap::Parser;
        #[derive(Parser)]
        struct Options {
            #[command(flatten)]
            blastn: BlastnArgs,
        }
        let mut state = 0x12345678_u32;
        let query: Vec<u8> = (0..512)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        // NCBI reference: c++/src/algo/blast/core/blast_gapalign.c:4170-4189
        // if (subject_length < MAX_SUBJECT_OFFSET) { *start_shift = 0; return; }
        // max_extension_left = query_offset + MAX_TOTAL_GAPS;
        // *start_shift = s_offset - max_extension_left;
        // Exercise a real subject shift and more than one 16-HSP batch.
        let mut subject = vec![b'N'; 91_000];
        for i in 0..20 {
            subject.extend_from_slice(&query[..240]);
            subject.extend_from_slice(b"GATTACA");
            subject.extend_from_slice(&query[240..512 - i * 13]);
            subject.extend_from_slice(&[b'N'; 64]);
        }
        let directory = std::env::temp_dir().join(format!(
            "losat-traceback-regression-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&directory).unwrap();
        let query_path = directory.join("query.fa");
        let subject_path = directory.join("subject.fa");
        let output_path = directory.join("output.txt");
        std::fs::write(
            &query_path,
            format!(">query\n{}\n", String::from_utf8(query).unwrap()),
        )
        .unwrap();
        std::fs::write(
            &subject_path,
            format!(">subject\n{}\n", String::from_utf8(subject).unwrap()),
        )
        .unwrap();
        let execute = |threads: usize, outfmt: &str| {
            let args = crate::cli::try_parse_from::<Options, _, _>([
                "test",
                "-query",
                query_path.to_str().unwrap(),
                "-subject",
                subject_path.to_str().unwrap(),
                "-task",
                "blastn",
                "-word_size",
                "11",
                "-num_threads",
                &threads.to_string(),
                "-outfmt",
                outfmt,
                "-out",
                output_path.to_str().unwrap(),
            ])
            .unwrap()
            .blastn;
            run(
                args,
                crate::blastinput::ncbi_environment::ApplicationSettings::default(),
            )
            .unwrap();
            std::fs::read(&output_path).unwrap()
        };
        let sequential = execute(1, "6");
        let text = std::str::from_utf8(&sequential).unwrap();
        assert!(
            text.lines().count() >= 20,
            "must cross the sixteen-HSP batch boundary"
        );
        assert!(text
            .lines()
            .any(|line| line.split('\t').nth(5).unwrap() != "0"));
        for threads in [4, 2, 4] {
            assert_eq!(execute(threads, "6"), sequential, "n{threads}");
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:509-513
        // ```c
        // BLAST_GappedAlignmentWithTraceback(program_number, query,
        //     adjusted_subject, gap_align, score_params, q_start, s_start,
        //     query_length, adjusted_s_length, fence_hit);
        // ```
        // Multiple seeds extend to the same gapped HSP. Compute in reverse input
        // order with independent scratch, then replay in the original seed order.
        let q = encode_iupac_to_blastna(
            &reader_records(
                &std::fs::read_to_string(&query_path).unwrap(),
                ReaderConfig::query("BLASTN", false, false),
            )[0]
            .seq(),
        );
        let subj = encode_iupac_to_blastna(
            &reader_records(
                &std::fs::read_to_string(&subject_path).unwrap(),
                ReaderConfig::subject("BLASTN", false, false),
            )[0]
            .seq(),
        );
        let matrix = build_blastna_matrix(2, -3);
        let seeds = [60, 200, 300, 400, 80, 180, 320, 420];
        let dp = |index: usize, scratch: &mut GapAlignScratch| {
            let qs = seeds[index];
            extend_gapped_heuristic_with_traceback_with_scratch(
                &q,
                &subj,
                qs,
                91_000 + qs + if qs >= 240 { 7 } else { 0 },
                1,
                2,
                -3,
                &matrix,
                5,
                2,
                60,
                scratch,
            )
        };
        let mut scratch = GapAlignScratch::new();
        let ordered: Vec<_> = (0..8).map(|i| dp(i, &mut scratch)).collect();
        // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
        // (*thread)->Run(); (*thread)->Join(&result);
        let reversed: Vec<_> = crate::utils::threading::with_search_pool(4, "blastn", |pool| {
            Ok(pool.install(|| {
                (0..8)
                    .into_par_iter()
                    .rev()
                    .map_init(GapAlignScratch::new, |scratch, i| dp(i, scratch))
                    .collect()
            }))
        })
        .unwrap();
        assert_eq!(ordered, reversed.into_iter().rev().collect::<Vec<_>>());
        assert!(ordered.iter().all(|hsp| hsp.7 > 0 && !hsp.9.is_empty()));
        std::fs::remove_dir_all(directory).unwrap();
    }

    fn query_prelim_hit(query_idx: u32, prelim_score: i32) -> PrelimHit {
        PrelimHit {
            query_idx,
            ..prelim_hit(query_idx * 2, 1, 10, 30, 20, 40, prelim_score)
        }
    }

    fn prelim_scores(hit_list: &HitList<PrelimHspList>) -> Vec<(u32, Vec<i32>)> {
        hit_list
            .hsplist_array
            .iter()
            .map(|list| {
                (
                    list.oid,
                    list.hsps.iter().map(|hsp| hsp.prelim_score).collect(),
                )
            })
            .collect()
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_collector.c:116-161
    // ```c
    //       for (index = 0; index < hsp_list->hspcnt; index++) {
    //          query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
    //          Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
    //       ...
    //             if (!results->hitlist_array[index]) {
    //                results->hitlist_array[index] =
    //                   Blast_HitListNew(params->prelim_hitlist_size);
    //             }
    //             Blast_HitListUpdate(results->hitlist_array[index],
    //                                 hsp_list_array[index]);
    // ```
    #[test]
    fn the_collector_keeps_prelim_hitlist_size_subjects_per_query() {
        let subject_hits = vec![
            vec![query_prelim_hit(0, 30), query_prelim_hit(1, 90)],
            vec![query_prelim_hit(0, 50)],
            vec![query_prelim_hit(0, 40), query_prelim_hit(0, 45)],
            Vec::new(),
        ];
        let mut hit_lists = collect_prelim_hit_lists(subject_hits, 3, 2);
        assert!(hit_lists[2].is_none());
        let mut query0 = hit_lists[0].take().unwrap();
        // Subject 2 replaces subject 0 (the worst first score); a list that enters the
        // heap is sorted by e-value, then score (`Blast_HSPListSortByEvalue`).
        query0.sort_by_evalue();
        assert_eq!(
            prelim_scores(&query0),
            vec![(1, vec![50]), (2, vec![45, 40])]
        );
        let mut query1 = hit_lists[1].take().unwrap();
        query1.sort_by_evalue();
        assert_eq!(prelim_scores(&query1), vec![(0, vec![90])]);
    }

    // `LOSAT_X_BNHITGROUP`: the lists made from a subject's HSPs alone equal the lists made with
    // one list per query of the batch, in every hit list, with the HSPs of a list in their order
    // (each HSP has its own score). Queries in random order, in runs and all in one query.
    #[test]
    fn grouped_hit_lists_equal_the_lists_per_query() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = move |below: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % below
        };
        for round in 0..200u64 {
            let num_queries = 1 + next(40) as usize;
            let size = 1 + next(5) as usize;
            let mut per_query = PrelimHitCollector::new(num_queries, size);
            let mut grouped = PrelimHitCollector::new(num_queries, size);
            let mut score = 0;
            for oid in 0..30 {
                let count = next(25);
                let run = 1 + next(4);
                let mut hits = Vec::new();
                for i in 0..count {
                    let query = match round % 3 {
                        0 => next(num_queries as u64),
                        1 => (i / run).min(num_queries as u64 - 1),
                        _ => num_queries as u64 - 1,
                    };
                    score += 1;
                    hits.push(query_prelim_hit(query as u32, score));
                }
                if hits.is_empty() {
                    continue;
                }
                per_query.add_per_query(oid, hits.clone());
                grouped.x_add_grouped(oid, hits);
            }
            for (a, b) in per_query.hit_lists.iter().zip(&grouped.hit_lists) {
                assert_eq!(a.as_ref().map(prelim_scores), b.as_ref().map(prelim_scores));
                assert_eq!(a.is_some(), b.is_some());
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2163-2210
    // ```c
    //     while (i < num_hsplists1 && j < num_hsplists2) {
    //         if (hsplist1->oid < hsplist2->oid) {
    //             Blast_HitListUpdate(new_hitlist, hsplist1);
    //         ...
    //             else {
    //                 Blast_HSPListAppend(hitlist1->hsplist_array + i,
    //                                     hitlist2->hsplist_array + j,
    //                                     hsplist2->hsp_max);
    //             }
    //             Blast_HitListUpdate(new_hitlist, hitlist2->hsplist_array[j]);
    // ```
    #[test]
    fn merging_prelim_hit_lists_appends_shared_subjects_and_keeps_the_size() {
        let list = |oid: u32, scores: &[i32]| PrelimHspList {
            oid,
            hsps: scores
                .iter()
                .map(|&score| query_prelim_hit(0, score))
                .collect(),
            best_evalue: 0.0,
        };
        let mut hitlist1 = HitList::new(3);
        hitlist1.update(list(4, &[50]));
        hitlist1.update(list(2, &[70]));
        let mut combined = HitList::new(3);
        combined.update(list(9, &[80]));
        combined.update(list(2, &[60]));
        combined.update(list(7, &[40]));
        let mut combined = Some(combined);
        merge_prelim_hit_list(hitlist1, &mut combined, [0, 0], 100);
        let mut merged = combined.unwrap();
        assert_eq!(merged.hsplist_max, 3);
        merged.sort_by_evalue();
        assert_eq!(
            prelim_scores(&merged),
            vec![(9, vec![80]), (2, vec![70, 60]), (4, vec![50])]
        );
        let mut empty = None;
        let mut single = HitList::new(3);
        single.update(list(5, &[10]));
        merge_prelim_hit_list(single, &mut empty, [0, 0], 100);
        assert_eq!(prelim_scores(empty.as_ref().unwrap()), vec![(5, vec![10])]);
    }

    /// The NCBI errors of BLASTN's `-outfmt` at NCBI's positions (blast_args.cpp:2745-2748,
    /// 2874-2886, 2975-2977; blastn_app.cpp:255-258).
    #[test]
    fn output_formats_raise_ncbis_errors_where_ncbi_does() {
        use clap::Parser;
        #[derive(Parser)]
        struct Options {
            #[command(flatten)]
            blastn: BlastnArgs,
        }
        let args = |words: &[&str]| {
            let mut argv = vec!["test", "-query", "q", "-subject", "s"];
            argv.extend_from_slice(words);
            crate::cli::try_parse_from::<Options, _, _>(argv)
                .unwrap()
                .blastn
        };
        let native = |error: anyhow::Error| {
            let native = error.downcast::<crate::cli::NativeError>().unwrap();
            (native.exit, native.message)
        };
        let options = |text: &str| {
            (
                1,
                format!(
                    "BLAST query/options error: {text}\nPlease refer to the BLAST+ user manual.\n"
                ),
            )
        };
        // The CLI: a format that NCBI rejects later has no report (`cli_output_formats`).
        let file = Some(std::path::Path::new("o.json"));
        for outfmt in ["13", "14", "19", "21"] {
            assert!(
                cli_output_formats(outfmt, None).unwrap().is_empty(),
                "{outfmt}"
            );
        }
        assert!(cli_output_formats("13", Some(std::path::Path::new("-")))
            .unwrap()
            .is_empty());
        // 13 and 14 to a file are written (OF5c part 1): the XML2 / JSON documents as files.
        for (outfmt, kind) in [
            ("13", crate::report::blast_xml2::Xml2Kind::Files13),
            ("14", crate::report::blast_xml2::Xml2Kind::Files14),
        ] {
            assert_eq!(
                cli_output_formats(outfmt, file).unwrap(),
                [BlastnOutputFormat::Xml2(kind)],
                "{outfmt}"
            );
        }
        for outfmt in ["1", "18"] {
            assert_eq!(
                cli_output_formats(outfmt, file).unwrap_err().to_string(),
                format!("output format {outfmt} is not supported by LOSAT's BLASTN")
            );
        }
        // SAM (OF6) reads its words and ignores the rest and a delimiter.
        assert_eq!(
            cli_output_formats("17", file).unwrap(),
            [BlastnOutputFormat::Sam(
                crate::report::sam::SamSpec::default()
            )]
        );
        assert_eq!(
            cli_output_formats(" 17 delim=, SR SQ foo", None).unwrap(),
            [BlastnOutputFormat::Sam(crate::report::sam::SamSpec {
                subject_as_reference: true,
                sequence: true
            })]
        );
        assert_eq!(
            cli_output_formats("0 qaccver", None).unwrap(),
            [BlastnOutputFormat::Pairwise]
        );
        // A field list, a delimiter, and outfmt 10 and 20 are the custom tabular report (line CLI
        // C2); the default 6 and 7 (`std` with a tab) keep their own writers.
        let custom = |outfmt: &str| match cli_output_formats(outfmt, None).unwrap().as_slice() {
            [BlastnOutputFormat::CustomTabular(spec)] => (*spec).clone(),
            other => panic!("{outfmt}: {other:?}"),
        };
        assert_eq!(
            custom("7 qaccver").fields,
            [crate::report::tabular::TabularField::QueryAccessionVersion]
        );
        assert_eq!(custom("6 delim=, std").delimiter, ",");
        assert_eq!(
            custom("20").kind,
            crate::report::tabular::TabularKind::CsvWithHeader
        );
        assert_eq!(
            cli_output_formats("6 std", None).unwrap(),
            [BlastnOutputFormat::Tabular]
        );
        // A tab given as `delim=` keeps the default writer; NCBI takes `delim=\t` as written, a
        // backslash and a `t`.
        assert_eq!(
            cli_output_formats("7 delim=\t std", None).unwrap(),
            [BlastnOutputFormat::TabularWithComments]
        );
        assert_eq!(custom(r"7 delim=\t").delimiter, r"\t");
        assert_eq!(
            native(cli_output_formats("22", None).unwrap_err()),
            (
                255,
                "Error: Formatting choice is out of range\n".to_string()
            )
        );
        // The formatting handler raises its errors before its warning; 13 passes it.
        for (outfmt, error) in [
            (
                "19",
                Some("AIRR rearrangement format is only applicable to igblastn"),
            ),
            (
                "21",
                Some("FASTA output format is only applicable to magicblast"),
            ),
            ("13", None),
        ] {
            let mut diagnostics = Vec::new();
            let mut sink = Vec::new();
            let mut outputs =
                ReportOutputs::single(outfmt, OutputSink::Writer(&mut sink), &mut diagnostics);
            let result = process_options(&mut args(&["-max_target_seqs", "1"]), &mut outputs);
            drop(outputs);
            match error {
                Some(error) => {
                    assert_eq!(native(result.unwrap_err()), options(error), "{outfmt}");
                    assert!(diagnostics.is_empty(), "{outfmt}");
                }
                None => {
                    result.unwrap();
                    assert_eq!(
                        String::from_utf8(diagnostics).unwrap(),
                        "Warning: [blastn] Examining 5 or more matches is recommended\n"
                    );
                }
            }
        }
        // `run_local` has no output file; its formats that NCBI rejects later fail at once.
        for (outfmt, error) in [
            ("13", "Please provide a file name for outfmt 13."),
            ("14", "Please provide a file name for outfmt 14."),
            (
                "19",
                "AIRR rearrangement format is only applicable to igblastn",
            ),
            ("21", "FASTA output format is only applicable to magicblast"),
        ] {
            assert_eq!(
                native(parse_output_formats(["6", outfmt]).unwrap_err()),
                options(error),
                "{outfmt}"
            );
        }
        assert_eq!(
            parse_output_formats(["0", "6", "07"]).unwrap(),
            [
                BlastnOutputFormat::Pairwise,
                BlastnOutputFormat::Tabular,
                BlastnOutputFormat::TabularWithComments
            ]
        );
        assert_eq!(
            parse_output_formats(["17", "17 SQ SR"]).unwrap(),
            [
                BlastnOutputFormat::Sam(crate::report::sam::SamSpec::default()),
                BlastnOutputFormat::Sam(crate::report::sam::SamSpec {
                    subject_as_reference: true,
                    sequence: true
                })
            ]
        );
    }
}
