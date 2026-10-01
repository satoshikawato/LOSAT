//! Main BLASTN run function
//!
//! This module contains the main `run()` function that coordinates the BLASTN
//! search process.

use crate::api::local_blast::{FormatProbe, OutputSink, ReportOutputs};
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
    encode_iupac_to_blastna, encode_subject_ncbi2na_packed, COMPRESSION_RATIO,
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
    greedy_gapped_alignment_score_only,
    greedy_gapped_alignment_with_traceback,
    stats_from_edit_ops,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    GapAlignScratch,
    GreedyAlignScratch,
};
use super::super::args::BlastnArgs;
use super::super::blast_extend::DiagStruct;
use super::super::constants::TWO_HIT_WINDOW;
use super::super::coordination::{
    build_lookup_tables, chunk_query_masks, collect_lowercase_masks, configure_task,
    finalize_task_config, prepare_sequence_data, query_masks, scan_subjects_metadata,
    subject_metadata_from_records, LookupTables, SubjectMetadata,
};
use super::super::extension::{
    build_compressed_query, build_nucl_score_table, build_query_four_base_bytes,
    extend_hit_ungapped_approx_ncbi, extend_hit_ungapped_exact_ncbi, type_of_word, SmallNaWord,
};
use super::super::filtering::{
    blast_hsp_test_identity_and_length, hsp_test, purge_hsps_with_common_endpoints,
    purge_hsps_with_common_endpoints_ex, reevaluate_hsp_with_ambiguities_gapped_ex,
    subject_best_hit, subject_best_hit_by, BestHitKey, ReevalParams,
};
use super::super::hsp::{
    evalue_comp, get_prelim_hitlist_size, parse_blastn_output_format, sort_hsplist_by_evalue,
    sort_hsps_by_score, trim_by_max_hsps, write_output_blastn_hitlists_to_writer, BlastnHitList,
    BlastnHsp, BlastnHspList, BlastnOutputFormat, HitList, HitListEntry, NCBI_BLASTN_VERSION,
};
use super::super::input::{
    check_deflines, check_records, check_records_have_residues, check_residues,
    check_sequence_lines, is_blank, with_u_as_t, write_title_warnings, UNREADABLE_FASTA,
};
use super::super::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use super::super::lookup::{build_unmasked_ranges, reverse_complement};
use super::super::pairwise::{pairwise_hits, DisplayMasks};
use super::super::query_split::{
    query_chunk_size, restrict_masks, split_query_batch, QueryChunk, QUERY_CHUNK_OVERLAP,
};
use super::super::scoring::{
    check_greedy_gap_costs, check_losat_limits, check_scoring_options, context_blocks,
    context_ungapped_blocks, gap_x_dropoffs, karlin_error, ContextKarlin,
};
use crate::blastinput::query_batch::{next_query_batch_end, BatchSizeMixer};
use crate::report::pairwise::{
    write_blastn_pairwise_prolog, write_blastn_pairwise_report, BlastnPairwiseQuery,
    BlastnPairwiseReport, PairwiseConfig,
};
use crate::report::query_warnings::{few_matches_warning, invalid_query_warning};
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

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/ncbi_math.h:160-161
// ```c
// #define NCBIMATH_LN2 0.69314718055994530941723212145818
// ```
const NCBIMATH_LN2: f64 = 0.69314718055994530941723212145818;

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showdefline.cpp:69-79
// ```c
// string CShowBlastDefline::GetDefline(const CBioseq_Handle& bioseq) const
// {
//     return bioseq.GetSeqId()->AsFastaString();
// }
// ```
#[inline]
fn fasta_defline(record: &bio::io::fasta::Record) -> Arc<str> {
    match record.desc() {
        Some(desc) => Arc::<str>::from(format!("{} {}", record.id(), desc)),
        None => Arc::<str>::from(record.id()),
    }
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
    hard_ranges: Vec<(i32, i32)>,
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
            hard_ranges: vec![(0, full_right)],
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
        Self {
            num_buckets: DIAGHASH_NUM_BUCKETS,
            occupancy: 1,
            capacity,
            backbone: vec![0; DIAGHASH_NUM_BUCKETS],
            chain: vec![DiagHashCell::default(); capacity as usize],
            offset: window,
            window,
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
    fn insert(
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

/// Structure to hold ungapped hit data for batch processing
/// NCBI reference: blast_gapalign.c - init_hsp_array is sorted by score descending
/// before gapped extension
#[derive(Clone)]
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
#[derive(Clone)]
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
    let mut kept: Vec<PrelimHit> = Vec::with_capacity(hsps.len());
    for hsp in hsps {
        let same_common_start = kept.last().is_some_and(|prev| {
            prev.context_idx == hsp.context_idx
                && prev.prelim_qs == hsp.prelim_qs
                && prev.prelim_ss == hsp.prelim_ss
        });
        if !same_common_start {
            kept.push(hsp);
        }
    }

    if kept.len() <= 1 {
        return kept;
    }

    qsort_prelim_hits_by(&mut kept, query_end_compare_prelim_hits);
    let mut out: Vec<PrelimHit> = Vec::with_capacity(kept.len());
    for hsp in kept {
        let same_common_end = out.last().is_some_and(|prev| {
            prev.context_idx == hsp.context_idx
                && prev.prelim_qe == hsp.prelim_qe
                && prev.prelim_se == hsp.prelim_se
        });
        if !same_common_end {
            out.push(hsp);
        }
    }

    out
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
    let mut hit_lists: PrelimHitLists = Vec::with_capacity(num_queries);
    hit_lists.resize_with(num_queries, || None);
    for (oid, hits) in subject_hits.into_iter().enumerate() {
        if hits.is_empty() {
            continue;
        }
        let mut per_query: Vec<Vec<PrelimHit>> = vec![Vec::new(); num_queries];
        for hit in hits {
            per_query[hit.query_idx as usize].push(hit);
        }
        for (query, hsps) in per_query.into_iter().enumerate() {
            if hsps.is_empty() {
                continue;
            }
            hit_lists[query]
                .get_or_insert_with(|| HitList::new(prelim_hitlist_size))
                .update(PrelimHspList {
                    oid: oid as u32,
                    hsps,
                    best_evalue: 0.0,
                });
        }
    }
    hit_lists
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
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
// ```c
// hsp_list = Blast_HSPListFree(hsp_list);
// BlastInitHitListReset(init_hitlist);
// ```
impl SubjectScratch {
    fn new(query_count: usize, offset_array_size: usize) -> Self {
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
            diag_hash: DiagHashTable::new(TWO_HIT_WINDOW as i32),
            diag_table_offset: TWO_HIT_WINDOW as i32,
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
        }
    }
}

#[derive(Clone)]
struct QueryContext {
    query_idx: u32,
    frame: i32,
    query_offset: i32,
    seq: Vec<u8>,
    masks: Vec<MaskedInterval>,
}

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
}

impl QueryContextIndex {
    fn new(contexts: &[QueryContext]) -> Self {
        let offsets: Vec<usize> = contexts
            .iter()
            .map(|ctx| ctx.query_offset.max(0) as usize)
            .collect();
        let min_length = contexts
            .iter()
            .map(|ctx| ctx.seq.len())
            .filter(|&len| len > 0)
            .min()
            .unwrap_or(0);
        let max_length = contexts.iter().map(|ctx| ctx.seq.len()).max().unwrap_or(0);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:219-238 (BSearchContextInfo)
        // ```c
        // if (A->contexts[m].query_offset > n) e = m;
        // else b = m;
        // ...
        // return b;
        // ```
        let total_len = contexts
            .iter()
            .map(|ctx| ctx.query_offset.max(0) as usize + ctx.seq.len())
            .max()
            .unwrap_or(0);
        let mut direct_map = vec![0u32; total_len];
        for (idx, ctx) in contexts.iter().enumerate() {
            let start = ctx.query_offset.max(0) as usize;
            let end = start.saturating_add(ctx.seq.len());
            if start < end && end <= direct_map.len() {
                direct_map[start..end].fill(idx as u32);
            }
        }
        Self {
            offsets,
            min_length,
            max_length,
            direct_map,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:219-238 (BSearchContextInfo)
    // ```c
    // if (A->contexts[m].query_offset > n) e = m;
    // else b = m;
    // ...
    // return b;
    // ```
    fn context_for_offset(&self, n: usize) -> usize {
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
fn build_query_blastna_concat_buffers(
    encoded_contexts: &[Vec<u8>],
    context_offsets: &[i32],
    query_concat_length: usize,
) -> (Vec<u8>, Vec<u8>) {
    debug_assert_eq!(encoded_contexts.len(), context_offsets.len());

    let mut logical = vec![NUCL_SENTINEL; query_concat_length];
    let mut with_sentinels = vec![NUCL_SENTINEL; query_concat_length + 2];

    for (encoded, &offset) in encoded_contexts.iter().zip(context_offsets.iter()) {
        let start = offset.max(0) as usize;
        let end = start.saturating_add(encoded.len());
        if start >= end || end > logical.len() {
            continue;
        }
        logical[start..end].copy_from_slice(encoded);
        with_sentinels[start + 1..end + 1].copy_from_slice(encoded);
    }

    (logical, with_sentinels)
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
        if ctx.seq.is_empty() {
            continue;
        }
        let ctx_offset = ctx.query_offset.max(0) as usize;
        let ranges = build_unmasked_ranges(ctx.seq.len(), masks);
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
            max_q_off = max_q_off.max(ctx_offset + ctx.seq.len() - 1);
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1173-1178
// ```c
// masks->ssr->left  = query_length - 1 - masks->ssr->right;
// masks->ssr->right = query_length - 1 - masks->ssr->left;
// ```
fn reverse_mask_intervals(masks: &[MaskedInterval], query_length: usize) -> Vec<MaskedInterval> {
    let mut out: Vec<MaskedInterval> = masks
        .iter()
        .map(|m| {
            let start = query_length.saturating_sub(m.end);
            let end = query_length.saturating_sub(m.start);
            MaskedInterval { start, end }
        })
        .collect();
    out.sort_by_key(|m| m.start);
    out
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

fn print_blastn_timing(
    timing: &BlastnTiming,
    t_search_start: Option<std::time::Instant>,
    t_total: Option<std::time::Instant>,
) {
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

    eprintln!("[TIMING] scan_lookup: {:.3}s (calls={})", scan_s, scan_n);
    eprintln!(
        "[TIMING] ungapped_extend: {:.3}s (calls={})",
        ungapped_s, ungapped_n
    );
    eprintln!(
        "[TIMING] gapped_extend: {:.3}s (calls={})",
        gapped_s, gapped_n
    );
    eprintln!(
        "[TIMING] traceback_prune: {:.3}s (calls={})",
        traceback_s, traceback_n
    );
    eprintln!(
        "[TIMING] traceback_sort_prelim: {:.3}s",
        timing_seconds(&timing.traceback_sort_prelim_ns)
    );
    eprintln!(
        "[TIMING] traceback_tree_precheck: {:.3}s",
        timing_seconds(&timing.traceback_tree_precheck_ns)
    );
    eprintln!(
        "[TIMING] traceback_start_offsets: {:.3}s",
        timing_seconds(&timing.traceback_start_offsets_ns)
    );
    eprintln!(
        "[TIMING] traceback_alignment: {:.3}s (dp={:.3}s greedy={:.3}s)",
        timing_seconds(&timing.traceback_alignment_ns),
        timing_seconds(&timing.traceback_alignment_dp_ns),
        timing_seconds(&timing.traceback_alignment_greedy_ns)
    );
    eprintln!(
        "[TIMING] traceback_hsp_build: {:.3}s",
        timing_seconds(&timing.traceback_hsp_build_ns)
    );
    eprintln!(
        "[TIMING] purge_common_endpoint_pass1: {:.3}s",
        timing_seconds(&timing.purge_common_endpoint_pass1_ns)
    );
    eprintln!(
        "[TIMING] reevaluate_ambiguities: {:.3}s",
        timing_seconds(&timing.reevaluate_ambiguities_ns)
    );
    eprintln!(
        "[TIMING] identity_length_test: {:.3}s",
        timing_seconds(&timing.identity_length_test_ns)
    );
    eprintln!(
        "[TIMING] purge_common_endpoint_pass2: {:.3}s",
        timing_seconds(&timing.purge_common_endpoint_pass2_ns)
    );
    eprintln!(
        "[TIMING] traceback_score_resort: {:.3}s",
        timing_seconds(&timing.traceback_score_resort_ns)
    );
    eprintln!(
        "[TIMING] traceback_final_tree_contains: {:.3}s",
        timing_seconds(&timing.traceback_final_tree_contains_ns)
    );
    eprintln!(
        "[TIMING] traceback_final_tree_add: {:.3}s",
        timing_seconds(&timing.traceback_final_tree_add_ns)
    );
    eprintln!(
        "[TIMING] traceback_postprocess: {:.3}s",
        timing_seconds(&timing.traceback_postprocess_ns)
    );
    eprintln!(
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
    eprintln!(
        "[TIMING] traceback_lengths: edit_script_samples={} edit_script_max={} edit_script_median={:.1} alignment_samples={} alignment_max={} alignment_median={:.1}",
        edit_count, edit_max, edit_median, aln_count, aln_max, aln_median
    );
    eprintln!(
        "[TIMING] format_output: {:.3}s (calls={})",
        format_s, format_n
    );
    if let Some(t_search_start) = t_search_start {
        eprintln!(
            "[TIMING] search_total: {:.3}s",
            t_search_start.elapsed().as_secs_f64()
        );
    }
    if let Some(t_total) = t_total {
        eprintln!("[TIMING] total: {:.3}s", t_total.elapsed().as_secs_f64());
    }
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

    let mut per_query: Vec<Vec<BlastnHsp>> = vec![Vec::new(); hit_lists.len()];
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

/// What the pairwise report, the hit records and the query warnings need besides the
/// final hit list.
struct BlastnReportInputs<'a> {
    queries: &'a [bio::io::fasta::Record],
    subjects: &'a [bio::io::fasta::Record],
    /// Per query, plus strand: DUST and, with `-lcase_masking`, the input lowercase.
    query_masks: &'a [Vec<MaskedInterval>],
    lcase_masking: bool,
    /// Per query context (the plus strand of query `q` is context `2 * q`).
    query_eff_searchsp: &'a [i64],
    megablast: bool,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    /// Per query context (the plus strand of query `q` is context `2 * q`); `None` for
    /// the contexts of an invalid query (not searched by NCBI).
    context_karlin: &'a [Option<ContextKarlin>],
    max_target_seqs: Option<usize>,
    db_num_seqs: usize,
    db_len_total: usize,
}

/// The per-query and run-level data of the pairwise report.
// Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
// shared post-processing makes every BLASTN run compile it in Wasm hosts.
#[inline(never)]
fn blastn_pairwise_report(
    report: &BlastnReportInputs<'_>,
    unsearched: Vec<bool>,
    query_titles: &[Arc<str>],
    subject_title: &str,
) -> Result<(Vec<BlastnPairwiseQuery>, BlastnPairwiseReport)> {
    let context_karlin = report.context_karlin;
    let queries = report
        .queries
        .iter()
        .enumerate()
        .map(|(q_idx, query)| BlastnPairwiseQuery {
            query_name: query_titles[q_idx].to_string(),
            query_length: query.seq().len(),
            karlin: context_karlin[2 * q_idx].map(|blocks| (blocks.ungapped, blocks.gapped)),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_results.cpp:82-104
            // ```c
            //     // find the first valid context corresponding to this query
            //     ...
            //     m_SearchSpace = ctx->eff_searchsp;
            // ```
            // An invalid query has no valid context, and its search space stays 0.
            effective_search_space: if context_karlin[2 * q_idx].is_some() {
                report.query_eff_searchsp[2 * q_idx]
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
    let (num_descriptions, num_alignments) = match report.max_target_seqs {
        Some(max_target_seqs) => (max_target_seqs, max_target_seqs),
        None => (500, 250),
    };
    let pairwise_report = BlastnPairwiseReport {
        version: NCBI_BLASTN_VERSION.to_string(),
        megablast: report.megablast,
        database_name: subject_title.to_string(),
        database_num_sequences: report.db_num_seqs,
        database_total_letters: report.db_len_total,
        reward: report.reward,
        penalty: report.penalty,
        gap_open: report.gap_open,
        gap_extend: report.gap_extend,
        num_descriptions,
        num_alignments,
        unsearched,
    };
    Ok((queries, pairwise_report))
}

fn post_process_hits_and_write(
    mut hit_lists: Vec<Option<BlastnHitList>>,
    hitlist_size: usize,
    max_hsps_per_subject: usize,
    subject_besthit: bool,
    query_lengths: &[usize],
    outputs: &mut ReportOutputs<'_>,
    verbose: bool,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    output_formats: &[BlastnOutputFormat],
    query_titles: &[Arc<str>],
    subject_title: &str,
    timing: Option<&BlastnTiming>,
    report: &BlastnReportInputs<'_>,
    unsearched: Vec<bool>,
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
    for (q_idx, hit_list_opt) in hit_lists.iter_mut().enumerate() {
        let hit_list = match hit_list_opt {
            Some(value) => value,
            None => continue,
        };
        for hsp_list in &mut hit_list.hsplist_array {
            if max_hsps_per_subject > 0 {
                trim_by_max_hsps(hsp_list, max_hsps_per_subject);
            }
            if subject_besthit {
                let query_len = query_lengths.get(q_idx).copied().unwrap_or_default();
                subject_best_hit(&mut hsp_list.hsps, query_len);
            }
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
    let pairwise_hits =
        if outputs.hits.is_some() || output_formats.contains(&BlastnOutputFormat::Pairwise) {
            let subject_masks = report.lcase_masking.then(|| {
                report
                    .subjects
                    .iter()
                    .map(|subject| collect_lowercase_masks(subject.seq()))
                    .collect::<Vec<_>>()
            });
            Some(pairwise_hits(
                &hit_lists,
                report.queries,
                report.subjects,
                &DisplayMasks {
                    query: report.query_masks,
                    subject: subject_masks.as_deref(),
                },
            )?)
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
    let observer = &mut outputs.observer;
    for (format_index, (format, &output_format)) in
        outputs.formats.iter_mut().zip(output_formats).enumerate()
    {
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:68-93
        // ```c
        // CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
        //     : m_FormatType(format_type), ..., m_Outfile(outfile),
        // ```
        let mut writer = format.sink.open()?;
        if output_format == BlastnOutputFormat::Pairwise {
            let hits = pairwise_hits
                .as_deref()
                .expect("pairwise hits are built for outfmt 0");
            let (queries, pairwise_report) = pairwise_data
                .as_ref()
                .expect("the pairwise report data is built for outfmt 0");
            write_blastn_pairwise_report(
                hits,
                &mut writer,
                &PairwiseConfig {
                    program: "blastn".to_string(),
                    show_frame: false,
                    ..PairwiseConfig::default()
                },
                queries,
                subject_ids,
                pairwise_report,
                probe.as_mut(),
            )?;
        } else {
            write_output_blastn_hitlists_to_writer(
                &hit_lists,
                &mut writer,
                query_ids,
                subject_ids,
                output_format,
                query_titles,
                subject_title,
                &unsearched,
                probe.as_mut(),
            )?;
        }
        writer.flush()?;
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1450-1452
    // ```c
    //     if (results.HasWarnings()) {
    //         ERR_POST(Warning << results.GetWarningStrings());
    //     }
    // ```
    // The warnings belong to the result, so they are written once, in query order.
    for (index, (query, karlin)) in report
        .queries
        .iter()
        .zip(report.context_karlin.iter().step_by(2))
        .enumerate()
    {
        if karlin.is_none() {
            outputs
                .diagnostics
                .write_all(&invalid_query_warning("blastn", index, query))?;
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

#[cfg(target_arch = "wasm32")]
fn fasta_records_from_bytes(bytes: &[u8]) -> Result<Vec<bio::io::fasta::Record>> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // void
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    bio::io::fasta::Reader::new(bytes)
        .records()
        // NCBI reference: c++/src/objtools/readers/fasta.cpp:428-431
        // FASTA_ERROR(LineNumber(), "CFastaReader: Expected defline around line " << LineNumber(), ...);
        .collect::<std::result::Result<Vec<_>, _>>()
        .context("failed to parse in-memory FASTA")
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
// ```c
// CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
// CRef<CBlastInput> input(new CBlastInput(fasta));
// CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
// sequences = input->GetAllSeqs(*scope);
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/readers/fasta.cpp:428-431
// ```c
// if (need_defline  &&  GetLineReader().AtEOF()) {
//     FASTA_ERROR(LineNumber(),
//         "CFastaReader: Expected defline around line " << LineNumber(),
//         CObjReaderParseException::eEOF);
// ```
fn read_blastn_fasta_records(
    path: &std::path::Path,
    role: &str,
) -> Result<Vec<bio::io::fasta::Record>> {
    let mut file = open_blastn_input(path, role)?;
    parse_blastn_fasta(&read_blastn_fasta_bytes(&mut file, path, role)?, path, role)
}

/// Opens a BLASTN input as NCBI's argument does when a handler asks for its stream: `-` is
/// standard input, and a file that does not open gets NCBI's error
/// (`crate::cli::inaccessible`).
///
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbiargs.cpp:717-735
/// ```c
///     if (AsString() == "-") {
/// #if defined(NCBI_OS_MSWIN)
///         NcbiSys_setmode(NcbiSys_fileno(stdin), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
/// #endif
///         m_Ios  = &cin;
///     } else if ( !AsString().empty() ) {
///         if (!fstrm) {
///             fstrm = new CNcbiIfstream;
///         }
///         if (fstrm) {
///             fstrm->open(AsString().c_str(),IOS_BASE::in | mode);
///             if ( !fstrm->is_open() ) {
///                 delete fstrm;
///                 fstrm = NULL;
///             } else {
///                 m_DeleteFlag = true;
///             }
///         }
///         m_Ios = fstrm;
///     }
/// ```
fn open_blastn_input(path: &std::path::Path, role: &str) -> Result<std::fs::File> {
    if path.as_os_str() == "-" {
        return standard_input().map_err(|_| {
            anyhow::anyhow!(
                "reading the {role} from standard input ('-') on this platform is not supported by LOSAT's BLASTN"
            )
        });
    }
    std::fs::File::open(path).map_err(|_| crate::cli::inaccessible(role, path))
}

/// Standard input as a file that shares its position, as `cin` does.
fn standard_input() -> std::io::Result<std::fs::File> {
    #[cfg(any(unix, target_os = "wasi"))]
    {
        use std::os::fd::AsFd;
        std::io::stdin()
            .as_fd()
            .try_clone_to_owned()
            .map(std::fs::File::from)
    }
    #[cfg(windows)]
    {
        use std::os::windows::io::AsHandle;
        std::io::stdin()
            .as_handle()
            .try_clone_to_owned()
            .map(std::fs::File::from)
    }
    #[cfg(not(any(unix, windows, target_os = "wasi")))]
    {
        Err(std::io::ErrorKind::Unsupported.into())
    }
}

/// The bytes of an opened FASTA file; a directory reads as no bytes, as NCBI's stream does.
fn read_blastn_fasta_bytes(
    file: &mut std::fs::File,
    path: &std::path::Path,
    role: &str,
) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    match std::io::Read::read_to_end(file, &mut bytes) {
        Ok(_) => Ok(bytes),
        Err(error) if error.kind() == std::io::ErrorKind::IsADirectory => Ok(Vec::new()),
        Err(error) => {
            Err(error).with_context(|| format!("failed to read {role} FASTA {}", path.display()))
        }
    }
}

/// The records of a FASTA file, after rejecting the residues that NCBI reads differently
/// (`input.rs`); a file of white space only has no record, as in NCBI. NCBI reads the
/// deflines that LOSAT rejects (`check_deflines`) without a message, so the callers check
/// them where the difference would change a result.
fn read_blastn_records(
    bytes: &[u8],
    path: &std::path::Path,
    role: &str,
) -> Result<Vec<bio::io::fasta::Record>> {
    if is_blank(bytes) {
        return Ok(Vec::new());
    }
    check_sequence_lines(bytes, role)?;
    let records = bio::io::fasta::Reader::new(bytes)
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| {
            format!(
                "failed to read {role} FASTA {} ({UNREADABLE_FASTA})",
                path.display()
            )
        })?;
    check_residues(&records, role)?;
    Ok(records)
}

/// The records of a FASTA file read where its deflines matter (the query), with the
/// deflines checked first, so that a defline that bio cannot read is named.
fn parse_blastn_fasta(
    bytes: &[u8],
    path: &std::path::Path,
    role: &str,
) -> Result<Vec<bio::io::fasta::Record>> {
    if !is_blank(bytes) {
        check_deflines(bytes, role)?;
    }
    let records = read_blastn_records(bytes, path, role)?;
    check_records_have_residues(&records, role)?;
    Ok(records)
}

/// ABI v1's reading of a BLASTN `-outfmt` value, which is frozen (plan TD-1): it keeps the
/// values that it accepted before `parse_blastn_output_format` followed NCBI.
#[cfg(target_arch = "wasm32")]
fn v1_output_format(spec: &str) -> Result<BlastnOutputFormat> {
    let mut parts = spec.split_whitespace();
    let format = parts.next().unwrap_or("6");
    if parts.next().is_some() {
        anyhow::bail!("unsupported BLASTN custom outfmt specification: {spec:?}");
    }
    match format {
        "0" => Ok(BlastnOutputFormat::Pairwise),
        "6" => Ok(BlastnOutputFormat::Tabular),
        "7" => Ok(BlastnOutputFormat::TabularWithComments),
        _ => anyhow::bail!("unsupported BLASTN output format: {format}"),
    }
}

#[cfg(target_arch = "wasm32")]
pub fn run_web_pair(args: BlastnArgs, query_fasta: &str, subject_fasta: &str) -> Result<Vec<u8>> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
    // ```c
    // BlastSeqBlkSetSequence(subj, sequence.data.release(),
    //    ((sentinels == eSentinels) ? sequence.length - 2 :
    //     sequence.length));
    // ...
    // SBlastSequence compressed_seq =
    //     subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
    //                               eNa_strand_plus, eNoSentinels);
    // ```
    let queries = fasta_records_from_bytes(query_fasta.as_bytes()).context("query FASTA")?;
    let subjects = fasta_records_from_bytes(subject_fasta.as_bytes()).context("subject FASTA")?;
    // NCBI reference: ncbi-blast/c++/include/algo/blast/api/local_blast.hpp:76-78
    // ```c
    // CLocalBlast(CRef<IQueryFactory> query_factory,
    //             CRef<CBlastOptionsHandle> opts_handle,
    //             CRef<CLocalDbAdapter> db);
    // ```
    // Web ABI v1 runs the same local search as the CLI and keeps the report in memory.
    let mut output = Vec::new();
    let outfmt = args.outfmt.clone();
    // Plan TD-1: ABI v1 is frozen except for fail-fast fixes, so its BLASTN keeps the
    // formats that it had (6 and 7) and rejects outfmt 0 with the error that it gave
    // before the engine implemented outfmt 0.
    let outfmt = match v1_output_format(&outfmt)? {
        BlastnOutputFormat::Pairwise => anyhow::bail!("unsupported BLASTN output format: 0"),
        BlastnOutputFormat::Tabular => "6".to_string(),
        BlastnOutputFormat::TabularWithComments => "7".to_string(),
    };
    // ABI v1 is frozen (plan TD-1) except fail-fast fixes. Before S07+ it checked the
    // thread count first (in its search pool) and gave the empty report of an empty
    // query, as NCBI does after reading the subjects and checking the options; S07+ keeps
    // NCBI's errors there (fail-fast) and adds its checks of the deflines, the records and
    // LOSAT's limits only for a search.
    crate::utils::threading::validate_threads(args.num_threads)?;
    if queries.is_empty() {
        check_subjects_not_empty(&subjects)?;
        check_scoring_options(&args)?;
        return Ok(output);
    }
    // The deflines that NCBI reads differently are rejected (a fail-fast fix, plan TD-1).
    check_deflines(subject_fasta.as_bytes(), "subject")?;
    check_deflines(query_fasta.as_bytes(), "query")?;
    check_sequence_lines(subject_fasta.as_bytes(), "subject")?;
    check_sequence_lines(query_fasta.as_bytes(), "query")?;
    if !subjects.is_empty() {
        check_scoring_options(&args)?;
        check_losat_limits(&args)?;
        check_records(&subjects, "subject")?;
        check_records(&queries, "query")?;
    }
    let mut stderr = std::io::stderr();
    run_local(
        args,
        &queries,
        &subjects,
        &mut ReportOutputs::single(&outfmt, OutputSink::Writer(&mut output), &mut stderr),
    )?;
    Ok(output)
}

pub fn run(args: BlastnArgs) -> Result<()> {
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
    // option handlers run.
    let output_formats = parse_output_formats([args.outfmt.as_str()])?;
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
    let mut subject_file = open_blastn_input(&args.subject, "subject")?;
    let subject_bytes = read_blastn_fasta_bytes(&mut subject_file, &args.subject, "subject")?;
    drop(subject_file);
    let subjects = read_blastn_records(&subject_bytes, &args.subject, "subject")?;
    write_title_warnings(&subjects, &mut std::io::stderr())?;
    // NCBI reads these deflines and records without a message; LOSAT rejects them where
    // the search would start.
    let subject_deflines = check_deflines(&subject_bytes, "subject")
        .and_then(|()| check_records_have_residues(&subjects, "subject"));
    drop(subject_bytes);
    check_subjects_not_empty(&subjects)?;
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
    let mut query_file = open_blastn_input(&args.query, "query")?;
    // The output file is opened once, as NCBI's stream (a named pipe gives one reader one
    // end of file).
    let mut out_file = match args.out.as_deref().filter(|path| path.as_os_str() != "-") {
        Some(path) => Some(std::io::BufWriter::new(
            std::fs::File::create(path).map_err(|_| crate::cli::inaccessible("out", path))?,
        )),
        None => None,
    };
    let outfmt = args.outfmt.clone();
    let mut stderr = std::io::stderr();
    let result = {
        let sink = match out_file.as_mut() {
            Some(file) => OutputSink::Writer(file),
            None => OutputSink::Stdout,
        };
        let mut outputs = ReportOutputs::single(&outfmt, sink, &mut stderr);
        search_cli(
            args,
            query_file,
            &subjects,
            subject_deflines,
            &mut outputs,
            output_formats,
        )
    };
    // What was written before an error (such as the outfmt 0 prolog) stays in the file.
    let flushed = out_file.as_mut().map_or(Ok(()), std::io::Write::flush);
    result?;
    flushed.context("failed to write the output")
}

/// The part of `run` after the output is opened: NCBI's filtering handler (`-dust`), the
/// processing of the options, `Query is Empty!`, LOSAT's deferred checks of the subjects,
/// and the search.
fn search_cli(
    mut args: BlastnArgs,
    mut query_file: std::fs::File,
    subjects: &[bio::io::fasta::Record],
    subject_deflines: Result<()>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
) -> Result<()> {
    args.resolve_dust()?;
    process_options(&args, outputs)?;
    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:209-212
    // ```c
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())) {
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    // ```
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:856-860
    // ```c
    // 	char c;
    // 	CNcbiStreampos orig_p = in.tellg();
    // 	// Piped input
    // 	if(orig_p < 0)
    // 		return false;
    // ```
    // The position is taken on the opened file, before it is read. When the subjects were
    // read from standard input too, `cin` has reached its end (a failed stream), so NCBI
    // gets no position for the query either.
    let seekable = !(args.query.as_os_str() == "-" && args.subject.as_os_str() == "-")
        && std::io::Seek::stream_position(&mut query_file).is_ok();
    let query_bytes = read_blastn_fasta_bytes(&mut query_file, &args.query, "query")?;
    drop(query_file);
    if is_blank(&query_bytes) {
        // NCBI reads a stream without a position (a pipe) as not empty, and then prints
        // the report of no query.
        if !seekable {
            anyhow::bail!(
                "an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's BLASTN"
            );
        }
        outputs
            .diagnostics
            .write_all(b"Warning: [blastn] Query is Empty!\n")?;
        return Ok(());
    }
    subject_deflines?;
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
    // ```c
    // sequences = input->GetAllSeqs(*scope);
    // ```
    let queries = parse_blastn_fasta(&query_bytes, &args.query, "query")?;
    search(args, &queries, subjects, outputs, output_formats)
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
fn check_subjects_not_empty(subjects: &[bio::io::fasta::Record]) -> Result<()> {
    if subjects.is_empty() {
        return Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
        }
        .into());
    }
    Ok(())
}

/// The requested output formats, each parsed exactly as the single CLI `-outfmt` is,
/// before NCBI's option handlers run (`run`).
fn parse_output_formats<'a>(
    formats: impl IntoIterator<Item = &'a str>,
) -> Result<Vec<BlastnOutputFormat>> {
    formats
        .into_iter()
        .map(parse_blastn_output_format)
        .collect()
}

/// NCBI's processing of the options after the query and the output are opened
/// (`SetOptions`): the few-matches warning of the formatting handler, and the check of the
/// options.
///
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
fn process_options(args: &BlastnArgs, outputs: &mut ReportOutputs<'_>) -> Result<()> {
    if args
        .max_target_seqs
        .is_some_and(|max_target_seqs| max_target_seqs < 5)
    {
        outputs
            .diagnostics
            .write_all(&few_matches_warning("blastn"))?;
    }
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
pub fn run_local(
    mut args: BlastnArgs,
    query_records: &[bio::io::fasta::Record],
    subject_records: &[bio::io::fasta::Record],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    let output_formats = parse_output_formats(outputs.formats.iter().map(|format| format.outfmt))?;
    write_title_warnings(subject_records, outputs.diagnostics)?;
    check_subjects_not_empty(subject_records)?;
    args.resolve_dust()?;
    process_options(&args, outputs)?;
    search(
        args,
        query_records,
        subject_records,
        outputs,
        output_formats,
    )
}

fn search(
    args: BlastnArgs,
    query_records: &[bio::io::fasta::Record],
    subject_records: &[bio::io::fasta::Record],
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
) -> Result<()> {
    // NCBI reads the subjects before the queries (`run`). Residues that NCBI reads
    // differently are rejected, and `U` is read as `T` (`input.rs`); the callers check the
    // deflines, in the bytes of the files, and an empty subject set. The CLI has already
    // reported an empty query file; records of the other callers get NCBI's warning.
    check_residues(subject_records, "subject")?;
    if query_records.is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [blastn] Query is Empty!\n")?;
        return Ok(());
    }
    check_records_have_residues(subject_records, "subject")?;
    // NCBI reads the queries after `Query is Empty!`, with its reader's warnings.
    write_title_warnings(query_records, outputs.diagnostics)?;
    // NCBI decodes HTML character references in the outfmt 0 titles of the subjects
    // (`NStr::HtmlDecode` in `CDeflineGenerator::GenerateDefline`, create_defline.cpp:4066),
    // which LOSAT does not reproduce (`report/defline.rs`).
    // NCBI's x_CleanAndCompress also reads past the end of some titles of punctuation
    // (NCBI crashes when it writes the title of such a subject with hits), which LOSAT
    // does not reproduce (`report/defline.rs`); LOSAT rejects the subject before the search.
    if output_formats.contains(&BlastnOutputFormat::Pairwise) {
        for (index, record) in subject_records.iter().enumerate() {
            let defline = match record.desc() {
                Some(desc) => format!("{} {desc}", record.id()),
                None => record.id().to_string(),
            };
            if crate::report::defline::has_html_character_reference(&defline) {
                anyhow::bail!(
                    "subject record {} has an HTML character reference (such as &amp;) in its defline, which NCBI BLAST+ decodes in the outfmt 0 titles; this is not supported by LOSAT's BLASTN",
                    index + 1
                );
            }
            if [false, true].into_iter().any(|leave_prefix| {
                crate::report::defline::ncbi_nucleotide_title(&defline, leave_prefix).is_none()
            }) {
                anyhow::bail!(
                    "subject record {} has a defline of punctuation that NCBI BLAST+ reads past its end when it writes the subject's outfmt 0 title (it crashes if the subject has hits), which LOSAT does not reproduce",
                    index + 1
                );
            }
        }
    }
    // LOSAT's limits come where NCBI starts the search, after its checks and its
    // `Query is Empty!` success.
    check_losat_limits(&args)?;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:86-91
    // ```c
    //     char* batch_sz_str = getenv("BATCH_SIZE");
    //     if (batch_sz_str) {
    //         retval = NStr::StringToInt(batch_sz_str);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:59-62
    // ```c
    //     char* chunk_sz_str = getenv("CHUNK_SIZE");
    //     if (chunk_sz_str && !NStr::IsBlank(chunk_sz_str)) {
    //         retval = NStr::StringToInt(chunk_sz_str);
    // ```
    // LOSAT follows NCBI's query batches without these variables (a blank CHUNK_SIZE is
    // none).
    for variable in ["BATCH_SIZE", "CHUNK_SIZE"] {
        if std::env::var_os(variable)
            .is_some_and(|value| variable == "BATCH_SIZE" || !is_blank(value.as_encoded_bytes()))
        {
            anyhow::bail!(
                "the environment variable {variable}, which changes NCBI BLAST+'s query batches, is not supported by LOSAT's BLASTN"
            );
        }
    }
    check_residues(query_records, "query")?;
    check_records_have_residues(query_records, "query")?;
    let subjects_read = with_u_as_t(subject_records);
    let queries_read = with_u_as_t(query_records);
    let subject_records = subjects_read.as_deref().unwrap_or(subject_records);
    let query_records = queries_read.as_deref().unwrap_or(query_records);
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
    // `int`, which LOSAT does not reproduce beyond its range.
    let subject_letters: usize = subject_records
        .iter()
        .map(|record| record.seq().len())
        .sum();
    if subject_letters > i32::MAX as usize {
        anyhow::bail!(
            "the subjects have {subject_letters} letters; NCBI BLAST+ reports a total of 2^31 letters or more through a 32-bit int, which is not supported by LOSAT's BLASTN"
        );
    }
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    crate::utils::threading::with_search_pool(args.num_threads, "blastn", |pool| {
        run_in_pool(
            args,
            query_records,
            subject_records,
            outputs,
            output_formats,
            pool,
        )
    })
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads());
// (*thread)->Run(); (*thread)->Join(&result);
/// The subjects of a search, read and encoded once for all query batches, as NCBI's
/// `CLocalDbAdapter` of the subject set serves every batch.
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
    records: &'a [bio::io::fasta::Record],
    metadata: SubjectMetadata,
    blastna: Vec<Vec<u8>>,
    packed: Vec<Vec<u8>>,
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
}

/// Which of NCBI's stages `search_query_batch` runs for its queries.
#[derive(Clone, Copy)]
enum BatchStage<'a> {
    /// The preliminary search and the traceback of a query batch (`CLocalBlast::Run`). A
    /// batch that NCBI splits has its preliminary search run in query chunks
    /// (`search_query_chunks`).
    Search,
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
/// `GetQueryBatchSize` is 0 for blastn (the environment variable `BATCH_SIZE` is rejected);
/// the chunk size is `query_chunk_size`.
fn run_in_pool(
    args: BlastnArgs,
    query_records: &[bio::io::fasta::Record],
    subject_records: &[bio::io::fasta::Record],
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
) -> Result<()> {
    if query_records.is_empty() {
        return Ok(());
    }
    let metadata = subject_metadata_from_records(subject_records);
    if metadata.db_num_seqs == 0 {
        return Ok(());
    }
    let subjects = PreparedSubjects {
        records: subject_records,
        blastna: subject_records
            .iter()
            .map(|record| encode_iupac_to_blastna(record.seq()))
            .collect(),
        packed: subject_records
            .iter()
            .map(|record| encode_subject_ncbi2na_packed(record.seq()))
            .collect(),
        metadata,
    };
    let megablast = args.task == "megablast";
    let chunk_size = query_chunk_size(megablast) as i32;
    let mut mixer = BatchSizeMixer::new(chunk_size - 1000);
    let total_length = subjects.metadata.db_len_total as i64;
    if total_length > 0 {
        mixer.set_target_hits((total_length / 3000) as i32);
    }
    let mut batch_size = mixer.batch_size(None);

    let lengths: Vec<usize> = query_records
        .iter()
        .map(|record| record.seq().len())
        .collect();
    let mut hit_lists: Vec<Option<BlastnHitList>> = Vec::with_capacity(lengths.len());
    let mut unsearched = Vec::with_capacity(lengths.len());
    let mut query_masks = Vec::with_capacity(lengths.len());
    let mut context_karlin = Vec::with_capacity(2 * lengths.len());
    let mut query_eff_searchsp = Vec::with_capacity(2 * lengths.len());
    let mut start = 0;
    while start < lengths.len() {
        let end = next_query_batch_end(&lengths, start, batch_size);
        let batch = search_query_batch(
            &args,
            &query_records[start..end],
            &subjects,
            outputs,
            &output_formats,
            parallel_pool,
            batch_size,
            start == 0,
            BatchStage::Search,
        )?;
        // The batch numbers its queries from 0.
        let first = start as u32;
        hit_lists.extend(batch.hit_lists.into_iter().map(|list| {
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
        unsearched.extend(std::iter::repeat_n(!batch.searched, end - start));
        query_masks.extend(batch.query_masks);
        context_karlin.extend(batch.context_karlin);
        query_eff_searchsp.extend(batch.query_eff_searchsp);
        // NCBI's `Int4` count of the batch's successful initial extensions.
        batch_size = mixer.batch_size(Some(batch.good_init_extends as u32 as i32));
        start = end;
    }

    let config = configure_task(&args);
    let report_inputs = BlastnReportInputs {
        queries: query_records,
        subjects: subject_records,
        query_masks: &query_masks,
        lcase_masking: args.lcase_masking,
        query_eff_searchsp: &query_eff_searchsp,
        megablast,
        reward: config.reward,
        penalty: config.penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
        context_karlin: &context_karlin,
        max_target_seqs: args.max_target_seqs,
        db_num_seqs: subjects.metadata.db_num_seqs,
        db_len_total: subjects.metadata.db_len_total,
    };
    let query_ids: Vec<Arc<str>> = query_records
        .iter()
        .map(|record| Arc::from(record.id().split_whitespace().next().unwrap_or("unknown")))
        .collect();
    let subject_ids: Vec<Arc<str>> = subjects
        .metadata
        .subject_ids
        .iter()
        .map(|id: &String| Arc::<str>::from(id.as_str()))
        .collect();
    let query_titles: Vec<Arc<str>> = query_records.iter().map(fasta_defline).collect();
    let subject_title = format!(
        "User specified sequence set (Input: {})",
        args.subject.display()
    );
    let hitlist_size = match args.max_target_seqs {
        Some(max_target_seqs) if max_target_seqs > 0 => max_target_seqs,
        _ => args.hitlist_size,
    };
    post_process_hits_and_write(
        hit_lists,
        hitlist_size,
        args.max_hsps_per_subject.unwrap_or(0),
        args.subject_besthit,
        &lengths,
        outputs,
        args.verbose,
        &query_ids,
        &subject_ids,
        &output_formats,
        &query_titles,
        &subject_title,
        None,
        &report_inputs,
        unsearched,
    )
}

/// Searches one query batch against the subjects, with the query block, lookup table and
/// diagonal table of the batch, as NCBI's `CLocalBlast` does for each batch that
/// `run_in_pool` reads.
fn search_query_batch(
    args: &BlastnArgs,
    query_records: &[bio::io::fasta::Record],
    subjects: &PreparedSubjects<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batch_size: i32,
    first_batch: bool,
    stage: BatchStage<'_>,
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
    let speculative_traceback = requested_parallel
        && config.use_dp
        && std::env::var_os("LOSAT_DEBUG_COORDS").is_none()
        && std::env::var_os("LOSAT_DEBUG_COORDS_START").is_none();
    // Read sequences
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    // The search owns its copy of the query records (`prepare_sequence_data`).
    let queries = query_records.to_vec();
    let query_ids = queries
        .iter()
        .map(|record| {
            record
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string()
        })
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
    let subject_records: Option<&[bio::io::fasta::Record]> = Some(subjects.records);
    let subject_metadata = subjects.metadata.clone();

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
    let use_parallel = requested_parallel && subject_metadata.db_num_seqs > 1;
    crate::utils::threading::report_stage(
        "blastn",
        "subjects",
        subject_metadata.db_num_seqs,
        use_parallel,
    );

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
    let subject_blastna_cache_ref = Some(subjects.blastna.as_slice());
    let subject_packed_cache_ref = Some(subjects.packed.as_slice());

    // Prepare sequence data (including DUST masking)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1401-1405
    // ```c
    // /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
    //   sequences are retieved in ncbistdaa/ncbi2na encodings respectively. */
    // seq_arg.encoding = eBlastEncodingProtein;
    // ```
    let masks = match stage {
        BatchStage::Search => query_masks(args, &queries),
        BatchStage::ChunkPrelim { query_masks, .. } => query_masks.to_vec(),
    };
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
            .map(fasta_defline)
            .collect::<Vec<Arc<str>>>(),
    );
    let subject_title = Arc::<str>::from(format!(
        "User specified sequence set (Input: {})",
        args.subject.display()
    ));

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2589-2593
    // ```c
    // curr_context = hsp_array[i]->context;
    // qlen = query_info->contexts[curr_context].query_length;
    // target_context = (hsp_array[i]->query.frame > 0) ?
    //                  curr_context + 1 : curr_context - 1;
    // ```
    let query_lengths = Arc::new(
        seq_data
            .queries
            .iter()
            .map(|record| record.seq().len())
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
    let subject_ids_arc = Arc::new(
        seq_data
            .subject_ids
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
    let mut query_concat_offset: usize = 0;
    for (q_idx, q_record) in seq_data.queries.iter().enumerate() {
        let seq = q_record.seq();
        let q_len = seq.len();
        let rc_seq = reverse_complement(seq);
        let plus_masks = seq_data.query_masks[q_idx].clone();
        let minus_masks = reverse_mask_intervals(&plus_masks, q_len);
        let plus_offset = query_concat_offset;
        let minus_offset = query_concat_offset + q_len + 1;

        query_base_offsets.push(plus_offset);
        query_contexts.push(QueryContext {
            query_idx: q_idx as u32,
            frame: 1,
            query_offset: plus_offset as i32,
            seq: seq.to_vec(),
            masks: plus_masks.clone(),
        });
        query_contexts.push(QueryContext {
            query_idx: q_idx as u32,
            frame: -1,
            query_offset: minus_offset as i32,
            seq: rc_seq.clone(),
            masks: minus_masks.clone(),
        });

        query_context_masks.push(plus_masks);
        query_context_masks.push(minus_masks);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:65-77
        // ```c
        // Uint4 prev_loc = qinfo->contexts[index-1].query_offset;
        // Uint4 prev_len = qinfo->contexts[index-1].query_length;
        // Uint4 shift = prev_len ? prev_len + 1 : 0;
        // qinfo->contexts[index].query_offset = prev_loc + shift;
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:246-249
        // ```c
        // BlastContextInfo * cinfo = & qinfo->contexts[qinfo->last_context];
        // return cinfo->query_offset + cinfo->query_length + (cinfo->query_length ? 2 : 1);
        // ```
        query_concat_offset += q_len * 2 + 1;
        if q_idx + 1 < seq_data.queries.len() {
            query_concat_offset += 1;
        }
    }
    let query_concat_length = query_concat_offset;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:166-174
    // ```c
    // query_info->contexts[0].query_offset = 0;
    // query_info->contexts[1].query_offset = kStrandLength + 1;
    // ```
    let query_context_offsets: Vec<i32> =
        query_contexts.iter().map(|ctx| ctx.query_offset).collect();
    let query_context_index = QueryContextIndex::new(&query_contexts);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:498-604
    // ```c
    // EBlastEncoding encoding = GetQueryEncoding(prog);
    // ...
    // sequence = queries.GetBlastSequence(index, encoding, strand, eSentinels);
    // ...
    // memcpy(&buf.get()[offset], sequence.data.get(),
    //        sequence.length);
    // ```
    let encoded_queries_blastna: Vec<Vec<u8>> = query_contexts
        .iter()
        .map(|ctx| encode_iupac_to_blastna(ctx.seq.as_slice()))
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:500-603
    // ```c
    // int buflen = QueryInfo_GetSeqBufLen(qinfo);
    // TAutoUint1Ptr buf((Uint1*) calloc(buflen+1, sizeof(Uint1)));
    // ...
    // memcpy(&buf.get()[offset], sequence.data.get(), sequence.length);
    // ```
    let (encoded_query_concat_blastna, encoded_query_concat_blastna_with_sentinels) =
        build_query_blastna_concat_buffers(
            &encoded_queries_blastna,
            &query_context_offsets,
            query_concat_length,
        );

    // NCBI reference: c++/src/algo/blast/core/na_ungapped.c:294,324,740-749
    // Uint1 q_byte = (q[0] << 6) | (q[1] << 4) | (q[2] << 2) | q[3];
    // if (word_params->matrix_only_scoring || word_length < 11)
    //     s_NuclUngappedExtendExact(...);
    // Use the same encoded concat/context slices as extension, including their
    // sentinels. Exact-only word sizes do not allocate this search-local table.
    let query_four_base = if config.effective_word_size >= 11 {
        build_query_four_base_bytes(&encoded_query_concat_blastna)
    } else {
        Vec::new()
    };

    // Finalize configuration with query-dependent parameters (adaptive lookup table selection)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:46-47
    // ```c
    // BlastChooseNaLookupTable(const LookupTableOptions* lookup_options,
    //                          Int4 approx_table_entries, Int4 max_q_off,
    //                          Int4 *lut_width)
    // ```
    let discontig_template = args.task == "dc-megablast";
    let scoring_spec = NuclScoringSpec {
        reward: config.reward,
        penalty: config.penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
    };
    let context_ungapped = context_ungapped_blocks(
        query_contexts.iter().map(|context| context.seq.as_slice()),
        &scoring_spec,
    );
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
    let diag_array_length = {
        let mut diag_array_length = 1usize;
        while diag_array_length < query_concat_length + TWO_HIT_WINDOW {
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
    let (context_karlin, round_down_evalue_score) = match context_blocks(
        &context_ungapped,
        &scoring_spec,
    ) {
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
            if !first_batch {
                anyhow::bail!(
                    "these scoring options have no Karlin-Altschul values, and the first query batch has only invalid queries; NCBI BLAST+ reports the error after the results of that batch, which LOSAT does not reproduce"
                );
            }
            if context_ungapped.iter().any(Option::is_none) {
                anyhow::bail!(
                    "these scoring options have no Karlin-Altschul values, and the first query batch (up to {} residues) has an invalid query; NCBI BLAST+ does not report the error then, which is not supported by LOSAT's BLASTN",
                    batch_size
                );
            }
            for (format, &output_format) in outputs.formats.iter_mut().zip(output_formats) {
                if output_format == BlastnOutputFormat::Pairwise {
                    let mut writer = format.sink.open()?;
                    write_blastn_pairwise_prolog(
                        &mut writer,
                        NCBI_BLASTN_VERSION,
                        megablast,
                        &subject_title,
                        seq_data.db_num_seqs,
                        seq_data.db_len_total,
                    )?;
                    writer.flush()?;
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
    // Gap costs that the tables accept but LOSAT's greedy extension does not (`scoring.rs`);
    // with no valid context, NCBI extends nothing.
    if any_valid_context {
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
    let db_len_total_i64 = seq_data.db_len_total as i64;
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
                BatchStage::Search => 0,
                BatchStage::ChunkPrelim {
                    batch_eff_searchsp, ..
                } => batch_eff_searchsp[index],
            };
            if !valid || options_searchsp != 0 {
                return options_searchsp;
            }
            let query_len = ctx.seq.len() as i64;
            let result = compute_length_adjustment_ncbi(
                query_len,
                db_len_total_i64,
                db_num_seqs_i64,
                &blocks.gapped,
            );
            let length_adjustment = result.length_adjustment;
            let effective_db_length =
                (db_len_total_i64 - db_num_seqs_i64 * length_adjustment).max(1);
            let effective_query_length = (query_len - length_adjustment).max(1);
            effective_db_length * effective_query_length
        })
        .collect();

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
        BatchStage::Search if any_valid_context => search_query_chunks(
            args,
            &seq_data.queries,
            &seq_data.query_masks,
            &query_eff_searchsp,
            &query_context_offsets,
            subjects,
            outputs,
            output_formats,
            parallel_pool,
            batch_size,
            first_batch,
        )?,
        _ => None,
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
    let subjects_for_lookup: &[bio::io::fasta::Record] = subject_records.as_deref().unwrap_or(&[]);
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
        build_lookup_tables(
            &config,
            &args,
            &encoded_queries_blastna,
            &query_context_masks,
            &query_context_offsets,
            subjects_for_lookup,
            subject_packed_cache_ref,
            approx_table_entries,
        )
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
    let progress_bar = if use_parallel {
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
    let hitlist_size = match args.max_target_seqs {
        Some(max_target_seqs) if max_target_seqs > 0 => max_target_seqs,
        _ => args.hitlist_size,
    };
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
        usize::try_from(get_prelim_hitlist_size(hitlist_size, false, true))
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

    // Pass lookup tables and queries for use in the closure
    let two_stage_lookup_ref = lookup_tables.two_stage_lookup.as_ref();
    let pv_direct_lookup_ref = lookup_tables.pv_direct_lookup.as_ref();
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:131-156
    // ```c
    // typedef struct BlastNaLookupTable {
    //     ...
    // } BlastNaLookupTable;
    // ```
    let na_lookup_ref = lookup_tables.na_lookup.as_ref();
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
    let subject_ids_ref = Arc::clone(&subject_ids_arc);
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
    let db_len_total = seq_data.db_len_total;
    let db_num_seqs = seq_data.db_num_seqs;
    let evalue_threshold = args.evalue;
    let subject_besthit = args.subject_besthit;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
    // ```c
    // return ((hsp->num_ident * 100.0 <
    //         align_length * hit_options->percent_identity) ||
    //         align_length < hit_options->min_hit_length) ;
    // ```
    let percent_identity = args.percent_identity;
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
    // Traditional megablast retains the zero-initialized X-drop option. Its
    // raw X-drop is the subject's word cutoff in ParametersUpdate below.
    let x_dropoff_init: Vec<i32> = search_karlin
        .iter()
        .map(|blocks| {
            if megablast {
                0
            } else {
                ((super::super::constants::X_DROP_UNGAPPED as f64 * NCBIMATH_LN2)
                    / blocks.ungapped.lambda)
                    .ceil() as i32
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
    let good_init_extends = std::sync::atomic::AtomicU64::new(0);
    let process_subject = |s_idx: usize,
                           s_record: &bio::io::fasta::Record,
                           gap_scratch: &mut GapAlignScratch,
                           subject_scratch: &mut SubjectScratch,
                           prelim_source: Option<&Vec<PrelimHit>>,
                           subject_hits: &mut Option<Vec<BlastnHsp>>,
                           prelim_out: &mut Vec<PrelimHit>| {
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
        let s_seq_full = s_record.seq();
        let s_len_full = s_seq_full.len();
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
        let s_id = subject_ids_ref
            .get(s_idx)
            .map(|id| id.as_ref())
            .unwrap_or("unknown");

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
        let subject_masks = if lcase_masking {
            collect_lowercase_masks(s_seq_full)
        } else {
            Vec::new()
        };
        let subject_masked = !subject_masks.is_empty();
        let soft_ranges = if subject_masked {
            build_subject_seq_ranges_from_masks(subject_masks.as_slice(), s_len_full)
        } else {
            vec![(0i32, s_len_full as i32)]
        };

        if s_len_full < effective_word_size {
            return;
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
        // ```c
        // BlastSeqBlkSetSequence(subj, sequence.data.release(), ...);
        // ...
        // BlastSeqBlkSetCompressedSequence(subj,
        //                                  compressed_seq.data.release());
        // ```
        let mut s_seq_blastna_full_vec: Option<Vec<u8>> = None;
        let s_seq_blastna_full = if let Some(blastna_cache) = subject_blastna_cache_ref {
            if let Some(blastna) = blastna_cache.get(s_idx) {
                blastna.as_slice()
            } else {
                s_seq_blastna_full_vec = Some(encode_iupac_to_blastna(s_seq_full));
                s_seq_blastna_full_vec.as_ref().unwrap().as_slice()
            }
        } else {
            s_seq_blastna_full_vec = Some(encode_iupac_to_blastna(s_seq_full));
            s_seq_blastna_full_vec.as_ref().unwrap().as_slice()
        };
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
        let mut s_seq_packed_full_vec: Option<Vec<u8>> = None;
        let s_seq_packed_full = if let Some(packed_cache) = subject_packed_cache_ref {
            if let Some(packed) = packed_cache.get(s_idx) {
                packed.as_slice()
            } else {
                s_seq_packed_full_vec = Some(encode_subject_ncbi2na_packed(s_seq_full));
                s_seq_packed_full_vec.as_ref().unwrap().as_slice()
            }
        } else {
            s_seq_packed_full_vec = Some(encode_subject_ncbi2na_packed(s_seq_full));
            s_seq_packed_full_vec.as_ref().unwrap().as_slice()
        };

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
            let hit_saving_cutoff = cutoff_score_max_from_evalue(
                evalue_threshold,
                query_eff_searchsp[context_idx],
                &search_karlin_ref[context_idx].gapped,
            );
            // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:342-374
            // gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda);
            // new_cutoff = gap_trigger;
            // new_cutoff *= (Int4)sbp->scale_factor;
            // new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
            let gap_trigger = gap_trigger_raw_score(
                GAP_TRIGGER_BIT_SCORE_NUCL,
                &search_karlin_ref[context_idx].ungapped,
            );
            let cutoff = cutoff_score_for_ungapped_extension(gap_trigger, hit_saving_cutoff, 1.0);
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
            let s_seq_blastna = &s_seq_blastna_full[chunk.offset..chunk.offset + chunk.length];
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
            let mut interval_tree = BlastIntervalTree::new(
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
            let use_diag_hash = query_concat_length > 8000;
            let use_array_indexing = !use_diag_hash;
            let diag_window = TWO_HIT_WINDOW as i32;

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
                    if TWO_HIT_WINDOW > 0 {
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
                        diag_hash_insert_window(TWO_HIT_WINDOW, scan_range, word_length)
                    } else {
                        0
                    };

                    let two_hits = TWO_HIT_WINDOW > 0;

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
                    let mut process_offset_pairs =
                        |offset_pairs: &[OffsetPair],
                         offset_pairs_len: &mut usize,
                         scan_pause_ns: &mut u64,
                         track_scan_pause: bool| {
                            let hitsfound = *offset_pairs_len;
                            for idx in 0..hitsfound {
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
                                let q_seq = ctx.seq.as_slice();
                                let q_context_start = ctx.query_offset as usize;
                                let q_context_end = q_context_start + q_seq.len();
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
                                        ctx.seq.len(),
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
                                let (last_hit, hit_saved) = if use_array_indexing
                                    && diag_idx < diag_array_size
                                {
                                    let diag_entry = &hit_level_array[diag_idx];
                                    (diag_entry.last_hit, diag_entry.flag != 0)
                                } else if !use_array_indexing {
                                    let (level, _hit_len, hit_saved) =
                                        diag_hash.retrieve(diag as i32).unwrap_or((0, 0, false));
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
                                        ctx.seq.len(),
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
                                let window_end = last_hit + TWO_HIT_WINDOW as i32;
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
                                        let window_size = TWO_HIT_WINDOW;
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
                                        ctx.seq.len(),
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
                                        diag_hash.insert(
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
                                        encoded_query_concat_blastna.as_slice(),
                                        search_seq_packed,
                                        q_off,
                                        s_off,
                                        s_len,
                                        x_dropoff,
                                        &score_matrix,
                                    )
                                } else {
                                    extend_hit_ungapped_approx_ncbi(
                                        encoded_query_concat_blastna.as_slice(),
                                        &query_four_base,
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
                                debug_assert!(ungapped.q_start >= q_context_start);
                                let qs = ungapped.q_start - q_context_start;
                                let qe = qs + ungapped.length;
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
                                        ctx.seq.len(),
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
                                            q_off.saturating_sub(q_context_start),
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
                                        diag_hash.insert(
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

                                if debug_enabled {
                                    dbg_gapped_attempted += 1;
                                }

                                if in_window && debug_mode {
                                    eprintln!("[DEBUG WINDOW] Seed at q={}, s={} -> COLLECT UNGAPPED (score={}, len={})", q_pos_usize, s_pos, ungapped_score, qe - qs);
                                }

                                // NCBI architecture: Collect ungapped hit for batch processing
                                // NCBI reference: blast_gapalign.c:3824 - init_hsp_array is sorted by score DESCENDING
                                // Gapped extension will be done later in score order with containment check
                                ungapped_hits.push(UngappedHit {
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
                                });

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
                                    diag_hash.insert(
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
                    let scan_start = if timing_enabled {
                        Some(std::time::Instant::now())
                    } else {
                        None
                    };
                    let mut scan_pause_ns: u64 = 0;
                    let mb_scan_kind = select_mb_scan_kind(safe_k, scan_step, subject_masked);
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
                                let q_seq = ctx.seq.as_slice();
                                let q_seq_blastna = encoded_queries_blastna[context_idx].as_slice();
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
                                        ctx.seq.len(),
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
                                let (last_hit, hit_saved) = if use_array_indexing
                                    && diag_idx < diag_array_size
                                {
                                    let diag_entry = &hit_level_array[diag_idx];
                                    (diag_entry.last_hit, diag_entry.flag != 0)
                                } else if !use_array_indexing {
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:833-836
                                    // ```c
                                    // rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
                                    // if(!rc)  last_hit = 0;
                                    // ```
                                    let (level, _hit_len, hit_saved) =
                                        diag_hash.retrieve(diag as i32).unwrap_or((0, 0, false));
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
                                let two_hits = TWO_HIT_WINDOW > 0;
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
                                    diag_hash_insert_window(TWO_HIT_WINDOW, scan_range, safe_k)
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
                                let window_end = last_hit + TWO_HIT_WINDOW as i32;
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
                                        q_seq.len(),
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
                                        let window_size = TWO_HIT_WINDOW;
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
                                        q_seq.len(),
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
                                        diag_hash.insert(
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
                                        &query_four_base[query_context_offsets[context_idx] as usize
                                            ..query_context_offsets[context_idx] as usize
                                                + q_seq_blastna.len().saturating_sub(3)],
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
                                        ctx.seq.len(),
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
                                        diag_hash.insert(
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
                                ungapped_hits.push(UngappedHit {
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
                                });

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
                                    diag_hash.insert(
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
            ungapped_hits.sort_by(score_compare_ungapped_hits);

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

            // Process each ungapped hit in score order
            for (idx, uh) in ungapped_hits.iter().enumerate() {
                // NCBI reference: blast_gapalign.c:3886-4090 (no per-HSP logging)
                if verbose && idx % 100 == 0 {
                    eprintln!("[INFO] Gapped extension: {}/{}", idx, total_ungapped);
                }
                // Get query sequence for this hit (context-specific)
                let ctx = &query_contexts[uh.context_idx as usize];
                let q_seq = ctx.seq.as_slice();

                let q_seq_blastna = encoded_queries_blastna[uh.context_idx as usize].as_slice();
                let subject_len = s_seq_blastna.len();

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
                let s_seq_score = s_seq_packed;
                let s_seq_trace = s_seq_blastna;

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
                    query_length: ctx.seq.len() as i32,
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
                        ctx.seq.len(),
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
                let (prelim_qs, prelim_qe, prelim_ss, prelim_se, prelim_score, seed_qs, seed_ss) =
                    if use_dp {
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
                            let mut score_q_anchor = seed_qs.saturating_add(offset_adjustment);
                            let mut score_s_anchor = seed_ss.saturating_add(offset_adjustment);
                            if score_q_anchor > q_seq_blastna.len() || score_s_anchor > subject_len
                            {
                                score_q_anchor = score_q_anchor.saturating_sub(COMPRESSION_RATIO);
                                score_s_anchor = score_s_anchor.saturating_sub(COMPRESSION_RATIO);
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

                        let prelim = match greedy_gapped_alignment_score_only(
                            q_seq_blastna,
                            s_seq_score,
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
                            ctx.seq.len(),
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
                    query_length: ctx.seq.len() as i32,
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
                chunk_prelim_hits = purge_prelim_hits_with_common_endpoints(chunk_prelim_hits);
                sort_prelim_hits_by_score(&mut chunk_prelim_hits);
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
                        crate::utils::threading::report_stage(
                            "blastn",
                            "subject_chunks",
                            chunk_batch.len(),
                            chunk_batch.len() > 1,
                        );
                        let soft_ranges = split_state.soft_ranges.as_slice();
                        if chunk_batch.len() == 1 {
                            let chunk = &chunk_batch[0];
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
                                true,
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
                                    true,
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
                                true,
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

            if !combined_prelim_hits.is_empty() {
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
        let s_seq_blastna = s_seq_blastna_full;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:371-373
        // ```c
        // tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
        //                               0, subject_length + 1);
        // ```
        let mut interval_tree = BlastIntervalTree::new(
            0,
            (query_concat_length + 1) as i32,
            0,
            (s_len_full + 1) as i32,
        );
        let mut traceback_edit_script_lengths: Vec<u32> = Vec::new();
        let mut traceback_alignment_lengths: Vec<u32> = Vec::new();

        let prelim_hits = &mut combined_prelim_hits;
        if !prelim_hits.is_empty() {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:358-365
            // ```c
            // /* Make sure the HSPs in the HSP list are sorted by score, as they should be. */
            // ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
            // ```
            if let Some(timing) = timing_ref {
                BlastnTiming::record_count(&timing.traceback_prelim_hsps, prelim_hits.len() as u64);
            }
            let sort_prelim_start = if timing_enabled {
                Some(std::time::Instant::now())
            } else {
                None
            };
            sort_prelim_hits_by_score(prelim_hits);
            if let (Some(timing), Some(sort_prelim_start)) = (timing_ref, sort_prelim_start) {
                BlastnTiming::record_duration(&timing.traceback_sort_prelim_ns, sort_prelim_start);
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
                let q_seq_blastna = encoded_queries_blastna[prelim.context_idx as usize].as_slice();
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
            let mut speculative_scratch: Vec<_> = if speculative_traceback {
                // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:509-513
                // BLAST_GappedAlignmentWithTraceback(...);
                // Each independent scratch slot retains its intermediate results.
                (0..num_threads)
                    .map(|_| (GapAlignScratch::new(), Vec::new()))
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
            const SPECULATIVE_TRACEBACK_BATCH_SIZE: usize =
                if cfg!(all(target_arch = "wasm32", not(feature = "wasm-threads"))) {
                    8
                } else {
                    16
                };
            for prelim_index in 0..prelim_hits.len() {
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
                if speculative_traceback && prelim_index % SPECULATIVE_TRACEBACK_BATCH_SIZE == 0 {
                    let batch_end =
                        (prelim_index + SPECULATIVE_TRACEBACK_BATCH_SIZE).min(prelim_hits.len());
                    speculative_jobs.clear();
                    speculative_jobs.extend((prelim_index..batch_end).filter(|&index| {
                        let p = &prelim_hits[index];
                        let h = TreeHsp {
                            query_offset: p.prelim_qs as i32,
                            query_end: p.prelim_qe as i32,
                            subject_offset: p.prelim_ss as i32,
                            subject_end: p.prelim_se as i32,
                            score: p.prelim_score,
                            query_frame: p.query_frame,
                            query_length: query_contexts[p.context_idx as usize].seq.len() as i32,
                            query_context_offset: p.query_context_offset,
                            subject_frame_sign: 1,
                        };
                        interval_tree
                            .containing_hsp(&h, p.query_context_offset, min_diag_separation)
                            .is_none()
                    }));
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
                            |(slot, (scratch, results))| {
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
                                            Some((index, (prepared, result)))
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
                    for (_, results) in &mut speculative_scratch {
                        for (index, result) in results.drain(..) {
                            speculative_results[index - prelim_index] = Some(result);
                        }
                    }
                }

                let prelim = &prelim_hits[prelim_index];
                let ctx = &query_contexts[prelim.context_idx as usize];
                let q_seq_blastna = encoded_queries_blastna[prelim.context_idx as usize].as_slice();
                let q_seq_nomask_blastna =
                    encoded_queries_blastna[prelim.context_idx as usize].as_slice();

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
                    query_length: ctx.seq.len() as i32,
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
                if let (Some(timing), Some(tree_precheck_start)) = (timing_ref, tree_precheck_start)
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
                            ctx.seq.len(),
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
                        BlastnTiming::record_count(&timing.traceback_tree_precheck_skipped_hsps, 1);
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
                let precomputed = if speculative_traceback {
                    speculative_results[prelim_index % SPECULATIVE_TRACEBACK_BATCH_SIZE].take()
                } else {
                    None
                };
                #[cfg(not(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                )))]
                let precomputed: Option<((usize, usize, usize, usize, usize), _)> = None;
                let start_offsets_start = timing_enabled.then(std::time::Instant::now);
                let prepared = precomputed
                    .as_ref()
                    .map(|(prepared, _)| *prepared)
                    .or_else(|| prepare_traceback(prelim));
                if let (Some(timing), Some(start)) = (timing_ref, start_offsets_start) {
                    BlastnTiming::record_duration(&timing.traceback_start_offsets_ns, start);
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
                let adjusted_subject = &s_seq_blastna[start_shift..start_shift + adjusted_s_len];
                let x_drop_trace = x_drop_final;
                if let Some(timing) = timing_ref {
                    BlastnTiming::record_count(&timing.traceback_full_traceback_hsps, 1);
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
                        precomputed.map(|(_, result)| result).unwrap_or_else(|| {
                            extend_gapped_heuristic_with_traceback_with_scratch(
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
                            )
                        })
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
                    match greedy_gapped_alignment_with_traceback(
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
                        &mut subject_scratch.greedy_align_scratch,
                    ) {
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
                            final_qs, final_qe, final_ss, final_se, score, 0usize, 0usize, 0usize,
                            0usize, aln_len, edit_ops,
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
                if let (Some(timing), Some(alignment_start)) = (timing_ref, alignment_start) {
                    let alignment_elapsed_ns = alignment_start.elapsed().as_nanos() as u64;
                    BlastnTiming::record_ns(&timing.traceback_alignment_ns, alignment_elapsed_ns);
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
                            ctx.seq.len(),
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
                    traceback_alignment_lengths.push(aln_len.min(u32::MAX as usize) as u32);
                }

                if use_dp && hsp_test(matches, aln_len, percent_identity, min_hit_length) {
                    if let Some(timing) = timing_ref {
                        BlastnTiming::record_count(
                            &timing.traceback_deleted_identity_length_hsps,
                            1,
                        );
                    }
                    if let (Some(timing), Some(hsp_build_start)) = (timing_ref, hsp_build_start) {
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
                    query_length: ctx.seq.len() as i32,
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

                let query_length = queries[prelim.query_idx as usize].seq().len();
                let (hit_q_start, hit_q_end, hit_s_start, hit_s_end) = adjust_blastn_offsets(
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
                if let (Some(timing), Some(hsp_build_start)) = (timing_ref, hsp_build_start) {
                    BlastnTiming::record_duration(&timing.traceback_hsp_build_ns, hsp_build_start);
                }
            }

            prelim_hits.clear();
        }

        // =================================================================
        // NCBI POST-GAPPED PROCESSING
        // Reference: blast_traceback.c:633-692
        // =================================================================

        // Step 1: Extract hits
        let local_hits: Vec<BlastnHsp> = hits.drain(..).collect();

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
        if blastn_trace_enabled && blastn_trace::should_trace_subject("purge", None, s_idx, s_id) {
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
            let context_idx = (hit.q_idx as usize) * 2 + if hit.query_frame < 0 { 1 } else { 0 };
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
            let q_seq_blastna = encoded_queries_blastna[context_idx].as_slice();
            let q_seq_nomask_blastna = encoded_queries_blastna[context_idx].as_slice();
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
            if let (Some(timing), Some(identity_length_start)) = (timing_ref, identity_length_start)
            {
                BlastnTiming::record_duration(
                    &timing.identity_length_test_ns,
                    identity_length_start,
                );
            }
            if delete_identity_length {
                if let Some(timing) = timing_ref {
                    BlastnTiming::record_count(&timing.traceback_deleted_identity_length_hsps, 1);
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
        if blastn_trace_enabled && blastn_trace::should_trace_subject("purge", None, s_idx, s_id) {
            blastn_trace::log(
                "purge",
                format!(
                    "subject={}({}) pass=post_reeval before={} after_reeval={} after_delete_purge={}",
                    s_id, s_idx, hits_step2, hits_step3, hits_step4
                ),
            );
        }

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
            BlastnTiming::record_duration(&timing.traceback_score_resort_ns, score_resort_start);
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
            let final_tree_containing_hsp =
                interval_tree.containing_hsp(&tree_hsp, query_context_offset, min_diag_separation);
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
                interval_tree.add_hsp(tree_hsp, query_context_offset, IndexMethod::QueryAndSubject);
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
            let context_idx = (hit.q_idx as usize) * 2 + if hit.query_frame < 0 { 1 } else { 0 };
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
    };
    let subject_records_ref = subjects.records;
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
    let search_subjects = |kept: Option<&[Vec<PrelimHit>]>| {
        let results: Vec<(Option<Vec<BlastnHsp>>, Vec<PrelimHit>)> = if use_parallel {
            #[cfg(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))]
            {
                parallel_pool.install(|| {
                    subject_records_ref
                        .par_iter()
                        .enumerate()
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
                                    SubjectScratch::new(queries_ref.len(), offset_array_size),
                                )
                            },
                            |state, (s_idx, s_record)| {
                                let (gap_scratch, subject_scratch) = state;
                                let mut subject_hits: Option<Vec<BlastnHsp>> = None;
                                let mut prelim_hits: Vec<PrelimHit> = Vec::new();
                                process_subject(
                                    s_idx,
                                    s_record,
                                    gap_scratch,
                                    subject_scratch,
                                    kept.map(|kept| &kept[s_idx]),
                                    &mut subject_hits,
                                    &mut prelim_hits,
                                );
                                (subject_hits, prelim_hits)
                            },
                        )
                        .collect()
                })
            }
            #[cfg(any(
                not(feature = "parallel"),
                all(target_arch = "wasm32", not(feature = "wasm-threads"))
            ))]
            {
                unreachable!("use_parallel is false when parallel threads are disabled");
            }
        } else {
            let mut gap_scratch = GapAlignScratch::new();
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:991-1041
            // ```c
            // Int4 offset_array_size = GetOffsetArraySize(lookup_wrap);
            // ...
            // aux_struct->offset_pairs =
            //   (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
            // ```
            let mut subject_scratch = SubjectScratch::new(queries_ref.len(), offset_array_size);
            subject_records_ref
                .iter()
                .enumerate()
                .map(|(s_idx, s_record)| {
                    let mut subject_hits: Option<Vec<BlastnHsp>> = None;
                    let mut prelim_hits: Vec<PrelimHit> = Vec::new();
                    process_subject(
                        s_idx,
                        s_record,
                        &mut gap_scratch,
                        &mut subject_scratch,
                        kept.map(|kept| &kept[s_idx]),
                        &mut subject_hits,
                        &mut prelim_hits,
                    );
                    (subject_hits, prelim_hits)
                })
                .collect()
        };
        let mut hit_lists: Vec<Option<BlastnHitList>> = Vec::with_capacity(queries_ref.len());
        hit_lists.resize_with(queries_ref.len(), || None);
        let mut prelim_hits = Vec::with_capacity(results.len());
        for (hits, prelim) in results {
            if let Some(hits) = hits {
                update_hitlists_with_subject_hits(&mut hit_lists, hits, hitlist_size);
            }
            prelim_hits.push(prelim);
        }
        (hit_lists, prelim_hits)
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
    let prelim_lists = match split_prelim_lists {
        Some(lists) => lists,
        None => {
            let (_, prelim_hits) = search_subjects(None);
            collect_prelim_hit_lists(prelim_hits, queries_ref.len(), prelim_hitlist_size)
        }
    };
    if let BatchStage::ChunkPrelim { .. } = stage {
        return Ok(batch_result(Vec::new(), prelim_lists));
    }
    let kept = prelim_hits_by_subject(prelim_lists, subject_records_ref.len());
    let (hit_lists, _) = search_subjects(Some(&kept));

    if let Some(bar) = progress_bar.as_ref() {
        bar.finish();
    }
    if let Some(timing) = timing.as_ref() {
        print_blastn_timing(timing.as_ref(), t_search_start, t_total);
    }
    Ok(batch_result(hit_lists, Vec::new()))
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
    queries: &[bio::io::fasta::Record],
    query_masks: &[Vec<MaskedInterval>],
    batch_eff_searchsp: &[i64],
    context_offsets: &[i32],
    subjects: &PreparedSubjects<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batch_size: i32,
    first_batch: bool,
) -> Result<Option<PrelimHitLists>> {
    let lengths: Vec<usize> = queries.iter().map(|record| record.seq().len()).collect();
    let Some(chunks) = split_query_batch(&lengths, args.task == "megablast") else {
        return Ok(None);
    };
    let mut merged: PrelimHitLists = Vec::with_capacity(queries.len());
    merged.resize_with(queries.len(), || None);
    for chunk in &chunks {
        // NCBI fails to make the query factory of a chunk without a query ("Empty
        // CBlastQueryVector"); `split_query_batch` makes none.
        anyhow::ensure!(
            !chunk.queries.is_empty(),
            "a query chunk without a query is not supported by LOSAT's BLASTN"
        );
        let parts: Vec<bio::io::fasta::Record> = chunk
            .queries
            .iter()
            .map(|part| {
                let record = &queries[part.query];
                bio::io::fasta::Record::with_attrs(
                    record.id(),
                    record.desc(),
                    &record.seq()[part.from..part.to],
                )
            })
            .collect();
        let restricted = chunk
            .queries
            .iter()
            .map(|part| restrict_masks(&query_masks[part.query], part))
            .collect();
        let part_masks = chunk_query_masks(args, &parts, restricted);
        let batch = search_query_batch(
            args,
            &parts,
            subjects,
            outputs,
            output_formats,
            parallel_pool,
            batch_size,
            first_batch,
            BatchStage::ChunkPrelim {
                query_masks: &part_masks,
                batch_eff_searchsp,
            },
        )?;
        merge_query_chunk(
            &mut merged,
            batch.chunk_prelim_lists,
            chunk,
            &lengths,
            context_offsets,
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
                let context = chunk.contexts[local_context];
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
        let split_points = [
            chunk.context_offsets[2 * local_query],
            chunk.context_offsets[2 * local_query + 1],
        ];
        merge_prelim_hit_list(hit_list, &mut merged[part.query], split_points);
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
                    QUERY_CHUNK_OVERLAP,
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
            let read_records = |fasta: &str| {
                bio::io::fasta::Reader::new(fasta.as_bytes())
                    .records()
                    .collect::<std::result::Result<Vec<_>, _>>()
                    .unwrap()
            };
            let mut output = Vec::new();
            let mut stderr = std::io::stderr();
            run_local(
                options.blastn,
                &read_records(&query),
                &read_records(&subject),
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
            run(args).unwrap();
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
            &read_blastn_fasta_records(&query_path, "query").unwrap()[0].seq(),
        );
        let subj = encode_iupac_to_blastna(
            &read_blastn_fasta_records(&subject_path, "subject").unwrap()[0].seq(),
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
        merge_prelim_hit_list(hitlist1, &mut combined, [0, 0]);
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
        merge_prelim_hit_list(single, &mut empty, [0, 0]);
        assert_eq!(prelim_scores(empty.as_ref().unwrap()), vec![(5, vec![10])]);
    }
}
